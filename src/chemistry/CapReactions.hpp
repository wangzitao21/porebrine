// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "BrineReactions.hpp"
#include "MineralInterface.hpp"

namespace saltchem {
struct CapReactionResult {
    IntegrationStats stats;
    std::array<double,49> fractionChange{}; // Positive for cap growth, per representative bulk volume.
};
inline CapReactionResult advanceCapPopulations(
    const Reactions& reactions,Inventory& cell,saltgeom::SurfaceState& surface,
    const RateParameters& input,const saltgeom::SurfaceMeasure& measured,
    const MineralVector& molarVolume,double cellVolume,double interval,double temperatureK,
    const IntegrationOptions& options={}){
    if(!input.geometry.enabled||!std::isfinite(input.geometry.bulkVolume)||input.geometry.bulkVolume<=0||
       !std::isfinite(cellVolume)||cellVolume<=0||!std::isfinite(interval)||interval<0||
       !std::isfinite(temperatureK)||temperatureK<=0)
        throw std::invalid_argument("Invalid cap-reaction geometry, interval or temperature");
    Reactions::validateParameters(input);
    if(!std::isfinite(cell.waterEquivalent)||cell.waterEquivalent<=0)throw std::invalid_argument("Invalid cap-reaction water inventory");
    for(double value:cell.aqueous)if(!std::isfinite(value)||value<0)throw std::invalid_argument("Invalid cap-reaction aqueous inventory");
    for(double value:cell.mineral)if(!std::isfinite(value)||value<0)throw std::invalid_argument("Invalid cap-reaction mineral inventory");
    const double bulkVolume=input.geometry.bulkVolume,epsilon=std::numeric_limits<double>::epsilon();
    const auto caps=saltgeom::subgrid(surface);
    std::array<bool,nMinerals> pureCaps{};
    for(int j=0;j<nMinerals;j++){
        if(!std::isfinite(molarVolume[j])||molarVolume[j]<=0)throw std::invalid_argument("Invalid cap mineral molar volume");
        const double amount=caps[j]*bulkVolume/molarVolume[j];
        pureCaps[j]=saltgeom::whollySubgrid(cell.mineral[j],amount);
        if(amount>cell.mineral[j]+64*epsilon*std::max(amount,cell.mineral[j])){
            std::ostringstream message;message<<std::setprecision(17)
                <<"Cap populations exceed the total mineral inventory: mineral="<<mineralNames[j]
                <<" caps_mol="<<amount<<" total_mol="<<cell.mineral[j]<<" excess_mol="<<amount-cell.mineral[j]
                <<" roundoff_budget_mol="<<64*epsilon*std::max(amount,cell.mineral[j])
                <<" cap_fraction="<<caps[j]<<" bulk_volume="<<bulkVolume;
            for(int k=0;k<7;k++)if(surface.population[7*j+k].fraction>0)
                message<<" population["<<k<<"]="<<surface.population[7*j+k].fraction<<" count="<<surface.population[7*j+k].count;
            throw std::invalid_argument(message.str());
        }
    }
    Inventory candidate=cell;auto populations=surface;CapReactionResult result;
    for(int pair=0;pair<49;pair++){
        auto& cap=populations.population[pair];const int mineral=pair/7;
        if(!std::isfinite(cap.fraction)||cap.fraction<0||!std::isfinite(cap.count)||cap.count<0||
           !std::isfinite(cap.heightFactor)||cap.heightFactor<=0||cap.heightFactor>2||
           !std::isfinite(cap.interfacialEnergy)||cap.interfacialEnergy<0||
           !std::isfinite(measured.populationArea[pair])||measured.populationArea[pair]<0)
            throw std::invalid_argument("Invalid cap population or reacting area");
        if(cap.fraction>0&&cap.count==0)throw std::invalid_argument("Nonzero cap volume without a population count");
        if(cap.fraction==0||measured.populationArea[pair]==0||interval==0||input.geometry.weights[mineral]==0||
           (input.neutral[mineral]==0&&input.acid[mineral]==0&&input.co2[mineral]==0))continue;
        const double amount=cap.fraction*bulkVolume/molarVolume[mineral];
        Inventory reactor=candidate;reactor.mineral.fill(0);reactor.mineral[mineral]=amount;
        RateParameters parameters=input;parameters.geometry.separateMineralAreas=true;
        parameters.area.fill(0);parameters.geometry.growthArea.fill(0);
        parameters.geometry.lockedMoles.fill(0);parameters.geometry.curvatureLogOmega.fill(0);
        parameters.geometry.surfaceExponent.fill(0);parameters.referenceMoles.fill(1);
        parameters.area[mineral]=parameters.geometry.growthArea[mineral]=measured.populationArea[pair]*bulkVolume/cellVolume;
        parameters.referenceMoles[mineral]=amount;parameters.geometry.surfaceExponent[mineral]=2./3;
        const double curvatureRadius=saltgeom::radius(cap,cellVolume);
        if(!(curvatureRadius>0)||!std::isfinite(curvatureRadius))throw std::invalid_argument("Invalid active cap curvature radius");
        if(cap.interfacialEnergy>0)parameters.geometry.curvatureLogOmega[mineral]=2*cap.interfacialEnergy*molarVolume[mineral]/
            (saltgeom::molarGasConstant*temperatureK*curvatureRadius);
        // Preserve the exact small active stock instead of subtracting two
        // nearly equal large mineral totals to recover a nucleus inventory.
        // Other caps and resolved minerals contribute only inert pore volume.
        long double inert=input.geometry.inertFraction;
        for(int j=0;j<nMinerals;j++)inert+=static_cast<long double>(candidate.mineral[j])*molarVolume[j]/bulkVolume;
        inert-=static_cast<long double>(amount)*molarVolume[mineral]/bulkVolume;
        if(inert<0&&inert>=-64*epsilon)inert=0;
        parameters.geometry.inertFraction=double(inert);
        parameters.geometry.initialPorosity=poreFraction(candidate,input.geometry,molarVolume);
        const auto stats=reactions.advance(reactor,interval,parameters,options);
        double extent=amount-reactor.mineral[mineral];
        const double activeRoundoff=64*epsilon*std::max(amount,reactor.mineral[mineral]);
        if(reactor.mineral[mineral]==0&&pureCaps[mineral]){
            bool last=true;for(int k=0;k<7;k++)if(7*mineral+k!=pair&&populations.population[7*mineral+k].fraction>0)last=false;
            // An exact final depletion consumes the actual remaining stock,
            // including its conversion ulp, with the matching aqueous/water
            // transfer. It must not leave a fake resolved mineral sliver.
            if(last){
                if(std::abs(extent-candidate.mineral[mineral])>activeRoundoff)
                    throw std::runtime_error("Final cap depletion differs from remaining mineral beyond roundoff");
                extent=candidate.mineral[mineral];
            }
        }
        if(extent>candidate.mineral[mineral]){
            if(extent-candidate.mineral[mineral]>activeRoundoff)throw std::runtime_error("Cap dissolution exceeds the total mineral inventory");
            extent=candidate.mineral[mineral];
        }
        if(extent<0){
            // Roundoff in the returned active stock must not turn an exactly
            // exhausted dissolved reagent into a negative inventory.
            double available=std::numeric_limits<double>::infinity();
            for(int j=0;j<nComponents;j++)if(stoichiometry[mineral][j]>0)
                available=std::min(available,candidate.aqueous[j]/stoichiometry[mineral][j]);
            if(hydration[mineral]>0)available=std::min(available,std::nextafter(candidate.waterEquivalent/(hydration[mineral]*waterMolarMass),0.));
            if(-extent>available){
                if(-extent-available>activeRoundoff)throw std::runtime_error("Cap growth exceeds shared aqueous inventory");
                extent=-available;
            }
        }
        MineralVector transfer{};transfer[mineral]=extent;
        const auto updated=Reactions::changed(candidate,transfer);
        const double represented=candidate.mineral[mineral]-updated.mineral[mineral];
        const double totalRoundoff=64*epsilon*std::max(candidate.mineral[mineral],updated.mineral[mineral]);
        if(std::abs(represented-extent)>totalRoundoff||updated.mineral[mineral]<0||updated.waterEquivalent<=0)
            throw std::runtime_error("Cap/total mineral transfer exceeds the representation roundoff budget");
        for(double value:updated.aqueous)if(value<0||!std::isfinite(value))throw std::runtime_error("Cap transfer produced an invalid aqueous inventory");
        const double oldFraction=cap.fraction;
        cap.fraction=(reactor.mineral[mineral]==0?0:std::max(0.,amount-extent))*molarVolume[mineral]/bulkVolume;
        result.fractionChange[pair]=-extent*molarVolume[mineral]/bulkVolume;
        if(std::abs((cap.fraction-oldFraction)-result.fractionChange[pair])>64*epsilon*std::max(oldFraction,cap.fraction))
            throw std::runtime_error("Cap volume transfer exceeds the representation roundoff budget");
        if(cap.fraction==0)saltgeom::clearPopulationMaterial(cap);
        candidate=updated;
        if(pureCaps[mineral])saltgeom::alignPureSubgridInventory(populations,mineral,candidate.mineral[mineral]*molarVolume[mineral]/bulkVolume);
        for(auto field:{&IntegrationStats::accepted,&IntegrationStats::rejected,&IntegrationStats::rateEvaluations,
            &IntegrationStats::speciationIterations,&IntegrationStats::newtonIterations,&IntegrationStats::depletionEvents,
            &IntegrationStats::historyRestarts,&IntegrationStats::integratorFallbacks,&IntegrationStats::analyticJacobians,&IntegrationStats::jacobianFallbacks})
            result.stats.*field+=stats.*field;
        result.stats.dissolved[mineral]+=std::max(extent,0.);result.stats.precipitated[mineral]+=std::max(-extent,0.);
    }
    cell=candidate;surface=populations;return result;
}
}
