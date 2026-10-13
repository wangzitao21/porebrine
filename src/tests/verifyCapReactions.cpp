// SPDX-License-Identifier: GPL-3.0-or-later
#include "CapReactions.hpp"
#include <iostream>
#include <iomanip>
using namespace saltchem;
using namespace saltgeom;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void close(double a,double b,double relative,const char* message){require(std::abs(a-b)<=relative*std::max({std::abs(a),std::abs(b),1e-100}),message);}
std::array<long double,8> conserved(const Inventory& state){
    std::array<long double,8> out{};for(int j=0;j<7;j++)out[j]=state.aqueous[j];out[7]=state.waterEquivalent;
    for(int j=0;j<7;j++){
        for(int k=0;k<7;k++)out[k]+=static_cast<long double>(state.mineral[j])*stoichiometry[j][k];
        out[7]+=static_cast<long double>(state.mineral[j])*hydration[j]*waterMolarMass;
    }return out;
}
void balance(const Inventory& before,const Inventory& after){
    const auto a=conserved(before),b=conserved(after);
    for(int j=0;j<8;j++)require(std::abs(a[j]-b[j])<=5e-14*std::max(std::abs(a[j]),1.L),"Independent cap reactions lost conserved mass");
    require(after.waterEquivalent>0,"Cap reactions exhausted the shared water inventory");
    for(double amount:after.aqueous)require(amount>=0,"Cap reactions exhausted a shared reagent twice");
}
int main(int argc,char** argv){try{
    require(argc==2,"Supply pitzer.dat");PitzerModel model(argv[1]);Speciation spec(model);Reactions reaction(model);
    constexpr double cellVolume=1e-15,T=298.15,rCritical=4e-8;
    auto liquid=spec.inventory({{6.3,0,0,0,6.3,0,0}},1000,{});
    double water;const auto chemistry=spec.conserved(liquid,water);
    const double logOmega=ln10*chemistry.saturationIndex[0];require(logOmega>0,"Ostwald fixture must be bulk-supersaturated");
    const double gamma=rCritical*molarGasConstant*T*logOmega/(2*model.molarVolume[0]);
    MineralVector bare{};bare[1]=bare[2]=1e-10;
    RateParameters parameters;parameters.geometry.enabled=true;parameters.geometry.inertFraction=.1;parameters.neutral[0]=.02;
    auto population=[&](double energy){SurfaceState out;for(int k=1;k<=2;k++){
        auto& p=out.population[k];const double r=k==1?rCritical/2:2*rCritical;
        p.count=1;p.fraction=capVolumeFactor(1)*r*r*r/cellVolume;p.supportArea=bare[k];p.interfacialEnergy=energy;
    }return out;};
    IntegrationOptions options;options.relativeTolerance=1e-9;options.absoluteTolerance=1e-16;
    std::cout<<std::setprecision(16)<<"test,value\n";
    {
        auto cell=liquid;cell.mineral[2]=.1/model.molarVolume[2];SurfaceState seed;
        NucleationLaw law{0,2,1e13,0,60,2*rCritical};law.physicalCNT=true;law.interfacialEnergy=gamma;
        require(addEmbryo(cell,seed,law,cellVolume,model.molarVolume,1e-5,parameters.geometry.inertFraction,bare[2],nullptr,chemistry.saturationIndex[0],T),
            "Postcritical physical nucleus could not be inserted");
        const double r=radius(seed.population[2],cellVolume);
        const auto updatedChemistry=spec.conserved(cell,water);
        const double effective=ln10*updatedChemistry.saturationIndex[0]-2*gamma*model.molarVolume[0]/(molarGasConstant*T*r);
        require(r>rCritical&&effective>0,"Postcritical fixture became subcritical after consuming its finite reservoir");
        const auto initial=cell;const auto measured=measure(seed,bare,cellVolume);
        const auto result=advanceCapPopulations(reaction,cell,seed,parameters,measured,model.molarVolume,cellVolume,.001,T,options);
        require(result.fractionChange[2]>0,"Inserted postcritical CNT nucleus failed to grow under its Gibbs driving force");
        balance(initial,cell);std::cout<<"cnt_postcritical_birth_grows_with_kelvin_law,1\n";
    }
    for(double energy:{0.,gamma}){
        Inventory answers[3];
        for(int method=0;method<
#ifdef SALT_HAS_CVODE
            3
#else
            2
#endif
            ;method++){
            auto caps=population(energy);auto cell=liquid;
            cell.mineral[0]=(.2+subgrid(caps)[0])/model.molarVolume[0];const auto oldCell=cell;
            const double resolved=cell.mineral[0]-subgrid(caps)[0]/model.molarVolume[0];
            const auto surface=measure(caps,bare,cellVolume);options.method=method;
            const auto result=advanceCapPopulations(reaction,cell,caps,parameters,surface,model.molarVolume,cellVolume,.001,T,options);
            balance(oldCell,cell);require(result.fractionChange[2]>0,"Supercritical cap failed to grow");
            require(energy>0?result.fractionChange[1]<0:result.fractionChange[1]>0,"Kelvin/legacy cap growth has the wrong sign");
            close(cell.mineral[0]-subgrid(caps)[0]/model.molarVolume[0],resolved,3e-14,"Cap reactions consumed resolved primary stock");
            answers[method]=cell;
            if(method)for(int j=0;j<7;j++)close(cell.aqueous[j],answers[0].aqueous[j],2e-10,"Cap integrators disagree on shared aqueous state");
            std::cout<<(energy>0?"kelvin_":"legacy_")<<"method_"<<method<<"_small_delta,"<<result.fractionChange[1]<<'\n';
        }
    }
    {
        auto caps=population(gamma);auto cell=liquid;cell.mineral[0]=subgrid(caps)[0]/model.molarVolume[0];const auto original=cell;
        const auto measured=measure(caps,bare,cellVolume);
        const auto result=advanceCapPopulations(reaction,cell,caps,parameters,measured,model.molarVolume,cellVolume,1.,T,options);
        require(caps.population[1].fraction==0&&caps.population[1].count==0,"Subcritical cap failed to dissolve completely");
        require(result.fractionChange[1]<0&&result.fractionChange[2]>0,"Independent populations failed to represent Ostwald directions");
        balance(original,cell);std::cout<<"ostwald_depletion_and_shared_pool,1\n";
    }
    {
        for(double concentration:{.3,6.3}){
            auto caps=population(0);auto cell=spec.inventory({{concentration,0,0,0,concentration,0,0}},1000,{});
            cell.mineral[0]=subgrid(caps)[0]/model.molarVolume[0];const auto initial=cell;
            for(int step=0;step<1000;step++){
                const auto measured=measure(caps,bare,cellVolume);
                advanceCapPopulations(reaction,cell,caps,parameters,measured,model.molarVolume,cellVolume,.001,T,options);
                close(subgrid(caps)[0],cell.mineral[0]*model.molarVolume[0],4*std::numeric_limits<double>::epsilon(),
                    "Pure subgrid fraction drift accumulated across chemical steps");
            }
            balance(initial,cell);
            if(concentration<1)require(cell.mineral[0]==0&&subgrid(caps)[0]==0,"Complete cap dissolution left a fake resolved sliver");
        }
        std::cout<<"thousand_step_pure_cap_geometry_alignment,1\n";
    }
    {
        // The same thermal/curvature expression is used by doubles and AD.
        auto caps=population(gamma);auto cell=liquid;cell.mineral[0]=caps.population[2].fraction/model.molarVolume[0];
        auto p=parameters;p.geometry.separateMineralAreas=true;p.area[0]=p.geometry.growthArea[0]=outerArea(caps.population[2],cellVolume)/cellVolume;
        p.referenceMoles[0]=cell.mineral[0];p.geometry.surfaceExponent[0]=2./3;
        for(double curvature:{0.,logOmega*.5,logOmega*2}){
            p.geometry.curvatureLogOmega[0]=curvature;
            const auto value=reaction.evaluate(cell,p);double jac[7][7];reaction.rateJacobian(cell,p,value,jac);
            const double delta=cell.mineral[0]*1e-4;MineralVector plus{},minus{};plus[0]=delta;minus[0]=-delta;
            const double numeric=(reaction.evaluate(Reactions::changed(cell,plus),p).rates[0]-reaction.evaluate(Reactions::changed(cell,minus),p).rates[0])/(2*delta);
            close(jac[0][0],numeric,3e-5,"Curvature analytic Jacobian disagrees with finite differences");
        }
        p.geometry.curvatureLogOmega[0]=logOmega;
        require(std::abs(reaction.evaluate(cell,p).rates[0])<1e-12,"Critical cap does not satisfy the Gibbs-Thomson stationary condition");
        cell.mineral[0]=0;require(reaction.evaluate(cell,p).rates[0]==0,"Exhausted curved cap regenerated without nucleation");
        std::cout<<"gibbs_thomson_stationarity_and_analytic_jacobian,1\n";
    }
    {
        auto quiet=spec.inventory({{.3,0,0,0,.3,0,0}},731.123456789,{});
        quiet.mineral[0]=.0123456789876543;const auto unchanged=quiet;
        auto growth=parameters;growth.geometry.separateMineralAreas=true;growth.area.fill(0);growth.geometry.growthArea[0]=1234;
        for(int method=0;method<
#ifdef SALT_HAS_CVODE
            3
#else
            2
#endif
            ;method++){
            options.method=method;
            for(int repetition=0;repetition<200;repetition++)reaction.advance(quiet,.01,growth,options);
            require(quiet.aqueous==unchanged.aqueous&&quiet.mineral==unchanged.mineral&&quiet.waterEquivalent==unchanged.waterEquivalent,
                "Inactive growth-only chemistry drifts under repeated normalization");
        }
        std::cout<<"growth_only_exact_noop_all_integrators,1\n";
        quiet=unchanged;quiet.mineral[0]=.1;quiet.mineral[4]=1.23456789876543;const auto initial=quiet;
        growth.area[0]=growth.geometry.growthArea[0]=10;
        growth.neutral[4]=.01;growth.geometry.growthArea[4]=1000;
        const auto statistics=reaction.advance(quiet,.01,growth,options);
        require(statistics.dissolved[0]>0&&statistics.dissolved[4]==0&&statistics.precipitated[4]==0,"Mixed active/inactive fixture did not isolate the growth-only phase");
        require(quiet.mineral[4]==initial.mineral[4]&&quiet.aqueous[Ca]==0&&quiet.aqueous[S]==0,
            "Inactive growth-only stock drifted or created phantom aqueous components in an active solve");
        balance(initial,quiet);std::cout<<"inactive_growth_stock_preserved_during_other_reactions,1\n";
    }
    {
        auto caps=population(0),untouched=caps;auto cell=liquid;
        cell.mineral[0]=(.2+subgrid(caps)[0])/model.molarVolume[0];const auto initial=cell;
        auto measured=measure(caps,bare,cellVolume);measured.populationArea[2]=-1;
        bool failed=false;try{advanceCapPopulations(reaction,cell,caps,parameters,measured,model.molarVolume,cellVolume,.001,T,options);}catch(const std::invalid_argument&){failed=true;}
        require(failed&&cell.aqueous==initial.aqueous&&cell.mineral==initial.mineral&&cell.waterEquivalent==initial.waterEquivalent,
            "Failed second cap left a partially committed liquid inventory");
        require(caps.population[1].fraction==untouched.population[1].fraction,"Failed cap batch partially changed its first population");
        std::cout<<"whole_cell_atomic_cap_commit,1\nall_passed,1\n";
    }
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
