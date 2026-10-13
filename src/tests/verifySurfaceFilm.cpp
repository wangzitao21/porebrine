// SPDX-License-Identifier: GPL-3.0-or-later
#include "SurfaceFilm.hpp"
#include "BrineInlet.hpp"
#include <iostream>
#include <iomanip>
#include <numeric>
#include <cfenv>

using namespace saltchem;
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
void close(double a,double b,double relative,const char* message){
    require(std::abs(a-b)<=relative*std::max({std::abs(a),std::abs(b),1e-100}),message);
}
void verify(const Reactions& reactions,const Inventory& bulk,int mineral,double k,double G,const FilmDissolution& film){
    require(film.surfaceStateAvailable&&film.rate>0&&film.rate<=k,"Film rate or surface state is invalid");
    require(std::abs(film.residual)<3e-10*k,"Film scalar residual is too large");
    close(std::accumulate(film.surfaceMassFractions.begin(),film.surfaceMassFractions.end(),0.),1.,4e-15,"Film mass fractions do not sum to one");
    Inventory surface;surface.waterEquivalent=film.surfaceMassFractions[7];
    for(int i=0;i<7;i++){require(film.surfaceMassFractions[i]>=0,"Negative film mass fraction");surface.aqueous[i]=film.surfaceMassFractions[i]/componentMass[i];}
    RateParameters zero;const auto aqueous=reactions.evaluate(surface,zero,nullptr,false).aqueous;
    close(aqueous.waterActivity,film.waterActivity,2e-10,"Reported surface water activity does not match Pitzer");
    require(std::abs(film.rate-k*(-std::expm1(ln10*aqueous.saturationIndex[mineral])))<3e-10*k,"Independent surface speciation fails nonlinear rate equation");
    const double pe=mineralMass[mineral]*film.rate/G,bulkMass=aqueousMass(bulk);
    const double transport=pe==0?G:G*pe/(-std::expm1(-pe));
    // J_i=G*B(-Pe)*(Ys-exp(-Pe)*Yb), checked independently against
    // each dissolved-mineral component including the hydrate-water flux.
    for(int i=0;i<8;i++){
        const double y=i<7?componentMass[i]*bulk.aqueous[i]/bulkMass:bulk.waterEquivalent/bulkMass;
        const double stoichMass=i<7?componentMass[i]*stoichiometry[mineral][i]:hydration[mineral]*waterMolarMass;
        const double flux=transport*(film.surfaceMassFractions[i]-std::exp(-pe)*y);
        const double expected=stoichMass*film.rate;
        require(std::abs(flux-expected)<3e-12*std::max(G,mineralMass[mineral]*film.rate),"Film flux is not stoichiometric or lost hydrate water");
    }
}
int main(int argc,char** argv){try{
    #if defined(__GLIBC__)
    // Production OpenFOAM enables these traps: invalid bracket trials must
    // fail by an exception rather than aborting the coupled calculation.
    feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);
    #endif
    require(argc==2,"Supply pitzer.dat");PitzerModel model(argv[1]);Speciation spec(model);Reactions reactions(model);
    constexpr double k=.02;
    const auto dilute=spec.inventory({{.05,0,0,0,.05,0,0}},1.);
    const auto brine=spec.inventory({{4,.2,.05,.5,5.28,.01,.005}},1.);
    RateParameters zero;const auto bulk=reactions.evaluate(dilute,zero,nullptr,false).aqueous;
    const auto bare=solveFilmDissolution(reactions,dilute,Halite,k,std::numeric_limits<double>::infinity(),&bulk);
    close(bare.rate,k*(-std::expm1(ln10*bulk.saturationIndex[Halite])),2e-15,"Bare surface does not recover the native kinetic law");
    require(bare.evaluations==0,"Cached bare bulk state caused a redundant speciation");
    const auto sealed=solveFilmDissolution(reactions,dilute,Halite,k,0.,&bulk);
    require(sealed.rate==0&&!sealed.surfaceStateAvailable,"Sealed film reacted or claimed a uniquely known surface composition");
    auto super=spec.inventory({{7,0,0,0,7,0,0}},1.);
    require(solveFilmDissolution(reactions,super,Halite,k,.01).rate==0,"Supersaturated primary mineral dissolved through film");
    require(solveFilmDissolution(reactions,dilute,Halite,0,.01).rate==0,"Zero-rate primary mineral reacted");
    std::cout<<std::setprecision(16)<<"test,value\nbare_rate,"<<bare.rate<<"\nsealed_rate,"<<sealed.rate<<'\n';
    double previous=0;
    for(double G:{1e-10,1e-6,.001,.01,.1,1.,1e3}){
        const auto film=solveFilmDissolution(reactions,dilute,Halite,k,G,&bulk);verify(reactions,dilute,Halite,k,G,film);
        require(film.rate>previous&&film.rate<bare.rate,"Film rate does not increase monotonically with conductance");previous=film.rate;
        std::cout<<"halite_rate_G_"<<G<<','<<film.rate<<"\nfilm_evaluations,"<<film.evaluations<<'\n';
    }
    const auto covered=solveFilmDissolution(reactions,brine,Halite,k,.01);
    verify(reactions,brine,Halite,k,.01,covered);
    require(covered.rate>0,"Complete permeable coverage incorrectly caused complete passivation");
    const double rho=1200,D1=2e-11,D2=8e-10,d1=1e-6,d2=20e-6;
    const double resistance=d1/(rho*D1)+d2/(rho*D2),series=1/resistance;
    const auto twoLayers=solveFilmDissolution(reactions,brine,Halite,k,series);
    const auto oneLayer=solveFilmDissolution(reactions,brine,Halite,k,rho*D2/d2);
    verify(reactions,brine,Halite,k,series,twoLayers);
    require(twoLayers.rate<oneLayer.rate,"Dense inner layer did not add series diffusion resistance");
    const double theta=.7,parallel=(1-theta)*solveFilmDissolution(reactions,brine,Halite,k,std::numeric_limits<double>::infinity()).rate+theta*covered.rate;
    require(parallel>covered.rate&&parallel<k,"Bare/coated parallel area weighting is invalid");
    std::cout<<"concentrated_brine_rate,"<<covered.rate<<"\ntwo_layer_rate,"<<twoLayers.rate<<"\nparallel_patch_rate,"<<parallel<<'\n';
    for(int mineral:{Sylvite,Carnallite,Bischofite,Gypsum,Polyhalite,Calcite}){
        const auto film=solveFilmDissolution(reactions,dilute,mineral,k,.02);
        verify(reactions,dilute,mineral,k,.02,film);
        std::cout<<"mineral_rate_"<<mineralNames[mineral]<<','<<film.rate<<'\n';
    }
    {
        // Extreme but positive conductance must approach saturation with a
        // nonnegative rate, not accept an invalid no-water upper trial.
        const auto tiny=solveFilmDissolution(reactions,dilute,Halite,k,1e-280);
        require(tiny.rate>0&&tiny.surfaceMassFractions[7]>0&&std::abs(tiny.saturationIndex)<1e-8,"Small-G film accepted an invalid surface or missed equilibrium");
        const auto high=solveFilmDissolution(reactions,dilute,Halite,k,1e280);
        close(high.rate,bare.rate,3e-9,"Large-G film does not approach bare kinetics");
        bool rejected=false;try{solveFilmDissolution(reactions,dilute,Halite,k,-1);}catch(const std::invalid_argument&){rejected=true;}
        require(rejected,"Negative film conductance was accepted");
        auto noWater=dilute;noWater.waterEquivalent=0;rejected=false;
        try{solveFilmDissolution(reactions,noWater,Halite,k,.01);}catch(const std::invalid_argument&){rejected=true;}
        require(rejected,"No-water film bulk state was accepted");
        std::cout<<"extreme_conductance_and_input_guards,1\n";
    }
    for(double temperature:{278.15,313.15,363.15}){
        model.setConditions(temperature,10.);
        const auto film=solveFilmDissolution(reactions,brine,Carnallite,k,.02);verify(reactions,brine,Carnallite,k,.02,film);
        std::cout<<"temperature_K_"<<temperature<<','<<film.rate<<'\n';
    }
    model.setConditions(298.15,1.);
    {
        // The production coupling freezes the film factor over each local
        // interval while integrating the native bulk chemical ODE. Verify
        // convergence under refinement, rather than claiming implicit film
        // dynamics from a single nonlinear surface solve.
        auto integrate=[&](int steps){
            auto cell=spec.inventory({{.5,0,0,0,.5,0,0}},.01,{});cell.mineral[Halite]=1;
            RateParameters parameters;parameters.geometry.enabled=true;parameters.geometry.separateMineralAreas=true;
            parameters.geometry.bulkVolume=.001;parameters.area[Halite]=.1;parameters.geometry.growthArea[Halite]=0;
            IntegrationOptions options;options.relativeTolerance=2e-10;options.absoluteTolerance=1e-14;
            for(int step=0;step<steps;step++){
                const auto aqueous=reactions.evaluate(cell,zero,nullptr,false).aqueous;
                const auto film=solveFilmDissolution(reactions,cell,Halite,k,.01,&aqueous);
                const double affinity=-std::expm1(ln10*aqueous.saturationIndex[Halite]);
                parameters.neutral[Halite]=affinity>0?film.rate/affinity:0;
                reactions.advance(cell,2./steps,parameters,options);
            }return cell.aqueous[Na];
        };
        const double reference=integrate(256),a=std::abs(integrate(4)-reference),b=std::abs(integrate(8)-reference),c=std::abs(integrate(16)-reference);
        require(a>b&&b>c&&a/b>1.5&&b/c>1.5,"Lagged Pitzer film coupling failed first-order time refinement");
        std::cout<<"lagged_error_4,"<<a<<"\nlagged_error_8,"<<b<<"\nlagged_error_16,"<<c<<'\n';
    }
    {
        const auto supersaturated=spec.inventory({{7,0,0,0,7,0,0}},1.);
        const auto hydratedBrine=spec.inventory({{1,4,0,7,19,0,0}},1.);
        for(const auto& sample:std::array<std::pair<int,Inventory>,2>{{{Halite,supersaturated},{Carnallite,hydratedBrine}}}){
            const int mineral=sample.first;const auto& state=sample.second;
            const auto aqueous=reactions.evaluate(state,zero,nullptr,false).aqueous;
            require(aqueous.saturationIndex[mineral]>0,"Precipitation fixture is not supersaturated");
            const auto bare=solveFilmPrecipitation(reactions,state,mineral,k,INFINITY,&aqueous);
            require(bare.rate<0&&bare.evaluations==0,"Infinite-conductance precipitation film did not recover bare TST");
            const auto sealed=solveFilmPrecipitation(reactions,state,mineral,k,0,&aqueous);
            require(sealed.rate==0&&sealed.evaluations==0&&!sealed.surfaceStateAvailable,"Sealed precipitation film reacted");
            double previous=0;
            for(double G:{1e-10,1e-6,.001,.01,.1,1.,1e3}){
                const auto film=solveFilmPrecipitation(reactions,state,mineral,k,G,&aqueous);
                require(film.rate<0&&-film.rate>previous&&-film.rate<=-bare.rate*(1+1e-8),"Precipitation film failed the signed conductance limit");previous=-film.rate;
                require(film.peclet<0&&film.surfaceMassFractions[7]>0,"Precipitation film lost its negative Stefan flux or positive solvent");
                Inventory surface;surface.waterEquivalent=film.surfaceMassFractions[7];
                for(int i=0;i<7;i++){require(film.surfaceMassFractions[i]>=0,"Precipitation film accepted a negative component");surface.aqueous[i]=film.surfaceMassFractions[i]/componentMass[i];}
                const auto check=reactions.evaluate(surface,zero,nullptr,false).aqueous;
                const double expected=k*(-std::expm1(ln10*check.saturationIndex[mineral]));
                require(std::abs(film.rate-expected)<3e-10*std::max(k,-bare.rate),"Precipitation-film surface Pitzer affinity fails its signed rate equation");
                const double pe=mineralMass[mineral]*film.rate/G,transport=G*pe/(-std::expm1(-pe));
                for(int i=0;i<8;i++){
                    const double y=i<7?componentMass[i]*state.aqueous[i]/aqueousMass(state):state.waterEquivalent/aqueousMass(state);
                    const double stoich=i<7?componentMass[i]*stoichiometry[mineral][i]:hydration[mineral]*waterMolarMass;
                    const double flux=transport*(film.surfaceMassFractions[i]-std::exp(-pe)*y);
                    require(std::abs(flux-stoich*film.rate)<3e-12*std::max(G,-mineralMass[mineral]*film.rate),"Negative Stefan flux lost component or hydrate-water stoichiometry");
                }
                close(std::accumulate(film.surfaceMassFractions.begin(),film.surfaceMassFractions.end(),0.),1.,4e-15,"Precipitation-film fractions lost total mass");
                std::cout<<"precipitation_"<<mineralNames[mineral]<<"_G_"<<G<<','<<film.rate<<'\n';
            }
            const auto tiny=solveFilmPrecipitation(reactions,state,mineral,k,1e-280,&aqueous);
            require(tiny.rate<0&&tiny.surfaceMassFractions[7]>0&&std::abs(tiny.saturationIndex)<1e-8,"Small-conductance precipitation failed its equilibrium limit");
            const auto large=solveFilmPrecipitation(reactions,state,mineral,k,1e280,&aqueous);
            close(large.rate,bare.rate,3e-9,"Large-conductance precipitation failed its bare limit");
        }
        require(solveFilmPrecipitation(reactions,dilute,Halite,k,.01).rate==0,"Undersaturated film precipitated primary mineral");
        std::cout<<"signed_precipitation_stefan_and_hydrate_checks,1\n";
    }
    std::cout<<"all_passed,1\n";
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
