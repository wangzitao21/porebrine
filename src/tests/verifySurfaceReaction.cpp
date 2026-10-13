// SPDX-License-Identifier: GPL-3.0-or-later
#include "SurfaceReaction.hpp"
#include "BrineInlet.hpp"
#include <iostream>
#include <iomanip>
#include <cfenv>

using namespace saltchem;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void close(double a,double b,double relative,const char* message){
    require(std::abs(a-b)<=relative*std::max({std::abs(a),std::abs(b),1e-100}),message);
}
std::array<double,8> conserved(const Inventory& state){
    std::array<double,8> out{};for(int j=0;j<7;j++)out[j]=state.aqueous[j];out[7]=state.waterEquivalent;
    for(int j=0;j<7;j++){
        for(int k=0;k<7;k++)out[k]+=state.mineral[j]*stoichiometry[j][k];
        out[7]+=state.mineral[j]*hydration[j]*waterMolarMass;
    }return out;
}
void balance(const Inventory& oldSolid,const Inventory& oldLiquid,const Inventory& solid,const Inventory& liquid){
    const auto a=conserved(oldSolid),b=conserved(oldLiquid),c=conserved(solid),d=conserved(liquid);
    for(int j=0;j<8;j++)close(a[j]+b[j],c[j]+d[j],3e-14,"Surface transfer lost a conserved component or water");
    require(solid.aqueous==oldSolid.aqueous&&solid.waterEquivalent==oldSolid.waterEquivalent,"Surface transfer changed solid-owner aqueous inventory");
    require(liquid.mineral==oldLiquid.mineral,"Surface transfer changed liquid-owner minerals");
    require(liquid.waterEquivalent>0,"Surface transfer exhausted liquid water");
    for(double n:liquid.aqueous)require(n>=0,"Surface transfer produced a negative aqueous inventory");
}
Inventory scaled(Inventory state,double factor){
    for(double& n:state.aqueous)n*=factor;for(double& n:state.mineral)n*=factor;state.waterEquivalent*=factor;return state;
}
int main(int argc,char** argv){try{
    #if defined(__GLIBC__)
    feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);
    #endif
    require(argc==2,"Supply pitzer.dat");PitzerModel model(argv[1]);Speciation spec(model);Reactions reaction(model);
    constexpr double volume=1e-6,area=1e-4;
    RateParameters rates;rates.geometry.enabled=true;rates.geometry.bulkVolume=volume;rates.geometry.inertFraction=.1;
    rates.neutral[0]=rates.neutral[2]=.02;
    IntegrationOptions options;options.relativeTolerance=2e-10;options.absoluteTolerance=1e-14;
    Inventory solid;solid.mineral[0]=.01;solid.waterEquivalent=.27;solid.aqueous[0]=solid.aqueous[4]=.013;
    const auto low=spec.inventory({{.01,0,0,0,.01,0,0}},.001,{});
    auto high=spec.inventory({{6,0,0,0,6,0,0}},.001,{});high.mineral[6]=1e-10;
    MineralVector face{},locked{};face[0]=area;locked[0]=.005;
    auto lowSolid=solid,lowFluid=low,highSolid=solid,highFluid=high;
    const auto a=advanceSurfaceDissolution(reaction,lowSolid,lowFluid,rates,face,locked,model.molarVolume,volume,1.,options);
    const auto b=advanceSurfaceDissolution(reaction,highSolid,highFluid,rates,face,locked,model.molarVolume,volume,1.,options);
    require(a.dissolved[0]>b.dissolved[0]&&b.dissolved[0]>0,"Surface reaction did not use each liquid owner's saturation state");
    balance(solid,low,lowSolid,lowFluid);balance(solid,high,highSolid,highFluid);
    close(a.retreat[0],a.dissolved[0]*model.molarVolume[0]/area,1e-14,"Surface recession has incorrect volume/area scaling");
    std::cout<<std::setprecision(16)<<"test,value\nlow_SI_dissolution,"<<a.dissolved[0]<<"\nhigh_SI_dissolution,"<<b.dissolved[0]<<'\n';
    for(double factor:{1e-12,1e6}){
        auto scaledSolid=scaled(solid,factor),scaledFluid=scaled(low,factor);auto scaledRates=rates;scaledRates.geometry.bulkVolume*=factor;
        auto scaledLocked=locked;for(double& n:scaledLocked)n*=factor;
        const auto result=advanceSurfaceDissolution(reaction,scaledSolid,scaledFluid,scaledRates,face,scaledLocked,model.molarVolume,volume,1.,options);
        close(result.dissolved[0]/factor,a.dissolved[0],2e-9,"Surface dissolution depends on inventory volume basis");
        close(result.retreat[0],a.retreat[0],2e-9,"Surface recession depends on inventory volume basis");
        for(int j=0;j<7;j++)close(scaledFluid.aqueous[j]/factor,lowFluid.aqueous[j],2e-10,"Liquid transfer depends on inventory volume basis");
    }
    std::cout<<"inventory_basis_invariance,1\n";
    {
        auto oversaturated=spec.inventory({{7,0,0,0,7,0,0}},.001,{});const auto initial=oversaturated;
        auto primary=solid;
        const auto result=advanceSurfaceDissolution(reaction,primary,oversaturated,rates,face,locked,model.molarVolume,volume,10.,options);
        require(result.dissolved==MineralVector{}&&primary.mineral==solid.mineral&&oversaturated.aqueous==initial.aqueous&&oversaturated.waterEquivalent==initial.waterEquivalent,
            "Supersaturated face precipitated primary mineral or dissolved against affinity");
        std::cout<<"supersaturated_primary_is_inactive,1\n";
    }
    for(int method=0;method<
#ifdef SALT_HAS_CVODE
        3
#else
        2
#endif
        ;method++){
        auto finiteSolid=solid;finiteSolid.mineral[0]=locked[0]+1e-6;const auto initialSolid=finiteSolid;
        auto first=low,second=high;options.method=method;
        const auto firstResult=advanceSurfaceDissolution(reaction,finiteSolid,first,rates,face,locked,model.molarVolume,volume,10.,options);
        close(finiteSolid.mineral[0],locked[0],2e-14,"Surface dissolution consumed a locked cap or failed to exhaust finite resolved mineral");
        balance(initialSolid,low,finiteSolid,first);
        const auto afterFirst=finiteSolid;
        const auto secondResult=advanceSurfaceDissolution(reaction,finiteSolid,second,rates,face,locked,model.molarVolume,volume,10.,options);
        require(secondResult.dissolved[0]==0&&second.aqueous==high.aqueous&&second.waterEquivalent==high.waterEquivalent,"Second fluid borrowed already exhausted mineral stock");
        require(finiteSolid.mineral==afterFirst.mineral,"Exhausted face changed the solid inventory");
        close(firstResult.dissolved[0]+secondResult.dissolved[0],initialSolid.mineral[0]-locked[0],2e-14,"Shared-face dissolution double-counted inventory");
        std::cout<<"shared_finite_stock_method_"<<method<<",1\n";
    }
    {
        auto firstSolid=solid,secondSolid=solid,fluid=low;MineralVector noCaps{};
        const auto first=advanceSurfaceDissolution(reaction,firstSolid,fluid,rates,face,noCaps,model.molarVolume,volume,1e6,options);
        const auto fluidAfterFirst=fluid;
        const auto second=advanceSurfaceDissolution(reaction,secondSolid,fluid,rates,face,noCaps,model.molarVolume,volume,1e6,options);
        require(first.dissolved[0]>1e-3&&second.dissolved[0]<first.dissolved[0]*1e-5,
            "Second solid did not see saturation established in the shared fluid");
        balance(solid,low,firstSolid,fluidAfterFirst);balance(solid,fluidAfterFirst,secondSolid,fluid);
        std::cout<<"shared_liquid_second_to_first_dissolution,"<<second.dissolved[0]/first.dissolved[0]<<'\n';
    }
    {
        Inventory hydrate=solid;hydrate.mineral.fill(0);hydrate.mineral[2]=.001001;
        auto fluid=low;const auto initialHydrate=hydrate;MineralVector hydrateFace{},hydrateLocked{};hydrateFace[2]=area;hydrateLocked[2]=.001;
        const auto result=advanceSurfaceDissolution(reaction,hydrate,fluid,rates,hydrateFace,hydrateLocked,model.molarVolume,volume,10.,options);
        balance(initialHydrate,low,hydrate,fluid);
        close(hydrate.mineral[2],hydrateLocked[2],2e-14,"Hydrated mineral did not stop at its locked inventory");
        close(fluid.waterEquivalent-low.waterEquivalent,6*waterMolarMass*result.dissolved[2],2e-12,"Hydrate water was not transferred to the liquid owner");
        std::cout<<"hydrate_water_transfer,"<<fluid.waterEquivalent-low.waterEquivalent<<'\n';
    }
    {
        auto invalidSolid=solid,invalidFluid=low;invalidFluid.waterEquivalent=0;const auto saved=invalidFluid;
        bool failed=false;try{advanceSurfaceDissolution(reaction,invalidSolid,invalidFluid,rates,face,locked,model.molarVolume,volume,1.,options);}catch(const std::invalid_argument&){failed=true;}
        require(failed&&invalidSolid.mineral==solid.mineral&&invalidFluid.aqueous==saved.aqueous&&invalidFluid.waterEquivalent==saved.waterEquivalent,"Failed surface reaction partially committed owners");
        failed=false;try{advanceSurfaceDissolution(reaction,invalidSolid,invalidSolid,rates,face,locked,model.molarVolume,volume,1.,options);}catch(const std::invalid_argument&){failed=true;}
        require(failed,"Aliased solid/liquid owners were accepted");
        auto closedFace=face;closedFace.fill(0);auto unchangedSolid=solid,unchangedFluid=low;
        const auto result=advanceSurfaceDissolution(reaction,unchangedSolid,unchangedFluid,rates,closedFace,locked,model.molarVolume,volume,10.,options);
        require(result.dissolved==MineralVector{}&&unchangedSolid.mineral==solid.mineral&&unchangedFluid.aqueous==low.aqueous&&unchangedFluid.waterEquivalent==low.waterEquivalent,"Closed face transferred mass");
        std::cout<<"positive_water_atomic_commit_and_closed_face,1\n";
    }
    {
        SurfaceFilmParameters film;film.coatedVelocity[0]=1e-6;
        auto primary=solid,fluid=low;
        const auto bare=advanceSurfaceDissolution(reaction,primary,fluid,rates,face,locked,model.molarVolume,volume,1.,options,&film);
        require(primary.mineral==lowSolid.mineral&&fluid.aqueous==lowFluid.aqueous&&fluid.waterEquivalent==lowFluid.waterEquivalent,
                "Zero-coverage film changed bare surface kinetics");
        require(bare.filmEvaluations==0&&bare.filmTransmission[0]==1,"Zero coverage evaluated an unused coated patch");
        film.coverage[0]=1;primary=solid;fluid=low;
        const auto covered=advanceSurfaceDissolution(reaction,primary,fluid,rates,face,locked,model.molarVolume,volume,1.,options,&film);
        require(covered.dissolved[0]>0&&covered.dissolved[0]<a.dissolved[0]&&covered.filmEvaluations>1,
                "Complete permeable coverage blocked dissolution or ignored film resistance");
        require(covered.filmTransmission[0]>0&&covered.filmTransmission[0]<1,"Film transmission diagnostic has the wrong range");
        balance(solid,low,primary,fluid);
        film.coatedVelocity[0]=0;primary=solid;fluid=low;
        const auto sealed=advanceSurfaceDissolution(reaction,primary,fluid,rates,face,locked,model.molarVolume,volume,1.,options,&film);
        require(sealed.dissolved==MineralVector{}&&sealed.filmEvaluations==0&&sealed.filmTransmission[0]==0&&primary.mineral==solid.mineral&&fluid.aqueous==low.aqueous,
                "Sealed coating reacted or performed a redundant film speciation");
        film.coverage[0]=.6;primary=solid;fluid=low;
        const auto partial=advanceSurfaceDissolution(reaction,primary,fluid,rates,face,locked,model.molarVolume,volume,1.,options,&film);
        auto referenceRates=rates;referenceRates.neutral[0]*=.4;auto referenceSolid=solid,referenceFluid=low;
        const auto reference=advanceSurfaceDissolution(reaction,referenceSolid,referenceFluid,referenceRates,face,locked,model.molarVolume,volume,1.,options);
        close(partial.dissolved[0],reference.dissolved[0],2e-11,"Bare/sealed parallel patches have incorrect area weighting");
        require(partial.filmEvaluations==0,"Algebraic bare/sealed mixture performed a redundant film speciation");
        std::cout<<"fully_covered_film_dissolution,"<<covered.dissolved[0]<<"\ncovered_film_transmission,"<<covered.filmTransmission[0]<<"\nfilm_fast_paths,1\n";
    }
    {
        RateParameters inactive;const auto bulk=reaction.evaluate(low,inactive,nullptr,true).aqueous;
        SurfaceFilmParameters film;film.coverage[0]=1;film.coatedVelocity[0]=1e-5;
        auto weightedRates=rates;weightedRates.geometry.weights[0]=.25;
        const auto patch=solveFilmDissolution(reaction,low,Halite,weightedRates.neutral[0]*weightedRates.geometry.weights[0],bulk.density*film.coatedVelocity[0],&bulk);
        const double affinity=-std::expm1(ln10*bulk.saturationIndex[0]);
        auto primary=solid,fluid=low;
        const auto weighted=advanceSurfaceDissolution(reaction,primary,fluid,weightedRates,face,locked,model.molarVolume,volume,1.,options,&film);
        auto equivalentRates=rates;equivalentRates.neutral[0]=patch.rate/affinity;auto equivalentSolid=solid,equivalentFluid=low;
        const auto equivalent=advanceSurfaceDissolution(reaction,equivalentSolid,equivalentFluid,equivalentRates,face,locked,model.molarVolume,volume,1.,options);
        close(weighted.dissolved[0],equivalent.dissolved[0],2e-9,"Surface weight was omitted or applied twice in film coupling");
        require(weighted.filmEvaluations==patch.evaluations+1,"Film wrapper repeated bulk speciation");
        film.coverage[0]=0;film.bareVelocity[0]=film.coatedVelocity[0];film.coatedVelocity[0]=1e-9;
        primary=solid;fluid=low;
        const auto finiteBare=advanceSurfaceDissolution(reaction,primary,fluid,weightedRates,face,locked,model.molarVolume,volume,1.,options,&film);
        close(finiteBare.dissolved[0],weighted.dissolved[0],2e-9,"Finite bare and coated films with equal conductance differ");
        require(finiteBare.filmEvaluations==patch.evaluations+1,"Zero coverage solved an unused coated film");
        film.coverage[0]=.7;film.coatedVelocity[0]=film.bareVelocity[0];primary=solid;fluid=low;
        const auto identical=advanceSurfaceDissolution(reaction,primary,fluid,weightedRates,face,locked,model.molarVolume,volume,1.,options,&film);
        require(identical.filmEvaluations==patch.evaluations+1&&identical.filmTransmission[0]==1,"Identical patch conductances caused duplicate film solves");
        std::cout<<"film_surface_weight_once,1\nfinite_bare_film_shared_speciation,1\n";
    }
    {
        Inventory hydrate=solid;hydrate.mineral.fill(0);hydrate.mineral[2]=.001;
        auto primary=hydrate,fluid=low;MineralVector hydrateFace{},noCaps{};hydrateFace[2]=area;
        SurfaceFilmParameters film;film.coverage[2]=1;film.coatedVelocity[2]=1e-5;
        const auto coated=advanceSurfaceDissolution(reaction,primary,fluid,rates,hydrateFace,noCaps,model.molarVolume,volume,1.,options,&film);
        auto bareSolid=hydrate,bareFluid=low;
        const auto bare=advanceSurfaceDissolution(reaction,bareSolid,bareFluid,rates,hydrateFace,noCaps,model.molarVolume,volume,1.,options);
        require(coated.dissolved[2]>0&&coated.dissolved[2]<bare.dissolved[2],"Hydrate film did not limit a positive dissolution flux");
        balance(hydrate,low,primary,fluid);
        close(fluid.waterEquivalent-low.waterEquivalent,6*waterMolarMass*coated.dissolved[2],3e-11,"Film hydrate water and mineral extent lost their stoichiometric pairing");
        std::cout<<"film_hydrate_water_transfer,"<<fluid.waterEquivalent-low.waterEquivalent<<'\n';
    }
    {
        SurfaceFilmParameters film;film.coverage[0]=1;film.coatedVelocity[0]=1e-5;
        for(int mode=0;mode<5;mode++){
            auto badFilm=film;auto badRates=rates;
            if(mode==0)badFilm.coverage[0]=1.01;if(mode==1)badFilm.coatedVelocity[0]=-1;
            if(mode==2)badRates.acid[0]=.01;if(mode==3)badRates.affinityExponent[0]=2;if(mode==4)badRates.co2[0]=.1;
            auto primary=solid,fluid=low;bool rejected=false;
            try{advanceSurfaceDissolution(reaction,primary,fluid,badRates,face,locked,model.molarVolume,volume,1.,options,&badFilm);}catch(const std::invalid_argument&){rejected=true;}
            require(rejected&&primary.mineral==solid.mineral&&fluid.aqueous==low.aqueous&&fluid.waterEquivalent==low.waterEquivalent,
                    "Invalid or unsupported finite film partially committed inventories");
        }
        auto nonlinear=rates;nonlinear.acid[0]=.01;nonlinear.affinityExponent[0]=2;film.coverage[0]=0;
        auto primary=solid,fluid=low,referenceSolid=solid,referenceFluid=low;
        const auto filmBare=advanceSurfaceDissolution(reaction,primary,fluid,nonlinear,face,locked,model.molarVolume,volume,1.,options,&film);
        advanceSurfaceDissolution(reaction,referenceSolid,referenceFluid,nonlinear,face,locked,model.molarVolume,volume,1.,options);
        require(primary.mineral==referenceSolid.mineral&&fluid.aqueous==referenceFluid.aqueous&&filmBare.filmEvaluations==0,"Uncoated surface lost existing nonlinear/acid kinetics");
        std::cout<<"finite_film_law_and_input_guards,1\n";
    }
    {
        auto integrate=[&](int steps){
            Inventory primary;primary.mineral[0]=1;
            auto fluid=spec.inventory({{.5,0,0,0,.5,0,0}},.01,{});
            RateParameters parameters;parameters.geometry.enabled=true;parameters.geometry.bulkVolume=.001;parameters.neutral[0]=.02;
            MineralVector exposed{},noCaps{};exposed[0]=.1;
            SurfaceFilmParameters film;film.coverage[0]=1;film.coatedVelocity[0]=1e-5;
            const auto initialSolid=primary,initialFluid=fluid;
            for(int step=0;step<steps;step++)advanceSurfaceDissolution(reaction,primary,fluid,parameters,exposed,noCaps,model.molarVolume,.001,2./steps,options,&film);
            balance(initialSolid,initialFluid,primary,fluid);return fluid.aqueous[Na];
        };
        const double reference=integrate(256),a=std::abs(integrate(4)-reference),b=std::abs(integrate(8)-reference),c=std::abs(integrate(16)-reference);
        require(a>b&&b>c&&a/b>1.5&&b/c>1.5,"Lagged two-owner film coupling failed first-order time refinement");
        std::cout<<"surface_film_error_4,"<<a<<"\nsurface_film_error_8,"<<b<<"\nsurface_film_error_16,"<<c<<'\n';
    }
    {
        for(int mineral:{Halite,Carnallite}){
            const auto initial=spec.inventory(mineral==Halite?Totals{{7,0,0,0,7,0,0}}:Totals{{1,4,0,7,19,0,0}},.001,{});
            RateParameters intrinsic;intrinsic.geometry.enabled=true;intrinsic.geometry.separateMineralAreas=true;
            intrinsic.geometry.bulkVolume=volume;intrinsic.geometry.growthArea[mineral]=area;intrinsic.neutral[mineral]=.02;
            SurfaceFilmParameters film;film.coverage[mineral]=1;film.coatedVelocity[mineral]=1e-5;
            auto coefficients=intrinsic;const auto diagnostic=applySurfaceFilmRates(reaction,initial,coefficients,film);
            require(initial.mineral[mineral]==0&&diagnostic.filmEvaluations>1&&coefficients.neutral[mineral]>0&&coefficients.neutral[mineral]<intrinsic.neutral[mineral],
                    "Film growth required fabricated parent stock or failed to attenuate the rate");
            auto coated=initial,bare=initial;reaction.advance(coated,1.,coefficients,options);reaction.advance(bare,1.,intrinsic,options);
            require(coated.mineral[mineral]>0&&coated.mineral[mineral]<bare.mineral[mineral],"Permeable film prevented incoming growth or increased its extent");
            const auto before=conserved(initial),after=conserved(coated);for(int j=0;j<8;j++)close(before[j],after[j],3e-14,"Film precipitation lost a component or water");
            if(hydration[mineral])close(initial.waterEquivalent-coated.waterEquivalent,hydration[mineral]*waterMolarMass*coated.mineral[mineral],3e-11,
                                      "Precipitation film did not consume stoichiometric hydrate water");
            film.coatedVelocity[mineral]=0;coefficients=intrinsic;const auto sealed=applySurfaceFilmRates(reaction,initial,coefficients,film);
            require(sealed.filmEvaluations==0&&sealed.filmTransmission[mineral]==0&&coefficients.geometry.growthArea[mineral]==0,"Sealed incoming growth retained a reactive area");
            std::cout<<"stock_free_film_growth_"<<mineralNames[mineral]<<','<<coated.mineral[mineral]<<'\n';
        }
        SurfaceFilmParameters film;film.coverage[0]=1;film.coatedVelocity[0]=1e-5;
        RateParameters growth;growth.geometry.enabled=true;growth.geometry.separateMineralAreas=true;growth.geometry.bulkVolume=volume;growth.geometry.growthArea[0]=area;growth.neutral[0]=.02;
        const auto closedGrowth=applySurfaceFilmRates(reaction,low,growth,film);
        require(growth.area[0]==0&&growth.geometry.growthArea[0]==0&&growth.neutral[0]==0,"Undersaturated incoming growth was converted into dissolution");
        auto probe=growth;probe.area[0]=area;probe.neutral[0]=.02;
        const auto probeResult=applySurfaceFilmRates(reaction,low,probe,film);
        require(low.mineral[0]==0&&probeResult.filmEvaluations>1&&probe.neutral[0]>0,"Stock-free surface diagnostic requires a fabricated mineral");
        auto curved=probe;curved.geometry.curvatureLogOmega[0]=.1;const auto saved=curved;bool rejected=false;
        try{applySurfaceFilmRates(reaction,low,curved,film);}catch(const std::invalid_argument&){rejected=true;}
        require(rejected&&curved.neutral==saved.neutral&&curved.area==saved.area,"Flat film accepted curved kinetics or partially changed coefficients");
        (void)closedGrowth;
        std::cout<<"film_direction_and_stock_free_probe_guards,1\n";
    }
    {
        auto integrate=[&](int steps){
            auto cell=spec.inventory({{7,0,0,0,7,0,0}},.001,{});
            RateParameters intrinsic;intrinsic.geometry.enabled=true;intrinsic.geometry.separateMineralAreas=true;
            intrinsic.geometry.bulkVolume=.001;intrinsic.geometry.growthArea[0]=.1;intrinsic.neutral[0]=.02;
            SurfaceFilmParameters film;film.coverage[0]=1;film.coatedVelocity[0]=1e-5;
            const auto initial=conserved(cell);
            for(int step=0;step<steps;step++){auto parameters=intrinsic;applySurfaceFilmRates(reaction,cell,parameters,film);reaction.advance(cell,2./steps,parameters,options);}
            const auto final=conserved(cell);for(int j=0;j<8;j++)close(initial[j],final[j],5e-14,"Refined precipitation film lost conserved inventory");
            return cell.aqueous[Na];
        };
        const double reference=integrate(256),a=std::abs(integrate(4)-reference),b=std::abs(integrate(8)-reference),c=std::abs(integrate(16)-reference);
        require(a>b&&b>c&&a/b>1.5&&b/c>1.5,"Lagged precipitation-film coupling failed first-order time refinement");
        std::cout<<"precipitation_film_error_4,"<<a<<"\nprecipitation_film_error_8,"<<b<<"\nprecipitation_film_error_16,"<<c<<'\n';
    }
    {
        SurfaceFilmParameters film;film.coverage[0]=.4;film.bareVelocity[0]=1e-5;film.coatedVelocity[0]=1e-7;
        RateParameters baseline;baseline.geometry.enabled=true;baseline.geometry.separateMineralAreas=true;baseline.area[0]=1;baseline.neutral[0]=.02;
        const auto intrinsic=baseline;applySurfaceFilmRates(reaction,low,baseline,film);
        for(double weight:{1e-300,1e-100,1.,1e100,1e300}){
            auto composite=intrinsic;composite.geometry.weights[0]=weight;composite.neutral[0]=.02/weight;
            applySurfaceFilmRates(reaction,low,composite,film);
            close(composite.neutral[0]*weight,baseline.neutral[0],2e-10,"Equivalent composite weights changed the film rate");
        }
        auto near=haliteSaturatedInlet(spec);near.aqueous[Na]*=.9999;near.aqueous[Cl]*=.9999;
        film.coverage[0]=1;film.bareVelocity[0]=INFINITY;film.coatedVelocity[0]=1e-5;
        auto tiny=intrinsic;tiny.geometry.weights[0]=1e-320;
        applySurfaceFilmRates(reaction,near,tiny,film);
        close(tiny.neutral[0],intrinsic.neutral[0],2e-8,"Near-equilibrium tiny composite rate underflow disabled a finite film multiplier");
        for(bool overflow:{false,true}){
            auto unsupported=intrinsic;
            unsupported.neutral[0]=unsupported.geometry.weights[0]=overflow?std::numeric_limits<double>::max():std::numeric_limits<double>::denorm_min();
            bool rejected=false;try{applySurfaceFilmRates(reaction,low,unsupported,film);}catch(const std::exception&){rejected=true;}
            require(rejected,"Unrepresentable composite film coefficient was silently accepted");
        }
        std::cout<<"extreme_composite_surface_weights,1\n";
    }
    std::cout<<"all_passed,1\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
