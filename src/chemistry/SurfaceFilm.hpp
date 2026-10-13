// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "BrineReactions.hpp"

namespace saltchem {
struct SurfaceFilmParameters {
    MineralVector coverage{};
    MineralVector bareVelocity{{INFINITY,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY}};
    MineralVector coatedVelocity{{INFINITY,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY}};
};
struct SurfaceFilmRateResult {
    int filmEvaluations=0;
    MineralVector filmTransmission{{1,1,1,1,1,1,1}};
};
struct FilmDissolution {
    double rate=0,kineticFraction=0,saturationIndex=0,waterActivity=1,residual=0,peclet=0;
    int evaluations=0;
    // A sealed film has no uniquely determined steady surface composition.
    // In that limit aqueous/fractions retain the bulk state as placeholders.
    bool surfaceStateAvailable=true;
    std::array<double,8> surfaceMassFractions{};
    AqueousState aqueous;
};

// One dissolution-only, neutral, linear-affinity mineral patch. k includes
// any surface weight, in mol/(m2 s); G=rho*D_eff/delta is in kg/(m2 s).
// The eight mass coordinates are the seven hydrogen-balanced components and
// waterEquivalent. For the outward total flux M*r, exact integration of a
// constant-conductance Stefan film gives Ys=Yb*exp(-Pe)+Ymin*(1-exp(-Pe)).
// Pitzer evaluates the surface affinity; no fixed affinity/concentration
// conversion or independently transported hydrate-water approximation occurs.
// This algebraic closure freezes geometry/conductance, not the bulk kinetics.
inline FilmDissolution solveFilmDissolution(
    const Reactions& reactions,const Inventory& bulk,int mineral,
    double neutralRate,double massConductance,const AqueousState* bulkAqueous=nullptr){
    if(mineral<0||mineral>=nMinerals||!std::isfinite(neutralRate)||neutralRate<0||
       std::isnan(massConductance)||massConductance<0||
       !std::isfinite(bulk.waterEquivalent)||bulk.waterEquivalent<=0)
        throw std::invalid_argument("Invalid Pitzer film mineral, rate, conductance or solvent");
    for(double n:bulk.aqueous)if(!std::isfinite(n)||n<0)
        throw std::invalid_argument("Invalid Pitzer film bulk composition");
    long double extendedMass=bulk.waterEquivalent;
    for(int i=0;i<nComponents;i++)extendedMass+=static_cast<long double>(componentMass[i])*bulk.aqueous[i];
    if(extendedMass>std::numeric_limits<double>::max())throw std::invalid_argument("Invalid Pitzer film bulk mass");
    const double mass=double(extendedMass);
    FilmDissolution result;
    std::array<double,8> bulkY{},mineralY{};
    for(int i=0;i<nComponents;i++){
        bulkY[i]=componentMass[i]*bulk.aqueous[i]/mass;
        mineralY[i]=componentMass[i]*stoichiometry[mineral][i]/mineralMass[mineral];
    }
    bulkY[7]=bulk.waterEquivalent/mass;
    mineralY[7]=hydration[mineral]*waterMolarMass/mineralMass[mineral];
    RateParameters inactive;
    if(bulkAqueous)result.aqueous=*bulkAqueous;
    else{result.aqueous=reactions.evaluate(bulk,inactive,nullptr,false).aqueous;result.evaluations++;}
    result.surfaceMassFractions=bulkY;
    result.saturationIndex=result.aqueous.saturationIndex[mineral];
    result.waterActivity=result.aqueous.waterActivity;
    if(!std::isfinite(result.saturationIndex)||!std::isfinite(result.waterActivity)||result.waterActivity<=0)
        throw std::invalid_argument("Invalid supplied Pitzer film bulk speciation");
    if(massConductance==0){result.surfaceStateAvailable=false;return result;}
    if(neutralRate==0||result.saturationIndex>=0)return result;
    const double bulkAffinity=-std::expm1(ln10*result.saturationIndex);
    if(std::isinf(massConductance)){result.kineticFraction=bulkAffinity;result.rate=neutralRate*bulkAffinity;return result;}
    // Pe, instead of r, resolves the near-equilibrium surface even as G->0.
    // Extended precision also keeps M*k/G finite over the double input range.
    const long double beta=static_cast<long double>(mineralMass[mineral])*neutralRate/massConductance;
    if(beta==0){result.kineticFraction=bulkAffinity;result.rate=neutralRate*bulkAffinity;return result;}
    struct Point {long double pe=0;double fraction=0,defect=0;std::array<double,8> y{};AqueousState aqueous;};
    Point low;low.y=bulkY;low.aqueous=result.aqueous;low.defect=-bulkAffinity;
    auto evaluate=[&](long double pe,const AqueousState& seed){
        Point point;point.pe=pe;point.fraction=double(pe/beta);
        const double argument=double(std::min(pe,static_cast<long double>(745)));
        const double bulkWeight=pe>=745?0:std::exp(-argument),mineralWeight=-std::expm1(-argument);
        Inventory surface;surface.mineral.fill(0);
        for(int i=0;i<8;i++)point.y[i]=bulkWeight*bulkY[i]+mineralWeight*mineralY[i];
        surface.waterEquivalent=point.y[7];
        if(!(surface.waterEquivalent>=std::numeric_limits<double>::min()))
            throw std::domain_error("Pitzer film surface has no numerically resolvable solvent");
        for(int i=0;i<nComponents;i++)surface.aqueous[i]=point.y[i]/componentMass[i];
        result.evaluations++;
        point.aqueous=reactions.evaluate(surface,inactive,&seed,false).aqueous;
        const double logOmega=ln10*point.aqueous.saturationIndex[mineral];
        if(!std::isfinite(logOmega))throw std::domain_error("Nonfinite Pitzer film surface affinity");
        // Only its positive sign is needed far beyond saturation; avoid
        // overflowing while constructing a bracket that will be reduced.
        point.defect=point.fraction+std::expm1(std::min(logOmega,50.));
        return point;
    };
    Point high;bool bracket=false;long double trial=std::min(beta,1.L);
    std::string lastFailure;
    for(int attempt=0;attempt<96&&!bracket;attempt++){
        try{
            auto point=evaluate(trial,low.aqueous);
            if(point.defect>=0){high=point;bracket=true;break;}
            low=point;
            if(trial==beta)throw std::runtime_error("Pitzer film failed to bracket its bounded kinetic rate");
            // At Pe>=64 the composition has reached its mineral endmember.
            // Jump to the kinetic bound only if that endmember was valid.
            trial=trial>=64?beta:std::min(beta,trial*2);
        }catch(const std::exception& error){
            lastFailure=error.what();long double invalid=trial;
            for(int backtrack=0;backtrack<96;backtrack++){
                const long double middle=(low.pe+invalid)/2;
                if(middle==low.pe||middle==invalid)break;
                try{
                    auto point=evaluate(middle,low.aqueous);
                    if(point.defect>=0){high=point;bracket=true;break;}
                    low=point;
                }catch(const std::exception& failure){lastFailure=failure.what();invalid=middle;}
            }
            if(!bracket)throw std::runtime_error("No valid Pitzer film equilibrium bracket: "+lastFailure);
        }
    }
    if(!bracket)throw std::runtime_error("Pitzer film bracket did not converge");
    const double tolerance=2e-12+2e-10*bulkAffinity;
    Point a=low,b=high,c=low;
    if(std::abs(a.defect)<std::abs(b.defect))std::swap(a,b);
    c=a;long double previous=c.pe;bool bisected=true;
    Point solution=b;bool converged=std::abs(b.defect)<=tolerance;
    for(int iteration=0;iteration<96;iteration++){
        if(converged)break;
        // Brent's safeguarded inverse interpolation reduces expensive Pitzer
        // calls while keeping a valid opposite-sign bracket at every step.
        long double candidate;
        if(a.defect!=c.defect&&b.defect!=c.defect){
            const long double fa=a.defect,fb=b.defect,fc=c.defect;
            candidate=a.pe*fb*fc/((fa-fb)*(fa-fc))+b.pe*fa*fc/((fb-fa)*(fb-fc))+c.pe*fa*fb/((fc-fa)*(fc-fb));
        }else candidate=b.pe-b.defect*(b.pe-a.pe)/(b.defect-a.defect);
        const long double inside=(3*a.pe+b.pe)/4;
        const long double precision=8*std::numeric_limits<long double>::epsilon()*std::max(std::abs(a.pe),std::abs(b.pe));
        if(!std::isfinite(candidate)||candidate<=std::min(inside,b.pe)||candidate>=std::max(inside,b.pe)||
           (bisected&&std::abs(candidate-b.pe)>=std::abs(b.pe-c.pe)/2)||
           (!bisected&&std::abs(candidate-b.pe)>=std::abs(c.pe-previous)/2)||
           (bisected&&std::abs(b.pe-c.pe)<=precision)||(!bisected&&std::abs(c.pe-previous)<=precision)){
            candidate=(a.pe+b.pe)/2;bisected=true;
        }else bisected=false;
        solution=evaluate(candidate,b.aqueous);
        if(std::abs(solution.defect)<=tolerance){converged=true;break;}
        previous=c.pe;c=b;
        if((a.defect<0)!=(solution.defect<0))b=solution;else a=solution;
        if(std::abs(a.defect)<std::abs(b.defect))std::swap(a,b);
    }
    if(!converged)throw std::runtime_error("Pitzer film nonlinear surface residual did not converge");
    result.kineticFraction=solution.fraction;result.rate=neutralRate*solution.fraction;
    result.peclet=double(std::min(solution.pe,static_cast<long double>(std::numeric_limits<double>::max())));
    result.surfaceMassFractions=solution.y;result.aqueous=solution.aqueous;
    result.saturationIndex=solution.aqueous.saturationIndex[mineral];result.waterActivity=solution.aqueous.waterActivity;
    result.residual=neutralRate*solution.defect;
    return result;
}

// The same conservative Stefan film on a precipitating flat patch. q=-Pe
// is positive and Ys=Yb+(Yb-Ymin)*expm1(q). Each consumed mass coordinate,
// including hydrate water, supplies an upper q bound. No negative surface
// coordinate is clipped into a purported thermodynamic solution.
inline FilmDissolution solveFilmPrecipitation(
    const Reactions& reactions,const Inventory& bulk,int mineral,double neutralRate,
    double massConductance,const AqueousState* bulkAqueous=nullptr){
    if(mineral<0||mineral>=nMinerals||!std::isfinite(neutralRate)||neutralRate<0||
       std::isnan(massConductance)||massConductance<0||!std::isfinite(bulk.waterEquivalent)||bulk.waterEquivalent<=0)
        throw std::invalid_argument("Invalid Pitzer precipitation-film input");
    long double mass=bulk.waterEquivalent;
    for(int i=0;i<nComponents;i++){
        if(!std::isfinite(bulk.aqueous[i])||bulk.aqueous[i]<0)throw std::invalid_argument("Invalid precipitation-film bulk composition");
        mass+=static_cast<long double>(componentMass[i])*bulk.aqueous[i];
    }
    if(mass>std::numeric_limits<double>::max())throw std::invalid_argument("Invalid precipitation-film bulk mass");
    FilmDissolution result;std::array<double,8> yb{},ym{};
    for(int i=0;i<nComponents;i++){yb[i]=componentMass[i]*bulk.aqueous[i]/double(mass);ym[i]=componentMass[i]*stoichiometry[mineral][i]/mineralMass[mineral];}
    yb[7]=bulk.waterEquivalent/double(mass);ym[7]=hydration[mineral]*waterMolarMass/mineralMass[mineral];result.surfaceMassFractions=yb;
    RateParameters inactive;
    if(bulkAqueous)result.aqueous=*bulkAqueous;
    else{result.aqueous=reactions.evaluate(bulk,inactive,nullptr,false).aqueous;result.evaluations++;}
    result.saturationIndex=result.aqueous.saturationIndex[mineral];result.waterActivity=result.aqueous.waterActivity;
    if(!std::isfinite(result.saturationIndex)||!std::isfinite(result.waterActivity)||result.waterActivity<=0)
        throw std::invalid_argument("Invalid precipitation-film bulk speciation");
    if(massConductance==0){result.surfaceStateAvailable=false;return result;}
    if(neutralRate==0||result.saturationIndex<=0)return result;
    if(result.saturationIndex>700/ln10)throw std::domain_error("Precipitation-film bulk saturation exceeds the supported kinetic range");
    const double logOmega=ln10*result.saturationIndex;
    const double supersaturation=std::expm1(logOmega);
    if(static_cast<long double>(neutralRate)*supersaturation>std::numeric_limits<double>::max())throw std::domain_error("Precipitation-film bare rate overflows");
    const double bareRate=-neutralRate*supersaturation;
    if(std::isinf(massConductance)){result.kineticFraction=-supersaturation;result.rate=bareRate;return result;}
    const long double beta=static_cast<long double>(mineralMass[mineral])*neutralRate/massConductance;
    long double upper=std::numeric_limits<long double>::infinity();bool same=true;
    for(int i=0;i<8;i++){
        same=same&&yb[i]==ym[i];
        if(ym[i]>yb[i])upper=std::min(upper,-std::log1p(-static_cast<long double>(yb[i])/ym[i]));
    }
    if(!(upper>0))throw std::domain_error("Supersaturated film lacks a required precipitating component");
    if(same||beta==0){result.kineticFraction=-supersaturation;result.rate=bareRate;return result;}
    struct Point {long double q=0;double fraction=0,defect=0;std::array<double,8> y{};AqueousState aqueous;};
    Point low;low.y=yb;low.aqueous=result.aqueous;low.defect=-supersaturation;
    auto evaluate=[&](long double q,const AqueousState& seed){
        if(!(q>=0&&q<upper)||q>700)throw std::domain_error("Precipitation film exceeds its nonnegative composition interval");
        const long double fraction=q/beta;
        if(fraction>std::numeric_limits<double>::max())throw std::domain_error("Precipitation-film kinetic ratio exceeds numerical range");
        Point p;p.q=q;p.fraction=double(fraction);const long double growth=std::expm1(q);
        Inventory surface;surface.mineral.fill(0);
        for(int i=0;i<8;i++){
            const long double value=static_cast<long double>(yb[i])+(static_cast<long double>(yb[i])-ym[i])*growth;
            if(!std::isfinite(value)||value<0||value>std::numeric_limits<double>::max())throw std::domain_error("Invalid precipitation-film surface mass fraction");
            p.y[i]=double(value);
        }
        surface.waterEquivalent=p.y[7];
        if(surface.waterEquivalent<std::numeric_limits<double>::min())throw std::domain_error("Precipitation film has no resolvable surface solvent");
        for(int i=0;i<nComponents;i++)surface.aqueous[i]=p.y[i]/componentMass[i];
        result.evaluations++;p.aqueous=reactions.evaluate(surface,inactive,&seed,false).aqueous;
        const double value=ln10*p.aqueous.saturationIndex[mineral];
        if(!std::isfinite(value)||value>700)throw std::domain_error("Invalid precipitation-film surface saturation");
        p.defect=p.fraction-std::expm1(value);return p;
    };
    Point high;bool bracket=false;long double trial=std::min(beta*supersaturation,1.L);
    if(std::isfinite(upper))trial=std::min(trial,upper/2);
    std::string failure;
    for(int attempt=0;attempt<128&&!bracket;attempt++){
        try{
            auto point=evaluate(trial,low.aqueous);
            if(point.defect>=0){high=point;bracket=true;break;}
            low=point;
            const long double next=std::isfinite(upper)&&2*trial>=upper?(trial+upper)/2:trial*2;
            if(next<=trial||next>=upper)throw std::runtime_error("Precipitation film reached its admissible composition boundary without a root");
            trial=next;
        }catch(const std::exception& error){
            failure=error.what();long double invalid=trial;
            for(int backtrack=0;backtrack<96;backtrack++){
                const long double middle=(low.q+invalid)/2;if(middle==low.q||middle==invalid)break;
                try{auto point=evaluate(middle,low.aqueous);if(point.defect>=0){high=point;bracket=true;break;}low=point;}
                catch(const std::exception& error){failure=error.what();invalid=middle;}
            }
            if(!bracket)throw std::runtime_error("No admissible Pitzer precipitation-film bracket: "+failure);
        }
    }
    if(!bracket)throw std::runtime_error("Precipitation-film bracket did not converge");
    // Scale by the local rate, not the possibly huge bulk supersaturation:
    // a nearly sealed film must still resolve its near-equilibrium surface.
    const auto tolerance=[](const Point& point){return 2e-12+2e-10*std::abs(point.fraction);};
    Point a=low,b=high,c=low;if(std::abs(a.defect)<std::abs(b.defect))std::swap(a,b);
    c=a;long double previous=c.q;bool bisected=true;Point solution=b;bool converged=std::abs(b.defect)<=tolerance(b);
    for(int iteration=0;iteration<128&&!converged;iteration++){
        long double candidate;
        if(a.defect!=c.defect&&b.defect!=c.defect){const long double fa=a.defect,fb=b.defect,fc=c.defect;candidate=a.q*fb*fc/((fa-fb)*(fa-fc))+b.q*fa*fc/((fb-fa)*(fb-fc))+c.q*fa*fb/((fc-fa)*(fc-fb));}
        else candidate=b.q-b.defect*(b.q-a.q)/(b.defect-a.defect);
        const long double inside=(3*a.q+b.q)/4,precision=8*std::numeric_limits<long double>::epsilon()*std::max(std::abs(a.q),std::abs(b.q));
        if(!std::isfinite(candidate)||candidate<=std::min(inside,b.q)||candidate>=std::max(inside,b.q)||
           (bisected&&std::abs(candidate-b.q)>=std::abs(b.q-c.q)/2)||(!bisected&&std::abs(candidate-b.q)>=std::abs(c.q-previous)/2)||
           (bisected&&std::abs(b.q-c.q)<=precision)||(!bisected&&std::abs(c.q-previous)<=precision)){
            candidate=(a.q+b.q)/2;bisected=true;
        }else bisected=false;
        solution=evaluate(candidate,b.aqueous);if(std::abs(solution.defect)<=tolerance(solution)){converged=true;break;}
        previous=c.q;c=b;if((a.defect<0)!=(solution.defect<0))b=solution;else a=solution;
        if(std::abs(a.defect)<std::abs(b.defect))std::swap(a,b);
    }
    if(!converged)throw std::runtime_error("Pitzer precipitation-film nonlinear residual did not converge");
    result.kineticFraction=-solution.fraction;result.rate=-neutralRate*solution.fraction;result.peclet=-double(solution.q);result.surfaceMassFractions=solution.y;result.aqueous=solution.aqueous;
    result.saturationIndex=solution.aqueous.saturationIndex[mineral];result.waterActivity=solution.aqueous.waterActivity;result.residual=-neutralRate*solution.defect;
    return result;
}

// Prepare a lagged signed film multiplier for the existing coupled ODE.
// Reactive area, not mineral stock, determines whether a hypothetical or
// incoming growth patch exists. Callers suppress unavailable solid stock.
// No-film paths preserve acid/nonlinear kinetics without another speciation.
inline SurfaceFilmRateResult applySurfaceFilmRates(const Reactions& reactions,const Inventory& bulkInventory,
                                                    RateParameters& parameters,const SurfaceFilmParameters& film){
    Reactions::validateParameters(parameters);RateParameters p=parameters;SurfaceFilmRateResult result;
    if(!p.geometry.enabled||!p.geometry.separateMineralAreas)throw std::invalid_argument("Surface film rates require explicit mineral areas");
    AqueousState bulk;bool haveBulk=false;
    for(int j=0;j<nMinerals;j++){
        const double theta=film.coverage[j],bareVelocity=film.bareVelocity[j],coatedVelocity=film.coatedVelocity[j];
        if(!std::isfinite(theta)||theta<0||theta>1||std::isnan(bareVelocity)||bareVelocity<0||std::isnan(coatedVelocity)||coatedVelocity<0)
            throw std::invalid_argument("Invalid surface film coverage or diffusion velocity");
        const bool dissolution=p.area[j]>0,precipitation=p.geometry.growthArea[j]>0;
        if((!dissolution&&!precipitation)||p.geometry.weights[j]==0)continue;
        const bool finiteBare=theta<1&&std::isfinite(bareVelocity),finiteCoated=theta>0&&std::isfinite(coatedVelocity);
        if(!finiteBare&&!finiteCoated)continue;
        if(p.acid[j]!=0||p.co2[j]!=0||p.affinityExponent[j]!=1||p.geometry.curvatureLogOmega[j]!=0)
            throw std::invalid_argument("Finite surface films require neutral, linear-affinity kinetics on a flat mineral patch");
        if(p.neutral[j]==0)continue;
        const bool needBare=theta<1&&bareVelocity>0&&std::isfinite(bareVelocity),needCoated=theta>0&&coatedVelocity>0&&std::isfinite(coatedVelocity);
        if(!needBare&&!needCoated){
            const double factor=(1-theta)*(bareVelocity>0)+theta*(coatedVelocity>0);p.neutral[j]*=factor;
            if(factor==0){p.area[j]=0;p.geometry.growthArea[j]=0;}
            if(theta>0)result.filmTransmission[j]=coatedVelocity==0?0:bareVelocity==0?INFINITY:1;
            continue;
        }
        if(!haveBulk){RateParameters inactive;bulk=reactions.evaluate(bulkInventory,inactive,nullptr,true).aqueous;
            if(!std::isfinite(bulk.density)||bulk.density<=0)throw std::runtime_error("Invalid bulk density for surface film");haveBulk=true;result.filmEvaluations++;}
        const double si=bulk.saturationIndex[j];
        if(!std::isfinite(si)||si>700/ln10)throw std::runtime_error("Invalid bulk saturation index for surface film");
        if(si==0||(si<0&&!dissolution)||(si>0&&!precipitation)){
            p.neutral[j]=0;p.area[j]=0;p.geometry.growthArea[j]=0;continue;
        }
        const double affinity=-std::expm1(ln10*si),weight=p.geometry.weights[j];
        if(static_cast<long double>(weight)*p.neutral[j]>std::numeric_limits<double>::max())throw std::invalid_argument("Surface-weighted neutral coefficient overflows");
        const double k=p.neutral[j]*weight;
        if(static_cast<long double>(k)*std::abs(affinity)>std::numeric_limits<double>::max())throw std::domain_error("Surface-film bulk kinetic rate overflows");
        if(k==0)throw std::domain_error("Surface-weighted neutral coefficient underflows the supported film range");
        auto patchFraction=[&](double velocity){
            if(velocity==0)return 0.;if(std::isinf(velocity))return affinity;
            const long double extendedConductance=static_cast<long double>(bulk.density)*velocity;
            const double conductance=extendedConductance>std::numeric_limits<double>::max()?INFINITY:double(extendedConductance);
            const auto solution=si<0?solveFilmDissolution(reactions,bulkInventory,j,k,conductance,&bulk):solveFilmPrecipitation(reactions,bulkInventory,j,k,conductance,&bulk);
            result.filmEvaluations+=solution.evaluations;return solution.kineticFraction;
        };
        // Use r/k instead of first multiplying by a very small k and then
        // dividing it out. The latter can silently turn a finite multiplier
        // into zero for small composite surface weights near saturation.
        const double bare=patchFraction(bareVelocity),covered=theta==0?bare:coatedVelocity==bareVelocity?bare:patchFraction(coatedVelocity);
        const long double ratio=((1-static_cast<long double>(theta))*bare+static_cast<long double>(theta)*covered)/affinity;
        if(!std::isfinite(ratio)||ratio<0)throw std::runtime_error("Invalid signed surface film rate multiplier");
        const long double effective=static_cast<long double>(p.neutral[j])*ratio;
        if(effective>std::numeric_limits<double>::max())throw std::domain_error("Surface-film effective coefficient overflows");
        if(effective>0&&double(effective)==0)throw std::domain_error("Surface-film effective coefficient underflows the supported kinetic representation");
        p.neutral[j]=double(effective);if(ratio==0){p.area[j]=0;p.geometry.growthArea[j]=0;}
        if(theta>0)result.filmTransmission[j]=bare!=0?covered/bare:covered!=0?INFINITY:1;
    }
    parameters=p;return result;
}
}
