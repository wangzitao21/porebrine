// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "BrineSpeciation.hpp"
#include <sstream>
#include <iomanip>
#include <boost/numeric/odeint/stepper/rosenbrock4.hpp>
#include "PhaseVector.hpp"
#ifdef SALT_HAS_CVODE
#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#endif
namespace saltchem {
struct PhaseGeometry {
    bool enabled=false;
    bool separateMineralAreas=false;
    MineralVector growthArea{},surfaceExponent{},lockedMoles{},curvatureLogOmega{};
    double bulkVolume=1,inertFraction=0,dx=5e-6,dy=5e-6,epsilonMin=1e-5,directionWidth=1e-9,closureWidth=1e-8;
    double initialPorosity=std::numeric_limits<double>::quiet_NaN();
    std::array<double,4> neighbours{{1,1,1,1}};
    MineralVector weights{{1,1,1,1,1,1,1}};
};
struct RateParameters {
    // k has units mol/(m2 s); H/CO2 terms may be used for acid-catalyzed TST.
    MineralVector neutral{},acid{},co2{},acidOrder{{1,1,1,1,1,1,1}};
    MineralVector area{},referenceMoles{{1,1,1,1,1,1,1}};
    MineralVector dissolutionExponent{{1,1,1,1,1,1,1}};
    MineralVector affinityExponent{{1,1,1,1,1,1,1}};
    PhaseGeometry geometry;
    // Externally specified reactive area is used for precipitation, including
    // a mineral initially absent on an available substrate. area=0 disables it.
};
inline double poreFraction(const Inventory& s,const PhaseGeometry& g,const MineralVector& volume){
    long double e=1-static_cast<long double>(g.inertFraction);for(int j=0;j<nMinerals;j++)e-=static_cast<long double>(s.mineral[j])*volume[j]/g.bulkVolume;return double(e);
}
// The same finite-rate surface law is used with native or PHREEQC activities.
// Geometry supplies area from the paper's phase-vector operator. Neighbouring
// pore fractions are held at the transport time level during each local ODE.
template<class Real> inline std::array<Real,nMinerals> surfaceRatesAt(const std::array<Real,nMinerals>& mineral,const RateParameters& p,const std::array<Real,nMinerals>& si,Real logH,Real logCO2,const MineralVector& volume,Real ep,Real* phaseQ=nullptr,Real* interfaceArea=nullptr,const std::array<Real,nMinerals>* saturationRatio=nullptr,const Real* co2Activity=nullptr){
    const auto& g=p.geometry;
    if(g.enabled&&(scalarValue(ep)<g.epsilonMin-1e-10||scalarValue(ep)>1+1e-12))throw std::domain_error("Reaction stage exceeds pore-volume bounds");
    std::array<Real,nMinerals> rates{};Real Q=0;
    for(int j=0;j<nMinerals;j++){
        if((!g.enabled&&p.area[j]==0)||(g.enabled&&g.weights[j]==0))continue;
        const Real logOmega=ln10*si[j];
        if(!std::isfinite(scalarValue(logOmega))||scalarValue(logOmega)>700)throw std::domain_error("Saturation ratio exceeds numerical range");
        const Real accessible=mineral[j]-g.lockedMoles[j];
        const Real fraction=scalarValue(accessible)>0?accessible/p.referenceMoles[j]:Real(0);
        const bool curved=g.curvatureLogOmega[j]>0;
        if(curved&&scalarValue(accessible)<=0)continue;
        Real curvature=0;
        if(curved){
            const Real radiusScale=pow(fraction,1./3);
            // Beyond this shift the largest supported bulk Omega multiplied
            // by exp(-curvature) is below exp(-700). Avoid a singular Kelvin
            // derivative as an exhausted cap reaches zero curvature radius.
            curvature=scalarValue(radiusScale)<=g.curvatureLogOmega[j]/1400?Real(1400):g.curvatureLogOmega[j]/radiusScale;
        }
        const Real effectiveLogOmega=logOmega-curvature;
        const Real affinity=saturationRatio?
            1-(curved?(*saturationRatio)[j]*exp(-curvature):(*saturationRatio)[j]):
            -expm1(scalarValue(effectiveLogOmega)>-700?effectiveLogOmega:Real(-700));
        const Real k=p.neutral[j]+p.acid[j]*pow10(p.acidOrder[j]*logH)+p.co2[j]*(co2Activity?*co2Activity:pow10(logCO2));
        Real area=g.separateMineralAreas?(scalarValue(affinity)>0?p.area[j]:g.growthArea[j]):(g.enabled?g.weights[j]:p.area[j]);
        if(g.separateMineralAreas){
            if(scalarValue(affinity)>0&&scalarValue(accessible)<=0)continue;
            const double exponent=curved?2./3:g.surfaceExponent[j];
            if(exponent>0){
                area*=pow(fraction,exponent);
            }
            if(scalarValue(affinity)<0){const Real gap=scalarValue(ep)>g.epsilonMin?ep-g.epsilonMin:Real(0);area*=gap/(gap+g.closureWidth);}
            area*=g.weights[j];
        }
        else if(scalarValue(affinity)>0){
            const Real fraction=g.enabled?mineral[j]*volume[j]/(g.bulkVolume*(scalarValue(ep)<1?1-ep+1e-6:Real(1e-6))):mineral[j]/p.referenceMoles[j];
            area*=pow(scalarValue(fraction)>=0?fraction:Real(0),p.dissolutionExponent[j]);
        }
        const Real signedAffinity=p.affinityExponent[j]==1?affinity:(scalarValue(affinity)<0?-pow(-affinity,p.affinityExponent[j]):pow(affinity,p.affinityExponent[j]));
        rates[j]=k*area*signedAffinity;
        if(!std::isfinite(scalarValue(rates[j])))throw std::domain_error("TST rate exceeds numerical range");
        Q+=rates[j]*volume[j];
    }
    if(g.enabled&&!g.separateMineralAreas){
        Real area=phaseVector(ep,g.neighbours,g.dx,g.dy,Q)[0];
        const Real positive=.5*(1+tanh(Q/g.directionWidth));
        if(std::abs(scalarValue(Q))<20*g.directionWidth)area=positive*phaseVector(ep,g.neighbours,g.dx,g.dy,Real(1))[0]+(1-positive)*phaseVector(ep,g.neighbours,g.dx,g.dy,Real(-1))[0];
        const Real gap=scalarValue(ep)>g.epsilonMin?ep-g.epsilonMin:Real(0);
        // Approach pore closure continuously. A hard stop at a positive gap
        // creates an ODE discontinuity and repeated step rejection there.
        if(scalarValue(Q)<0)area*=gap/(gap+g.closureWidth);
        for(Real& r:rates)r*=area*g.bulkVolume;
        if(interfaceArea)*interfaceArea=area;
    }
    if(g.separateMineralAreas&&interfaceArea){Real area=0;for(double a:p.area)area+=a;*interfaceArea=area/g.bulkVolume;}
    if(phaseQ)*phaseQ=Q;return rates;
}
inline MineralVector surfaceRates(const Inventory& state,const RateParameters& p,const MineralVector& si,double logH,double logCO2,const MineralVector& volume,double* phaseQ=nullptr,double* interfaceArea=nullptr,double porosity=std::numeric_limits<double>::quiet_NaN()){
    const double ep=p.geometry.enabled?(std::isfinite(porosity)?porosity:poreFraction(state,p.geometry,volume)):1;
    return surfaceRatesAt(state.mineral,p,si,logH,logCO2,volume,ep,phaseQ,interfaceArea);
}
// Absolute tolerance is mol per initial kg water-equivalent, independent of cell size.
struct IntegrationOptions {
    double relativeTolerance=2e-12,absoluteTolerance=1e-15;int maximumAttempts=100000;
    bool analyticJacobian=true;
    // 0: TR-BDF2; 1: Rosenbrock4; 2: variable-order CVODE BDF.
#ifdef SALT_HAS_CVODE
    int method=2;
#else
    int method=1;
#endif
};
struct IntegrationStats {int accepted=0,rejected=0,rateEvaluations=0,speciationIterations=0,newtonIterations=0,depletionEvents=0,historyRestarts=0,integratorFallbacks=0,analyticJacobians=0,jacobianFallbacks=0;MineralVector dissolved{},precipitated{};};
struct ReactionOutput {AqueousState aqueous;double freeWaterKg=1;MineralVector rates{};};
class Reactions {
    const PitzerModel& model;Speciation spec;
    static bool linear(double matrix[nMinerals][nMinerals],MineralVector b,MineralVector& x){
        for(int j=0;j<nMinerals;j++){int p=j;for(int i=j+1;i<nMinerals;i++)if(std::abs(matrix[i][j])>std::abs(matrix[p][j]))p=i;
            if(!std::isfinite(matrix[p][j])||std::abs(matrix[p][j])<1e-25)return false;
            if(p!=j){for(int k=j;k<nMinerals;k++)std::swap(matrix[j][k],matrix[p][k]);std::swap(b[j],b[p]);}
            for(int i=j+1;i<nMinerals;i++){double factor=matrix[i][j]/matrix[j][j];for(int k=j+1;k<nMinerals;k++)matrix[i][k]-=factor*matrix[j][k];b[i]-=factor*b[j];}}
        for(int i=nMinerals-1;i>=0;i--){x[i]=b[i];for(int j=i+1;j<nMinerals;j++)x[i]-=matrix[i][j]*x[j];x[i]/=matrix[i][i];}return true;
    }
    static bool valid(const Inventory& s){
        if(s.waterEquivalent<=0||!std::isfinite(s.waterEquivalent))return false;
        for(double x:s.aqueous)if(x<0||!std::isfinite(x))return false;
        for(double x:s.mineral)if(x<0||!std::isfinite(x))return false;return true;
    }
    bool admissible(const Inventory& s,const RateParameters& p)const{if(!valid(s))return false;if(!p.geometry.enabled)return true;const double ep=poreFraction(s,p.geometry,model.molarVolume);return ep>=p.geometry.epsilonMin-1e-10&&ep<=1+1e-12;}
    double incrementedPorosity(const Inventory& base,const MineralVector& extent,const RateParameters& p)const{
        if(!p.geometry.enabled)return std::numeric_limits<double>::quiet_NaN();
        long double ep=std::isfinite(p.geometry.initialPorosity)?p.geometry.initialPorosity:poreFraction(base,p.geometry,model.molarVolume);for(int j=0;j<nMinerals;j++)ep+=static_cast<long double>(extent[j])*model.molarVolume[j]/p.geometry.bulkVolume;return double(ep);
    }
    bool implicit(const Inventory& base,const RateParameters& p,const MineralVector& rhs,double coefficient,MineralVector& x,const AqueousState* seed,IntegrationStats& stats,const IntegrationOptions& options)const{
        x.fill(0);AqueousState last;if(seed)last=*seed;bool haveSeed=seed!=nullptr;
        double derivative[nMinerals][nMinerals]={};bool haveJacobian=false,haveEvaluation=false;ReactionOutput cached;
        auto norm=[&](const MineralVector& r,const MineralVector& y){double v=0;for(int j=0;j<nMinerals;j++){double scale=.03*(options.absoluteTolerance+options.relativeTolerance*std::abs(y[j]));v=std::max(v,std::abs(r[j])/scale);}return v;};
        for(int it=0;it<24;it++){
            stats.newtonIterations++;Inventory state=changed(base,x);if(!valid(state))return false;
            ReactionOutput eval;
            if(haveEvaluation)eval=cached;
            else {try{eval=evaluate(state,p,haveSeed?&last:nullptr,false,incrementedPorosity(base,x,p));}catch(const std::exception&){return false;}stats.rateEvaluations++;stats.speciationIterations+=eval.aqueous.iterations;}
            last=eval.aqueous;haveSeed=true;
            MineralVector residual{};for(int j=0;j<nMinerals;j++)residual[j]=x[j]-rhs[j]-coefficient*eval.rates[j];double oldNorm=norm(residual,x);if(oldNorm<=1)return true;
            double jac[nMinerals][nMinerals]={};
            if(!haveJacobian)for(int j=0;j<nMinerals;j++){
                double delta=1e-6*std::max(std::abs(x[j]),1e-3*base.waterEquivalent);MineralVector trial=x;bool ok=false;
                for(int k=0;k<16;k++){trial=x;trial[j]+=delta;if(valid(changed(base,trial))){ok=true;break;}trial=x;trial[j]-=delta;if(valid(changed(base,trial))){delta=-delta;ok=true;break;}delta*=.1;}
                for(int i=0;i<nMinerals;i++)derivative[i][j]=0;
                if(ok){try{auto e=evaluate(changed(base,trial),p,&last,false,incrementedPorosity(base,trial,p));stats.rateEvaluations++;stats.speciationIterations+=e.aqueous.iterations;for(int i=0;i<nMinerals;i++)derivative[i][j]=(e.rates[i]-eval.rates[i])/delta;}catch(const std::exception&){return false;}}
            }
            haveJacobian=true;
            for(int i=0;i<nMinerals;i++)for(int j=0;j<nMinerals;j++)jac[i][j]=(i==j?1.:0.)-coefficient*derivative[i][j];
            MineralVector step{},rhsNewton{};for(int j=0;j<nMinerals;j++)rhsNewton[j]=-residual[j];if(!linear(jac,rhsNewton,step))return false;
            bool accepted=false;
            for(double alpha=1;alpha>=1./16384;alpha*=.5){
                MineralVector y=x;for(int j=0;j<nMinerals;j++)y[j]+=alpha*step[j];auto s=changed(base,y);if(!valid(s))continue;
                try{auto e=evaluate(s,p,&last,false,incrementedPorosity(base,y,p));stats.rateEvaluations++;stats.speciationIterations+=e.aqueous.iterations;MineralVector r{};for(int j=0;j<nMinerals;j++)r[j]=y[j]-rhs[j]-coefficient*e.rates[j];
                    double newNorm=norm(r,y);
                    if(newNorm<oldNorm*(1-1e-4*alpha)||newNorm<1){
                        // Broyden's rank-one secant update saves seven speciation solves.
                        // The nonlinear residual is still checked against the same tolerance.
                        double denominator=0;MineralVector dx{};for(int j=0;j<nMinerals;j++){dx[j]=y[j]-x[j];denominator+=dx[j]*dx[j];}
                        if(denominator>1e-60)for(int i=0;i<nMinerals;i++){double correction=e.rates[i]-eval.rates[i];for(int j=0;j<nMinerals;j++)correction-=derivative[i][j]*dx[j];for(int j=0;j<nMinerals;j++)derivative[i][j]+=correction*dx[j]/denominator;}
                        haveJacobian=alpha>.125&&newNorm<oldNorm*.8;x=y;last=e.aqueous;cached=e;haveEvaluation=true;accepted=true;break;
                    }
                }catch(const std::exception&){}
            }
            if(!accepted)return false;
        }return false;
    }
    bool step(const Inventory& s,const RateParameters& p,double dt,MineralVector& extent,const AqueousState* seed,IntegrationStats& stats,const IntegrationOptions& options)const{
        constexpr double gamma=.5857864376269049512,a=.2928932188134524756,b=1.2071067811865475244;
        ReactionOutput initial;try{initial=evaluate(s,p,seed,false,incrementedPorosity(s,{},p));}catch(const std::exception&){return false;}stats.rateEvaluations++;stats.speciationIterations+=initial.aqueous.iterations;
        MineralVector rhs{},middle{};for(int i=0;i<nMinerals;i++)rhs[i]=gamma*dt*.5*initial.rates[i];
        if(!implicit(s,p,rhs,gamma*dt*.5,middle,&initial.aqueous,stats,options))return false;
        for(int i=0;i<nMinerals;i++)rhs[i]=b*middle[i];return implicit(s,p,rhs,a*dt,extent,&initial.aqueous,stats,options);
    }
public:
    explicit Reactions(const PitzerModel& m):model(m),spec(m){}
    void rateJacobian(const Inventory& state,const RateParameters& p,const ReactionOutput& value,double derivative[nMinerals][nMinerals],double porosity=std::numeric_limits<double>::quiet_NaN())const;
    static Inventory changed(const Inventory& initial,const MineralVector& extent){
        Inventory s=initial;
        long double water=initial.waterEquivalent;std::array<long double,nComponents> aqueous{};for(int i=0;i<nComponents;i++)aqueous[i]=initial.aqueous[i];
        for(int j=0;j<nMinerals;j++){s.mineral[j]-=extent[j];water+=static_cast<long double>(hydration[j])*waterMolarMass*extent[j];for(int i=0;i<nComponents;i++)aqueous[i]+=static_cast<long double>(stoichiometry[j][i])*extent[j];}
        s.waterEquivalent=double(water);for(int i=0;i<nComponents;i++)s.aqueous[i]=double(aqueous[i]);return s;
    }
    static void validateParameters(const RateParameters& p){
        if(p.geometry.separateMineralAreas&&!p.geometry.enabled)throw std::invalid_argument("Mineral interfaces require pore-volume constraints");
        for(int j=0;j<nMinerals;j++){
            if(!std::isfinite(p.geometry.growthArea[j])||p.geometry.growthArea[j]<0||!std::isfinite(p.geometry.surfaceExponent[j])||p.geometry.surfaceExponent[j]<0)throw std::invalid_argument("Invalid mineral interface area/exponent");
            if(!std::isfinite(p.geometry.lockedMoles[j])||p.geometry.lockedMoles[j]<0||(!p.geometry.separateMineralAreas&&p.geometry.lockedMoles[j]!=0))throw std::invalid_argument("Invalid inaccessible mineral inventory");
            if(!std::isfinite(p.geometry.curvatureLogOmega[j])||p.geometry.curvatureLogOmega[j]<0||
               (!p.geometry.separateMineralAreas&&p.geometry.curvatureLogOmega[j]!=0))throw std::invalid_argument("Invalid spherical-cap curvature correction");
            for(double v:{p.neutral[j],p.acid[j],p.co2[j],p.acidOrder[j],p.area[j],p.referenceMoles[j],p.dissolutionExponent[j],p.affinityExponent[j]})if(!std::isfinite(v))throw std::invalid_argument("Nonfinite kinetic parameter");
            if(p.neutral[j]<0||p.acid[j]<0||p.co2[j]<0||p.acidOrder[j]<0||p.area[j]<0||p.referenceMoles[j]<=0||p.dissolutionExponent[j]<=0||p.affinityExponent[j]<=0)throw std::invalid_argument("Invalid kinetic parameter for "+std::string(mineralNames[j]));
        }
        if(p.geometry.enabled){const auto& g=p.geometry;for(double v:{g.bulkVolume,g.dx,g.dy,g.directionWidth,g.closureWidth,g.epsilonMin,g.inertFraction})if(!std::isfinite(v))throw std::invalid_argument("Nonfinite phase geometry");if(!(g.bulkVolume>0&&g.dx>0&&g.dy>0&&g.directionWidth>0&&g.closureWidth>0&&g.epsilonMin>0&&g.epsilonMin<1&&g.inertFraction>=0&&g.inertFraction<1))throw std::invalid_argument("Invalid phase geometry");for(double w:g.weights)if(!std::isfinite(w)||w<0)throw std::invalid_argument("Invalid reactive surface weight");for(double e:g.neighbours)if(!std::isfinite(e)||e<0||e>1)throw std::invalid_argument("Invalid neighbouring porosity");}
    }
    ReactionOutput evaluate(const Inventory& state,const RateParameters& p,const AqueousState* seed=nullptr,bool computeDensity=true,double porosity=std::numeric_limits<double>::quiet_NaN())const{
        ReactionOutput r;r.aqueous=spec.conserved(state,r.freeWaterKg,seed,computeDensity);
        r.rates=surfaceRates(state,p,r.aqueous.saturationIndex,r.aqueous.logActivity[H],r.aqueous.logActivity[CO2],model.molarVolume,nullptr,nullptr,porosity);return r;
    }
    IntegrationStats advance(Inventory& state,double interval,const RateParameters& p,const IntegrationOptions& options={})const{
        if(!valid(state)||interval<0||!std::isfinite(interval)||!std::isfinite(options.absoluteTolerance)||!std::isfinite(options.relativeTolerance)||options.absoluteTolerance<=0||options.relativeTolerance<=0||options.maximumAttempts<=0)throw std::invalid_argument("Invalid kinetic initial state or tolerance");validateParameters(p);
        const double scale=state.waterEquivalent;Inventory normalized=state;RateParameters parameters=p;
        if(p.geometry.enabled&&!std::isfinite(parameters.geometry.initialPorosity))parameters.geometry.initialPorosity=poreFraction(state,p.geometry,model.molarVolume);
        MineralVector locked{};
        for(int j=0;j<nMinerals;j++){
            if(p.geometry.lockedMoles[j]>state.mineral[j])throw std::invalid_argument("Inaccessible mineral exceeds inventory");
            locked[j]=p.geometry.lockedMoles[j];normalized.mineral[j]-=locked[j];
            // All integrators already enforce positivity and handle phase
            // depletion. Integrate only accessible inventory, retaining the
            // buried volume as inert geometry until this local step ends.
            // No aqueous or hydrate-water transfer accompanies this shift.
            if(locked[j]>0)parameters.geometry.inertFraction+=locked[j]*model.molarVolume[j]/p.geometry.bulkVolume;
            parameters.geometry.lockedMoles[j]=0;
        }
        for(int j=0;j<nComponents;j++)normalized.aqueous[j]/=scale;
        for(int j=0;j<nMinerals;j++){normalized.mineral[j]/=scale;parameters.area[j]/=scale;parameters.geometry.growthArea[j]/=scale;parameters.referenceMoles[j]/=scale;}
        parameters.geometry.bulkVolume/=scale;
        normalized.waterEquivalent=1;IntegrationStats stats;
        if(options.method==0)stats=advanceNormalized(normalized,interval,parameters,options);
        else if(options.method==1)stats=advanceRosenbrock(normalized,interval,parameters,options);
#ifdef SALT_HAS_CVODE
        else if(options.method==2)stats=advanceCvode(normalized,interval,parameters,options);
#endif
        else throw std::invalid_argument("Unsupported chemical integration method");
        // An exact no-reaction solve must remain an exact inventory no-op.
        // Repeated divide/multiply normalization otherwise accumulates ulps
        // even when a growth-only phase stays undersaturated throughout.
        if(std::all_of(stats.dissolved.begin(),stats.dissolved.end(),[](double amount){return amount==0;})&&
           std::all_of(stats.precipitated.begin(),stats.precipitated.end(),[](double amount){return amount==0;}))return stats;
        for(int j=0;j<nComponents;j++)normalized.aqueous[j]*=scale;
        for(int j=0;j<nMinerals;j++){normalized.mineral[j]=normalized.mineral[j]*scale+locked[j];stats.dissolved[j]*=scale;stats.precipitated[j]*=scale;}
        normalized.waterEquivalent*=scale;
        for(int j=0;j<nMinerals;j++)if(stats.dissolved[j]==0&&stats.precipitated[j]==0){
            const double roundoff=64*std::numeric_limits<double>::epsilon()*std::max(state.mineral[j],normalized.mineral[j]);
            if(std::abs(normalized.mineral[j]-state.mineral[j])>roundoff)
                throw std::runtime_error("Unreacted mineral changed beyond normalization roundoff");
            // No reaction extent was applied to this phase. Its normalization
            // roundoff did not transfer ions/water, so restore only its stock.
            normalized.mineral[j]=state.mineral[j];
        }
        if(p.geometry.separateMineralAreas){
            // Dense implicit solves can leave roundoff in an exactly inactive
            // extent. Enforce the exact zero-area invariant, with the paired
            // aqueous/water correction, before it can seed a new interface.
            MineralVector correction{};bool repair=false;
            for(int j=0;j<nMinerals;j++)if(p.area[j]==0&&p.geometry.growthArea[j]==0){
                correction[j]=normalized.mineral[j]-state.mineral[j];repair|=correction[j]!=0;
                if(std::abs(correction[j])>64*(scale*options.absoluteTolerance+options.relativeTolerance*state.mineral[j]))throw std::runtime_error("Inactive mineral changed beyond roundoff");
                stats.dissolved[j]=stats.precipitated[j]=0;
            }
            if(repair)normalized=changed(normalized,correction);
        }
        state=normalized;return stats;
    }
private:
#ifdef SALT_HAS_CVODE
    IntegrationStats advanceCvode(Inventory& state,double interval,const RateParameters& p,const IntegrationOptions& options)const;
#endif
    // Boost.Odeint's embedded fourth-order Rosenbrock integrator uses one local
    // Jacobian factorization for all stages. Thermodynamics and TST remain ours.
    // Stage positivity is enforced before evaluating any chemical properties.
    IntegrationStats advanceRosenbrock(Inventory& state,double interval,const RateParameters& inputParameters,const IntegrationOptions& options)const{
        using Stepper=boost::numeric::odeint::rosenbrock4<double>;
        using Vector=Stepper::state_type;using Matrix=Stepper::matrix_type;
        IntegrationStats stats;if(interval==0)return stats;
        RateParameters p=inputParameters;const Inventory origin=state;MineralVector cumulative{},compensation{};
        const double originPorosity=incrementedPorosity(origin,{},p);
        auto commit=[&](const MineralVector& extent){for(int j=0;j<nMinerals;j++){const double v=extent[j]-compensation[j],t=cumulative[j]+v;compensation[j]=(t-cumulative[j])-v;cumulative[j]=t;}state=changed(origin,cumulative);auto base=inputParameters;base.geometry.initialPorosity=originPorosity;p.geometry.initialPorosity=incrementedPorosity(origin,cumulative,base);};
        Stepper stepper;Vector zero(nMinerals,0),extent(nMinerals),estimate(nMinerals);
        double time=0,dt=interval;int attempts=0,consecutiveFailures=0;
        AqueousState seed=evaluate(state,p,nullptr,false,originPorosity).aqueous;
        while(time<interval*(1-1e-13)){
            if(++attempts>options.maximumAttempts){std::ostringstream message;message<<std::setprecision(15)<<"Rosenbrock attempt limit: t="<<time<<" dt="<<dt<<" accepted="<<stats.accepted<<" rejected="<<stats.rejected<<" epsilon="<<(p.geometry.enabled?poreFraction(state,p.geometry,model.molarVolume):1)<<" water="<<state.waterEquivalent;throw std::runtime_error(message.str());}
            auto eval=[&](const MineralVector& x){auto candidate=changed(state,x);if(!valid(candidate))throw std::domain_error("Chemical stage violates positivity");auto output=evaluate(candidate,p,&seed,false,incrementedPorosity(state,x,p));seed=output.aqueous;stats.rateEvaluations++;stats.speciationIterations+=output.aqueous.iterations;return output;};
            auto initial=eval({});MineralVector exhausted{};bool event=false;
            for(int j=0;j<nMinerals;j++)if(state.mineral[j]>0&&state.mineral[j]<=options.absoluteTolerance&&initial.rates[j]>0){exhausted[j]=state.mineral[j];stats.dissolved[j]+=exhausted[j];stats.depletionEvents++;event=true;}
            if(event){commit(exhausted);initial=eval({});}
            // With no local forcing and every rate exactly zero, the autonomous
            // ODE has a constant solution over the remainder of this interval.
            if(std::all_of(initial.rates.begin(),initial.rates.end(),[](double r){return r==0;}))return stats;
            dt=std::min(dt,interval-time);
            auto rhs=[&](const Vector& x,Vector& out,double){MineralVector v{};for(int j=0;j<nMinerals;j++)v[j]=x[j];auto result=eval(v);for(int j=0;j<nMinerals;j++)out[j]=result.rates[j];};
            auto jac=[&](const Vector&,Matrix& out,double,Vector& dfdt){
                for(int j=0;j<nMinerals;j++){
                    dfdt[j]=0;for(int i=0;i<nMinerals;i++)out(i,j)=0;
                    if((!p.geometry.enabled&&p.area[j]==0)||(p.geometry.enabled&&p.geometry.weights[j]==0))continue;
                    // The derivative is with respect to dissolved reaction
                    // extent. A large solid inventory in a nearly closed pore
                    // must not produce an order-one aqueous perturbation.
                    double delta=1e-6*std::max(std::min(state.mineral[j],state.waterEquivalent),1e-6*state.waterEquivalent);MineralVector plus{},minus{};bool ok=false;
                    for(int k=0;k<24;k++){
                        plus[j]=delta;minus[j]=-delta;bool vp=admissible(changed(state,plus),p),vm=admissible(changed(state,minus),p);
                        if(vp||vm){auto a=vp?eval(plus):initial,b=vm?eval(minus):initial;double denom=(vp&&vm?2:1)*delta;
                            for(int i=0;i<nMinerals;i++)out(i,j)=(a.rates[i]-b.rates[i])/denom;ok=true;break;}
                        delta*=.1;
                    }
                    if(!ok&&state.mineral[j]>0)throw std::runtime_error("No feasible chemical Jacobian perturbation");
                }
            };
            double error=0;bool ok=true;MineralVector change{};
            try{
                stepper.do_step(std::make_pair(rhs,jac),zero,time,extent,dt,estimate);
                for(int j=0;j<nMinerals;j++){
                    change[j]=extent[j];double scale=options.absoluteTolerance+options.relativeTolerance*std::max(std::abs(state.mineral[j]),std::abs(change[j]));
                    error=std::max(error,std::abs(estimate[j])/scale);
                }
                if(!admissible(changed(state,change),p))ok=false;
                for(int i=0;i<nComponents;i++){double dc=0,de=0;for(int j=0;j<nMinerals;j++){dc+=stoichiometry[j][i]*change[j];de+=stoichiometry[j][i]*estimate[j];}double scale=options.absoluteTolerance+options.relativeTolerance*std::max(std::abs(state.aqueous[i]),std::abs(dc));error=std::max(error,std::abs(de)/scale);}
                double dw=0,ew=0;for(int j=0;j<nMinerals;j++){dw+=hydration[j]*waterMolarMass*change[j];ew+=hydration[j]*waterMolarMass*estimate[j];}error=std::max(error,std::abs(ew)/(waterMolarMass*options.absoluteTolerance+options.relativeTolerance*std::max(state.waterEquivalent,std::abs(dw))));
                if(!std::isfinite(error))ok=false;
            }catch(const std::exception&){ok=false;}
            if(!ok||error>1){
                seed=initial.aqueous;
                stats.rejected++;consecutiveFailures++;dt*=ok?std::clamp(.85*std::pow(error,-.25),.1,.8):.25;
                if(consecutiveFailures>=10||dt<interval*1e-14){
                    auto tail=advanceNormalized(state,interval-time,p,options);stats.accepted+=tail.accepted;stats.rejected+=tail.rejected;stats.rateEvaluations+=tail.rateEvaluations;stats.speciationIterations+=tail.speciationIterations;stats.newtonIterations+=tail.newtonIterations;stats.depletionEvents+=tail.depletionEvents;
                    for(int j=0;j<nMinerals;j++){stats.dissolved[j]+=tail.dissolved[j];stats.precipitated[j]+=tail.precipitated[j];}return stats;
                }continue;
            }
            commit(change);time+=dt;stats.accepted++;consecutiveFailures=0;
            for(int j=0;j<nMinerals;j++){stats.dissolved[j]+=std::max(change[j],0.);stats.precipitated[j]+=std::max(-change[j],0.);}
            dt*=error>1e-16?std::clamp(.9*std::pow(error,-.25),.3,4.):4.;
        }return stats;
    }
    IntegrationStats advanceNormalized(Inventory& state,double interval,const RateParameters& inputParameters,const IntegrationOptions& options)const{
        IntegrationStats stats;if(interval==0)return stats;double time=0,dt=interval;int attempts=0;
        RateParameters p=inputParameters;const Inventory origin=state;MineralVector cumulative{},compensation{};
        const double originPorosity=incrementedPorosity(origin,{},p);
        auto commit=[&](const MineralVector& extent){for(int j=0;j<nMinerals;j++){const double v=extent[j]-compensation[j],t=cumulative[j]+v;compensation[j]=(t-cumulative[j])-v;cumulative[j]=t;}state=changed(origin,cumulative);auto base=inputParameters;base.geometry.initialPorosity=originPorosity;p.geometry.initialPorosity=incrementedPorosity(origin,cumulative,base);};
        AqueousState seed=evaluate(state,p,nullptr,false,originPorosity).aqueous;
        while(time<interval*(1-1e-13)){
            // Resolve a disappearing phase within the requested absolute accuracy.
            // Transfer its residual inventory stoichiometrically; never discard mass.
            // This avoids asymptotic step collapse for area exponents below one.
            auto endpoint=evaluate(state,p,&seed,false,incrementedPorosity(state,{},p));stats.rateEvaluations++;stats.speciationIterations+=endpoint.aqueous.iterations;
            MineralVector exhausted{};bool event=false;
            for(int j=0;j<nMinerals;j++)if(state.mineral[j]>0&&state.mineral[j]<=options.absoluteTolerance&&endpoint.rates[j]>0){exhausted[j]=state.mineral[j];stats.dissolved[j]+=exhausted[j];stats.depletionEvents++;event=true;}
            if(event){commit(exhausted);seed=evaluate(state,p,&seed,false,incrementedPorosity(state,{},p)).aqueous;}
            if(++attempts>options.maximumAttempts)throw std::runtime_error("Seven-mineral integration exceeded attempt limit");dt=std::min(dt,interval-time);
            MineralVector coarse{},first{},second{};bool ok=step(state,p,dt,coarse,&seed,stats,options)&&step(state,p,dt*.5,first,&seed,stats,options);
            auto midpoint=changed(state,first);auto midpointParameters=p;midpointParameters.geometry.initialPorosity=incrementedPorosity(state,first,p);if(ok)ok=step(midpoint,midpointParameters,dt*.5,second,&seed,stats,options);
            MineralVector fine{};double error=0;
            if(ok)for(int j=0;j<nMinerals;j++){
                fine[j]=first[j]+second[j];double scale=options.absoluteTolerance+options.relativeTolerance*std::max(std::abs(state.mineral[j]),std::abs(fine[j]));
                error=std::max(error,std::abs(fine[j]-coarse[j])/(3*scale));
            }
            // Also constrain errors in aqueous inventories, which may be much smaller than solids.
            if(ok)for(int i=0;i<nComponents;i++){double difference=0,change=0;for(int j=0;j<nMinerals;j++){difference+=stoichiometry[j][i]*(fine[j]-coarse[j]);change+=stoichiometry[j][i]*fine[j];}double scale=options.absoluteTolerance+options.relativeTolerance*std::max(std::abs(state.aqueous[i]),std::abs(change));error=std::max(error,std::abs(difference)/(3*scale));}
            if(!ok||error>1){stats.rejected++;dt*=ok?std::max(.1,.8*std::pow(error,-1./3)):.25;if(dt<interval*1e-14)throw std::runtime_error("Seven-mineral implicit solve failed near t="+std::to_string(time));continue;}
            commit(fine);seed=evaluate(state,p,&seed,false,incrementedPorosity(state,{},p)).aqueous;time+=dt;stats.accepted++;
            for(int j=0;j<nMinerals;j++){stats.dissolved[j]+=std::max(first[j],0.)+std::max(second[j],0.);stats.precipitated[j]+=std::max(-first[j],0.)+std::max(-second[j],0.);}
            dt*=error>1e-12?std::clamp(.9*std::pow(error,-1./3),.5,2.5):2.5;
        }return stats;
    }
};
}
#include "BrineJacobian.hpp"
#ifdef SALT_HAS_CVODE
#include "BrineCvode.hpp"
#endif
