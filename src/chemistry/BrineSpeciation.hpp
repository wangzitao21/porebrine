// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "BrinePitzer.hpp"
namespace saltchem {
class Speciation {
    const PitzerModel& model;
    template<class Real> struct Distribution {
        std::array<Real,nComponents> totals;std::array<Real,nSpecies> gamma;Real kw,bicarbonate,co2,bisulfate,mgoh,mgco3;
        explicit Distribution(const PitzerModel& p,const std::array<Real,nComponents>& t,const std::array<Real,nSpecies>& g,Real lnAw,Real I):totals(t),gamma(g){
            const auto logK=p.aqueousLogK(I,checkedExp(lnAw));
            const Real la=lnAw/ln10;
            kw=pow10(logK[OH]+la-g[OH]);
            bicarbonate=pow10(logK[Bicarbonate]+g[Carbonate]-g[Bicarbonate]);
            co2=pow10(logK[CO2]+g[Carbonate]-la-g[CO2]);
            bisulfate=pow10(logK[Bisulfate]+g[Sulfate]-g[Bisulfate]);
            mgoh=pow10(logK[MgOH]+g[MgIon]+la-g[MgOH]);
            mgco3=pow10(logK[MgCO3]+g[MgIon]+g[Carbonate]-g[MgCO3]);
        }
        std::array<Real,nSpecies> at(Real logH)const{
            std::array<Real,nSpecies> m{};const Real ah=pow10(logH),rH=bicarbonate*ah,rC=co2*ah*ah,rS=bisulfate*ah,rM=mgoh/ah;
            m[NaIon]=totals[Na];m[KIon]=totals[K];m[CaIon]=totals[Ca];m[ClIon]=totals[Cl];
            m[H]=pow10(logH-gamma[H]);m[OH]=kw/ah;
            m[Sulfate]=totals[S]/(1+rS);m[Bisulfate]=m[Sulfate]*rS;
            Real assoc=0;
            if(scalarValue(totals[C])>=0&&scalarValue(totals[Mg])>=0){
                const Real t=mgco3/((1+rM)*(1+rH+rC)),sum=totals[Mg]+totals[C];
                if(scalarValue(t)>0){
                    // Equivalent roots scaled for weak and strong association.
                    // Avoid squaring a huge 1/t at small association constants.
                    if(scalarValue(t)*scalarValue(sum)<1){
                        const Real d=t*(totals[Mg]-totals[C]);
                        assoc=2*t*totals[Mg]*totals[C]/(1+t*sum+sqrt(1+2*t*sum+d*d));
                    }else{
                        const Real inv=1/t;assoc=2*totals[Mg]*totals[C]/(sum+inv+hypot(totals[Mg]-totals[C],sqrt(inv*(2*sum+inv))));
                    }
                }
            }
            m[MgCO3]=assoc;m[MgIon]=(totals[Mg]-assoc)/(1+rM);m[MgOH]=m[MgIon]*rM;
            m[Carbonate]=(totals[C]-assoc)/(1+rH+rC);m[Bicarbonate]=m[Carbonate]*rH;m[CO2]=m[Carbonate]*rC;return m;
        }
        double residual(const SpeciesVector& m)const{
            // Work with conserved strong-ion difference to avoid repeated cancellation.
            long double q=static_cast<long double>(totals[Na])+totals[K]+2*static_cast<long double>(totals[Ca])+2*static_cast<long double>(totals[Mg])-totals[Cl]-2*static_cast<long double>(totals[S])-2*static_cast<long double>(totals[C]);
            return static_cast<double>(q+m[H]-m[OH]+m[Bicarbonate]+2*m[CO2]+m[Bisulfate]-m[MgOH]);
        }
        double derivative(const SpeciesVector& m,double logH)const{
            const double ah=pow10(logH),rH=bicarbonate*ah,rC=co2*ah*ah,rS=bisulfate*ah,rM=mgoh/ah;
            const double u=1+rM,v=1+rH+rC,dc=(rH+2*rC)/v;
            double da=0,freeMg=totals[Mg]-m[MgCO3],freeC=totals[C]-m[MgCO3];
            if(m[MgCO3]>0&&freeMg>0&&freeC>0)da=m[MgCO3]*(rM/u-dc)/(1+m[MgCO3]/freeMg+m[MgCO3]/freeC);
            const double dCarb=-da/v-m[Carbonate]*dc;
            const double dMgOH=-da*rM/u-m[MgOH]/u;
            return ln10*(m[H]+m[OH]+dCarb*(rH+2*rC)+m[Bicarbonate]+4*m[CO2]+m[Bisulfate]/(1+rS)-dMgOH);
        }
        SpeciesVector electroneutral(double& logH)const{
            double lo=-30,hi=10;logH=std::clamp(logH,lo,hi);SpeciesVector m{};
            for(int it=0;it<120;it++){
                m=at(logH);double r=residual(m);if(r>0)hi=logH;else lo=logH;
                if(hi-lo<2e-14||r==0)return m;
                double slope=derivative(m,logH);
                double next=logH-r/slope;
                if(std::isfinite(next)&&std::abs(next-logH)<=2*std::numeric_limits<double>::epsilon()*std::max(std::abs(logH),1.))return m;
                if(!std::isfinite(next)||next<=lo||next>=hi)next=.5*(lo+hi);
                // Force bracket contraction periodically even if Newton repeatedly
                // takes valid but tiny steps along one edge of the bracket.
                if(it%4==3)next=.5*(lo+hi);
                logH=next;
            }
            throw std::runtime_error("Acid-base charge solve did not converge");
        }
    };
public:
    explicit Speciation(const PitzerModel& p):model(p){}
    template<class Real> static std::array<Real,nSpecies> distribute(const PitzerModel& p,const std::array<Real,nComponents>& t,const std::array<Real,nSpecies>& g,Real lnAw,Real I,Real logH){return Distribution<Real>(p,t,g,lnAw,I).at(logH);}
    AqueousState solve(const Totals& totals,const AqueousState* seed=nullptr,double tolerance=2e-13,bool computeDensity=true)const{
        for(double x:totals)if(x<0||!std::isfinite(x))throw std::invalid_argument("Component totals must be finite nonnegative mol/kg free water");
        SpeciesVector lg{};double lnAw=0,logH=-7,I=0;for(int i=0;i<nComponents;i++)I+=.5*charge[i]*charge[i]*totals[i];
        if(seed){lg=seed->logGamma;lnAw=std::log(seed->waterActivity);logH=-seed->pH;I=seed->ionicStrength;}
        AqueousState out;double last=1e99,relax=1;
        for(int it=0;it<300;it++){
            Distribution<double> distribution(model,totals,lg,lnAw,I);out.molality=distribution.electroneutral(logH);
            SpeciesVector next{};double nextAw,nextI;model.activities(out.molality,next,nextAw,nextI,checkedExp(lnAw));
            double error=std::abs(nextAw-lnAw);for(int j=0;j<nSpecies;j++)error=std::max(error,std::abs(next[j]-lg[j]));
            if(model.pressureAtm()>1)error=std::max(error,std::abs(nextI-I)/(1+nextI));
            if(error<tolerance){out.iterations=it+1;out.logGamma=next;out.ionicStrength=nextI;model.finish(out,nextAw,computeDensity);auto t=componentTotals(out.molality);for(int j=0;j<nComponents;j++)out.massResidual=std::max(out.massResidual,std::abs(t[j]-totals[j])/(1+totals[j]));out.equilibriumResidual=error;return out;}
            if(it>2&&error>last*1.05)relax=std::max(.08,relax*.5);else if(error<last*.7)relax=std::min(1.,relax*1.1);
            for(int j=0;j<nSpecies;j++)lg[j]+=relax*(next[j]-lg[j]);lnAw+=relax*(nextAw-lnAw);I+=relax*(nextI-I);last=error;
        }
        std::ostringstream msg;msg<<"Pitzer speciation did not converge; last log-activity correction="<<std::scientific<<last;throw std::runtime_error(msg.str());
    }
    // Free-water based initial state -> conserved O-equivalent water inventory.
    Inventory inventory(const Totals& totals,double freeWaterKg,const MineralVector& minerals={})const{
        if(freeWaterKg<=0)throw std::invalid_argument("Water must be positive");
        const auto state=solve(totals);Inventory v;v.mineral=minerals;
        for(int i=0;i<nComponents;i++)v.aqueous[i]=totals[i]*freeWaterKg;
        v.waterEquivalent=freeWaterKg*(1+waterMolarMass*(state.molality[OH]+state.molality[MgOH]-state.molality[CO2]));return v;
    }
    AqueousState conserved(const Inventory& v,double& freeWaterKg,const AqueousState* seed=nullptr,bool computeDensity=true)const{
        if(v.waterEquivalent<=0||!std::isfinite(v.waterEquivalent))throw std::invalid_argument("Nonpositive water-equivalent inventory");
        freeWaterKg=v.waterEquivalent;
        if(seed)freeWaterKg/=1+waterMolarMass*(seed->molality[OH]+seed->molality[MgOH]-seed->molality[CO2]);
        AqueousState out;bool have=false;
        for(int it=0;it<40;it++){
            Totals t{};for(int j=0;j<nComponents;j++)t[j]=v.aqueous[j]/freeWaterKg;
            out=solve(t,have?&out:seed,2e-13,computeDensity);have=true;
            double corrected=v.waterEquivalent/(1+waterMolarMass*(out.molality[OH]+out.molality[MgOH]-out.molality[CO2]));
            if(std::abs(corrected/freeWaterKg-1)<2e-12)return out;
            freeWaterKg=corrected;
        }
        throw std::runtime_error("Solvent mass balance did not converge");
    }
};
}
