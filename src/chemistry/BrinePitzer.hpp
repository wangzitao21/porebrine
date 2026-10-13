// SPDX-License-Identifier: GPL-3.0-or-later
// Fixed chemical system; PHREEQC-compatible temperature/pressure dependence.
// Reference: USGS public-domain phreeqcpp/pitzer.cpp, prep.cpp, utilities.cpp.
#pragma once
#include "BrineSystem.hpp"
#include "BrineThermal.hpp"
#include "DatabaseContract.hpp"
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <limits>
namespace saltchem {
class PitzerModel {
    enum Kind {B0,B1,B2,C0,Theta,Lambda,Psi,Zeta};
    struct Term {Kind kind;int i,j,k=-1,group=-1;double value=0,alpha=0,normalized=0;std::array<double,6> coefficients{};};
    std::vector<Term> terms;
    std::vector<double> activityAlphas;
    std::vector<std::array<int,2>> differentChargePairs;
    std::array<std::array<double,10>,nSpecies> vm{};
    SpeciesVector v0{},bi{};
    WaterProperties water;
    std::array<EquilibriumConstant,nSpecies> aqueousConstants{};
    std::array<EquilibriumConstant,nMinerals> phaseConstants{};
    double temperature=298.15,pressure=1;
    double macB0=0,macB1=0,macC0=0;
    bool strictWater=false;
    template<class Real> static std::array<Real,2> jay(Real x){
        static constexpr double a[42]={1.925154014814667,-.060076477753119,-.029779077456514,-.007299499690937,.000388260636404,.000636874599598,.000036583601823,-.000045036975204,-.000004537895710,.000002937706971,.000000396566462,-.000000202099617,-.000000025267769,.000000013522610,.000000001229405,-.000000000821969,-.000000000050847,.000000000046333,.000000000001943,-.000000000002563,-.000000000010991,.628023320520852,.462762985338493,.150044637187895,-.028796057604906,-.036552745910311,-.001668087945272,.006519840398744,.001130378079086,-.000887171310131,-.000242107641309,.000087294451594,.000034682122751,-.000004583768938,-.000003548684306,-.000000250453880,.000000216991779,.000000080779570,.000000004558555,-.000000006944757,-.000000002849257,.000000000237816};
        const bool low=scalarValue(x)<=1;Real u=low?pow(x,.2):pow(x,-.1),z=low?4*u-2:(40*u-22)/9,dz=low?.4*u:-4*u/18;
        const double* p=a+(low?0:21);Real b[23]={},d[23]={};b[20]=p[20];b[19]=z*p[20]+p[19];d[19]=p[20];
        for(int i=18;i>=0;i--){b[i]=z*b[i+1]-b[i+2]+p[i];d[i]=b[i+1]+z*d[i+1]-d[i+2];}
        return {{x/4-1+.5*(b[0]-b[2]),x*.25+dz*(d[0]-d[2])}};
    }
    template<class Real> static Real G(Real x){return scalarValue(x)<1e-3?1-2*x/3+x*x/4-x*x*x/15:2*(1-(1+x)*exp(-x))/(x*x);}
    template<class Real> static Real GP(Real x){return scalarValue(x)<1e-3?-x/3+x*x/4-x*x*x/10:-2*(1-(1+x+x*x/2)*exp(-x))/(x*x);}
    static double numeric(std::istringstream& in){double value;if(!(in>>value)||!std::isfinite(value))throw std::runtime_error("Malformed or nonfinite thermodynamic parameter");return value;}
    template<std::size_t N> static void coefficients(std::istringstream& in,std::array<double,N>& values){
        std::size_t count=0;in>>std::ws;
        while(!in.eof()){if(count==N)throw std::runtime_error("Too many thermodynamic coefficients");values[count++]=numeric(in);in>>std::ws;}
        if(count==0)throw std::runtime_error("Missing thermodynamic coefficient");
    }
    static void constantInput(EquilibriumConstant& k,const std::string& name,std::istringstream& in){
        if(name=="log_k"||name=="-log_k")k.reference=numeric(in);
        if(name=="delta_h"||name=="-delta_h"){
            std::string unit;k.enthalpyKJ=numeric(in);in>>unit;
            if(unit=="kcal")k.enthalpyKJ*=4.184;
            else if(unit=="cal")k.enthalpyKJ*=.004184;
            else if(unit=="J"||unit=="joules")k.enthalpyKJ*=.001;
            else if(!unit.empty()&&unit!="kJ"&&unit!="kj")throw std::runtime_error("Unsupported enthalpy unit: "+unit);
        }
        if(name.rfind("-analytic",0)==0){k.analytic=true;coefficients(in,k.coefficients);}
    }
    template<class Real> std::array<Real,nSpecies> volumes(Real I,const WaterProperties::StateType<Real>& w)const {
        std::array<Real,nSpecies> v;for(int i=0;i<nSpecies;i++)v[i]=v0[i];const Real positiveI=scalarValue(I)>0?I:Real(0),root=sqrt(positiveI);
        for(int i=0;i<nSpecies;i++){
            const auto& a=vm[i];
            if(scalarValue(w.pressure)!=water.dielectricPressure)v[i]=a[0]*4.184004+a[1]*4184.004/(2600+w.pressure*1.01325)+(a[2]*41.84004+a[3]*418400.4/(2600+w.pressure*1.01325))/(temperature-228)-a[4]*1e5*w.qBorn;
            if(a[0]!=0&&charge[i]!=0)v[i]+=.5*charge[i]*charge[i]*w.av*root/(1+a[5]*w.dhb*root)+bi[i]*pow(positiveI,a[9]);
        }return v;
    }
    template<class Real> Real pressureFactor(const WaterProperties::StateType<Real>& w)const {return (scalarValue(w.pressure)>1?w.pressure-1:Real(0))*1.01325e-4/(ln10*.00831470*temperature);}
public:
    SpeciesVector logK{};
    MineralVector mineralLogK{},molarVolume{};
    explicit PitzerModel(const std::string& path,double temperatureK=298.15,double pressureAtm=1){
        std::ifstream f(path);if(!f)throw std::runtime_error("Cannot read Pitzer database: "+path);
        std::string line,section,kind;int species=-1,mineral=-1;std::array<bool,nSpecies> hasK{};std::array<bool,nMinerals> hasPhaseK{};
        struct Alpha {int i,j;double a,b;};std::vector<Alpha> alphas;
        std::string database((std::istreambuf_iterator<char>(f)),{});std::replace(database.begin(),database.end(),';','\n');DatabaseContract::validate(database);std::istringstream lines(database);
        while(std::getline(lines,line)){
            line=line.substr(0,line.find('#'));std::istringstream in(line);std::string first;if(!(in>>first))continue;
            if(first=="SOLUTION_MASTER_SPECIES"||first=="SOLUTION_SPECIES"||first=="PHASES"||first=="PITZER"||first=="EXCHANGE_MASTER_SPECIES"||first=="EXCHANGE_SPECIES"||first=="SURFACE_MASTER_SPECIES"||first=="SURFACE_SPECIES"||first=="RATES"){section=first;continue;}
            if(section=="SOLUTION_SPECIES"){
                const auto eq=line.find('=');if(eq!=std::string::npos){species=DatabaseContract::aqueousIndex(line);continue;}
                if(species<0)continue;
                if(first=="-Vm"){vm[species][9]=1;coefficients(in,vm[species]);}
                constantInput(aqueousConstants[species],first,in);
                if(first=="log_k"||first=="-log_k"||first.rfind("-analytic",0)==0)hasK[species]=true;
            }else if(section=="PHASES"){
                if(!std::isspace(static_cast<unsigned char>(line[0]))){mineral=mineralIndex(first);continue;}
                if(mineral<0)continue;
                if(first=="-Vm"||first=="Vm"){molarVolume[mineral]=numeric(in)*1e-6;}
                constantInput(phaseConstants[mineral],first,in);
                if(first=="log_k"||first=="-log_k"||first.rfind("-analytic",0)==0)hasPhaseK[mineral]=true;
            }else if(section=="PITZER"){
                if(first[0]=='-'){kind=first;continue;}
                int i=speciesIndex(first);std::string second,third;in>>second;int j=speciesIndex(second);if(i<0||j<0)continue;
                Term t{};t.i=i;t.j=j;
                if(kind=="-PSI"||kind=="-ZETA"||kind=="-MU"||kind=="-ETA"){in>>third;t.k=speciesIndex(third);if(t.k<0)continue;}
                if(kind=="-alphas"||kind=="-ALPHAS"){double a=numeric(in),b=numeric(in);if(a<=0||b<=0)throw std::runtime_error("Nonpositive Pitzer alpha");alphas.push_back({i,j,a,b});continue;}
                coefficients(in,t.coefficients);
                if(kind=="-B0")t.kind=B0;else if(kind=="-B1")t.kind=B1;else if(kind=="-B2")t.kind=B2;else if(kind=="-C0")t.kind=C0;else if(kind=="-THETA")t.kind=Theta;else if(kind=="-LAMBDA")t.kind=Lambda;else if(kind=="-PSI")t.kind=Psi;else if(kind=="-ZETA")t.kind=Zeta;else throw std::runtime_error("Unsupported active Pitzer parameter: "+kind);
                if(t.kind==B1)t.alpha=std::abs(charge[i])==2&&std::abs(charge[j])==2?1.4:2;
                if(t.kind==B2)t.alpha=12;
                // PHREEQC replaces duplicate parameter definitions, independent of ordering.
                auto key=[](const Term& p){std::array<int,3> k{{p.i,p.j,p.k}};std::sort(k.begin(),k.end());return k;};bool replaced=false;
                for(auto& old:terms)if(old.kind==t.kind&&key(old)==key(t)){old=t;replaced=true;break;}
                if(!replaced)terms.push_back(t);
            }
        }
        for(auto a:alphas)for(auto& t:terms)if((a.i==t.i&&a.j==t.j)||(a.i==t.j&&a.j==t.i)){if(t.kind==B1)t.alpha=a.a;if(t.kind==B2)t.alpha=a.b;}
        for(int i:{OH,Bicarbonate,CO2,Bisulfate,MgOH,MgCO3})if(!hasK[i])throw std::runtime_error("Missing aqueous equilibrium: "+std::string(speciesNames[i]));
        for(int i=0;i<nMinerals;i++)if(!hasPhaseK[i]||molarVolume[i]<=0)throw std::runtime_error("Missing mineral data: "+std::string(mineralNames[i]));
        for(auto& t:terms){
            if(t.kind==B1||t.kind==B2){auto p=std::find(activityAlphas.begin(),activityAlphas.end(),t.alpha);t.group=int(p-activityAlphas.begin());if(p==activityAlphas.end())activityAlphas.push_back(t.alpha);}
        }
        if(activityAlphas.size()>128)throw std::runtime_error("Too many Pitzer alpha groups");
        for(int i=0;i<nSpecies;i++)for(int j=i+1;j<nSpecies;j++)if(charge[i]*charge[j]>0&&charge[i]!=charge[j])differentChargePairs.push_back({i,j});
        setConditions(temperatureK,pressureAtm,true);
    }
    double temperatureK()const{return temperature;}
    double pressureAtm()const{return pressure;}
    struct WaterState {double requestedTemperatureK,waterTemperatureK,requestedPressureAtm,effectivePressureAtm;bool temperatureCapped,pressureRaised;};
    WaterState waterState(double aw)const{
        const auto w=water.at(aw);return {temperature,water.waterT,pressure,w.pressure,water.waterT!=temperature,w.pressure>pressure};
    }
    void setStrictWater(bool strict){
        if(strict&&water.waterT!=temperature)throw std::domain_error("PHREEQC water-property temperature cap would be used");strictWater=strict;
    }
    // Mutable conditions belong to one worker. Reuse all temperature-only
    // coefficients when visiting cells with different pressures at the same T.
    bool setConditions(double T,double pa,bool force=false){
        if(!force&&T==temperature&&pa==pressure)return false;
        const bool thermalChanged=force||T!=temperature;
        WaterProperties next=water;next.set(T,pa);if(strictWater&&next.waterT!=T)throw std::domain_error("PHREEQC water-property temperature cap would be used");water=next;temperature=T;pressure=pa;
        const double pb=water.dielectricPressure*1.01325,dT=T-228;
        for(int i=0;i<nSpecies;i++){
            const auto& a=vm[i];v0[i]=a[0]*4.184004+a[1]*4184.004/(2600+pb)+(a[2]*41.84004+a[3]*418400.4/(2600+pb))/dT-a[4]*1e5*water.qBorn;
            if(thermalChanged){bi[i]=a[6]+a[7]/dT+a[8]*dT;logK[i]=aqueousConstants[i].at(T);}
        }
        if(!thermalChanged)return true;
        for(int j=0;j<nMinerals;j++)mineralLogK[j]=phaseConstants[j].at(T);
        const double inv=1/T-1/298.15,l=std::log(T/298.15),dt=T-298.15,dt2=T*T-298.15*298.15,inv2=1/(T*T)-1/(298.15*298.15);
        for(auto& t:terms){const auto& a=t.coefficients;t.value=std::abs(dt)<.001?a[0]:a[0]+a[1]*inv+a[2]*l+a[3]*dt+a[4]*dt2+a[5]*inv2;
            if(t.kind==C0)t.normalized=t.value/(2*std::sqrt(std::abs(charge[t.i]*charge[t.j])));
            if((t.i==KIon&&t.j==ClIon)||(t.i==ClIon&&t.j==KIon)){if(t.kind==B0)macB0=t.value;if(t.kind==B1)macB1=t.value;if(t.kind==C0)macC0=t.value;}
        }return true;
    }
    template<class Real> std::array<Real,nSpecies> aqueousLogK(Real I,Real aw)const {
        std::array<Real,nSpecies> k;for(int i=0;i<nSpecies;i++)k[i]=logK[i];if(pressure<=1&&water.pSat<=1)return k;
        const auto w=water.at(aw);const Real f=pressureFactor(w);if(scalarValue(f)==0)return k;
        const auto v=volumes(I,w);const Real vw=18.016/w.rho;
        k[OH]-=f*(v[OH]-vw);
        k[Bicarbonate]-=f*(v[Bicarbonate]-v[Carbonate]);
        k[CO2]-=f*(v[CO2]+vw-v[Carbonate]);
        k[Bisulfate]-=f*(v[Bisulfate]-v[Sulfate]);
        k[MgOH]-=f*(v[MgOH]-v[MgIon]-vw);
        k[MgCO3]-=f*(v[MgCO3]-v[MgIon]-v[Carbonate]);return k;
    }
    template<class Real> std::array<Real,nMinerals> phaseLogK(Real I,Real aw)const {
        std::array<Real,nMinerals> k;for(int j=0;j<nMinerals;j++)k[j]=mineralLogK[j];if(pressure<=1&&water.pSat<=1)return k;
        const auto w=water.at(aw);const Real f=pressureFactor(w);if(scalarValue(f)==0)return k;
        const auto v=volumes(I,w);const Real vw=18.016/w.rho;
        for(int j=0;j<nMinerals;j++){Real dv=hydration[j]*vw-molarVolume[j]*1e6;for(int i=0;i<nComponents;i++)dv+=stoichiometry[j][i]*v[i];k[j]-=f*dv;}return k;
    }
    std::size_t interactionCount()const{return terms.size();}
    template<class Real> void activities(const std::array<Real,nSpecies>& m,std::array<Real,nSpecies>& lg,Real& lnAw,Real& I,Real solventActivity=Real(1))const {
        lg.fill(Real(0));I=0;Real sum=0,Z=0;
        for(int i=0;i<nSpecies;i++){if(scalarValue(m[i])<0||!std::isfinite(scalarValue(m[i])))throw std::domain_error("Invalid species molality");I+=.5*charge[i]*charge[i]*m[i];Z+=std::abs(charge[i])*m[i];sum+=m[i];}
        if(scalarValue(I)<1e-30)I=Real(1e-30);const Real s=sqrt(I);
        Real gb[128],gp[128],ex[128];for(size_t k=0;k<activityAlphas.size();k++){Real x=activityAlphas[k]*s;gb[k]=G(x);gp[k]=GP(x);ex[k]=exp(-x);}
        const auto solvent=water.at(solventActivity);const Real aphi=solvent.aphi;
        const Real Fdh=-aphi*(s/(1+1.2*s)+2*log1p(1.2*s)/1.2);Real F=Fdh,cs=0,osm=-aphi*I*s/(1+1.2*s);
        for(const auto& t:terms){if(t.value==0)continue;const int i=t.i,j=t.j,k=t.k;const double p=t.value;const Real ij=m[i]*m[j];
            switch(t.kind){
            case B0:lg[i]+=2*m[j]*p;lg[j]+=2*m[i]*p;osm+=ij*p;break;
            case B1:case B2:{Real b=p*gb[t.group];F+=ij*p*gp[t.group]/I;lg[i]+=2*m[j]*b;lg[j]+=2*m[i]*b;osm+=ij*p*ex[t.group];break;}
            case C0:{double c=t.normalized;cs+=ij*c;lg[i]+=m[j]*Z*c;lg[j]+=m[i]*Z*c;osm+=ij*Z*c;break;}
            case Theta:lg[i]+=2*m[j]*p;lg[j]+=2*m[i]*p;osm+=ij*p;break;
            case Lambda:{double a=i==j?1:2;lg[i]+=a*m[j]*p;lg[j]+=a*m[i]*p;osm+=ij*p*(i==j?.5:1);break;}
            case Psi:case Zeta:lg[i]+=m[j]*m[k]*p;lg[j]+=m[i]*m[k]*p;lg[k]+=ij*p;osm+=ij*m[k]*p;break;
            }
        }
        // All supported charged species have |z| = 1 or 2. Compute E-theta once.
        if(scalarValue(I)>1e-14){Real x=6*aphi*s;auto a=jay(x),b=jay(2*x),c=jay(4*x);Real e=(b[0]-.5*a[0]-.5*c[0])/(2*I),ep=(b[1]-.5*a[1]-.5*c[1])/(4*I*I)-e/I;
            for(auto pair:differentChargePairs){int i=pair[0],j=pair[1];F+=m[i]*m[j]*ep;lg[i]+=2*m[j]*e;lg[j]+=2*m[i]*e;osm+=m[i]*m[j]*(e+I*ep);}
        }
        Real f1=Fdh,f2=Fdh;
        if(scalarValue(solvent.pressure)>1){
            Real pap=(7e-5+1.93e-9*std::pow(temperature-250.,2))*solvent.pressure;
            const Real b1=1.2-(scalarValue(pap)<.2?pap:Real(.2));
            if(temperature>263)pap=9.65e-10*std::pow(temperature-263.,2.773)*pow(solvent.pressure,.623);
            const Real b2=1.2-(scalarValue(pap)<.2?pap:Real(.2));
            f1=-aphi*(s/(1+b1*s)+2*log1p(b1*s)/b1);
            f2=-aphi*(s/(1+b2*s)+2*log1p(b2*s)/b2);
        }
        for(int i=0;i<nSpecies;i++)lg[i]+=charge[i]*charge[i]*(F+(std::abs(charge[i])==1?f1:std::abs(charge[i])==2?f2:Fdh)-Fdh)+std::abs(charge[i])*cs;
        Real x=2*s;Real gx=scalarValue(x)<1e-3?1-5*x/6+3*x*x/8: (1-(1+x-x*x/2)*exp(-x))/(x*x);
        const Real mac=lg[ClIon]-(f1+2*I*(macB0+macB1*gx)+1.5*macC0*I*I);
        for(int i=0;i<nSpecies;i++)lg[i]=(lg[i]+charge[i]*mac)/ln10;
        lnAw=-(sum+2*osm)/55.50837;
    }
    double density(const SpeciesVector& m,double I,double waterActivity=1)const{
        const auto w=water.at(waterActivity);const auto v=volumes(I,w);double volume=1/w.rho,mass=1;
        for(int i=0;i<nSpecies;i++){volume+=m[i]*v[i]/1000;mass+=m[i]*speciesMass[i];}
        if(volume<=0||!std::isfinite(volume))throw std::domain_error("Invalid solution volume");return 1000*mass/volume;
    }
    void finish(AqueousState& s,double aw,bool computeDensity=true)const{
        if(strictWater&&waterState(checkedExp(aw)).pressureRaised)throw std::domain_error("PHREEQC effective water pressure exceeds configured pressure; use compatible water policy or raise pressure");
        s.waterActivity=checkedExp(aw);if(s.waterActivity<=0)throw std::domain_error("Water activity below floating-point range");if(computeDensity)s.density=density(s.molality,s.ionicStrength,s.waterActivity);s.chargeResidual=chargeBalance(s.molality);
        for(int i=0;i<nSpecies;i++)s.logActivity[i]=std::log10(std::max(s.molality[i],1e-300))+s.logGamma[i];s.pH=-s.logActivity[H];
        const auto k=phaseLogK(s.ionicStrength,s.waterActivity);
        for(int j=0;j<nMinerals;j++){
            double si=-k[j]+hydration[j]*aw/ln10;bool present=true;
            for(int i=0;i<nComponents;i++)if(stoichiometry[j][i]){si+=stoichiometry[j][i]*s.logActivity[i];if(s.molality[i]==0)present=false;}
            // Match PHREEQC's SI sentinel when a required aqueous component
            // is absent; a floor inside log(activity) is not a physical SI.
            s.saturationIndex[j]=present?si:-99.99;
        }
    }
    void complete(AqueousState& s)const{double aw;activities(s.molality,s.logGamma,aw,s.ionicStrength,s.waterActivity);finish(s,aw);}
};
}
