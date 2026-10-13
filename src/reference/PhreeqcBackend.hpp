// SPDX-License-Identifier: GPL-3.0-or-later
// Reference backend only. The production chemistry kernel does not link it.
// PhreeqcRM receives conserved inventories through its array API; rate
// parameters and geometric areas use BASIC callbacks. Rebase amounts to the
// current water inventory before each call, keeping the physical system fixed.
#pragma once
#include "BrineReactions.hpp"
#include "PhreeqcRM.h"
#include "IPhreeqc.hpp"
#include <memory>
#include <iomanip>
#include <cstring>
#include <atomic>
namespace saltchem {
class SevenPhreeqcBackend {
    PhreeqcRM rm;
    std::vector<double> scales;
    std::vector<double> originalScales;
    std::vector<RateParameters> current;
    std::vector<std::string> components;
    const bool phaseGeometry;
    MineralVector molarVolume;
    struct RateCache {Inventory state;MineralVector si{},rates{};double logH=0,logCO2=0;bool dirty=true;};
    std::vector<RateCache> rateCache;
    std::vector<std::string> callbackErrors;
    double time=0;
    std::atomic<bool> invalidCallback{false};
    static double callback(double cell,double phase,const char* field,void* ptr){
        auto& self=*static_cast<SevenPhreeqcBackend*>(ptr);const int i=int(std::llround(cell));
        if(i<0||i>=int(self.current.size())){self.invalidCallback=true;return !std::strcmp(field,"reference")?1:0;}
        const auto& p=self.current[i];
        if(self.phaseGeometry){
            auto& c=self.rateCache[i];
            auto set=[&](double& destination,double value){if(destination!=value){destination=value;c.dirty=true;}return 0.;};
            if(field[0]=='M'&&field[1]>='0'&&field[1]<='6')return set(c.state.mineral[field[1]-'0'],phase*self.scales[i]);
            if(field[0]=='S'&&field[1]>='0'&&field[1]<='6')return set(c.si[field[1]-'0'],phase);
            if(!std::strcmp(field,"logH"))return set(c.logH,phase);
            if(!std::strcmp(field,"logCO2"))return set(c.logCO2,phase);
            if(!std::strcmp(field,"rate")){
                const int j=int(std::llround(phase));if(j<0||j>=nMinerals){self.invalidCallback=true;return 0;}
                try{if(c.dirty){c.rates=surfaceRates(c.state,p,c.si,c.logH,c.logCO2,self.molarVolume);c.dirty=false;}return c.rates[j]/self.scales[i];}
                catch(const std::exception& e){self.callbackErrors[i]=e.what();self.invalidCallback=true;return 0;}
            }
        }
        const int j=int(std::llround(phase));if(j<0||j>=nMinerals){self.invalidCallback=true;return 0;}
        if(!std::strcmp(field,"area"))return p.area[j]/self.scales[i];
        if(!std::strcmp(field,"reference"))return p.referenceMoles[j]/self.scales[i];
        if(!std::strcmp(field,"neutral"))return p.neutral[j];
        if(!std::strcmp(field,"acid"))return p.acid[j];
        if(!std::strcmp(field,"co2"))return p.co2[j];
        if(!std::strcmp(field,"acidOrder"))return p.acidOrder[j];
        if(!std::strcmp(field,"surfaceExponent"))return p.dissolutionExponent[j];
        if(!std::strcmp(field,"affinityExponent"))return p.affinityExponent[j];
        self.invalidCallback=true;return 0;
    }
    std::string errors() {
        std::string message=rm.GetErrorString();
        for(int worker=0;worker<rm.GetThreadCount();worker++)
            message+=rm.GetIPhreeqcPointer(worker)->GetErrorString();
        return message;
    }
    void check(int status){if(status<0)throw std::runtime_error(errors());}
    void rebase(const std::vector<Inventory>& state,const std::vector<RateParameters>& rates){
        for(size_t i=0;i<state.size();i++){
            if(!(state[i].waterEquivalent>0)||!std::isfinite(state[i].waterEquivalent))
                throw std::invalid_argument("Nonpositive reference normalization water");
            scales[i]=state[i].waterEquivalent;
        }
        std::vector<std::string> failures(rm.GetThreadCount());
        const auto& first=rm.GetStartCell();const auto& last=rm.GetEndCell();
        #pragma omp parallel for num_threads(rm.GetThreadCount()) schedule(static)
        for(int worker=0;worker<rm.GetThreadCount();worker++){
            std::ostringstream input;input<<std::setprecision(17);
            for(int i=first[worker];i<=last[worker];i++){
                input<<"KINETICS_MODIFY "<<i<<'\n';
                for(int j=0;j<nMinerals;j++)input<<"-component "<<mineralNames[j]
                    <<"\n-m "<<(state[i].mineral[j]-rates[i].geometry.lockedMoles[j])/scales[i]
                    <<"\n-m0 "<<originalScales[i]/scales[i]<<'\n';
            }
            input<<"USE solution none\nUSE kinetics none\nEND\n";
            auto* interpreter=rm.GetIPhreeqcPointer(worker);
            if(interpreter->RunString(input.str().c_str())!=0)failures[worker]=interpreter->GetErrorString();
        }
        for(const auto& failure:failures)if(!failure.empty())throw std::runtime_error(failure);
    }
    std::string definitions()const{
        std::ostringstream s;s<<"RATES\n";
        for(int j=0;j<nMinerals;j++){
            if(phaseGeometry){
                s<<mineralNames[j]<<"\n-start\n";int line=10;
                for(int k=0;k<nMinerals;k++)s<<line++<<" z = CALLBACK(PARM(1),KIN(\""<<mineralNames[k]<<"\"),\"M"<<k<<"\")\n";
                for(int k=0;k<nMinerals;k++)s<<line++<<" z = CALLBACK(PARM(1),SI(\""<<mineralNames[k]<<"\"),\"S"<<k<<"\")\n";
                s<<line++<<" z = CALLBACK(PARM(1),LA(\"H+\"),\"logH\")\n";
                s<<line++<<" z = CALLBACK(PARM(1),LA(\"CO2\"),\"logCO2\")\n";
                s<<line++<<" SAVE CALLBACK(PARM(1),"<<j<<",\"rate\")*TIME\n-end\n";continue;
            }
            s<<mineralNames[j]<<"\n-start\n10 a = 1 - SR(\""<<mineralNames[j]<<"\")\n"
             <<"20 area = CALLBACK(PARM(1),"<<j<<",\"area\")\n"
             <<"30 r = CALLBACK(PARM(1),"<<j<<",\"neutral\") + CALLBACK(PARM(1),"<<j<<",\"acid\")*ACT(\"H+\")^CALLBACK(PARM(1),"<<j<<",\"acidOrder\") + CALLBACK(PARM(1),"<<j<<",\"co2\")*ACT(\"CO2\")\n"
             <<"40 IF a > 0 AND M <= 0 THEN r = 0 : GOTO 70\n"
             <<"50 IF a > 0 THEN area = area*(M/CALLBACK(PARM(1),"<<j<<",\"reference\"))^CALLBACK(PARM(1),"<<j<<",\"surfaceExponent\")\n"
             <<"60 r = r*area*ABS(a)^CALLBACK(PARM(1),"<<j<<",\"affinityExponent\")\n"
             <<"70 IF a < 0 THEN r = -r\n80 SAVE r*TIME\n-end\n";
        }
        s<<"SELECTED_OUTPUT 1\n-reset false\n-high_precision true\nUSER_PUNCH 1\n-headings H S C P G B Ca Na K Ca_aq Mg Cl S_aq C_aq water_equivalent rho pH water\n-start\n10 PUNCH ";
        for(int j=0;j<nMinerals;j++)s<<(j?", ":"")<<"KIN(\""<<mineralNames[j]<<"\")";
        s<<"\n20 PUNCH ";for(int j=0;j<nComponents;j++)s<<(j?", ":"")<<"TOT(\""<<componentNames[j]<<"\")*TOT(\"water\")";
        s<<"\n30 PUNCH TOT(\"water\")*(1+0.018016*(MOL(\"OH-\")+MOL(\"MgOH+\")-MOL(\"CO2\"))), RHO, -LA(\"H+\"), TOT(\"water\")\n-end\nEND\n";return s.str();
    }
public:
    SevenPhreeqcBackend(const std::string& db,const std::vector<Inventory>& initial,int threads=1,bool geometry=false,double kineticTolerance=1e-11,double temperatureC=25,double pressureAtm=1):rm(int(initial.size()),threads),scales(initial.size()),current(initial.size()),phaseGeometry(geometry),molarVolume(PitzerModel(db).molarVolume),rateCache(initial.size()),callbackErrors(initial.size()){
        if(initial.empty())throw std::invalid_argument("Empty reference batch");
        if(!(kineticTolerance>0)||!std::isfinite(kineticTolerance))throw std::invalid_argument("Invalid reference kinetic tolerance");
        rm.SetErrorHandlerMode(1);rm.SetScreenOn(false);rm.SetPrintChemistryOn(false,false,false);rm.SetSelectedOutputOn(true);rm.SetUnitsSolution(2);rm.SetUnitsKinetics(0);rm.SetComponentH2O(true);rm.UseSolutionDensityVolume(true);
        const int n=int(initial.size());std::vector<double> one(n,1);rm.SetRepresentativeVolume(one);rm.SetPorosity(one);rm.SetSaturationUser(one);
        check(rm.LoadDatabase(db));
        // Loading a database recreates the internal PHREEQC interpreter.
        for(int worker=0;worker<rm.GetThreadCount()+2;worker++)rm.GetIPhreeqcPointer(worker)->SetBasicCallback(callback,this);
        check(rm.RunString(true,true,true,definitions()));
        std::ostringstream s;s<<std::setprecision(17)<<"SOLUTION 1\n-units mol/kgw\n-temp "<<temperatureC<<"\n-pressure "<<pressureAtm<<"\n-pH 7 charge\n-water 1\nNa 1\nK .1\nCa .01\nMg .1\nCl 1.3\nS(6) .01\nC(4) .001\nEND\n";
        std::vector<int> mapping(n*7,-1);
        for(int i=0;i<n;i++){
            scales[i]=initial[i].waterEquivalent;if(scales[i]<=0)throw std::invalid_argument("Reference cell needs positive water");
            s<<"KINETICS "<<i+1<<'\n';
            // CELL_NO may identify PHREEQC's temporary CVODE solution; keep
            // the transport index explicitly with each kinetic reactant.
            for(int j=0;j<nMinerals;j++)s<<mineralNames[j]<<"\n-formula "<<mineralNames[j]<<" 1\n-m "<<initial[i].mineral[j]/scales[i]<<"\n-m0 1\n-tol "<<kineticTolerance<<"\n-parms "<<i<<'\n';
            s<<"-cvode true\n-cvode_steps 10000\n-steps 0\nUSE solution none\nEND\n";mapping[i]=1;mapping[6*n+i]=i+1;
        }
        check(rm.RunString(false,true,false,s.str()));check(rm.FindComponents());components=rm.GetComponents();check(rm.InitialPhreeqc2Module(mapping));rm.SetCurrentSelectedOutputUserNumber(1);setConditions(temperatureC,pressureAtm);originalScales=scales;
    }
    void setConditions(double temperatureC,double pressureAtm){
        if(!std::isfinite(temperatureC)||!std::isfinite(pressureAtm)||pressureAtm<=0)throw std::invalid_argument("Invalid reference temperature/absolute pressure");
        check(rm.SetTemperature(std::vector<double>(scales.size(),temperatureC)));
        check(rm.SetPressure(std::vector<double>(scales.size(),pressureAtm)));
    }
    void setPressures(const std::vector<double>& pressureAtm){
        if(pressureAtm.size()!=scales.size())throw std::invalid_argument("Reference pressure field size changed");
        for(double p:pressureAtm)if(!std::isfinite(p)||p<=0)throw std::invalid_argument("Nonpositive reference absolute pressure");
        check(rm.SetPressure(pressureAtm));
    }
    void setTime(double value){time=value;}
    void advance(std::vector<Inventory>& state,const std::vector<RateParameters>& rates,double dt,std::vector<AqueousState>* properties=nullptr){
        const int n=int(scales.size());if(state.size()!=scales.size()||rates.size()!=scales.size())throw std::invalid_argument("Reference batch size changed");
        current=rates;
        for(int i=0;i<n;i++){
            Reactions::validateParameters(rates[i]);
            for(int j=0;j<nMinerals;j++){
                const double locked=rates[i].geometry.lockedMoles[j];
                if(locked>state[i].mineral[j])throw std::invalid_argument("Reference inaccessible mineral exceeds inventory");
                if(locked>0){
                    if(!phaseGeometry)throw std::invalid_argument("Reference inaccessible minerals require geometry callbacks");
                    current[i].geometry.inertFraction+=locked*molarVolume[j]/rates[i].geometry.bulkVolume;
                }
                current[i].geometry.lockedMoles[j]=0;
            }
        }
        // PHREEQC's kinetic positivity bound must be the same exposed
        // inventory as in the native solver, rather than the buried reserve.
        rebase(state,rates);
        for(int i=0;i<n;i++){rateCache[i].dirty=true;callbackErrors[i].clear();}std::vector<double> concentrations(components.size()*n);
        for(int i=0;i<n;i++){
            Totals a=state[i].aqueous;for(double& v:a)v/=scales[i];const double water=state[i].waterEquivalent/(waterMolarMass*scales[i]);
            for(size_t k=0;k<components.size();k++){
                const auto& name=components[k];double value=0;
                if(name=="H2O")value=water;else if(name=="H")value=-strongCharge(a);else if(name=="O")value=4*a[S]+3*a[C];else if(name=="Charge")value=0;
                else {bool found=false;for(int j=0;j<nComponents;j++)if(name==componentNames[j]||(j==S&&name=="S")||(j==C&&name=="C")){value=a[j];found=true;break;}if(!found)throw std::runtime_error("Unexpected reference component: "+name);}
                concentrations[k*n+i]=value;
            }
        }
        check(rm.SetConcentrations(concentrations));rm.SetTimeStep(dt);time+=dt;rm.SetTime(time);invalidCallback=false;
        try{check(rm.RunCells());}catch(const std::exception& e){throw std::runtime_error(std::string(e.what())+"\n"+errors());}
        if(invalidCallback){for(int i=0;i<n;i++)if(!callbackErrors[i].empty())throw std::runtime_error("Reference cell "+std::to_string(i)+": "+callbackErrors[i]);throw std::runtime_error("Invalid PhreeqcRM rate callback index or field");}
        std::vector<double> out;check(rm.GetSelectedOutput(out));if(out.size()!=size_t(n*18))throw std::runtime_error("Unexpected reference output width");
        if(properties)properties->resize(n);
        for(int i=0;i<n;i++){
            for(int j=0;j<nMinerals;j++)state[i].mineral[j]=out[j*n+i]*scales[i]+rates[i].geometry.lockedMoles[j];
            for(int j=0;j<nComponents;j++)state[i].aqueous[j]=out[(7+j)*n+i]*scales[i];
            state[i].waterEquivalent=out[14*n+i]*scales[i];
            if(properties){(*properties)[i].density=out[15*n+i]*1000;(*properties)[i].pH=out[16*n+i];}
        }
    }
};
}
