// SPDX-License-Identifier: GPL-3.0-or-later
#include "BrineReactions.hpp"
#include "IPhreeqc.hpp"
#include "PhreeqcBackend.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <filesystem>
using namespace saltchem;
void run(IPhreeqc& p,const std::string& s){if(p.RunString(s.c_str()))throw std::runtime_error(p.GetErrorString());}
double value(IPhreeqc& p,int col){VAR v;VarInit(&v);if(p.GetSelectedOutputValue(p.GetSelectedOutputRowCount()-1,col,&v))throw std::runtime_error("Missing reference output");double x=v.type==TT_DOUBLE?v.dVal:v.lVal;VarClear(&v);return x;}
int main(int argc,char** argv){try{
    const std::string db=argc>1?argv[1]:"src/chemistry/pitzer.dat";
    const std::filesystem::path destination=argc>2?argv[2]:"results/checks/thermochemistry";std::filesystem::create_directories(destination);
    IPhreeqc ref;ref.SetOutputFileOn(false);ref.SetLogFileOn(false);ref.SetSelectedOutputFileOn(false);if(ref.LoadDatabase(db.c_str()))throw std::runtime_error(ref.GetErrorString());
    std::ostringstream setup;setup<<"SELECTED_OUTPUT 1\n-reset false\n-high_precision true\nUSER_PUNCH 1\n-start\n10 PUNCH -LA(\"H+\"), ACT(\"H2O\"), RHO, MU\n20 PUNCH ";
    for(int j=0;j<7;j++)setup<<(j?", ":"")<<"SI(\""<<mineralNames[j]<<"\")";
    setup<<"\n30 PUNCH ";for(int j=0;j<14;j++)setup<<(j?", ":"")<<"MOL(\""<<speciesNames[j]<<"\"), LA(\""<<speciesNames[j]<<"\"), LG(\""<<speciesNames[j]<<"\")";
    setup<<"\n-end\nEND\n";run(ref,setup.str());
    PitzerModel model(db);Speciation spec(model);int cases=0;double worstSI=0,worstRho=0,worstSpecies=0,worstGamma=0,worstAw=0,worstPH=0;
    std::ofstream csv(destination/"thermodynamics.csv");csv<<std::setprecision(16)<<"temperature_C,pressure_atm,brine,max_SI_error,max_log_activity_error,max_log_gamma_error,density_relative_error,aw_relative_error,pH_error\n";
    const Totals states[]={{{6,0,0,0,6,0,0}},{{4,.2,.1,.5,5.328,.035,.002}},{{.2,.02,.01,.03,.30,.001,.002}},{{1,.1,.03,.1,1.4,.002,.001}}};
    for(auto tp:std::initializer_list<std::array<double,2>>{{5,1},{25,1},{25,50},{40,10},{60,100},{90,100},{10,.8},{120,50},{150,500}}){
        model.setConditions(tp[0]+273.15,tp[1]);
        for(int b=0;b<4;b++){
            auto s=spec.solve(states[b]);std::ostringstream in;in<<std::setprecision(17)<<"SOLUTION 1\n-units mol/kgw\n-temp "<<tp[0]<<"\n-pressure "<<tp[1]<<"\n-pH 7 charge\n-water 1\n";
            for(int j=0;j<7;j++)in<<componentNames[j]<<' '<<states[b][j]<<'\n';in<<"END\n";run(ref,in.str());
            double si=0,act=0,gamma=0;
            for(int j=0;j<7;j++)si=std::max(si,std::abs(s.saturationIndex[j]-value(ref,4+j)));
            for(int j=0;j<14;j++)if(value(ref,11+3*j)>1e-25){act=std::max(act,std::abs(s.logActivity[j]-value(ref,12+3*j)));gamma=std::max(gamma,std::abs(s.logGamma[j]-value(ref,13+3*j)));}
            const double rho=std::abs(s.density/(1000*value(ref,2))-1),aw=std::abs(s.waterActivity/value(ref,1)-1),ph=std::abs(s.pH-value(ref,0));
            csv<<tp[0]<<','<<tp[1]<<','<<b<<','<<si<<','<<act<<','<<gamma<<','<<rho<<','<<aw<<','<<ph<<'\n';
            worstSI=std::max(worstSI,si);worstSpecies=std::max(worstSpecies,act);worstGamma=std::max(worstGamma,gamma);worstRho=std::max(worstRho,rho);worstAw=std::max(worstAw,aw);worstPH=std::max(worstPH,ph);cases++;
            if(si>1e-6||act>1e-6||gamma>1e-6||rho>1e-7||aw>1e-7||ph>1e-6)std::cout<<std::setprecision(10)<<"MISMATCH T="<<tp[0]<<" P="<<tp[1]<<" brine="<<b<<" SI="<<si<<" loga="<<act<<" gamma="<<gamma<<" rho="<<rho<<" aw="<<aw<<" pH="<<ph<<'\n';
        }
    }
    std::cout<<std::setprecision(10)<<"TP cases="<<cases<<" SI="<<worstSI<<" loga="<<worstSpecies<<" gamma="<<worstGamma<<" rho="<<worstRho<<" aw="<<worstAw<<" pH="<<worstPH<<std::endl;
    bool failed=worstSI>1e-6||worstSpecies>1e-6||worstGamma>1e-6||worstRho>1e-7||worstAw>1e-7||worstPH>1e-6;
    std::ofstream kinetic(destination/"kinetics.csv");
    kinetic<<std::setprecision(16)<<"step,temperature_C,pressure_atm,state_error,water_error,density_relative_error,native_seconds,reference_seconds\n";
    model.setConditions(298.15,1);Reactions reaction(model);MineralVector minerals;minerals.fill(.2);
    std::vector<Inventory> native(4);std::vector<RateParameters> rates(4);
    for(int b=0;b<4;b++){
        native[b]=spec.inventory(states[b],1,minerals);
        rates[b].neutral={{.02,.002,.02,1e-5,1e-4,.02,1e-5}};
        rates[b].area.fill(.1);rates[b].referenceMoles=minerals;
    }
    auto reference=native;SevenPhreeqcBackend rm(db,reference,1,false,1e-12,25,1);
    IntegrationOptions options;
    int step=0;double maxKinetic=0,maxWater=0,maxDensity=0,totalNative=0,totalReference=0;
    using Clock=std::chrono::steady_clock;
    for(auto tp:std::initializer_list<std::array<double,2>>{{25,1},{40,10},{60,100},{10,50},{25,1}}){
        model.setConditions(tp[0]+273.15,tp[1]);rm.setConditions(tp[0],tp[1]);
        auto start=Clock::now();for(int b=0;b<4;b++)reaction.advance(native[b],.5,rates[b],options);
        double nt=std::chrono::duration<double>(Clock::now()-start).count();start=Clock::now();rm.advance(reference,rates,.5);
        double rt=std::chrono::duration<double>(Clock::now()-start).count(),err=0,we=0,re=0;
        for(int b=0;b<4;b++){
            for(int j=0;j<7;j++){err=std::max(err,std::abs(native[b].mineral[j]-reference[b].mineral[j]));err=std::max(err,std::abs(native[b].aqueous[j]-reference[b].aqueous[j]));}
            we=std::max(we,std::abs(native[b].waterEquivalent-reference[b].waterEquivalent));
            double w1,w2;const auto a=spec.conserved(native[b],w1),c=spec.conserved(reference[b],w2);re=std::max(re,std::abs(a.density/c.density-1));
        }
        maxKinetic=std::max(maxKinetic,err);maxWater=std::max(maxWater,we);maxDensity=std::max(maxDensity,re);totalNative+=nt;totalReference+=rt;
        kinetic<<step++<<','<<tp[0]<<','<<tp[1]<<','<<err<<','<<we<<','<<re<<','<<nt<<','<<rt<<'\n';
    }
    std::ofstream metrics(destination/"metrics.csv");
    metrics<<std::setprecision(16)<<"metric,value\nthermodynamic_states,"<<cases<<"\nmax_SI_error,"<<worstSI<<"\nmax_log_activity_error,"<<worstSpecies<<"\nmax_log_gamma_error,"<<worstGamma<<"\nmax_density_relative_error,"<<worstRho<<"\nkinetic_cells,4\nkinetic_steps,"<<step<<"\nkinetic_state_error_mol,"<<maxKinetic<<"\nkinetic_water_error_kg,"<<maxWater<<"\nkinetic_density_relative_error,"<<maxDensity<<"\nnative_kinetics_seconds,"<<totalNative<<"\nreference_kinetics_seconds,"<<totalReference<<'\n';
    std::cout<<"TP KINETICS state="<<maxKinetic<<" water="<<maxWater<<" rho="<<maxDensity<<" native_s="<<totalNative<<" reference_s="<<totalReference<<std::endl;
    {
        MineralVector mixedMinerals{};mixedMinerals[1]=.5;
        std::vector<Inventory> exposed(1,spec.inventory({{6,.001,0,0,6.001,0,0}},1000,mixedMinerals));
        auto compared=exposed;std::vector<RateParameters> surface(1);auto& p=surface[0];
        p.geometry.enabled=true;p.geometry.separateMineralAreas=true;p.geometry.lockedMoles[1]=.4;
        p.geometry.surfaceExponent[1]=2./3;p.referenceMoles[1]=.1;p.area[1]=p.geometry.growthArea[1]=1e4;p.neutral[1]=.02;
        SevenPhreeqcBackend interfaceReference(db,compared,1,true,1e-14);
        for(int iteration=0;iteration<2;iteration++){
            reaction.advance(exposed[0],.01,p,options);interfaceReference.advance(compared,surface,.01);
            if(exposed[0].mineral[1]<.4||compared[0].mineral[1]<.4||exposed[0].mineral[1]>.40000001||compared[0].mineral[1]>.40000001)
                throw std::runtime_error("Reference/native solver consumed buried mineral or failed to exhaust exposed inventory");
            for(int j=0;j<7;j++)if(std::abs(exposed[0].aqueous[j]-compared[0].aqueous[j])>1e-7)
                throw std::runtime_error("Reference/native locked-inventory aqueous balance mismatch");
            if(std::abs(exposed[0].waterEquivalent-compared[0].waterEquivalent)>1e-8)
                throw std::runtime_error("Reference/native locked-inventory water balance mismatch");
        }
        std::cout<<"REFERENCE_LOCKED_INVENTORY native="<<exposed[0].mineral[1]<<" reference="<<compared[0].mineral[1]<<std::endl;
    }
    return failed||maxKinetic>1e-7||maxWater>1e-8||maxDensity>1e-7;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 2;}}
