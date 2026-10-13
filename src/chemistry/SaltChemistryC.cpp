// SPDX-License-Identifier: GPL-3.0-or-later
#include "SaltChemistryC.h"
#include "BrineReactions.hpp"
#include <omp.h>
#include <memory>
#include <cstring>
using namespace saltchem;
struct SaltChemistryHandle {PitzerModel model;Speciation spec;Reactions reaction;SaltChemistryHandle(const char* path,double t,double p):model(path,t,p),spec(model),reaction(model){}};
static thread_local std::string lastError;
static Inventory inventory(const SaltCell& c){Inventory v;std::copy_n(c.aqueous_mol,7,v.aqueous.begin());std::copy_n(c.mineral_mol,7,v.mineral.begin());v.waterEquivalent=c.water_equivalent_kg;return v;}
static SaltCell cell(const Inventory& v){SaltCell c{};std::copy(v.aqueous.begin(),v.aqueous.end(),c.aqueous_mol);std::copy(v.mineral.begin(),v.mineral.end(),c.mineral_mol);c.water_equivalent_kg=v.waterEquivalent;return c;}
static RateParameters parameters(const SaltRates& c){RateParameters p;std::copy_n(c.k_neutral,7,p.neutral.begin());std::copy_n(c.k_acid,7,p.acid.begin());std::copy_n(c.k_co2,7,p.co2.begin());std::copy_n(c.acid_order,7,p.acidOrder.begin());std::copy_n(c.area_m2,7,p.area.begin());std::copy_n(c.reference_mol,7,p.referenceMoles.begin());std::copy_n(c.dissolution_exponent,7,p.dissolutionExponent.begin());std::copy_n(c.affinity_exponent,7,p.affinityExponent.begin());return p;}
static SaltResult properties(const SaltChemistryHandle& h,const Inventory& v){
    SaltResult r{};auto s=h.spec.conserved(v,r.free_water_kg);r.density_kg_m3=s.density;r.pH=s.pH;r.water_activity=s.waterActivity;r.aqueous_mass_kg=aqueousMass(v);
    r.aqueous_volume_m3=r.aqueous_mass_kg/r.density_kg_m3;
    std::copy(s.saturationIndex.begin(),s.saturationIndex.end(),r.saturation_index);return r;
}
extern "C" int salt_chemistry_create(const char* database,double t,double p,SaltChemistryHandle** out){
    lastError.clear();if(!out){lastError="Null output handle";return 1;}*out=nullptr;
    try{if(!database)throw std::invalid_argument("Null database path");*out=new SaltChemistryHandle(database,t,p);return 0;}catch(const std::exception& e){lastError=e.what();return 1;}
}
extern "C" void salt_chemistry_destroy(SaltChemistryHandle* h){delete h;}
extern "C" int salt_chemistry_set_conditions(SaltChemistryHandle* h,double t,double p){
    lastError.clear();try{
        if(!h)throw std::invalid_argument("Null chemistry handle");
        if(t==h->model.temperatureK()&&p==h->model.pressureAtm())return 0;
        auto next=h->model;next.setConditions(t,p);h->model=std::move(next);return 0;
    }catch(const std::exception& e){lastError=e.what();return 1;}
}
extern "C" const char* salt_chemistry_last_error(){return lastError.c_str();}
extern "C" int salt_chemistry_set_water_policy(SaltChemistryHandle* h,int strict){
    lastError.clear();try{if(!h||(strict!=0&&strict!=1))throw std::invalid_argument("Invalid water policy");h->model.setStrictWater(strict!=0);return 0;}catch(const std::exception& e){lastError=e.what();return 1;}
}
extern "C" int salt_chemistry_water_state(const SaltChemistryHandle* h,double aw,SaltWaterState* result){
    lastError.clear();try{if(!h||!result||!(aw>0)||!std::isfinite(aw))throw std::invalid_argument("Invalid water diagnostic argument");const auto s=h->model.waterState(aw);*result={s.requestedTemperatureK,s.waterTemperatureK,s.requestedPressureAtm,s.effectivePressureAtm,s.temperatureCapped,s.pressureRaised};return 0;}catch(const std::exception& e){lastError=e.what();return 1;}
}
extern "C" int salt_chemistry_initialize(const SaltChemistryHandle* h,const double t[7],double w,const double mineral[7],SaltCell* out,SaltResult* result){
    lastError.clear();try{
        if(!h||!t||!mineral||!out||!result)throw std::invalid_argument("Null initialization argument");Totals totals;MineralVector minerals;std::copy_n(t,7,totals.begin());std::copy_n(mineral,7,minerals.begin());
        for(double n:minerals)if(n<0||!std::isfinite(n))throw std::invalid_argument("Invalid mineral inventory");
        auto v=h->spec.inventory(totals,w,minerals);auto r=properties(*h,v);*out=cell(v);*result=r;return 0;
    }catch(const std::exception& e){lastError=e.what();return 1;}
}
extern "C" int salt_chemistry_advance(const SaltChemistryHandle* h,size_t count,SaltCell* cells,const SaltRates* rates,double dt,const SaltOptions* opt,SaltResult* results,int* status){
    lastError.clear();try{
        if(!h||!opt)throw std::invalid_argument("Null chemistry handle or batch options");if(dt<0||!std::isfinite(dt)||opt->threads<1)throw std::invalid_argument("Invalid interval or thread count");
        if(!std::isfinite(opt->relative_tolerance)||!std::isfinite(opt->absolute_tolerance_mol_per_kg)||opt->relative_tolerance<=0||opt->absolute_tolerance_mol_per_kg<=0||opt->maximum_attempts<=0)throw std::invalid_argument("Invalid integration tolerance or attempt limit");
        if(count==0)return 0; // An MPI rank may own no wet cells.
        if(!cells||!rates||!results||!status)throw std::invalid_argument("Null batch array");
        IntegrationOptions options;options.relativeTolerance=opt->relative_tolerance;options.absoluteTolerance=opt->absolute_tolerance_mol_per_kg;options.maximumAttempts=opt->maximum_attempts;
        std::vector<SaltCell> next(count);std::vector<SaltResult> values(count);std::vector<std::string> errors(count);
        #pragma omp parallel for schedule(dynamic,8) num_threads(opt->threads)
        for(size_t i=0;i<count;i++){
            status[i]=0;
            try{
                auto state=inventory(cells[i]);const auto before=state;const double oldMass=aqueousMass(state);auto stats=h->reaction.advance(state,dt,parameters(rates[i]),options);auto r=properties(*h,state);
                if(dt>0){for(int j=0;j<7;j++){r.component_source_mol_s[j]=(state.aqueous[j]-before.aqueous[j])/dt;r.mineral_change_mol[j]=state.mineral[j]-before.mineral[j];}r.aqueous_mass_source_kg_s=(r.aqueous_mass_kg-oldMass)/dt;}
                r.accepted_steps=stats.accepted;r.rejected_steps=stats.rejected;r.rate_evaluations=stats.rateEvaluations;next[i]=cell(state);values[i]=r;
            }catch(const std::exception& e){status[i]=1;errors[i]=e.what();}
        }
        for(size_t i=0;i<count;i++)if(status[i]){lastError="Cell "+std::to_string(i)+": "+errors[i];return 1;}
        std::copy(next.begin(),next.end(),cells);std::copy(values.begin(),values.end(),results);return 0;
    }catch(const std::exception& e){lastError=e.what();return 1;}
}
