// SPDX-License-Identifier: GPL-3.0-or-later
#include "SaltChemistryC.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <vector>
#include <chrono>
#include <stdexcept>
#include <iomanip>
#include <filesystem>
void require(bool yes,const char* message){if(!yes)throw std::runtime_error(message);}
void check(int status){if(status)throw std::runtime_error(salt_chemistry_last_error());}
SaltRates rate(double scale){SaltRates p{};const double k[7]={.02,.002,.02,1e-5,1e-4,.02,1e-5};for(int j=0;j<7;j++){p.k_neutral[j]=k[j];p.area_m2[j]=scale;p.reference_mol[j]=scale;p.acid_order[j]=p.dissolution_exponent[j]=p.affinity_exponent[j]=1;}return p;}
int main(int argc,char** argv){try{
    const char* db=argc>1?argv[1]:"src/chemistry/pitzer.dat";std::string dest=argc>2?argv[2]:"results/checks/api";std::filesystem::create_directories(dest);
    SaltChemistryHandle* h=nullptr;check(salt_chemistry_create(db,298.15,1,&h));
    const double t[7]={4,.2,.1,.5,5.328,.035,.002};double mineral[7]={1,1,1,1,1,1,1};SaltCell base;SaltResult initial;
    check(salt_chemistry_initialize(h,t,1,mineral,&base,&initial));require(std::abs(initial.free_water_kg-1)<1e-10,"Free-water initialization mismatch");
    SaltOptions options{2e-10,1e-13,10000,1};const double dt=.5;
    check(salt_chemistry_advance(h,0,nullptr,nullptr,dt,&options,nullptr,nullptr));
    std::vector<SaltCell> scaleCells;std::vector<SaltRates> scaleRates;std::vector<double> scales;
    for(int e:{0,-3,-9,-15,3}){double s=std::pow(10.,e);SaltCell c=base;for(int j=0;j<7;j++){c.aqueous_mol[j]*=s;c.mineral_mol[j]*=s;}c.water_equivalent_kg*=s;scaleCells.push_back(c);scaleRates.push_back(rate(s));scales.push_back(s);}
    std::vector<SaltResult> output(scales.size());std::vector<int> status(scales.size());check(salt_chemistry_advance(h,scales.size(),scaleCells.data(),scaleRates.data(),dt,&options,output.data(),status.data()));
    double scalingError=0;
    for(size_t i=0;i<scales.size();i++){for(int j=0;j<7;j++){scalingError=std::max(scalingError,std::abs(scaleCells[i].aqueous_mol[j]/scales[i]-scaleCells[0].aqueous_mol[j]));scalingError=std::max(scalingError,std::abs(scaleCells[i].mineral_mol[j]/scales[i]-scaleCells[0].mineral_mol[j]));}scalingError=std::max(scalingError,std::abs(output[i].density_kg_m3/output[0].density_kg_m3-1));}
    require(scalingError<1e-9,"Cell-size scaling changed the chemistry");
    // Invalid cells must abort the whole batch without any partially advanced output.
    auto invalid=scaleCells;invalid[2].aqueous_mol[0]=-1;auto saved=invalid;auto savedOutput=output;
    require(salt_chemistry_advance(h,scales.size(),invalid.data(),scaleRates.data(),dt,&options,output.data(),status.data())!=0,"Invalid batch was accepted");
    require(std::memcmp(invalid.data(),saved.data(),saved.size()*sizeof(SaltCell))==0,"Failed batch partially changed cell states");
    require(std::memcmp(output.data(),savedOutput.data(),savedOutput.size()*sizeof(SaltResult))==0,"Failed batch partially changed outputs");
    // Same varied batch, serial and parallel: test every conserved state and density.
    const size_t count=256;std::vector<SaltCell> cells(count);std::vector<SaltRates> rates(count);std::vector<SaltResult> result(count);status.resize(count);
    for(size_t i=0;i<count;i++){double total[7];std::copy_n(t,7,total);double factor=.5+(i%31)/30.;total[0]*=factor;total[4]+=total[0]-t[0];double m[7];for(int j=0;j<7;j++)m[j]=.01+((i+3*j)%17)/16.;SaltResult r;check(salt_chemistry_initialize(h,total,1,m,&cells[i],&r));rates[i]=rate(1);}
    const auto start=cells;const auto tick=std::chrono::steady_clock::now();check(salt_chemistry_advance(h,count,cells.data(),rates.data(),dt,&options,result.data(),status.data()));double serial=std::chrono::duration<double>(std::chrono::steady_clock::now()-tick).count();
    const auto serialCells=cells;const auto serialResults=result;cells=start;options.threads=8;auto parallelStart=std::chrono::steady_clock::now();check(salt_chemistry_advance(h,count,cells.data(),rates.data(),dt,&options,result.data(),status.data()));double parallel=std::chrono::duration<double>(std::chrono::steady_clock::now()-parallelStart).count();
    double parallelError=0,massSourceError=0;const double mineralMass[7]={.0584428,.0745513,.2778583,.6029496,.172176,.203307,.1000911};
    for(size_t i=0;i<count;i++){
        for(int j=0;j<7;j++){parallelError=std::max(parallelError,std::abs(cells[i].aqueous_mol[j]-serialCells[i].aqueous_mol[j]));parallelError=std::max(parallelError,std::abs(cells[i].mineral_mol[j]-serialCells[i].mineral_mol[j]));}
        parallelError=std::max(parallelError,std::abs(result[i].density_kg_m3-serialResults[i].density_kg_m3));double source=result[i].aqueous_mass_source_kg_s;for(int j=0;j<7;j++)source+=mineralMass[j]*result[i].mineral_change_mol[j]/dt;massSourceError=std::max(massSourceError,std::abs(source));
    }
    require(parallelError<1e-12,"Parallel execution changed the chemistry");require(massSourceError<1e-10,"Aqueous mass source disagrees with mineral mass transfer");
    salt_chemistry_destroy(h);
    std::ofstream f(dest+"/api_metrics.csv");f<<std::setprecision(16)<<"metric,value\nscale_invariance_error,"<<scalingError<<"\nparallel_error,"<<parallelError<<"\nmass_source_error_kg_s,"<<massSourceError<<"\nbatch_cells,"<<count<<"\nserial_seconds,"<<serial<<"\neight_threads_seconds,"<<parallel<<"\ntransaction_rollback_pass,1\n";
    std::cout<<"API PASS scale_error="<<scalingError<<" parallel_error="<<parallelError<<" source_error="<<massSourceError<<" serial_s="<<serial<<" parallel_s="<<parallel<<std::endl;return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
