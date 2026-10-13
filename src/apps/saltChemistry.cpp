// SPDX-License-Identifier: GPL-3.0-or-later
#include "SaltChemistryC.h"
#include "BrineSystem.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <algorithm>
using namespace saltchem;
std::vector<std::string> fields(std::string line){std::vector<std::string> out;std::istringstream in(line);std::string s;while(std::getline(in,s,',')){if(!s.empty()&&s.back()=='\r')s.pop_back();out.push_back(s);}return out;}
double number(const std::string& s){size_t pos=0;double v=std::stod(s,&pos);if(pos!=s.size()||!std::isfinite(v))throw std::runtime_error("Invalid numeric input: "+s);return v;}
void check(int code){if(code)throw std::runtime_error(salt_chemistry_last_error());}
int main(int argc,char** argv){try{
    if(argc<5){std::cerr<<"Usage: saltChemistry database initial.csv rates.csv output-directory [duration_s=20] [output_interval_s=0.5] [threads=1] [temperature_C=25] [absolute_pressure_atm=1]\n";return 1;}
    double duration=argc>5?number(argv[5]):20,interval=argc>6?number(argv[6]):.5;int threads=argc>7?std::stoi(argv[7]):1;if(duration<0||interval<=0||threads<1)throw std::runtime_error("Invalid duration, output interval or thread count");
    std::ifstream initial(argv[2]),inputRates(argv[3]);if(!initial||!inputRates)throw std::runtime_error("Cannot read chemical input files");std::string line;std::getline(initial,line);
    const auto headings=fields(line);std::vector<std::string> expected{"water_kg","Na","K","Ca","Mg","Cl","SO4","C"};for(auto n:mineralNames)expected.emplace_back(n);if(headings!=expected)throw std::runtime_error("Incorrect initial.csv column order");
    std::getline(initial,line);auto values=fields(line);if(values.size()!=15)throw std::runtime_error("initial.csv requires 15 values");double t[7],mineral[7];double w=number(values[0]);for(int i=0;i<7;i++){t[i]=number(values[1+i]);mineral[i]=number(values[8+i]);}
    SaltRates rates{};std::array<bool,7> present{};std::getline(inputRates,line);
    const std::vector<std::string> rateHeadings{"mineral","k_neutral_mol_m2_s","k_acid","k_co2","acid_order","area_m2","reference_mol","dissolution_exponent","affinity_exponent"};if(fields(line)!=rateHeadings)throw std::runtime_error("Incorrect rates.csv column order");
    while(std::getline(inputRates,line)){if(line.empty()||line=="\r")continue;auto v=fields(line);if(v.size()!=9)throw std::runtime_error("Each rate row requires nine columns");int i=mineralIndex(v[0]);if(i<0||present[i])throw std::runtime_error("Unknown or repeated mineral: "+v[0]);present[i]=true;rates.k_neutral[i]=number(v[1]);rates.k_acid[i]=number(v[2]);rates.k_co2[i]=number(v[3]);rates.acid_order[i]=number(v[4]);rates.area_m2[i]=number(v[5]);rates.reference_mol[i]=number(v[6]);rates.dissolution_exponent[i]=number(v[7]);rates.affinity_exponent[i]=number(v[8]);}
    for(bool x:present)if(!x)throw std::runtime_error("All seven minerals must have explicit rate parameters; use zero area to disable a mineral");
    const double temperatureC=argc>8?number(argv[8]):25,pressureAtm=argc>9?number(argv[9]):1;
    SaltChemistryHandle* h=nullptr;check(salt_chemistry_create(argv[1],temperatureC+273.15,pressureAtm,&h));SaltCell cell;SaltResult result;check(salt_chemistry_initialize(h,t,w,mineral,&cell,&result));
    std::filesystem::path dest(argv[4]);std::filesystem::create_directories(dest);std::ofstream history(dest/"history.csv");history<<std::setprecision(16)<<"time_s";for(auto n:mineralNames)history<<','<<n;for(auto n:componentNames)history<<",m_"<<n;history<<",water_kg,density_kg_m3,pH,water_activity,aqueous_mass_kg";for(auto n:mineralNames)history<<",SI_"<<n;history<<'\n';
    auto write=[&](double time){history<<time;for(double x:cell.mineral_mol)history<<','<<x;for(double x:cell.aqueous_mol)history<<','<<x/result.free_water_kg;history<<','<<result.free_water_kg<<','<<result.density_kg_m3<<','<<result.pH<<','<<result.water_activity<<','<<result.aqueous_mass_kg;for(double x:result.saturation_index)history<<','<<x;history<<'\n';};
    SaltOptions options{2e-12,1e-15,100000,threads};int status=0;check(salt_chemistry_advance(h,1,&cell,&rates,0,&options,&result,&status));write(0);
    auto start=std::chrono::steady_clock::now();double time=0;int accepted=0,rejected=0;
    while(time<duration*(1-1e-13)){double dt=std::min(interval,duration-time);check(salt_chemistry_advance(h,1,&cell,&rates,dt,&options,&result,&status));time+=dt;accepted+=result.accepted_steps;rejected+=result.rejected_steps;write(time);}
    const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();std::ofstream metrics(dest/"run_metrics.csv");metrics<<std::setprecision(16)<<"metric,value\nphysical_time_s,"<<time<<"\nwall_seconds,"<<elapsed<<"\naccepted_steps,"<<accepted<<"\nrejected_steps,"<<rejected<<"\ntemperature_K,"<<temperatureC+273.15<<"\npressure_atm,"<<pressureAtm<<'\n';
    salt_chemistry_destroy(h);std::cout<<"Seven-mineral batch complete: physical_time_s="<<time<<" wall_s="<<elapsed<<" density_kg_m3="<<result.density_kg_m3<<" pH="<<result.pH<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
