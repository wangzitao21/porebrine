// SPDX-License-Identifier: GPL-3.0-or-later
#include "LocalWaterReaction.hpp"
#include "LocalChemistry.H"
#include <iostream>
#include <chrono>
int main(int argc,char** argv){try{
    if(argc!=4)throw std::invalid_argument("Supply database, geometry and chemical-state binaries");
    auto g=saltdual::Geometry::read(argv[2]);std::ifstream f(argv[3],std::ios::binary);std::uint32_t h[8]{};double v[8]{};f.read(reinterpret_cast<char*>(h),sizeof(h));f.read(reinterpret_cast<char*>(v),sizeof(v));
    if(h[0]!=0x57415443||h[1]!=1||h[2]!=g.nx*g.ny)throw std::runtime_error("Invalid chemical replay header");
    saltchem::RateParameters rates;rates.geometry.enabled=true;rates.geometry.epsilonMin=v[3];rates.geometry.closureWidth=v[4];rates.geometry.directionWidth=v[5];
    saltchem::IntegrationOptions options;options.maximumAttempts=h[4];options.method=h[5];options.analyticJacobian=h[6];options.relativeTolerance=v[6];options.absoluteTolerance=v[7];
    for(auto* vector:{&rates.neutral,&rates.acid,&rates.co2,&rates.acidOrder,&rates.dissolutionExponent,&rates.affinityExponent,&rates.geometry.weights})f.read(reinterpret_cast<char*>(vector->data()),sizeof(*vector));
    std::vector<saltchem::Inventory> cells(h[2]);std::vector<double> pressure(h[2]);for(std::uint32_t p=0;p<h[2];p++){f.read(reinterpret_cast<char*>(&pressure[p]),8);saltdual::WaterStock stock{};f.read(reinterpret_cast<char*>(stock.data()),sizeof(stock));saltdual::setWater(cells[p],stock);f.read(reinterpret_cast<char*>(cells[p].mineral.data()),sizeof(cells[p].mineral));}
    if(!f)throw std::runtime_error("Incomplete chemical replay data");
    saltchem::PitzerModel model(argv[1]);model.setConditions(v[1],1);saltflow::LocalChemistry workers(model,1);saltdual::ReactionMap map;map.rebuild(g);auto tasks=saltdual::waterTasks(g,map);if(h[3]>=tasks.size())throw std::runtime_error("Replay task no longer matches geometry");
    std::cout<<"regions="<<tasks[h[3]].regions.size()<<" links="<<tasks[h[3]].links.size()<<" dt="<<v[0]<<std::endl;auto start=std::chrono::steady_clock::now();
    for(int r:tasks[h[3]].regions){const auto& region=g.waterRegions[r];std::cout<<"initial_region="<<r<<" parent="<<region.parent<<" bulk="<<region.bulk<<" connected="<<region.connected<<" volume="<<region.fraction<<" mass="<<saltdual::waterMass(saltdual::regionStock(g,cells,r))<<std::endl;
        for(const auto& edge:g.waterLinks)if(edge.a==r||edge.b==r){int other=edge.a==r?edge.b:edge.a;if(std::find(tasks[h[3]].regions.begin(),tasks[h[3]].regions.end(),other)==tasks[h[3]].regions.end())std::cout<<"omitted_neighbour="<<other<<" volume="<<g.waterRegions[other].fraction<<" bulk="<<g.waterRegions[other].bulk<<std::endl;}}
    auto result=saltdual::reactWaterTask(tasks[h[3]],g,map,cells,rates,model.molarVolume,v[2],v[0],options,[&](int p)->saltflow::ChemicalWorker&{return workers.at(v[1],pressure[p]);});
    double residual=0,error=0;long evaluations=0;for(const auto& r:result){residual=std::max(residual,r.nonlinearResidualRatio);error=std::max(error,r.nonlinearStateErrorRatio);evaluations+=r.stats.rateEvaluations;}
    std::cout<<"PASS chemical replay state_error_ratio="<<error<<" raw_residual_ratio="<<residual<<" evaluations="<<evaluations<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
    for(const auto& r:result)if(!g.waterRegions[r.region].bulk){double volume=g.waterRegions[r.region].fraction;for(int m=0;m<3;m++)volume+=(r.stats.dissolved[m]-r.stats.precipitated[m])*model.molarVolume[m];saltchem::Inventory s;saltdual::setWater(s,r.water);double water;auto density=workers.at(v[1],pressure[g.waterRegions[r.region].parent]).speciation.conserved(s,water).density;std::cout<<"local_region="<<r.region<<" inventory_density="<<saltdual::waterMass(r.water)/volume<<" eos_density="<<density<<'\n';}
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
