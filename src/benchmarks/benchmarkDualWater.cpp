// SPDX-License-Identifier: GPL-3.0-or-later
#include "LocalWaterReaction.hpp"
#include "LocalChemistry.H"
#include <chrono>
#include <iostream>
using namespace saltdual;
using namespace saltchem;
int main(int argc,char** argv){try{
    if(argc!=3)throw std::invalid_argument("Supply pitzer.dat and a dual geometry checkpoint");
    PitzerModel model(argv[1]);Speciation spec(model);saltflow::LocalChemistry workers(model,1);auto g=Geometry::read(argv[2]);
    auto liquid=spec.inventory(Totals{{6.1292275648037755,0,0,0,6.1292275648037755,0,0}},1,{});double w,density=spec.conserved(liquid,w).density;
    std::vector<Inventory> cells(g.nx*g.ny);
    for(int p=0;p<g.nx*g.ny;p++){cells[p]=liquid;double factor=(1-std::accumulate(g.fractions[p].begin(),g.fractions[p].end(),0.))*density/aqueousMass(liquid);for(double& a:cells[p].aqueous)a*=factor;cells[p].waterEquivalent*=factor;for(int m=0;m<3;m++)cells[p].mineral[m]=g.fractions[p][m]/model.molarVolume[m];}
    initializeWater(g,cells,model.molarVolume);ReactionMap map;map.rebuild(g);auto tasks=waterTasks(g,map);
    RateParameters rates;rates.geometry.enabled=true;rates.geometry.epsilonMin=1e-5;rates.neutral[0]=.02;rates.neutral[1]=.01;rates.neutral[2]=.05;IntegrationOptions options;options.relativeTolerance=2e-8;options.absoluteTolerance=1e-11;
    std::cout<<"regions="<<g.waterRegions.size()<<" local_pools="<<g.pockets.size()<<" tasks="<<tasks.size()<<std::endl;
    for(std::size_t k=0;k<tasks.size();k++){
        const auto& task=tasks[k];std::cout<<"task="<<k<<" regions="<<task.regions.size()<<" links="<<task.links.size()<<std::flush;
        auto start=std::chrono::steady_clock::now();auto result=reactWaterTask(task,g,map,cells,rates,model.molarVolume,1e-9,.01,options,[&](int)->saltflow::ChemicalWorker&{return workers.at(298.15,1);});
        double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();int evaluations=0;for(const auto& r:result)evaluations+=r.stats.rateEvaluations;
        std::cout<<" seconds="<<seconds<<" evaluations="<<evaluations<<std::endl;
    }
}catch(const std::exception& e){std::cerr<<'\n'<<e.what()<<'\n';return 1;}}
