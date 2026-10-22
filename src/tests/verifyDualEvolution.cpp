// SPDX-License-Identifier: GPL-3.0-or-later
#include "LocalWaterReaction.hpp"
#include "LocalChemistry.H"
#include <iostream>
#include <chrono>
using namespace saltdual;
using namespace saltchem;
static void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
static MineralVector mineralVolumes;
static std::array<long double,8> total(const Geometry& g,const std::vector<Inventory>& cells){
    std::array<long double,8> sum{};for(const auto& cell:cells){for(int j=0;j<7;j++)sum[j]+=cell.aqueous[j];sum[7]+=cell.waterEquivalent;for(int m=0;m<3;m++){for(int j=0;j<7;j++)sum[j]+=stoichiometry[m][j]*cell.mineral[m];sum[7]+=hydration[m]*waterMolarMass*cell.mineral[m];}}
    for(const auto& p:g.pockets)for(int j=0;j<8;j++)sum[j]+=p.inventory[j];auto debt=mineralDebtComponents(g,mineralVolumes);for(int j=0;j<8;j++)sum[j]+=debt[j];return sum;
}
int main(int argc,char** argv){try{
    require(argc==2,"Supply pitzer.dat");PitzerModel model(argv[1]);mineralVolumes=model.molarVolume;Speciation spec(model);saltflow::LocalChemistry workers(model,1);omp_set_num_threads(1);
    Geometry g(16,12,4,5e-6,20e-6,1e-5);g.circle(2,40e-6,30e-6,12.5e-6,128);g.rebuild();
    auto solution=spec.inventory(Totals{{6,1.5,0,1.5,10.5,0,0}},1,{});double water,density=spec.conserved(solution,water).density;
    std::vector<Inventory> cells(g.nx*g.ny);
    for(int p=0;p<g.nx*g.ny;p++){cells[p]=solution;double factor=(1-g.fractions[p][2])*density/aqueousMass(solution);for(double& a:cells[p].aqueous)a*=factor;cells[p].waterEquivalent*=factor;cells[p].mineral[2]=g.fractions[p][2]/model.molarVolume[2];}
    initializeWater(g,cells,model.molarVolume);RateParameters rates;rates.geometry.enabled=true;rates.geometry.epsilonMin=1e-5;rates.neutral[0]=.02;rates.neutral[1]=.01;rates.neutral[2]=.05;
    IntegrationOptions options;options.relativeTolerance=2e-8;options.absoluteTolerance=1e-11;
    std::vector<saltgeom::NucleationLaw> laws;for(int substrate=0;substrate<3;substrate++)for(int product=0;product<3;product++)if(substrate!=product){saltgeom::NucleationLaw law;law.product=product;law.substrate=substrate;law.prefactor=1e11;law.barrier=.02;law.radius=1e-7;law.contactAngle=60;laws.push_back(law);}
    int rejected=0;double maxBalance=0,coveredBeforeFlush=0;std::uint64_t nuclei=0;auto begin=std::chrono::steady_clock::now();
    for(int frame=0;frame<20;frame++){
        if(frame==10){
            for(const auto& a:g.covered)for(double area:a)coveredBeforeFlush+=area;
            auto fresh=spec.inventory(Totals{{.1,0,0,0,.1,0,0}},1,{});double mass=aqueousMass(fresh);
            for(auto& cell:cells){double factor=aqueousMass(cell)/mass;auto incoming=waterStock(fresh);for(double& a:incoming)a*=factor;setWater(cell,incoming);}
        }
        const auto before=total(g,cells);double remaining=.01,step=.01;int attempts=0;
        while(remaining>1e-13){
            if(++attempts>256)throw std::runtime_error("Fixed-grid evolution exceeded its event budget");step=std::min(step,remaining);const auto original=g;const auto stock=cells;
            std::vector<SurfaceBudget> attempted;std::vector<std::array<double,3>> desired;
            try{
                applyMineralDebt(g,cells,model.molarVolume);
                ReactionMap map;map.rebuild(g);std::vector<RegionReaction> results;
                for(const auto& task:waterTasks(g,map)){auto part=reactWaterTask(task,g,map,cells,rates,model.molarVolume,1e-9,step,options,[&](int)->saltflow::ChemicalWorker&{return workers.at(298.15,1);});results.insert(results.end(),part.begin(),part.end());}
                attempted=commitRegions(g,cells,results,model.molarVolume);
                desired=g.fractions;for(int p=0;p<g.nx*g.ny;p++)for(int m=0;m<3;m++)desired[p][m]=cells[p].mineral[m]*model.molarVolume[m];
                g.advance(desired,false,.2,&attempted);for(int p=0;p<g.nx*g.ny;p++)for(int m=0;m<3;m++)cells[p].mineral[m]=g.fractions[p][m]/model.molarVolume[m];remapWater(original,g,cells,&g.pockets);
                auto born=nucleateSurface(g,map,cells,laws,model.molarVolume,step,298.15,2024,[&](int)->saltflow::ChemicalWorker&{return workers.at(298.15,1);});for(auto n:born.count)nuclei+=n;
                for(int r=0;r<int(g.waterRegions.size());r++)if(resolvedWater(g,r))require(regionStock(g,cells,r)[7]>0,"Newly dissolved volume received no solvent");
                remaining-=step;step=std::min(.01,step*1.5);
            }catch(const StepFailure& e){
                if(std::getenv("SALT_DUAL_CAPTURE")&&std::string(e.what()).find("Fine geometry")!=std::string::npos){
                    original.write("bin/dual_evolution_failure.bin");g.write("bin/dual_evolution_after.bin");std::ofstream f("bin/dual_evolution_budget.bin",std::ios::binary);
                    for(const auto& a:desired)f.write(reinterpret_cast<const char*>(a.data()),sizeof(a));std::uint32_t count=attempted.size();f.write(reinterpret_cast<const char*>(&count),4);
                    for(const auto& b:attempted){std::int32_t values[]{b.parent,b.mineral,b.liquidRegion,int(b.surfaces.size())};f.write(reinterpret_cast<const char*>(values),sizeof(values));f.write(reinterpret_cast<const char*>(&b.change),8);f.write(reinterpret_cast<const char*>(b.surfaces.data()),4*b.surfaces.size());}throw;
                }
                g=original;cells=stock;step*=std::clamp(e.scale,1e-6,.999999);rejected++;if(rejected<=5||rejected%32==0)std::cerr<<"event step="<<step<<" reason="<<e.what()<<std::endl;if(step<1e-12)throw;
            }
        }
        auto after=total(g,cells);for(int j=0;j<8;j++){double error=double(std::abs(after[j]-before[j])/std::max(1.L,std::abs(before[j])));maxBalance=std::max(maxBalance,error);require(error<1e-10,"Evolving local water lost conserved mass");}
        std::cout<<"time="<<(frame+1)*.01<<" local_water="<<g.pockets.size()<<" nuclei="<<nuclei<<" retries="<<rejected<<std::endl;
    }
    require(nuclei>0,"Coupled brine evolution produced no new mineral sites");
    double finalCovered=0;for(const auto& a:g.covered)for(double area:a)finalCovered+=area;
    require(coveredBeforeFlush>0&&finalCovered<coveredBeforeFlush,"Changing the bath failed to dissolve coverage and restore exposure");
    std::cout<<"PASS evolving local water and three-mineral interfaces: balance="<<maxBalance<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
