// SPDX-License-Identifier: GPL-3.0-or-later
#include "LocalWaterReaction.hpp"
#include "LocalChemistry.H"
#include <iostream>
using namespace saltdual;
using namespace saltchem;
static void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
static std::array<long double,8> total(const Geometry& g,const std::vector<Inventory>& cells){
    std::array<long double,8> result{};
    for(const auto& c:cells){for(int j=0;j<7;j++)result[j]+=c.aqueous[j];result[7]+=c.waterEquivalent;
        for(int m=0;m<7;m++){for(int j=0;j<7;j++)result[j]+=stoichiometry[m][j]*c.mineral[m];result[7]+=hydration[m]*waterMolarMass*c.mineral[m];}}
    for(const auto& p:g.pockets)for(int j=0;j<8;j++)result[j]+=p.inventory[j];return result;
}
int main(int argc,char** argv){try{
    require(argc==2,"Supply pitzer.dat");PitzerModel model(argv[1]);Speciation spec(model);saltflow::LocalChemistry worker(model,1);
    {
        Geometry graph;graph.waterRegions.resize(3);graph.waterLinks={{0,1,1,1},{1,2,1,3}};
        std::vector<Inventory> nodes(3);for(int i=0;i<3;i++){nodes[i].aqueous[0]=nodes[i].aqueous[4]=i+1;nodes[i].waterEquivalent=1-(i+1)*(componentMass[0]+componentMass[4]);}
        std::array<long double,8> initial{};for(const auto& node:nodes)for(int j=0;j<8;j++)initial[j]+=waterStock(node)[j];
        redistributeReactionWater(graph,{0,1,0},{1,.4,1},nodes);
        require(std::abs(aqueousMass(nodes[0])-1.15)<1e-13&&std::abs(aqueousMass(nodes[1])-.4)<1e-13&&std::abs(aqueousMass(nodes[2])-1.45)<1e-13,"Pore-volume displacement ignored hydraulic throat competition");
        redistributeReactionWater(graph,{0,1,0},{1,.8,1},nodes);
        for(int j=0;j<8;j++){long double after=0;for(const auto& node:nodes)after+=waterStock(node)[j];require(std::abs(after-initial[j])<1e-13,"Pore-volume displacement lost a component");}
        require(std::abs(aqueousMass(nodes[1])-.8)<1e-13,"Pore-volume refill failed");
    }
    {
        // An internally connected sealed pore needs displacement even with
        // zero diffusion. The pressure gauge must not select its chemistry.
        Geometry graph;graph.waterRegions.resize(3);graph.waterLinks={{0,1,1,1},{1,2,1,3}};
        std::vector<Inventory> original(3);for(int i=0;i<3;i++){original[i].aqueous[0]=original[i].aqueous[4]=i+.5;original[i].waterEquivalent=1-(i+.5)*(componentMass[0]+componentMass[4]);}
        std::array<long double,8> initial{};for(const auto& node:original)for(int j=0;j<8;j++)initial[j]+=waterStock(node)[j];
        auto forward=original;redistributeReactionWater(graph,{1,1,1},{.5,1.5,2.},forward);
        for(int i=0;i<3;i++)require(std::abs(aqueousMass(forward[i])-.75*std::array<double,3>{.5,1.5,2.}[i])<1e-13,"Closed water failed to redistribute finite mass by liquid volume");
        auto reverse=original;std::reverse(reverse.begin(),reverse.end());Geometry reversed;reversed.waterRegions.resize(3);reversed.waterLinks={{2,1,1,1},{1,0,1,3}};
        redistributeReactionWater(reversed,{1,1,1},{2.,1.5,.5},reverse);
        for(int j=0;j<8;j++){long double after=0;for(int i=0;i<3;i++){after+=waterStock(forward[i])[j];require(std::abs(waterStock(forward[i])[j]-waterStock(reverse[2-i])[j])<1e-13,"Closed-water displacement depends on node order");}require(std::abs(after-initial[j])<1e-13,"Closed-water displacement lost a component");}
    }
    {
        // A just-created connected sliver can have only roundoff inventory.
        // It shares the finite closed pore's solution, rather than presenting
        // an almost massless reactor to the precipitation rate law.
        Geometry g(1,1,4,5e-6,20e-6,1e-5);g.waterRegions={{0,false,false,6e-6,{},{}},{0,false,false,1.04e-10,{},{}}};
        g.waterLinks={{0,1,1e-8,1e-8}};g.bulkRegion={-1};g.pocketRegion={0,1};g.pockets={{0,6e-6,{}},{0,1.04e-10,{}}};g.waterInitialized=true;
        auto brine=spec.inventory(Totals{{6,0,0,0,6,0,0}},1,{});double water;double density=spec.conserved(brine,water).density;
        auto stock=waterStock(brine);for(double& value:stock)value*=density*6e-6/aqueousMass(brine);g.pockets[0].inventory=stock;for(double& value:stock)value*=1e-21;g.pockets[1].inventory=stock;
        std::vector<Inventory> cells(1);ReactionMap map;map.patches.resize(2);RateParameters rates;IntegrationOptions options;const auto before=total(g,cells);
        auto result=reactWaterTask(WaterTask{{0,1},{0}},g,map,cells,rates,model.molarVolume,0,.01,options,[&](int)->saltflow::ChemicalWorker&{return worker.at(298.15,1);});
        require(std::abs(waterMass(result[0].water)/6e-6-waterMass(result[1].water)/1.04e-10)<1e-8,"Closed connected sliver stayed almost massless");
        commitRegions(g,cells,result,model.molarVolume);const auto after=total(g,cells);for(int j=0;j<8;j++)require(std::abs(after[j]-before[j])<1e-13,"Closed sliver refill lost a component");
    }
    Geometry initial(3,1,4,5e-6,20e-6,1e-5);
    for(int i=0;i<int(initial.cells.size());i++){auto o=initial.origin(i);paint(initial.cells[i],-1,2,{{{1,0},3.9-o.x}});merge(initial.cells[i]);}initial.rebuild();
    auto brine=spec.inventory(Totals{{6.1292275648,0,0,0,6.1292275648,0,0}},1,{});double water;double rho=spec.conserved(brine,water).density;
    std::vector<Inventory> starting(3);
    for(int p=0;p<3;p++){starting[p]=brine;double factor=(1-initial.fractions[p][2])*rho/aqueousMass(brine);for(double& a:starting[p].aqueous)a*=factor;starting[p].waterEquivalent*=factor;starting[p].mineral[2]=initial.fractions[p][2]/model.molarVolume[2];}
    RateParameters rates;rates.geometry.enabled=true;rates.geometry.epsilonMin=1e-5;rates.neutral[2]=.05;IntegrationOptions options;
    std::vector<double> dissolved;
    // Fixed geometry and interval; vary only the physical exchange coefficient.
    for(double diffusion:{0.,1e-9,1e-6}){
        auto g=initial;auto cells=starting;initializeWater(g,cells,model.molarVolume);ReactionMap map;map.rebuild(g);const auto before=total(g,cells);
        std::vector<RegionReaction> results;
        for(const auto& task:waterTasks(g,map)){
            auto part=reactWaterTask(task,g,map,cells,rates,model.molarVolume,diffusion,.01,options,[&](int)->saltflow::ChemicalWorker&{return worker.at(298.15,1);});
            results.insert(results.end(),part.begin(),part.end());
        }
        commitRegions(g,cells,results,model.molarVolume);auto after=total(g,cells);
        for(int j=0;j<8;j++)require(std::abs(after[j]-before[j])<1e-11*std::max(1.L,std::abs(before[j])),"Coupled water exchange and reaction lost a component");
        dissolved.push_back(starting[0].mineral[2]-cells[0].mineral[2]);
        require(dissolved.back()>0,"Thin-water mineral reaction stalled");
        std::cout<<"diffusivity="<<diffusion<<" dissolved="<<dissolved.back()<<'\n';
    }
    require(dissolved[1]>dissolved[0]&&dissolved[2]>dissolved[1],"Opening water exchange failed to increase sustained dissolution");
    {
        Geometry g(3,1,4,5e-6,20e-6,1e-5);g.waterCouplingFraction=1;
        g.circle(2,2.5e-6,2.5e-6,1.5e-6);g.circle(1,7.5e-6,2.5e-6,1.5e-6);g.circle(0,12.5e-6,2.5e-6,1.5e-6);g.rebuild();
        auto brine=spec.inventory(Totals{{6,3,0,.3,9.6,0,0}},1,{});double rho=spec.conserved(brine,water).density;std::vector<Inventory> original(3);
        for(int p=0;p<3;p++){original[p]=brine;double scale=(1-std::accumulate(g.fractions[p].begin(),g.fractions[p].end(),0.))*rho/aqueousMass(brine);for(double& a:original[p].aqueous)a*=scale;original[p].waterEquivalent*=scale;for(int m=0;m<3;m++)original[p].mineral[m]=g.fractions[p][m]/model.molarVolume[m];}
        std::array<double,2> retreat{},potassium{};
        // A fixed geometry and interval isolate the chemical effect of KCl
        // precipitation on the neighbouring carnallite reaction.
        for(int active=0;active<2;active++){
            auto geometry=g;auto cells=original;initializeWater(geometry,cells,model.molarVolume);ReactionMap map;map.rebuild(geometry);auto localRates=rates;localRates.neutral[1]=active?1.:0.;const auto before=total(geometry,cells);
            std::vector<RegionReaction> results;for(const auto& task:waterTasks(geometry,map)){auto part=reactWaterTask(task,geometry,map,cells,localRates,model.molarVolume,1e-9,.01,options,[&](int)->saltflow::ChemicalWorker&{return worker.at(298.15,1);});results.insert(results.end(),part.begin(),part.end());}
            double crystallized=0;for(const auto& node:results){retreat[active]+=node.stats.dissolved[2]-node.stats.precipitated[2];crystallized+=node.stats.precipitated[1];if(node.region==geometry.bulkRegion[0]){Inventory state;setWater(state,node.water);double freeWater;auto solution=spec.conserved(state,freeWater);potassium[active]=componentTotals(solution.molality)[1];}}
            if(active)require(crystallized>0,"Feedback fixture failed to precipitate Sylvite");commitRegions(geometry,cells,results,model.molarVolume);const auto after=total(geometry,cells);
            for(int j=0;j<8;j++)require(std::abs(after[j]-before[j])<1e-11*std::max(1.L,std::abs(before[j])),"Precipitation feedback lost a conserved component");
        }
        require(potassium[1]<potassium[0]&&retreat[1]>retreat[0],"K consumption failed to accelerate the coupled dissolution reaction");
    }
    for(int mineral=0;mineral<3;mineral++){
        Geometry g(2,1,4,4e-6,20e-6,1e-5);g.waterCouplingFraction=1;g.circle(mineral,2e-6,2e-6,.1e-6,32);g.rebuild();
        auto fresh=spec.inventory(Totals{{.1,0,0,0,.1,0,0}},1,{});double rho=spec.conserved(fresh,water).density;std::vector<Inventory> cells(2);
        for(int p=0;p<2;p++){cells[p]=fresh;double scale=(1-g.fractions[p][mineral])*rho/aqueousMass(fresh);for(double& a:cells[p].aqueous)a*=scale;cells[p].waterEquivalent*=scale;cells[p].mineral[mineral]=g.fractions[p][mineral]/model.molarVolume[mineral];}
        initializeWater(g,cells,model.molarVolume);ReactionMap map;map.rebuild(g);const auto before=total(g,cells);auto localRates=rates;localRates.neutral={};localRates.neutral[mineral]=1;
        std::vector<RegionReaction> results;for(const auto& task:waterTasks(g,map)){auto part=reactWaterTask(task,g,map,cells,localRates,model.molarVolume,1e-9,.01,options,[&](int)->saltflow::ChemicalWorker&{return worker.at(298.15,1);});results.insert(results.end(),part.begin(),part.end());}
        commitRegions(g,cells,results,model.molarVolume);const auto after=total(g,cells);require(cells[0].mineral[mineral]<1e-10,"Implicit local reaction failed to exhaust a small mineral");
        for(int j=0;j<8;j++)require(std::abs(after[j]-before[j])<1e-11*std::max(1.L,std::abs(before[j])),"Local depletion lost a conserved component");
    }
    {
        Geometry g(2,1,4,4e-6,20e-6,1e-5);g.waterCouplingFraction=1;
        for(int cell=0;cell<int(g.cells.size());cell++){auto o=g.origin(cell);paint(g.cells[cell],-1,0,{{{0,1},2-o.y}});merge(g.cells[cell]);}g.rebuild();
        auto fresh=spec.inventory(Totals{{.1,0,0,0,.1,0,0}},1,{});double rho=spec.conserved(fresh,water).density;std::vector<Inventory> cells(2);
        for(int p=0;p<2;p++){cells[p]=fresh;double scale=(1-g.fractions[p][0])*rho/aqueousMass(fresh);for(double& a:cells[p].aqueous)a*=scale;cells[p].waterEquivalent*=scale;cells[p].mineral[0]=g.fractions[p][0]/model.molarVolume[0];}
        initializeWater(g,cells,model.molarVolume);ReactionMap map;map.rebuild(g);RateParameters stiff=rates;stiff.neutral={};stiff.neutral[0]=1000;
        auto tasks=waterTasks(g,map);require(tasks.size()==1&&tasks[0].regions.size()==2,"Stiff symmetric-pool fixture is disconnected");
        const int r=tasks[0].regions[0];auto problem=regionProblem(r,g,map,cells,stiff,model.molarVolume);auto& chemistry=worker.at(298.15,1).reactions;
        // Independent scalar bracketing gives the symmetric implicit solution.
        double lo=0,hi=problem.initial.mineral[0];for(int iteration=0;iteration<90;iteration++){
            double extent=.5*(lo+hi);MineralVector change{};change[0]=extent;auto state=Reactions::changed(problem.initial,change);
            double ep=(g.waterRegions[r].fraction+g.epsilon+extent*model.molarVolume[0])/problem.rates.geometry.bulkVolume;
            double residual=extent-.01*chemistry.evaluate(state,problem.rates,nullptr,false,ep).rates[0];if(residual>0)hi=extent;else lo=extent;
        }
        auto result=reactWaterTask(tasks[0],g,map,cells,stiff,model.molarVolume,1e-9,.01,options,[&](int)->saltflow::ChemicalWorker&{return worker.at(298.15,1);});
        for(const auto& node:result)require(std::abs(node.stats.dissolved[0]-.5*(lo+hi))<1e-8,"Stiff coupled chemistry disagrees with independently bracketed equilibrium approach");
    }
    {
        // A fixed two-reservoir problem: two salts share a nearly filled
        // retained pore and draw solution from the same finite bulk water.
        Geometry g(2,1,4,5e-6,20e-6,1e-5);const double pore=6.86e-7;
        g.waterRegions={{0,true,true,1,{},{}},{1,false,true,pore,{},{}}};g.waterLinks={{0,1,1e-8,1e-8}};
        g.bulkRegion={0,-1};g.pocketRegion={-1,0};g.pockets={{1,pore,{}}};g.waterInitialized=true;
        auto brine=spec.inventory(Totals{{5.2,1.6,0,1.6,10,0,0}},1,{});double density=spec.conserved(brine,water).density;
        auto stock=waterStock(brine);for(double& v:stock)v*=density/aqueousMass(brine);
        std::vector<Inventory> cells(2);setWater(cells[0],stock);setWater(cells[1],{});
        auto trapped=stock;for(double& v:trapped)v*=pore;g.pockets[0].inventory=trapped;waterLedger(g);
        cells[1].mineral[0]=1.;cells[1].mineral[1]=.1;
        ReactionMap map;map.patches.resize(2);map.patches[1]={{1,0,1e-12,1e-12,1,{},{}},{1,1,5e-13,5e-13,1,{},{}}};map.active={1};
        auto localRates=rates;localRates.neutral={};localRates.neutral[0]=.02;localRates.neutral[1]=.01;
        auto localOptions=options;localOptions.relativeTolerance=2e-8;localOptions.absoluteTolerance=1e-11;
        const auto before=total(g,cells);auto result=reactWaterTask(WaterTask{{0,1},{0}},g,map,cells,localRates,model.molarVolume,1e-9,.01,localOptions,[&](int)->saltflow::ChemicalWorker&{return worker.at(298.15,1);});
        double remaining=pore;std::array<double,2> grown{};
        for(const auto& node:result)if(node.region==1)for(int m=0;m<2;m++){grown[m]=node.stats.precipitated[m];remaining+=(node.stats.dissolved[m]-node.stats.precipitated[m])*model.molarVolume[m];}
        require(grown[0]>0&&grown[1]>0&&remaining>0&&remaining<.1*pore,"Competing salts failed to fill a retained pore while preserving its remaining water");
        commitRegions(g,cells,result,model.molarVolume);const auto after=total(g,cells);
        for(int j=0;j<8;j++)require(std::abs(after[j]-before[j])<1e-11*std::max(1.L,std::abs(before[j])),"Competing retained-pore precipitation lost a component");
    }
    std::cout<<"PASS coupled water-reaction solve: finite exchange, thin-water limit, competing pore closure, three-mineral depletion and all-component conservation\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
