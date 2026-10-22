// SPDX-License-Identifier: GPL-3.0-or-later
#include "DualSurfaceChemistry.hpp"
#include "LocalChemistry.H"
#include <iostream>
#include <cstdio>
using namespace saltdual;
using namespace saltchem;
static MineralVector inventoryVolumes;
void require(bool v,const char* message){if(!v)throw std::runtime_error(message);}
Geometry slab(double edge){
    Geometry g(3,1,4,5e-6,20e-6,1e-5);
    for(int k=0;k<int(g.cells.size());k++){auto o=g.origin(k);paint(g.cells[k],-1,2,{{{1,0},edge-o.x}});merge(g.cells[k]);}g.rebuild();return g;
}
std::vector<Inventory> initialize(const Geometry& g,const Speciation& spec,const MineralVector& vm,bool supersaturated=false){
    auto solution=spec.inventory(supersaturated?Totals{{6,1.5,0,1.5,10.5,0,0}}:Totals{{6.1292275648,0,0,0,6.1292275648,0,0}},1,{});
    double water;double density=spec.conserved(solution,water).density;std::vector<Inventory> cells(g.nx*g.ny);
    for(int p=0;p<g.nx*g.ny;p++){cells[p]=solution;double factor=(1-std::accumulate(g.fractions[p].begin(),g.fractions[p].end(),0.))*density/aqueousMass(solution);
        for(auto& a:cells[p].aqueous)a*=factor;cells[p].waterEquivalent*=factor;
        for(int m=0;m<3;m++)cells[p].mineral[m]=g.fractions[p][m]/vm[m];
    }return cells;
}
std::array<long double,8> inventory(const std::vector<Inventory>& cells,const Geometry* geometry=nullptr){
    std::array<long double,8> out{};for(const auto& c:cells){for(int j=0;j<7;j++)out[j]+=c.aqueous[j];out[7]+=c.waterEquivalent;
        for(int m=0;m<7;m++){for(int j=0;j<7;j++)out[j]+=c.mineral[m]*stoichiometry[m][j];out[7]+=c.mineral[m]*hydration[m]*waterMolarMass;}}
    if(geometry){for(const auto& pocket:geometry->pockets)for(int j=0;j<8;j++)out[j]+=pocket.inventory[j];auto debt=mineralDebtComponents(*geometry,inventoryVolumes);for(int j=0;j<8;j++)out[j]+=debt[j];}
    return out;
}
void balance(const std::vector<Inventory>& a,const std::vector<Inventory>& b){auto before=inventory(a),after=inventory(b);for(int j=0;j<8;j++)require(std::abs(before[j]-after[j])<3e-13*std::max(1.L,std::abs(before[j])),"Receiver/solid split lost a component or hydrate water");}
int main(int argc,char** argv){try{
    require(argc==2,"Supply pitzer.dat");PitzerModel model(argv[1]);inventoryVolumes=model.molarVolume;Speciation spec(model);Reactions reactions(model);saltflow::LocalChemistry workers(model,1);
    RateParameters rates;rates.geometry.enabled=true;rates.geometry.epsilonMin=1e-5;rates.neutral[2]=.05;
    IntegrationOptions options;options.method=2;
    for(int mineral=0;mineral<3;mineral++){
        Geometry hidden(1,1,4,4e-6,20e-6,1e-5);int cover=(mineral+1)%3;
        for(int cell=0;cell<int(hidden.cells.size());cell++){auto o=hidden.origin(cell);paint(hidden.cells[cell],-1,cover,{{{1,0},2-o.x},{{0,1},2-o.y}});paint(hidden.cells[cell],cover,mineral,{{{1,0},.75-o.x},{{-1,0},o.x-.25},{{0,1},.75-o.y},{{0,-1},o.y-.25}});merge(hidden.cells[cell]);}
        hidden.rebuild();double buried=hidden.fractions[0][mineral];hidden.circle(mineral,3e-6,3e-6,.02e-6,32);hidden.rebuild();auto initial=initialize(hidden,spec,model.molarVolume);initializeWater(hidden,initial,model.molarVolume);ReactionMap access;access.rebuild(hidden);
        double available=0;for(int r:access.active)available+=regionProblem(r,hidden,access,initial,rates,model.molarVolume).initial.mineral[mineral];
        require(std::abs(available*model.molarVolume[mineral]-(hidden.fractions[0][mineral]-buried))<1e-12,"A wet crystal borrowed mineral stock from a buried crystal in the same coarse cell");
    }
    for(double edge:{4.,3.999996,3.6}){
        std::cout<<"local slab edge="<<edge<<std::endl;
        auto g=slab(edge);ReactionMap map;map.rebuild(g);auto cells=initialize(g,spec,model.molarVolume);initializeWater(g,cells,model.molarVolume);
        require(g.openings[0][0]==0&&std::abs(g.openings[2][1]-1)<1e-12,"Boundary face aperture incorrect");
        const auto before=inventory(cells,&g);std::vector<RegionReaction> results;
        for(int r:map.active)results.push_back(reactRegion(r,g,map,cells,reactions,rates,model.molarVolume,.001,options));
        auto budgets=commitRegions(g,cells,results,model.molarVolume);auto after=inventory(cells,&g);
        for(int j=0;j<8;j++)require(std::abs(before[j]-after[j])<3e-12*std::max(1.L,std::abs(before[j])),"Local aqueous reaction lost component or hydrate water");
        require(!budgets.empty(),"Exposed mineral failed to react with its adjacent liquid");
    }
    {auto share=boundedPartition(1.5,{1,1},{.2,2});require(std::abs(share[0]-.2)<1e-14&&std::abs(share[1]-1.3)<1e-14,"Bounded area partition failed at mineral exhaustion");}
    std::vector<saltgeom::NucleationLaw> laws;
    for(int m=0;m<2;m++){saltgeom::NucleationLaw law;law.product=m;law.substrate=2;law.prefactor=1e14;law.barrier=.02;law.radius=1e-7;law.contactAngle=60;laws.push_back(law);}
    {
        Geometry rare(128,1,1,5e-6,20e-6,1e-5);for(auto& cell:rare.cells){paint(cell,-1,2,{{{0,1},.5}});merge(cell);}rare.rebuild();
        auto liquid=initialize(rare,spec,model.molarVolume,true);initializeWater(rare,liquid,model.molarVolume);ReactionMap faces;faces.rebuild(rare);const auto conserved=inventory(liquid,&rare);
        double water;auto solution=spec.conserved(liquid[0],water);auto law=laws[0];law.prefactor=1;
        const double interval=.01,hazard=.2;law.prefactor=hazard/(interval*rare.dissolution[0][2]*law.frequency(solution.saturationIndex[0],model.molarVolume[0],298.15));
        std::uint64_t expected=0;for(int p=0;p<rare.nx;p++)expected+=saltgeom::exponential(2024,p,2,0)<=hazard;
        require(expected>0&&expected<128,"Rare-event fixture has no statistical contrast");
        auto births=nucleateSurface(rare,faces,liquid,{law},model.molarVolume,interval,298.15,2024,[&](int)->saltflow::ChemicalWorker&{return workers.at(298.15,1);});
        require(births.count[0]==expected,"Sub-unit surface exposure silently suppressed rare nucleation events");
        for(int p=0;p<rare.nx;p++){double threshold=saltgeom::exponential(2024,p,2,0),remaining=hazard-(threshold<=hazard?threshold:0.);require(std::abs(rare.nucleationHazard[p][2]-remaining)<1e-12,"Nucleation failed to preserve its unused exposure history");}
        const auto after=inventory(liquid,&rare);for(int j=0;j<8;j++)require(std::abs(after[j]-conserved[j])<3e-12*std::max(1.L,std::abs(conserved[j])),"Rare nucleation lost aqueous or mineral inventory");
    }
    auto g=slab(3.75);ReactionMap map;map.rebuild(g);auto cells=initialize(g,spec,model.molarVolume,true);initializeWater(g,cells,model.molarVolume);const auto before=inventory(cells,&g);
    std::array<std::uint64_t,3> counts{};
    for(int k=0;k<6;k++){auto born=nucleateSurface(g,map,cells,laws,model.molarVolume,.01,298.15,2024,[&](int)->saltflow::ChemicalWorker&{return workers.at(298.15,1);});for(int m=0;m<3;m++)counts[m]+=born.count[m];map.rebuild(g);}
    if(!(counts[0]>0&&counts[1]>0))std::cerr<<"births H="<<counts[0]<<" S="<<counts[1]<<" hazard="<<g.nucleationHazard[0][2]<<','<<g.nucleationHazard[0][5]<<'\n';
    require(counts[0]>0&&counts[1]>0,"Supersaturated bare substrate failed to nucleate both products");{auto after=inventory(cells,&g);for(int j=0;j<8;j++)require(std::abs(before[j]-after[j])<3e-12*std::max(1.L,std::abs(before[j])),"Nucleation lost aqueous or mineral inventory");}
    require(g.fractions[0][0]>0&&g.fractions[0][1]>0,"New nuclei were not placed at their geometric owner");
    for(int p=0;p<g.nx*g.ny;p++)for(int m=0;m<3;m++)require(std::abs(g.fractions[p][m]-cells[p].mineral[m]*model.molarVolume[m])<1e-9,"Nucleus geometry and stock disagree");
    const std::string path="bin/tests/dual_nucleation_restart.bin";g.write(path);auto loaded=Geometry::read(path);std::remove(path.c_str());
    require(g.nucleationEvents==loaded.nucleationEvents&&g.nucleationHazard==loaded.nucleationHazard,"Restart changed nucleation state");
    {auto a=slab(3.75),b=a;auto left=initialize(a,spec,model.molarVolume,true),right=left;initializeWater(a,left,model.molarVolume);initializeWater(b,right,model.molarVolume);ReactionMap ma,mb;ma.rebuild(a);mb.rebuild(b);auto reversed=laws;std::reverse(reversed.begin(),reversed.end());
     nucleateSurface(a,ma,left,laws,model.molarVolume,.01,298.15,2024,[&](int)->saltflow::ChemicalWorker&{return workers.at(298.15,1);});
     nucleateSurface(b,mb,right,reversed,model.molarVolume,.01,298.15,2024,[&](int)->saltflow::ChemicalWorker&{return workers.at(298.15,1);});
     require(a.nucleationEvents==b.nucleationEvents&&a.nucleationHazard==b.nucleationHazard,"Nucleation depends on mineral-law list order");for(int p=0;p<a.nx*a.ny;p++)for(int m=0;m<3;m++)require(std::abs(a.fractions[p][m]-b.fractions[p][m])<1e-12,"Mineral-law ordering changed new-site geometry");}
    auto low=slab(3.75);map.rebuild(low);auto lowCells=initialize(low,spec,model.molarVolume);initializeWater(low,lowCells,model.molarVolume);auto noBirth=nucleateSurface(low,map,lowCells,laws,model.molarVolume,.01,298.15,2024,[&](int)->saltflow::ChemicalWorker&{return workers.at(298.15,1);});
    require(noBirth.count[0]==0&&noBirth.count[1]==0,"Undersaturated/saturated brine created nuclei");
    // Pre-existing product grains must not monopolise later precipitation.
    // Retreat the substrate, then require both products to appear and grow
    // at the newly exposed inner face while the old outer grains still exist.
    {
        auto moving=slab(5.75);const double h=moving.h();
        moving.circle(0,5.9*h,.75*h,.18*h,32);moving.circle(1,5.9*h,3.25*h,.18*h,32);moving.rebuild();
        auto liquid=initialize(moving,spec,model.molarVolume,true);initializeWater(moving,liquid,model.molarVolume);
        const auto original=moving;const auto conserved=inventory(liquid,&moving);auto target=moving.fractions;target[1][2]-=.25;
        moving.advance(target,false,1.2);
        for(int p=0;p<moving.nx*moving.ny;p++){
            MineralVector extent{};extent[2]=(original.fractions[p][2]-moving.fractions[p][2])/model.molarVolume[2];
            liquid[p]=Reactions::changed(liquid[p],extent);
        }
        remapWater(original,moving,liquid,&moving.pockets);ReactionMap current;current.rebuild(moving);
        auto innerAmount=[&](int mineral){double amount=0;for(int cell=0;cell<int(moving.cells.size());cell++)for(const auto& piece:moving.cells[cell])if(piece.phase==mineral){
            const double edge=5.3-moving.origin(cell).x;amount+=area(clip(piece.polygon,{1,0},edge));for(const auto& hole:piece.holes)amount-=area(clip(hole,{1,0},edge));
        }return amount;};
        require(innerAmount(0)<1e-12&&innerAmount(1)<1e-12,"Moving-face fixture contains an initial inner precipitate");
        for(int k=0;k<10;k++){
            nucleateSurface(moving,current,liquid,laws,model.molarVolume,.01,298.15,2024,[&](int)->saltflow::ChemicalWorker&{return workers.at(298.15,1);});
            current.rebuild(moving);
        }
        require(innerAmount(0)>1e-6&&innerAmount(1)>1e-6,"Precipitation stayed at the old outer grains after substrate retreat");
        const std::array<double,2> innerBefore{innerAmount(0),innerAmount(1)};const auto beforeGrowth=moving;
        applyMineralDebt(moving,liquid,model.molarVolume);
        auto localRates=rates;localRates.neutral[0]=.02;localRates.neutral[1]=.01;std::vector<RegionReaction> grown;
        std::array<double,2> deposited{};
        for(int r:current.active){auto result=reactRegion(r,moving,current,liquid,reactions,localRates,model.molarVolume,1e-4,options);for(int m=0;m<2;m++)deposited[m]+=result.stats.precipitated[m];grown.push_back(std::move(result));}
        auto growthBudgets=commitRegions(moving,liquid,grown,model.molarVolume);target=moving.fractions;
        for(int p=0;p<moving.nx*moving.ny;p++)for(int m=0;m<3;m++)target[p][m]=liquid[p].mineral[m]*model.molarVolume[m];
        moving.advance(target,false,.2,&growthBudgets);
        for(int p=0;p<moving.nx*moving.ny;p++)for(int m=0;m<3;m++)liquid[p].mineral[m]=moving.fractions[p][m]/model.molarVolume[m];
        remapWater(beforeGrowth,moving,liquid,&moving.pockets);
        require(deposited[0]>0&&deposited[1]>0,"New product interfaces failed to support subsequent growth");
        require(innerAmount(0)>innerBefore[0]&&innerAmount(1)>innerBefore[1],"Subsequent geometric growth was assigned only to old outer grains");
        const auto after=inventory(liquid,&moving);for(int j=0;j<8;j++)require(std::abs(after[j]-conserved[j])<3e-11*std::max(1.L,std::abs(conserved[j])),"Retreat, inner-surface nucleation and growth lost a conserved component");
    }
    // A closed water pocket continues reacting with every mineral type.
    for(int mineral=0;mineral<3;mineral++){
        std::cout<<"closed mineral="<<mineral<<std::endl;
        Geometry pocket(1,1,4,5e-6,20e-6,1e-5);
        auto box=[&](int from,int to,double lo,double hi){for(int k=0;k<int(pocket.cells.size());k++){auto o=pocket.origin(k);paint(pocket.cells[k],from,to,{{{1,0},hi-o.x},{{-1,0},o.x-lo},{{0,1},hi-o.y},{{0,-1},o.y-lo}});merge(pocket.cells[k]);}pocket.rebuild();};
        box(-1,mineral,.5,3.5);box(mineral,-1,1.5,2.5);
        auto mobile=initialize(pocket,spec,model.molarVolume);auto fresh=spec.inventory(Totals{{.1,0,0,0,.1,0,0}},1,{});double initialWater;double density=spec.conserved(fresh,initialWater).density;
        double factor=(1-pocket.fractions[0][mineral])*density/aqueousMass(fresh);for(double& a:fresh.aqueous)a*=factor;fresh.waterEquivalent*=factor;setWater(mobile[0],waterStock(fresh));
        initializeWater(pocket,mobile,model.molarVolume);map.rebuild(pocket);
        require(pocket.pockets.size()==1,"Independent closed water state missing");auto initial=inventory(mobile,&pocket);const auto trappedBefore=pocket.pockets[0].inventory;
        auto localRates=rates;localRates.neutral={};localRates.neutral[mineral]=.001;
        std::vector<RegionReaction> reactionsByWater;for(int r:map.active)reactionsByWater.push_back(reactRegion(r,pocket,map,mobile,reactions,localRates,model.molarVolume,1e-4,options));
        auto changes=commitRegions(pocket,mobile,reactionsByWater,model.molarVolume);auto after=inventory(mobile,&pocket);
        for(int j=0;j<8;j++)require(std::abs(after[j]-initial[j])<3e-12*std::max(1.L,std::abs(initial[j])),"Closed-water reaction lost total inventory");
        require(pocket.pockets[0].inventory[4]>trappedBefore[4],"Closed water stopped reacting with its mineral boundary");
        require(!changes.empty(),"Closed-water reaction supplied no interface budget");
    }
    std::cout<<"PASS dual surface chemistry: local liquid states, shared stock conservation, new-site growth after substrate retreat with pre-existing products, symmetric nucleation, closed-pocket reactions for all three minerals and restart\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
