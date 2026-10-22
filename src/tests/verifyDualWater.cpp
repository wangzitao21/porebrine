// SPDX-License-Identifier: GPL-3.0-or-later
#include "DualWater.hpp"
#include <iostream>
#include <cstdio>
using namespace saltdual;
static void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
static void rectangle(Geometry& g,int from,int to,double x0,double x1,double y0,double y1){
    for(int i=0;i<int(g.cells.size());i++){auto o=g.origin(i);paint(g.cells[i],from,to,{{{1,0},x1-o.x},{{-1,0},o.x-x0},{{0,1},y1-o.y},{{0,-1},o.y-y0}});merge(g.cells[i]);}g.rebuild();
}
static WaterStock total(const Geometry& g,const std::vector<saltchem::Inventory>& mobile){WaterStock t{};for(const auto& c:mobile)addWater(t,waterStock(c));for(const auto& p:g.pockets)addWater(t,p.inventory);return t;}
static void balance(const WaterStock& a,const WaterStock& b){for(int j=0;j<8;j++)require(std::abs(a[j]-b[j])<2e-12*std::max(1.,std::abs(a[j])),"Local water lost a conserved component");}
static WaterStock composition(double mass,double fraction){WaterStock s{};s[0]=mass*fraction/saltchem::componentMass[0];s[4]=mass*fraction/saltchem::componentMass[4];s[7]=mass*(1-2*fraction);return s;}
int main(){try{
    {
        Geometry graph;graph.waterRegions.resize(3);graph.waterRegions[0].bulk=graph.waterRegions[2].bulk=true;graph.waterLinks={{0,1,1,1},{1,2,1,3}};
        std::vector<std::map<int,double>> overlap(3);overlap[0][0]=overlap[2][2]=1;
        auto destinations=closingWaterRecipients(graph,1,overlap,0);require(destinations.size()==2&&destinations[0].first==0&&destinations[1].first==2&&std::abs(destinations[0].second-.25)<1e-13&&std::abs(destinations[1].second-.75)<1e-13,"Closing-water residue ignored hydraulic outlet resistance");
        Geometry chain;chain.waterRegions.resize(72);chain.waterRegions.front().bulk=chain.waterRegions.back().bulk=true;for(int i=0;i<71;i++)chain.waterLinks.push_back({i,i+1,1,1});
        overlap.assign(72,{});overlap.front()[0]=overlap.back()[71]=1;destinations=closingWaterRecipients(chain,36,overlap,0);
        require(destinations.size()==2&&std::abs(destinations[0].second-35./71)<1e-11&&std::abs(destinations[1].second-36./71)<1e-11,"Large closing-water network disagrees with the analytical resistance split");
    }
    // The analytical backward-Euler solution of two exchanging liquid pools.
    Geometry pair(1,1,1,1,1,0);pair.waterRegions={{0,true,true,.4,{},{}},{0,false,true,.6,{},{}}};pair.bulkRegion={0};pair.pocketRegion={-1,0};pair.pockets.resize(1);pair.pockets[0].parent=0;pair.pockets[0].fraction=.6;pair.waterInitialized=true;
    std::vector<saltchem::Inventory> mobile(1);setWater(mobile[0],composition(1,.05));pair.pockets[0].inventory=composition(2,.15);pair.waterLinks={{0,1,.02}};waterLedger(pair);
    const auto initial=total(pair,mobile);double density=2/(.4/1+.6/2),k=.5*.02*density;
    auto stats=exchangeWater(pair,mobile,1,.5);double difference=pair.pockets[0].inventory[0]*saltchem::componentMass[0]/2-mobile[0].aqueous[0]*saltchem::componentMass[0];
    require(std::abs(difference-.1/(1+k*1.5))<1e-13,"Finite water exchange disagrees with two-pool solution");balance(initial,total(pair,mobile));require(stats.links==1&&stats.residual<1e-12,"Water exchange matrix residual");
    {
        Geometry tinyPair;tinyPair.nx=2;tinyPair.ny=tinyPair.ratio=1;tinyPair.spacing=tinyPair.depth=1;tinyPair.epsilon=0;
        tinyPair.waterRegions={{0,true,true,1e-12,{},{}},{1,true,true,1,{},{}}};tinyPair.waterLinks={{0,1,.1}};
        std::vector<saltchem::Inventory> nodes(2);auto tiny=composition(1e-12,.1);tiny[0]-=.001/saltchem::componentMass[0];tiny[4]-=.001/saltchem::componentMass[4];tiny[7]+=.002;setWater(nodes[0],tiny);setWater(nodes[1],composition(1,.1));
        auto before=total(tinyPair,nodes);auto solved=exchangeWater(tinyPair,nodes,1,1,0);
        require(solved.links==1&&nodes[0].aqueous[0]>0&&nodes[0].aqueous[4]>0,"Closing-water exchange switched off or lost positivity");balance(before,total(tinyPair,nodes));
    }
    auto sealed=pair;auto closed=mobile;sealed.waterLinks.clear();auto before=total(sealed,closed);auto pocketBefore=sealed.pockets[0].inventory;exchangeWater(sealed,closed,1,1e6);require(sealed.pockets[0].inventory==pocketBefore,"Sealed water exchanged through solid");balance(before,total(sealed,closed));
    for(int mineral=0;mineral<3;mineral++){
        // Opposite faces retreat into the same grain while a solid wall
        // still separates their water. A tracer must remain on its own side.
        Geometry g(1,1,1,1,1,0);rectangle(g,-1,mineral,.2,.8,0,1);std::vector<saltchem::Inventory> stock(1);
        setWater(stock[0],composition(400,.05));saltchem::MineralVector volumes{};initializeWater(g,stock,volumes);
        std::vector<SurfaceBudget> budgets;for(int r=0;r<int(g.waterRegions.size());r++){auto fluid=regionStock(g,stock,r);if(g.waterRegions[r].centre.x<.5){fluid[2]=1;fluid[4]+=2;}setRegionStock(g,stock,r,fluid);budgets.push_back({0,mineral,g.waterRegions[r].centre.x<.5?-.02:-.04,{},r});}
        waterLedger(g);const auto original=g;const auto conserved=total(g,stock);for(int k=0;k<int(g.wet.size());k++)budgets[g.wet[k].liquidRegion].surfaces.push_back(k);
        auto desired=g.fractions;desired[0][mineral]-=.06;g.advance(desired,false,.2,&budgets);remapWater(original,g,stock,&g.pockets);
        require(g.waterRegions.size()==2&&g.waterLinks.empty(),"Opposing retreat removed the separating wall");
        for(int r=0;r<int(g.waterRegions.size());r++){double tracer=regionStock(g,stock,r)[2];require(std::abs(tracer-(g.waterRegions[r].centre.x<.5?1.:0.))<1e-12,"Newly dissolved space mixed water through the remaining solid wall");}
        balance(conserved,total(g,stock));
    }
    {
        auto geometry=[](double shift){Geometry g(1,1,1,1,1,0);g.cells[0]={{2,square(),{}}};
            for(int y=0;y<20;y++)for(int x=0;x<20;x++){double a=.02+.047*x+shift*(x%2?1:-1),b=.02+.047*y;Poly p{{a,b},{a+.018,b},{a+.018,b+.018},{a,b+.018}};g.cells[0].push_back({-1,p,{}});std::reverse(p.begin(),p.end());g.cells[0][0].holes.push_back(std::move(p));}g.rebuild();return g;};
        auto old=geometry(0),next=geometry(.001);makePockets(old);old.waterInitialized=true;
        for(std::size_t i=0;i<old.pockets.size();i++){old.pockets[i].inventory=composition(1000*old.pockets[i].fraction,.05);old.pockets[i].inventory[2]=(i+1)*old.pockets[i].fraction;}waterLedger(old);
        std::vector<saltchem::Inventory> cells(1);const auto before=total(old,cells);remapWater(old,next,cells);
        require(next.pockets.size()==old.pockets.size(),"Dense water remap changed the number of separated pools");
        for(std::size_t i=0;i<old.pockets.size();i++)require(next.pockets[i].inventory==old.pockets[i].inventory,"Dense water remap transferred inventory through a separating mineral wall");
        balance(before,total(next,cells));
    }
    // Geometry alone determines the throat, independent of mineral identity.
    for(int m=0;m<3;m++){
        Geometry g(1,1,4,4e-6,20e-6,1e-5);rectangle(g,-1,m,1.8,2.2,0,4);rectangle(g,m,-1,1.8,2.2,1.96,2.04);
        require(g.globalLiquidComponents==1&&g.waterRegions.size()>1&&!g.waterLinks.empty(),"Narrow water throat was homogenized or disconnected");
        std::vector<saltchem::Inventory> cells(1);double liquid=1-g.fractions[0][m];setWater(cells[0],composition(1000*liquid,.05));saltchem::MineralVector volumes{};
        initializeWater(g,cells,volumes);require(!g.pockets.empty(),"Narrow throat retained no local water");
        for(auto& p:g.pockets)p.inventory=composition(waterMass(p.inventory),.1);waterLedger(g);auto originalTotal=total(g,cells);
        double conductance=0;for(const auto& link:g.waterLinks)conductance+=link.conductance;
        auto old=g;rectangle(g,m,-1,1.8,2.2,1.8,2.2);remapWater(old,g,cells);balance(originalTotal,total(g,cells));
        require(g.pockets.empty(),"Opening a broad water connection failed to merge local inventories");
        require(conductance>0,"Narrow throat has zero exchange conductance");
    }
    {
        Geometry g(2,1,1,1,1,0);rectangle(g,-1,0,0,1,0,.5);rectangle(g,0,-1,.75,1,.1,.2);
        std::vector<saltchem::Inventory> liquid(2);for(int p=0;p<2;p++)setWater(liquid[p],composition(1000*(1-g.fractions[p][0]),.05));saltchem::MineralVector volumes{};initializeWater(g,liquid,volumes);
        require(g.pockets.size()==1,"Closing-pore fixture is missing its retained water");g.pockets[0].inventory=composition(1000*g.pockets[0].fraction,.15);waterLedger(g);
        const auto original=g;const auto conserved=total(g,liquid);const auto untouched=waterStock(liquid[0]);auto expelled=g.pockets[0].inventory;
        double fraction=1e-10/g.pockets[0].fraction;for(double& v:g.pockets[0].inventory)v*=fraction;addWater(expelled,g.pockets[0].inventory,-1);auto receiver=waterStock(liquid[1]);addWater(receiver,expelled);setWater(liquid[1],receiver);
        rectangle(g,-1,0,.75,1,.1,.2);remapWater(original,g,liquid,&g.pockets);
        require(g.pockets.empty()&&waterStock(liquid[0])==untouched,"Closing-water residue crossed solid into the unrelated parent reservoir");balance(conserved,total(g,liquid));
    }
    Geometry pocket(1,1,4,4e-6,20e-6,1e-5);rectangle(pocket,-1,0,.5,3.5,.5,3.5);rectangle(pocket,0,-1,1.5,2.5,1.5,2.5);
    std::vector<saltchem::Inventory> cells(1);setWater(cells[0],composition(1000*(1-pocket.fractions[0][0]),.1));saltchem::MineralVector volumes{};initializeWater(pocket,cells,volumes);
    require(pocket.pockets.size()==1&&pocket.waterLinks.empty(),"Closed cavity topology");
    auto initialTotal=total(pocket,cells);pocket.chemicalVolumeDebt[0]=1e-11;const std::string path="bin/tests/dual_water_restart.bin";pocket.write(path);auto restored=Geometry::read(path);std::remove(path.c_str());
    require(restored.waterInitialized&&restored.pockets.size()==1&&restored.pockets[0].inventory==pocket.pockets[0].inventory,"Local water checkpoint changed composition");balance(initialTotal,total(restored,cells));
    require(restored.chemicalVolumeDebt==pocket.chemicalVolumeDebt,"Checkpoint lost deferred chemical mineral volume");
    {auto original=pocket,closing=pocket;auto liquid=cells;for(double& v:closing.pockets[0].inventory)v*=1e-10;rectangle(closing,-1,0,1.5,2.5,1.5,2.5);bool rejected=false;try{remapWater(original,closing,liquid,&closing.pockets);}catch(const StepFailure&){rejected=true;}require(rejected,"A resolved isolated pore was drained through solid during closure");}
    auto old=pocket;rectangle(pocket,0,-1,0,2,.95,2.05);remapWater(old,pocket,cells);balance(initialTotal,total(pocket,cells));
    require(pocket.pockets.empty(),"Reopened cavity retained an isolated inventory");
    std::cout<<"PASS local water: finite throat exchange, analytical balance, sealed water, symmetric mineral handling, conservative reopening and restart\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
