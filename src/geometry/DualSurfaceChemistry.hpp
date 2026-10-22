// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "DualWater.hpp"
#include <tuple>

namespace saltdual {
using Areas=std::array<double,3>;
struct ReactionPatch {
    int owner=-1,mineral=-1;
    double dissolution=0,growth=0,stockShare=0;
    std::vector<int> dissolveFaces,growFaces;
};
struct ReactionMap {
    std::vector<int> waterReceiver,receiver,active;
    std::vector<Areas> dissolution,growth;
    std::vector<std::array<Point,3>> dissolveNormal,growNormal;
    std::vector<double> pore,closed;
    std::vector<std::vector<ReactionPatch>> patches;
    void rebuild(const Geometry& g,bool unblocked=false,double minimumLiquid=.5){
        const int n=g.nx*g.ny;
        waterReceiver.assign(n,-1);receiver.assign(n,-1);active.clear();patches.assign(g.waterRegions.size(),{});
        dissolution.assign(n,{});growth.assign(n,{});dissolveNormal.assign(n,{});growNormal.assign(n,{});pore.resize(n);closed=g.closedFraction;
        std::vector<std::map<std::pair<int,int>,ReactionPatch>> entries(g.waterRegions.size());
        auto add=[&](int r,int owner,int m,int face,double area,bool growing,Point normal){
            auto& patch=entries[r][{owner,m}];patch.owner=owner;patch.mineral=m;
            if(growing){patch.growth+=area;patch.growFaces.push_back(face);growth[owner][m]+=area;growNormal[owner][m]=growNormal[owner][m]+normal*area;}
            else{patch.dissolution+=area;patch.dissolveFaces.push_back(face);dissolution[owner][m]+=area;dissolveNormal[owner][m]=dissolveNormal[owner][m]+normal*area;}
            receiver[owner]=g.waterRegions[r].parent;
        };
        std::map<std::pair<int,int>,double> pieceContact;
        for(const auto& s:g.wet)if(resolvedWater(g,s.liquidRegion)&&s.solidPiece>=0)pieceContact[{s.solid,s.solidPiece}]+=length(s.b-s.a);
        for(int k=0;k<int(g.wet.size());k++){const auto& s=g.wet[k];if(!resolvedWater(g,s.liquidRegion))continue;double a=length(s.b-s.a)*g.h()*g.depth;
            add(s.liquidRegion,g.parent(s.solid),s.mineral,k,a,false,s.normal);add(s.liquidRegion,g.parent(s.water),s.mineral,k,a,true,s.normal);}
        // Only physically exposed fine pieces supply dissolution stock.
        // A small wet crystal cannot borrow a buried crystal's inventory
        // merely because both share the same coarse cell and mineral label.
        for(const auto& s:g.wet)if(resolvedWater(g,s.liquidRegion)&&s.solidPiece>=0){int p=g.parent(s.solid);double phase=g.fractions[p][s.mineral];if(phase>0)entries[s.liquidRegion][{p,s.mineral}].stockShare+=pieceArea(g.cells[s.solid][s.solidPiece])*(1-g.epsilon)/(g.ratio*g.ratio*phase)*length(s.b-s.a)/pieceContact.at({s.solid,s.solidPiece});}
        if(unblocked)for(int k=0;k<int(g.buried.size());k++){const auto& s=g.buried[k];int p=g.parent(s.solid),r=g.bulkRegion[p];
            if(r>=0)add(r,p,s.mineral,-k-1,length(s.b-s.a)*g.h()*g.depth,false,s.normal);}
        for(int p=0;p<n;p++){
            pore[p]=std::max(0.,1-std::accumulate(g.fractions[p].begin(),g.fractions[p].end(),0.)-closed[p]);
            if(g.bulkRegion[p]>=0)waterReceiver[p]=p;
            for(int m=0;m<3;m++){if(dissolution[p][m]>0)dissolveNormal[p][m]=dissolveNormal[p][m]*(1/dissolution[p][m]);if(growth[p][m]>0)growNormal[p][m]=growNormal[p][m]*(1/growth[p][m]);}
        }
        for(int r=0;r<int(entries.size());r++)if(!entries[r].empty()){
            for(auto& entry:entries[r]){auto patch=std::move(entry.second);if(unblocked&&patch.dissolution>0)patch.stockShare=patch.dissolution/dissolution[patch.owner][patch.mineral];patches[r].push_back(std::move(patch));}
            active.push_back(r);
        }
        (void)minimumLiquid;
    }
};

inline std::vector<double> boundedPartition(double amount,const std::vector<double>& weights,const std::vector<double>& capacity){
    std::vector<double> result(weights.size());double remaining=amount;
    const double tolerance=128*std::numeric_limits<double>::epsilon()*std::max(1.,amount);
    for(std::size_t pass=0;pass<=weights.size()&&remaining>0;pass++){
        double sum=0;int last=-1;for(std::size_t k=0;k<weights.size();k++)if(capacity[k]>result[k]&&weights[k]>0){sum+=weights[k];last=int(k);}
        if(last<0){if(remaining<=tolerance)break;throw StepFailure("Reaction extent exceeds accessible geometric stock/pore volume");}
        const double before=remaining;bool limited=false;
        for(std::size_t k=0;k<weights.size();k++)if(weights[k]>0&&capacity[k]>result[k]){
            const double room=capacity[k]-result[k];if(before*weights[k]/sum>=room){result[k]+=room;remaining-=room;limited=true;}
        }
        if(limited){remaining=std::max(0.,remaining);continue;}
        for(std::size_t k=0;k<weights.size();k++)if(int(k)!=last&&weights[k]>0&&capacity[k]>result[k]){double a=before*weights[k]/sum;result[k]+=a;remaining-=a;}
        result[last]+=remaining;remaining=0;
    }
    // Restore only floating-point summation residue, never clip an extent.
    double residue=amount-std::accumulate(result.begin(),result.end(),0.);
    if(residue!=0){int best=-1;for(std::size_t k=0;k<weights.size();k++)if(weights[k]>0&&result[k]+residue>=-tolerance&&result[k]+residue<=capacity[k]+tolerance)best=int(k);if(best<0)throw StepFailure("Conservative extent partition failed");result[best]+=residue;}
    return result;
}

struct RegionReaction {
    int region=-1;
    WaterStock water{};
    saltchem::IntegrationStats stats;
    double exchangeResidual=0,nonlinearResidualRatio=0;
    double nonlinearStateErrorRatio=0;
    std::vector<SurfaceBudget> budgets;
};
struct RegionProblem {saltchem::Inventory initial;saltchem::RateParameters rates;};
inline RegionProblem regionProblem(int r,const Geometry& g,const ReactionMap& map,const std::vector<saltchem::Inventory>& cells,
    const saltchem::RateParameters& base,const saltchem::MineralVector& volumes){
    saltchem::Inventory initial;setWater(initial,regionStock(g,cells,r));initial.mineral={};auto rates=base;
    rates.area={};rates.geometry.growthArea={};rates.geometry.lockedMoles={};rates.geometry.surfaceExponent={};rates.geometry.curvatureLogOmega={};
    rates.geometry.separateMineralAreas=true;rates.geometry.enabled=true;rates.geometry.inertFraction=0;rates.geometry.initialPorosity=NAN;
    rates.geometry.bulkVolume=g.waterRegions[r].fraction+(g.waterRegions[r].bulk?g.epsilon:0);
    const double liquidVolume=rates.geometry.bulkVolume;
    for(const auto& patch:map.patches[r]){
        int m=patch.mineral;double stock=cells[patch.owner].mineral[m]*patch.stockShare;
        initial.mineral[m]+=stock;rates.area[m]+=patch.dissolution/g.coarseVolume();rates.geometry.growthArea[m]+=patch.growth/g.coarseVolume();
        rates.geometry.bulkVolume+=stock*volumes[m];
    }
    for(int m=0;m<3;m++)rates.referenceMoles[m]=std::max(initial.mineral[m],1e-30);
    rates.geometry.initialPorosity=liquidVolume/rates.geometry.bulkVolume;
    // The coarse-cell pore floor is a regularisation of the main flow field.
    // A real thin local water region has its own smaller positive volume.
    rates.geometry.epsilonMin=std::min(base.geometry.epsilonMin,.1*liquidVolume/rates.geometry.bulkVolume);
    return {initial,rates};
}
inline RegionReaction regionResult(int r,const Geometry& g,const ReactionMap& map,const std::vector<saltchem::Inventory>& cells,
    const saltchem::MineralVector& volumes,const saltchem::Inventory& initial,const saltchem::Inventory& reacted,const saltchem::IntegrationStats& stats={},const saltchem::MineralVector* integrated=nullptr){
    RegionReaction result;result.region=r;result.water=waterStock(reacted);result.stats=stats;
    for(int m=0;m<3;m++){
        double extent=integrated?(*integrated)[m]:initial.mineral[m]-reacted.mineral[m];if(extent==0)continue;
        std::vector<double> weights,capacity;std::vector<const ReactionPatch*> selected;
        for(const auto& patch:map.patches[r])if(patch.mineral==m&&(extent>0?patch.dissolution:patch.growth)>0){
            selected.push_back(&patch);weights.push_back(extent>0?patch.dissolution:patch.growth);
            capacity.push_back(extent>0?cells[patch.owner].mineral[m]*patch.stockShare:std::abs(extent));
        }
        auto shares=boundedPartition(std::abs(extent),weights,capacity);
        for(std::size_t k=0;k<shares.size();k++)if(shares[k]>0){const auto& patch=*selected[k];result.budgets.push_back({patch.owner,m,(extent>0?-1:1)*shares[k]*volumes[m],extent>0?patch.dissolveFaces:patch.growFaces,r});}
    }
    return result;
}
inline RegionReaction reactRegion(int r,const Geometry& g,const ReactionMap& map,const std::vector<saltchem::Inventory>& cells,
    const saltchem::Reactions& chemistry,const saltchem::RateParameters& base,const saltchem::MineralVector& volumes,
    double dt,const saltchem::IntegrationOptions& options){
    auto problem=regionProblem(r,g,map,cells,base,volumes);auto reacted=problem.initial;saltchem::IntegrationStats stats;
    bool active=false;for(int m=0;m<3;m++)active|=base.neutral[m]>0&&(problem.rates.area[m]>0||problem.rates.geometry.growthArea[m]>0);
    if(active&&reacted.waterEquivalent>1e-30)stats=chemistry.advance(reacted,dt,problem.rates,options);
    saltchem::MineralVector extent{};for(int m=0;m<3;m++)extent[m]=stats.dissolved[m]-stats.precipitated[m];
    return regionResult(r,g,map,cells,volumes,problem.initial,reacted,stats,&extent);
}
inline std::vector<SurfaceBudget> commitRegions(Geometry& g,std::vector<saltchem::Inventory>& cells,const std::vector<RegionReaction>& results,const saltchem::MineralVector& volumes){
    std::vector<SurfaceBudget> budgets;std::vector<std::array<long double,3>> changes(cells.size());
    for(const auto& result:results)if(result.region>=0){setRegionStock(g,cells,result.region,result.water);for(const auto& b:result.budgets){changes[b.parent][b.mineral]+=static_cast<long double>(b.change)/volumes[b.mineral];budgets.push_back(b);}}
    for(int p=0;p<int(cells.size());p++)for(int m=0;m<3;m++){
        double value=double(static_cast<long double>(cells[p].mineral[m])+changes[p][m]);
        if(value< -1e-11*std::max(1.,cells[p].mineral[m]))throw StepFailure("Shared mineral stock exceeded by local water reactions");
        cells[p].mineral[m]=std::max(0.,value);
    }
    waterLedger(g);return budgets;
}

struct NucleationResult {std::array<std::uint64_t,3> count{};std::vector<Point> volumeVectors;};
template<class Worker>
NucleationResult nucleateSurface(Geometry& g,const ReactionMap& map,std::vector<saltchem::Inventory>& cells,
    const std::vector<saltgeom::NucleationLaw>& laws,const saltchem::MineralVector& volume,double dt,double temperature,
    std::uint64_t seed,Worker worker){
    NucleationResult result;result.volumeVectors.resize(cells.size());if(laws.empty())return result;
    struct Candidate {int parent,law;double time,total,threshold;std::uint64_t key;std::vector<std::pair<int,double>> faces;};
    std::vector<Candidate> candidates;std::vector<double> saturation(g.waterRegions.size()*3,NAN);
    std::vector<std::array<std::vector<int>,3>> surfaceByParent(g.nx*g.ny);
    for(int i=0;i<int(g.wet.size());i++)if(resolvedWater(g,g.wet[i].liquidRegion))surfaceByParent[g.parent(g.wet[i].water)][g.wet[i].mineral].push_back(i);
    for(int p=0;p<g.nx*g.ny;p++)for(int k=0;k<int(laws.size());k++){
        const auto& law=laws[k];Candidate event{p,k,0,0,0,0,{}};int pair=3*law.product+law.substrate;
        for(int i:surfaceByParent[p][law.substrate]){const auto& s=g.wet[i];
            double& si=saturation[3*s.liquidRegion+law.product];
            if(std::isnan(si)){auto stock=regionStock(g,cells,s.liquidRegion);if(stock[7]<=1e-30){si=-99;continue;}saltchem::Inventory inventory;setWater(inventory,stock);double water;si=worker(p).speciation.conserved(inventory,water,nullptr,false).saturationIndex[law.product];}
            double frequency=law.frequency(si,volume[law.product],temperature)*length(s.b-s.a)*g.h()*g.depth;
            if(frequency>0){event.faces.push_back({i,frequency});event.total+=frequency;}
        }
        double& hazard=g.nucleationHazard[p][pair];double previous=hazard;hazard+=event.total*dt;
        // Independent exponential thresholds retain rare events on a moving
        // surface. A fixed threshold of one suppresses every cell whose finite
        // exposure ends below that threshold, even with a large total hazard.
        // The saved counter and seed reproduce the threshold after a restart.
        event.threshold=saltgeom::exponential(seed,std::uint64_t(p),pair,g.nucleationEvents[p][pair]);
        if(hazard<event.threshold||event.total==0)continue;
        event.time=std::max(0.,event.threshold-previous)/event.total;event.key=saltgeom::mix(seed^std::uint64_t(p)*0x9e3779b97f4a7c15ULL^std::uint64_t(pair)*0xbf58476d1ce4e5b9ULL^g.nucleationEvents[p][pair]);candidates.push_back(std::move(event));
    }
    if(candidates.empty())return result;
    std::sort(candidates.begin(),candidates.end(),[](const auto& a,const auto& b){return std::tie(a.time,a.key)<std::tie(b.time,b.key);});
    const auto original=g;std::set<int> used;bool changed=false;
    for(const auto& event:candidates){
        const auto& law=laws[event.law];int p=event.parent,pair=3*law.product+law.substrate;
        double selection=double(event.key>>11)*0x1.0p-53*event.total;int selected=event.faces.back().first;
        for(auto [face,weight]:event.faces){if(selection<=weight){selected=face;break;}selection-=weight;}
        const auto& face=g.wet[selected];if(used.count(face.water))continue;int r=face.liquidRegion;
        saltchem::Inventory stock;setWater(stock,regionStock(g,cells,r));stock.mineral=cells[p].mineral;double water;
        auto& chemical=worker(p);double si=chemical.speciation.conserved(stock,water,nullptr,false).saturationIndex[law.product];if(si<=0)continue;
        double radius=law.embryoRadius(si,volume[law.product],temperature)/g.h();if(!(radius>0&&radius<=.2))throw std::invalid_argument("Dual-grid embryo radius must not exceed 0.2 h");
        double position=.15+.7*double(saltgeom::mix(event.key)>>11)*0x1.0p-53;
        auto centre=face.a+(face.b-face.a)*position-face.normal*(radius*std::cos(law.contactAngle*saltgeom::pi/180))-g.origin(face.water);
        std::vector<std::pair<Point,double>> planes;for(int k=0;k<16;k++){double a=2*saltgeom::pi*(k+.5)/16;Point n{std::cos(a),std::sin(a)};planes.push_back({n,radius*std::cos(saltgeom::pi/16)+dot(n,centre)});}
        Cell localWater{g.cells[face.water][face.liquidPiece]};double areaAdded=paint(localWater,-1,law.product,planes);if(areaAdded<1e-12)continue;
        auto candidate=g.cells[face.water];candidate.erase(candidate.begin()+face.liquidPiece);candidate.insert(candidate.end(),localWater.begin(),localWater.end());
        double fraction=areaAdded*(1-g.epsilon)/(g.ratio*g.ratio),moles=fraction/volume[law.product];
        saltchem::MineralVector extent{};extent[law.product]=-moles;auto after=saltchem::Reactions::changed(stock,extent);
        if(!(after.waterEquivalent>0)||std::any_of(after.aqueous.begin(),after.aqueous.end(),[](double a){return !std::isfinite(a)||a<0;}))continue;
        if(chemical.speciation.conserved(after,water,nullptr,false).saturationIndex[law.product]<0)continue;
        merge(candidate);g.cells[face.water]=std::move(candidate);cells[p].mineral=after.mineral;setRegionStock(g,cells,r,waterStock(after));
        result.volumeVectors[p]=result.volumeVectors[p]-face.normal*fraction;g.nucleationHazard[p][pair]-=event.threshold;g.nucleationEvents[p][pair]++;result.count[law.product]++;used.insert(face.water);changed=true;
    }
    if(changed){
        g.createdWater.clear();g.rebuild();
        for(int m=0;m<3;m++){long double difference=0;for(int p=0;p<g.nx*g.ny;p++){difference+=static_cast<long double>(cells[p].mineral[m])*volume[m]-g.fractions[p][m];cells[p].mineral[m]=g.fractions[p][m]/volume[m];}g.chemicalVolumeDebt[m]+=double(difference);}
        remapWater(original,g,cells,&g.pockets);
    }else waterLedger(g);
    (void)map;return result;
}
}
