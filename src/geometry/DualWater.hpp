// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "DualMineralGeometry.hpp"
#include "BrineReactions.hpp"

namespace saltdual {
using WaterStock=std::array<double,8>;
inline bool resolvedWater(const Geometry& g,int r){return g.waterRegions[r].fraction+(g.waterRegions[r].bulk?g.epsilon:0)>1e-10;}
inline WaterStock waterStock(const saltchem::Inventory& s){WaterStock a;for(int j=0;j<7;j++)a[j]=s.aqueous[j];a[7]=s.waterEquivalent;return a;}
inline void setWater(saltchem::Inventory& s,const WaterStock& a){for(int j=0;j<7;j++)s.aqueous[j]=a[j];s.waterEquivalent=a[7];}
inline void addWater(WaterStock& a,const WaterStock& b,double factor=1){for(int j=0;j<8;j++)a[j]+=factor*b[j];}
inline double waterMass(const WaterStock& a){double m=a[7];for(int j=0;j<7;j++)m+=saltchem::componentMass[j]*a[j];return m;}
inline double massFactor(int j){return j==7?1:saltchem::componentMass[j];}
// Sub-resolution interface motion remains in the chemical mineral ledger.
// It is reintroduced before the next solve and included in the next volume
// budget, rather than repeatedly forcing subatomic polygon operations.
inline void applyMineralDebt(const Geometry& g,std::vector<saltchem::Inventory>& cells,const saltchem::MineralVector& volumes){
    for(int m=0;m<3;m++)if(g.chemicalVolumeDebt[m]!=0){
        long double total=0;int largest=0;for(int p=0;p<int(cells.size());p++){total+=cells[p].mineral[m];if(cells[p].mineral[m]>cells[largest].mineral[m])largest=p;}
        const long double target=total+static_cast<long double>(g.chemicalVolumeDebt[m])/volumes[m];
        if(target< -1e-12)throw std::runtime_error("Negative deferred chemical mineral inventory");
        if(total==0){cells[largest].mineral[m]=std::max(0.,double(target));continue;}
        const long double desired=std::max(0.L,target);const double factor=double(desired/total);long double after=0;
        for(auto& cell:cells){cell.mineral[m]*=factor;after+=cell.mineral[m];}cells[largest].mineral[m]+=double(desired-after);
    }
}
inline WaterStock mineralDebtComponents(const Geometry& g,const saltchem::MineralVector& volumes){
    WaterStock result{};for(int m=0;m<3;m++){double n=g.chemicalVolumeDebt[m]/volumes[m];for(int j=0;j<7;j++)result[j]+=saltchem::stoichiometry[m][j]*n;result[7]+=saltchem::hydration[m]*saltchem::waterMolarMass*n;}return result;
}
inline WaterStock regionStock(const Geometry& g,const std::vector<saltchem::Inventory>& mobile,int r){
    return g.waterRegions[r].bulk?waterStock(mobile[g.waterRegions[r].parent]):g.pockets.at(g.pocketRegion.at(r)).inventory;
}
inline void setRegionStock(Geometry& g,std::vector<saltchem::Inventory>& mobile,int r,const WaterStock& stock){
    if(g.waterRegions[r].bulk)setWater(mobile[g.waterRegions[r].parent],stock);else g.pockets.at(g.pocketRegion.at(r)).inventory=stock;
}
inline void waterLedger(Geometry& g){
    g.storedClosedFraction.assign(g.nx*g.ny,0);g.closedInventory.assign(g.nx*g.ny,{});
    for(const auto& p:g.pockets){g.storedClosedFraction[p.parent]+=p.fraction;addWater(g.closedInventory[p.parent],p.inventory);}
}
inline void makePockets(Geometry& g){
    g.pockets.clear();g.pocketRegion.assign(g.waterRegions.size(),-1);
    for(int r=0;r<int(g.waterRegions.size());r++)if(!g.waterRegions[r].bulk){
        const auto& region=g.waterRegions[r];WaterPocket p;p.parent=region.parent;p.fraction=region.fraction;
        g.pocketRegion[r]=g.pockets.size();g.pockets.push_back(std::move(p));
    }
}
inline void initializeWater(Geometry& g,std::vector<saltchem::Inventory>& mobile,const saltchem::MineralVector& volumes){
    if(g.waterInitialized)return;
    std::vector<WaterStock> total(mobile.size());std::vector<double> pore(mobile.size(),g.epsilon);
    for(int p=0;p<int(mobile.size());p++){total[p]=waterStock(mobile[p]);addWater(total[p],g.closedInventory[p]);}
    for(const auto& r:g.waterRegions)pore[r.parent]+=r.fraction;
    makePockets(g);
    for(auto& p:g.pockets){p.inventory=total[p.parent];for(double& a:p.inventory)a*=p.fraction/pore[p.parent];}
    waterLedger(g);
    for(int p=0;p<int(mobile.size());p++){auto a=total[p];addWater(a,g.closedInventory[p],-1);setWater(mobile[p],a);}
    g.waterInitialized=true;
    (void)volumes;
}
inline bool sameRing(const Poly& a,const Poly& b){if(a.size()!=b.size())return false;for(std::size_t i=0;i<a.size();i++)if(a[i].x!=b[i].x||a[i].y!=b[i].y)return false;return true;}
inline bool samePiece(const Piece& a,const Piece& b){if(a.phase!=b.phase||a.holes.size()!=b.holes.size()||!sameRing(a.polygon,b.polygon))return false;for(std::size_t i=0;i<a.holes.size();i++)if(!sameRing(a.holes[i],b.holes[i]))return false;return true;}

inline std::vector<std::pair<int,double>> closingWaterRecipients(const Geometry& old,int source,
    const std::vector<std::map<int,double>>& overlap,double mobileFloor){
    struct Edge{int to;double conductance;};std::vector<std::vector<Edge>> links(old.waterRegions.size());
    for(const auto& e:old.waterLinks){double k=e.hydraulicConductance>0?e.hydraulicConductance:e.conductance;if(k>0){links[e.a].push_back({e.b,k});links[e.b].push_back({e.a,k});}}
    auto solve=[&](bool bulkOnly){
        std::vector<std::pair<int,double>> result;std::vector<int> index(links.size(),-1),nodes;std::vector<unsigned char> outlet(links.size());
        index[source]=0;nodes.push_back(source);int outlets=0;
        for(std::size_t head=0;head<nodes.size();head++)for(const auto& edge:links[nodes[head]]){
            int r=edge.to;bool receives=!overlap[r].empty()||(old.waterRegions[r].bulk&&mobileFloor>0);
            if(r!=source&&receives&&(!bulkOnly||old.waterRegions[r].bulk)){if(!outlet[r]){outlet[r]=1;outlets++;}}
            else if(index[r]<0){index[r]=nodes.size();nodes.push_back(r);}
        }
        if(!outlets)return result;
        const int n=nodes.size();double scale=0;for(int r:nodes)for(const auto& e:links[r])scale=std::max(scale,e.conductance);
        if(!(scale>0))return result;
        struct Internal{int a,b;long double k;};std::vector<Internal> internal;std::vector<long double> diagonal(n),pressure(n),rhs(n);rhs[0]=1;
        for(int i=0;i<n;i++)for(const auto& edge:links[nodes[i]]){long double k=edge.conductance/scale;diagonal[i]+=k;int j=index[edge.to];if(j>=0&&i<j)internal.push_back({i,j,k});}
        if(n<=64){
            std::vector<long double> factor(n*n);for(int i=0;i<n;i++)factor[i*n+i]=diagonal[i];for(const auto& e:internal){factor[e.a*n+e.b]-=e.k;factor[e.b*n+e.a]-=e.k;}
            for(int i=0;i<n;i++)for(int j=0;j<=i;j++){long double value=factor[i*n+j];for(int k=0;k<j;k++)value-=factor[i*n+k]*factor[j*n+k];if(i==j){if(!(value>0))return result;factor[i*n+j]=std::sqrt(value);}else factor[i*n+j]=value/factor[j*n+j];}
            for(int i=0;i<n;i++){long double value=rhs[i];for(int j=0;j<i;j++)value-=factor[i*n+j]*pressure[j];pressure[i]=value/factor[i*n+i];}
            for(int i=n-1;i>=0;i--){long double value=pressure[i];for(int j=i+1;j<n;j++)value-=factor[j*n+i]*pressure[j];pressure[i]=value/factor[i*n+i];}
        }else{
            auto multiply=[&](const auto& x,auto& y){for(int i=0;i<n;i++)y[i]=diagonal[i]*x[i];for(const auto&e:internal){y[e.a]-=e.k*x[e.b];y[e.b]-=e.k*x[e.a];}};
            auto residual=rhs;std::vector<long double> z(n),direction(n),product(n);long double rz=0;
            for(int i=0;i<n;i++){if(!(diagonal[i]>0))return result;z[i]=residual[i]/diagonal[i];direction[i]=z[i];rz+=residual[i]*z[i];}
            bool converged=false;
            for(int iteration=0;iteration<std::max(128,8*n);iteration++){
                long double error=0;for(auto v:residual)error=std::max(error,std::abs(v));if(error<1e-13L){converged=true;break;}
                multiply(direction,product);long double pap=0;for(int i=0;i<n;i++)pap+=direction[i]*product[i];if(!(pap>0&&rz>0))return result;
                long double alpha=rz/pap,next=0;for(int i=0;i<n;i++){pressure[i]+=alpha*direction[i];residual[i]-=alpha*product[i];z[i]=residual[i]/diagonal[i];next+=residual[i]*z[i];}
                long double beta=next/rz;for(int i=0;i<n;i++)direction[i]=z[i]+beta*direction[i];rz=next;
            }
            if(!converged)return result;
        }
        std::map<int,long double> discharge;long double sum=0;
        for(int i=0;i<n;i++)for(const auto& edge:links[nodes[i]])if(outlet[edge.to]){long double q=edge.conductance/scale*pressure[i];if(q< -1e-12L)return result;q=std::max(0.L,q);discharge[edge.to]+=q;sum+=q;}
        if(sum>0)for(auto [r,q]:discharge)if(q>0)result.push_back({r,double(q/sum)});
        return result;
    };
    auto result=solve(true);return result.empty()?solve(false):result;
}

// Remap the complete old water inventories using actual polygon overlap.
// Normalisation keeps solute and solvent when a reacting boundary moves;
// splitting, merging and reopening all use the same conservative operation.
inline void remapWater(const Geometry& old,Geometry& g,std::vector<saltchem::Inventory>& mobile,const std::vector<WaterPocket>* reactedPockets=nullptr){
    if(!old.waterInitialized)throw std::logic_error("Initialize local water before moving geometry");
    using Ref=std::pair<int,int>;
    std::vector<std::vector<Ref>> prior(g.cells.size()),next(g.cells.size());
    for(int r=0;r<int(old.waterRegions.size());r++)for(auto [cell,k]:old.waterRegions[r].pieces)prior[cell].push_back({r,k});
    for(int r=0;r<int(g.waterRegions.size());r++)for(auto [cell,k]:g.waterRegions[r].pieces)next[cell].push_back({r,k});
    std::vector<std::vector<double>> inherited(g.cells.size());for(int cell=0;cell<int(g.cells.size());cell++)if(!next[cell].empty())inherited[cell].resize(g.cells[cell].size());
    std::vector<std::map<int,double>> weights(old.waterRegions.size());
    BooleanEnvironment environment;auto context=booleanContext().handle;
    for(int cell=0;cell<int(g.cells.size());cell++){
        const auto& sources=prior[cell];const auto& targets=next[cell];if(sources.empty()||targets.empty())continue;
        using Box=std::array<double,4>;
        auto bounds=[](const Piece& piece){Box b{INFINITY,INFINITY,-INFINITY,-INFINITY};for(auto p:piece.polygon){b[0]=std::min(b[0],p.x);b[1]=std::min(b[1],p.y);b[2]=std::max(b[2],p.x);b[3]=std::max(b[3],p.y);}return b;};
        auto overlaps=[](const Box& a,const Box& b){constexpr double pad=4e-10;return a[0]<=b[2]+pad&&b[0]<=a[2]+pad&&a[1]<=b[3]+pad&&b[1]<=a[3]+pad;};
        std::vector<Box> targetBounds;std::vector<std::pair<Point,Point>> boxes;
        for(auto [r,k]:targets){auto b=bounds(g.cells[cell][k]);targetBounds.push_back(b);boxes.push_back({{b[0],b[1]},{b[2],b[3]}});}
        std::unique_ptr<BoundaryEdgeIndex> lookup;if(sources.size()*targets.size()>64)lookup=std::make_unique<BoundaryEdgeIndex>(boxes);
        std::vector<Geos> beforeCache(sources.size()),afterCache(targets.size());std::vector<std::size_t> candidates;
        for(std::size_t u=0;u<sources.size();u++){
            auto [a,i]=sources[u];const auto& before=old.cells[cell][i];auto box=bounds(before);
            if(lookup)lookup->candidates({box[0],box[1]},{box[2],box[3]},candidates);
            else{candidates.resize(targets.size());std::iota(candidates.begin(),candidates.end(),0);}
            for(auto v:candidates){
                if(!overlaps(box,targetBounds[v]))continue;auto [b,j]=targets[v];const auto& after=g.cells[cell][j];double overlap=0;
                if(samePiece(before,after))overlap=pieceArea(before);
                else{
                    if(!beforeCache[u])beforeCache[u]=region(before);if(!afterCache[v])afterCache[v]=region(after);
                    auto common=intersection(beforeCache[u].get(),afterCache[v].get());GEOSArea_r(context,common.get(),&overlap);
                }
                if(overlap>0){weights[a][b]+=overlap*(1-g.epsilon)/(g.ratio*g.ratio);inherited[cell][j]+=overlap;}
            }
        }
    }
    // Newly dissolved space belongs to the liquid that supplied that reaction.
    // Include its swept volume even when it crosses a main-cell boundary.
    std::map<int,std::vector<const CreatedWater*>> creators;
    for(const auto& water:g.createdWater){
        if(water.region<0||water.region>=int(weights.size()))throw std::logic_error("Invalid water provenance for interface motion");
        creators[water.cell].push_back(&water);
    }
    // Attribute only space actually reached by each dissolving interface.
    // Normalise overlapping fronts within each final connected liquid piece,
    // never across disconnected pieces merely sharing the same fine cell.
    for(const auto& entry:creators){int cell=entry.first;std::vector<Geos> sources;for(const auto* water:entry.second)sources.push_back(unionPieces(water->shape));
        for(auto [r,k]:next[cell]){
            double created=std::max(0.,pieceArea(g.cells[cell][k])-inherited[cell][k]);if(created==0)continue;
            auto target=region(g.cells[cell][k]);std::map<int,double> supplied;double sum=0;
            for(std::size_t i=0;i<sources.size();i++){auto overlap=intersection(sources[i].get(),target.get());double amount=0;GEOSArea_r(context,overlap.get(),&amount);amount*=entry.second[i]->weight;if(amount>0){supplied[entry.second[i]->region]+=amount;sum+=amount;}}
            if(sum>0)for(auto [source,amount]:supplied)weights[source][r]+=created*amount/sum*(1-g.epsilon)/(g.ratio*g.ratio);
        }
    }
    std::vector<WaterStock> incoming(g.waterRegions.size()),bulk(g.nx*g.ny);
    for(int p=0;p<g.nx*g.ny;p++)if(old.bulkRegion[p]<0)bulk[p]=waterStock(mobile[p]);
    auto project=[&](int r,const WaterStock& stock){
        const auto& source=old.waterRegions[r];double volume=source.bulk?g.epsilon:0;for(auto [id,a]:weights[r])volume+=a;if(volume<=0)return false;
        WaterStock assigned{};int last=-1;for(auto [id,a]:weights[r]){addWater(incoming[id],stock,a/volume);addWater(assigned,stock,a/volume);last=id;}
        WaterStock residual=stock;addWater(residual,assigned,-1);if(source.bulk||last<0)addWater(bulk[source.parent],residual);else addWater(incoming[last],residual);return true;
    };
    for(int r=0;r<int(old.waterRegions.size());r++){
        const auto& source=old.waterRegions[r];
        const auto stock=source.bulk?waterStock(mobile[source.parent]):(reactedPockets?*reactedPockets:old.pockets).at(old.pocketRegion.at(r)).inventory;
        if(!project(r,stock)){
            // A finite isolated water volume remains present through an event.
            // A precipitating pool may already have expelled almost all its
            // water. Send only its numerical remnant through the old hydraulic
            // connections, rather than through solid into its parent cell.
            double remaining=source.fraction;
            if(!source.bulk){double priorMass=waterMass(old.pockets.at(old.pocketRegion.at(r)).inventory);if(priorMass>0)remaining*=waterMass(stock)/priorMass;}
            if(remaining<=1e-9){
                auto recipients=closingWaterRecipients(old,r,weights,g.epsilon);
                if(!recipients.empty()){
                    WaterStock left=stock;for(std::size_t k=0;k<recipients.size();k++){auto [receiver,fraction]=recipients[k];auto portion=stock;if(k+1==recipients.size())portion=left;else for(double& v:portion)v*=fraction;addWater(left,portion,-1);if(!project(receiver,portion))throw std::logic_error("Closing-water receiver has no surviving volume");}
                    continue;
                }
            }
            // Retain the host-ledger rule only for water that was already
            // below geometric resolution before the reaction.
            if((source.fraction>1e-9||remaining>1e-9)&&!source.bulk){std::ostringstream message;message<<std::setprecision(17)<<"A reacting local water pool reached its volume floor: region="<<r<<" old_fraction="<<source.fraction<<" remaining_fraction="<<remaining<<" mass="<<waterMass(stock)<<" connected="<<source.connected;throw StepFailure(message.str());}
            addWater(bulk[source.parent],stock);continue;
        }
    }
    makePockets(g);
    for(int r=0;r<int(g.waterRegions.size());r++){
        const auto& region=g.waterRegions[r];
        if(region.bulk)addWater(bulk[region.parent],incoming[r]);else g.pockets[g.pocketRegion[r]].inventory=incoming[r];
    }
    for(int p=0;p<g.nx*g.ny;p++)setWater(mobile[p],bulk[p]);
    g.waterInitialized=true;
    // Fill newly represented, connected void slivers conservatively from
    // their actual neighbours. This also handles microscopic polygon cleanup.
    for(std::size_t pass=0;pass<g.waterRegions.size();pass++){
        bool changed=false;
        for(const auto& edge:g.waterLinks){
            auto a=regionStock(g,mobile,edge.a),b=regionStock(g,mobile,edge.b);double ma=waterMass(a),mb=waterMass(b);
            if((ma>0)==(mb>0))continue;
            int donor=ma>0?edge.a:edge.b,receiver=ma>0?edge.b:edge.a;auto stock=ma>0?a:b;
            double vd=g.waterRegions[donor].fraction+(g.waterRegions[donor].bulk?g.epsilon:0),vr=g.waterRegions[receiver].fraction;
            if(!(vd>0&&vr>0))continue;auto transfer=stock;double share=vr/(vd+vr);for(double& v:transfer)v*=share;addWater(stock,transfer,-1);
            setRegionStock(g,mobile,donor,stock);setRegionStock(g,mobile,receiver,transfer);changed=true;
        }
        if(!changed)break;
    }
    waterLedger(g);g.createdWater.clear();
}

struct WaterExchangeStats {int nodes=0,links=0,iterations=0;double residual=0;};
// A common M-matrix for all seven conserved components and solvent.
// Only links involving a retained local pool are included: the main-grid
// finite-volume operator owns bulk-to-bulk transport exactly once.
inline WaterExchangeStats exchangeWater(Geometry& g,std::vector<saltchem::Inventory>& mobile,double diffusivity,double dt,double minimumVolume=1e-10){
    if(!(diffusivity>=0&&dt>=0))throw std::invalid_argument("Invalid local water diffusivity/interval");
    WaterExchangeStats stats;if(dt==0||diffusivity==0||g.waterLinks.empty())return stats;
    std::vector<int> indices(g.waterRegions.size(),-1),regions;
    auto index=[&](int r){if(indices[r]<0){indices[r]=regions.size();regions.push_back(r);}return indices[r];};
    struct Link{int a,b;double k;};std::vector<Link> links;
    for(const auto& edge:g.waterLinks){
        auto volume=[&](int r){return g.waterRegions[r].fraction+(g.waterRegions[r].bulk?g.epsilon:0);};
        if(volume(edge.a)<=minimumVolume||volume(edge.b)<=minimumVolume)continue;
        auto a=regionStock(g,mobile,edge.a),b=regionStock(g,mobile,edge.b);
        double ma=waterMass(a),mb=waterMass(b);if(!(ma>0&&mb>0))continue;
        double va=g.waterRegions[edge.a].fraction+(g.waterRegions[edge.a].bulk?g.epsilon:0);
        double vb=g.waterRegions[edge.b].fraction+(g.waterRegions[edge.b].bulk?g.epsilon:0);
        double density=2/(va/ma+vb/mb);
        links.push_back({index(edge.a),index(edge.b),dt*diffusivity*edge.conductance*density/g.coarseVolume()});
    }
    int n=regions.size();if(n==0)return stats;stats.nodes=n;stats.links=links.size();
    std::vector<WaterStock> stocks(n);std::vector<double> mass(n),diagonal(n);
    for(int i=0;i<n;i++){stocks[i]=regionStock(g,mobile,regions[i]);mass[i]=waterMass(stocks[i]);diagonal[i]=mass[i];}
    for(const auto& l:links){diagonal[l.a]+=l.k;diagonal[l.b]+=l.k;}
    auto multiply=[&](const std::vector<double>& x,std::vector<double>& y){for(int i=0;i<n;i++)y[i]=diagonal[i]*x[i];for(const auto& l:links){y[l.a]-=l.k*x[l.b];y[l.b]-=l.k*x[l.a];}};
    // Tiny closing reservoirs can have signed pre-diffusion sources and
    // extreme storage ratios. Factor their small SPD network once in extended
    // precision, shared by all components, without a cancellation-prone
    // residual tolerance based on the provisional concentration.
    std::vector<long double> factor;
    if(n<=64){
        factor.assign(n*n,0);for(int i=0;i<n;i++)factor[i*n+i]=mass[i];
        for(const auto& l:links){factor[l.a*n+l.a]+=l.k;factor[l.b*n+l.b]+=l.k;factor[l.a*n+l.b]-=l.k;factor[l.b*n+l.a]-=l.k;}
        for(int i=0;i<n;i++)for(int j=0;j<=i;j++){long double value=factor[i*n+j];for(int k=0;k<j;k++)value-=factor[i*n+k]*factor[j*n+k];if(i==j){if(!(value>0))throw StepFailure("Local water exchange factor lost positivity");factor[i*n+j]=std::sqrt(value);}else factor[i*n+j]=value/factor[j*n+j];}
    }
    const int largest=std::max_element(mass.begin(),mass.end())-mass.begin();
    for(int component=0;component<8;component++){
        std::vector<double> x(n),rhs(n),r(n),z(n),p(n),ap(n);double scale=0;
        long double initial=0;
        for(int i=0;i<n;i++){rhs[i]=stocks[i][component]*massFactor(component);x[i]=rhs[i]/mass[i];scale=std::max(scale,x[i]);initial+=rhs[i];}
        if(scale==0)continue;scale=std::min(1.,scale);
        if(!factor.empty()){
            std::vector<long double> answer(n);for(int i=0;i<n;i++){long double value=rhs[i];for(int j=0;j<i;j++)value-=factor[i*n+j]*answer[j];answer[i]=value/factor[i*n+i];}
            for(int i=n-1;i>=0;i--){long double value=answer[i];for(int j=i+1;j<n;j++)value-=factor[j*n+i]*answer[j];answer[i]=value/factor[i*n+i];x[i]=double(answer[i]);}
            stats.iterations++;
        }
        multiply(x,ap);double rz=0;
        for(int i=0;i<n;i++){r[i]=rhs[i]-ap[i];z[i]=r[i]/diagonal[i];p[i]=z[i];rz+=r[i]*z[i];}
        double error=0;
        for(int iteration=0;factor.empty()&&iteration<std::max(64,4*n);iteration++){
            error=0;for(int i=0;i<n;i++)error=std::max(error,std::abs(r[i])/diagonal[i]);
            if(error<=2e-13*scale+1e-16)break;
            multiply(p,ap);double pap=0;for(int i=0;i<n;i++)pap+=p[i]*ap[i];
            if(!(pap>0&&rz>0))throw StepFailure("Local water exchange lost its positive matrix");
            double alpha=rz/pap,newRz=0;
            for(int i=0;i<n;i++){x[i]+=alpha*p[i];r[i]-=alpha*ap[i];z[i]=r[i]/diagonal[i];newRz+=r[i]*z[i];}
            double beta=newRz/rz;for(int i=0;i<n;i++)p[i]=z[i]+beta*p[i];rz=newRz;stats.iterations++;
        }
        multiply(x,ap);error=0;long double final=0;
        for(int i=0;i<n;i++){
            error=std::max(error,std::abs(ap[i]-rhs[i])/diagonal[i]);
            if(x[i]<-1e-12*scale)throw StepFailure("Local water exchange produced a negative component");
            stocks[i][component]=std::max(0.,x[i])*mass[i]/massFactor(component);final+=stocks[i][component]*massFactor(component);
        }
        stats.residual=std::max(stats.residual,error);
        if(error>2e-11*scale+1e-15)throw StepFailure("Local water exchange residual exceeded tolerance");
        stocks[largest][component]+=double(initial-final)/massFactor(component);
        if(stocks[largest][component]<0)throw StepFailure("Local water exchange conservation correction exceeded inventory");
    }
    for(int i=0;i<n;i++)setRegionStock(g,mobile,regions[i],stocks[i]);waterLedger(g);return stats;
}
}
