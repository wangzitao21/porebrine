// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "MineralInterface.hpp"
#include <fstream>
#include <map>
#include <numeric>
#include <set>
#include <queue>
#include <tuple>

namespace saltdual {
using saltgeom::Point;using saltgeom::dot;using saltgeom::cross;using saltgeom::length;
using saltgeom::operator+;using saltgeom::operator-;using saltgeom::operator*;
using Poly=std::vector<Point>;
constexpr double tiny=2e-12;
constexpr std::size_t maxFineCellPieces=4096;
inline double area(const Poly& p){double v=0;for(std::size_t k=1;k+1<p.size();k++)v+=cross(p[k]-p[0],p[k+1]-p[0]);return std::abs(v)*.5;}
inline Poly clean(Poly p){
    Poly q;for(auto v:p)if(q.empty()||length(v-q.back())>tiny)q.push_back(v);
    if(q.size()>1&&length(q.front()-q.back())<=tiny)q.pop_back();
    bool changed=true;while(changed&&q.size()>2){changed=false;for(std::size_t i=0;i<q.size();i++){
        const auto a=q[i]-q[(i+q.size()-1)%q.size()],b=q[(i+1)%q.size()]-q[i];
        if(std::abs(cross(a,b))<=tiny*(length(a)+length(b))&&dot(a,b)>=0){q.erase(q.begin()+i);changed=true;break;}
    }}
    return area(q)>1e-16?q:Poly{};
}
inline Poly clip(const Poly& p,Point n,double d){
    Poly q;if(p.empty())return q;
    for(std::size_t i=0;i<p.size();i++){
        auto a=p[i],b=p[(i+1)%p.size()];double da=dot(n,a)-d,db=dot(n,b)-d;
        if(da<=0)q.push_back(a);
        if((da<0&&db>0)||(da>0&&db<0))q.push_back(a+(b-a)*(da/(da-db)));
    }return clean(std::move(q));
}
inline Poly square(){return {{0,0},{1,0},{1,1},{0,1}};}
inline Poly hull(Poly p){
    std::sort(p.begin(),p.end(),[](Point a,Point b){return a.x<b.x||(a.x==b.x&&a.y<b.y);});
    p.erase(std::unique(p.begin(),p.end(),[](Point a,Point b){return length(a-b)<tiny;}),p.end());
    if(p.size()<3)return {};
    Poly h(2*p.size());std::size_t k=0;
    for(auto v:p){while(k>=2&&cross(h[k-1]-h[k-2],v-h[k-1])<=0)--k;h[k++]=v;}
    const auto lower=k+1;for(int i=int(p.size())-2;i>=0;i--){while(k>=lower&&cross(h[k-1]-h[k-2],p[i]-h[k-1])<=0)--k;h[k++]=p[i];}
    h.resize(k-1);return clean(std::move(h));
}
inline double common(Point a,Point b,Point c,Point d,Point& start,Point& end){
    // A positive common segment must have overlapping coordinate bounds.
    // Keep a margin larger than the collinearity tolerance below so this
    // broad-phase rejection leaves all accepted contacts bit-for-bit intact.
    constexpr double margin=4e-10;
    if(std::max(a.x,b.x)+margin<std::min(c.x,d.x)||std::max(c.x,d.x)+margin<std::min(a.x,b.x)||
       std::max(a.y,b.y)+margin<std::min(c.y,d.y)||std::max(c.y,d.y)+margin<std::min(a.y,b.y))return 0;
    double l=length(b-a);if(l<tiny)return 0;auto t=(b-a)*(1/l);
    if(std::abs(cross(t,c-a))>2e-10||std::abs(cross(t,d-a))>2e-10)return 0;
    double lo=std::max(0.,std::min(dot(c-a,t),dot(d-a,t))),hi=std::min(l,std::max(dot(c-a,t),dot(d-a,t)));
    if(hi-lo<tiny)return 0;start=a+t*lo;end=a+t*hi;return hi-lo;
}
struct Piece {int phase=-1;Poly polygon;std::vector<Poly> holes;};
inline double pieceArea(const Piece& p){long double a=area(p.polygon);for(const auto& h:p.holes)a-=area(h);return std::max(0.,double(a));}
inline Point pieceCentre(const Piece& p){
    double a=0;Point centre{};
    auto add=[&](const Poly& ring){for(std::size_t i=0;i<ring.size();i++){auto u=ring[i],v=ring[(i+1)%ring.size()];double w=cross(u,v);a+=w;centre=centre+(u+v)*w;}};
    add(p.polygon);for(const auto& h:p.holes)add(h);
    return std::abs(a)>1e-16?centre*(1/(3*a)):p.polygon.front();
}
using Cell=std::vector<Piece>;
inline double paint(Cell& cell,int from,int to,const std::vector<std::pair<Point,double>>& planes);
inline void merge(Cell& cell,std::array<double,3>* removedArea=nullptr);

struct Surface {int solid,water,mineral;Point a,b,normal;int liquidPiece=-1;bool connected=true;int liquidRegion=-1,solidPiece=-1;};
struct StepFailure:std::runtime_error {
    double scale=.5;
    explicit StepFailure(const std::string& message,double fraction=.5):std::runtime_error(message),scale(fraction){}
};
struct WaterRegion {
    int parent=-1;bool bulk=false,connected=false;
    double fraction=0;Point centre{};
    std::vector<std::pair<int,int>> pieces;
};
struct WaterLink {int a=-1,b=-1;double conductance=0,hydraulicConductance=0;};
struct WaterPocket {
    int parent=-1;
    double fraction=0;
    std::array<double,8> inventory{};
};
struct SurfaceBudget {
    int parent=-1,mineral=-1;
    double change=0;
    std::vector<int> surfaces;
    int liquidRegion=-1;
};
struct CreatedWater {int region=-1,cell=-1;Cell shape;double weight=1;};

// Query only edges whose bounding rectangles can contain a shared segment.
// Candidate indices retain the original ring order before the exact predicate
// is called, preserving contact summation and liquid-region numbering.
class BoundaryEdgeIndex {
    using Box=std::array<double,4>;
    struct Node {Box box;std::size_t begin,end;int left=-1,right=-1;};
    std::vector<Box> boxes;
    std::vector<std::size_t> order;
    std::vector<Node> nodes;
    static bool overlap(const Box& a,const Box& b){
        constexpr double margin=4e-10;
        return a[0]<=b[2]+margin&&b[0]<=a[2]+margin&&a[1]<=b[3]+margin&&b[1]<=a[3]+margin;
    }
    int build(std::size_t begin,std::size_t end){
        Box box{INFINITY,INFINITY,-INFINITY,-INFINITY};
        for(std::size_t k=begin;k<end;k++){const auto& b=boxes[order[k]];for(int d=0;d<2;d++){box[d]=std::min(box[d],b[d]);box[d+2]=std::max(box[d+2],b[d+2]);}}
        const int id=int(nodes.size());nodes.push_back({box,begin,end});
        if(end-begin>8){
            const int axis=box[2]-box[0]>=box[3]-box[1]?0:1;const auto middle=begin+(end-begin)/2;
            std::nth_element(order.begin()+begin,order.begin()+middle,order.begin()+end,[&](std::size_t a,std::size_t b){double x=boxes[a][axis]+boxes[a][axis+2],y=boxes[b][axis]+boxes[b][axis+2];return x<y||(x==y&&a<b);});
            const int left=build(begin,middle),right=build(middle,end);nodes[id].left=left;nodes[id].right=right;
        }
        return id;
    }
    void query(int id,const Box& box,std::vector<std::size_t>& found)const{
        const auto& node=nodes[id];if(!overlap(node.box,box))return;
        if(node.left<0){for(auto k=node.begin;k<node.end;k++)if(overlap(boxes[order[k]],box))found.push_back(order[k]);}
        else{query(node.left,box,found);query(node.right,box,found);}
    }
    static long double lowerDistance(Point p,const Box& box){
        // Outward padding keeps the lower bound conservative at roundoff.
        constexpr long double margin=4e-10L;
        long double x=std::max({0.L,static_cast<long double>(box[0])-p.x-margin,static_cast<long double>(p.x)-box[2]-margin});
        long double y=std::max({0.L,static_cast<long double>(box[1])-p.y-margin,static_cast<long double>(p.y)-box[3]-margin});return x*x+y*y;
    }
    template<class Distance>
    void nearest(int id,Point p,const Distance& distance,long double& best)const{
        const auto& node=nodes[id];if(lowerDistance(p,node.box)>best)return;
        if(node.left<0){for(auto k=node.begin;k<node.end;k++)if(lowerDistance(p,boxes[order[k]])<=best)best=std::min(best,distance(order[k]));}
        else{int a=node.left,b=node.right;if(lowerDistance(p,nodes[b].box)<lowerDistance(p,nodes[a].box))std::swap(a,b);nearest(a,p,distance,best);nearest(b,p,distance,best);}
    }
    void initialize(){order.resize(boxes.size());std::iota(order.begin(),order.end(),0);nodes.reserve(2*boxes.size());if(!boxes.empty())build(0,boxes.size());}
public:
    explicit BoundaryEdgeIndex(const Poly& ring){
        boxes.reserve(ring.size());for(std::size_t k=0;k<ring.size();k++){auto a=ring[k],b=ring[(k+1)%ring.size()];boxes.push_back({std::min(a.x,b.x),std::min(a.y,b.y),std::max(a.x,b.x),std::max(a.y,b.y)});}
        initialize();
    }
    explicit BoundaryEdgeIndex(const std::vector<std::pair<Point,Point>>& segments){
        boxes.reserve(segments.size());for(auto [a,b]:segments)boxes.push_back({std::min(a.x,b.x),std::min(a.y,b.y),std::max(a.x,b.x),std::max(a.y,b.y)});initialize();
    }
    template<class Distance>
    long double nearestSquared(Point p,const Distance& distance)const{
        long double best=std::numeric_limits<long double>::infinity();if(!nodes.empty())nearest(0,p,distance,best);return best;
    }
    void candidates(Point a,Point b,std::vector<std::size_t>& found)const{
        found.clear();if(nodes.empty())return;query(0,{std::min(a.x,b.x),std::min(a.y,b.y),std::max(a.x,b.x),std::max(a.y,b.y)},found);std::sort(found.begin(),found.end());
    }
};

}
#include "GeometryBoolean.hpp"
namespace saltdual {
class Geometry {
public:
    int nx=0,ny=0,ratio=0;double spacing=0,depth=0,epsilon=0;
    std::vector<Cell> cells;
    std::vector<Surface> wet,buried;
    bool storeBuriedInterfaces=true;
    std::vector<std::array<double,3>> fractions,dissolution,growth,covered;
    std::vector<int> liquidComponents;
    int globalLiquidComponents=0;
    std::vector<double> globalLiquidAreas;
    double isolatedLiquidArea=0;
    std::uint64_t booleanFallbacks=0;
    double chemicalVolumeError=0;
    std::array<double,3> chemicalVolumeDebt{};
    std::vector<std::array<double,4>> openings;
    std::vector<std::array<double,4>> diffusiveOpenings;
    std::vector<std::array<double,9>> nucleationHazard;
    std::vector<std::array<std::uint64_t,9>> nucleationEvents;
    std::vector<double> closedFraction,storedClosedFraction;
    std::vector<std::array<double,8>> closedInventory;
    std::map<int,Cell> closedWater;
    // Geometry-derived aqueous control volumes. Broad fine-grid connections
    // are combined within a main cell; narrow throats retain a finite link.
    double waterNeckFraction=.1;
    double waterCouplingFraction=.5;
    std::vector<WaterRegion> waterRegions;
    std::vector<WaterLink> waterLinks;
    std::vector<int> bulkRegion;
    std::vector<WaterPocket> pockets;
    std::vector<int> pocketRegion;
    bool waterInitialized=false;
    std::vector<CreatedWater> createdWater;
    std::vector<Point> sweptVolumeVectors;
    bool lastMovementChanged=false;
    int fx()const{return nx*ratio;}int fy()const{return ny*ratio;}
    int parent(int g)const{return (g/fx()/ratio)*nx+(g%fx()/ratio);}
    Point origin(int g)const{return {double(g%fx()),double(g/fx())};}
    double h()const{return spacing/ratio;}
    double unitVolume()const{return h()*h()*depth*(1-epsilon);}
    double coarseVolume()const{return spacing*spacing*depth;}
    Geometry()=default;
    Geometry(int x,int y,int r,double H,double z,double ep):nx(x),ny(y),ratio(r),spacing(H),depth(z),epsilon(ep),cells(x*y*r*r,Cell{{-1,square()}}),nucleationHazard(x*y),nucleationEvents(x*y),closedFraction(x*y),storedClosedFraction(x*y),closedInventory(x*y){validateHeader();}
    void validateHeader()const{
        if(nx<1||ny<1||ratio<1||ratio>32||std::uint64_t(nx)*ny*ratio*ratio>10000000||!(spacing>0&&depth>0&&epsilon>=0&&epsilon<.01))throw std::invalid_argument("Invalid dual geometry dimensions");
    }
    void circle(int phase,double x,double y,double radius,int facets=256){
        if(phase<0||phase>2||!(radius>0)||facets<16)throw std::invalid_argument("Invalid mineral circle");
        x/=h();y/=h();radius/=h();
        for(int iy=std::max(0,int(std::floor(y-radius)));iy<std::min(fy(),int(std::ceil(y+radius)));iy++)
        for(int ix=std::max(0,int(std::floor(x-radius)));ix<std::min(fx(),int(std::ceil(x+radius)));ix++){
            std::vector<std::pair<Point,double>> planes;
            for(int k=0;k<facets;k++){double angle=2*saltgeom::pi*(k+.5)/facets;Point n{std::cos(angle),std::sin(angle)};planes.push_back({n,radius*std::cos(saltgeom::pi/facets)+dot(n,{x-ix,y-iy})});}
            paint(cells[iy*fx()+ix],-1,phase,planes);merge(cells[iy*fx()+ix]);
        }
    }
    void rebuild(bool topology=true){
        wet.clear();buried.clear();fractions.assign(nx*ny,{});dissolution.assign(nx*ny,{});growth.assign(nx*ny,{});covered.assign(nx*ny,{});openings.assign(nx*ny,{});
        std::vector<std::vector<int>> ids(cells.size());std::vector<int> roots;
        struct Contact {int a,b;double length,distance;};std::vector<Contact> contacts;
        std::vector<int> nodeCells,nodePieces;std::vector<Point> centres;
        struct Aperture {int parent,direction,liquid;double amount;};std::vector<Aperture> apertures;
        for(int g=0;g<int(cells.size());g++){
            double sum=0;for(const auto& p:cells[g]){double a=pieceArea(p);sum+=a;if(p.phase>=0)fractions[parent(g)][p.phase]+=a*(1-epsilon)/(ratio*ratio);}
            if(std::abs(sum-1)>2e-9)throw std::runtime_error("Dual geometry lost its material partition");
            if(topology)for(std::size_t k=0;k<cells[g].size();k++){const auto& p=cells[g][k];int id=-1;if(p.phase<0){id=roots.size();roots.push_back(id);nodeCells.push_back(g);nodePieces.push_back(k);centres.push_back(pieceCentre(p)+origin(g));}ids[g].push_back(id);}
        }
        auto globalRoots=roots,mixingRoots=roots;
        auto root=[&](int a){while(roots[a]!=a){roots[a]=roots[roots[a]];a=roots[a];}return a;};
        auto globalRoot=[&](int a){while(globalRoots[a]!=a){globalRoots[a]=globalRoots[globalRoots[a]];a=globalRoots[a];}return a;};
        auto mixingRoot=[&](int a){while(mixingRoots[a]!=a){mixingRoots[a]=mixingRoots[mixingRoots[a]];a=mixingRoots[a];}return a;};
        auto connect=[&](int ga,int gb,int a,int b,double l){
            if(parent(ga)==parent(gb)){
                roots[root(a)]=root(b);
                if(l>=waterNeckFraction)mixingRoots[mixingRoot(a)]=mixingRoot(b);
            }
            globalRoots[globalRoot(a)]=globalRoot(b);
            contacts.push_back({a,b,l,std::max(.05,length(centres[a]-centres[b]))});
        };
        std::map<const Poly*,BoundaryEdgeIndex> edgeIndices;
        auto pair=[&](int ga,int gb,bool same){
            Point shift=origin(gb)-origin(ga);
            // Most of the uniform fine grid is single-phase and stationary.
            // Its exact shared face needs no polygon-edge intersection.
            if(!same&&cells[ga].size()==1&&cells[gb].size()==1&&cells[ga][0].phase==cells[gb][0].phase){
                if(cells[ga][0].phase>=0)return;
                if(topology)connect(ga,gb,ids[ga][0],ids[gb][0],1.);
                if(parent(ga)!=parent(gb)){int dir=shift.x!=0?1:3;openings[parent(ga)][dir]+=1./ratio;openings[parent(gb)][dir^1]+=1./ratio;if(topology){apertures.push_back({parent(ga),dir,ids[ga][0],1./ratio});apertures.push_back({parent(gb),dir^1,ids[gb][0],1./ratio});}}
                return;
            }
            for(std::size_t ia=0;ia<cells[ga].size();ia++)for(std::size_t ib=0;ib<cells[gb].size();ib++){
                if(same&&ia>=ib)continue;const auto& a=cells[ga][ia];const auto& b=cells[gb][ib];
                if(a.phase==b.phase&&a.phase>=0)continue;
                for(std::size_t ar=0;ar<=a.holes.size();ar++)for(std::size_t br=0;br<=b.holes.size();br++){
                const auto& pa=ar?a.holes[ar-1]:a.polygon;const auto& pb=br?b.holes[br-1]:b.polygon;
                auto contact=[&](std::size_t u,std::size_t v){
                    Point start,end;double l=common(pa[u],pa[(u+1)%pa.size()],pb[v]+shift,pb[(v+1)%pb.size()]+shift,start,end);if(l==0)return;
                    auto tangent=(end-start)*(1/l);Point normal{tangent.y,-tangent.x};
                    if(a.phase<0&&b.phase<0){
                        if(topology)connect(ga,gb,ids[ga][ia],ids[gb][ib],l);
                        if(parent(ga)!=parent(gb)){int dir=shift.x!=0?1:3;openings[parent(ga)][dir]+=l/ratio;openings[parent(gb)][dir^1]+=l/ratio;if(topology){apertures.push_back({parent(ga),dir,ids[ga][ia],l/ratio});apertures.push_back({parent(gb),dir^1,ids[gb][ib],l/ratio});}}
                    }else if(a.phase>=0&&b.phase>=0){
                        if(storeBuriedInterfaces){
                            buried.push_back({ga,gb,a.phase,start+origin(ga),end+origin(ga),normal});
                            buried.push_back({gb,ga,b.phase,end+origin(ga),start+origin(ga),normal*(-1)});
                        }
                        covered[parent(ga)][a.phase]+=l*h()*depth;covered[parent(gb)][b.phase]+=l*h()*depth;
                    }else{
                        Surface s=a.phase>=0?Surface{ga,gb,a.phase,start+origin(ga),end+origin(ga),normal}:Surface{gb,ga,b.phase,end+origin(ga),start+origin(ga),normal*(-1)};
                        s.liquidPiece=a.phase>=0?int(ib):int(ia);
                        s.solidPiece=a.phase>=0?int(ia):int(ib);
                        wet.push_back(s);dissolution[parent(s.solid)][s.mineral]+=l*h()*depth;growth[parent(s.water)][s.mineral]+=l*h()*depth;
                    }
                };
                if(pb.size()>=32){
                    const auto& lookup=edgeIndices.try_emplace(&pb,pb).first->second;std::vector<std::size_t> candidates;
                    for(std::size_t u=0;u<pa.size();u++){lookup.candidates(pa[u]-shift,pa[(u+1)%pa.size()]-shift,candidates);for(auto v:candidates)contact(u,v);}
                }else for(std::size_t u=0;u<pa.size();u++)for(std::size_t v=0;v<pb.size();v++)contact(u,v);
                }
            }
        };
        for(int g=0;g<int(cells.size());g++){
            pair(g,g,true);if(g%fx()+1<fx())pair(g,g+1,false);if(g/fx()+1<fy())pair(g,g+fx(),false);
            // Boundary apertures use the same physical liquid cross section
            // as internal faces, including grains touching an inlet/wall.
            if(g%fx()==0||g%fx()+1==fx()||g/fx()==0||g/fx()+1==fy())
            for(std::size_t ip=0;ip<cells[g].size();ip++){const auto& piece=cells[g][ip];if(piece.phase<0)for(std::size_t k=0;k<piece.polygon.size();k++){
                auto a=piece.polygon[k],b=piece.polygon[(k+1)%piece.polygon.size()];
                auto add=[&](int d){double l=length(b-a)/ratio;openings[parent(g)][d]+=l;if(topology)apertures.push_back({parent(g),d,ids[g][ip],l});};
                if(g%fx()==0&&std::abs(a.x)+std::abs(b.x)<tiny)add(0);
                if(g%fx()+1==fx()&&std::abs(a.x-1)+std::abs(b.x-1)<tiny)add(1);
                if(g/fx()==0&&std::abs(a.y)+std::abs(b.y)<tiny)add(2);
                if(g/fx()+1==fy()&&std::abs(a.y-1)+std::abs(b.y-1)<tiny)add(3);
            }}
        }
        if(topology){
            std::vector<std::map<int,double>> regions(nx*ny);
            for(int g=0;g<int(cells.size());g++)for(std::size_t k=0;k<ids[g].size();k++)if(ids[g][k]>=0)regions[parent(g)][root(ids[g][k])]+=pieceArea(cells[g][k]);
            liquidComponents.assign(nx*ny,0);
            // Connectivity has the same absolute area resolution as the
            // polygon predicates. Roundoff fragments retain their volume in
            // the ledger but do not create new physical aqueous reservoirs.
            for(int p=0;p<nx*ny;p++)for(auto region:regions[p])if(region.second>1e-8)liquidComponents[p]++;
            std::map<int,double> globalAreas;
            for(int g=0;g<int(cells.size());g++)for(std::size_t k=0;k<ids[g].size();k++)if(ids[g][k]>=0)globalAreas[globalRoot(ids[g][k])]+=pieceArea(cells[g][k]);
            std::set<int> exterior;
            for(const auto& a:apertures)if((a.direction==0&&a.parent%nx==0)||(a.direction==1&&a.parent%nx==nx-1)||(a.direction==2&&a.parent<nx)||(a.direction==3&&a.parent>=(ny-1)*nx))exterior.insert(globalRoot(a.liquid));
            waterRegions.clear();waterLinks.clear();bulkRegion.assign(nx*ny,-1);
            std::map<int,int> regionIds;std::vector<int> nodeRegion(roots.size());
            for(int id=0;id<int(roots.size());id++){
                int key=mixingRoot(id);auto found=regionIds.emplace(key,waterRegions.size());
                if(found.second){WaterRegion r;r.parent=parent(nodeCells[id]);r.connected=exterior.count(globalRoot(id));waterRegions.push_back(r);}
                int r=found.first->second;nodeRegion[id]=r;auto& region=waterRegions[r];
                double a=pieceArea(cells[nodeCells[id]][nodePieces[id]])*(1-epsilon)/(ratio*ratio);
                region.fraction+=a;region.centre=region.centre+centres[id]*a;region.pieces.push_back({nodeCells[id],nodePieces[id]});
            }
            for(int r=0;r<int(waterRegions.size());r++){
                auto& region=waterRegions[r];if(region.fraction>0)region.centre=region.centre*(1/region.fraction);
                int& chosen=bulkRegion[region.parent];
                if(region.connected&&(chosen<0||region.fraction>waterRegions[chosen].fraction))chosen=r;
            }
            for(int r:bulkRegion)if(r>=0)waterRegions[r].bulk=true;
            for(auto& s:wet){s.liquidRegion=nodeRegion[ids[s.water][s.liquidPiece]];s.connected=waterRegions[s.liquidRegion].connected;}
            closedFraction.assign(nx*ny,0);closedWater.clear();
            for(const auto& r:waterRegions)if(!r.bulk){closedFraction[r.parent]+=r.fraction;for(auto [g,k]:r.pieces)closedWater[g].push_back(cells[g][k]);}
            openings.assign(nx*ny,{});
            for(const auto& a:apertures)if(nodeRegion[a.liquid]==bulkRegion[a.parent]&&((a.direction==0&&a.parent%nx==0)||(a.direction==1&&a.parent%nx==nx-1)||(a.direction==2&&a.parent<nx)||(a.direction==3&&a.parent>=(ny-1)*nx)))openings[a.parent][a.direction]+=a.amount;
            diffusiveOpenings=openings;std::map<std::pair<int,int>,double> links,hydraulicLinks;
            auto addLink=[&](int a,int b,double aperture,double distance){
                auto pair=std::minmax(a,b);double diffusion=aperture*depth/distance;
                double small=std::min(aperture*h(),depth),large=std::max(aperture*h(),depth);
                links[pair]+=diffusion;hydraulicLinks[pair]+=diffusion*small*small*(1-.63*small/large)/12;
            };
            for(const auto& c:contacts){
                int a=nodeRegion[c.a],b=nodeRegion[c.b];if(a==b)continue;
                const auto& ra=waterRegions[a];const auto& rb=waterRegions[b];
                if(ra.bulk&&rb.bulk){int d=rb.parent-ra.parent==1?1:rb.parent-ra.parent==-1?0:rb.parent>ra.parent?3:2;openings[ra.parent][d]+=c.length/ratio;openings[rb.parent][d^1]+=c.length/ratio;
                    if(std::min(ra.fraction,rb.fraction)>=waterCouplingFraction){diffusiveOpenings[ra.parent][d]+=c.length/ratio;diffusiveOpenings[rb.parent][d^1]+=c.length/ratio;}
                    else addLink(a,b,c.length,std::max(.05,length(ra.centre-rb.centre)));
                }else addLink(a,b,c.length,c.distance);
            }
            for(auto [pair,value]:links)waterLinks.push_back({pair.first,pair.second,value,hydraulicLinks[pair]});
            // Pockets below the already configured numerical pore-volume
            // floor retain their inventory in the main-cell ledger, but
            // cannot supply a reactive connection to the external liquid.
            const double poreResolution=std::max(1e-8,epsilon*ratio*ratio);
            globalLiquidComponents=0;globalLiquidAreas.clear();isolatedLiquidArea=0;
            for(auto region:globalAreas){if(!exterior.count(region.first))isolatedLiquidArea+=region.second;if(region.second>poreResolution)globalLiquidComponents++;if(region.second>1e-8)globalLiquidAreas.push_back(region.second);}
        }
    }
    // Exact volume correction of the normal swept polygons, with a fixed
    // displacement ceiling. Chemistry is retried from its initial state if
    // any budget cannot fit; no unpaired clipping of reaction extents occurs.
    void advance(const std::vector<std::array<double,3>>& target,bool unblocked=false,double maxMove=.2,const std::vector<SurfaceBudget>* budgets=nullptr);
    void advanceImpl(const std::vector<std::array<double,3>>& target,bool unblocked,double maxMove,const std::vector<SurfaceBudget>* budgets,bool fullGrowthBand);
    void write(const std::string& path)const{
        std::ofstream f(path,std::ios::binary);const std::uint32_t head[]{0x44554752,6,std::uint32_t(nx),std::uint32_t(ny),std::uint32_t(ratio)};
        const double values[]{spacing,depth,epsilon};f.write(reinterpret_cast<const char*>(head),sizeof(head));f.write(reinterpret_cast<const char*>(values),sizeof(values));
        auto ring=[&](const Poly& p){std::uint32_t n=p.size();f.write(reinterpret_cast<const char*>(&n),4);for(auto v:p){f.write(reinterpret_cast<const char*>(&v.x),8);f.write(reinterpret_cast<const char*>(&v.y),8);}};
        for(const auto& c:cells){std::uint32_t n=c.size();f.write(reinterpret_cast<const char*>(&n),4);for(const auto& p:c){std::int32_t phase=p.phase;std::uint32_t holes=p.holes.size();f.write(reinterpret_cast<const char*>(&phase),4);ring(p.polygon);f.write(reinterpret_cast<const char*>(&holes),4);for(const auto& h:p.holes)ring(h);}}
        for(int p=0;p<nx*ny;p++){f.write(reinterpret_cast<const char*>(nucleationHazard[p].data()),9*sizeof(double));f.write(reinterpret_cast<const char*>(nucleationEvents[p].data()),9*sizeof(std::uint64_t));}
        for(int p=0;p<nx*ny;p++){f.write(reinterpret_cast<const char*>(&storedClosedFraction[p]),sizeof(double));f.write(reinterpret_cast<const char*>(closedInventory[p].data()),8*sizeof(double));}
        f.write(reinterpret_cast<const char*>(&waterNeckFraction),sizeof(double));
        f.write(reinterpret_cast<const char*>(&waterCouplingFraction),sizeof(double));
        const std::uint32_t initialized=waterInitialized,count=pockets.size();f.write(reinterpret_cast<const char*>(&initialized),4);f.write(reinterpret_cast<const char*>(&count),4);
        for(const auto& p:pockets){std::int32_t parent=p.parent;f.write(reinterpret_cast<const char*>(&parent),4);f.write(reinterpret_cast<const char*>(&p.fraction),8);f.write(reinterpret_cast<const char*>(p.inventory.data()),8*sizeof(double));}
        f.write(reinterpret_cast<const char*>(chemicalVolumeDebt.data()),sizeof(chemicalVolumeDebt));
        if(!f)throw std::runtime_error("Dual geometry write failed");
    }
    static Geometry read(const std::string& path,bool retainBuriedInterfaces=true){
        std::ifstream f(path,std::ios::binary);f.seekg(0,std::ios::end);const auto fileBytes=f.tellg();if(fileBytes<44)throw std::runtime_error("Invalid dual geometry checkpoint length");f.seekg(0);
        // A valid saved polygon can exceed an arbitrary vertex-count cutoff.
        // Bound allocation by the actual input size and reject truncation at
        // the first failed coordinate read instead.
        const auto vertexLimit=static_cast<std::uint64_t>(fileBytes)/16;
        std::uint32_t head[5]{};double values[3]{};f.read(reinterpret_cast<char*>(head),sizeof(head));f.read(reinterpret_cast<char*>(values),sizeof(values));
        if(!f||head[0]!=0x44554752||head[1]<1||head[1]>6)throw std::runtime_error("Invalid dual geometry checkpoint");
        Geometry g(head[2],head[3],head[4],values[0],values[1],values[2]);
        g.storeBuriedInterfaces=retainBuriedInterfaces;
        auto ring=[&](){std::uint32_t nv=0;f.read(reinterpret_cast<char*>(&nv),4);if(!f||nv<3||nv>vertexLimit)throw std::runtime_error("Invalid geometry ring");Poly p(nv);for(auto& v:p){f.read(reinterpret_cast<char*>(&v.x),8);f.read(reinterpret_cast<char*>(&v.y),8);if(!f||!std::isfinite(v.x)||!std::isfinite(v.y)||v.x< -tiny||v.x>1+tiny||v.y< -tiny||v.y>1+tiny)throw std::runtime_error("Invalid fine geometry coordinate");}return p;};
        for(auto& c:g.cells){std::uint32_t n=0;f.read(reinterpret_cast<char*>(&n),4);if(!n||n>maxFineCellPieces)throw std::runtime_error("Invalid fine-cell region count");c.clear();
            for(std::uint32_t k=0;k<n;k++){std::int32_t phase=-2;f.read(reinterpret_cast<char*>(&phase),4);if(phase< -1||phase>2)throw std::runtime_error("Invalid mineral phase");Piece p{phase,ring(),{}};
                if(head[1]>=2){std::uint32_t holes=0;f.read(reinterpret_cast<char*>(&holes),4);if(holes>4096)throw std::runtime_error("Invalid geometry hole count");for(std::uint32_t j=0;j<holes;j++)p.holes.push_back(ring());}
                c.push_back(std::move(p));
            }
        }
        if(head[1]>=3)for(int p=0;p<g.nx*g.ny;p++){f.read(reinterpret_cast<char*>(g.nucleationHazard[p].data()),9*sizeof(double));f.read(reinterpret_cast<char*>(g.nucleationEvents[p].data()),9*sizeof(std::uint64_t));for(double h:g.nucleationHazard[p])if(!std::isfinite(h)||h<0)throw std::runtime_error("Invalid nucleation exposure state");}
        if(head[1]>=4)for(int p=0;p<g.nx*g.ny;p++){f.read(reinterpret_cast<char*>(&g.storedClosedFraction[p]),sizeof(double));f.read(reinterpret_cast<char*>(g.closedInventory[p].data()),8*sizeof(double));if(!std::isfinite(g.storedClosedFraction[p])||g.storedClosedFraction[p]<0||g.storedClosedFraction[p]>1)throw std::runtime_error("Invalid closed pore fraction");for(double a:g.closedInventory[p])if(!std::isfinite(a)||a<0)throw std::runtime_error("Invalid closed-water inventory");}
        if(head[1]>=5){
            std::uint32_t initialized=0,count=0;f.read(reinterpret_cast<char*>(&g.waterNeckFraction),8);f.read(reinterpret_cast<char*>(&g.waterCouplingFraction),8);f.read(reinterpret_cast<char*>(&initialized),4);f.read(reinterpret_cast<char*>(&count),4);
            if(initialized>1||count>g.cells.size()*maxFineCellPieces||!(g.waterNeckFraction>=0&&g.waterNeckFraction<=1)||!(g.waterCouplingFraction>0&&g.waterCouplingFraction<=1))throw std::runtime_error("Invalid local water checkpoint header");
            g.waterInitialized=initialized;g.pockets.resize(count);
            for(auto& p:g.pockets){std::int32_t parent;f.read(reinterpret_cast<char*>(&parent),4);p.parent=parent;f.read(reinterpret_cast<char*>(&p.fraction),8);f.read(reinterpret_cast<char*>(p.inventory.data()),8*sizeof(double));
                if(parent<0||parent>=g.nx*g.ny||!std::isfinite(p.fraction)||p.fraction<0||p.fraction>1)throw std::runtime_error("Invalid local water volume");for(double a:p.inventory)if(!std::isfinite(a)||a<0)throw std::runtime_error("Invalid local water inventory");}
        }
        if(head[1]>=6){f.read(reinterpret_cast<char*>(g.chemicalVolumeDebt.data()),sizeof(g.chemicalVolumeDebt));for(double v:g.chemicalVolumeDebt)if(!std::isfinite(v))throw std::runtime_error("Invalid deferred mineral volume");}
        if(!f||f.peek()!=std::char_traits<char>::eof())throw std::runtime_error("Truncated/trailing dual geometry state");g.rebuild();
        g.pocketRegion.assign(g.waterRegions.size(),-1);
        if(g.waterInitialized){int k=0;for(int r=0;r<int(g.waterRegions.size());r++)if(!g.waterRegions[r].bulk){const auto& region=g.waterRegions[r];
                if(k>=int(g.pockets.size())||g.pockets[k].parent!=region.parent||std::abs(g.pockets[k].fraction-region.fraction)>1e-9)throw std::runtime_error("Local water checkpoint geometry mismatch");
                g.pocketRegion[r]=k;k++;}
            if(k!=int(g.pockets.size()))throw std::runtime_error("Local water checkpoint count mismatch");}
        return g;
    }
};
}


#include "GeometryAdvance.hpp"
