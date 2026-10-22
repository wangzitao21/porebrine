// SPDX-License-Identifier: GPL-3.0-or-later
#include "DualMineralGeometry.hpp"
#include <iostream>
using namespace saltdual;
static void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
static void rectangle(Geometry& g,int from,int to,double x0,double x1,double y0,double y1){
    for(int i=0;i<int(g.cells.size());i++){auto o=g.origin(i);paint(g.cells[i],from,to,{{{1,0},x1-o.x},{{-1,0},o.x-x0},{{0,1},y1-o.y},{{0,-1},o.y-y0}});merge(g.cells[i]);}g.rebuild();
}
int main(){try{
    {std::vector<std::pair<Point,Point>> segments;for(int k=0;k<1000;k++){double a=k*.023,r=.4+.1*std::sin(7*a);Point p{r*std::cos(a),r*std::sin(a)},q{p.x+.01*std::cos(17*a),p.y+.01*std::sin(17*a)};segments.push_back({p,k%19?q:p});}BoundaryEdgeIndex index(segments);
     for(int k=0;k<1000;k++){Point p{k%2?.8*std::sin(k*.019):segments[k].first.x+1e-12,.8*std::cos(k*.017)};auto nearest=index.nearestSquared(p,[&](std::size_t i){return segmentDistanceSquared(p,segments[i].first,segments[i].second);});require(nearest==frontDistanceSquared(p,segments),"Indexed arrival distance changed the competition partition");}}
    for(int mineral=0;mineral<3;mineral++){
        Geometry g(2,1,4,4,1,0);g.circle(mineral,4,2,.1,32);g.circle(mineral,2,2,.5,64);g.rebuild();auto target=g.fractions;
        std::vector<SurfaceBudget> budgets;for(int p=0;p<2;p++){budgets.push_back({p,mineral,-g.fractions[p][mineral],{},g.bulkRegion[p]});target[p][mineral]=0;}
        for(int k=0;k<int(g.wet.size());k++)budgets[g.parent(g.wet[k].water)].surfaces.push_back(k);
        g.advance(target,false,1.,&budgets);
        require(g.fractions[0][mineral]+g.fractions[1][mineral]<1e-10,"Shared-crystal allocation stranded a water region despite sufficient exposed stock");
    }
    std::vector<std::array<double,2>> reference;
    for(int a=0;a<3;a++)for(int b=0;b<3;b++)if(a!=b){
        Geometry g(1,1,4,4,1,0);rectangle(g,-1,a,0,1.9,0,4);rectangle(g,-1,b,2.1,4,0,4);
        const auto original=g;auto target=g.fractions;target[0][a]+=.035;target[0][b]+=.035;double event=0;
        try{g.advance(target);}catch(const StepFailure& e){event=e.scale;}
        require(event>.7&&event<.73,"Competing growth failed to locate first contact");
        g=original;target=g.fractions;target[0][a]+=.035*event;target[0][b]+=.035*event;g.advance(target);
        require(std::abs(g.fractions[0][a]-.5)<2e-9&&std::abs(g.fractions[0][b]-.5)<2e-9,"Mineral ordering changed the competition boundary");
        require(g.covered[0][a]>3.999&&g.covered[0][b]>3.999,"Touching growth fronts did not form a solid-solid interface");
        // The covering mineral can be any of the three salts.
        Geometry film(1,1,4,4,1,0);rectangle(film,-1,a,0,2,0,4);rectangle(film,-1,b,2,2.1,0,4);
        require(film.dissolution[0][a]<1e-10,"Coverage left an exposed substrate area");auto desired=film.fractions;desired[0][b]=0;film.advance(desired);
        require(std::abs(film.dissolution[0][a]-4)<1e-8,"Dissolving the covering mineral failed to expose the other mineral");
        Geometry junction(1,1,4,4,1,0);rectangle(junction,-1,a,0,2,0,2);rectangle(junction,-1,b,0,2,2,4);
        desired=junction.fractions;desired[0][a]+=.01;desired[0][b]+=.02;junction.advance(desired);
        require(std::abs(junction.fractions[0][a]-desired[0][a])<1e-9&&std::abs(junction.fractions[0][b]-desired[0][b])<1e-9,"Triple-junction growth lost a mineral volume");
        Geometry narrow(1,1,4,4,1,0);rectangle(narrow,-1,a,0,2,0,1.9998);rectangle(narrow,-1,b,0,2,2.0002,4);auto closing=narrow.fractions;closing[0][a]+=.01;closing[0][b]+=.02;narrow.advance(closing);
        require(std::abs(narrow.fractions[0][a]-closing[0][a])<1e-9&&std::abs(narrow.fractions[0][b]-closing[0][b])<1e-9&&narrow.covered[0][a]>0,"Microscopic contact failed its simultaneous conservative update");
        std::vector<std::array<double,2>> shape(junction.cells.size());for(std::size_t k=0;k<shape.size();k++)for(const auto& piece:junction.cells[k]){if(piece.phase==a)shape[k][0]+=pieceArea(piece);if(piece.phase==b)shape[k][1]+=pieceArea(piece);}
        if(reference.empty())reference=shape;else for(std::size_t k=0;k<shape.size();k++)for(int m=0;m<2;m++)require(std::abs(reference[k][m]-shape[k][m])<1e-8,"Mineral labels changed a competing triple-junction shape");
    }
    {Geometry g(1,1,4,4,1,0);rectangle(g,-1,1,0,2,0,4);std::vector<SurfaceBudget> budgets{{0,1,-.01,{}},{0,1,-.02,{}}};
     for(int k=0;k<int(g.wet.size());k++){const auto& s=g.wet[k];budgets[(s.a.y+s.b.y)*.5<2?0:1].surfaces.push_back(k);}
     auto desired=g.fractions;desired[0][1]-=.03;g.advance(desired,false,.2,&budgets);require(std::abs(g.fractions[0][1]-desired[0][1])<1e-9,"Different local water compositions double-counted their shared dissolution corner");}
    {Geometry compact(1,1,4,4,1,0);rectangle(compact,-1,0,0,2,0,2);rectangle(compact,-1,1,0,2,2,4);auto full=compact;auto target=compact.fractions;target[0][0]+=.01;target[0][1]+=.02;
     std::vector<SurfaceBudget> budgets{{0,0,.01,{},compact.bulkRegion[0]},{0,1,.02,{},compact.bulkRegion[0]}};for(int k=0;k<int(compact.wet.size());k++)budgets[compact.wet[k].mineral].surfaces.push_back(k);
     compact.advance(target,false,.2,&budgets);full.advanceImpl(target,false,.2,&budgets,true);
     double mismatch=0;for(int cell=0;cell<int(compact.cells.size());cell++)for(int m=0;m<3;m++){Cell a,b;for(const auto& p:compact.cells[cell])if(p.phase==m)a.push_back(p);for(const auto& p:full.cells[cell])if(p.phase==m)b.push_back(p);auto x=unionPieces(a),y=unionPieces(b),left=difference(x.get(),y.get()),right=difference(y.get(),x.get());double area=0;GEOSArea_r(booleanContext().handle,left.get(),&area);mismatch+=area;GEOSArea_r(booleanContext().handle,right.get(),&area);mismatch+=area;}
     require(mismatch<1e-7,"Local growth search changed the conservative competing boundary");}
    {Geometry g(1,1,4,4,1,0);rectangle(g,-1,2,0,2,0,4);auto face=g.wet.front();g.wet.clear();constexpr double width=1e-4,retreat=.05;
     for(auto ends:{std::pair<double,double>{0,2-width/2},{2-width/2,2+width/2},{2+width/2,4}}){auto s=face;s.a={2,ends.first};s.b={2,ends.second};g.wet.push_back(s);}
     std::vector<SurfaceBudget> budgets{{0,2,-retreat*(4-width)/16,{0,2}},{0,2,-retreat*width/16,{1}}};auto desired=g.fractions;desired[0][2]-=retreat/4;g.advance(desired,false,.2,&budgets);
     auto context=booleanContext().handle;auto point=owned(GEOSGeom_createPointFromXY_r(context,.98,width*.1));bool water=false;
     for(const auto& piece:g.cells[9])if(piece.phase<0){auto polygon=region(piece);water|=GEOSCovers_r(context,polygon.get(),point.get())==1;}
     require(water,"A short interface between two long patches left an artificial mineral ridge");require(g.globalLiquidComponents==1,"Uniform retreat created a trapped water artifact");}
    for(int mineral=0;mineral<3;mineral++){
        Geometry g(1,1,8,8,1,0);g.circle(mineral,2,4,1,128);g.circle(mineral,6,4,.01,32);g.rebuild();
        double perimeter=0;for(const auto& face:g.wet)perimeter+=length(face.b-face.a);SurfaceBudget budget{0,mineral,-.1*perimeter/64,{},g.wet.front().liquidRegion};for(int k=0;k<int(g.wet.size());k++)budget.surfaces.push_back(k);
        auto target=g.fractions;target[0][mineral]+=budget.change;std::vector<SurfaceBudget> budgets{budget};g.advance(target,false,.2,&budgets);
        require(std::abs(g.fractions[0][mineral]-target[0][mineral])<1e-9,"Finite exposed-crystal stocks lost the water-region reaction amount");
        for(int cell=0;cell<int(g.cells.size());cell++)if(g.origin(cell).x>=5)for(const auto& piece:g.cells[cell])require(piece.phase!=mineral,"A depleted small crystal retained artificial surface stock");
    }
    for(int mineral=0;mineral<3;mineral++){
        Geometry pocket(1,1,4,4,1,0);rectangle(pocket,-1,mineral,0,2,0,4);rectangle(pocket,mineral,-1,1.95,1.99,1,2);rectangle(pocket,mineral,-1,1.99,2,1.49,1.51);
        Cell oldLocal;double localVolume=0;std::map<int,SurfaceBudget> grouped;
        for(int r=0;r<int(pocket.waterRegions.size());r++)if(!pocket.waterRegions[r].bulk){const auto& region=pocket.waterRegions[r];localVolume+=region.fraction;for(auto [cell,k]:region.pieces){auto piece=pocket.cells[cell][k];auto o=pocket.origin(cell);for(auto& p:piece.polygon)p=p+o;for(auto& ring:piece.holes)for(auto& p:ring)p=p+o;oldLocal.push_back(std::move(piece));}}
        require(localVolume>0,"Retained growth fixture has no local reservoir");
        for(int k=0;k<int(pocket.wet.size());k++){const auto& face=pocket.wet[k];int r=face.liquidRegion;auto entry=grouped.emplace(r,SurfaceBudget{0,mineral,pocket.waterRegions[r].bulk?.001:.8*pocket.waterRegions[r].fraction,{},r});entry.first->second.surfaces.push_back(k);}
        std::vector<SurfaceBudget> localBudgets;auto target=pocket.fractions;for(auto& entry:grouped){target[0][mineral]+=entry.second.change;localBudgets.push_back(std::move(entry.second));}pocket.advance(target,false,.2,&localBudgets);
        Cell newWater;for(int cell=0;cell<int(pocket.cells.size());cell++)for(auto piece:pocket.cells[cell])if(piece.phase<0){auto o=pocket.origin(cell);for(auto& p:piece.polygon)p=p+o;for(auto& ring:piece.holes)for(auto& p:ring)p=p+o;newWater.push_back(std::move(piece));}
        auto local=unionPieces(oldLocal),water=unionPieces(newWater),remaining=intersection(local.get(),water.get());double area=0;GEOSArea_r(booleanContext().handle,remaining.get(),&area);
        require(std::abs(area/16-.2*localVolume)<1e-9,"Bulk growth stole a retained reservoir's competing-growth volume");
    }
    for(int mineral=0;mineral<3;mineral++){
        Geometry g(8,8,4,4,1,0);g.circle(mineral,16,16,8,256);g.rebuild();double totalArea=0,totalLength=0;for(const auto& a:g.fractions)totalArea+=a[mineral]*16;for(const auto& face:g.wet)totalLength+=length(face.b-face.a);
        const double apothem=8*std::cos(saltgeom::pi/256),distance=.1,factor=std::pow((apothem-distance)/apothem,2),volume=totalArea*(1-factor);
        auto desired=g.fractions;for(auto& a:desired)a[mineral]*=factor;
        std::map<int,SurfaceBudget> byWater;for(int k=0;k<int(g.wet.size());k++){const auto& face=g.wet[k];auto inserted=byWater.emplace(face.liquidRegion,SurfaceBudget{g.parent(face.solid),mineral,0,{},face.liquidRegion});auto& b=inserted.first->second;b.change-=volume*length(face.b-face.a)/totalLength/16;b.surfaces.push_back(k);}
        std::vector<SurfaceBudget> budgets;for(auto& entry:byWater)budgets.push_back(std::move(entry.second));g.advance(desired,false,.2,&budgets);
        const double radius=8-distance/std::cos(saltgeom::pi/256);for(const auto& face:g.wet)for(auto p:{face.a,face.b})require(std::abs(length(p-Point{16,16})-radius)<1e-3,"Splitting water regions distorted uniform circular retreat");
        require(g.globalLiquidComponents==1,"Uniform mineral dissolution created an isolated cavity");
        const double advance=.05,growthFactor=std::pow((apothem-distance+advance)/(apothem-distance),2);double waterArea=0,frontLength=0;for(const auto& a:g.fractions)waterArea+=a[mineral]*16;for(const auto& face:g.wet)frontLength+=length(face.b-face.a);
        auto growing=g.fractions;for(auto& a:growing)a[mineral]*=growthFactor;byWater.clear();
        for(int k=0;k<int(g.wet.size());k++){const auto& face=g.wet[k];auto entry=byWater.emplace(face.liquidRegion,SurfaceBudget{g.parent(face.water),mineral,0,{},face.liquidRegion});entry.first->second.change+=waterArea*(growthFactor-1)*length(face.b-face.a)/frontLength/16;entry.first->second.surfaces.push_back(k);}
        budgets.clear();for(auto& entry:byWater)budgets.push_back(std::move(entry.second));g.advance(growing,false,.2,&budgets);
        const double grownRadius=radius+advance/std::cos(saltgeom::pi/256);for(const auto& face:g.wet)for(auto point:{face.a,face.b})require(std::abs(length(point-Point{16,16})-grownRadius)<1e-3,"Splitting water regions distorted continuous mineral growth");
    }
    std::cout<<"PASS competing mineral fronts: first contact, all six mineral pairs, label-invariant triple junctions, local dissolution and re-exposure\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
