// SPDX-License-Identifier: GPL-3.0-or-later
#include "DualMineralGeometry.hpp"
#include <iostream>
#include <cstdio>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace saltdual;
static void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
static void rectangle(Geometry& g,int from,int to,double x0,double x1,double y0,double y1){
    for(int i=0;i<int(g.cells.size());i++){auto o=g.origin(i);paint(g.cells[i],from,to,{{{1,0},x1-o.x},{{-1,0},o.x-x0},{{0,1},y1-o.y},{{0,-1},o.y-y0}});merge(g.cells[i]);}g.rebuild();
}
static double total(const Geometry& g,int m,int kind){double v=0;for(const auto& a:kind==0?g.fractions:kind==1?g.dissolution:g.covered)v+=a[m];return v;}
int main(){
    {Poly ring;for(int k=0;k<512;k++){double angle=2*saltgeom::pi*k/512,r=.4+.015*std::cos(7*angle);ring.push_back({.5+r*std::cos(angle),.5+r*std::sin(angle)});}BoundaryEdgeIndex lookup(ring);
     for(std::size_t i=0;i<ring.size();i+=17)for(double offset:{0.,1e-10,3e-10}){auto a=ring[i],b=ring[(i+1)%ring.size()],t=(b-a)*(1/length(b-a));Point normal{t.y,-t.x};a=a+normal*offset;b=b+normal*offset;
        std::vector<std::size_t> expected,actual,candidates;Point start,end;
        for(std::size_t k=0;k<ring.size();k++)if(common(a,b,ring[k],ring[(k+1)%ring.size()],start,end)>0)expected.push_back(k);
        lookup.candidates(a,b,candidates);for(auto k:candidates)if(common(a,b,ring[k],ring[(k+1)%ring.size()],start,end)>0)actual.push_back(k);
        require(actual==expected,"Spatial edge lookup changed contact detection or ordering near the geometric tolerance");
     }}
    {Cell cell{{2,{{0,0},{1,0},{1,.6},{.75,.55004},{.5,.5},{.25,.44997},{0,.4}}}};
     double initial=area(cell[0].polygon);merge(cell);double final=0,left=0,right=0;
     for(const auto& piece:cell)if(piece.phase==2){final+=pieceArea(piece);for(std::size_t j=0;j<piece.polygon.size();j++){auto a=piece.polygon[j],b=piece.polygon[(j+1)%piece.polygon.size()];if(std::abs(a.x)<1e-10&&std::abs(b.x)<1e-10)left+=std::abs(a.y-b.y);if(std::abs(a.x-1)<1e-10&&std::abs(b.x-1)<1e-10)right+=std::abs(a.y-b.y);}}
     require(std::abs(final-initial)<1e-9,"Boundary reconstruction changed mineral volume");require(std::abs(left-.4)<1e-9&&std::abs(right-.6)<1e-9,"Boundary reconstruction changed a shared face aperture");}
    {auto c=booleanContext().handle;auto a=unionPolygons({square(),{{2,0},{3,0},{3,1},{2,1}}}),b=polygon({{1,0},{2.5,0},{2.5,1},{1,1}});auto v=intersection(a.get(),b.get());double amount=0;GEOSArea_r(c,v.get(),&amount);require(std::abs(amount-.5)<1e-12,"Polygon/line intersection changed area");require(GEOSGeomTypeId_r(c,v.get())==GEOS_POLYGON||GEOSGeomTypeId_r(c,v.get())==GEOS_MULTIPOLYGON,"Zero-volume boundary remnants reached an area overlay");}
    {auto c=booleanContext().handle;auto outer=polygon({{100,100},{104,100},{100,104}}),inner=polygon({{100.1,100.1},{100.2,100.1},{100.1,100.2}});
     auto ring=owned(GEOSDifferencePrec_r(c,outer.get(),inner.get(),1e-12));Cell parts;unpackPolygon(ring.get(),0,parts);double expected=0,actual=0;GEOSArea_r(c,ring.get(),&expected);for(auto& p:parts)actual+=pieceArea(p);require(std::abs(actual-expected)<1e-11,"Triangular exterior lost its interior liquid hole");}
    Geometry g(1,1,4,4,1,0);rectangle(g,-1,2,0,2,0,4);
    require(std::abs(total(g,2,1)-4)<1e-11,"Planar contact area");
    rectangle(g,-1,1,2,2.2,0,1.6);
    require(std::abs(total(g,2,1)-2.4)<1e-11,"Solid-solid contact was counted as reactive");
    {auto compact=g;compact.storeBuriedInterfaces=false;compact.rebuild();require(compact.buried.empty()&&compact.covered==g.covered&&compact.dissolution==g.dissolution&&compact.growth==g.growth&&compact.wet.size()==g.wet.size(),"Omitting unused buried records changed physical coverage or reaction area");}
    rectangle(g,1,-1,2,2.2,.8,1.6);
    require(std::abs(total(g,2,1)-3.2)<1e-11,"Re-exposure did not restore reactive area");
    auto target=g.fractions;target[0][2]-=.01;double product=total(g,1,0);g.advance(target);
    require(std::abs(g.fractions[0][2]-target[0][2])<2e-12,"Conservative dissolution budget");
    require(std::abs(total(g,1,0)-product)<2e-12,"Substrate retreat moved product inventory");
    require(total(g,2,2)>.7&&total(g,2,2)<=.8+1e-9,"Free-edge undercutting exceeded the interface displacement");
    Geometry exhausted(1,1,4,4,1,0);rectangle(exhausted,-1,2,0,1e-3,0,4);target=exhausted.fractions;target[0][2]=0;exhausted.advance(target);
    require(total(exhausted,2,0)<1e-11&&total(exhausted,2,1)<1e-10,"Depleted mineral retains reactive area");
    Geometry birth(1,1,4,4,1,0);rectangle(birth,-1,0,0,1,0,4);target=birth.fractions;target[0][0]+=.025;birth.advance(target);
    require(std::abs(birth.fractions[0][0]-target[0][0])<2e-12,"Growth failed across fine-cell boundary");
    Geometry film(1,1,4,4,1,0);rectangle(film,-1,2,0,2,0,4);rectangle(film,-1,0,2,2.2,0,4);
    require(total(film,2,1)<1e-11,"Complete coating failed to mask primary");
    target=film.fractions;target[0][0]-=.01;film.advance(target);
    require(total(film,2,1)<1e-11,"Thin surviving film spuriously exposes primary");
    rectangle(film,0,-1,0,4,0,4);require(std::abs(total(film,2,1)-4)<1e-11,"Removing film did not restore surface");
    Geometry crystal(1,1,4,4,1,0);crystal.circle(0,2,2,1,32);crystal.rebuild();
    for(int step=0;step<12;step++){target=crystal.fractions;target[0][0]+=.001;crystal.advance(target);require(crystal.liquidComponents[0]==1,"Growing crystal creates artificial closed water at its vertices");}
    Geometry split(1,1,4,4,1,0);rectangle(split,-1,2,1.8,2.2,0,4);require(split.liquidComponents[0]==2,"Disconnected water was silently merged");
    Geometry corner(2,2,4,4,1,0);rectangle(corner,-1,2,0,8,0,8);rectangle(corner,2,-1,0,4.012,0,4);rectangle(corner,2,-1,0,4,0,4.000001);
    target=corner.fractions;target[3][2]-=.0002;corner.advance(target);
    require(corner.globalLiquidComponents==1&&corner.isolatedLiquidArea<1e-9,"Short wet corner created a detached dissolution cavity");
    Geometry circle(8,8,4,4,1,1e-5);circle.circle(2,16,16,8);circle.rebuild();
    require(std::abs(total(circle,2,0)*16/(1-1e-5)/(saltgeom::pi*64)-1)<1.1e-4,"Circle volume");
    require(std::abs(total(circle,2,1)/(2*saltgeom::pi*8)-1)<4e-5,"Reconstructed circle perimeter");
#ifdef _OPENMP
    {auto a=circle,b=circle;auto desired=circle.fractions;for(std::size_t p=0;p<desired.size();p++)desired[p][2]=std::max(0.,desired[p][2]-.01*circle.dissolution[p][2]/circle.coarseVolume());
     omp_set_num_threads(1);a.advance(desired);omp_set_num_threads(4);b.advance(desired);
     for(std::size_t p=0;p<desired.size();p++)require(std::abs(a.fractions[p][2]-b.fractions[p][2])<1e-9&&std::abs(a.dissolution[p][2]-b.dissolution[p][2])<1e-8,"Parallel geometry changed volume or exposed area");}
#endif
    const std::string file="bin/tests/dual_geometry_roundtrip.bin";circle.write(file);auto loaded=Geometry::read(file);std::remove(file.c_str());
    require(total(circle,2,0)==total(loaded,2,0)&&total(circle,2,1)==total(loaded,2,1),"Geometry restart changed volume or area");
    Geometry holeState(1,1,1,1,1,0);rectangle(holeState,-1,2,.1,.9,.1,.9);rectangle(holeState,2,-1,.3,.7,.3,.7);holeState.write(file);auto holeLoaded=Geometry::read(file);std::remove(file.c_str());
    require(holeLoaded.globalLiquidComponents==2&&std::abs(total(holeLoaded,2,0)-.48)<1e-11&&std::abs(total(holeLoaded,2,1)-4.8)<1e-9,"Polygon-hole restart lost volume, area or liquid separation");
    {Geometry dense(1,1,1,1,1,0);dense.cells[0].clear();for(int i=0;i<600;i++){double x=double(i)/601,y=double(i+1)/601;dense.cells[0].push_back({i%2,{{x,0},{y,0},{y,1},{x,1}},{}});}dense.cells[0].push_back({-1,{{600./601,0},{1,0},{1,1},{600./601,1}},{}});dense.rebuild();dense.write(file);auto restored=Geometry::read(file,false);std::remove(file.c_str());
     require(restored.cells[0].size()==601&&restored.fractions==dense.fractions&&restored.covered==dense.covered&&restored.buried.empty(),"Dense multi-mineral geometry could not be resumed with its coverage intact");}
    {Geometry detailed(1,1,1,1,1,0);Poly solid;for(int i=0;i<20000;i++){double angle=2*saltgeom::pi*i/20000;solid.push_back({.5+.2*std::cos(angle),.5+.2*std::sin(angle)});}auto hole=solid;std::reverse(hole.begin(),hole.end());detailed.cells[0]={{0,solid,{}},{-1,square(),{hole}}};detailed.rebuild();detailed.write(file);auto restored=Geometry::read(file,false);std::remove(file.c_str());
     require(restored.cells[0][0].polygon.size()==20000&&restored.fractions==detailed.fractions&&restored.dissolution==detailed.dissolution,"A valid large saved interface failed its checkpoint round trip");}
    std::cout<<"PASS dual geometry: exact contacts, coating, re-exposure, depletion, growth, connectivity, conservation and restart\n";
}
