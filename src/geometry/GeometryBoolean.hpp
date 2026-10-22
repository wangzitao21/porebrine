// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <geos_c.h>
#include <memory>
#include <cfenv>
#include <exception>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace saltdual {
constexpr double overlayPrecision=1e-11;
constexpr double wetFragmentAreaFloor=1e-10;
struct BooleanContext {
    GEOSContextHandle_t handle=GEOS_init_r();std::string error;std::uint64_t fallbacks=0;
    BooleanContext(){GEOSContext_setErrorMessageHandler_r(handle,[](const char* text,void* data){static_cast<BooleanContext*>(data)->error=text;},this);}
    ~BooleanContext(){GEOS_finish_r(handle);}
};
inline BooleanContext& booleanContext(){static thread_local BooleanContext c;return c;}
struct GeosDelete {GEOSContextHandle_t context;void operator()(GEOSGeometry* g)const{if(g)GEOSGeom_destroy_r(context,g);}};
using Geos=std::unique_ptr<GEOSGeometry,GeosDelete>;
struct BooleanEnvironment {fenv_t previous;BooleanEnvironment(){feholdexcept(&previous);}~BooleanEnvironment(){fesetenv(&previous);}};
inline Geos owned(GEOSGeometry* g){auto& c=booleanContext();if(!g)throw StepFailure("Mineral polygon operation failed: "+c.error);return Geos(g,{c.handle});}
inline void collectPolygons(const GEOSGeometry* value,std::vector<GEOSGeometry*>& output){
    auto c=booleanContext().handle;const int type=GEOSGeomTypeId_r(c,value);
    if(type==GEOS_POLYGON){if(!GEOSisEmpty_r(c,value))output.push_back(GEOSGeom_clone_r(c,value));}
    else if(type==GEOS_MULTIPOLYGON||type==GEOS_GEOMETRYCOLLECTION)for(int k=0;k<GEOSGetNumGeometries_r(c,value);k++)collectPolygons(GEOSGetGeometryN_r(c,value,k),output);
}
inline Geos polygonOnly(Geos value){
    auto c=booleanContext().handle;const int type=GEOSGeomTypeId_r(c,value.get());
    if(type==GEOS_POLYGON||type==GEOS_MULTIPOLYGON)return value;
    // An area intersection can also contain a touching line/point. Those
    // zero-volume remnants are not mineral stock or additional interfaces.
    std::vector<GEOSGeometry*> polygons;collectPolygons(value.get(),polygons);
    if(polygons.empty())return owned(GEOSGeom_createEmptyPolygon_r(c));
    return owned(GEOSGeom_createCollection_r(c,GEOS_MULTIPOLYGON,polygons.data(),polygons.size()));
}
inline GEOSGeometry* linearRing(const Poly& points){
    auto c=booleanContext().handle;auto* seq=GEOSCoordSeq_create_r(c,points.size()+1,2);
    for(std::size_t k=0;k<=points.size();k++){const auto& p=points[k%points.size()];GEOSCoordSeq_setX_r(c,seq,k,p.x);GEOSCoordSeq_setY_r(c,seq,k,p.y);}
    return GEOSGeom_createLinearRing_r(c,seq);
}
inline Geos rawRegion(const Piece& piece){
    auto c=booleanContext().handle;if(piece.polygon.empty())return owned(GEOSGeom_createEmptyPolygon_r(c));
    std::vector<GEOSGeometry*> holes;for(const auto& h:piece.holes)holes.push_back(linearRing(h));
    return owned(GEOSGeom_createPolygon_r(c,linearRing(piece.polygon),holes.data(),holes.size()));
}
inline Geos region(const Piece& piece){
    auto c=booleanContext().handle;auto result=rawRegion(piece);
    if(GEOSisValid_r(c,result.get())!=1){booleanContext().fallbacks++;result=polygonOnly(owned(GEOSMakeValid_r(c,result.get())));}
    auto snapped=GEOSGeom_setPrecision_r(c,result.get(),overlayPrecision,GEOS_PREC_VALID_OUTPUT);
    if(snapped)return polygonOnly(owned(snapped));
    booleanContext().fallbacks++;return polygonOnly(std::move(result));
}
inline Geos polygon(const Poly& points){return region(Piece{0,points,{}});}
inline Geos unionPieces(const Cell& pieces){
    auto c=booleanContext().handle;std::vector<GEOSGeometry*> input;for(const auto& p:pieces)input.push_back(region(p).release());
    auto collection=owned(GEOSGeom_createCollection_r(c,GEOS_GEOMETRYCOLLECTION,input.data(),input.size()));
    auto* result=GEOSUnaryUnionPrec_r(c,collection.get(),overlayPrecision);if(!result){booleanContext().fallbacks++;result=GEOSUnaryUnion_r(c,collection.get());}
    return polygonOnly(owned(result));
}
inline Geos unionPolygons(const std::vector<Poly>& polygons){
    auto c=booleanContext().handle;std::vector<GEOSGeometry*> input;
    for(const auto& p:polygons)input.push_back(polygon(p).release());
    auto collection=owned(GEOSGeom_createCollection_r(c,GEOS_GEOMETRYCOLLECTION,input.data(),input.size()));
    auto* result=GEOSUnaryUnionPrec_r(c,collection.get(),overlayPrecision);
    if(!result){booleanContext().fallbacks++;result=GEOSUnaryUnion_r(c,collection.get());}
    return polygonOnly(owned(result));
}
inline Geos intersection(const GEOSGeometry* a,const GEOSGeometry* b){auto c=booleanContext().handle;auto* r=GEOSIntersectionPrec_r(c,a,b,overlayPrecision);if(!r){booleanContext().fallbacks++;auto x=polygonOnly(owned(GEOSMakeValid_r(c,a))),y=polygonOnly(owned(GEOSMakeValid_r(c,b)));r=GEOSIntersectionPrec_r(c,x.get(),y.get(),overlayPrecision);if(!r)r=GEOSIntersection_r(c,x.get(),y.get());}return polygonOnly(owned(r));}
inline Geos difference(const GEOSGeometry* a,const GEOSGeometry* b){auto c=booleanContext().handle;auto* r=GEOSDifferencePrec_r(c,a,b,overlayPrecision);if(!r){booleanContext().fallbacks++;auto x=polygonOnly(owned(GEOSMakeValid_r(c,a))),y=polygonOnly(owned(GEOSMakeValid_r(c,b)));r=GEOSDifferencePrec_r(c,x.get(),y.get(),overlayPrecision);if(!r)r=GEOSDifference_r(c,x.get(),y.get());}return polygonOnly(owned(r));}
inline Geos unite(const GEOSGeometry* a,const GEOSGeometry* b){auto c=booleanContext().handle;auto* r=GEOSUnionPrec_r(c,a,b,overlayPrecision);if(!r){booleanContext().fallbacks++;auto x=polygonOnly(owned(GEOSMakeValid_r(c,a))),y=polygonOnly(owned(GEOSMakeValid_r(c,b)));r=GEOSUnionPrec_r(c,x.get(),y.get(),overlayPrecision);if(!r)r=GEOSUnion_r(c,x.get(),y.get());}return polygonOnly(owned(r));}
inline Poly readRing(const GEOSGeometry* ring,bool hole,bool cleanup=true){
    auto c=booleanContext().handle;auto seq=GEOSGeom_getCoordSeq_r(c,ring);unsigned count=0;GEOSCoordSeq_getSize_r(c,seq,&count);Poly p;
    for(unsigned k=0;k+1<count;k++){Point v;GEOSCoordSeq_getX_r(c,seq,k,&v.x);GEOSCoordSeq_getY_r(c,seq,k,&v.y);p.push_back(v);}
    if(cleanup)p=clean(std::move(p));if(p.empty())return p;
    double a=0;for(std::size_t k=1;k+1<p.size();k++)a+=cross(p[k]-p[0],p[k+1]-p[0]);
    if((a>0)==hole)std::reverse(p.begin(),p.end());return p;
}
inline void unpackPolygon(const GEOSGeometry* g,int phase,Cell& result,bool cleanup=true){
    auto c=booleanContext().handle;int type=GEOSGeomTypeId_r(c,g);
    if(type==GEOS_POLYGON){
        auto ring=GEOSGetExteriorRing_r(c,g);if(!ring)return;Piece p{phase,readRing(ring,false,cleanup),{}};if(p.polygon.empty())return;
        for(int k=0;k<GEOSGetNumInteriorRings_r(c,g);k++){auto h=readRing(GEOSGetInteriorRingN_r(c,g,k),true,cleanup);if(!h.empty())p.holes.push_back(std::move(h));}
        result.push_back(std::move(p));
    }else if(type==GEOS_MULTIPOLYGON||type==GEOS_GEOMETRYCOLLECTION){
        const int n=GEOSGetNumGeometries_r(c,g);for(int k=0;k<n;k++)unpackPolygon(GEOSGetGeometryN_r(c,g,k),phase,result,cleanup);
    }
}
inline double paint(Cell& cell,int from,int to,const std::vector<std::pair<Point,double>>& planes){
    BooleanEnvironment environment;auto c=booleanContext().handle;auto mask=square();
    for(const auto& plane:planes){mask=clip(mask,plane.first,plane.second);if(mask.empty())return 0;}
    Cell source,destination,result;for(const auto& p:cell){if(p.phase==from)source.push_back(p);else if(p.phase==to)destination.push_back(p);else result.push_back(p);}
    if(source.empty())return 0;
    auto a=unionPieces(source),b=polygon(mask),moved=intersection(a.get(),b.get());double amount=0;GEOSArea_r(c,moved.get(),&amount);if(amount<1e-16)return 0;
    auto remainder=difference(a.get(),b.get());unpackPolygon(remainder.get(),from,result);
    if(!destination.empty()){auto old=unionPieces(destination);moved=unite(old.get(),moved.get());}
    unpackPolygon(moved.get(),to,result);cell=std::move(result);return amount;
}
inline void merge(Cell& cell,std::array<double,3>* removedArea){
    // GEOS internally uses non-finite sentinels. Isolate its exception flags
    // from OpenFOAM's optional SIGFPE traps; validate the returned partition
    // and conservative volume explicitly in the geometry layer.
    BooleanEnvironment environment;
    auto c=booleanContext().handle;Cell result;
    std::array<Geos,3> phases{{Geos(nullptr,{c}),Geos(nullptr,{c}),Geos(nullptr,{c})}};
    for(int phase=0;phase<3;phase++){
        std::vector<GEOSGeometry*> input;
        for(const auto& p:cell)if(p.phase==phase&&!p.polygon.empty())input.push_back(region(p).release());
        if(input.empty())continue;
        auto collection=owned(GEOSGeom_createCollection_r(c,GEOS_GEOMETRYCOLLECTION,input.data(),input.size()));
        auto* raw=GEOSUnaryUnionPrec_r(c,collection.get(),overlayPrecision);
        if(!raw){booleanContext().fallbacks++;raw=GEOSUnaryUnion_r(c,collection.get());}
        auto united=polygonOnly(owned(raw));
        phases[phase]=std::move(united);
    }
    // Construct the exact swept phase partition shared by all minerals.
    auto liquid=polygon(square());
    for(int phase=0;phase<3;phase++)if(phases[phase])liquid=difference(liquid.get(),phases[phase].get());
    Cell unresolved;
    for(int phase=0;phase<3;phase++)if(phases[phase]){
        Cell pieces;unpackPolygon(phases[phase].get(),phase,pieces);
        for(auto& piece:pieces){double amount=pieceArea(piece);bool discard=false;
            if(amount<wetFragmentAreaFloor){
                auto solid=region(piece);double distance=0;GEOSDistance_r(c,solid.get(),liquid.get(),&distance);
                if(distance<=tiny){auto raw=GEOSIntersection_r(c,solid.get(),liquid.get());if(raw){auto contact=owned(raw);double length=0;GEOSLength_r(c,contact.get(),&length);discard=length>tiny;}}
            }
            if(discard){if(removedArea)(*removedArea)[phase]+=amount;unresolved.push_back(std::move(piece));}
            else result.push_back(std::move(piece));
        }
    }
    if(!unresolved.empty()){auto extra=unionPieces(unresolved);liquid=unite(liquid.get(),extra.get());}
    unpackPolygon(liquid.get(),-1,result);cell=std::move(result);
}
struct CellEdit {int cell;std::array<std::vector<Cell>,3> changes;};
inline std::vector<std::array<double,3>> commitCellEdits(std::vector<Cell>& cells,const std::vector<CellEdit>& edits,bool growing){
    std::vector<std::array<double,3>> removed(edits.size());std::vector<std::exception_ptr> errors(edits.size());std::vector<std::uint64_t> repairs(edits.size());
    const auto before=booleanContext().fallbacks;
    bool parallel=false;
#ifdef _OPENMP
    parallel=!omp_in_parallel()&&edits.size()>16;
#endif
    // Each worker owns its GEOS context and all geometry objects. Only plain
    // polygon coordinates cross the worker boundary; cell edits are disjoint.
    #pragma omp parallel for schedule(dynamic,1) if(parallel)
    for(int i=0;i<int(edits.size());i++){
        auto& context=booleanContext();const auto previous=context.fallbacks;
        try{
            BooleanEnvironment environment;const auto& edit=edits[i];Cell revised;
            for(int m=0;m<3;m++){
                Cell old;for(const auto& piece:cells[edit.cell])if(piece.phase==m)old.push_back(piece);
                if(edit.changes[m].empty()){revised.insert(revised.end(),old.begin(),old.end());continue;}
                auto material=unionPieces(old);
                for(const auto& pieces:edit.changes[m]){
                    std::vector<GEOSGeometry*> parts;for(const auto& piece:pieces)parts.push_back(rawRegion(piece).release());
                    auto cut=owned(GEOSGeom_createCollection_r(context.handle,GEOS_MULTIPOLYGON,parts.data(),parts.size()));
                    material=growing?unite(material.get(),cut.get()):difference(material.get(),cut.get());
                }
                unpackPolygon(material.get(),m,revised);
            }
            merge(revised,&removed[i]);
            if(revised.size()>maxFineCellPieces)throw std::runtime_error("Fine-cell polygon storage capacity exceeded");
            cells[edit.cell]=std::move(revised);
        }catch(...){errors[i]=std::current_exception();}
        repairs[i]=context.fallbacks-previous;
    }
    booleanContext().fallbacks=before+std::accumulate(repairs.begin(),repairs.end(),std::uint64_t(0));
    for(const auto& error:errors)if(error)std::rethrow_exception(error);
    return removed;
}

}
