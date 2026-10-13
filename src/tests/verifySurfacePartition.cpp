// SPDX-License-Identifier: GPL-3.0-or-later
#include "MineralInterface.hpp"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cfenv>
using namespace saltgeom;
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
void close(double a,double b,double tolerance,const char* message){
    require(std::abs(a-b)<=tolerance*std::max({std::abs(a),std::abs(b),1e-100}),message);
}
int main(){try{
    #if defined(__GLIBC__)
    feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
    #endif
    constexpr double volume=1e-15,area=1e-10;
    MineralVector support{};support[2]=area;std::array<Point,7> normals{};normals[2]={1,0};
    ProductLayer layer;layer.enabled=layer.crossGrid=true;
    SurfaceState initial;auto& p=initial.population[2];p.fraction=.1;p.count=400;p.supportArea=area;
    updatePartition(initial,support,normals);const auto original=measure(initial,support,volume,true,layer);
    auto state=initial;state.population[2].retreat=100*capHeight(state.population[2],volume);
    updatePartition(state,support,normals);auto m=measure(state,support,volume,true,layer);
    close(m.coverage[2],original.coverage[2],1e-14,"Recession erased local film coverage");
    close(m.attachedFraction[0],original.attachedFraction[0],1e-14,"Recession reassigned local film inventory");
    auto smaller=support;smaller[2]*=.25;updatePartition(state,smaller,normals);m=measure(state,smaller,volume,true,layer);
    close(m.coverage[2],original.coverage[2],1e-14,"Support contraction changed footprint density");
    close(m.layerThickness[2],original.layerThickness[2],1e-14,"Support contraction squeezed the film");
    close(m.attachedFraction[0],.025,1e-14,"Support contraction did not split local film inventory");
    close(m.residueFraction[0],.075,1e-14,"Support contraction lost residual inventory");
    updatePartition(state,support,normals);m=measure(state,support,volume,true,layer);
    close(m.coverage[2],original.coverage[2]*.25,1e-14,"Old film spread onto newly exposed support");
    close(m.attachedFraction[0],.025,1e-14,"Support regrowth reattached old residue");
    close(m.layerThickness[2],original.layerThickness[2],1e-14,"Support regrowth diluted the old film");
    updatePartition(state,smaller,normals);m=measure(state,smaller,volume,true,layer);
    close(m.attachedFraction[0],.025,1e-14,"Revisiting a support minimum retired the same material twice");
    auto enlarged=state;enlarged.population[2].fraction*=8;
    const auto grown=measure(enlarged,smaller,volume,true,layer);
    require(grown.coverage[2]>m.coverage[2],"Actual cap growth did not enlarge its footprint");
    close(grown.attachedFraction[0],8*m.attachedFraction[0],1e-14,"Material growth broke the local/residual partition");
    updatePartition(state,MineralVector{},normals);updatePartition(state,support,normals);m=measure(state,support,volume,true,layer);
    require(m.coverage[2]==0&&m.attachedFraction[0]==0&&m.residueFraction[0]==.1,"Lost support reattached old material");
    require(m.capArea[0]>0,"Retained product lost its independent chemistry");
    // Only newly born mass can enter the local account after support returns.
    const double oldFraction=state.population[2].fraction;addPartitionVolume(state.population[2],.01,area);
    state.population[2].fraction+=.01;state.population[2].count+=1;
    updatePartition(state,support,normals);m=measure(state,support,volume,true,layer);
    close(m.attachedFraction[0],.01,1e-14,"New nucleation reattached old residual mass");
    close(m.residueFraction[0],oldFraction,1e-14,"New nucleation changed old residual mass");
    auto reverse=initial;auto reversed=normals;reversed[2]={-1,0};updatePartition(reverse,support,reversed);
    updatePartition(reverse,support,normals);require(reverse.population[2].localShare==0,"Direction reversal reattached a different interface");
    // A raw state copy is the same data carried by checkpoints and rollback.
    SurfaceState restored;std::memcpy(&restored,&state,sizeof(state));
    const auto read=measure(restored,support,volume,true,layer);
    require(read.attachedFraction==m.attachedFraction&&read.residueFraction==m.residueFraction,"Partition state was not checkpoint complete");
    auto multiple=initial;multiple.population[9]=multiple.population[2];multiple.population[9].supportArea=area*.5;
    multiple.population[9].partitionSupport=-1;multiple.population[9].fraction*=.5;
    updatePartition(multiple,support,normals);const auto before=measure(multiple,support,volume,true,layer);
    auto twice=support;twice[2]*=2;updatePartition(multiple,twice,normals);const auto after=measure(multiple,twice,volume,true,layer);
    close(before.coverage[2]*area,after.coverage[2]*2*area,1e-14,"Multiple population footprints spread when support expanded");
    close(before.layerThickness[2],after.layerThickness[2],1e-14,"Multiple population film thickness depends on fresh bare area");
    for(bool porous:{false,true}){
        auto selected=layer;selected.enabled=porous;const auto v=measure(state,support,volume,true,selected);
        close(v.attachedFraction[0]+v.residueFraction[0],state.population[2].fraction,1e-14,"Layer mode changed material partition");
    }
    auto tangent=initial;tangent.population[2].heightFactor=2;updatePartition(tangent,support,normals);
    const auto sphere=measure(tangent,support,volume,true,layer);
    require(sphere.coverage[2]==0&&sphere.attachedFraction[0]==0&&sphere.capArea[0]>0,"Zero-footprint sphere acquired local-film coverage");
    std::cout<<std::setprecision(16)<<"test,value\npartition_recession_area_history_new_birth_normals_and_restart,1\n";
    layer.diffusivity=2e-10;layer.porosity=.4;layer.tortuosity=2;const double waterD=1e-9;
    const auto pure=retainedConductivity(1,0,0,waterD,layer);require(pure.relative==1,"No product changed water conductance");
    for(double theta:{0.,.1,.4,1.}){
        const auto c=retainedConductivity(.8,.2,theta,waterD,layer);
        close(c.coverage*c.thicknessFraction*(1-c.porosity),.2,1e-14,"Retained conductivity discarded product volume");
        require(c.relative>0&&c.relative<=1&&c.coverage>=.2,"Invalid porous retained-layer conductance");
    }
    auto sealed=layer;sealed.porosity=0;
    const auto partial=retainedConductivity(.8,.2,.4,waterD,sealed),complete=retainedConductivity(.8,.2,1,waterD,sealed);
    close(partial.relative,.6,1e-14,"Sealed partial layer erased the open parallel path");
    require(complete.relative==0,"Complete impermeable layer leaked");
    const auto sparse=retainedConductivity(.999,.001,.01,waterD,sealed);
    close(sparse.relative,.99,1e-14,"Sparse impermeable particles sealed a whole cell");
    const auto filled=retainedConductivity(.2,.8,.8,waterD,layer);
    require(filled.porosity==0&&filled.thicknessFraction==1,"Crowded layer did not reduce porosity to fit its mass");
    const auto full=retainedConductivity(.8,.2,1,waterD,layer);
    const double length=.2/(1-layer.porosity),effective=layer.diffusivity*layer.porosity/layer.tortuosity;
    close(full.relative,1/(1-length+length*waterD/effective),1e-14,"Full layer failed the independent series-resistance solution");
    close(full.relativeNormal,full.relative,1e-15,"Normal conductance alias mismatch");
    close(full.relativeParallel,1-length+length*effective/waterD,1e-14,"Tangential layer failed the independent parallel-conductance solution");
    require(full.relativeParallel>=full.relativeNormal,"Laminate normal conductance exceeded its tangential conductance");
    auto legacy=layer;legacy.enabled=false;
    close(retainedConductivity(.8,.2,.4,waterD,legacy).relative,.6,1e-14,"Impermeable mode lost its parallel openings");
    require(retainedConductivity(0,1,1,waterD,layer).relative==0,"No-water cell transmitted solute");
    std::cout<<"partial_sealed_relative_conductance,"<<partial.relative<<"\nfull_porous_relative_conductance,"<<full.relative<<"\nretained_series_parallel_volume_limits,1\nall_passed,1\n";
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
