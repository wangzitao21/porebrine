// SPDX-License-Identifier: GPL-3.0-or-later
#include "CapReactions.hpp"
#include "SurfaceFilm.hpp"
#include <iostream>
#include <iomanip>
using namespace saltgeom;
using namespace saltchem;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void close(double a,double b,double tol,const char* message){require(std::abs(a-b)<=tol*std::max({std::abs(a),std::abs(b),1e-100}),message);}
int main(int argc,char** argv){try{
    require(argc==2,"Supply pitzer.dat");PitzerModel model(argv[1]);Speciation spec(model);Reactions chemistry(model);
    // Recover straight interfaces from independently integrated cell volumes,
    // over all quadrants and nearly empty/full central cells. This catches
    // spurious stair-step area without assuming the implementation's normal.
    for(int angle=0;angle<360;angle+=7)for(double target:{.03,.2,.5,.8,.97}){
        const Point normal{std::cos(angle*pi/180),std::sin(angle*pi/180)};const auto square=box(1,1);
        double lo=std::min(0.,normal.x)+std::min(0.,normal.y),hi=std::max(0.,normal.x)+std::max(0.,normal.y);
        for(int k=0;k<55;k++){const double mid=(lo+hi)/2;if(clip(square,normal,mid).area()<target)lo=mid;else hi=mid;}
        const double offset=(lo+hi)/2;std::array<double,9> f{};
        for(int y=0;y<3;y++)for(int x=0;x<3;x++)f[3*y+x]=clip(square,normal,offset-normal.x*(x-1)-normal.y*(y-1)).area();
        const auto fitted=heightFitNormal(f,1,1);
        require(dot(fitted,normal)>1-1e-11,"Planar normal depends on interface angle or volume fraction");
    }
    constexpr double A=1e-8,cellVolume=1e-13,delta=2e-7,r=1.5e-6;
    auto seeded=[&](bool film,double scale=1.){
        SurfaceState s;
        for(int j=0;j<2;j++){auto& p=s.population[7*j+2];const double share=j==0?.95:.05;
            p.fraction=A*delta*share/cellVolume;p.count=A*delta*share*scale/(capVolumeFactor(1)*r*r*r);
            p.supportArea=p.partitionSupport=A*scale;p.filmArea=film?A*share*scale:0;p.birthNormal={1,0};
        }return s;
    };
    MineralVector bare{};bare[2]=A;ProductLayer layer;layer.enabled=true;layer.porosity=0;
    auto particles=seeded(false),sheet=seeded(true);
    const auto grains=measure(particles,bare,cellVolume,true,layer),dense=measure(sheet,bare,cellVolume,true,layer);
    close(grains.coverage[2],-std::expm1(-1.5*delta/r),1e-14,"Particle coverage disagrees with Poisson sphere-cap footprint");
    close(dense.coverage[2],1,1e-14,"Tiled planar film does not fully cover its substrate");
    close(dense.layerThickness[2],delta,1e-14,"Dense film thickness is not V/A");
    layer.porosity=.35;const auto porous=measure(sheet,bare,cellVolume,true,layer);
    close(porous.layerThickness[2],delta/.65,1e-14,"Porous film thickness does not preserve solid volume");
    // Partition a physical patch without changing its volume, number density,
    // contact angle or coverage. Expected area is independent of cell count.
    for(int n:{1,2,4,16,64})for(bool film:{false,true}){
        const double scale=1./n;auto sub=seeded(film,scale);auto patch=bare;patch[2]*=scale;
        const auto m=measure(sub,patch,cellVolume*scale,true,layer);
        const auto full=measure(film?sheet:particles,bare,cellVolume,true,layer);
        close(m.coverage[2],full.coverage[2],2e-14,"Mesh subdivision changed physical coverage");
        close(m.layerThickness[2],full.layerThickness[2],2e-14,"Mesh subdivision changed physical layer thickness");
        close(m.capArea[0]*n,full.capArea[0],2e-14,"Mesh subdivision changed total product reactive area");
    }
    // At unchanged support, cap area scales as V^(2/3), planar area remains
    // constant and its thickness scales as V, until exact exhaustion.
    auto shrinking=sheet;shrinking.population[2].fraction*=.5;shrinking.population[9].fraction*=.5;
    auto thin=measure(shrinking,bare,cellVolume,true,layer);
    close(thin.coverage[2],1,1e-14,"Normal film dissolution spuriously opened holes");
    close(thin.layerThickness[2],porous.layerThickness[2]/2,1e-14,"Film thinning violates volume conservation");
    clearPopulationMaterial(shrinking.population[9]);thin=measure(shrinking,bare,cellVolume,true,layer);
    close(thin.coverage[2],.95,1e-14,"Depleted mineral patch did not open the correct area");
    auto normal=std::array<Point,7>{};normal[2]={1,0};auto half=bare;half[2]*=.5;
    updatePartition(sheet,half,normal);const auto detached=measure(sheet,half,cellVolume,true,layer);
    close(detached.attachedFraction[0]+detached.residueFraction[0],sheet.population[2].fraction,1e-14,"Support retirement lost product volume");
    updatePartition(sheet,bare,normal);close(sheet.population[2].localShare,.5,1e-14,"Growing substrate area reattached old residue");
    sheet=seeded(true);
    // All shapes below see exactly the same solution and mineral amounts.
    auto liquid=spec.inventory({{6.1292275648038,.1,0,.1,6.4292275648038,0,0}},1000,{});
    liquid.mineral[2]=.2/model.molarVolume[2];
    RateParameters rate;rate.geometry.enabled=true;rate.geometry.separateMineralAreas=true;
    rate.neutral[2]=.05;rate.area[2]=A/cellVolume;rate.referenceMoles[2]=liquid.mineral[2];
    const double bareRate=chemistry.evaluate(liquid,rate).rates[2];require(bareRate>0,"Coating fixture is not dissolving");
    std::cout<<std::setprecision(16)<<"morphology,solid_volume_m3,coverage,thickness_m,product_area_m2,carnallite_flux_ratio\n";
    for(int mode=0;mode<4;mode++){
        auto l=layer;l.porosity=mode==3?.35:0;
        const SurfaceState empty;const auto& state=mode==0?empty:mode==1?particles:sheet;
        const auto m=measure(state,bare,cellVolume,true,l);
        SurfaceFilmParameters film;film.coverage=m.coverage;
        if(m.coverage[2]>0)film.coatedVelocity[2]=l.porosity>0?l.diffusivity*l.porosity/(l.tortuosity*m.layerThickness[2]):0;
        auto p=rate;applySurfaceFilmRates(chemistry,liquid,p,film);
        const double ratio=chemistry.evaluate(liquid,p).rates[2]/bareRate;
        if(mode==0)close(ratio,1,1e-14,"Uncoated limit changed kinetics");
        if(mode==1)close(ratio,1-m.coverage[2],1e-13,"Impermeable particle contact limit is wrong");
        if(mode==2)require(ratio<1e-14,"Fully sealed film leaked substrate reaction beyond roundoff");
        if(mode==3)require(ratio>0&&ratio<1,"Porous film lost its finite transport limit");
        std::cout<<(mode==0?"baseline":mode==1?"particles":mode==2?"dense":"porous")<<','<<(mode?A*delta:0)<<','
            <<m.coverage[2]<<','<<m.layerThickness[2]<<','<<m.capArea[0]+m.capArea[1]<<','<<ratio<<'\n';
    }
    // Finite-volume chemistry must deplete a planar patch without consuming
    // the substrate, and transfer its ions through the ordinary mass ledger.
    auto film=seeded(true);Inventory cell=liquid;
    for(int j=0;j<2;j++)cell.mineral[j]=film.population[7*j+2].fraction/model.molarVolume[j];
    const auto before=cell;RateParameters p;p.geometry.enabled=true;p.neutral[1]=.01;
    IntegrationOptions opts;opts.relativeTolerance=1e-9;opts.absoluteTolerance=1e-13;
    advanceCapPopulations(chemistry,cell,film,p,measure(film,bare,cellVolume,true,layer),model.molarVolume,cellVolume,2.,298.15,opts);
    require(cell.mineral[1]<before.mineral[1]&&cell.mineral[2]==before.mineral[2],"Planar product reaction did not preserve primary stock");
    for(int i=0;i<7;i++){
        long double a=before.aqueous[i],b=cell.aqueous[i];
        for(int j=0;j<7;j++){a+=before.mineral[j]*stoichiometry[j][i];b+=cell.mineral[j]*stoichiometry[j][i];}
        close(double(a),double(b),1e-13,"Planar product reaction lost an element");
    }
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
