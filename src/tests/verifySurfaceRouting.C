// SPDX-License-Identifier: GPL-3.0-or-later
#include "MineralSurface.H"
#include <fstream>
#include <iomanip>

int main(int argc,char** argv){try{
    using namespace Foam;
    argList args(argc,argv);if(!args.checkRootCase())return 1;
    Time runTime(Time::controlDictName,args);fvMesh mesh(IOobject(fvMesh::defaultRegion,runTime.timeName(),runTime,IOobject::MUST_READ));
    IOdictionary cfg(IOobject("reactiveProperties",runTime.constant(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE));
    const double h=readScalar(cfg.lookup("dx")),ep=readScalar(cfg.lookup("epsilonMin"));
    PhaseStencil stencil(mesh,h,h);
    volScalarField inert(IOobject("inert",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE),mesh);
    volScalarField normalX(IOobject("normalX",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE),mesh);
    volScalarField normalY(IOobject("normalY",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE),mesh);
    PtrList<volScalarField> mineral(7);for(int j=0;j<7;j++)
        mineral.set(j,new volScalarField(IOobject("mineral"+name(j),runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE),mesh));
    saltflow::MineralSurface surface(cfg,mesh,stencil,h,h,ep);
    forAll(inert,i)surface.state[i].inheritedNormal[0]={normalX[i],normalY[i]};
    saltchem::RateParameters rates;rates.geometry.enabled=true;rates.neutral[0]=.02;
    surface.reconstruct(mineral,inert,rates);
    fileName database(cfg.lookup("database"));saltchem::PitzerModel model(database.c_str());saltchem::Speciation spec(model);saltchem::Reactions reactions(model);
    const auto bulk=spec.inventory({{.05,0,0,0,.05,0,0}},1.);double water;
    const auto aqueous=spec.conserved(bulk,water);saltchem::SurfaceFilmParameters raw;
    std::ofstream output((runTime.path()/"surface_routes.csv").c_str());
    output<<std::setprecision(17)<<"sx,sy,rx,ry,local,patch,area,velocity,molar_flux\n";
    auto report=[&](label sx,label sy,label rx,label ry,int local,std::size_t patch,double area,const saltchem::SurfaceFilmParameters& film){
        if(area<=0)return;
        const double velocity=film.bareVelocity[0];
        if(std::isnan(velocity)||velocity<0)throw std::runtime_error("Invalid constructed route conductance");
        const auto solved=saltchem::solveFilmDissolution(reactions,bulk,0,.02,aqueous.density*velocity,&aqueous);
        if(solved.rate<0||!std::isfinite(solved.rate)||std::abs(solved.residual)>1e-10)
            throw std::runtime_error("Constructed route fails Pitzer surface-flux closure");
        output<<sx<<','<<sy<<','<<rx<<','<<ry<<','<<local<<','<<patch<<','<<area<<','<<(std::isinf(velocity)?-1:velocity)<<','<<area*solved.rate<<'\n';
    };
    const auto& links=surface.surfaceExchange().links();
    for(std::size_t link=0;link<links.size();link++)for(std::size_t patch=0;patch<surface.facePatchCount(link);patch++){
        const auto& l=links[link];report(l.gx,l.gy,l.rx,l.ry,0,patch,surface.facePatchArea(link,patch)[0],surface.faceFilm(link,raw,patch));
    }
    forAll(inert,i)for(std::size_t patch=0;patch<surface.localPatchCount(i);patch++){
        const label x=std::llround(mesh.C()[i].x()/h-.5),y=std::llround(mesh.C()[i].y()/h-.5);
        report(x,y,x,y,1,patch,surface.localPatchArea(i,patch)[0],surface.localPatchFilm(i,raw,patch));
    }
    if(!output)throw std::runtime_error("Failed to write surface route verification output");
    Info<<"SURFACE_ROUTING_CHECK completed"<<endl;
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
