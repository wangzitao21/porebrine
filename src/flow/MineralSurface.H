// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "MineralInterface.hpp"
#include "SurfaceReaction.hpp"
#include "CapReactions.hpp"
#include "PhaseTransform.H"
#include "DiffusiveConductance.H"
#include "FaceReactionExchange.H"
#include "SurfaceLinkExchange.H"
#include "labelIOList.H"
#include <fstream>
#include <iomanip>
#include <memory>

namespace saltflow {
class SurfaceTransportModelError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};
class MineralSurface {
    const Foam::fvMesh& mesh;const PhaseStencil& stencil;
    double dx,dy,epsilonMin,resolvedRadius,surfaceTimeTolerance,molecularDiffusivity;
    bool blocking;std::uint64_t seed;
    saltgeom::ProductLayer productLayer;
    std::vector<saltgeom::NucleationLaw> laws;
    std::vector<std::uint64_t> keys;
    std::vector<saltgeom::Reconstruction> shape;
    std::vector<saltchem::MineralVector> bare,incoming;
    std::vector<saltgeom::SurfaceMeasure> measures;
    Foam::volScalarField exchange;
    std::vector<std::array<double,2>> transportResistance;
    std::vector<std::array<double,2>> transportRelative;
    bool retainedTransportActive=false;
    std::vector<double> transportFraction;
    std::vector<std::array<bool,4>> connected;
    std::vector<std::array<Foam::label,4>> directionalFace;
    std::vector<std::array<Foam::label,2>> cellIndex;
    std::array<Foam::label,2> gridExtent{};
    std::vector<Foam::label> coupledOffsets;
    struct RouteEdge {std::array<Foam::label,3> key{};double resistance=0;};
    struct ReactionPatch {
        saltchem::MineralVector area{};
        std::array<saltgeom::Point,7> normal{};
        std::array<RouteEdge,2> route{};
        int hops=0;
        double waterResistance=0; // Unit-patch-area resistance, before division by D.
    };
    struct FaceContact {
        saltchem::MineralVector area{},support{};
        std::array<saltgeom::Point,7> normal{};
        std::vector<ReactionPatch> patches;
    };
    std::array<bool,7> routeReactive{};
    std::vector<FaceContact> contacts;
    std::vector<SurfaceLink> surfaceLinks;
    Foam::scalarField faceOpening;
    std::vector<saltchem::MineralVector> localArea,support;
    std::vector<saltchem::MineralVector> directLocalArea;
    std::vector<std::vector<ReactionPatch>> localRoutes;
    std::vector<std::array<saltgeom::Point,7>> reactionNormal;
    std::unique_ptr<FaceReactionExchange> meshExchange;
    std::unique_ptr<SurfaceLinkExchange> faceExchange;
    Foam::word solidOwnerSurfacePolicy;
    Foam::label sharpSurfaceCells_=0;
    double sharpSurfaceArea_=0;
    bool reportedSolidOwnerSurface=false;
public:
    bool enabled;
    std::vector<saltgeom::SurfaceState> state;
    std::vector<saltchem::MineralVector> reactionRates;
    std::vector<std::array<double,2>> reactionVectors;
    MineralSurface(const Foam::dictionary& cfg,const Foam::fvMesh& m,const PhaseStencil& st,double x,double y,double ep):
        mesh(m),stencil(st),dx(x),dy(y),epsilonMin(ep),resolvedRadius(cfg.lookupOrDefault<Foam::scalar>("resolvedNucleusRadius",.35*std::min(x,y))),
        surfaceTimeTolerance(cfg.lookupOrDefault<Foam::scalar>("timeSurfaceTolerance",1e-3)),
        molecularDiffusivity(Foam::dimensionedScalar("diffusivity",Foam::dimViscosity,cfg).value()),
        blocking(cfg.lookupOrDefault<bool>("surfaceBlocking",true)),seed(cfg.lookupOrDefault<Foam::label>("nucleationSeed",2024)),
        exchange(Foam::IOobject("surfaceExchange",m.time().timeName(),m),m,Foam::dimensionedScalar(Foam::dimless,0),Foam::zeroGradientFvPatchScalarField::typeName),
        solidOwnerSurfacePolicy(cfg.lookupOrDefault<Foam::word>("solidOwnerSurfacePolicy","fluidSide")),
        enabled(cfg.lookupOrDefault<Foam::word>("mineralInterfaceModel","phaseVector")=="heterogeneous") {
        const auto mode=cfg.lookupOrDefault<Foam::word>("mineralInterfaceModel","phaseVector");
        if(mode!="phaseVector"&&mode!="heterogeneous")throw std::invalid_argument("Unknown mineralInterfaceModel");
        if(solidOwnerSurfacePolicy!="fluidSide"&&solidOwnerSurfacePolicy!="volumeAveraged"&&solidOwnerSurfacePolicy!="reject")
            throw std::invalid_argument("solidOwnerSurfacePolicy must be fluidSide, volumeAveraged or reject");
        if(!enabled)return;
        const auto layerMode=cfg.lookupOrDefault<Foam::word>("productLayerModel","impermeable");
        if(layerMode!="impermeable"&&layerMode!="porous"&&layerMode!="crossGrid")throw std::invalid_argument("Unknown productLayerModel");
        productLayer.enabled=layerMode!="impermeable";productLayer.crossGrid=layerMode=="crossGrid";
        productLayer.diffusivity=cfg.lookupOrDefault<Foam::scalar>("productLayerWaterDiffusivity",1e-9);
        productLayer.porosity=cfg.lookupOrDefault<Foam::scalar>("productLayerPorosity",.35);
        productLayer.tortuosity=cfg.lookupOrDefault<Foam::scalar>("productLayerTortuosity",2);
        productLayer.affinityConcentration=cfg.lookupOrDefault<Foam::scalar>("productLayerAffinityConcentration",1000);
        productLayer.validate();
        if(!(resolvedRadius>0)||!std::isfinite(resolvedRadius))throw std::invalid_argument("Invalid resolvedNucleusRadius");
        if(!std::isfinite(surfaceTimeTolerance)||surfaceTimeTolerance<=0||surfaceTimeTolerance>1)throw std::invalid_argument("Invalid timeSurfaceTolerance");
        state.resize(m.nCells());shape.resize(m.nCells());bare.resize(m.nCells());incoming.resize(m.nCells());measures.resize(m.nCells());keys.resize(m.nCells());
        transportResistance.resize(m.nCells());transportRelative.resize(m.nCells(),{{1,1}});transportFraction.resize(m.nCells());
        // Geometry changes phase fractions, never Cartesian connectivity.
        connected.resize(m.nCells());
        directionalFace.resize(m.nCells());for(auto& f:directionalFace)f.fill(-1);
        forAll(m.neighbour(),f){const int d=stencil.directions[f];connected[m.owner()[f]][d]=true;connected[m.neighbour()[f]][d^1]=true;}
        forAll(m.neighbour(),f){const int d=stencil.directions[f];directionalFace[m.owner()[f]][d]=f;directionalFace[m.neighbour()[f]][d^1]=f;}
        forAll(m.boundary(),p)if(m.boundary()[p].coupled())forAll(m.boundary()[p],f)
        {
            connected[m.boundary()[p].faceCells()[f]][stencil.patchDirections[p][f]]=true;
            directionalFace[m.boundary()[p].faceCells()[f]][stencil.patchDirections[p][f]]=m.boundary()[p].start()+f;
        }
        if(fluidSide()){
            reactionRates.resize(m.nCells());reactionVectors.resize(m.nCells());
            faceOpening.setSize(m.nFaces(),0);localArea.resize(m.nCells());directLocalArea.resize(m.nCells());localRoutes.resize(m.nCells());support.resize(m.nCells());reactionNormal.resize(m.nCells());
            meshExchange=std::make_unique<FaceReactionExchange>(m,stencil,dx,dy);
            cellIndex.resize(m.nCells());
            const Foam::scalarField x(m.C().primitiveField().component(0)),y(m.C().primitiveField().component(1));
            const double x0=Foam::gMin(x),y0=Foam::gMin(y);
            gridExtent={{Foam::label(std::llround((Foam::gMax(x)-x0)/dx))+1,Foam::label(std::llround((Foam::gMax(y)-y0)/dy))+1}};
            forAll(m.C(),i)cellIndex[i]={{Foam::label(std::llround((m.C()[i].x()-x0)/dx)),Foam::label(std::llround((m.C()[i].y()-y0)/dy))}};
            coupledOffsets.resize(m.boundary().size(),-1);Foam::label offset=0;
            forAll(m.boundary(),p)if(m.boundary()[p].coupled()){coupledOffsets[p]=offset;offset+=m.boundary()[p].size();}
        }
        // Coordinates, not local cell order or rank, key the random stream.
        forAll(m.C(),i){const auto c=m.C()[i];keys[i]=saltgeom::mix(std::uint64_t(std::llround(c.x()/dx-.5)))^saltgeom::mix(std::uint64_t(std::llround(c.y()/dy-.5))+0x19f3);}
        if(cfg.found("heterogeneousNucleation"))for(const auto& entry:cfg.subDict("heterogeneousNucleation")){
            const auto& d=entry.dict();saltgeom::NucleationLaw law;
            law.product=saltchem::mineralIndex(Foam::word(d.lookup("product")).c_str());law.substrate=saltchem::mineralIndex(Foam::word(d.lookup("substrate")).c_str());
            law.prefactor=Foam::readScalar(d.lookup("prefactor"));
            const auto nucleationModel=d.lookupOrDefault<Foam::word>("nucleationModel",d.found("interfacialEnergyJm2")?"classical":"effectiveBarrier");
            if(nucleationModel!="classical"&&nucleationModel!="effectiveBarrier")throw std::invalid_argument("nucleationModel must be classical or effectiveBarrier");
            law.physicalCNT=nucleationModel=="classical";
            if(law.physicalCNT){
                if(!fluidSide())throw std::invalid_argument("Classical nucleation and curved cap growth require solidOwnerSurfacePolicy fluidSide");
                law.interfacialEnergy=Foam::readScalar(d.lookup("interfacialEnergyJm2"));
                if(d.found("barrier"))throw std::invalid_argument("Classical nucleation derives its barrier from interfacialEnergyJm2; do not also supply barrier");
            }else{
                if(d.found("interfacialEnergyJm2"))throw std::invalid_argument("effectiveBarrier mode does not infer a Kelvin energy; use classical with interfacialEnergyJm2");
                law.barrier=Foam::readScalar(d.lookup("barrier"));
            }
            law.contactAngle=d.lookupOrDefault<Foam::scalar>("contactAngle",90);law.radius=d.lookupOrDefault<Foam::scalar>("embryoRadius",2e-8);law.validate();
            for(const auto& other:laws)if(other.product==law.product&&other.substrate==law.substrate)throw std::invalid_argument("Duplicate nucleation pair");
            laws.push_back(law);
        }
        if(m.time().value()!=0)checkpoint(false);
    }
    Foam::label sharpSurfaceCells()const {return sharpSurfaceCells_;}
    double sharpSurfaceArea()const {return sharpSurfaceArea_;}
    bool fluidSide()const {return enabled&&solidOwnerSurfacePolicy=="fluidSide";}
    SurfaceLinkExchange& surfaceExchange(){return *faceExchange;}
    const saltchem::MineralVector& supportArea(int i)const {return fluidSide()?support[i]:bare[i];}
    saltchem::MineralVector lockedCaps(int i,const saltchem::MineralVector& volume)const {
        auto amount=saltgeom::subgrid(state[i]);for(int j=0;j<7;j++)amount[j]/=volume[j];return amount;
    }
    const saltchem::MineralVector& localResolvedArea(int i)const {return localArea[i];}
    saltgeom::SurfaceMeasure currentMeasure(int i)const {
        return saltgeom::measure(state[i],supportArea(i),mesh.V()[i],blocking,productLayer);
    }
    saltgeom::ProductLayer localFilmLayer(int i,double waterFraction)const {
        auto layer=productLayer;
        if(layer.enabled){
            const auto caps=saltgeom::subgrid(state[i]);double product=0;for(double f:caps)product+=f;
            // Allocate available pore water to product populations in
            // proportion to their volume. Local and retained material cannot
            // each claim the same full pore-water volume for its porosity.
            if(product>0)layer.porosity=std::min(layer.porosity,std::max(0.,waterFraction)/(std::max(0.,waterFraction)+product));
        }else layer.porosity=0;
        return layer;
    }
    saltchem::SurfaceFilmParameters filmParameters(int i,double waterFraction)const {
        saltchem::SurfaceFilmParameters film;if(!blocking)return film;
        const auto layer=localFilmLayer(i,waterFraction);
        const auto measured=saltgeom::measure(state[i],supportArea(i),mesh.V()[i],blocking,layer);film.coverage=measured.layerCoverage;
        for(int j=0;j<7;j++)if(film.coverage[j]>0){
            const double thickness=measured.layerThickness[j];
            film.coatedVelocity[j]=!layer.enabled||layer.diffusivity==0||layer.porosity==0?0:
                thickness>0?layer.diffusivity*layer.porosity/(layer.tortuosity*thickness):INFINITY;
        }
        return film;
    }
    void commitPrimary(int i,const saltchem::SurfaceReactionResult& result,const saltchem::MineralVector& volume,std::size_t patch){
        saltchem::MineralVector retreat{};
        for(int j=0;j<7;j++){
            const double fraction=result.dissolved[j]*volume[j];
            const auto normal=localPatchNormal(i,patch,j);
            reactionVectors[i][0]+=fraction*normal.x;
            reactionVectors[i][1]+=fraction*normal.y;
            if(support[i][j]>0)retreat[j]=fraction*mesh.V()[i]/(support[i][j]*(1-epsilonMin));
        }
        retreatExternal(i,retreat);
    }
    void commitCaps(int i,const saltchem::CapReactionResult& result){
        for(int j=0;j<7;j++)for(int k=0;k<7;k++){
            auto normal=reactionNormal[i][k];if(saltgeom::length(normal)==0)normal=state[i].inheritedNormal[j];
            const double change=result.fractionChange[7*j+k];
            reactionVectors[i][0]-=change*normal.x;reactionVectors[i][1]-=change*normal.y;
        }
        saltgeom::promoteResolvedPopulations(state[i],mesh.V()[i],resolvedRadius,reactionNormal[i]);
    }
    void configureResolvedGrowth(int i,const saltchem::Inventory& cell,saltchem::RateParameters& r,std::size_t patch)const {
        r.geometry.separateMineralAreas=true;r.area.fill(0);
        r.geometry.lockedMoles=cell.mineral;r.geometry.surfaceExponent.fill(0);r.geometry.curvatureLogOmega.fill(0);
        for(int j=0;j<7;j++){
            r.geometry.growthArea[j]=localPatchArea(i,patch)[j]/mesh.V()[i];r.referenceMoles[j]=std::max(cell.mineral[j],1e-30);
        }
    }
    void commitResolvedGrowth(int i,const saltchem::Inventory& before,const saltchem::Inventory& after,const saltchem::MineralVector& volume,std::size_t patch){
        for(int j=0;j<7;j++){
            const double change=(after.mineral[j]-before.mineral[j])*volume[j];
            const auto normal=localPatchNormal(i,patch,j);
            reactionVectors[i][0]-=change*normal.x;reactionVectors[i][1]-=change*normal.y;
            if(change>0&&saltgeom::length(normal)>0)state[i].inheritedNormal[j]=normal;
        }
    }
    void updateFilmDiagnostics(int i,const saltchem::Reactions& chemistry,const saltchem::Inventory& cell,const saltchem::RateParameters& kinetics,double waterFraction){
        if(!fluidSide())return;
        auto parameters=kinetics;parameters.geometry.separateMineralAreas=true;
        parameters.geometry.lockedMoles.fill(0);parameters.geometry.curvatureLogOmega.fill(0);
        for(int j=0;j<7;j++)parameters.area[j]=parameters.geometry.growthArea[j]=support[i][j]/mesh.V()[i];
        const auto film=filmParameters(i,waterFraction);
        const auto result=saltchem::applySurfaceFilmRates(chemistry,cell,parameters,film);
        auto& measured=measures[i];measured=saltgeom::measure(state[i],supportArea(i),mesh.V()[i],blocking,localFilmLayer(i,waterFraction));
        for(int j=0;j<7;j++){
            measured.layerTransmission[j]=result.filmTransmission[j];
            measured.accessibility[j]=1-film.coverage[j]+film.coverage[j]*result.filmTransmission[j];
            measured.exposed[j]=support[i][j]*measured.accessibility[j]+measured.capArea[j];
        }
    }
    bool crossGridTransport()const {return enabled&&(productLayer.crossGrid||fluidSide());}
    // Refresh after a chemical segment and after restart/reconstruction. A
    // retained population stays in its original cell even with no substrate.
    // Its disjoint attached part is handled by the local interface film.
    void updateTransport(const Foam::volScalarField& eps){
        if(!crossGridTransport())return;
        int changedTransport=0;
        forAll(eps,i){
            double retained=0,totalProduct=0,xx=0,yy=0,opticalVolume=0;
            for(int j=0;j<7;j++)for(int k=0;k<7;k++){
                const auto& p=state[i].population[7*j+k];
                totalProduct+=p.fraction;
                const double residue=1-saltgeom::attachedShare(p,supportArea(i)[k],mesh.V()[i]);
                const double f=p.fraction*residue;
                if(f<=0)continue;
                auto n=p.birthNormal;if(saltgeom::length(n)==0)n=state[i].inheritedNormal[j];const double norm=saltgeom::length(n);
                retained+=f;
                if(norm>0){
                    xx+=f*n.x*n.x/(norm*norm);yy+=f*n.y*n.y/(norm*norm);
                    opticalVolume+=saltgeom::footprintArea(p,mesh.V()[i])*residue*(std::abs(n.x)*dx+std::abs(n.y)*dy)/norm;
                }else{
                    xx+=.5*f;yy+=.5*f;opticalVolume+=saltgeom::footprintArea(p,mesh.V()[i])*residue*.5*(dx+dy);
                }
            }
            transportFraction[i]=retained;
            double water=eps[i],product=retained;
            if(fluidSide()){
                const auto layerProperties=localFilmLayer(i,water);
                const double localWater=(totalProduct-retained)*layerProperties.porosity/(1-layerProperties.porosity);
                water=std::max(0.,water-localWater);
                const double region=water+retained;
                const double coverage=region>0?-std::expm1(-opticalVolume/(region*mesh.V()[i])):0;
                const auto layer=saltgeom::retainedConductivity(water,retained,coverage,molecularDiffusivity,layerProperties);
                // Diagonal tensor closure: exact for grid-aligned laminates.
                // For oblique populations the off-diagonal transport term is
                // unresolved; normals never cancel in these second moments.
                transportResistance[i]={};transportRelative[i]={{1,1}};
                if(blocking&&retained>0)for(int axis=0;axis<2;axis++){
                    const double normalWeight=(axis==0?xx:yy)/retained;
                    const double relative=normalWeight*layer.relativeNormal+(1-normalWeight)*layer.relativeParallel;
                    transportRelative[i][axis]=relative;
                    changedTransport|=relative!=1;
                    transportResistance[i][axis]=relative>1e-100?1/relative-1:1e100;
                }
                continue;
            }
            const double resistance=blocking?saltgeom::retainedLayerResistance(water,product,productLayer):0;
            transportResistance[i]=retained>0?std::array<double,2>{{resistance*xx/retained,resistance*yy/retained}}:std::array<double,2>{};
        }
        if(fluidSide()){Foam::reduce(changedTransport,Foam::maxOp<int>());retainedTransportActive=changedTransport;}
    }
    // One shared positive coefficient per face is used by the solve, residual,
    // boundary conditions and mass ledger. Half-cell excess resistances add in
    // series; the untouched FV water cells transmit through gaps and around
    // openings. There is no nearest-parent remapping or nonlocal mass source.
    Foam::tmp<Foam::surfaceScalarField> diffusivity(const Foam::dimensionedScalar& D,const Foam::volScalarField& eps,const Foam::volScalarField& rho){
        if(fluidSide()){
            if(retainedTransportActive){forAll(exchange,i)exchange[i]=transportRelative[i][0];exchange.correctBoundaryConditions();}
            auto result=retainedTransportActive?diffusiveConductance((D*rho*exchange)()):diffusiveConductance((D*rho)());auto& gamma=result.ref();
            auto setDirection=[&](int axis,const Foam::surfaceScalarField& directional){
                forAll(mesh.neighbour(),f)if(stencil.directions[f]/2==axis)gamma[f]=directional[f]*faceOpening[f];
                forAll(mesh.boundary(),patch)forAll(gamma.boundaryField()[patch],f)
                    if(mesh.boundary()[patch].type()=="empty")gamma.boundaryFieldRef()[patch][f]=0;
                    else if(stencil.patchDirections[patch][f]/2==axis)
                        gamma.boundaryFieldRef()[patch][f]=directional.boundaryField()[patch][f]*faceOpening[mesh.boundary()[patch].start()+f];
            };
            setDirection(0,gamma);
            if(retainedTransportActive){
                forAll(exchange,i)exchange[i]=transportRelative[i][1];exchange.correctBoundaryConditions();
                const auto directional=diffusiveConductance((D*rho*exchange)());setDirection(1,directional());
            }else setDirection(1,gamma);
            return result;
        }
        auto result=diffusiveConductance((D*Foam::sqr(eps)*rho)());
        auto& gamma=result.ref();
        if(crossGridTransport()&&blocking)for(int axis=0;axis<2;axis++){
            forAll(exchange,i)exchange[i]=transportResistance[i][axis]>0?
                transportResistance[i][axis]/(std::max(productLayer.diffusivity,1e-100)*rho[i]):0;
            exchange.correctBoundaryConditions();
            const Foam::surfaceScalarField additional(Foam::fvc::interpolate(exchange));
            forAll(mesh.neighbour(),f)if(stencil.directions[f]/2==axis)
                gamma[f]=saltgeom::faceDiffusion(gamma[f],additional[f],additional[f]);
            forAll(mesh.boundary(),patch)if(mesh.boundary()[patch].type()!="empty")forAll(gamma.boundaryField()[patch],f)if(stencil.patchDirections[patch][f]/2==axis){
                auto& value=gamma.boundaryFieldRef()[patch][f];
                value=saltgeom::faceDiffusion(value,additional.boundaryField()[patch][f],additional.boundaryField()[patch][f]);
            }
        }
        return result;
    }
    // Associate an interface with a resolved liquid control volume. A cut
    // cell whose centre is still solid cannot use epsilonMin water as a
    // reaction reservoir. Its internal interface is conservatively assigned
    // through open liquid faces to adjacent liquid-centred cells instead.
    template<unsigned N>
    std::vector<std::vector<Foam::FixedList<Foam::scalar,N>>> exchangeRoutePackets(
        const std::vector<std::vector<Foam::FixedList<Foam::scalar,N>>>& outgoing)const {
        using namespace Foam;using Packet=FixedList<scalar,N>;
        labelList send(Pstream::nProcs(),0),receive(Pstream::nProcs(),0),neighbours;
        for(label rank=0;rank<Pstream::nProcs();rank++)send[rank]=outgoing[rank].size();
        UPstream::allToAll(send,receive);
        DynamicList<label> peers;for(label rank=0;rank<Pstream::nProcs();rank++)
            if(rank!=Pstream::myProcNo()&&(send[rank]||receive[rank]))peers.append(rank);
        neighbours=peers;
        std::vector<std::vector<Packet>> incoming(Pstream::nProcs());
        incoming[Pstream::myProcNo()]=outgoing[Pstream::myProcNo()];
        if(Pstream::parRun()){
            PstreamBuffers buffers(UPstream::commsTypes::nonBlocking);
            for(label rank:neighbours)if(send[rank]){
                List<Packet> values(send[rank]);forAll(values,i)values[i]=outgoing[rank][i];
                UOPstream stream(rank,buffers);stream<<values;
            }
            buffers.finishedNeighbourSends(neighbours);
            for(label rank:neighbours)if(receive[rank]){
                List<Packet> values;UIPstream stream(rank,buffers);stream>>values;
                if(values.size()!=receive[rank])FatalErrorInFunction<<"Surface route packet count mismatch"<<abort(FatalError);
                incoming[rank].assign(values.begin(),values.end());
            }
        }
        return incoming;
    }
    // Share every unresolved liquid edge's finite conductance among all
    // patches using it. A sparse hash owner sums only active interface-edge
    // records; no rank allocates a global cell or interface array.
    void updateRouteConductances(){
        using Packet=Foam::FixedList<Foam::scalar,11>;
        std::vector<std::vector<Packet>> outgoing(Foam::Pstream::nProcs());
        auto enqueue=[&](std::size_t link,std::size_t patch,ReactionPatch& p,Foam::label gx,Foam::label gy,Foam::label rx,Foam::label ry){
            p.waterResistance=0;double area=0;
            for(int j=0;j<7;j++)if(routeReactive[j])area+=p.area[j];
            if(!(area>0))return;
            for(int e=0;e<p.hops;e++){
                const auto& key=p.route[e].key;
                const auto hash=saltgeom::mix(std::uint64_t(key[0]))^saltgeom::mix(std::uint64_t(key[1])+0x7319)^saltgeom::mix(std::uint64_t(key[2])+0x82a1);
                const auto rank=hash%Foam::Pstream::nProcs();
                outgoing[rank].push_back(Packet({Foam::scalar(key[0]),Foam::scalar(key[1]),Foam::scalar(key[2]),Foam::scalar(link),Foam::scalar(patch),area,
                    Foam::scalar(gx),Foam::scalar(gy),Foam::scalar(rx),Foam::scalar(ry),Foam::scalar(e)}));
            }
        };
        for(std::size_t link=0;link<contacts.size();link++)for(std::size_t patch=0;patch<contacts[link].patches.size();patch++){
            const auto& geometry=surfaceLinks[link];enqueue(link,patch,contacts[link].patches[patch],geometry.gx,geometry.gy,geometry.rx,geometry.ry);
        }
        for(std::size_t cell=0;cell<localRoutes.size();cell++)for(std::size_t patch=0;patch<localRoutes[cell].size();patch++)
            enqueue(contacts.size()+cell,patch,localRoutes[cell][patch],cellIndex[cell][0],cellIndex[cell][1],cellIndex[cell][0],cellIndex[cell][1]);
        int any=0;for(const auto& rank:outgoing)any=any||!rank.empty();Foam::reduce(any,Foam::maxOp<int>());
        if(!any)return;
        const auto incoming=exchangeRoutePackets(outgoing);
        struct Entry {Packet value;int rank;};std::vector<Entry> entries;
        for(int rank=0;rank<Foam::Pstream::nProcs();rank++)for(const auto& value:incoming[rank])entries.push_back({value,rank});
        std::sort(entries.begin(),entries.end(),[](const Entry& a,const Entry& b){
            for(int field:{0,1,2,6,7,8,9,4}){if(a.value[field]<b.value[field])return true;if(a.value[field]>b.value[field])return false;}return false;
        });
        using Reply=Foam::FixedList<Foam::scalar,4>;
        std::vector<std::vector<Reply>> replies(Foam::Pstream::nProcs());
        for(std::size_t first=0;first<entries.size();){
            std::size_t end=first+1;const auto& key=entries[first].value;
            while(end<entries.size()&&entries[end].value[0]==key[0]&&entries[end].value[1]==key[1]&&entries[end].value[2]==key[2])end++;
            long double area=0;for(std::size_t n=first;n<end;n++)area+=entries[n].value[5];
            for(std::size_t n=first;n<end;n++){const auto& entry=entries[n];replies[entry.rank].push_back(Reply({entry.value[3],entry.value[4],entry.value[10],double(area)}));}
            first=end;
        }
        const auto returned=exchangeRoutePackets(replies);
        for(const auto& rank:returned)for(const auto& reply:rank){
            const std::size_t id=reply[0],p=reply[1];
            auto& patch=id<contacts.size()?contacts[id].patches[p]:localRoutes[id-contacts.size()][p];
            patch.waterResistance+=patch.route[std::size_t(reply[2])].resistance*reply[3];
        }
    }
    void collectFluidSupport(){
        using Packet=Foam::FixedList<Foam::scalar,21>;
        std::vector<Packet> amount(contacts.size()),received(mesh.nCells());
        for(auto& p:received)p=0.;
        for(std::size_t link=0;link<contacts.size();link++)for(int j=0;j<7;j++){
            amount[link][j]=contacts[link].support[j];amount[link][7+j]=contacts[link].normal[j].x;amount[link][14+j]=contacts[link].normal[j].y;
        }
        faceExchange->scatterAddTargets(amount,received);
        forAll(exchange,i){
            for(int j=0;j<7;j++){
                incoming[i][j]=received[i][j];support[i][j]=localArea[i][j]+incoming[i][j];
                auto normal=shape[i].normal[j]*directLocalArea[i][j]+saltgeom::Point{received[i][7+j],received[i][14+j]};
                for(const auto& patch:localRoutes[i])normal=normal+patch.normal[j];
                const double norm=saltgeom::length(normal);reactionNormal[i][j]=norm>0?normal*(1/norm):saltgeom::Point{};
            }
            saltgeom::updatePartition(state[i],support[i],reactionNormal[i]);
            measures[i]=saltgeom::measure(state[i],support[i],mesh.V()[i],blocking,productLayer);
        }
        // Physical area and the receiver's film conductance are separate.
        // Do not attenuate area and then apply the same film a second time.
        for(std::size_t link=0;link<contacts.size();link++)for(int j=0;j<7;j++)
            contacts[link].area[j]=contacts[link].support[j];
        updateRouteConductances();
    }
    void reconstructFluidContacts(const saltchem::RateParameters& kinetics){
        contacts.clear();surfaceLinks.clear();faceOpening=0;
        auto containsCentre=[&](const saltgeom::Polygon& p){
            if(p.n<3)return false;const saltgeom::Point centre{dx*.5,dy*.5};
            for(int k=0;k<p.n;k++)if(saltgeom::cross(p.vertex[(k+1)%p.n]-p.vertex[k],centre-p.vertex[k])< -1e-12*dx*dy)return false;
            return true;
        };
        struct Receiver {Foam::label rank=-1,cell=-1,x=0,y=0,hops=3;double resistance=0;std::array<RouteEdge,2> route{};};
        std::vector<Receiver> receiver(mesh.nCells());
        forAll(exchange,i){
            localArea[i]={};directLocalArea[i]={};localRoutes[i].clear();support[i]={};reactionNormal[i]={};incoming[i]={};
            if(containsCentre(shape[i].liquid))receiver[i]={Foam::Pstream::myProcNo(),i,cellIndex[i][0],cellIndex[i][1],0,0};
        }
        struct BoundaryContact {Foam::label source;int direction;saltchem::MineralVector area;};
        std::vector<BoundaryContact> boundaryContacts;
        std::vector<std::array<double,4>> openArea(mesh.nCells());
        for(int d=0;d<4;d++){
            std::vector<std::array<double,4>> lower,upper;
            for(int end=0;end<2;end++){
                forAll(exchange,i){const auto a=saltgeom::aperture(shape[i].liquid,d^1,dx,dy);exchange[i]=end?a.hi:a.lo;}
                exchange.correctBoundaryConditions();(end?upper:lower)=stencil.values(exchange);
            }
            forAll(exchange,i)if(connected[i][d]){
                const saltgeom::Interval other{lower[i][d],upper[i][d]};const double depth=mesh.V()[i]/(dx*dy);
                openArea[i][d]=saltgeom::overlap(saltgeom::aperture(shape[i].liquid,d,dx,dy),other)*depth;
                faceOpening[directionalFace[i][d]]=openArea[i][d]/(d<2?dy*depth:dx*depth);
                BoundaryContact contact{i,d,{}};bool present=false;
                for(int j=0;j<7;j++){
                    contact.area[j]=saltgeom::overlap(saltgeom::aperture(shape[i].solid[j],d,dx,dy),other)*depth;
                    present=present||contact.area[j]>0;
                }
                if(present)boundaryContacts.push_back(contact);
            }
        }
        forAll(mesh.boundary(),patch)if(!mesh.boundary()[patch].coupled()&&mesh.boundary()[patch].type()!="empty"){
            const auto& cells=mesh.boundary()[patch].faceCells();
            forAll(cells,f){const int d=stencil.patchDirections[patch][f];const auto aperture=saltgeom::aperture(shape[cells[f]].liquid,d,dx,dy);
                faceOpening[mesh.boundary()[patch].start()+f]=std::max(0.,aperture.hi-aperture.lo)/(d<2?dy:dx);}
        }
        const auto& coupled=meshExchange->coupledFaces();
        auto haloReceivers=[&](const std::vector<Receiver>& values){
            std::vector<Foam::FixedList<Foam::scalar,14>> halo(coupled.size());
            for(std::size_t k=0;k<coupled.size();k++){
                const auto& r=values[mesh.faceOwner()[coupled[k]]];halo[k][0]=r.rank;halo[k][1]=r.cell;halo[k][2]=r.x;halo[k][3]=r.y;halo[k][4]=r.hops;halo[k][5]=r.resistance;
                for(int e=0;e<2;e++){for(int axis=0;axis<3;axis++)halo[k][6+4*e+axis]=r.route[e].key[axis];halo[k][9+4*e]=r.route[e].resistance;}
            }
            meshExchange->swapCoupled(halo);return halo;
        };
        auto neighbourReceiver=[&](int i,int d,const std::vector<Receiver>& values,const std::vector<Foam::FixedList<Foam::scalar,14>>& halo){
            const auto f=directionalFace[i][d];
            if(f<mesh.nInternalFaces())return values[mesh.owner()[f]==i?mesh.neighbour()[f]:mesh.owner()[f]];
            const auto patch=mesh.boundaryMesh().whichPatch(f),k=coupledOffsets[patch]+f-mesh.boundary()[patch].start();const auto& r=halo[k];
            Receiver result{Foam::label(r[0]),Foam::label(r[1]),Foam::label(r[2]),Foam::label(r[3]),Foam::label(r[4]),r[5]};
            for(int e=0;e<2;e++){for(int axis=0;axis<3;axis++)result.route[e].key[axis]=Foam::label(r[6+4*e+axis]);result.route[e].resistance=r[9+4*e];}
            return result;
        };
        // A bounded cut-cell agglomeration follows actual liquid openings.
        // It never connects two pores through the intervening mineral. Limit
        // the surface-to-reservoir path to two Cartesian links so the boundary
        // approximation remains local and converges with grid refinement.
        for(int hop=1;hop<=2;hop++){
            const auto previous=receiver;const auto halo=haloReceivers(previous);
            forAll(exchange,i)if(previous[i].hops!=0&&shape[i].liquid.area()>0){
                Receiver best=previous[i];
                for(int d=0;d<4;d++)if(connected[i][d]&&openArea[i][d]>0){
                    auto candidate=neighbourReceiver(i,d,previous,halo);if(candidate.rank<0||candidate.hops>=2)continue;
                    const int axis=d/2;auto coordinate=cellIndex[i];coordinate[axis]+=(d&1);
                    if(meshExchange->periodicAxes()[axis]&&coordinate[axis]==gridExtent[axis])coordinate[axis]=0;
                    const double resistance=(axis==0?dx:dy)/openArea[i][d];
                    candidate.route[candidate.hops]={{{coordinate[0],coordinate[1],axis}},resistance};
                    candidate.hops++;candidate.resistance+=resistance;
                    if(best.rank<0||candidate.hops<best.hops||(candidate.hops==best.hops&&
                       (candidate.resistance<best.resistance||(candidate.resistance==best.resistance&&
                        (candidate.x<best.x||(candidate.x==best.x&&candidate.y<best.y))))))best=candidate;
                }
                receiver[i]=best;
            }
        }
        const auto halo=haloReceivers(receiver);
        std::array<bool,7> relevant{};
        for(int j=0;j<7;j++)relevant[j]=kinetics.geometry.weights[j]>0&&(kinetics.neutral[j]>0||kinetics.acid[j]>0||kinetics.co2[j]>0);
        routeReactive=relevant;
        for(const auto& law:laws)relevant[law.substrate]=relevant[law.substrate]||law.prefactor>0;
        std::vector<std::vector<std::size_t>> bySource(mesh.nCells());
        int inaccessible=0,reported=0;
        auto assign=[&](int source,const Receiver& target,int j,double area,saltgeom::Point normal,int maximumHops){
            if(area<=0)return;
            if(target.rank<0||target.hops>maximumHops){
                inaccessible|=relevant[j];
                if(relevant[j]&&reported++<4)Foam::Pout<<"SURFACE_ROUTE no_liquid_reservoir centre="<<mesh.C()[source]<<" mineral="<<j<<" area="<<area<<" hops="<<target.hops<<Foam::endl;
                return;
            }
            auto addPatch=[&](std::vector<ReactionPatch>& patches){
                std::size_t patch=patches.size();
                for(std::size_t k=0;k<patches.size();k++){
                    bool same=patches[k].hops==target.hops;
                    for(int e=0;e<target.hops&&same;e++)same=patches[k].route[e].key==target.route[e].key&&patches[k].route[e].resistance==target.route[e].resistance;
                    if(same){patch=k;break;}
                }
                if(patch==patches.size()){patches.push_back(ReactionPatch{});patches.back().hops=target.hops;patches.back().route=target.route;}
                patches[patch].area[j]+=area;patches[patch].normal[j]=patches[patch].normal[j]+normal*area;
            };
            if(target.rank==Foam::Pstream::myProcNo()&&target.cell==source){
                localArea[source][j]+=area;
                if(target.hops==0)directLocalArea[source][j]+=area;else addPatch(localRoutes[source]);
                return;
            }
            std::size_t link=contacts.size();
            for(auto k:bySource[source])if(surfaceLinks[k].receiverRank==target.rank&&surfaceLinks[k].receiverCell==target.cell){link=k;break;}
            if(link==contacts.size()){
                surfaceLinks.push_back({source,target.rank,target.cell,cellIndex[source][0],cellIndex[source][1],target.x,target.y});
                contacts.push_back(FaceContact{});bySource[source].push_back(link);
            }
            contacts[link].support[j]+=area;contacts[link].normal[j]=contacts[link].normal[j]+normal*area;
            addPatch(contacts[link].patches);
        };
        forAll(exchange,i)for(int j=0;j<7;j++)assign(i,receiver[i],j,shape[i].internalArea[j],shape[i].normal[j],2);
        for(const auto& contact:boundaryContacts){
            const int d=contact.direction;const auto target=neighbourReceiver(contact.source,d,receiver,halo);
            const saltgeom::Point normal=d<2?saltgeom::Point{d==1?1.:-1.,0}:saltgeom::Point{0,d==3?1.:-1.};
            for(int j=0;j<7;j++)assign(contact.source,target,j,contact.area[j],normal,1);
        }
        Foam::reduce(inaccessible,Foam::maxOp<int>());
        if(inaccessible)throw SurfaceTransportModelError("fluidSide surface has no liquid-centred reservoir within two open Cartesian links; resolve the liquid pore or refine the mesh");
        int changed=!faceExchange||!faceExchange->matches(surfaceLinks);
        Foam::reduce(changed,Foam::maxOp<int>());
        if(changed)faceExchange=std::make_unique<SurfaceLinkExchange>(mesh,surfaceLinks,dx,dy,meshExchange->periodicAxes());
        collectFluidSupport();
    }
    const saltchem::MineralVector& faceArea(Foam::label link,bool own)const {
        static const saltchem::MineralVector empty{};return own?contacts[link].area:empty;
    }
    std::size_t facePatchCount(Foam::label link)const {return contacts[link].patches.size();}
    const saltchem::MineralVector& facePatchArea(Foam::label link,std::size_t patch)const {return contacts[link].patches[patch].area;}
    saltgeom::Point facePatchNormal(Foam::label link,std::size_t patch,int mineral)const {
        const auto& p=contacts[link].patches[patch];return p.area[mineral]>0?p.normal[mineral]*(1/p.area[mineral]):saltgeom::Point{};
    }
    saltchem::SurfaceFilmParameters withRouteWater(const saltchem::SurfaceFilmParameters& receiverFilm,double resistance)const {
        auto film=receiverFilm;
        if(resistance==0)return film;
        const double velocity=molecularDiffusivity/resistance;
        auto series=[](double a,double b){if(a==0||b==0)return 0.;if(std::isinf(a))return b;if(std::isinf(b))return a;
            const double small=std::min(a,b),large=std::max(a,b);return small/(1+small/large);};
        for(int j=0;j<7;j++){film.bareVelocity[j]=series(film.bareVelocity[j],velocity);film.coatedVelocity[j]=series(film.coatedVelocity[j],velocity);}
        return film;
    }
    saltchem::SurfaceFilmParameters faceFilm(Foam::label link,const saltchem::SurfaceFilmParameters& receiverFilm,std::size_t patch=0)const {
        return withRouteWater(receiverFilm,contacts[link].patches[patch].waterResistance);
    }
    std::size_t localPatchCount(int cell)const {return 1+localRoutes[cell].size();}
    const saltchem::MineralVector& localPatchArea(int cell,std::size_t patch)const {return patch==0?directLocalArea[cell]:localRoutes[cell][patch-1].area;}
    saltgeom::Point localPatchNormal(int cell,std::size_t patch,int mineral)const {
        if(patch==0)return shape[cell].normal[mineral];const auto& p=localRoutes[cell][patch-1];
        return p.area[mineral]>0?p.normal[mineral]*(1/p.area[mineral]):saltgeom::Point{};
    }
    saltchem::SurfaceFilmParameters localPatchFilm(int cell,const saltchem::SurfaceFilmParameters& receiverFilm,std::size_t patch=0)const {
        return patch==0?receiverFilm:withRouteWater(receiverFilm,localRoutes[cell][patch-1].waterResistance);
    }
    saltgeom::Point faceNormal(Foam::label link,bool own,int mineral)const {
        if(!own)return {};const auto& c=contacts[link];const double area=c.support[mineral];
        return area>0?c.normal[mineral]*(1/area):saltgeom::Point{};
    }
    Foam::boolList activeReactionFaces(const saltchem::RateParameters& kinetics)const {
        Foam::boolList active(contacts.size(),false);
        for(std::size_t link=0;link<contacts.size();link++)for(int j=0;j<7;j++)if(kinetics.geometry.weights[j]>0&&
            (kinetics.neutral[j]>0||kinetics.acid[j]>0||kinetics.co2[j]>0)&&contacts[link].area[j]>0)active[link]=true;
        return active;
    }
    void retreatExternal(int i,const saltchem::MineralVector& distance){
        for(int j=0;j<7;j++)for(int k=0;k<7;k++)if(state[i].population[7*j+k].fraction>0)
            state[i].population[7*j+k].retreat+=distance[k];
    }
    // Retire exhausted substrates without changing reactor ownership inside
    // a chemical half interval. Sparse graph origins own their mineral stock;
    // support is re-summed on the receiving rank after deleting dead links.
    void retireExternalSupports(const std::vector<saltchem::Inventory>& cells,const saltchem::MineralVector& volume){
        forAll(exchange,i)retireDepletedSurfaces(i,cells[i],volume);
        for(std::size_t link=0;link<contacts.size();link++)for(int j=0;j<7;j++){
            const auto i=surfaceLinks[link].source;
            if(!saltgeom::hasResolvedSubstrate(cells[i],state[i],j,volume)){
                contacts[link].area[j]=contacts[link].support[j]=0;contacts[link].normal[j]={};
                for(auto& patch:contacts[link].patches){patch.area[j]=0;patch.normal[j]={};}
            }
        }
        collectFluidSupport();
    }
    void reconstruct(const Foam::PtrList<Foam::volScalarField>& solid,const Foam::volScalarField& inert,const saltchem::RateParameters& kinetics){
        if(!enabled)return;
        if(productLayer.enabled){
            // This inexpensive film closure uses the project's neutral TST
            // laws. Reject unsupported affinity/catalytic laws explicitly.
            for(int j=0;j<7;j++)if(kinetics.acid[j]!=0||kinetics.co2[j]!=0||kinetics.affinityExponent[j]!=1)
                throw std::invalid_argument("Porous product layers require neutral, linear-affinity TST kinetics");
            productLayer.kinetic=kinetics.neutral;productLayer.validate();
        }
        std::array<std::vector<std::array<double,4>>,7> neighbours;
        std::vector<saltchem::MineralVector> resolved(mesh.nCells());
        forAll(solid[0],i){const auto sub=saltgeom::subgrid(state[i]);for(int j=0;j<7;j++){
            const double difference=solid[j][i]-sub[j];
            const double roundoff=64*std::numeric_limits<double>::epsilon()*std::max(solid[j][i],sub[j]);
            resolved[i][j]=sub[j]>0&&difference<=roundoff?0:std::max(0.,difference);
        }}
        for(int j=0;j<7;j++){
            forAll(exchange,i)exchange[i]=resolved[i][j];exchange.correctBoundaryConditions();neighbours[j]=stencil.values(exchange);
        }
        forAll(exchange,i){
            std::array<std::array<double,4>,7> n;for(int j=0;j<7;j++)n[j]=neighbours[j][i];
            shape[i]=saltgeom::reconstruct(resolved[i],inert[i],n,state[i].inheritedNormal,dx,dy,mesh.V()[i]/(dx*dy),epsilonMin);
            bare[i]=shape[i].internalArea;incoming[i]={};
        }
        if(fluidSide()){reconstructFluidContacts(kinetics);return;}
        // Only actual mesh connections participate: no artificial substrate
        // on walls/inlets, and processor interfaces have the same apertures.
        for(int phase=-1;phase<7;phase++){
            double total=0;if(phase>=0){for(const auto& r:resolved)total+=r[phase];Foam::reduce(total,Foam::sumOp<double>());if(total==0)continue;}
            std::vector<std::array<double,4>> accessible;
            if(phase>=0){
                forAll(exchange,i)exchange[i]=measures[i].accessibility[phase];
                exchange.correctBoundaryConditions();accessible=stencil.values(exchange);
            }
            for(int d=0;d<4;d++){
                std::vector<std::array<double,4>> lower,upper;
                for(int end=0;end<2;end++){
                    forAll(exchange,i){const auto a=saltgeom::aperture(phase<0?shape[i].liquid:shape[i].solid[phase],d^1,dx,dy);exchange[i]=end?a.hi:a.lo;}
                    exchange.correctBoundaryConditions();(end?upper:lower)=stencil.values(exchange);
                }
                forAll(exchange,i)if(connected[i][d]){
                    const saltgeom::Interval other{lower[i][d],upper[i][d]};const double depth=mesh.V()[i]/(dx*dy);
                    if(phase<0){for(int j=0;j<7;j++)bare[i][j]+=saltgeom::overlap(saltgeom::aperture(shape[i].solid[j],d,dx,dy),other)*depth;}
                    else incoming[i][phase]+=saltgeom::overlap(saltgeom::aperture(shape[i].liquid,d,dx,dy),other)*depth*accessible[i][d];
                }
            }
            if(phase<0)forAll(exchange,i)measures[i]=saltgeom::measure(state[i],bare[i],mesh.V()[i],blocking,productLayer);
        }
        forAll(exchange,i)measures[i]=saltgeom::measure(state[i],bare[i],mesh.V()[i],blocking,productLayer);
        // A closed resolved solid still has epsilonMin numerical pore water.
        // Applying reactions there plus bulk harmonic diffusion introduces a
        // half-solid-cell resistance; it is not a liquid-side sharp-interface
        // law. Diagnose geometry and enabled kinetics, independent of SI or
        // a concentration threshold, without changing the reaction equations.
        Foam::vector solidOwnerSurface(Foam::vector::zero);
        forAll(exchange,i)if(shape[i].liquid.area()==0){
            double area=0;
            for(int j=0;j<7;j++)if(kinetics.geometry.weights[j]>0&&
                (kinetics.neutral[j]>0||kinetics.acid[j]>0||kinetics.co2[j]>0))
                area+=std::max(0.,bare[i][j]-shape[i].internalArea[j])*measures[i].accessibility[j];
            if(area>0){solidOwnerSurface.x()+=1;solidOwnerSurface.y()+=area;}
        }
        // One collective carries both count and area. Every rank reaches the
        // same policy decision, including ranks with no affected local cells.
        Foam::reduce(solidOwnerSurface,Foam::sumOp<Foam::vector>());
        sharpSurfaceCells_=std::max(sharpSurfaceCells_,Foam::label(solidOwnerSurface.x()));
        sharpSurfaceArea_=std::max(sharpSurfaceArea_,double(solidOwnerSurface.y()));
        if(solidOwnerSurface.x()>0){
            if(solidOwnerSurfacePolicy=="reject")
                throw SurfaceTransportModelError("solidOwnerSurfacePolicy reject: grid-face mineral reactions in fully solid-owner cells require a liquid-side interface transport model");
            if(!reportedSolidOwnerSurface){
                if(Foam::Pstream::master())Foam::Info
                    <<"SURFACE_TRANSPORT_MODEL volumeAveraged: cells="<<Foam::label(solidOwnerSurface.x())
                    <<" exposed_area_m2="<<solidOwnerSurface.y()<<Foam::nl
                    <<"Reactions use numerical pore water in solid-owner cells. Harmonic bulk diffusion does not provide a liquid-side sharp-interface reaction flux; rates can depend on epsilonMin. Set solidOwnerSurfacePolicy reject to reject this geometry."
                    <<Foam::endl;
                reportedSolidOwnerSurface=true;
            }
        }
    }
    void retireDepletedSurfaces(int i,const saltchem::Inventory& cell,const saltchem::MineralVector& volume){
        for(int j=0;j<7;j++)if(!saltgeom::hasResolvedSubstrate(cell,state[i],j,volume)){
            // A local chemical half-step can remove the complete resolved
            // substrate before the next global geometric reconstruction.
            // Neighbour-owned growth faces stay available; the vanished
            // local solid cannot supply a nucleation/dissolution surface.
            bare[i][j]=0;shape[i].internalArea[j]=0;
            if(fluidSide()){
                support[i][j]=std::max(0.,support[i][j]-localArea[i][j]);localArea[i][j]=directLocalArea[i][j]=0;
                for(auto& patch:localRoutes[i]){patch.area[j]=0;patch.normal[j]={};}
            }
        }
    }
    saltgeom::SurfaceMeasure localMeasure(int i)const {
        auto value=measures[i];
        if(fluidSide())for(int j=0;j<7;j++)value.exposed[j]=localArea[i][j]*value.accessibility[j]+value.capArea[j];
        return value;
    }
    void configure(int i,const saltchem::Inventory& cell,saltchem::RateParameters& r,const saltchem::MineralVector& volume){
        if(!enabled)return;
        retireDepletedSurfaces(i,cell,volume);
        measures[i]=saltgeom::measure(state[i],supportArea(i),mesh.V()[i],blocking,productLayer);
        const auto local=localMeasure(i);
        saltchem::MineralVector fractions{};for(int j=0;j<7;j++)fractions[j]=cell.mineral[j]*volume[j];
        const auto locked=saltgeom::inaccessibleFractions(state[i],fractions,local);
        r.geometry.separateMineralAreas=true;
        for(int j=0;j<7;j++){
            r.area[j]=local.exposed[j]/mesh.V()[i];r.geometry.growthArea[j]=growthArea(i,j)/mesh.V()[i];
            r.geometry.lockedMoles[j]=std::min(cell.mineral[j],locked[j]/volume[j]);
            r.referenceMoles[j]=std::max(cell.mineral[j]-r.geometry.lockedMoles[j],1e-30);
            // Cap area grows with V^(2/3); resolved interfaces stay frozen
            // during this first-order geometry/chemistry split step.
            r.geometry.surfaceExponent[j]=local.capArea[j]>0&&local.exposed[j]==local.capArea[j]&&incoming[i][j]==0?2./3:0;
        }
    }
    double growthArea(int i,int j)const {
        // A grid-face interface dissolves in its solid owner and grows in its
        // liquid owner. Counting both owners for precipitation doubles area.
        if(fluidSide())return (localArea[i][j]+incoming[i][j])*measures[i].accessibility[j]+measures[i].capArea[j];
        return shape[i].internalArea[j]*measures[i].accessibility[j]+measures[i].capArea[j]+incoming[i][j];
    }
    void seedCell(int i,saltchem::Inventory& cell,const saltchem::MineralVector& si,double dt,const saltchem::MineralVector& volume,double inert,double temperatureK){
        if(!enabled)return;
        if(fluidSide()){
            saltchem::MineralVector available{};
            std::array<double,49> prior{};
            for(int pair=0;pair<49;pair++)prior[pair]=state[i].population[pair].fraction;
            for(int j=0;j<7;j++)available[j]=(incoming[i][j]>0||saltgeom::hasResolvedSubstrate(cell,state[i],j,volume))?1:0;
            saltgeom::nucleate(cell,state[i],laws,si,support[i],dt,mesh.V()[i],volume,epsilonMin,inert,seed,keys[i],productLayer,&available,temperatureK);
            for(int j=0;j<7;j++)for(int k=0;k<7;k++){
                const double added=state[i].population[7*j+k].fraction-prior[7*j+k];
                reactionVectors[i][0]-=added*reactionNormal[i][k].x;
                reactionVectors[i][1]-=added*reactionNormal[i][k].y;
                if(state[i].population[7*j+k].fraction>0&&saltgeom::length(reactionNormal[i][k])>0)
                    state[i].inheritedNormal[j]=reactionNormal[i][k];
            }
        }else saltgeom::nucleate(cell,state[i],laws,si,bare[i],dt,mesh.V()[i],volume,epsilonMin,inert,seed,keys[i],productLayer,nullptr,temperatureK);
    }
    void commit(int i,const saltchem::Inventory& start,const saltchem::Inventory& end,const saltchem::MineralVector& volume){
        if(!enabled)return;
        saltchem::MineralVector a{},b{},growth{};
        for(int j=0;j<7;j++){a[j]=start.mineral[j]*volume[j];b[j]=end.mineral[j]*volume[j];growth[j]=growthArea(i,j);}
        if(fluidSide())for(int j=0;j<7;j++){
            // Capture the direction while the reacting interface still
            // exists; a depleted phase has no normal after reconstruction.
            auto normal=shape[i].normal[j];
            if(saltgeom::length(normal)==0)normal=state[i].inheritedNormal[j];
            if(saltgeom::length(normal)==0)normal=reactionNormal[i][j];
            reactionVectors[i][0]+=(a[j]-b[j])*normal.x;
            reactionVectors[i][1]+=(a[j]-b[j])*normal.y;
        }
        const auto oldCaps=saltgeom::subgrid(state[i]);
        saltgeom::grow(state[i],a,b,localMeasure(i),growth,mesh.V()[i],resolvedRadius,
            fluidSide()?reactionNormal[i]:shape[i].normal,supportArea(i),productLayer,epsilonMin,nullptr,!fluidSide());
        if(fluidSide()){
            const auto newCaps=saltgeom::subgrid(state[i]);saltchem::MineralVector retreat{};
            for(int k=0;k<7;k++)if(support[i][k]>0)
                retreat[k]=std::max(0.,(a[k]-oldCaps[k])-(b[k]-newCaps[k]))*mesh.V()[i]/(support[i][k]*(1-epsilonMin));
            retreatExternal(i,retreat);
        }
        retireDepletedSurfaces(i,end,volume);
    }
    bool active(int i)const {const auto local=localMeasure(i);for(int j=0;j<7;j++)if(local.exposed[j]>0||incoming[i][j]>0)return true;return false;}
    double timeError(const std::vector<saltgeom::SurfaceState>& coarse)const {
        if(!enabled)return 0;double error=0;
        for(std::size_t i=0;i<state.size();i++){
            const auto a=saltgeom::measure(coarse.at(i),supportArea(i),mesh.V()[i],blocking,productLayer);
            const auto b=saltgeom::measure(state[i],supportArea(i),mesh.V()[i],blocking,productLayer);
            for(int j=0;j<7;j++){
                error=std::max(error,std::abs(a.coverage[j]-b.coverage[j])/surfaceTimeTolerance);
                if(fluidSide()){
                    error=std::max(error,std::abs(a.attachedFraction[j]-b.attachedFraction[j])/surfaceTimeTolerance);
                    error=std::max(error,std::abs(a.layerThickness[j]-b.layerThickness[j])/(std::min(dx,dy)*surfaceTimeTolerance));
                }else error=std::max(error,std::abs(a.accessibility[j]-b.accessibility[j])/surfaceTimeTolerance);
            }
        }
        Foam::reduce(error,Foam::maxOp<double>());return error;
    }
    void beginReactions(){
        if(!fluidSide())return;
        std::fill(reactionRates.begin(),reactionRates.end(),saltchem::MineralVector{});
        std::fill(reactionVectors.begin(),reactionVectors.end(),std::array<double,2>{});
    }
    void accumulateFaceVector(int i,const std::array<double,2>& integrated){
        if(!fluidSide())return;
        for(int axis=0;axis<2;axis++)reactionVectors[i][axis]+=integrated[axis];
    }
    void recordReactions(const std::vector<saltchem::Inventory>& before,
                         const std::vector<saltchem::Inventory>& after,double dt){
        if(!fluidSide())return;
        if(before.size()!=state.size()||after.size()!=state.size()||!std::isfinite(dt)||dt<=0)
            throw std::invalid_argument("Invalid surface reaction diagnostic state or interval");
        for(std::size_t i=0;i<state.size();i++)for(int j=0;j<7;j++){
            const double rate=(before[i].mineral[j]-after[i].mineral[j])/dt;
            if(!std::isfinite(rate))throw std::runtime_error("Nonfinite surface reaction diagnostic rate");
            reactionRates[i][j]=rate;
        }
        for(auto& vector:reactionVectors)for(double& value:vector){
            value/=dt;
            if(!std::isfinite(value))throw std::runtime_error("Nonfinite surface reaction diagnostic vector");
        }
    }
    std::array<double,2> transformation(int i,const saltchem::MineralVector& rates,const saltchem::MineralVector& volumes)const {
        if(fluidSide())return reactionVectors[i];
        std::array<double,2> q{};
        for(int j=0;j<7;j++){
            auto normal=shape[i].normal[j];
            if(saltgeom::length(normal)==0)for(int k=0;k<7;k++)if(state[i].population[7*j+k].fraction>0){normal=shape[i].normal[k];break;}
            const double rate=rates[j];
            q[0]+=rate*volumes[j]*normal.x;q[1]+=rate*volumes[j]*normal.y;
        }return q;
    }
    void checkpoint(bool write){
        if(!enabled)return;
        const auto path=mesh.time().timePath()/"mineralSurface.bin";
        const std::uint64_t recordBytes=sizeof(std::uint64_t)+sizeof(saltgeom::SurfaceState)+(fluidSide()?sizeof(saltchem::MineralVector)+sizeof(std::array<double,2>):0);
        if(write){
            std::ofstream f(path.c_str(),std::ios::binary);std::uint64_t h[6]={0x53414c5453555246ULL,6,state.size(),sizeof(saltgeom::SurfaceState),fluidSide()?1ULL:0ULL,recordBytes};
            f.write(reinterpret_cast<const char*>(h),sizeof(h));
            for(std::size_t i=0;i<state.size();i++){
                f.write(reinterpret_cast<const char*>(&keys[i]),sizeof(keys[i]));f.write(reinterpret_cast<const char*>(&state[i]),sizeof(state[i]));
                if(fluidSide()){
                    f.write(reinterpret_cast<const char*>(&reactionRates[i]),sizeof(reactionRates[i]));
                    f.write(reinterpret_cast<const char*>(&reactionVectors[i]),sizeof(reactionVectors[i]));
                }
            }
            if(!f)throw std::runtime_error("Failed to write mineral surface checkpoint");
        }else{
            std::ifstream f(path.c_str(),std::ios::binary);std::uint64_t h[6]{};f.read(reinterpret_cast<char*>(h),sizeof(h));
            if(!f||h[0]!=0x53414c5453555246ULL||h[1]!=6||h[2]!=state.size()||h[3]!=sizeof(saltgeom::SurfaceState)||h[4]!=(fluidSide()?1ULL:0ULL)||h[5]!=recordBytes)
                throw std::runtime_error("Missing/incompatible mineralSurface.bin restart; version 6 requires saved curvature/film-partition history, reaction rates/vectors, the same surface transport model and mesh/partition");
            f.seekg(0,std::ios::end);
            if(f.tellg()!=std::streamoff(sizeof(h)+state.size()*recordBytes))throw std::runtime_error("Truncated/oversized mineralSurface.bin restart");
            f.seekg(sizeof(h));
            for(std::size_t i=0;i<state.size();i++){
                std::uint64_t key;f.read(reinterpret_cast<char*>(&key),sizeof(key));f.read(reinterpret_cast<char*>(&state[i]),sizeof(state[i]));
                if(fluidSide()){
                    f.read(reinterpret_cast<char*>(&reactionRates[i]),sizeof(reactionRates[i]));
                    f.read(reinterpret_cast<char*>(&reactionVectors[i]),sizeof(reactionVectors[i]));
                    for(double rate:reactionRates[i])if(!std::isfinite(rate))throw std::runtime_error("Nonfinite mineral surface reaction checkpoint");
                    for(double value:reactionVectors[i])if(!std::isfinite(value))throw std::runtime_error("Nonfinite mineral surface reaction vector checkpoint");
                }
                if(!f||key!=keys[i])throw std::runtime_error("Mineral surface checkpoint mesh mismatch");
            }
        }
    }
    void save(const Foam::fileName& output)const {
        if(!enabled)return;
        Foam::mkDir(output/"interfaces");std::ofstream f((output/"interfaces"/(mesh.time().timeName()+".csv")).c_str());
        f<<std::setprecision(16)<<"x_m,y_m";
        for(const char* prefix:{"bare_m2_","exposed_m2_","coverage_","nuclei_","subgrid_fraction_","births_","promoted_","growth_m2_","layer_coverage_","layer_thickness_m_","layer_transmission_","accessibility_","attached_nuclei_","attached_fraction_","residue_fraction_"})for(const char* name:saltchem::mineralNames)f<<','<<prefix<<name;
        f<<",transport_fraction,transport_resistance_x,transport_resistance_y\n";
        forAll(exchange,i){f<<mesh.C()[i].x()<<','<<mesh.C()[i].y();const auto sub=saltgeom::subgrid(state[i]);saltchem::MineralVector counts{};for(int j=0;j<7;j++)for(int k=0;k<7;k++)counts[j]+=state[i].population[7*j+k].count;
            saltchem::MineralVector growth{};for(int j=0;j<7;j++)growth[j]=growthArea(i,j);
            for(const auto& values:{supportArea(i),measures[i].exposed,measures[i].coverage,counts,sub,state[i].births,state[i].promotions,growth,
                                  measures[i].layerCoverage,measures[i].layerThickness,measures[i].layerTransmission,measures[i].accessibility,
                                  measures[i].attachedCount,measures[i].attachedFraction,measures[i].residueFraction})for(double v:values)f<<','<<v;
            f<<','<<transportFraction[i]<<','<<transportResistance[i][0]<<','<<transportResistance[i][1]<<'\n';
        }if(!f)throw std::runtime_error("Failed to write mineral interface diagnostics");
    }
};
}
