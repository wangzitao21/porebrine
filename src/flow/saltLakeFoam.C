// SPDX-License-Identifier: GPL-3.0-or-later
// Seven-component conservative, variable-density DBS reactive transport.
// Native and reference executables share every spatial and temporal equation.
#include "BrineReactions.hpp"
#include "BrineInlet.hpp"
#include "SurfaceReaction.hpp"
#ifdef SALT_REFERENCE_BACKEND
#include "PhreeqcBackend.hpp"
#endif
#include "fvCFD.H"
#include "simpleControl.H"
#include "PhaseTransform.H"
#include "MineralSurface.H"
#include "ThermodynamicSchedule.H"
#include "ComponentTransport.H"
#include "LocalChemistry.H"
#include "HydraulicBoundaries.H"
#include "CouplingAcceleration.H"
#include "TimeStepControl.H"
#include "constrainPressure.H"
#include "labelIOList.H"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <memory>
#include <numeric>
#include <omp.h>
using namespace Foam;
using Timer=std::chrono::steady_clock;
double elapsed(Timer::time_point t){return std::chrono::duration<double>(Timer::now()-t).count();}
// Euler needs the accepted cell state and its conservative face flux. The
// controller and ledger are also continued, so a restart cannot reset a
// conservation error or change the next accepted step/output event.
struct FlowCheckpoint {
    enum Control {time,offset,nextFrame,maxBalance,nominalDt,proposedDt,writeOrigin,
        minAcceptedDt,maxAcceptedDt,deltaT,maximumDt,acceptedSteps,rejectedSteps,
        trialSegments,lastTimeError,maxTimeError,count};
    std::array<double,count> control{};
    std::array<double,8> initial{},boundary{},balance{},inventory{};
    saltchem::MineralVector dissolved{},precipitated{};
    scalarField deviation;
    std::vector<double> fields;
    std::uint64_t meshHash=0,index=0;
    static std::uint64_t geometryHash(const fvMesh& mesh){
        std::uint64_t hash=1469598103934665603ULL;
        auto add=[&](const void* data,std::size_t size){const auto* bytes=static_cast<const unsigned char*>(data);for(std::size_t i=0;i<size;i++){hash^=bytes[i];hash*=1099511628211ULL;}};
        for(const auto& point:mesh.points())for(int j=0;j<3;j++){const double value=point[j];add(&value,sizeof(value));}
        for(const auto& face:mesh.faces()){const label size=face.size();add(&size,sizeof(size));for(label vertex:face)add(&vertex,sizeof(vertex));}
        for(label cell:mesh.faceOwner())add(&cell,sizeof(cell));
        for(label cell:mesh.faceNeighbour())add(&cell,sizeof(cell));
        for(const auto& patch:mesh.boundary()){const std::string name(patch.name()),type(patch.type());add(name.data(),name.size());add(type.data(),type.size());const label size=patch.size();add(&size,sizeof(size));}
        return hash;
    }
    void read(const fileName& path,const fvMesh& mesh){
        std::ifstream stream(path.c_str(),std::ios::binary);std::array<std::uint64_t,9> header{};
        stream.read(reinterpret_cast<char*>(header.data()),sizeof(header));
        if(!stream||header[0]!=0x53414c54464c4f57ULL||header[1]!=2||header[2]!=std::uint64_t(mesh.nCells())||header[3]!=std::uint64_t(mesh.nFaces())||header[4]!=geometryHash(mesh)||header[6]!=std::uint64_t(Pstream::nProcs())||header[7]!=std::uint64_t(Pstream::myProcNo())||header[8]>23*std::uint64_t(mesh.nCells()+mesh.nFaces()))
            throw std::runtime_error("Missing/incompatible flowState.bin restart; a complete checkpoint on the same mesh and partition is required");
        meshHash=header[4];index=header[5];deviation.setSize(mesh.nCells());fields.resize(header[8]);
        auto get=[&](auto& values){stream.read(reinterpret_cast<char*>(values.data()),sizeof(values[0])*values.size());};
        get(control);get(initial);get(boundary);get(balance);get(inventory);get(dissolved);get(precipitated);get(deviation);get(fields);
        bool valid=bool(stream)&&stream.peek()==std::char_traits<char>::eof();
        auto finite=[&](const auto& values){for(double x:values)valid=valid&&std::isfinite(x);};
        finite(control);finite(initial);finite(boundary);finite(balance);finite(inventory);finite(dissolved);finite(precipitated);finite(deviation);finite(fields);
        if(!valid||control[nextFrame]<=control[time]||control[nominalDt]<=0||control[proposedDt]<=0||control[deltaT]<=0||control[maximumDt]<=0)
            throw std::runtime_error("Corrupt flowState.bin restart");
    }
    void write(const fileName& path,const fvMesh& mesh)const{
        const fileName temporary(path+".tmp");std::ofstream stream(temporary.c_str(),std::ios::binary|std::ios::trunc);
        const std::array<std::uint64_t,9> header{{0x53414c54464c4f57ULL,2,std::uint64_t(mesh.nCells()),std::uint64_t(mesh.nFaces()),geometryHash(mesh),index,std::uint64_t(Pstream::nProcs()),std::uint64_t(Pstream::myProcNo()),fields.size()}};
        stream.write(reinterpret_cast<const char*>(header.data()),sizeof(header));
        auto put=[&](const auto& values){if(!values.empty())stream.write(reinterpret_cast<const char*>(&values[0]),sizeof(values[0])*values.size());};
        put(control);put(initial);put(boundary);put(balance);put(inventory);put(dissolved);put(precipitated);put(deviation);put(fields);stream.close();
        if(!stream)throw std::runtime_error("Failed to write flowState.bin checkpoint");
        if(!Foam::mv(temporary,path))throw std::runtime_error("Failed to finalize flowState.bin checkpoint");
    }
    // Native field files preserve patch types, while this binary payload
    // preserves numerical state independently of writeFormat/writePrecision.
    template<class Field> void transfer(Field& field,bool loading,std::size_t& position){
        auto values=[&](auto& list){for(auto& value:list)for(direction c=0;c<pTraits<typename Field::value_type>::nComponents;c++){
            if(loading){if(position>=fields.size())throw std::runtime_error("Incomplete checkpoint field payload");setComponent(value,c)=fields[position++];}
            else fields.push_back(component(value,c));
        }};
        values(field.primitiveFieldRef());for(auto& patch:field.boundaryFieldRef())values(patch);
    }
    void transferFields(bool loading,volScalarField& p,volScalarField& eps,volScalarField& rho,volScalarField& inert,
                        volVectorField& U,surfaceScalarField& phi,surfaceScalarField& massPhi,PtrList<volScalarField>& Y,PtrList<volScalarField>& solid){
        std::size_t position=0;if(!loading)fields.clear();
        for(auto* field:{&p,&eps,&rho,&inert})transfer(*field,loading,position);
        transfer(U,loading,position);transfer(phi,loading,position);transfer(massPhi,loading,position);
        for(int j=0;j<7;j++){transfer(Y[j],loading,position);transfer(solid[j],loading,position);}
        if(loading&&position!=fields.size())throw std::runtime_error("Checkpoint field layout mismatch");
    }
};
int main(int argc,char** argv){
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    const bool restarting=runTime.value()!=0;
    const double restartTime=runTime.value();FlowCheckpoint restart;
    if(restarting){
        restart.read(runTime.timePath()/"flowState.bin",mesh);
        const double tolerance=1e-8*max(std::abs(restartTime),1.);
        if(std::abs(restart.control[FlowCheckpoint::time]-restartTime)>tolerance)throw std::runtime_error("Checkpoint time does not match the selected restart directory");
        runTime.setTime(restart.control[FlowCheckpoint::time],label(restart.index));runTime.setDeltaT(restart.control[FlowCheckpoint::deltaT],false);
    }
    const auto wallStart=Timer::now();simpleControl simple(mesh);
    IOdictionary cfg(IOobject("reactiveProperties",runTime.constant(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE));
    const ThermodynamicSchedule thermalSchedule(cfg);const auto initialTP=thermalSchedule.at(runTime.value());
    fileName db(cfg.lookup("database"));db.expand();if(!db.isAbsolute())db=runTime.globalPath()/db;
    saltchem::PitzerModel model(db.c_str(),initialTP[0]+273.15,initialTP[1]);saltchem::Speciation spec(model);saltchem::Reactions chemistry(model);
    const word waterPolicy=cfg.lookupOrDefault<word>("waterPropertyPolicy","phreeqcCompatible");
    if(waterPolicy!="strict"&&waterPolicy!="phreeqcCompatible")throw std::invalid_argument("Unknown waterPropertyPolicy");model.setStrictWater(waterPolicy=="strict");
    const double dx=readScalar(cfg.lookup("dx")),dy=readScalar(cfg.lookup("dy")),epMin=cfg.lookupOrDefault<scalar>("epsilonMin",1e-5);
    const double flowTol=cfg.lookupOrDefault<scalar>("flowTolerance",1e-5),couplingTol=cfg.lookupOrDefault<scalar>("couplingTolerance",1e-8),frame=cfg.lookupOrDefault<scalar>("frameInterval",1);
    const double conservationTolerance=cfg.lookupOrDefault<scalar>("conservationTolerance",1e-7);
    if(!(conservationTolerance>0&&conservationTolerance<=1e-7))throw std::invalid_argument("conservationTolerance must be in (0, 1e-7]");
    const double anchorThreshold=cfg.lookupOrDefault<scalar>("anchorThreshold",.01),anchorStrength=cfg.lookupOrDefault<scalar>("anchorStrength",1e12);
    const double anchorWidth=cfg.lookupOrDefault<scalar>("anchorTransitionWidth",0);
    if(!std::isfinite(anchorThreshold)||!std::isfinite(anchorStrength)||!std::isfinite(anchorWidth)||anchorThreshold<0||anchorStrength<0||anchorWidth<0||anchorWidth>2*anchorThreshold)throw std::invalid_argument("Invalid solid anchoring controls");
    const int maxIterations=cfg.lookupOrDefault<label>("maxFlowIterations",1500),threads=cfg.lookupOrDefault<label>("chemistryThreads",2);
    const int maxCouplingIterations=cfg.lookupOrDefault<label>("maxCouplingIterations",min(200,maxIterations));
    if(maxCouplingIterations<6||maxCouplingIterations>maxIterations)throw std::invalid_argument("maxCouplingIterations must be between 6 and maxFlowIterations");
    const int flowCorrectors=cfg.lookupOrDefault<label>("flowCorrectors",5);
    const bool adaptiveFlowCorrectors=cfg.lookupOrDefault<bool>("adaptiveFlowCorrectors",false);
    const bool conservativeFlowClosure=cfg.lookupOrDefault<bool>("conservativeFlowClosure",true);
    const bool accelerateFaceFlux=cfg.lookupOrDefault<bool>("accelerateFaceFlux",false);
    const bool consistentMomentum=cfg.lookupOrDefault<bool>("consistentMomentum",true);
    const word faceMomentumInterpolation=cfg.lookupOrDefault<word>("faceMomentumInterpolation","linear");
    if(faceMomentumInterpolation!="linear"&&faceMomentumInterpolation!="harmonic")throw std::invalid_argument("faceMomentumInterpolation must be linear or harmonic");
    const bool detectCouplingStagnation=cfg.lookupOrDefault<bool>("detectCouplingStagnation",true);
    const bool traceCoupling=cfg.lookupOrDefault<bool>("traceCoupling",false);
    if(flowCorrectors<1||flowCorrectors>50)throw std::invalid_argument("flowCorrectors must be between 1 and 50");
    if(!std::isfinite(flowTol)||!std::isfinite(couplingTol)||!std::isfinite(frame)||flowTol<=0||couplingTol<=0||frame<=0||maxIterations<1||threads<1)throw std::invalid_argument("Invalid flow, coupling, output or thread controls");
    if(!saltflow::schemeIs(mesh.ddtScheme("ddt(rho,U)"),{"Euler"}))throw std::invalid_argument("Momentum/face-flux time integration requires Euler discretization");
    omp_set_num_threads(threads);
    const bool reuseTransport=cfg.lookupOrDefault<bool>("reuseTransportMatrix",true),accelerateCoupling=cfg.lookupOrDefault<bool>("accelerateCoupling",true);
    const word couplingAcceleration=cfg.lookupOrDefault<word>("couplingAcceleration",accelerateCoupling?"anderson":"none");
    if(couplingAcceleration!="anderson"&&couplingAcceleration!="aitken"&&couplingAcceleration!="none")throw std::invalid_argument("Unknown couplingAcceleration");
    const double minimumCouplingRelaxation=cfg.lookupOrDefault<scalar>("minimumCouplingRelaxation",.001);
    const double maximumCouplingRelaxation=cfg.lookupOrDefault<scalar>("maximumCouplingRelaxation",1.);
    if(!(minimumCouplingRelaxation>0&&minimumCouplingRelaxation<=maximumCouplingRelaxation&&maximumCouplingRelaxation<=1))throw std::invalid_argument("Coupling relaxation bounds must satisfy 0 < minimum <= maximum <= 1");
    const bool flowPressure=cfg.lookupOrDefault<word>("chemistryPressureMode","flow")=="flow";
    const word densityPressureMode=cfg.lookupOrDefault<word>("densityPressureMode","local");
    if(densityPressureMode!="local"&&densityPressureMode!="thermodynamic")throw std::invalid_argument("densityPressureMode must be local or thermodynamic");
    const bool pressureCompressible=flowPressure&&densityPressureMode=="local";
    const double densityPressureTolerance=cfg.lookupOrDefault<scalar>("densityPressureTolerance",1e-4);
    if(!(densityPressureTolerance>0&&densityPressureTolerance<=.01))throw std::invalid_argument("Invalid densityPressureTolerance");
    const bool linearizeEOS=cfg.lookupOrDefault<bool>("linearizePressureEOS",true);
    if(!flowPressure&&cfg.lookupOrDefault<word>("chemistryPressureMode","flow")!="uniform")throw std::invalid_argument("chemistryPressureMode must be flow or uniform");
    saltflow::LocalChemistry localChemistry(model,threads);
    const double pressureRelaxation=cfg.lookupOrDefault<scalar>("pressureRelaxation",1),positivityTolerance=cfg.lookupOrDefault<scalar>("transportRoundoffTolerance",1e-12);
    if(!(pressureRelaxation>0&&pressureRelaxation<=1&&positivityTolerance>0&&positivityTolerance<=1e-10))throw std::invalid_argument("Invalid pressure relaxation or round-off bound");
    saltchem::IntegrationOptions options;options.relativeTolerance=cfg.lookupOrDefault<scalar>("kineticTolerance",options.relativeTolerance);options.absoluteTolerance=cfg.lookupOrDefault<scalar>("kineticAbsoluteTolerance",options.absoluteTolerance);options.maximumAttempts=cfg.lookupOrDefault<label>("kineticMaximumAttempts",options.maximumAttempts);options.method=cfg.lookupOrDefault<label>("kineticIntegrator",options.method);
    options.analyticJacobian=cfg.lookupOrDefault<bool>("analyticChemicalJacobian",true);
    List<scalar> coefficients(cfg.lookup("rateCoefficients")),inletFractions(cfg.lookup("inletFractions"));if(coefficients.size()!=7||inletFractions.size()!=7)FatalErrorInFunction<<"Seven rates and inlet fractions required"<<exit(FatalError);
    saltchem::RateParameters rateTemplate;
    for(int j=0;j<7;j++)rateTemplate.neutral[j]=coefficients[j];
    auto optionalRates=[&](const word& key,saltchem::MineralVector& values){
        if(!cfg.found(key))return;const List<scalar> input(cfg.lookup(key));
        if(input.size()!=7)FatalErrorInFunction<<key<<" requires seven mineral values"<<exit(FatalError);
        for(int j=0;j<7;j++)values[j]=input[j];
    };
    optionalRates("acidCoefficients",rateTemplate.acid);
    optionalRates("co2Coefficients",rateTemplate.co2);
    optionalRates("acidOrders",rateTemplate.acidOrder);
    optionalRates("dissolutionExponents",rateTemplate.dissolutionExponent);
    optionalRates("affinityExponents",rateTemplate.affinityExponent);
    optionalRates("surfaceWeights",rateTemplate.geometry.weights);
    saltchem::MineralVector neutralEnergy{},acidEnergy{},co2Energy{};
    optionalRates("rateActivationEnergiesKJ",neutralEnergy);
    optionalRates("acidActivationEnergiesKJ",acidEnergy);
    optionalRates("co2ActivationEnergiesKJ",co2Energy);
    const double rateReferenceK=cfg.lookupOrDefault<scalar>("rateReferenceTemperatureC",25)+273.15;
    if(!(rateReferenceK>0)||!std::isfinite(rateReferenceK))throw std::invalid_argument("Invalid rate reference temperature");
    const auto neutralReference=rateTemplate.neutral,acidReference=rateTemplate.acid,co2Reference=rateTemplate.co2;
    auto updateRateTemperature=[&](){
        const double inverse=(1/model.temperatureK()-1/rateReferenceK)/.00831470;
        for(int j=0;j<7;j++){
            rateTemplate.neutral[j]=neutralReference[j]*saltchem::checkedExp(-neutralEnergy[j]*inverse);
            rateTemplate.acid[j]=acidReference[j]*saltchem::checkedExp(-acidEnergy[j]*inverse);
            rateTemplate.co2[j]=co2Reference[j]*saltchem::checkedExp(-co2Energy[j]*inverse);
        }
    };updateRateTemperature();
    rateTemplate.geometry.enabled=true;rateTemplate.geometry.bulkVolume=1;
    rateTemplate.geometry.dx=dx;rateTemplate.geometry.dy=dy;rateTemplate.geometry.epsilonMin=epMin;
    rateTemplate.geometry.directionWidth=cfg.lookupOrDefault<scalar>("phaseDirectionWidth",rateTemplate.geometry.directionWidth);
    rateTemplate.geometry.closureWidth=cfg.lookupOrDefault<scalar>("poreClosureWidth",rateTemplate.geometry.closureWidth);
    saltchem::Reactions::validateParameters(rateTemplate);
    saltchem::Inventory inletState;inletState.waterEquivalent=1;for(int j=0;j<7;j++){inletState.aqueous[j]=inletFractions[j]/saltchem::componentMass[j];inletState.waterEquivalent-=inletFractions[j];}
    const word inletChemistry=cfg.lookupOrDefault<word>("inletChemistry","fixedFractions");
    if(inletChemistry!="fixedFractions"&&inletChemistry!="haliteSaturated")throw std::invalid_argument("Unknown inletChemistry");
    auto updateInlet=[&](){
        if(inletChemistry=="haliteSaturated"){
            inletState=saltchem::haliteSaturatedInlet(spec);const double mass=saltchem::aqueousMass(inletState);
            for(int j=0;j<7;j++)inletFractions[j]=saltchem::componentMass[j]*inletState.aqueous[j]/mass;
        }
    };updateInlet();
    double inletWater;const auto inletProperties=spec.conserved(inletState,inletWater);
    dimensionedScalar D("diffusivity",dimViscosity,cfg),mu("mu",dimDynamicViscosity,readScalar(cfg.lookup("viscosity")));
    dimensionedScalar refDensity("rhoReference",dimDensity,cfg.lookupOrDefault<scalar>("hydraulicReferenceDensity",inletProperties.density));dimensionedVector gravity("gravity",dimAcceleration,cfg.lookupOrDefault<vector>("gravity",vector(0,-9.80665,0)));
    if(!(refDensity.value()>0)||!std::isfinite(refDensity.value()))throw std::invalid_argument("hydraulicReferenceDensity must be positive");
    const vector pressureReferencePoint=cfg.lookupOrDefault<vector>("pressureReferencePoint",vector::zero);
    volScalarField p(IOobject("p",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh),eps(IOobject("epsilon",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh),inert(IOobject("inert",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh);
    volVectorField U(IOobject("U",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh);
    PtrList<volScalarField> solid(7),Y(7),source(7);
    for(int j=0;j<7;j++){
        solid.set(j,new volScalarField(IOobject("mineral"+name(j),runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh));
        Y.set(j,new volScalarField(IOobject("Y"+name(j),runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh));
        source.set(j,new volScalarField(IOobject("source"+name(j),runTime.timeName(),mesh),mesh,dimensionedScalar(dimDensity/dimTime,0)));mesh.setFluxRequired(Y[j].name());
    }
    volScalarField rho(IOobject("rho",runTime.timeName(),mesh,restarting?IOobject::MUST_READ:IOobject::NO_READ,IOobject::AUTO_WRITE),eps*refDensity);if(!restarting)rho=refDensity;
    volScalarField totalSource(IOobject("totalSource",runTime.timeName(),mesh),mesh,dimensionedScalar(dimDensity/dimTime,0)),storage(IOobject("massStorage",runTime.timeName(),mesh),mesh,dimensionedScalar(dimDensity/dimTime,0)),drag(IOobject("drag",runTime.timeName(),mesh),mesh,dimensionedScalar(dimless/dimTime,0));
    surfaceScalarField phi(IOobject("phi",runTime.timeName(),mesh,restarting?IOobject::MUST_READ:IOobject::NO_READ,IOobject::AUTO_WRITE),fvc::flux(U)),massPhi(IOobject("massPhi",runTime.timeName(),mesh,restarting?IOobject::MUST_READ:IOobject::NO_READ,IOobject::AUTO_WRITE),fvc::interpolate(rho)*phi);mesh.setFluxRequired(p.name());
    // Hydraulic port construction replaces patch objects. Retain the accepted
    // numerical values until the first ordinary transient pressure correction.
    saltflow::HydraulicBoundaries hydraulics(cfg,mesh,inletFractions,inletChemistry,p,U,Y);
    if(restarting)restart.transferFields(true,p,eps,rho,inert,U,phi,massPhi,Y,solid);
    if(flowPressure&&densityPressureMode=="thermodynamic"&&p.needReference())
        throw std::invalid_argument("Thermodynamic density pressure requires an imposed pressure/head boundary; use local pressure for closed domains");
    label pressureReferenceCell=-1;scalar pressureReferenceValue=0;
    setRefCell(p,mesh.solutionDict().subDict("SIMPLE"),pressureReferenceCell,pressureReferenceValue);
    volScalarField chemicalPressure(IOobject("chemicalPressureAtm",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),mesh,dimensionedScalar(dimless,initialTP[1]));
    volScalarField densityPressureDerivative(IOobject("densityPressureDerivative",runTime.timeName(),mesh),mesh,dimensionedScalar(dimDensity/dimPressure,0));
    const auto& centres=mesh.C();const double bulkVolume=gSum(mesh.V());PhaseStencil stencil(mesh,dx,dy);
    saltflow::MineralSurface mineralSurface(cfg,mesh,stencil,dx,dy,epMin);
#ifdef SALT_REFERENCE_BACKEND
    if(mineralSurface.enabled)throw std::invalid_argument("Heterogeneous interfaces currently require the native Pitzer backend");
#endif
    int invalidInitial=0;
    forAll(eps,i){
        double phase=eps[i]+inert[i],sumY=0;
        bool valid=std::isfinite(eps[i])&&eps[i]>=epMin-1e-12&&eps[i]<=1&&std::isfinite(inert[i])&&inert[i]>=0&&std::isfinite(p[i]);
        for(int j=0;j<3;j++)valid=valid&&std::isfinite(U[i][j]);
        for(int j=0;j<7;j++){
            phase+=solid[j][i];sumY+=Y[j][i];
            valid=valid&&std::isfinite(solid[j][i])&&solid[j][i]>=0&&std::isfinite(Y[j][i])&&Y[j][i]>=0;
        }
        invalidInitial|=!valid||std::abs(phase-1)>1e-10||!(sumY<1);
    }
    reduce(invalidInitial,maxOp<int>());
    if(invalidInitial)throw std::invalid_argument("Initial fields require finite nonnegative inventories, positive solvent, and epsilon + inert + minerals = 1");
    double chemistryTime=0,flowTime=0,transportTime=0,outputTime=0,initializationTime=0,flowChange=0,couplingChange=0,massResidual=0,momentumResidual=0,maxBalance=0;
    double chemicalAccepted=0,chemicalEvaluations=0,chemicalRestarts=0,chemicalFallbacks=0;
    double chemicalJacobians=0,jacobianFallbacks=0,transportAssemblies=0,positivityMass=0,equationResidual=0;
    double propertyTime=0,flowSolves=0,interfaceGeometryTime=0;
    std::array<double,8> boundaryIntegral{},initialInventory{},balance{};saltchem::MineralVector dissolved{},precipitated{};
    std::vector<saltchem::AqueousState> properties(mesh.nCells());std::vector<saltchem::Totals> cachedY(mesh.nCells());std::vector<unsigned char> cached(mesh.nCells(),0);
    struct BoundaryDensity {saltchem::Totals fractions{};double pressure=0,temperature=0,density=0;bool valid=false;};
    const bool cacheBoundaryProperties=cfg.lookupOrDefault<bool>("cacheBoundaryProperties",true);
    std::vector<std::vector<BoundaryDensity>> boundaryDensity(mesh.boundary().size());
    forAll(mesh.boundary(),k)if(!mesh.boundary()[k].coupled()&&mesh.boundary()[k].type()!="empty")boundaryDensity[k].resize(mesh.boundary()[k].size());
    double boundaryPropertyEvaluations=0,boundaryPropertyHits=0;
    std::vector<saltchem::Inventory> cells(mesh.nCells());std::vector<saltchem::RateParameters> rates(mesh.nCells());
    std::vector<double> cachedPressure(mesh.nCells(),0),chemistryCost(mesh.nCells(),1),cellSeconds(mesh.nCells(),0),cellEvaluations(mesh.nCells(),0),cellRejected(mesh.nCells(),0),cellFallbacks(mesh.nCells(),0);
    std::vector<label> workOrder(mesh.nCells());std::iota(workOrder.begin(),workOrder.end(),0);
    double skippedChemistry=0,couplingIterations=0,pressureOffset=0,eosLinearized=0,maxDensityPressureError=0,minimumUsedCouplingRelaxation=1,maxReservoirFluxError=0;
    // Preserve the small mechanical pressure across iterations. Repeatedly
    // adding/subtracting a large imposed pressure destroys its low bits.
    scalarField pressureDeviation(p.primitiveField());
    auto absolutePressure=[&](double reduced,const vector& x){return model.pressureAtm()+(flowPressure?(reduced+refDensity.value()*(gravity.value()&(x-pressureReferencePoint)))/101325.:0);};
    // Low-Mach option: retain the imposed thermodynamic pressure and its
    // hydrostatic part in the EOS, while mechanical pressure enforces mass
    // conservation. Chemical equilibria and kinetics still use local total
    // absolute pressure. This removes acoustic storage, not solutal density.
    auto densityPressure=[&](double reduced,const vector& x){return absolutePressure(densityPressureMode=="local"?reduced:pressureOffset,x);};
    auto updateHydraulics=[&](){
        const double previousOffset=pressureOffset;
        hydraulics.update(runTime.value(),eps,refDensity.value(),gravity.value(),pressureReferencePoint,p,U,Y,spec);
        double area=0,sumPressure=0;forAll(p.boundaryField(),k)if(p.boundaryField()[k].fixesValue()){
            area+=sum(mesh.magSf().boundaryField()[k]);sumPressure+=sum(p.boundaryField()[k]*mesh.magSf().boundaryField()[k]);
        }
        reduce(area,sumOp<double>());reduce(sumPressure,sumOp<double>());pressureOffset=area>0?sumPressure/area:pressureReferenceValue;
        pressureDeviation+=previousOffset-pressureOffset;
    };
    auto failCollectively=[&](const std::vector<std::string>& errors,const char* stage){int failed=0;for(label i=0;i<label(errors.size());i++)if(!errors[i].empty()){if(failed++<3)Pout<<stage<<" index="<<i<<" "<<errors[i]<<endl;}reduce(failed,sumOp<int>());if(failed)throw std::runtime_error(std::string(stage)+" failed in "+std::to_string(failed)+" cells/faces");};
    auto updateReservoirs=[&](){hydraulics.updateTransport(massPhi,mineralSurface.diffusivity(D,eps,rho)(),Y);for(auto& y:Y)y.correctBoundaryConditions();};
    auto updateProperties=[&](bool updateCompressibility=false,bool exact=true){
        const auto propertyStart=Timer::now();
        updateReservoirs();
        std::vector<std::string> errors(mesh.nCells());
        double linearized=0;
        #pragma omp parallel for schedule(static) reduction(+:linearized)
        for(label i=0;i<mesh.nCells();i++)try{
            chemicalPressure[i]=absolutePressure(p[i],centres[i]);const double pa=densityPressure(p[i],centres[i]);
            bool changed=!cached[i]||pa!=cachedPressure[i];double compositionChange=0;saltchem::Inventory state;state.waterEquivalent=1;
            for(int j=0;j<7;j++){changed=changed||Y[j][i]!=cachedY[i][j];compositionChange=max(compositionChange,std::abs(Y[j][i]-cachedY[i][j]));state.waterEquivalent-=Y[j][i];state.aqueous[j]=Y[j][i]/saltchem::componentMass[j];}
            if(!exact&&linearizeEOS&&cached[i]&&pa==cachedPressure[i]&&compositionChange<=1e-12){
                // The same inexact-iterate allowance also applies when
                // thermodynamic pressure is fixed. Always restore the full
                // composition-dependent EOS before accepting a time segment.
                rho[i]=properties[i].density;linearized++;continue;
            }
            if(!exact&&linearizeEOS&&pressureCompressible&&cached[i]&&pa>0&&compositionChange<=1e-12&&std::abs(pa-cachedPressure[i])*101325.<=1&&densityPressureDerivative[i]>0){
                // Inexact nonlinear iterate only. Every acceptance candidate
                // restores the full EOS and rechecks mass and transport.
                rho[i]=properties[i].density+densityPressureDerivative[i]*(pa-cachedPressure[i])*101325.;linearized++;continue;
            }
            if(changed){double water;auto& worker=localChemistry.at(model.temperatureK(),pa);properties[i]=worker.speciation.conserved(state,water,cached[i]?&properties[i]:nullptr);for(int j=0;j<7;j++)cachedY[i][j]=Y[j][i];cachedPressure[i]=pa;cached[i]=true;}
            rho[i]=properties[i].density;
            if(pressureCompressible&&updateCompressibility){
                // A lagged EOS Jacobian stabilizes compressible storage. It is
                // only a nonlinear preconditioner: the final mass residual
                // still uses the unmodified EOS and conservative FV fluxes.
                // Refresh once per segment, not in every coupling iteration.
                const double dp=max(1e-3,pa*1e-5);double water;
                auto& perturbed=localChemistry.at(model.temperatureK(),pa+dp);
                const double high=perturbed.speciation.conserved(state,water,&properties[i]).density;
                densityPressureDerivative[i]=max(0.,(high-rho[i])/(dp*101325.));
            }
        }catch(const std::exception& e){errors[i]=e.what();}
        eosLinearized+=linearized;
        failCollectively(errors,"EOS");rho.correctBoundaryConditions();
        // Algebraic construction eps*rhoReference creates calculated physical
        // patches. correctBoundaryConditions() alone leaves their old density
        // unchanged. Evaluate the EOS at the actual boundary composition and
        // absolute pressure, including outflow and stagnant wall faces.
        std::vector<std::string> boundaryErrors(mesh.boundary().size());
        forAll(mesh.boundary(),patch)if(!mesh.boundary()[patch].coupled()&&mesh.boundary()[patch].type()!="empty")try{
            auto& rb=rho.boundaryFieldRef()[patch];auto& cp=chemicalPressure.boundaryFieldRef()[patch];const auto& cf=mesh.Cf().boundaryField()[patch];
            forAll(rb,f){
                cp[f]=absolutePressure(p.boundaryField()[patch][f],cf[f]);saltchem::Inventory state;state.waterEquivalent=1;
                const double pa=densityPressure(p.boundaryField()[patch][f],cf[f]);auto& cachedBoundary=boundaryDensity[patch][f];
                bool hit=cacheBoundaryProperties&&cachedBoundary.valid&&cachedBoundary.pressure==pa&&cachedBoundary.temperature==model.temperatureK();
                for(int j=0;j<7;j++){const double fraction=Y[j].boundaryField()[patch][f];hit=hit&&cachedBoundary.fractions[j]==fraction;state.waterEquivalent-=fraction;state.aqueous[j]=fraction/saltchem::componentMass[j];}
                if(hit){rb[f]=cachedBoundary.density;boundaryPropertyHits++;continue;}
                double water;auto& worker=localChemistry.at(model.temperatureK(),pa);rb[f]=worker.speciation.conserved(state,water,&properties[mesh.boundary()[patch].faceCells()[f]]).density;
                boundaryPropertyEvaluations++;cachedBoundary.valid=true;cachedBoundary.pressure=pa;cachedBoundary.temperature=model.temperatureK();cachedBoundary.density=rb[f];
                for(int j=0;j<7;j++)cachedBoundary.fractions[j]=Y[j].boundaryField()[patch][f];
            }
        }catch(const std::exception& e){boundaryErrors[patch]=e.what();}
        failCollectively(boundaryErrors,"BOUNDARY_EOS");
        propertyTime+=elapsed(propertyStart);
    };
    auto loadCells=[&](){
        const auto geometryStart=Timer::now();mineralSurface.reconstruct(solid,inert,rateTemplate);mineralSurface.updateTransport(eps);if(mineralSurface.enabled)interfaceGeometryTime+=elapsed(geometryStart);
        const auto neighbours=stencil.values(eps);
        forAll(eps,i){auto& s=cells[i];s.waterEquivalent=eps[i]*rho[i];for(int j=0;j<7;j++){s.aqueous[j]=eps[i]*rho[i]*Y[j][i]/saltchem::componentMass[j];s.waterEquivalent-=eps[i]*rho[i]*Y[j][i];s.mineral[j]=solid[j][i]/model.molarVolume[j];}
            auto& r=rates[i];r=rateTemplate;r.geometry.inertFraction=inert[i];r.geometry.neighbours=neighbours[i];mineralSurface.configure(i,s,r,model.molarVolume);}
    };
    auto inventory=[&](){std::array<double,8> a{};forAll(eps,i){const double cap=eps[i]*rho[i],v=mesh.V()[i];double sumY=0;for(int j=0;j<7;j++){sumY+=Y[j][i];a[j]+=cap*Y[j][i]*v/saltchem::componentMass[j];}a[7]+=cap*(1-sumY)*v;for(int k=0;k<7;k++){const double m=solid[k][i]*v/model.molarVolume[k];for(int j=0;j<7;j++)a[j]+=saltchem::stoichiometry[k][j]*m;a[7]+=saltchem::hydration[k]*saltchem::waterMolarMass*m;}}for(double& x:a)reduce(x,sumOp<double>());return a;};
    auto flowIteration=[&](bool transient,const scalarField& oldMass){
        flowSolves++;
        auto tick=Timer::now();vectorField previous(U.primitiveField());volVectorField previousVelocity("previousVelocity",U);
        // The pressure equation and gradients use a small gauge field. A large
        // uniform imposed pressure belongs in chemistry, not in subtraction of
        // nearly equal FV matrix terms. Put the datum back before EOS/output.
        auto shiftPressure=[&](double delta){
            p.primitiveFieldRef()+=delta;
            // fixedValue deliberately ignores operator += / -=; force only
            // the numerical boundary values, preserving each boundary type.
            for(auto& b:p.boundaryFieldRef())b==scalarField(b)+delta;
        };
        shiftPressure(-pressureOffset);p.primitiveFieldRef()=pressureDeviation;p.storePrevIter();
        if(transient)storage.primitiveFieldRef()=(eps.primitiveField()*rho.primitiveField()-oldMass)/runTime.deltaTValue();else storage.primitiveFieldRef()=0;
        forAll(eps,i){
            if(anchorWidth==0)drag[i]=eps[i]<=anchorThreshold?anchorStrength:0;
            else{
                // Interpolate mobility across the anchoring threshold. A hard
                // 1e12 -> 0 resistance jump launches a numerical transient
                // every time an interface cell opens. Harmonic blending
                // keeps the solid endpoint anchored and the open endpoint
                // unpenalized, with continuous cell-scale mobility between.
                const double s=std::clamp((eps[i]-anchorThreshold)/anchorWidth+.5,0.,1.);
                const double open=s*s*(3-2*s),cellDrag=mu.value()/(rho[i]*eps[i]*min(dx*dx,dy*dy));
                drag[i]=anchorStrength*(1-open)/(1+anchorStrength*open/cellDrag);
            }
        }drag.correctBoundaryConditions();
        fvVectorMatrix ue((1./eps)*fvm::div(massPhi/fvc::interpolate(eps),U)-fvm::laplacian(mu/eps,U)+fvm::Sp(rho*drag,U));if(transient)ue+=(1./eps)*fvm::ddt(rho,U);
        // Preserve the matching face source when momentum is under-relaxed.
        // Otherwise its artificial diagonal makes the converged Rhie-Chow
        // flux depend on the relaxation factor and the time-step size.
        volScalarField unrelaxedA(ue.A());
        ue.relax();volVectorField buoyancy((rho-refDensity)*gravity);solve(ue==-fvc::grad(p)+buoyancy);
        volScalarField rAU(1./ue.A()),rAtU(rAU);
        // SIMPLE and SIMPLEC solve the same converged equations. The latter
        // is an optional pressure preconditioner, not a required model term.
        if(consistentMomentum)rAtU=1./max(1./rAU-ue.H1(),.1/rAU);
        volVectorField H(constrainHbyA(rAU*(ue.H()+buoyancy),U,p));surfaceScalarField rhoFace(fvc::interpolate(rho));
        surfaceScalarField rAUf(1./fvc::interpolate(1./rAU)),rAtUf(1./fvc::interpolate(1./rAtU));
        volScalarField relaxationA(1./rAU-unrelaxedA);
        volVectorField spatialH(H/rAU-relaxationA*previousVelocity);
        volScalarField timeA("timeDiagonal",0.*unrelaxedA);
        volScalarField oldTimeA("oldTimeDiagonal",0.*unrelaxedA);
        if(transient){
            timeA=rho/eps/runTime.deltaT();oldTimeA=rho.oldTime()/eps/runTime.deltaT();
            spatialH-=oldTimeA*U.oldTime();
        }
        surfaceScalarField historyFlux(fvc::interpolate(relaxationA)*massPhi);
        if(transient)historyFlux+=rhoFace*fvc::interpolate(oldTimeA)*massPhi.oldTime()/fvc::interpolate(rho.oldTime());
        surfaceScalarField ph("predictedMassFlux",rAUf*(rhoFace*fvc::flux(spatialH)+historyFlux));
        if(faceMomentumInterpolation=="linear"){
            // Keep the original linear interpolation of spatial mobility.
            // Interpolate the temporal and relaxation terms after division
            // by that same spatial diagonal. At a stationary fixed point,
            // face histories cancel the denominator exactly, independent
            // of dt and relaxation, without changing the spatial stencil.
            volScalarField spatialA(unrelaxedA-timeA);
            if(gMin(spatialA.primitiveField())<=0)throw std::runtime_error("Non-positive spatial momentum diagonal");
            surfaceScalarField denominator(1.+fvc::interpolate((relaxationA+timeA)/spatialA));
            rAUf=fvc::interpolate(1./spatialA)/denominator;
            rAtUf=rAUf*fvc::interpolate(rAtU)/fvc::interpolate(rAU);
            ph=(rhoFace*fvc::flux(spatialH/spatialA)+fvc::interpolate(relaxationA/spatialA)*massPhi)/denominator;
            if(transient)ph+=rhoFace*fvc::interpolate(oldTimeA/spatialA)*massPhi.oldTime()/fvc::interpolate(rho.oldTime())/denominator;
        }
        ph+=rhoFace*(rAtUf-rAUf)*fvc::snGrad(p)*mesh.magSf();H-=(rAU-rAtU)*fvc::grad(p);
        forAll(U.boundaryField(),patch)if(U.boundaryField()[patch].fixesValue()&&!p.boundaryField()[patch].fixesValue())ph.boundaryFieldRef()[patch]==rho.boundaryField()[patch]*(U.boundaryField()[patch]&mesh.Sf().boundaryField()[patch]);
        surfaceScalarField rf(rhoFace*rAtUf);constrainPressure(p,rho,U,ph,rf);
        volScalarField pressureStorage(eps*densityPressureDerivative/runTime.deltaT());
        if(!transient)pressureStorage=dimensionedScalar(pressureStorage.dimensions(),0);
        fvScalarMatrix pe(fvm::laplacian(rf,p)-fvm::Sp(pressureStorage,p)==fvc::div(ph)+storage-totalSource-pressureStorage*p);
        // Compressible storage fixes absolute pressure even in a closed
        // domain. Pinning one cell there would create a spurious mass source.
        if(!transient||!pressureCompressible)pe.setReference(pressureReferenceCell,pressureReferenceValue-pressureOffset);
        pe.solve();massPhi=ph-pe.flux();phi=massPhi/rhoFace;
        p.relax(pressureRelaxation);U=H-rAtU*fvc::grad(p);pressureDeviation=p.primitiveField();shiftPressure(pressureOffset);U.correctBoundaryConditions();flowChange=gMax(mag(U.primitiveField()-previous))/max(gMax(mag(U.primitiveField())),1e-12);
        massResidual=gMax(mag((fvc::div(massPhi)+storage-totalSource)().primitiveField()));flowTime+=elapsed(tick);
    };
    auto checkPressureApproximation=[&](){
        if(densityPressureMode=="local"||!flowPressure)return;
        double error=0;std::vector<std::string> errors(mesh.nCells());
        #pragma omp parallel for schedule(static) reduction(max:error)
        for(label i=0;i<mesh.nCells();i++)try{
            saltchem::Inventory state;state.waterEquivalent=1;
            for(int j=0;j<7;j++){state.waterEquivalent-=Y[j][i];state.aqueous[j]=Y[j][i]/saltchem::componentMass[j];}
            double water;auto& worker=localChemistry.at(model.temperatureK(),chemicalPressure[i]);
            const auto full=worker.speciation.conserved(state,water,&properties[i]);
            error=std::max(error,std::abs(full.density/rho[i]-1));
        }catch(const std::exception& e){errors[i]=e.what();}
        failCollectively(errors,"LOCAL_PRESSURE_CHECK");
        std::vector<std::string> boundaryErrors(mesh.boundary().size());
        forAll(mesh.boundary(),patch)if(!mesh.boundary()[patch].coupled()&&mesh.boundary()[patch].type()!="empty")try{
            forAll(rho.boundaryField()[patch],f){
                saltchem::Inventory state;state.waterEquivalent=1;
                for(int j=0;j<7;j++){const double y=Y[j].boundaryField()[patch][f];state.waterEquivalent-=y;state.aqueous[j]=y/saltchem::componentMass[j];}
                double water;auto& worker=localChemistry.at(model.temperatureK(),chemicalPressure.boundaryField()[patch][f]);
                const auto full=worker.speciation.conserved(state,water,&properties[mesh.boundary()[patch].faceCells()[f]]);
                error=std::max(error,std::abs(full.density/rho.boundaryField()[patch][f]-1));
            }
        }catch(const std::exception& e){boundaryErrors[patch]=e.what();}
        failCollectively(boundaryErrors,"BOUNDARY_PRESSURE_CHECK");reduce(error,maxOp<double>());
        maxDensityPressureError=max(maxDensityPressureError,error);
        if(error>densityPressureTolerance)throw std::runtime_error("Thermodynamic-pressure EOS exceeds densityPressureTolerance; use densityPressureMode local");
    };
    auto checkMomentum=[&](bool transient){
        // A small iterate change alone can be caused by an under-solved
        // linear system. Test the unrelaxed momentum equation at the final
        // pressure, density and velocity, including processor boundaries.
        auto tick=Timer::now();
        volScalarField gaugePressure(p-dimensionedScalar("pressureDatum",dimPressure,pressureOffset));
        gaugePressure.primitiveFieldRef()=pressureDeviation;
        fvVectorMatrix equation((1./eps)*fvm::div(massPhi/fvc::interpolate(eps),U)-fvm::laplacian(mu/eps,U)+fvm::Sp(rho*drag,U));
        if(transient)equation+=(1./eps)*fvm::ddt(rho,U);
        equation+=fvc::grad(gaugePressure)-(rho-refDensity)*gravity;
        // Count each processor-face term once. Do not add the coupled
        // boundary source and then apply an LDU interface residual again.
        scalarField squaredError(mesh.nCells(),0);
        for(direction component=0;component<3;component++){
            scalarField diagonal(equation.diag());
            scalarField r(equation.source().component(component)-diagonal*U.primitiveField().component(component));
            const auto& lower=equation.lduAddr().lowerAddr();const auto& upper=equation.lduAddr().upperAddr();
            forAll(lower,f){r[lower[f]]-=equation.upper()[f]*U[upper[f]][component];r[upper[f]]-=equation.lower()[f]*U[lower[f]][component];}
            forAll(U.boundaryField(),patch){
                const auto& b=U.boundaryField()[patch];const auto& faceCells=mesh.boundary()[patch].faceCells();
                const auto& ic=equation.internalCoeffs()[patch];const auto& bc=equation.boundaryCoeffs()[patch];
                vectorField other;if(b.coupled())other=b.patchNeighbourField();
                forAll(faceCells,f){const label c=faceCells[f];diagonal[c]+=ic[f][component];r[c]-=ic[f][component]*U[c][component];r[c]+=bc[f][component]*(b.coupled()?other[f][component]:1);}
            }
            squaredError+=sqr(r/diagonal);
        }
        momentumResidual=std::sqrt(gMax(squaredError))/max(gMax(mag(U.primitiveField())),1e-12);
        flowTime+=elapsed(tick);return momentumResidual;
    };
    auto ledger=[&](double dt){
        auto tick=Timer::now();surfaceScalarField deff(mineralSurface.diffusivity(D,eps,rho));double soluteSum=0,totalFlux=0;
        for(int j=0;j<7;j++){fvScalarMatrix matrix(fvm::div(massPhi,Y[j])-fvm::laplacian(deff,Y[j]));surfaceScalarField flux(matrix.flux());maxReservoirFluxError=max(maxReservoirFluxError,hydraulics.fluxError(j,massPhi,flux,deff,Y[j]));double componentFlux=0;forAll(flux.boundaryField(),patch)if(!flux.boundaryField()[patch].coupled()&&mesh.boundary()[patch].type()!="empty")componentFlux+=sum(flux.boundaryField()[patch]);reduce(componentFlux,sumOp<double>());boundaryIntegral[j]+=dt*componentFlux/saltchem::componentMass[j];soluteSum+=componentFlux;}
        forAll(massPhi.boundaryField(),patch)if(!massPhi.boundaryField()[patch].coupled()&&mesh.boundary()[patch].type()!="empty")totalFlux+=sum(massPhi.boundaryField()[patch]);reduce(totalFlux,sumOp<double>());boundaryIntegral[7]+=dt*(totalFlux-soluteSum);
        auto now=inventory();for(int j=0;j<8;j++){const double scale=max(std::abs(initialInventory[j]),j==7?bulkVolume*1e-6:bulkVolume*1e-3);balance[j]=(now[j]-initialInventory[j]+boundaryIntegral[j])/scale;maxBalance=max(maxBalance,std::abs(balance[j]));}transportTime+=elapsed(tick);
    };
    fileName out=runTime.path()/"output";mkDir(out);mkDir(out/"frames");std::ofstream hist((out/"history.csv").c_str());hist<<std::setprecision(14)<<"time_s,porosity,Qin_m3_s,Qout_m3_s";
    std::ofstream boundaryHistory((out/"boundary_history.csv").c_str());boundaryHistory<<std::setprecision(16)<<"time_s,patch,Qin_m3_s,Qout_m3_s,net_mass_kg_s,mean_absolute_pressure_atm\n";
    for(const char* prefix:{"min_","prec_","diss_"})for(const char* n:saltchem::mineralNames)hist<<','<<prefix<<n;hist<<",meanRho,minRho,maxRho,maxP_Pa,maxU_m_s";for(const char* n:saltchem::componentNames)hist<<",balance_"<<n;hist<<",balance_water,mass_residual,flow_change,coupling_change,chemistry_s,flow_s,transport_s\n";
    auto save=[&](){
        if(mineralSurface.fluidSide()){
            std::vector<std::string> diagnosticErrors(mesh.nCells());
            #pragma omp parallel for schedule(static)
            for(label i=0;i<mesh.nCells();i++)try{
                auto& worker=localChemistry.at(model.temperatureK(),absolutePressure(p[i],centres[i]));
                mineralSurface.updateFilmDiagnostics(i,worker.reactions,cells[i],rateTemplate,eps[i]);
            }catch(const std::exception& error){diagnosticErrors[i]=error.what();}
            failCollectively(diagnosticErrors,"SURFACE_FILM_DIAGNOSTIC");
        }
        mineralSurface.save(out);
        auto tick=Timer::now();saltchem::MineralVector minerals{};for(int j=0;j<7;j++)minerals[j]=gSum(solid[j].primitiveField()*mesh.V())/model.molarVolume[j];const double por=gSum(eps.primitiveField()*mesh.V())/bulkVolume;
        auto globalPrecipitated=precipitated,globalDissolved=dissolved;for(int j=0;j<7;j++){reduce(globalPrecipitated[j],sumOp<double>());reduce(globalDissolved[j],sumOp<double>());}
        double qin=0,qout=0;forAll(phi.boundaryField(),patch)if(!mesh.boundary()[patch].coupled()&&mesh.boundary()[patch].type()!="empty")for(double q:phi.boundaryField()[patch]){qin-=min(q,0.);qout+=max(q,0.);}reduce(qin,sumOp<double>());reduce(qout,sumOp<double>());
        forAll(mesh.boundary(),patch)if(!mesh.boundary()[patch].coupled()&&mesh.boundary()[patch].type()!="empty"){
            double incoming=0,outgoing=0,mass=sum(massPhi.boundaryField()[patch]);
            for(double q:phi.boundaryField()[patch]){incoming-=min(q,0.);outgoing+=max(q,0.);}
            double area=sum(mesh.magSf().boundaryField()[patch]),pressure=sum(chemicalPressure.boundaryField()[patch]*mesh.magSf().boundaryField()[patch]);
            for(double* v:{&incoming,&outgoing,&mass,&area,&pressure})reduce(*v,sumOp<double>());
            if(Pstream::master()&&area>0)boundaryHistory<<runTime.value()<<','<<mesh.boundary()[patch].name()<<','<<incoming<<','<<outgoing<<','<<mass<<','<<pressure/area<<'\n';
        }
        boundaryHistory.flush();
        hist<<runTime.value()<<','<<por<<','<<qin<<','<<qout;for(const auto& a:{minerals,globalPrecipitated,globalDissolved})for(double v:a)hist<<','<<v;
        hist<<','<<gSum(eps.primitiveField()*rho.primitiveField()*mesh.V())/(por*bulkVolume)<<','<<gMin(rho.primitiveField())<<','<<gMax(rho.primitiveField())<<','<<gMax(mag(p.primitiveField()))<<','<<gMax(mag(U.primitiveField()));for(double b:balance)hist<<','<<b;
        hist<<','<<massResidual<<','<<flowChange<<','<<couplingChange<<','<<chemistryTime<<','<<flowTime<<','<<transportTime<<std::endl;
        // Fixed 33-column little-endian float32 snapshots; schema in reader script.
        std::ofstream f((out/"frames"/(runTime.timeName()+".f32")).c_str(),std::ios::binary);std::uint32_t header[4]={0x53414C54,1,std::uint32_t(mesh.nCells()),33};f.write(reinterpret_cast<char*>(header),sizeof(header));
        forAll(eps,i){auto chemical=properties[i];if(densityPressureMode!="local"&&flowPressure){double water;auto& worker=localChemistry.at(model.temperatureK(),chemicalPressure[i]);chemical=worker.speciation.conserved(cells[i],water,&properties[i]);}std::array<float,33> row{};int n=0;for(double v:{centres[i].x(),centres[i].y(),eps[i],U[i].x(),U[i].y(),p[i]})row[n++]=float(v);for(int j=0;j<7;j++)row[n++]=float(solid[j][i]);row[n++]=float(inert[i]);row[n++]=float(rho[i]);for(double v:saltchem::componentTotals(chemical.molality))row[n++]=float(v);for(double v:chemical.saturationIndex)row[n++]=float(v);row[n++]=float(chemical.pH);row[n++]=float(chemical.waterActivity);
            double Q=0;saltchem::MineralVector mineralRates{};
            if(!mineralSurface.fluidSide())mineralRates=saltchem::surfaceRates(cells[i],rates[i],chemical.saturationIndex,chemical.logActivity[saltchem::H],chemical.logActivity[saltchem::CO2],model.molarVolume,&Q);
            if(mineralSurface.enabled){const auto q=mineralSurface.transformation(i,mineralRates,model.molarVolume);row[n++]=float(q[0]);row[n++]=float(q[1]);}
            else{const auto qv=phaseVector(eps[i],rates[i].geometry.neighbours,dx,dy,Q);row[n++]=float(qv[1]);row[n++]=float(qv[2]);}
            f.write(reinterpret_cast<const char*>(row.data()),sizeof(row));}
        if(!f)FatalErrorInFunction<<"Snapshot write failed"<<exit(FatalError);
        if(runTime.value()+1e-10>=runTime.endTime().value()){
            std::ofstream exact((out/"final_state.f64").c_str(),std::ios::binary);std::uint32_t h[4]={0x53414C54,2,std::uint32_t(mesh.nCells()),22};exact.write(reinterpret_cast<char*>(h),sizeof(h));
            forAll(eps,i){std::array<double,22> row{};int j=0;for(double v:{centres[i].x(),centres[i].y(),eps[i],U[i].x(),U[i].y(),p[i],rho[i],cells[i].waterEquivalent})row[j++]=v;for(double v:cells[i].aqueous)row[j++]=v;for(double v:cells[i].mineral)row[j++]=v;exact.write(reinterpret_cast<const char*>(row.data()),sizeof(row));}
            if(!exact)FatalErrorInFunction<<"Final-state write failed"<<exit(FatalError);
        }
        outputTime+=elapsed(tick);Info<<"FRAME "<<runTime.value()<<" porosity="<<por<<" carnallite="<<minerals[2]<<" balance="<<maxBalance<<endl;
    };
    updateHydraulics();
    // Lift the initial pressure guess to the imposed pressure level. Starting
    // an elevated-pressure domain at p=0 would launch a spurious huge velocity
    // before the first pressure correction. Users can retain their own guess.
    if(!restarting&&cfg.lookupOrDefault<bool>("initializePressure",true)){
        double area=0,pressureSum=0;forAll(p.boundaryField(),k)if(p.boundaryField()[k].fixesValue()){
            area+=sum(mesh.magSf().boundaryField()[k]);pressureSum+=sum(p.boundaryField()[k]*mesh.magSf().boundaryField()[k]);
        }
        reduce(area,sumOp<double>());reduce(pressureSum,sumOp<double>());
        if(area>0)p.primitiveFieldRef()=pressureSum/area;
        forAll(p.boundaryField(),k)if(isA<fixedFluxPressureFvPatchScalarField>(p.boundaryField()[k]))
            refCast<fixedFluxPressureFvPatchScalarField>(p.boundaryFieldRef()[k]).updateSnGrad(scalarField(mesh.boundary()[k].size(),0));
        p.correctBoundaryConditions();
    }
    if(restarting){
        restart.transferFields(true,p,eps,rho,inert,U,phi,massPhi,Y,solid);
        pressureDeviation=restart.deviation;pressureOffset=restart.control[FlowCheckpoint::offset];
        // Rebuild chemistry caches without replacing the accepted density or
        // boundary values. Their next ordinary update occurs inside a step.
        loadCells();updateProperties();
        restart.transferFields(true,p,eps,rho,inert,U,phi,massPhi,Y,solid);
        initialInventory=restart.initial;boundaryIntegral=restart.boundary;balance=restart.balance;
        maxBalance=restart.control[FlowCheckpoint::maxBalance];dissolved=restart.dissolved;precipitated=restart.precipitated;
        const auto restored=inventory();for(int j=0;j<8;j++)if(restored[j]!=restart.inventory[j])
            throw std::runtime_error("Checkpoint inventory mismatch");
        Info<<"RESTART restored accepted transient state at time="<<runTime.value()<<endl;
    }else{
        pressureDeviation=p.primitiveField()-pressureOffset;
        U.correctBoundaryConditions();phi=fvc::flux(U);massPhi=fvc::interpolate(rho)*phi;for(auto& y:Y)y.correctBoundaryConditions();updateProperties();massPhi=fvc::interpolate(rho)*phi;loadCells();
    }
#ifdef SALT_REFERENCE_BACKEND
    const double referenceTolerance=cfg.lookupOrDefault<scalar>("referenceKineticTolerance",1e-13);
    const auto initStart=Timer::now();saltchem::SevenPhreeqcBackend reference(db.c_str(),cells,threads,true,referenceTolerance,initialTP[0],initialTP[1]);initializationTime=elapsed(initStart);Info<<"CHEMISTRY_BACKEND PhreeqcRM kinetic_tolerance="<<referenceTolerance<<endl;
#else
    Info<<"CHEMISTRY_BACKEND native LOCAL_INTEGRATOR "<<(options.method==2?"CVODE_BDF":options.method==1?"ROS4":"TR_BDF2")<<endl;
#endif
    Info<<"THERMODYNAMIC_CONDITIONS temperature_C="<<model.temperatureK()-273.15<<" absolute_pressure_atm="<<model.pressureAtm()<<endl;
    Info<<"PRESSURE_CONVENTION mode="<<(flowPressure?"flow":"uniform")<<" reference_point="<<pressureReferencePoint<<" gravity="<<gravity.value()<<" reference_density="<<refDensity.value()<<endl;
    Info<<"NUMERICAL_CONTROLS epsilonMin="<<epMin<<" poreClosureWidth="<<rateTemplate.geometry.closureWidth<<" phaseDirectionWidth="<<rateTemplate.geometry.directionWidth<<" waterPolicy="<<waterPolicy<<" reuseTransportMatrix="<<reuseTransport<<" analyticJacobian="<<options.analyticJacobian<<endl;
    scalarField initialMass(eps.primitiveField()*rho.primitiveField());int initialIterations=0;
    saltflow::CouplingAcceleration initialAccelerator(6,.8);
    double initialVelocityScale=1e-6;for(const auto& b:U.boundaryField())initialVelocityScale=max(initialVelocityScale,max(mag(b)));reduce(initialVelocityScale,maxOp<double>());
    const double initialPressureScale=max(gMax(mag(pressureDeviation)),1.);
    auto initialPack=[&](){scalarField x(4*mesh.nCells());forAll(eps,i){for(int j=0;j<3;j++)x[4*i+j]=U[i][j]/initialVelocityScale;x[4*i+3]=pressureDeviation[i]/initialPressureScale;}return x;};
    for(;!restarting&&initialIterations<maxIterations;initialIterations++){
        const scalarField previous(initialPack());flowIteration(false,initialMass);for(auto& y:Y)y.correctBoundaryConditions();const scalarField previousRho(rho.primitiveField());updateProperties();
        const double densityChange=gMax(mag(rho.primitiveField()-previousRho)/rho.primitiveField());
        if(initialIterations%100==0)Info<<"INITIAL_FLOW "<<initialIterations<<" change="<<flowChange<<endl;
        if(initialIterations>5&&flowChange<flowTol&&densityChange<couplingTol&&massResidual<1e-5&&checkMomentum(false)<flowTol)break;
        if(cfg.lookupOrDefault<bool>("accelerateInitialFlow",true)){
            const scalarField proposed(initialAccelerator.next(previous,initialPack()));
            forAll(eps,i){for(int j=0;j<3;j++)U[i][j]=proposed[4*i+j]*initialVelocityScale;pressureDeviation[i]=proposed[4*i+3]*initialPressureScale;p[i]=pressureOffset+pressureDeviation[i];}
            U.correctBoundaryConditions();p.correctBoundaryConditions();updateProperties();
        }
    }
    if(initialIterations==maxIterations)FatalErrorInFunction<<"Initial flow did not converge"<<exit(FatalError);checkPressureApproximation();if(!restarting)initialInventory=inventory();loadCells();save();
    double next=restarting?restart.control[FlowCheckpoint::nextFrame]:runTime.value()+frame;
    const saltflow::TimeStepControl timeControl(cfg,restarting?restart.control[FlowCheckpoint::maximumDt]:runTime.deltaTValue());
    double acceptedSteps=0,rejectedSteps=0,trialSegments=0,lastTimeError=0,maxTimeError=0,minAcceptedDt=GREAT,maxAcceptedDt=0;
    if(restarting){acceptedSteps=restart.control[FlowCheckpoint::acceptedSteps];rejectedSteps=restart.control[FlowCheckpoint::rejectedSteps];trialSegments=restart.control[FlowCheckpoint::trialSegments];lastTimeError=restart.control[FlowCheckpoint::lastTimeError];maxTimeError=restart.control[FlowCheckpoint::maxTimeError];minAcceptedDt=restart.control[FlowCheckpoint::minAcceptedDt];maxAcceptedDt=restart.control[FlowCheckpoint::maxAcceptedDt];}
    auto capture=[&](){
        saltflow::FlowState state(runTime,p,eps,rho,U,phi,massPhi,Y,solid);
        state.pressureDeviation=pressureDeviation;state.pressureOffset=pressureOffset;
        state.boundaryIntegral=boundaryIntegral;state.balance=balance;state.maxBalance=maxBalance;state.dissolved=dissolved;state.precipitated=precipitated;state.mineralSurface=mineralSurface.state;state.mineralSurfaceRates=mineralSurface.reactionRates;state.mineralSurfaceVectors=mineralSurface.reactionVectors;return state;
    };
    auto restore=[&](const saltflow::FlowState& state){
        static_cast<TimeState&>(runTime)=state.nativeTime;
        runTime.setTime(state.time,state.index);state.p.restore(p);state.eps.restore(eps);state.rho.restore(rho);state.U.restore(U);state.phi.restore(phi);state.massPhi.restore(massPhi);
        pressureDeviation=state.pressureDeviation;pressureOffset=state.pressureOffset;
        for(int j=0;j<7;j++){state.Y[j].restore(Y[j]);state.solid[j].restore(solid[j]);}
        boundaryIntegral=state.boundaryIntegral;balance=state.balance;maxBalance=state.maxBalance;dissolved=state.dissolved;precipitated=state.precipitated;
        mineralSurface.state=state.mineralSurface;
        mineralSurface.reactionRates=state.mineralSurfaceRates;
        mineralSurface.reactionVectors=state.mineralSurfaceVectors;
        std::fill(cached.begin(),cached.end(),0);
    };
    auto advanceSegment=[&](double dt){
        trialSegments++;scalarField oldMass(eps.primitiveField()*rho.primitiveField());
        // Subtracting two nearly equal mass densities then dividing by a very
        // small dt has a predictable double-precision floor. This allowance
        // is roundoff-scaled, not a relaxation of the time-error tolerances.
        const double massTolerance=max(1e-5,64*std::numeric_limits<double>::epsilon()*gMax(mag(oldMass))/dt);
        runTime.setDeltaT(dt,false);runTime.setTime(runTime.value()+dt,runTime.timeIndex()+1);U.oldTime()==U;rho.oldTime()==rho;massPhi.oldTime()==massPhi;
        // Preserve previous mass inventories. The new EOS enters the coupled
        // storage/pressure solve after the local reaction stage.
        const auto tp=thermalSchedule.at(runTime.value());
        if(model.setConditions(tp[0]+273.15,tp[1])){
            std::fill(cached.begin(),cached.end(),0);updateRateTemperature();
            updateInlet();
            double w;if(!cfg.found("hydraulicReferenceDensity"))refDensity.value()=spec.conserved(inletState,w).density;
#ifdef SALT_REFERENCE_BACKEND
            reference.setConditions(tp[0],tp[1]);
#endif
            Info<<"THERMODYNAMIC_CONDITIONS time="<<runTime.value()<<" temperature_C="<<tp[0]<<" absolute_pressure_atm="<<tp[1]<<endl;
        }
        updateHydraulics();
        std::array<scalarField,7> oldY;for(int j=0;j<7;j++)oldY[j]=Y[j].primitiveField();auto tick=Timer::now();loadCells();const auto before=cells;mineralSurface.beginReactions();std::vector<std::string> errors(mesh.nCells());
#ifdef SALT_REFERENCE_BACKEND
        try{std::vector<double> pa(mesh.nCells());forAll(p,i)pa[i]=absolutePressure(p[i],centres[i]);reference.setConditions(tp[0],tp[1]);reference.setPressures(pa);reference.setTime(runTime.value()-dt);reference.advance(cells,rates,dt);}catch(const std::exception& e){errors[0]=e.what();}
#else
        double stepAccepted=0,stepEvaluations=0,stepRestarts=0,stepFallbacks=0,stepJacobians=0,stepJacobianFallbacks=0,stepSkipped=0;
        auto advanceFaces=[&](double faceDt,bool precipitation){
            // Face matchings serialize every shared liquid/solid inventory.
            // Amounts use mol or kg per bulk m3, with uniform physical cell
            // volumes guaranteed by PhaseStencil. No solid numerical pore
            // water is borrowed by this surface reaction stage.
            using Packet=FixedList<scalar,58>;
            std::vector<Packet> packets(mesh.nCells());
            auto putInventory=[](Packet& packet,const saltchem::Inventory& cell){
                for(int j=0;j<7;j++){packet[j]=cell.aqueous[j];packet[7+j]=cell.mineral[j];}packet[14]=cell.waterEquivalent;
            };
            auto getInventory=[](const Packet& packet){
                saltchem::Inventory cell;for(int j=0;j<7;j++){cell.aqueous[j]=packet[j];cell.mineral[j]=packet[7+j];}cell.waterEquivalent=packet[14];return cell;
            };
            forAll(eps,i){
                auto& packet=packets[i];packet=0.;putInventory(packet,cells[i]);
                const auto caps=mineralSurface.lockedCaps(i,model.molarVolume);
                for(int j=0;j<7;j++){packet[15+j]=min(caps[j],cells[i].mineral[j]);packet[31+j]=mineralSurface.supportArea(i)[j];}
                const auto film=mineralSurface.filmParameters(i,eps[i]);
                for(int j=0;j<7;j++){packet[44+j]=film.coverage[j];packet[51+j]=std::isinf(film.coatedVelocity[j])?-1:film.coatedVelocity[j];}
                packet[22]=inert[i];packet[23]=absolutePressure(p[i],centres[i]);
            }
            const auto activeFaces=mineralSurface.activeReactionFaces(rateTemplate);
            auto& exchange=mineralSurface.surfaceExchange();
            const double cellVolume=gMin(mesh.V());
            exchange.execute(packets,activeFaces,[&](label face,Packet& own,Packet& other){
                auto react=[&](Packet& solidPacket,Packet& liquidPacket,const saltchem::MineralVector& area,std::size_t patch){
                    bool active=false;for(int j=0;j<7;j++)active=active||(area[j]>0&&rateTemplate.geometry.weights[j]>0&&
                        (rateTemplate.neutral[j]>0||rateTemplate.acid[j]>0||rateTemplate.co2[j]>0));if(!active)return;
                    const auto started=Timer::now();
                    auto solidCell=getInventory(solidPacket),liquidCell=getInventory(liquidPacket);
                    auto rates=rateTemplate;rates.geometry.inertFraction=solidPacket[22];
                    saltchem::MineralVector locked{};for(int j=0;j<7;j++)locked[j]=min(solidPacket[15+j],solidCell.mineral[j]);
                    auto& worker=localChemistry.at(model.temperatureK(),liquidPacket[23]);
                    saltchem::SurfaceFilmParameters film;
                    for(int j=0;j<7;j++){film.coverage[j]=liquidPacket[44+j];film.coatedVelocity[j]=liquidPacket[51+j]<0?INFINITY:liquidPacket[51+j];}
                    film=mineralSurface.faceFilm(face,film,patch);
                    saltchem::SurfaceReactionResult result;
                    if(precipitation){
                        rates.geometry.inertFraction=liquidPacket[22];rates.geometry.separateMineralAreas=true;
                        rates.area.fill(0);rates.geometry.lockedMoles=liquidCell.mineral;rates.geometry.surfaceExponent.fill(0);rates.geometry.curvatureLogOmega.fill(0);
                        rates.geometry.initialPorosity=saltchem::poreFraction(liquidCell,rates.geometry,model.molarVolume);
                        for(int j=0;j<7;j++){rates.geometry.growthArea[j]=area[j]/cellVolume;rates.referenceMoles[j]=std::max(liquidCell.mineral[j],1e-30);}
                        const auto adjusted=saltchem::applySurfaceFilmRates(worker.reactions,liquidCell,rates,film);
                        const auto initial=liquidCell;
                        result.stats=worker.reactions.advance(liquidCell,faceDt,rates,options);result.filmEvaluations=adjusted.filmEvaluations;
                        for(int j=0;j<7;j++)result.dissolved[j]=initial.mineral[j]-liquidCell.mineral[j];
                    }else result=saltchem::advanceSurfaceDissolution(worker.reactions,solidCell,liquidCell,rates,area,locked,model.molarVolume,cellVolume,faceDt,options,&film);
                    putInventory(solidPacket,solidCell);putInventory(liquidPacket,liquidCell);
                    for(int j=0;j<7;j++)if(!precipitation&&liquidPacket[31+j]>0)
                        liquidPacket[24+j]+=result.dissolved[j]*model.molarVolume[j]*cellVolume/(liquidPacket[31+j]*(1-epMin));
                    for(int j=0;j<7;j++){
                        const auto normal=mineralSurface.facePatchNormal(face,patch,j);const double change=result.dissolved[j]*model.molarVolume[j];
                        auto& location=precipitation?liquidPacket:solidPacket;
                        location[42]+=change*normal.x;location[43]+=change*normal.y;
                    }
                    const auto& stats=result.stats;stepAccepted+=stats.accepted;stepEvaluations+=stats.rateEvaluations+result.filmEvaluations;
                    stepRestarts+=stats.historyRestarts;stepFallbacks+=stats.integratorFallbacks;
                    stepJacobians+=stats.analyticJacobians;stepJacobianFallbacks+=stats.jacobianFallbacks;
                    auto& cost=own;
                    cost[38]+=stats.rateEvaluations+result.filmEvaluations;cost[39]+=stats.rejected;cost[40]+=stats.integratorFallbacks;cost[41]+=elapsed(started);
                };
                // Each directed graph link has one solid origin and one
                // geometrically connected liquid reservoir, even across MPI.
                for(std::size_t patch=0;patch<mineralSurface.facePatchCount(face);patch++)
                    react(own,other,mineralSurface.facePatchArea(face,patch),patch);
            });
            forAll(eps,i){
                cells[i]=getInventory(packets[i]);saltchem::MineralVector retreat{};
                for(int j=0;j<7;j++){retreat[j]=packets[i][24+j];solid[j][i]=cells[i].mineral[j]*model.molarVolume[j];}
                mineralSurface.retreatExternal(i,retreat);
                mineralSurface.accumulateFaceVector(i,{{packets[i][42],packets[i][43]}});
                cellEvaluations[i]+=packets[i][38];cellRejected[i]+=packets[i][39];cellFallbacks[i]+=packets[i][40];cellSeconds[i]+=packets[i][41];
                auto geometry=rateTemplate.geometry;geometry.inertFraction=inert[i];
                eps[i]=saltchem::poreFraction(cells[i],geometry,model.molarVolume);
            }
            mineralSurface.retireExternalSupports(cells,model.molarVolume);
        };
        auto refreshSurfaceGeometry=[&](){
            forAll(eps,i){for(int j=0;j<7;j++)solid[j][i]=cells[i].mineral[j]*model.molarVolume[j];auto geometry=rateTemplate.geometry;geometry.inertFraction=inert[i];eps[i]=saltchem::poreFraction(cells[i],geometry,model.molarVolume);}
            eps.correctBoundaryConditions();for(auto& mineral:solid)mineral.correctBoundaryConditions();
            mineralSurface.reconstruct(solid,inert,rateTemplate);mineralSurface.updateTransport(eps);
            const auto neighbours=stencil.values(eps);
            forAll(eps,i){rates[i]=rateTemplate;rates[i].geometry.inertFraction=inert[i];rates[i].geometry.neighbours=neighbours[i];mineralSurface.configure(i,cells[i],rates[i],model.molarVolume);}
        };
        if(mineralSurface.fluidSide()){
            auto advanceLocalPart=[&](double interval,bool primary){
                std::fill(errors.begin(),errors.end(),std::string{});
                std::stable_sort(workOrder.begin(),workOrder.end(),[&](label a,label b){return chemistryCost[a]>chemistryCost[b];});
                #pragma omp parallel for schedule(dynamic,1) reduction(+:stepAccepted,stepEvaluations,stepRestarts,stepFallbacks,stepJacobians,stepJacobianFallbacks,stepSkipped)
                for(label task=0;task<mesh.nCells();task++){
                    const label i=workOrder[task];const auto started=Timer::now();
                    const double evaluationsBefore=cellEvaluations[i];
                    try{
                        if(!mineralSurface.active(i)){stepSkipped++;chemistryCost[i]=1;continue;}
                        auto& worker=localChemistry.at(model.temperatureK(),absolutePressure(p[i],centres[i]));
                        auto record=[&](const saltchem::IntegrationStats& stats,int filmEvaluations=0){
                            stepAccepted+=stats.accepted;stepEvaluations+=stats.rateEvaluations+filmEvaluations;stepRestarts+=stats.historyRestarts;
                            stepFallbacks+=stats.integratorFallbacks;stepJacobians+=stats.analyticJacobians;stepJacobianFallbacks+=stats.jacobianFallbacks;
                            cellEvaluations[i]+=stats.rateEvaluations+filmEvaluations;cellRejected[i]+=stats.rejected;cellFallbacks[i]+=stats.integratorFallbacks;
                        };
                        auto parameters=rateTemplate;parameters.geometry.inertFraction=inert[i];
                        if(primary){
                            for(std::size_t patch=0;patch<mineralSurface.localPatchCount(i);patch++){
                            const auto& area=mineralSurface.localPatchArea(i,patch);
                            bool any=false;for(int j=0;j<7;j++)any=any||(area[j]>0&&parameters.geometry.weights[j]>0&&
                                (parameters.neutral[j]>0||parameters.acid[j]>0||parameters.co2[j]>0));
                            if(any){
                                auto solidCell=cells[i],liquidCell=cells[i];
                                auto locked=mineralSurface.lockedCaps(i,model.molarVolume);
                                for(int j=0;j<7;j++)locked[j]=min(locked[j],solidCell.mineral[j]);
                                const auto film=mineralSurface.localPatchFilm(i,mineralSurface.filmParameters(i,saltchem::poreFraction(cells[i],parameters.geometry,model.molarVolume)),patch);
                                const auto result=saltchem::advanceSurfaceDissolution(worker.reactions,solidCell,liquidCell,parameters,area,locked,model.molarVolume,mesh.V()[i],interval,options,&film);
                                cells[i].mineral=solidCell.mineral;cells[i].aqueous=liquidCell.aqueous;cells[i].waterEquivalent=liquidCell.waterEquivalent;
                                mineralSurface.commitPrimary(i,result,model.molarVolume,patch);record(result.stats,result.filmEvaluations);
                            }
                            }
                        }else{
                            // Every population owns its curvature and stock.
                            // Opposite cap rates must not cancel before their
                            // individual volumes and shared solution change.
                            const auto measured=mineralSurface.currentMeasure(i);
                            const auto caps=saltchem::advanceCapPopulations(worker.reactions,cells[i],mineralSurface.state[i],parameters,measured,model.molarVolume,mesh.V()[i],interval,model.temperatureK(),options);
                            mineralSurface.commitCaps(i,caps);record(caps.stats);
                            for(std::size_t patch=0;patch<mineralSurface.localPatchCount(i);patch++){
                                const auto& area=mineralSurface.localPatchArea(i,patch);bool any=false;
                                for(int j=0;j<7;j++)any=any||(area[j]>0&&parameters.geometry.weights[j]>0&&
                                    (parameters.neutral[j]>0||parameters.acid[j]>0||parameters.co2[j]>0));
                                if(!any)continue;
                                const auto initial=cells[i];auto growthParameters=parameters;
                                mineralSurface.configureResolvedGrowth(i,cells[i],growthParameters,patch);
                                const auto film=mineralSurface.localPatchFilm(i,mineralSurface.filmParameters(i,saltchem::poreFraction(cells[i],parameters.geometry,model.molarVolume)),patch);
                                const auto filmResult=saltchem::applySurfaceFilmRates(worker.reactions,cells[i],growthParameters,film);
                                const auto stats=worker.reactions.advance(cells[i],interval,growthParameters,options);
                                mineralSurface.commitResolvedGrowth(i,initial,cells[i],model.molarVolume,patch);record(stats,filmResult.filmEvaluations);
                            }
                        }
                        if(!primary)chemistryCost[i]=std::max(1.,cellEvaluations[i]-evaluationsBefore);
                    }catch(const std::exception& error){errors[i]=error.what();}
                    cellSeconds[i]+=elapsed(started);
                }
                failCollectively(errors,"LOCAL_SURFACE_REACTION");
            };
            // Every physical surface belongs to either a face reaction or a
            // local liquid-centred reaction during each half interval. Freeze
            // that allocation within the half: crossing a cell centre cannot
            // make the same interface react for a second full interval.
            advanceFaces(dt/2,false);advanceLocalPart(dt/2,true);
            mineralSurface.retireExternalSupports(cells,model.molarVolume);advanceLocalPart(dt/2,false);advanceFaces(dt/2,true);
            refreshSurfaceGeometry();
            std::fill(errors.begin(),errors.end(),std::string{});
            #pragma omp parallel for schedule(static)
            for(label i=0;i<mesh.nCells();i++)try{
                if(!mineralSurface.active(i))continue;
                auto& worker=localChemistry.at(model.temperatureK(),absolutePressure(p[i],centres[i]));
                double water;const auto aqueous=worker.speciation.conserved(cells[i],water,&properties[i],false);
                mineralSurface.seedCell(i,cells[i],aqueous.saturationIndex,dt,model.molarVolume,inert[i],model.temperatureK());
            }catch(const std::exception& error){errors[i]=error.what();}
            failCollectively(errors,"SURFACE_NUCLEATION");
            refreshSurfaceGeometry();
            advanceFaces(dt/2,false);advanceLocalPart(dt/2,true);
            mineralSurface.retireExternalSupports(cells,model.molarVolume);advanceLocalPart(dt/2,false);advanceFaces(dt/2,true);
            // Transport must use the final geometry, including any newly
            // blocked opening. This refresh performs no additional reaction.
            refreshSurfaceGeometry();
        }else{
        std::stable_sort(workOrder.begin(),workOrder.end(),[&](label a,label b){return chemistryCost[a]>chemistryCost[b];});
        #pragma omp parallel for schedule(dynamic,1) reduction(+:stepAccepted,stepEvaluations,stepRestarts,stepFallbacks,stepJacobians,stepJacobianFallbacks,stepSkipped)
        for(label task=0;task<mesh.nCells();task++){
            const label i=workOrder[task];const auto cellStart=Timer::now();
            try{
                // Exactly zero geometric area for either sign of the rate.
                // No threshold is used: small but nonzero reactions are kept.
                bool flat=true;for(double e:rates[i].geometry.neighbours)flat=flat&&e==eps[i];
                if(mineralSurface.enabled?!mineralSurface.active(i):flat){stepSkipped++;chemistryCost[i]=1;continue;}
                auto& worker=localChemistry.at(model.temperatureK(),absolutePressure(p[i],centres[i]));
                if(mineralSurface.enabled){
                    const auto firstHalf=cells[i];
                    const auto firstStats=worker.reactions.advance(cells[i],dt/2,rates[i],options);
                    mineralSurface.commit(i,firstHalf,cells[i],model.molarVolume);
                    stepAccepted+=firstStats.accepted;stepEvaluations+=firstStats.rateEvaluations;stepRestarts+=firstStats.historyRestarts;stepFallbacks+=firstStats.integratorFallbacks;stepJacobians+=firstStats.analyticJacobians;stepJacobianFallbacks+=firstStats.jacobianFallbacks;
                    cellEvaluations[i]+=firstStats.rateEvaluations;cellRejected[i]+=firstStats.rejected;cellFallbacks[i]+=firstStats.integratorFallbacks;
                    // Sample nucleation after the first chemical half-step:
                    // the initially K-free inlet cannot predict the local
                    // supersaturation created by carnallite dissolution.
                    double water;const auto aqueous=worker.speciation.conserved(cells[i],water,&properties[i],false);
                    mineralSurface.seedCell(i,cells[i],aqueous.saturationIndex,dt,model.molarVolume,inert[i],model.temperatureK());
                    mineralSurface.configure(i,cells[i],rates[i],model.molarVolume);
                }
                const auto seeded=cells[i];
                const auto st=worker.reactions.advance(cells[i],mineralSurface.enabled?dt/2:dt,rates[i],options);
                mineralSurface.commit(i,seeded,cells[i],model.molarVolume);
                stepAccepted+=st.accepted;stepEvaluations+=st.rateEvaluations;stepRestarts+=st.historyRestarts;stepFallbacks+=st.integratorFallbacks;stepJacobians+=st.analyticJacobians;stepJacobianFallbacks+=st.jacobianFallbacks;
                cellEvaluations[i]+=st.rateEvaluations;cellRejected[i]+=st.rejected;cellFallbacks[i]+=st.integratorFallbacks;chemistryCost[i]=std::max(1.,double(st.rateEvaluations));
            }catch(const std::exception& e){errors[i]=e.what();}
            cellSeconds[i]+=elapsed(cellStart);
        }
        }
        skippedChemistry+=stepSkipped;
        chemicalAccepted+=stepAccepted;chemicalEvaluations+=stepEvaluations;chemicalRestarts+=stepRestarts;chemicalFallbacks+=stepFallbacks;chemicalJacobians+=stepJacobians;jacobianFallbacks+=stepJacobianFallbacks;
#endif
        failCollectively(errors,"KINETICS");
        mineralSurface.recordReactions(before,cells,dt);
        forAll(eps,i){const auto& s=cells[i];double mass=saltchem::aqueousMass(s);totalSource[i]=0;for(int j=0;j<7;j++){const double extent=before[i].mineral[j]-s.mineral[j];solid[j][i]=s.mineral[j]*model.molarVolume[j];totalSource[i]+=saltchem::mineralMass[j]*extent/dt;source[j][i]=saltchem::componentMass[j]*(s.aqueous[j]-before[i].aqueous[j])/dt;Y[j][i]=saltchem::componentMass[j]*s.aqueous[j]/mass;if(extent>=0)dissolved[j]+=extent*mesh.V()[i];else precipitated[j]-=extent*mesh.V()[i];}if(mineralSurface.fluidSide()){totalSource[i]=(s.waterEquivalent-before[i].waterEquivalent)/dt;for(int j=0;j<7;j++)totalSource[i]+=source[j][i];}eps[i]=saltchem::poreFraction(s,rates[i].geometry,model.molarVolume);}
        // Accumulated dissolution/precipitation below is resolved per global
        // step. Local sign reversals within a step are not counted twice.
        for(int j=0;j<7;j++){solid[j].correctBoundaryConditions();Y[j].correctBoundaryConditions();}eps.correctBoundaryConditions();mineralSurface.updateTransport(eps);updateProperties(true);updateHydraulics();chemistryTime+=elapsed(tick);
        std::array<bool,7> active{};for(int j=0;j<7;j++){double boundaryMaximum=0;for(const auto& b:Y[j].boundaryField())for(double v:b)boundaryMaximum=max(boundaryMaximum,std::abs(v));reduce(boundaryMaximum,maxOp<double>());active[j]=gMax(mag(Y[j].primitiveField()))>0||gMax(mag(source[j].primitiveField()))>0||hydraulics.hasComponent(j)||boundaryMaximum>0;}
        int iteration=0;double omega=.7,windowResidual=GREAT;std::array<scalarField,7> lastDefect;bool haveDefect=false;
        saltflow::CouplingAcceleration accelerator(cfg.lookupOrDefault<label>("andersonDepth",5),cfg.lookupOrDefault<scalar>("andersonDamping",.8));
        // Include momentum/pressure in the fixed point, with fixed scales per
        // segment. This also accelerates the SIMPLE part of the coupling.
        const double velocityScale=max(gMax(mag(U.primitiveField())),1e-6),pressureScale=max(gMax(mag(pressureDeviation)),1.);
        std::array<double,7> componentScale;for(int j=0;j<7;j++)componentScale[j]=max(gMax(mag(Y[j].primitiveField())),1e-4);
        label fluxCount=massPhi.size();for(const auto& boundary:massPhi.boundaryField())fluxCount+=boundary.size();
        const double fluxScale=refDensity.value()*velocityScale;
        auto pack=[&](){scalarField x(mesh.nCells()*11+(accelerateFaceFlux?fluxCount:0));forAll(eps,i){for(int j=0;j<7;j++)x[11*i+j]=Y[j][i]/componentScale[j];for(int j=0;j<3;j++)x[11*i+7+j]=U[i][j]/velocityScale;x[11*i+10]=pressureDeviation[i]/pressureScale;}
            // The collocated momentum iteration also depends on its face
            // flux history. Include both sides of processor faces: the same
            // global Anderson coefficients preserve their opposite signs.
            if(accelerateFaceFlux){label n=mesh.nCells()*11;forAll(massPhi,f)x[n++]=massPhi[f]/(fluxScale*mesh.magSf()[f]);forAll(massPhi.boundaryField(),k)forAll(massPhi.boundaryField()[k],f)x[n++]=massPhi.boundaryField()[k][f]/(fluxScale*mesh.magSf().boundaryField()[k][f]);}return x;};
        for(;iteration<maxCouplingIterations;iteration++){
            const scalarField previousPacked(pack());const vectorField previousOuterVelocity(U.primitiveField());
            for(int corrector=0;corrector<flowCorrectors;corrector++){
                flowIteration(true,oldMass);
                // Only shorten the inner iteration. Exact EOS, continuity,
                // transport and unrelaxed momentum still decide acceptance.
                if(adaptiveFlowCorrectors&&flowChange<flowTol&&massResidual<massTolerance)break;
            }
            flowChange=max(flowChange,gMax(mag(U.primitiveField()-previousOuterVelocity))/max(gMax(mag(U.primitiveField())),1e-12));
            updateReservoirs();tick=Timer::now();scalarField previousRho(rho.primitiveField());std::array<scalarField,7> previousY;for(int j=0;j<7;j++)previousY[j]=Y[j].primitiveField();volScalarField cap(eps*rho);surfaceScalarField deff(mineralSurface.diffusivity(D,eps,rho));couplingChange=0;
            saltflow::transport(Y,active,cap,deff,massPhi,oldMass,oldY,source,dt,reuseTransport,positivityTolerance,positivityMass,transportAssemblies);
            std::array<scalarField,7> defect;double dot=0,norm=0;
            for(int j=0;j<7;j++){
                defect[j]=Y[j].primitiveField()-previousY[j];couplingChange=max(couplingChange,gMax(mag(defect[j])));
                if(haveDefect){const scalarField change(defect[j]-lastDefect[j]);dot+=sum(lastDefect[j]*change);norm+=sum(sqr(change));}
            }
            if(couplingAcceleration=="anderson"&&(flowChange>=flowTol||couplingChange>=couplingTol)){
                const scalarField proposed(accelerator.next(previousPacked,pack()));double feasible=1;
                forAll(eps,i){double sum=0,change=0;for(int j=0;j<7;j++){const double d=proposed[11*i+j]*componentScale[j]-previousY[j][i];if(d<0)feasible=min(feasible,.99*previousY[j][i]/(-d));sum+=previousY[j][i];change+=d;}if(change>0)feasible=min(feasible,.99*(1-sum)/change);}
                reduce(feasible,minOp<double>());
                if(feasible>.05){
                    forAll(eps,i){for(int j=0;j<7;j++)defect[j][i]=proposed[11*i+j]*componentScale[j]-previousY[j][i];for(int j=0;j<3;j++)U[i][j]=velocityScale*(previousPacked[11*i+7+j]+feasible*(proposed[11*i+7+j]-previousPacked[11*i+7+j]));pressureDeviation[i]=pressureScale*(previousPacked[11*i+10]+feasible*(proposed[11*i+10]-previousPacked[11*i+10]));p[i]=pressureOffset+pressureDeviation[i];}
                    if(accelerateFaceFlux){label n=mesh.nCells()*11;forAll(massPhi,f){massPhi[f]=fluxScale*mesh.magSf()[f]*(previousPacked[n]+feasible*(proposed[n]-previousPacked[n]));n++;}forAll(massPhi.boundaryField(),k)forAll(massPhi.boundaryField()[k],f){massPhi.boundaryFieldRef()[k][f]=fluxScale*mesh.magSf().boundaryField()[k][f]*(previousPacked[n]+feasible*(proposed[n]-previousPacked[n]));n++;}phi=massPhi/fvc::interpolate(rho);}
                    omega=feasible;U.correctBoundaryConditions();p.correctBoundaryConditions();
                }else omega=.7;
            }else if(couplingAcceleration=="aitken"&&haveDefect){reduce(dot,sumOp<double>());reduce(norm,sumOp<double>());if(norm>1e-40)omega=std::clamp(-omega*dot/norm,minimumCouplingRelaxation,maximumCouplingRelaxation);}
            else omega=couplingAcceleration=="anderson"?1.:.7;
            double feasible=omega;forAll(eps,i){double a=0,b=0;for(int j=0;j<7;j++){a+=previousY[j][i];b+=Y[j][i];}if(b>a)feasible=min(feasible,.5*(1-a)/(b-a));}reduce(feasible,minOp<double>());
            if(!(feasible>0))throw std::runtime_error("No positive solvent fraction during coupling");
            minimumUsedCouplingRelaxation=min(minimumUsedCouplingRelaxation,feasible);
            for(int j=0;j<7;j++){Y[j].primitiveFieldRef()=previousY[j]+feasible*defect[j];Y[j].correctBoundaryConditions();}
            lastDefect=std::move(defect);haveDefect=true;updateProperties(false,false);couplingChange=max(couplingChange,gMax(mag(rho.primitiveField()-previousRho)/rho.primitiveField()));transportTime+=elapsed(tick);
            storage.primitiveFieldRef()=(eps.primitiveField()*rho.primitiveField()-oldMass)/dt;
            massResidual=gMax(mag((fvc::div(massPhi)+storage-totalSource)().primitiveField()));
            if(traceCoupling&&iteration>=99&&iteration<=103){
                const scalarField du(mag(U.primitiveField()-previousOuterVelocity)),dr(mag(rho.primitiveField()-previousRho));
                const double maxDu=gMax(du),maxDr=gMax(dr);
                forAll(eps,i)if(du[i]>=maxDu||dr[i]>=maxDr){
                    Pout<<"COUPLING_CELL time="<<runTime.value()<<" iteration="<<iteration<<" x="<<centres[i]<<" epsilon="<<eps[i]<<" pressure="<<p[i]<<" Uold="<<previousOuterVelocity[i]<<" Unew="<<U[i]<<" rhoOld="<<previousRho[i]<<" rhoNew="<<rho[i]<<" Yold=(";
                    for(int j=0;j<7;j++)Pout<<previousY[j][i]<<' ';Pout<<") Ynew=(";for(int j=0;j<7;j++)Pout<<Y[j][i]<<' ';Pout<<')'<<endl;
                }
            }
            if(iteration%100==0)Info<<"COUPLING "<<runTime.value()<<" iteration="<<iteration<<" flow="<<flowChange<<" transport="<<couplingChange<<" mass="<<massResidual<<" relaxation="<<feasible<<endl;
            if((iteration+1)%25==0){
                const double residual=max(flowChange/flowTol,max(couplingChange/couplingTol,massResidual/massTolerance));
                if(detectCouplingStagnation&&iteration>=49&&residual>10&&residual>.9*windowResidual)
                    throw saltflow::StepConvergenceFailure("Flow/transport iteration stagnated; retry a smaller time step");
                windowResidual=residual;
            }
            if(iteration>=5&&flowChange<flowTol&&couplingChange<couplingTol){
                tick=Timer::now();const scalarField approximateRho(rho.primitiveField());updateProperties();
                couplingChange=max(couplingChange,gMax(mag(rho.primitiveField()-approximateRho)/rho.primitiveField()));
                storage.primitiveFieldRef()=(eps.primitiveField()*rho.primitiveField()-oldMass)/dt;
                massResidual=gMax(mag((fvc::div(massPhi)+storage-totalSource)().primitiveField()));
                transportTime+=elapsed(tick);
                // Transport changes density after the pressure projection.
                // Close continuity at that final density before deciding to
                // repeat all component solves. The unmodified momentum and
                // transport equations below still decide acceptance.
                if(conservativeFlowClosure&&couplingChange<couplingTol&&massResidual>=massTolerance){
                    for(int corrector=0;corrector<flowCorrectors;corrector++){
                        flowIteration(true,oldMass);for(auto& y:Y)y.correctBoundaryConditions();
                        auto propertyTick=Timer::now();updateProperties();transportTime+=elapsed(propertyTick);
                        storage.primitiveFieldRef()=(eps.primitiveField()*rho.primitiveField()-oldMass)/dt;
                        massResidual=gMax(mag((fvc::div(massPhi)+storage-totalSource)().primitiveField()));
                        if(massResidual<massTolerance)break;
                    }
                }
                if(couplingChange>=couplingTol||massResidual>=massTolerance)continue;
                tick=Timer::now();
                volScalarField finalCapacity(eps*rho);surfaceScalarField finalDiffusivity(mineralSurface.diffusivity(D,eps,rho));
                equationResidual=saltflow::residual(Y,active,finalCapacity,finalDiffusivity,massPhi,oldMass,oldY,source,dt);transportTime+=elapsed(tick);
                if(iteration%25==0)Info<<"TRANSPORT_RESIDUAL "<<runTime.value()<<" iteration="<<iteration<<" value="<<equationResidual<<endl;
                if(equationResidual<couplingTol&&checkMomentum(true)<flowTol)break;
            }
        }
        couplingIterations+=iteration+1;
        if(iteration==maxCouplingIterations)throw saltflow::StepConvergenceFailure("Flow/transport iteration did not converge");ledger(dt);if(maxBalance>conservationTolerance)throw std::runtime_error("Conservation check failed");
        Info<<"STEP "<<runTime.value()<<" iterations="<<iteration+1<<" balance="<<maxBalance<<" equation_residual="<<equationResidual<<" chemistry_s="<<chemistryTime<<" flow_s="<<flowTime<<" transport_s="<<transportTime<<endl;
    };
    const double nominalDt=restarting?restart.control[FlowCheckpoint::nominalDt]:min(timeControl.maximum,runTime.deltaTValue());
    double proposedDt=restarting?restart.control[FlowCheckpoint::proposedDt]:nominalDt;
    const double writeOrigin=restarting?restart.control[FlowCheckpoint::writeOrigin]:runTime.startTime().value();
    std::vector<double>().swap(restart.fields);restart.deviation.clear();
    auto checkpoint=[&](){
        restart.control={{runTime.value(),pressureOffset,next,maxBalance,nominalDt,proposedDt,writeOrigin,
            minAcceptedDt,maxAcceptedDt,runTime.deltaTValue(),timeControl.maximum,acceptedSteps,rejectedSteps,trialSegments,lastTimeError,maxTimeError}};
        restart.index=runTime.timeIndex();restart.initial=initialInventory;restart.boundary=boundaryIntegral;restart.balance=balance;restart.inventory=inventory();
        restart.dissolved=dissolved;restart.precipitated=precipitated;restart.deviation=pressureDeviation;
        restart.transferFields(false,p,eps,rho,inert,U,phi,massPhi,Y,solid);
        restart.write(runTime.timePath()/"flowState.bin",mesh);
        std::vector<double>().swap(restart.fields);restart.deviation.clear();
    };
    std::ofstream timeLog((out/"time_steps.csv").c_str());timeLog<<std::setprecision(16)<<"time_s,dt_s,error,accepted,reason\n";
    while(runTime.value()<runTime.endTime().value()-1e-12){
        const auto start=capture();double stop=min(next,runTime.endTime().value());stop=thermalSchedule.nextKnot(start.time,stop);stop=hydraulics.nextKnot(start.time,stop);
        const word writeControl=runTime.controlDict().lookupOrDefault<word>("writeControl","timeStep");
        const bool physicalWrite=writeControl=="runTime"||writeControl=="adjustableRunTime";double nextNativeWrite=GREAT,writeTolerance=0;
        if(physicalWrite){
            const double interval=readScalar(runTime.controlDict().lookup("writeInterval")),origin=writeOrigin;
            if(!(interval>0)||!std::isfinite(interval))throw std::invalid_argument("writeInterval must be finite and positive");
            const double index=std::floor((start.time-origin)/interval+1e-10)+1;
            nextNativeWrite=origin+index*interval;writeTolerance=interval*1e-10;stop=min(stop,nextNativeWrite);
        }
        double dt=min(proposedDt,stop-start.time);int attempts=0;
        for(;;){
            std::string failure;bool convergenceFailure=false;lastTimeError=0;
            try{
                advanceSegment(dt);
                if(timeControl.enabled){const auto coarse=capture();restore(start);advanceSegment(dt/2);advanceSegment(dt/2);lastTimeError=max(timeControl.error(coarse,eps,rho,U,Y,solid,model.molarVolume),mineralSurface.timeError(coarse.mineralSurface));}
            }catch(const saltflow::SurfaceTransportModelError&){throw;}
            catch(const saltflow::StepConvergenceFailure& e){failure=e.what();convergenceFailure=true;lastTimeError=GREAT;}
            catch(const std::exception& e){failure=e.what();lastTimeError=GREAT;}
            const bool accepted=failure.empty()&&(!timeControl.enabled||lastTimeError<=1);
            timeLog<<start.time<<','<<dt<<','<<lastTimeError<<','<<accepted<<','<<(failure.empty()?(timeControl.enabled?"estimate":"fixed"):failure)<<std::endl;
            if(accepted){acceptedSteps++;minAcceptedDt=min(minAcceptedDt,dt);maxAcceptedDt=max(maxAcceptedDt,dt);maxTimeError=max(maxTimeError,lastTimeError);proposedDt=timeControl.enabled?min(timeControl.maximum,dt*timeControl.factor(lastTimeError)):min(nominalDt,1.25*(attempts?dt:proposedDt));break;}
            rejectedSteps++;restore(start);Info<<"REJECT time="<<start.time<<" dt="<<dt<<" error="<<lastTimeError<<" reason="<<failure<<endl;
            if((!timeControl.enabled&&!(timeControl.recoverFailedSteps&&convergenceFailure))||++attempts>timeControl.maximumRejections||dt<=timeControl.minimum*(1+1e-10))throw std::runtime_error("Time step could not satisfy accuracy/convergence bounds: "+failure);
            dt=max(timeControl.minimum,dt*(convergenceFailure?.5:min(.5,timeControl.factor(lastTimeError))));
        }
        // Trial segments move physical time only. Commit the native clock
        // once per accepted global step, so retries cannot consume output
        // events or leave deltaT/writeTime state from a rejected coarse step.
        static_cast<TimeState&>(runTime)=start.nativeTime;
        runTime.setTime(start.time,start.index);runTime.setDeltaT(dt,false);++runTime;
        checkPressureApproximation();
        Info<<"ACCEPT time="<<runTime.value()<<" dt="<<dt<<" error="<<lastTimeError<<" balance="<<maxBalance<<endl;
        if(runTime.value()+1e-10>=next||runTime.value()+1e-10>=runTime.endTime().value()){loadCells();save();while(next<=runTime.value()+1e-10)next+=frame;}
        // Native runTime output rounds to the closest step (t + dt/2).
        // We already clip to exact write events, so suppress early rounding
        // and force writing at the actual boundary instead.
        const bool writeDue=physicalWrite?runTime.value()+writeTolerance>=nextNativeWrite:runTime.writeTime();
        auto writeStart=Timer::now();if(writeDue){chemicalPressure.correctBoundaryConditions();runTime.writeNow();mineralSurface.checkpoint(true);checkpoint();}outputTime+=elapsed(writeStart);
    }
    // Compact per-cell profiling survives export and can drive a weighted
    // decomposition for a subsequent run. Counts include rejected trial work.
    {labelList globalCells(mesh.nCells());std::iota(globalCells.begin(),globalCells.end(),0);
        if(Pstream::parRun()){labelIOList addressing(IOobject("cellProcAddressing",mesh.facesInstance(),polyMesh::meshSubDir,mesh,IOobject::MUST_READ,IOobject::NO_WRITE));globalCells=addressing;}
        std::ofstream f((out/"chemistry_cells.csv").c_str());f<<std::setprecision(16)<<"global_cell,x_m,y_m,z_m,wall_s,rate_evaluations,rejected_substeps,fallbacks\n";
        forAll(eps,i)f<<globalCells[i]<<','<<centres[i].x()<<','<<centres[i].y()<<','<<centres[i].z()<<','<<cellSeconds[i]<<','<<cellEvaluations[i]<<','<<cellRejected[i]<<','<<cellFallbacks[i]<<'\n';}
    volScalarField chemistryWork(IOobject("chemistryWork",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),mesh,dimensionedScalar(dimless,1));
    forAll(eps,i)chemistryWork[i]=1+cellEvaluations[i]/max(trialSegments,1.);chemistryWork.correctBoundaryConditions();runTime.writeNow();mineralSurface.checkpoint(true);checkpoint();
    double rankWork=std::accumulate(cellSeconds.begin(),cellSeconds.end(),0.),maximumRankWork=rankWork,totalRankWork=rankWork;
    reduce(maximumRankWork,maxOp<double>());reduce(totalRankWork,sumOp<double>());reduce(skippedChemistry,sumOp<double>());reduce(eosLinearized,sumOp<double>());reduce(boundaryPropertyEvaluations,sumOp<double>());reduce(boundaryPropertyHits,sumOp<double>());
    const double chemistryImbalance=totalRankWork>0?maximumRankWork*Pstream::nProcs()/totalRankWork:1.;
    double wall=elapsed(wallStart);reduce(wall,maxOp<double>());for(double* v:{&chemistryTime,&flowTime,&transportTime,&outputTime,&initializationTime,&propertyTime,&interfaceGeometryTime})reduce(*v,maxOp<double>());
    for(double* v:{&chemicalAccepted,&chemicalEvaluations,&chemicalRestarts,&chemicalFallbacks,&chemicalJacobians,&jacobianFallbacks})reduce(*v,sumOp<double>());
    reduce(transportAssemblies,maxOp<double>());reduce(positivityMass,maxOp<double>());
    double effectivePressure=0,pressureRaisedCells=0;forAll(eps,i){const auto& local=localChemistry.at(model.temperatureK(),chemicalPressure[i]);const auto w=local.model.waterState(properties[i].waterActivity);effectivePressure=max(effectivePressure,w.effectivePressureAtm);pressureRaisedCells+=w.pressureRaised;}
    reduce(effectivePressure,maxOp<double>());reduce(pressureRaisedCells,sumOp<double>());
    const double minimumPorosity=gMin(eps.primitiveField());
    const double minimumChemicalPressure=gMin(chemicalPressure.primitiveField()),maximumChemicalPressure=gMax(chemicalPressure.primitiveField());
    if(Pstream::master()){std::ofstream f((out/"summary.csv").c_str());f<<std::setprecision(16)<<"metric,value\nphysical_time_s,"<<runTime.value()<<"\nsolver_wall_s,"<<wall<<"\nchemistry_s,"<<chemistryTime<<"\nflow_s,"<<flowTime<<"\ntransport_s,"<<transportTime<<"\noutput_s,"<<outputTime<<"\nreference_initialization_s,"<<initializationTime<<"\nmax_component_balance_error,"<<maxBalance<<"\nnative_accepted_substeps,"<<chemicalAccepted<<"\nnative_rate_evaluations,"<<chemicalEvaluations<<"\nnative_history_restarts,"<<chemicalRestarts<<"\nnative_integrator_fallbacks,"<<chemicalFallbacks<<"\nanalytic_jacobians,"<<chemicalJacobians<<"\njacobian_fallbacks,"<<jacobianFallbacks<<"\ntransport_matrix_assemblies_max_rank,"<<transportAssemblies<<"\npositivity_redistributed_kg,"<<positivityMass<<"\nfinal_transport_equation_residual,"<<equationResidual<<"\nwater_temperature_K,"<<model.waterState(1).waterTemperatureK<<"\nmax_effective_pressure_atm,"<<effectivePressure<<"\npressure_raised_cells,"<<pressureRaisedCells<<"\nepsilon_min_configured,"<<epMin<<"\nminimum_final_porosity,"<<minimumPorosity<<"\npore_closure_width,"<<rateTemplate.geometry.closureWidth<<"\nphase_direction_width,"<<rateTemplate.geometry.directionWidth<<'\n';}
    if(Pstream::master()){
        std::ofstream f((out/"summary.csv").c_str(),std::ios::app);f<<std::setprecision(16)
            <<"min_chemical_pressure_atm,"<<minimumChemicalPressure<<"\nmax_chemical_pressure_atm,"<<maximumChemicalPressure<<"\nhydraulic_reference_density_kg_m3,"<<refDensity.value()
            <<"\nchemistry_pressure_feedback,"<<flowPressure<<"\nadaptive_time_step,"<<timeControl.enabled
            <<"\naccepted_global_steps,"<<acceptedSteps<<"\nrejected_global_steps,"<<rejectedSteps<<"\ntrial_segments,"<<trialSegments
            <<"\nmax_accepted_time_error,"<<maxTimeError<<"\nminimum_accepted_dt_s,"<<(acceptedSteps?minAcceptedDt:0)<<"\nmaximum_accepted_dt_s,"<<maxAcceptedDt
            <<"\ncoupling_iterations,"<<couplingIterations<<"\nzero_area_chemistry_skips,"<<skippedChemistry<<"\nchemistry_rank_imbalance,"<<chemistryImbalance<<'\n';
        f<<"anchor_transition_width,"<<anchorWidth<<'\n';
        f<<"eos_linearized_evaluations,"<<eosLinearized<<'\n';
        f<<"density_pressure_local,"<<(densityPressureMode=="local")<<"\nmaximum_density_pressure_approximation_error,"<<maxDensityPressureError<<"\ndensity_pressure_tolerance,"<<densityPressureTolerance<<'\n';
        f<<"flow_tolerance,"<<flowTol<<"\ncoupling_tolerance,"<<couplingTol<<"\nconservation_tolerance,"<<conservationTolerance<<"\nflow_correctors,"<<flowCorrectors<<'\n';
        f<<"final_momentum_equation_residual,"<<momentumResidual<<'\n';
        f<<"adaptive_flow_correctors,"<<adaptiveFlowCorrectors<<"\nflow_solves,"<<flowSolves<<"\nproperties_s_inclusive,"<<propertyTime<<'\n';
        f<<"conservative_flow_closure,"<<conservativeFlowClosure<<'\n';
        f<<"recover_failed_time_steps,"<<timeControl.recoverFailedSteps<<"\nmax_coupling_iterations,"<<maxCouplingIterations<<'\n';
        f<<"consistent_momentum,"<<consistentMomentum<<"\ndetect_coupling_stagnation,"<<detectCouplingStagnation<<'\n';
        f<<"linear_face_momentum,"<<(faceMomentumInterpolation=="linear")<<"\nmaximum_reservoir_flux_error,"<<maxReservoirFluxError<<'\n';
        f<<"minimum_coupling_relaxation,"<<minimumCouplingRelaxation<<"\nminimum_used_coupling_relaxation,"<<minimumUsedCouplingRelaxation<<'\n';
        f<<"accelerate_face_flux,"<<accelerateFaceFlux<<'\n';
        f<<"heterogeneous_interface_model,"<<mineralSurface.enabled<<"\ninterface_geometry_s_inclusive,"<<interfaceGeometryTime<<'\n';
        f<<"solid_owned_sharp_surface_peak_cells,"<<mineralSurface.sharpSurfaceCells()
         <<"\nsolid_owned_sharp_surface_peak_area_m2,"<<mineralSurface.sharpSurfaceArea()<<'\n';
        f<<"cache_boundary_properties,"<<cacheBoundaryProperties<<"\nboundary_property_evaluations,"<<boundaryPropertyEvaluations<<"\nboundary_property_cache_hits,"<<boundaryPropertyHits<<'\n';
    }
    Info<<"COMPLETE wall_s="<<wall<<" chemistry_s="<<chemistryTime<<" flow_s="<<flowTime<<" transport_s="<<transportTime<<" balance="<<maxBalance<<endl;return 0;
}
