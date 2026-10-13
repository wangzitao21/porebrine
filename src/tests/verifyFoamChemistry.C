// SPDX-License-Identifier: GPL-3.0-or-later
#include "FoamSaltChemistry.H"
#include "ComponentTransport.H"
#include <cmath>
#include <algorithm>
int main(int argc,char** argv){
    using namespace Foam;argList::addOption("database","path","Pitzer database path");argList args(argc,argv);if(!args.checkRootCase())return 1;
    Time runTime(Time::controlDictName,args);fvMesh mesh(IOobject(fvMesh::defaultRegion,runTime.timeName(),runTime,IOobject::MUST_READ));
    const fileName db=args.getOrDefault<fileName>("database","src/chemistry/pitzer.dat");SaltChemistryHandle* handle=nullptr;
    if(salt_chemistry_create(db.c_str(),298.15,1,&handle))FatalErrorInFunction<<salt_chemistry_last_error()<<exit(FatalError);
    double totals[7]={4,.2,.1,.5,5.328,.035,.002},minerals[7]={1,1,1,1,1,1,1};SaltCell initial;SaltResult initResult;
    if(salt_chemistry_initialize(handle,totals,1,minerals,&initial,&initResult))FatalErrorInFunction<<salt_chemistry_last_error()<<exit(FatalError);
    const double waterPerBulk=400;const auto makeField=[&](const word& n,const dimensionSet& dim,double value){return new volScalarField(IOobject(n,runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedScalar(dim,value),"zeroGradient");};
    PtrList<volScalarField> aqueous(7),solid(7),area(7),reference(7),sources(7);
    for(int j=0;j<7;j++){aqueous.set(j,makeField("testAqueous"+name(j),dimMoles/dimVolume,waterPerBulk*initial.aqueous_mol[j]));solid.set(j,makeField("testMineral"+name(j),dimMoles/dimVolume,waterPerBulk*initial.mineral_mol[j]));area.set(j,makeField("testArea"+name(j),dimless/dimLength,0));reference.set(j,makeField("testReference"+name(j),dimMoles/dimVolume,waterPerBulk));sources.set(j,makeField("testSource"+name(j),dimMoles/dimVolume/dimTime,0));}
    auto waterPtr=makeField("testWater",dimDensity,waterPerBulk*initial.water_equivalent_kg),rhoPtr=makeField("testRho",dimDensity,0),sourcePtr=makeField("testMassSource",dimDensity/dimTime,0);
    SaltRates parameters{};double coeff[7]={.02,.002,.02,1e-5,1e-4,.02,1e-5};for(int j=0;j<7;j++){parameters.k_neutral[j]=coeff[j];parameters.acid_order[j]=parameters.dissolution_exponent[j]=parameters.affinity_exponent[j]=1;for(label i=0;i<std::min<label>(mesh.nCells(),16);i++)area[j][i]=waterPerBulk;}
    // Fully dry solid cells must remain unchanged and receive zero fluid sources.
    for(label i=16;i<mesh.nCells();i+=2){(*waterPtr)[i]=0;(*rhoPtr)[i]=1234;for(int j=0;j<7;j++)aqueous[j][i]=0;}
    SaltOptions options{2e-7,1e-12,10000,Pstream::parRun()?2:8};advanceSaltChemistry(mesh,handle,aqueous,solid,*waterPtr,area,reference,parameters,.5,options,*rhoPtr,*sourcePtr,sources);
    SaltCell direct=initial;SaltResult directResult;int status=0;for(int j=0;j<7;j++)parameters.area_m2[j]=parameters.reference_mol[j]=1;options.threads=1;
    if(salt_chemistry_advance(handle,1,&direct,&parameters,.5,&options,&directResult,&status))FatalErrorInFunction<<salt_chemistry_last_error()<<exit(FatalError);
    auto testFraction=makeField("testFraction",dimless,1e-4),testCapacity=makeField("testCapacity",dimDensity,700);
    if(mesh.nCells())(*testFraction)[0]=-1e-15;
    const double before=gSum(testFraction->primitiveField()*testCapacity->primitiveField()*mesh.V());
    const double redistributed=saltflow::positiveConservative(*testFraction,*testCapacity,1e-12);
    const double after=gSum(testFraction->primitiveField()*testCapacity->primitiveField()*mesh.V());
    if(gMin(testFraction->primitiveField())<0||std::abs(after-before)>1e-14*std::abs(before)||redistributed<=0)FatalErrorInFunction<<"Nonconservative positivity repair"<<exit(FatalError);
    delete testFraction;delete testCapacity;
    double error=0;for(label i=0;i<std::min<label>(mesh.nCells(),16);i++){error=max(error,std::abs((*rhoPtr)[i]-directResult.density_kg_m3));error=max(error,std::abs((*sourcePtr)[i]/waterPerBulk-directResult.aqueous_mass_source_kg_s));for(int j=0;j<7;j++)error=max(error,std::abs(aqueous[j][i]/waterPerBulk-direct.aqueous_mol[j]));}
    for(label i=16;i<mesh.nCells();i++){bool dry=i%2==0;error=max(error,std::abs((*rhoPtr)[i]-(dry?1234:initResult.density_kg_m3)));error=max(error,std::abs((*sourcePtr)[i]));if(dry){error=max(error,std::abs((*waterPtr)[i]));for(int j=0;j<7;j++){error=max(error,std::abs(aqueous[j][i]));error=max(error,std::abs(solid[j][i]/waterPerBulk-initial.mineral_mol[j]));error=max(error,std::abs(sources[j][i]));}}}
    reduce(error,maxOp<double>());salt_chemistry_destroy(handle);delete waterPtr;delete rhoPtr;delete sourcePtr;
    label totalCells=mesh.nCells(),reactiveCells=std::min<label>(mesh.nCells(),16);reduce(totalCells,sumOp<label>());reduce(reactiveCells,sumOp<label>());
    Info<<"FOAM_CHEMISTRY_CHECK cells="<<totalCells<<" reactive_cells="<<reactiveCells<<" ranks="<<Pstream::nProcs()<<" max_error="<<error<<" density_and_sources="<<(error<1e-8?"PASS":"FAIL")<<endl;return error<1e-8?0:1;
}
