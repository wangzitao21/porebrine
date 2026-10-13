// SPDX-License-Identifier: GPL-3.0-or-later
#include "ComponentTransport.H"
#include "CouplingAcceleration.H"
#include "DiffusiveConductance.H"
#include <fstream>
#include <iomanip>

int main(int argc,char** argv){
    using namespace Foam;
    argList args(argc,argv);if(!args.checkRootCase())return 1;
    Time runTime(Time::controlDictName,args);
    fvMesh mesh(IOobject(fvMesh::defaultRegion,runTime.timeName(),runTime,IOobject::MUST_READ));
    IOdictionary config(IOobject("transportTest",runTime.constant(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE));
    const word mode(config.lookup("mode"));
    const scalar speed=config.lookupOrDefault<scalar>("speed",0),D=config.lookupOrDefault<scalar>("diffusivity",.1);
    const scalar rightD=config.lookupOrDefault<scalar>("rightDiffusivity",D),reaction=config.lookupOrDefault<scalar>("reaction",.7);
    const scalar dt=readScalar(config.lookup("dt"));const label steps=readLabel(config.lookup("steps"));
    const bool reuse=config.lookupOrDefault<bool>("reuse",true);
    const scalar pi=constant::mathematical::pi;
    const bool stationary=mode!="decay"&&!(mode=="barrier"&&rightD==0);
    // Independent series-resistance values also exercise unequal distances
    // and coefficients whose reciprocals would overflow in double precision.
    for(const auto& pair:std::initializer_list<std::array<double,3>>{
        {1,1,.5},{1,1e-12,.25},{1e-310,1e-300,.75},{1e-310,1e300,.5},
        {1e300,1e300,.4},{0,1,.5},{2,3,0},{2,3,1}}){
        const long double a=pair[0],b=pair[1],w=pair[2];
        const double expected=a==0||b==0?0:double(1/((1-w)/a+w/b));
        const double actual=saltflow::harmonicConductance(pair[0],pair[1],pair[2]);
        if(!std::isfinite(actual)||std::abs(actual-expected)>1e-13*std::abs(expected))return 2;
    }
    bool invalidRejected=false;
    try{saltflow::harmonicConductance(-1,1);}catch(const std::invalid_argument&){invalidRejected=true;}
    if(!invalidRejected)return 2;
    PtrList<volScalarField> Y(7),source(7);
    Y.set(0,new volScalarField(IOobject("Y0",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE),mesh));
    for(int j=1;j<7;j++)Y.set(j,new volScalarField(IOobject("Y"+name(j),runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),Y[0]));
    const auto exact=[&](scalar x,scalar time){
        if(mode=="layer"){
            const scalar flux=1/(.5/D+.5/rightD+1/reaction);
            return 1-flux*(x<.5?x/D:.5/D+(x-.5)/rightD);
        }
        if(mode=="barrier"){
            if(rightD==0)return x<1./3?1.:0.;
            const scalar resistance=x<1./3?x/D:x<2./3?1/(3*D)+(x-1./3)/rightD:
                1/(3*D)+1/(3*rightD)+(x-2./3)/D;
            return 1-resistance/(2/(3*D)+1/(3*rightD));
        }
        return 1+.2*std::sin(pi*x)*(mode=="decay"?std::exp(-D*pi*pi*time):1);
    };
    for(int j=0;j<7;j++){
        forAll(Y[j],i)Y[j][i]=(j+1)*(stationary?1:exact(mesh.C()[i].x(),0));
        forAll(Y[j].boundaryField(),patch){
            auto& boundary=Y[j].boundaryFieldRef()[patch];
            if(boundary.type()=="fixedValue")boundary==scalar(mode=="barrier"&&mesh.boundary()[patch].name()=="outlet"?0:j+1);
            if(boundary.type()=="mixed"){
                auto& mixed=refCast<mixedFvPatchScalarField>(boundary);
                mixed.refValue()=0;mixed.refGrad()=0;
                mixed.valueFraction()=reaction/(reaction+rightD*mesh.boundary()[patch].deltaCoeffs());
            }
        }
        Y[j].correctBoundaryConditions();
        source.set(j,new volScalarField(IOobject("source"+name(j),runTime.timeName(),mesh),mesh,dimensionedScalar(dimDensity/dimTime,0),"zeroGradient"));
    }
    volScalarField capacity(IOobject("capacity",runTime.timeName(),mesh),mesh,dimensionedScalar(dimDensity,stationary?0:1),"zeroGradient");
    volScalarField cellDiffusivity(IOobject("cellDiffusivity",runTime.timeName(),mesh),mesh,dimensionedScalar(dimDensity*dimArea/dimTime,D),"zeroGradient");
    forAll(cellDiffusivity,i){
        const scalar x=mesh.C()[i].x();
        cellDiffusivity[i]=(mode=="layer"&&x>=.5)||(mode=="barrier"&&x>=1./3&&x<2./3)?rightD:D;
    }
    cellDiffusivity.correctBoundaryConditions();
    surfaceScalarField diffusion(IOobject("diffusion",runTime.timeName(),mesh),saltflow::diffusiveConductance(cellDiffusivity));
    surfaceScalarField flux(IOobject("massPhi",runTime.timeName(),mesh),mesh,dimensionedScalar(dimMass/dimTime,0));
    forAll(flux,f)flux[f]=speed*mesh.Sf()[f].x();
    forAll(flux.boundaryField(),patch)forAll(flux.boundaryField()[patch],f)flux.boundaryFieldRef()[patch][f]=speed*mesh.Sf().boundaryField()[patch][f].x();
    const std::array<bool,7> active{{true,true,true,true,true,true,true}};
    const scalarField oldMass(capacity.primitiveField());
    double corrected=0,assemblies=0,worstResidual=0;label nonlinearIterations=0;
    const bool linearOperator=saltflow::reusable(Y,flux,diffusion),canReuse=reuse&&linearOperator;
    const auto pack=[&](){
        scalarField values(7*mesh.nCells());
        forAll(Y[0],i)for(int j=0;j<7;j++)values[7*i+j]=Y[j][i]/(j+1);
        return values;
    };
    for(label step=0;step<steps;step++){
        std::array<scalarField,7> oldY;
        for(int j=0;j<7;j++){
            oldY[j]=Y[j].primitiveField();
            forAll(Y[j],i){
                const scalar x=mesh.C()[i].x();
                source[j][i]=mode=="steady"?(j+1)*.2*(speed*pi*std::cos(pi*x)+D*pi*pi*std::sin(pi*x)):0;
            }
        }
        double finalResidual=GREAT;
        // The limiter depends on the current component field. Keep the time
        // level fixed while converging its coefficients with the same
        // safeguarded Anderson implementation used by the production loop.
        saltflow::CouplingAcceleration accelerator;
        for(label iteration=0;iteration<(linearOperator?1:2000);iteration++){
            const scalarField previous(linearOperator?scalarField():pack());
            saltflow::transport(Y,active,capacity,diffusion,flux,oldMass,oldY,source,dt,reuse,1e-12,corrected,assemblies);
            if(!linearOperator){
                const scalarField proposed(accelerator.next(previous,pack()));double feasible=1;
                forAll(proposed,i)if(proposed[i]<0)feasible=min(feasible,.99*previous[i]/(previous[i]-proposed[i]));
                reduce(feasible,minOp<double>());
                forAll(Y[0],i)for(int j=0;j<7;j++)Y[j][i]=(j+1)*(previous[7*i+j]+feasible*(proposed[7*i+j]-previous[7*i+j]));
            }
            for(int j=0;j<7;j++)Y[j].correctBoundaryConditions();
            nonlinearIterations++;
            finalResidual=saltflow::residual(Y,active,capacity,diffusion,flux,oldMass,oldY,source,dt);
            if(finalResidual<1e-11)break;
        }
        worstResidual=max(worstResidual,finalResidual);
    }
    double squareError=0,volume=0,maxError=0,minimum=GREAT,maxComponentDifference=0;
    forAll(Y[0],i){
        const double reference=exact(mesh.C()[i].x(),steps*dt),error=Y[0][i]-reference;
        squareError+=error*error*mesh.V()[i];volume+=mesh.V()[i];maxError=max(maxError,std::abs(error));minimum=min(minimum,Y[0][i]);
        for(int j=1;j<7;j++)maxComponentDifference=max(maxComponentDifference,std::abs(Y[j][i]/(j+1)-Y[0][i]));
    }
    reduce(squareError,sumOp<double>());reduce(volume,sumOp<double>());reduce(maxError,maxOp<double>());reduce(minimum,minOp<double>());reduce(maxComponentDifference,maxOp<double>());
    double leftFlux=0,rightFlux=0;
    const surfaceScalarField diffusiveFlux(-diffusion*fvc::snGrad(Y[0])*mesh.magSf());
    forAll(mesh.boundary(),patch){
        if(mesh.boundary()[patch].name()=="inlet")leftFlux=-gSum(diffusiveFlux.boundaryField()[patch])/gSum(mesh.magSf().boundaryField()[patch]);
        if(mesh.boundary()[patch].name()=="outlet")rightFlux=gSum(diffusiveFlux.boundaryField()[patch])/gSum(mesh.magSf().boundaryField()[patch]);
    }
    const scalar l2=std::sqrt(squareError/volume);
    std::ofstream cells((runTime.path()/"transport_cells.csv").c_str());
    cells<<std::setprecision(17)<<"x,y,Y0,Y1,Y2,Y3,Y4,Y5,Y6\n";
    forAll(Y[0],i){
        cells<<mesh.C()[i].x()<<','<<mesh.C()[i].y();
        for(int j=0;j<7;j++)cells<<','<<Y[j][i];cells<<'\n';
    }
    if(Pstream::master()){
        std::ofstream out((runTime.path()/"transport_metrics.csv").c_str());
        out<<std::setprecision(17)<<"metric,value\nl2_error,"<<l2<<"\nmax_error,"<<maxError<<"\nminimum,"<<minimum
           <<"\ncomponent_scaling_error,"<<maxComponentDifference<<"\nresidual,"<<worstResidual<<"\ncorrected_mass,"<<corrected
           <<"\nassemblies,"<<assemblies<<"\nnonlinear_iterations,"<<nonlinearIterations<<"\nmatrix_reused,"<<canReuse
           <<"\nleft_flux,"<<leftFlux<<"\nright_flux,"<<rightFlux<<'\n';
        Info<<"TRANSPORT_ANALYTIC_CHECK mode="<<mode<<" speed="<<speed<<" L2="<<l2<<" residual="<<worstResidual<<endl;
    }
    return minimum<0||!std::isfinite(l2)||worstResidual>1e-9||maxComponentDifference>1e-9;
}
