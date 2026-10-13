// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "fvCFD.H"
#include "mixedFvPatchFields.H"
namespace saltflow {
using namespace Foam;
inline bool mixedReservoir(const fvPatchScalarField& b){return b.type()=="inletOutlet"||b.type()=="mixed";}
inline bool schemeIs(ITstream& stream,std::initializer_list<const char*> expected){
    stream.rewind();bool same=stream.size()==label(expected.size());
    for(const char* s:expected){if(!same)break;word v;stream>>v;same=v==s;}stream.rewind();return same;
}
inline bool reusable(const PtrList<volScalarField>& Y,const surfaceScalarField& flux,const surfaceScalarField& diffusivity){
    const auto& mesh=Y[0].mesh();
    forAll(Y,j){
        if(!schemeIs(mesh.divScheme("div("+flux.name()+','+Y[j].name()+')'),{"Gauss","upwind"})||
           !schemeIs(mesh.laplacianScheme("laplacian("+diffusivity.name()+','+Y[j].name()+')'),{"Gauss","linear","orthogonal"}))return false;
        forAll(Y[j].boundaryField(),patch){
            const auto& b=Y[j].boundaryField()[patch];
            if(b.type()!=Y[0].boundaryField()[patch].type())return false;
            if(b.type()!="fixedValue"&&b.type()!="zeroGradient"&&b.type()!="empty"&&b.type()!="processor"&&b.type()!="cyclic"&&!mixedReservoir(b))return false;
            if(mixedReservoir(b)){
                const auto& mixed=refCast<const mixedFvPatchScalarField>(b);
                const auto& first=refCast<const mixedFvPatchScalarField>(Y[0].boundaryField()[patch]);
                if(mixed.valueFraction()!=first.valueFraction()||sum(mag(mixed.refGrad()))!=0)return false;
            }
        }
    }return true;
}
// Redistribute only round-off-sized negative mass, with one global inventory
// constraint per component. A subsequent equation residual check is mandatory.
inline double positiveConservative(volScalarField& Y,const volScalarField& capacity,double tolerance){
    const double lowest=gMin(Y.primitiveField());if(lowest>=0)return 0;
    if(lowest < -tolerance)throw std::runtime_error("Transport positivity failure in "+std::string(Y.name()));
    double negative=0,positive=0;const auto& V=Y.mesh().V();
    forAll(Y,i){const double amount=capacity[i]*V[i]*Y[i];if(amount<0)negative-=amount;else positive+=amount;}
    reduce(negative,sumOp<double>());reduce(positive,sumOp<double>());
    if(negative>positive||positive<=0)throw std::runtime_error("Negative total component inventory in "+std::string(Y.name()));
    const double factor=1-negative/positive;
    forAll(Y,i)Y[i]=Y[i]<0?0:Y[i]*factor;
    return negative;
}
inline void transport(PtrList<volScalarField>& Y,const std::array<bool,7>& active,
    const volScalarField& capacity,const surfaceScalarField& diffusivity,const surfaceScalarField& flux,
    const scalarField& oldMass,const std::array<scalarField,7>& oldY,const PtrList<volScalarField>& source,
    double dt,bool reuse,double positivityTolerance,double& correctedMass,double& assemblies){
    const fvMesh& mesh=Y[0].mesh();const dimensionedScalar dtime("dt",dimTime,dt);
    reuse=reuse&&reusable(Y,flux,diffusivity);
    if(!reuse){
        for(int j=0;j<7;j++)if(active[j]){
            fvScalarMatrix equation(fvm::Sp(capacity/dtime,Y[j])+fvm::div(flux,Y[j])-fvm::laplacian(diffusivity,Y[j]));
            equation.source()+=(oldMass*oldY[j]/dt+source[j].primitiveField())*mesh.V();equation.solve();assemblies++;
            correctedMass+=positiveConservative(Y[j],capacity,positivityTolerance);
        }return;
    }
    volScalarField unknown(IOobject("Ytransport",mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE,false),Y[0]);
    forAll(unknown.boundaryField(),patch){
        auto& b=unknown.boundaryFieldRef()[patch];
        if(b.type()=="fixedValue")b==1;
        if(mixedReservoir(b)){refCast<mixedFvPatchScalarField>(b).refValue()=1;b.evaluate();}
    }
    fvScalarMatrix equation(fvm::Sp(capacity/dtime,unknown)+fvm::div(flux,unknown,"div("+flux.name()+','+Y[0].name()+')')-fvm::laplacian(diffusivity,unknown,"laplacian("+diffusivity.name()+','+Y[0].name()+')'));assemblies++;
    const auto unitBoundary=equation.boundaryCoeffs();const scalarField explicitSource(equation.source());
    auto solver=equation.solver(mesh.solver(Y[0].name()));
    for(int j=0;j<7;j++)if(active[j]){
        unknown.primitiveFieldRef()=Y[j].primitiveField();
        forAll(unknown.boundaryField(),patch)if(unknown.boundaryField()[patch].type()=="fixedValue"||mixedReservoir(unknown.boundaryField()[patch])){
            unknown.boundaryFieldRef()[patch]==Y[j].boundaryField()[patch];
            equation.boundaryCoeffs()[patch]=unitBoundary[patch]*(mixedReservoir(Y[j].boundaryField()[patch])?
                refCast<const mixedFvPatchScalarField>(Y[j].boundaryField()[patch]).refValue():scalarField(Y[j].boundaryField()[patch]));
        }
        equation.source()=explicitSource+(oldMass*oldY[j]/dt+source[j].primitiveField())*mesh.V();
        solver->solve(mesh.solver(Y[j].name()));Y[j].primitiveFieldRef()=unknown.primitiveField();
        correctedMass+=positiveConservative(Y[j],capacity,positivityTolerance);
    }
}
// Check the actual, unrelaxed finite-volume equations at the final EOS state.
inline double residual(const PtrList<volScalarField>& Y,const std::array<bool,7>& active,
    const volScalarField& capacity,const surfaceScalarField& diffusivity,const surfaceScalarField& flux,
    const scalarField& oldMass,const std::array<scalarField,7>& oldY,const PtrList<volScalarField>& source,double dt){
    const dimensionedScalar dtime("dt",dimTime,dt);double maximum=0;
    for(int j=0;j<7;j++)if(active[j]){
        fvScalarMatrix equation(fvm::Sp(capacity/dtime,Y[j])+fvm::div(flux,Y[j])-fvm::laplacian(diffusivity,Y[j]));
        equation.source()+=(oldMass*oldY[j]/dt+source[j].primitiveField())*Y[j].mesh().V();
        scalarField diagonal(equation.diag()),r(equation.source()-equation.diag()*Y[j].primitiveField());
        const auto& lower=equation.lduAddr().lowerAddr();const auto& upper=equation.lduAddr().upperAddr();
        forAll(lower,f){r[lower[f]]-=equation.upper()[f]*Y[j][upper[f]];r[upper[f]]-=equation.lower()[f]*Y[j][lower[f]];}
        // Each coupled-face contribution is included exactly once, just as
        // in the linear solve. Do not combine an LDU interface residual with
        // an additional fvMatrix coupled boundary source.
        forAll(Y[j].boundaryField(),patch){
            const auto& b=Y[j].boundaryField()[patch];const auto& cells=Y[j].mesh().boundary()[patch].faceCells();
            const auto& ic=equation.internalCoeffs()[patch];const auto& bc=equation.boundaryCoeffs()[patch];
            scalarField other;if(b.coupled())other=b.patchNeighbourField();
            forAll(cells,f){const label c=cells[f];diagonal[c]+=ic[f];r[c]-=ic[f]*Y[j][c];r[c]+=bc[f]*(b.coupled()?other[f]:1);}
        }
        maximum=max(maximum,gMax(mag(r)/mag(diagonal)));
    }return maximum;
}
}
