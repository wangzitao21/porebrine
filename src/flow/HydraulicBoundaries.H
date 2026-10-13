// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "fvCFD.H"
#include "Tuple2.H"
#include "inletOutletFvPatchFields.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "pressureInletOutletVelocityFvPatchVectorField.H"
#include "BrineInlet.hpp"

namespace saltflow {
using namespace Foam;
struct HydraulicPort {
    label patch;
    word type,chemistry,concentrationBoundary,velocityBoundary;
    double value=0,minimumPorosity=0;
    List<Tuple2<scalar,scalar>> schedule;
    std::array<double,7> fractions{};
    double at(double t)const {
        if(schedule.empty())return value;
        if(t<=schedule.first().first())return schedule.first().second();
        for(label i=1;i<schedule.size();i++)if(t<=schedule[i].first()){
            const double w=(t-schedule[i-1].first())/(schedule[i].first()-schedule[i-1].first());
            return (1-w)*schedule[i-1].second()+w*schedule[i].second();
        }
        return schedule.last().second();
    }
};

// p is the reduced pressure [Pa]: Pabs = Pbackground + p + rhoRef*g.(x-xRef).
// Ports accept inward flowRate [m3/s], gauge pressure [Pa], or hydraulic head
// [m] relative to xRef. Reservoir fractions apply only on inflowing faces.
// Example: hydraulicBoundaries { west { type head; value 0.01;
//   inletFractions (0.1 0 0 0 0.16 0 0); } east { type head; value 0; } }
// Optional per-port schedule ((time value) ...) is linearly interpolated.
class HydraulicBoundaries {
    const fvMesh& mesh;
    std::vector<HydraulicPort> ports;
    template<class Type> void replace(GeometricField<Type,fvPatchField,volMesh>& field,
                                     label patch,const dictionary& dict) {
        field.boundaryFieldRef().set(patch,
            fvPatchField<Type>::New(mesh.boundary()[patch],field.internalField(),dict).ptr());
    }
    void add(const word& name,const dictionary& d,const List<scalar>& fallback,
             const word& chemistry,volScalarField& p,volVectorField& U,PtrList<volScalarField>& Y) {
        HydraulicPort port;port.patch=mesh.boundaryMesh().findPatchID(name);
        if(port.patch<0)throw std::invalid_argument("Unknown hydraulic patch: "+std::string(name));
        const auto& patch=mesh.boundary()[port.patch];
        if(patch.coupled()||patch.type()=="empty")throw std::invalid_argument("Hydraulic port must be an external physical patch");
        port.type=d.lookupOrDefault<word>("type","field");
        if(port.type!="flowRate"&&port.type!="pressure"&&port.type!="head"&&port.type!="field")
            throw std::invalid_argument("Hydraulic type must be flowRate, pressure, head or field");
        port.velocityBoundary=d.lookupOrDefault<word>("velocityBoundary","zeroGradient");
        if(port.velocityBoundary!="zeroGradient"&&port.velocityBoundary!="pressureInletOutletVelocity"&&port.velocityBoundary!="pressureNormalInletOutletVelocity"&&port.velocityBoundary!="field")
            throw std::invalid_argument("Unknown hydraulic velocityBoundary");
        port.value=d.lookupOrDefault<scalar>("value",0);
        port.minimumPorosity=d.lookupOrDefault<scalar>("minimumPorosity",0);
        if(!std::isfinite(port.value)||port.minimumPorosity<0||port.minimumPorosity>=1)
            throw std::invalid_argument("Invalid hydraulic port value or minimumPorosity");
        if(d.found("schedule"))d.lookup("schedule")>>port.schedule;
        for(label i=0;i<port.schedule.size();i++)if(!std::isfinite(port.schedule[i].first())||
            !std::isfinite(port.schedule[i].second())||(i&&port.schedule[i].first()<=port.schedule[i-1].first()))
            throw std::invalid_argument("Hydraulic schedule times must increase strictly");
        port.chemistry=d.lookupOrDefault<word>("inletChemistry",chemistry);
        // A supplied incoming concentration specifies total incoming solute
        // flux. Danckwerts' mixed condition tends continuously to zero flux
        // at stagnation. Dirichlet/zeroGradient switching would abruptly
        // switch a finite diffusive flux when a face changes direction.
        port.concentrationBoundary=d.lookupOrDefault<word>("concentrationBoundary","advectiveReservoir");
        if(port.concentrationBoundary!="advectiveReservoir"&&port.concentrationBoundary!="inletOutlet"&&port.concentrationBoundary!="fixedValue"&&port.concentrationBoundary!="field")
            throw std::invalid_argument("concentrationBoundary must be advectiveReservoir, inletOutlet, fixedValue or field");
        if(port.chemistry!="fixedFractions"&&port.chemistry!="haliteSaturated")throw std::invalid_argument("Unknown port inletChemistry");
        List<scalar> fractions(fallback);if(d.found("inletFractions"))d.lookup("inletFractions")>>fractions;
        if(fractions.size()!=7)throw std::invalid_argument("A reservoir requires seven inletFractions");
        double sum=0;for(int j=0;j<7;j++){port.fractions[j]=fractions[j];sum+=fractions[j];if(!std::isfinite(fractions[j])||fractions[j]<0)throw std::invalid_argument("Invalid reservoir fraction");}
        if(sum>=1)throw std::invalid_argument("Reservoir must have positive water fraction");
        if(port.type=="flowRate"){
            // Construct from the type, then explicitly assign patch values.
            p.boundaryFieldRef().set(port.patch,fvPatchField<scalar>::New("fixedFluxPressure",patch,p.internalField()).ptr());
            U.boundaryFieldRef().set(port.patch,fvPatchField<vector>::New("fixedValue",patch,U.internalField()).ptr());
        }else if(port.type=="pressure"||port.type=="head"){
            p.boundaryFieldRef().set(port.patch,fvPatchField<scalar>::New("fixedValue",patch,p.internalField()).ptr());
            // A pressure/head condition does not prescribe incoming
            // tangential velocity. Keep the natural velocity gradient by
            // default. Imposing zero tangent only on backflow introduces a
            // discontinuity at stagnation; normal-entry models remain an
            // explicit physical boundary choice instead of an automatic one.
            if(port.velocityBoundary!="field")U.boundaryFieldRef().set(port.patch,fvPatchField<vector>::New(port.velocityBoundary,patch,U.internalField()).ptr());
        }
        for(int j=0;j<7;j++){
            if(port.concentrationBoundary=="field")continue;
            if(port.concentrationBoundary=="fixedValue"){
                Y[j].boundaryFieldRef().set(port.patch,fvPatchField<scalar>::New("fixedValue",patch,Y[j].internalField()).ptr());
                Y[j].boundaryFieldRef()[port.patch]==port.fractions[j];continue;
            }
            if(port.concentrationBoundary=="advectiveReservoir"){
                OStringStream text;text.precision(17);text<<"type mixed; refValue uniform "<<port.fractions[j]<<"; refGradient uniform 0; valueFraction uniform 0; value uniform "<<port.fractions[j]<<';';
                IStringStream stream(text.str());dictionary yd(stream);replace(Y[j],port.patch,yd);continue;
            }
            OStringStream text;text.precision(17);text<<"type inletOutlet; phi massPhi; inletValue uniform "<<port.fractions[j]<<"; value uniform "<<port.fractions[j]<<';';
            IStringStream stream(text.str());dictionary yd(stream);
            replace(Y[j],port.patch,yd);
        }
        ports.push_back(port);
    }
public:
    HydraulicBoundaries(const dictionary& cfg,const fvMesh& m,const List<scalar>& fallback,
                        const word& chemistry,volScalarField& p,volVectorField& U,PtrList<volScalarField>& Y):mesh(m) {
        if(cfg.found("hydraulicBoundaries")){
            const auto& all=cfg.subDict("hydraulicBoundaries");
            for(const word& name:all.toc())add(name,all.subDict(name),fallback,chemistry,p,U,Y);
        }else if(cfg.found("inletFlow")){
            dictionary d;d.add("type",word("flowRate"));d.add("value",readScalar(cfg.lookup("inletFlow")));d.add("minimumPorosity",scalar(.9));
            add("inlet",d,fallback,chemistry,p,U,Y);
            // Keep the existing outlet pressure value (including spatial data).
            dictionary out;out.add("type",word("field"));add("outlet",out,fallback,chemistry,p,U,Y);
            const label k=mesh.boundaryMesh().findPatchID("outlet");
            U.boundaryFieldRef().set(k,fvPatchField<vector>::New("zeroGradient",mesh.boundary()[k],U.internalField()).ptr());
        }
    }
    double nextKnot(double time,double end)const {
        for(const auto& port:ports)for(const auto& row:port.schedule)if(row.first()>time+1e-12)end=min(end,row.first());
        return end;
    }
    bool hasComponent(int j)const {for(const auto& port:ports)if(port.fractions[j]>0)return true;return false;}
    void updateTransport(const surfaceScalarField& flux,const surfaceScalarField& diffusivity,PtrList<volScalarField>& Y)const {
        for(const auto& port:ports)if(port.concentrationBoundary=="advectiveReservoir"){
            const label k=port.patch;
            const scalarField incoming(max(-flux.boundaryField()[k],scalar(0)));
            const scalarField conductance(diffusivity.boundaryField()[k]*mesh.magSf().boundaryField()[k]*mesh.deltaCoeffs().boundaryField()[k]);
            scalarField fraction(incoming.size(),0);
            forAll(fraction,f)if(incoming[f]>0)fraction[f]=incoming[f]/(incoming[f]+conductance[f]);
            for(int j=0;j<7;j++){
                auto& boundary=refCast<mixedFvPatchScalarField>(Y[j].boundaryFieldRef()[k]);
                boundary.refValue()=port.fractions[j];boundary.refGrad()=0;boundary.valueFraction()=fraction;
            }
        }
    }
    double fluxError(int component,const surfaceScalarField& massFlux,const surfaceScalarField& componentFlux,const surfaceScalarField& diffusivity,const volScalarField& Y)const {
        double error=0;
        for(const auto& port:ports)if(port.concentrationBoundary=="advectiveReservoir"){
            const label k=port.patch;const auto& cells=mesh.boundary()[k].faceCells();
            forAll(cells,f){
                const double q=massFlux.boundaryField()[k][f],yi=Y[cells[f]],yr=port.fractions[component];
                const double expected=q*(q<0?yr:yi);
                const double conductance=diffusivity.boundaryField()[k][f]*mesh.magSf().boundaryField()[k][f]*mesh.deltaCoeffs().boundaryField()[k][f];
                const double scale=(std::abs(q)+conductance)*max(max(std::abs(yi),std::abs(yr)),1e-12);
                error=max(error,std::abs(componentFlux.boundaryField()[k][f]-expected)/max(scale,1e-300));
            }
        }
        reduce(error,maxOp<double>());return error;
    }
    void update(double time,const volScalarField& eps,double rhoRef,const vector& gravity,
                const vector& referencePoint,volScalarField& p,volVectorField& U,
                PtrList<volScalarField>& Y,const saltchem::Speciation& spec) {
        for(auto& port:ports){
            const label k=port.patch;const auto& patch=mesh.boundary()[k];const auto& sf=mesh.Sf().boundaryField()[k];
            const double value=port.at(time);
            if(port.type=="flowRate"){
                double area=0;forAll(sf,f)if(eps[patch.faceCells()[f]]>port.minimumPorosity)area+=mag(sf[f]);reduce(area,sumOp<double>());
                if(area<=0&&value!=0)throw std::runtime_error("Hydraulic port is blocked: "+std::string(patch.name()));
                auto& ub=U.boundaryFieldRef()[k];forAll(ub,f)ub[f]=area>0&&eps[patch.faceCells()[f]]>port.minimumPorosity?-value/area*sf[f]/mag(sf[f]):vector::zero;
            }else if(port.type=="pressure"||port.type=="head"){
                if(port.type=="head"&&mag(gravity)==0)throw std::invalid_argument("Head boundary requires nonzero gravity");
                auto& pb=p.boundaryFieldRef()[k];const auto& cf=mesh.Cf().boundaryField()[k];
                forAll(pb,f)pb[f]=port.type=="head"?rhoRef*mag(gravity)*value:value-rhoRef*(gravity&(cf[f]-referencePoint));
            }
            if(port.chemistry=="haliteSaturated"){
                const auto state=saltchem::haliteSaturatedInlet(spec);const double mass=saltchem::aqueousMass(state);
                for(int j=0;j<7;j++)port.fractions[j]=saltchem::componentMass[j]*state.aqueous[j]/mass;
            }
            for(int j=0;j<7;j++){
                if(port.concentrationBoundary=="inletOutlet")refCast<inletOutletFvPatchScalarField>(Y[j].boundaryFieldRef()[k]).refValue()=port.fractions[j];
                else if(port.concentrationBoundary=="fixedValue")Y[j].boundaryFieldRef()[k]==port.fractions[j];
            }
        }
    }
};
}
