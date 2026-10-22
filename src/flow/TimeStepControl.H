// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "fvCFD.H"
#include "BrineSystem.hpp"
#include "MineralInterface.hpp"
#include "DualMineralGeometry.hpp"

namespace saltflow {
struct StepConvergenceFailure:std::runtime_error {using std::runtime_error::runtime_error;};
template<class Type,template<class> class PatchField,class GeoMesh> struct FieldState {
    Foam::Field<Type> internal;
    Foam::List<Foam::Field<Type>> boundary;
    using Field=Foam::GeometricField<Type,PatchField,GeoMesh>;
    FieldState()=default;
    explicit FieldState(const Field& field):internal(field.primitiveField()),boundary(field.boundaryField().size()) {
        forAll(boundary,k)boundary[k]=field.boundaryField()[k];
    }
    void restore(Field& field)const {
        field.primitiveFieldRef()=internal;
        forAll(boundary,k)field.boundaryFieldRef()[k]==boundary[k];
    }
};
using ScalarState=FieldState<Foam::scalar,Foam::fvPatchField,Foam::volMesh>;
using VectorState=FieldState<Foam::vector,Foam::fvPatchField,Foam::volMesh>;
using FluxState=FieldState<Foam::scalar,Foam::fvsPatchField,Foam::surfaceMesh>;
struct FlowState {
    double time,maxBalance;
    double pressureOffset=0;
    Foam::scalarField pressureDeviation;
    Foam::label index;
    Foam::TimeState nativeTime;
    ScalarState p,eps,rho;
    VectorState U;
    FluxState phi,massPhi;
    std::array<ScalarState,7> Y,solid;
    std::array<double,8> boundaryIntegral,balance;
    saltchem::MineralVector dissolved,precipitated;
    std::vector<saltgeom::SurfaceState> mineralSurface;
    std::vector<saltchem::MineralVector> mineralSurfaceRates;
    std::vector<std::array<double,2>> mineralSurfaceVectors;
    std::shared_ptr<saltdual::Geometry> dualGeometry;
    std::vector<std::array<double,2>> dualGeometryVectors;
    FlowState(const Foam::Time& time,const Foam::volScalarField& pressure,
              const Foam::volScalarField& porosity,const Foam::volScalarField& density,
              const Foam::volVectorField& velocity,const Foam::surfaceScalarField& volumeFlux,
              const Foam::surfaceScalarField& massFlux,const Foam::PtrList<Foam::volScalarField>& components,
              const Foam::PtrList<Foam::volScalarField>& minerals):
        time(time.value()),index(time.timeIndex()),nativeTime(time),p(pressure),eps(porosity),rho(density),U(velocity),phi(volumeFlux),massPhi(massFlux) {
        for(int j=0;j<7;j++){Y[j]=ScalarState(components[j]);solid[j]=ScalarState(minerals[j]);}
    }
};

class TimeStepControl {
public:
    bool enabled,recoverFailedSteps;
    double relative,porosityAbsolute,inventoryAbsolute,waterAbsolute,velocityAbsolute,minimum,maximum;
    int maximumRejections;
    explicit TimeStepControl(const Foam::dictionary& cfg,double initialDt):
        enabled(cfg.lookupOrDefault<bool>("adaptiveTimeStep",false)),
        recoverFailedSteps(cfg.lookupOrDefault<bool>("recoverFailedTimeSteps",true)),
        relative(cfg.lookupOrDefault<Foam::scalar>("timeRelativeTolerance",1e-3)),
        porosityAbsolute(cfg.lookupOrDefault<Foam::scalar>("timePorosityTolerance",1e-4)),
        inventoryAbsolute(cfg.lookupOrDefault<Foam::scalar>("timeInventoryTolerance",1e-3)),
        waterAbsolute(cfg.lookupOrDefault<Foam::scalar>("timeWaterTolerance",1e-4)),
        velocityAbsolute(cfg.lookupOrDefault<Foam::scalar>("timeVelocityTolerance",1e-8)),
        minimum(cfg.lookupOrDefault<Foam::scalar>("minDeltaT",1e-8)),
        maximum(cfg.lookupOrDefault<Foam::scalar>("maxDeltaT",initialDt)),
        maximumRejections(cfg.lookupOrDefault<Foam::label>("maximumTimeStepRejections",20)) {
        for(double v:{relative,porosityAbsolute,inventoryAbsolute,waterAbsolute,velocityAbsolute,minimum,maximum})
            if(!std::isfinite(v)||v<=0)throw std::invalid_argument("Time tolerances and step bounds must be finite and positive");
        if(maximum<minimum||maximumRejections<1)throw std::invalid_argument("Invalid time-step bounds/rejection limit");
    }
    // Euler + sequential reaction/transport is first order. The difference
    // between one full step and two half steps estimates the fine-step error.
    // Use a global MAX over cells/components: a thin moving interface cannot
    // disappear in an RMS average. Accept the two conservative half steps;
    // never extrapolate inventories or accept an over-tolerance minimum step.
    double error(const FlowState& coarse,const Foam::volScalarField& eps,
                 const Foam::volScalarField& rho,const Foam::volVectorField& U,
                 const Foam::PtrList<Foam::volScalarField>& Y,
                 const Foam::PtrList<Foam::volScalarField>& solid,const saltchem::MineralVector& volume)const {
        std::array<double,4> errors{};
        auto add=[&](double fine,double rough,double absolute,int category){
            if(!std::isfinite(fine)||!std::isfinite(rough)){errors[category]=Foam::GREAT;return;}
            errors[category]=std::max(errors[category],std::abs(fine-rough)/(absolute+relative*std::max(std::abs(fine),std::abs(rough))));
        };
        const double velocityScale=std::max(Foam::gMax(Foam::mag(U.primitiveField())),Foam::gMax(Foam::mag(coarse.U.internal)));
        forAll(eps,i){
            add(eps[i],coarse.eps.internal[i],porosityAbsolute,0);
            const double velocityError=Foam::mag(U[i]-coarse.U.internal[i])/(velocityAbsolute+relative*velocityScale);
            errors[1]=std::isfinite(velocityError)?std::max(errors[1],velocityError):Foam::GREAT;
            const double cap=eps[i]*rho[i],oldCap=coarse.eps.internal[i]*coarse.rho.internal[i];double water=cap,oldWater=oldCap;
            for(int j=0;j<7;j++){
                add(cap*Y[j][i]/saltchem::componentMass[j],oldCap*coarse.Y[j].internal[i]/saltchem::componentMass[j],inventoryAbsolute,2);
                add(solid[j][i]/volume[j],coarse.solid[j].internal[i]/volume[j],inventoryAbsolute,2);
                water-=cap*Y[j][i];oldWater-=oldCap*coarse.Y[j].internal[i];
            }
            add(water,oldWater,waterAbsolute,3);
        }
        for(double& e:errors)Foam::reduce(e,Foam::maxOp<double>());
        Foam::Info<<"TIME_ERROR porosity="<<errors[0]<<" velocity="<<errors[1]<<" inventory="<<errors[2]<<" water="<<errors[3]<<Foam::endl;
        return *std::max_element(errors.begin(),errors.end());
    }
    double factor(double error)const {return std::clamp(.9/std::sqrt(std::max(error,1e-12)),.2,2.);}
};
}
