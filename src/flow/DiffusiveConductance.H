// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "fvCFD.H"
#include <algorithm>
#include <cmath>

namespace saltflow {
// Two half-cell diffusion resistances in series. The mesh interpolation
// weight is the neighbour distance / centre distance, hence the opposite
// weights on the two resistivities. Scaling by the smaller coefficient
// avoids overflowing 1/k at a nearly impermeable solid.
inline double harmonicConductance(double owner,double neighbour,double ownerWeight=.5){
    if(!std::isfinite(owner)||!std::isfinite(neighbour)||owner<0||neighbour<0||
       !std::isfinite(ownerWeight)||ownerWeight<0||ownerWeight>1)
        throw std::invalid_argument("Invalid diffusion coefficient or face weight");
    if(ownerWeight==0)return owner;
    if(ownerWeight==1)return neighbour;
    if(owner==0||neighbour==0)return 0;
    const double scale=std::min(owner,neighbour);
    return scale/((1-ownerWeight)*(scale/owner)+ownerWeight*(scale/neighbour));
}

inline Foam::tmp<Foam::surfaceScalarField> diffusiveConductance(const Foam::volScalarField& coefficient){
    const auto& mesh=coefficient.mesh();
    // Keep physical-patch coefficients supplied by the field. Coupled
    // patches use exactly the same two-point operator as internal faces.
    auto result=Foam::fvc::interpolate(coefficient);
    auto& face=result.ref();const auto& weights=mesh.surfaceInterpolation::weights();
    forAll(mesh.neighbour(),f)
        face[f]=harmonicConductance(coefficient[mesh.owner()[f]],coefficient[mesh.neighbour()[f]],weights[f]);
    forAll(coefficient.boundaryField(),p){
        const auto& patch=coefficient.boundaryField()[p];
        if(patch.coupled()){
            const auto other=patch.patchNeighbourField();const auto& cells=mesh.boundary()[p].faceCells();
            forAll(cells,f)face.boundaryFieldRef()[p][f]=harmonicConductance(coefficient[cells[f]],other()[f],weights.boundaryField()[p][f]);
        }
    }
    return result;
}
}
