// SPDX-License-Identifier: GPL-3.0-or-later
// Wang, Hu & Steefel (2024), doi:10.1029/2023WR036122, equations 14,17,18.
#pragma once
#include "fvCFD.H"
#include <array>
#include <vector>
#include "PhaseVector.hpp"

struct PhaseStencil {
    const Foam::fvMesh& mesh;
    std::vector<unsigned char> directions;
    std::vector<std::vector<unsigned char>> patchDirections;
    PhaseStencil(const Foam::fvMesh& m,double dx,double dy):mesh(m),directions(m.nInternalFaces()),patchDirections(m.boundary().size()) {
        // The phase-vector model assumes one uniform Cartesian layer in XY.
        // Reject unsupported geometry before solving; cache only topology,
        // never the porosity values exchanged across processor boundaries.
        int invalid=!(std::isfinite(dx)&&std::isfinite(dy)&&dx>0&&dy>0);
        const double scale=std::max(dx,dy),tol=scale*1e-7;
        const auto z=mesh.C().primitiveField().component(Foam::vector::Z);
        const double zmin=Foam::gMin(z()),zmax=Foam::gMax(z());
        const double vmin=Foam::gMin(mesh.V()),vmax=Foam::gMax(mesh.V());
        invalid|=zmax-zmin>tol||!(vmin>0)||vmax-vmin>vmax*1e-7;
        forAll(mesh.cells(),i)invalid|=mesh.cells()[i].size()!=6;
        const auto& own=mesh.owner();const auto& nei=mesh.neighbour();
        forAll(nei,f){
            const auto d=mesh.C()[nei[f]]-mesh.C()[own[f]];
            const bool x=std::abs(d.x())>std::abs(d.y());
            invalid|=std::abs(d.z())>tol||std::abs(x?d.y():d.x())>tol||std::abs(std::abs(x?d.x():d.y())-(x?dx:dy))>tol;
            directions[f]=(x?0:2)+((x?d.x():d.y())>0);
        }
        auto axisAligned=[&](const Foam::vector& s){
            const double a=Foam::mag(s),largest=std::max({std::abs(s.x()),std::abs(s.y()),std::abs(s.z())});
            return a>0&&a-largest<=a*1e-12;
        };
        for(const auto& s:mesh.Sf().primitiveField())invalid|=!axisAligned(s)||std::abs(s.z())>Foam::mag(s)*1e-7;
        forAll(mesh.boundary(),p){
            const auto& patch=mesh.boundary()[p];const auto& sf=mesh.Sf().boundaryField()[p];
            patchDirections[p].resize(patch.size());
            forAll(sf,f){
                const bool z=std::abs(sf[f].z())>Foam::mag(sf[f])*0.5;
                invalid|=!axisAligned(sf[f])||(z!=(patch.type()=="empty"));
                const bool x=std::abs(sf[f].x())>std::abs(sf[f].y());
                patchDirections[p][f]=(x?0:2)+((x?sf[f].x():sf[f].y())>0);
            }
        }
        Foam::reduce(invalid,Foam::maxOp<int>());
        if(invalid)throw std::invalid_argument("Phase geometry requires a uniform, axis-aligned XY mesh with one empty Z layer and matching dx/dy");
    }
    std::vector<std::array<double,4>> values(const Foam::volScalarField& e) const {
        std::vector<std::array<double,4>> v(mesh.nCells());
        forAll(e,i) v[i]={{e[i],e[i],e[i],e[i]}};
        const auto& own=mesh.owner(); const auto& nei=mesh.neighbour();
        forAll(nei,f) {
            const int direction=directions[f];
            v[own[f]][direction]=e[nei[f]];v[nei[f]][direction^1]=e[own[f]];
        }
        forAll(e.boundaryField(),p) if(e.boundaryField()[p].coupled()) {
            Foam::tmp<Foam::scalarField> other=e.boundaryField()[p].patchNeighbourField();
            const auto& cells=mesh.boundary()[p].faceCells();
            forAll(cells,f) {
                v[cells[f]][patchDirections[p][f]]=other()[f];
            }
        }
        return v;
    }
    static std::array<double,3> evaluate(double e,const std::array<double,4>& v,
                                        double dx,double dy,double Q) {
        return phaseVector(e,v,dx,dy,Q);
    }
};
