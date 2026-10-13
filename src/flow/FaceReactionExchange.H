// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "PhaseTransform.H"
#include "FixedList.H"
#include "PstreamBuffers.H"
#include "UOPstream.H"
#include "UIPstream.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace saltflow {

// A face sweep is a sequence of matchings: a cell belongs to at most one
// face per colour. Every callback therefore consumes the current endpoint
// inventory, including changes made by earlier colours on another rank.
class FaceReactionExchange {
    const Foam::fvMesh& mesh_;
    struct CoupledFace {Foam::label face,cell,slot;unsigned char direction,colour;};
    std::vector<CoupledFace> coupled_;
    std::vector<Foam::label> coupledFaces_;
    std::array<std::vector<Foam::label>,6> internalGroups_,coupledGroups_;
    std::array<bool,6> globalGroups_{};
    Foam::labelList offsets_;
    Foam::label coupledSize_=0;
    std::vector<unsigned char> directions_,colours_;
    std::array<Foam::label,2> extent_{};
    std::array<bool,2> periodicAxes_{};

    // OpenFOAM's syncTools boundary list also allocates values for both
    // empty extrusion patches. Use the same processor/cyclic pairing with
    // compact coupled-patch storage instead; scalar packets need no spatial
    // transformation, and rotating/nonmatching patches are rejected below.
    template<class T> void swapPacked(Foam::UList<T>& values)const {
        using namespace Foam;
        if(Pstream::parRun()){
            PstreamBuffers buffers(UPstream::commsTypes::nonBlocking);
            DynamicList<label> neighbours;
            forAll(mesh_.boundaryMesh(),k){
                const auto* patch=isA<processorPolyPatch>(mesh_.boundaryMesh()[k]);
                if(!patch||!patch->size())continue;
                neighbours.push_uniq(patch->neighbProcNo());
                const SubList<T> outgoing(values,patch->size(),offsets_[k]);
                UOPstream toNeighbour(patch->neighbProcNo(),buffers);toNeighbour<<outgoing;
            }
            buffers.finishedNeighbourSends(neighbours);
            forAll(mesh_.boundaryMesh(),k){
                const auto* patch=isA<processorPolyPatch>(mesh_.boundaryMesh()[k]);
                if(!patch||!patch->size())continue;
                List<T> received;UIPstream fromNeighbour(patch->neighbProcNo(),buffers);fromNeighbour>>received;
                if(received.size()!=patch->size())FatalErrorInFunction<<"Coupled reaction packet size mismatch"<<abort(FatalError);
                forAll(received,f)values[offsets_[k]+f]=received[f];
            }
        }
        forAll(mesh_.boundaryMesh(),k){
            const auto* patch=Foam::isA<Foam::cyclicPolyPatch>(mesh_.boundaryMesh()[k]);
            if(!patch||!patch->owner())continue;
            const Foam::label other=offsets_[patch->neighbPatchID()];
            forAll(*patch,f)std::swap(values[offsets_[k]+f],values[other+f]);
        }
    }

    static void requireCollectively(int invalid,const char* message){
        Foam::reduce(invalid,Foam::maxOp<int>());
        if(invalid)throw std::invalid_argument(message);
    }
public:
    FaceReactionExchange(const Foam::fvMesh& mesh,const PhaseStencil& stencil,double dx,double dy):
        mesh_(mesh),offsets_(mesh.boundaryMesh().size(),-1),
        directions_(mesh.nFaces(),0),colours_(mesh.nFaces(),6) {
        using namespace Foam;
        int invalid=&stencil.mesh!=&mesh||!(std::isfinite(dx)&&std::isfinite(dy)&&dx>0&&dy>0);
        forAll(mesh.boundaryMesh(),k){
            const auto& patch=mesh.boundaryMesh()[k];if(!patch.coupled())continue;
            const auto* coupled=isA<coupledPolyPatch>(patch);
            invalid|=!(isA<processorPolyPatch>(patch)||isA<cyclicPolyPatch>(patch))||!coupled||!coupled->parallel();
            if(const auto* cyclic=isA<cyclicPolyPatch>(patch))invalid|=cyclic->size()!=cyclic->neighbPatch().size();
            offsets_[k]=coupledSize_;coupledSize_+=patch.size();
        }
        requireCollectively(invalid,"Face reactions require matching processor or translational cyclic patches on the PhaseStencil mesh");
        const std::array<double,2> spacing{{dx,dy}};
        std::array<double,2> origin{};auto& extent=extent_;
        for(int axis=0;axis<2;axis++){
            const scalarField coordinate(mesh.C().primitiveField().component(axis));
            origin[axis]=gMin(coordinate);const double high=gMax(coordinate);
            extent[axis]=label(std::llround((high-origin[axis])/spacing[axis]))+1;
        }
        std::vector<std::array<label,2>> index(mesh.nCells());
        forAll(mesh.C(),cell)for(int axis=0;axis<2;axis++){
            const double coordinate=(mesh.C()[cell][axis]-origin[axis])/spacing[axis];
            index[cell][axis]=label(std::llround(coordinate));invalid|=std::abs(coordinate-index[cell][axis])>1e-7;
        }
        forAll(mesh.neighbour(),face){
            const label own=mesh.owner()[face],other=mesh.neighbour()[face];
            const unsigned char direction=stencil.directions[face];const int axis=direction/2;
            const unsigned char colour=3*axis+(std::min(index[own][axis],index[other][axis])&1);
            directions_[face]=direction;colours_[face]=colour;internalGroups_[colour].push_back(face);
        }
        List<FixedList<label,3>> remote(coupledSize_);
        forAll(mesh.boundaryMesh(),patch)if(offsets_[patch]>=0){
            const auto& faces=mesh.boundary()[patch].faceCells();
            forAll(faces,f){
                const label cell=faces[f],face=mesh.boundaryMesh()[patch].start()+f,slot=offsets_[patch]+f;
                const unsigned char direction=stencil.patchDirections[patch][f];directions_[face]=direction;
                remote[slot][0]=index[cell][0];remote[slot][1]=index[cell][1];remote[slot][2]=direction;
                coupled_.push_back({face,cell,slot,direction,6});
                coupledFaces_.push_back(face);
            }
        }
        swapPacked(remote);
        for(label i=0;i<label(coupled_.size());i++){
            auto& face=coupled_[i];const auto& other=remote[face.slot];const int axis=face.direction/2,sign=(face.direction&1)?1:-1;
            const label delta=other[axis]-index[face.cell][axis],n=extent[axis];
            invalid|=other[2]!=(face.direction^1)||other[axis^1]!=index[face.cell][axis^1];
            unsigned char parity=0;
            if(delta==sign)parity=std::min(index[face.cell][axis],other[axis])&1;
            else if(n>1&&delta==-sign*(n-1)){parity=(n&1)?2:1;periodicAxes_[axis]=true;}
            else invalid=1;
            face.colour=3*axis+parity;colours_[face.face]=face.colour;coupledGroups_[face.colour].push_back(i);
        }
        requireCollectively(invalid,"Coupled face reactions require adjacent Cartesian cells or matching whole-axis periodic endpoints (at least two cells)");
        for(int axis=0;axis<2;axis++){int periodic=periodicAxes_[axis];reduce(periodic,maxOp<int>());periodicAxes_[axis]=periodic;}
        // A topology check, independent of active-face flags, also catches
        // unsupported self-coupling or duplicate connections before stepping.
        for(int colour=0;colour<6;colour++){
            std::vector<unsigned char> used(mesh.nCells(),0);
            for(label face:internalGroups_[colour]){invalid|=used[mesh.owner()[face]]++!=0;invalid|=used[mesh.neighbour()[face]]++!=0;}
            for(label f:coupledGroups_[colour])invalid|=used[coupled_[f].cell]++!=0;
            int present=!internalGroups_[colour].empty()||!coupledGroups_[colour].empty();reduce(present,maxOp<int>());globalGroups_[colour]=present;
        }
        requireCollectively(invalid,"Face reaction colouring is not a matching; shared endpoint inventories would be consumed twice");
    }

    unsigned char direction(Foam::label face)const {return directions_.at(face);}
    unsigned char colour(Foam::label face)const {return colours_.at(face);}
    Foam::label cell(Foam::label face)const {return mesh_.faceOwner()[face];}
    Foam::label coupledFaceCount()const {return coupledSize_;}
    const std::array<bool,2>& periodicAxes()const {return periodicAxes_;}
    const std::array<Foam::label,2>& gridExtent()const {return extent_;}
    const std::vector<Foam::label>& coupledFaces()const {return coupledFaces_;}
    template<unsigned N> void swapCoupled(std::vector<Foam::FixedList<Foam::scalar,N>>& faceValues)const {
        requireCollectively(faceValues.size()!=std::size_t(coupledSize_),"Face reaction halo requires one packet per compact coupled face");
        Foam::UList<Foam::FixedList<Foam::scalar,N>> values(faceValues.data(),faceValues.size());swapPacked(values);
    }

    // The caller owns the packet layout. The callback may only mutate its
    // two packets; stage-wide rollback of earlier successful faces belongs
    // to the caller. A failing face commits neither endpoint, and all ranks
    // complete the return exchange before reporting a callback exception.
    template<unsigned N,class Callback>
    void execute(std::vector<Foam::FixedList<Foam::scalar,N>>& cells,
                 const Foam::boolList& activeFaces,Callback&& callback)const {
        using namespace Foam;using Packet=FixedList<scalar,N>;
        requireCollectively(cells.size()!=std::size_t(mesh_.nCells())||activeFaces.size()!=mesh_.nFaces(),
                            "Face reaction executor requires one packet per cell and one activity flag per mesh face");
        List<unsigned char> active(coupledSize_,0);
        for(const auto& face:coupled_)active[face.slot]=activeFaces[face.face];
        swapPacked(active);
        for(const auto& face:coupled_)active[face.slot]=active[face.slot]||activeFaces[face.face];
        List<Packet> halo(coupledSize_);
        for(int colour=0;colour<6;colour++){
            if(!globalGroups_[colour])continue;
            int activeCount=0;
            for(label face:internalGroups_[colour])activeCount+=activeFaces[face];
            for(label f:coupledGroups_[colour])activeCount+=active[coupled_[f].slot];
            reduce(activeCount,sumOp<int>());if(!activeCount)continue;
            for(const auto& face:coupled_)halo[face.slot]=cells[face.cell];
            swapPacked(halo);
            int failed=0;std::string firstError;
            auto apply=[&](label face,bool boundary,Packet& owner,Packet& neighbour){
                Packet a(owner),b(neighbour);
                try{
                    callback(face,boundary,a,b);
                    for(unsigned k=0;k<N;k++)if(!std::isfinite(a[k])||!std::isfinite(b[k]))throw std::runtime_error("Nonfinite face reaction packet");
                    owner=a;neighbour=b;
                }catch(const std::exception& error){if(!failed)firstError=error.what();failed++;}
                catch(...){if(!failed)firstError="Unknown face reaction callback failure";failed++;}
            };
            for(label face:internalGroups_[colour])if(activeFaces[face])apply(face,false,cells[mesh_.owner()[face]],cells[mesh_.neighbour()[face]]);
            for(label f:coupledGroups_[colour]){
                const auto& face=coupled_[f];
                if(active[face.slot]&&(face.direction&1))apply(face.face,true,cells[face.cell],halo[face.slot]);
            }
            swapPacked(halo);
            for(label f:coupledGroups_[colour]){
                const auto& face=coupled_[f];if(active[face.slot]&&!(face.direction&1))cells[face.cell]=halo[face.slot];
            }
            if(failed)Pout<<"FACE_REACTION colour="<<colour<<" failures="<<failed<<" "<<firstError<<endl;
            reduce(failed,sumOp<int>());
            if(failed)throw std::runtime_error("Face reaction callback failed in "+std::to_string(failed)+" faces after collective exchange");
        }
    }
};
}
