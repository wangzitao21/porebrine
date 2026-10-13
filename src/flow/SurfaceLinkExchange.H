// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "PhaseTransform.H"
#include "FixedList.H"
#include "PstreamBuffers.H"
#include "UOPstream.H"
#include "UIPstream.H"
#include <array>
#include <cmath>
#include <string>
#include <tuple>
#include <vector>

namespace saltflow {
struct SurfaceLink {
    Foam::label source,receiverRank,receiverCell,gx,gy,rx,ry;
};

// Links describe an already validated liquid path supplied by the geometry
// layer. This executor only validates endpoints and schedules conservative
// updates; proximity alone never establishes liquid connectivity.
class SurfaceLinkExchange {
    using Label=Foam::label;
    using Request=Foam::FixedList<Label,6>;
    static constexpr int colourCount=36;
    struct Partner {
        Label rank;
        Label incomingOffset=0;
        std::vector<Label> outgoing;
        std::vector<Request> incoming;
        std::array<std::vector<Label>,colourCount> outColour,inColour;
    };
    const Foam::fvMesh& mesh_;
    const std::array<bool,2> periodic_;
    std::array<Label,2> extent_{};
    std::vector<SurfaceLink> links_;
    std::vector<unsigned char> colours_;
    std::array<std::vector<Label>,colourCount> groups_;
    std::array<bool,colourCount> globalGroups_{};
    std::vector<Partner> partners_;
    Foam::labelList neighbours_;
    std::vector<Label> remoteSlot_;
    Label remoteCount_=0;
    struct Contribution {Label cell,gx,gy,index;bool local;};
    std::vector<Contribution> scatterOrder_;
    Label incomingCount_=0;

    static void requireCollectively(int invalid,const char* message){
        Foam::reduce(invalid,Foam::maxOp<int>());
        if(invalid)throw std::invalid_argument(message);
    }
    int colourOf(Label gx,Label gy,Label rx,Label ry)const {
        static constexpr std::array<std::array<int,2>,12> displacement{{
            {{1,0}},{{-1,0}},{{0,1}},{{0,-1}},{{2,0}},{{-2,0}},
            {{0,2}},{{0,-2}},{{1,1}},{{1,-1}},{{-1,1}},{{-1,-1}}
        }};
        const std::array<Label,2> raw{{rx-gx,ry-gy}};auto canonical=raw;
        for(int axis=0;axis<2;axis++)if(periodic_[axis]&&std::abs(canonical[axis])>extent_[axis]/2)
            canonical[axis]+=canonical[axis]>0?-extent_[axis]:extent_[axis];
        for(int direction=0;direction<12;direction++){
            const auto& d=displacement[direction];
            if(canonical[0]==d[0]&&canonical[1]==d[1]){
                const int axis=d[0]?0:1;
                const Label origin=d[0]?std::min(gx,rx):std::min(gy,ry);
                const Label stride=std::abs(d[0]?d[0]:d[1]);
                // Half-domain ties retain the raw direction. A shortest
                // periodic edge crossing the selected axis gets its own
                // matching; wrapping on the other axis cannot change which
                // endpoints share the selected-axis parity.
                return 3*direction+(raw[axis]!=canonical[axis]?2:((origin/stride)&1));
            }
        }
        return -1;
    }

    // The neighbour graph and request order are fixed at construction.
    // Each later exchange contains only endpoint packets, without global
    // cell gathers or repeated discovery of MPI partners.
    template<unsigned N,class Fill,class Receive>
    void communicate(int colour,bool toSources,Fill&& fill,Receive&& receive)const {
        using namespace Foam;using Packet=FixedList<scalar,N>;
        if(!Pstream::parRun())return;
        PstreamBuffers buffers(UPstream::commsTypes::nonBlocking);
        for(const auto& partner:partners_){
            const auto& indices=toSources?(colour<0?partner.inColour[0]:partner.inColour[colour]):(colour<0?partner.outgoing:partner.outColour[colour]);
            const Label count=colour<0&&toSources?partner.incoming.size():indices.size();
            if(!count)continue;
            List<Packet> values(count);
            for(Label k=0;k<count;k++)fill(partner,colour<0&&toSources?k:indices[k],values[k]);
            UOPstream stream(partner.rank,buffers);stream<<values;
        }
        buffers.finishedNeighbourSends(neighbours_);
        for(const auto& partner:partners_){
            const auto& indices=toSources?(colour<0?partner.outgoing:partner.outColour[colour]):(colour<0?partner.inColour[0]:partner.inColour[colour]);
            const Label count=colour<0&&!toSources?partner.incoming.size():indices.size();
            if(!count)continue;
            List<Packet> values;UIPstream stream(partner.rank,buffers);stream>>values;
            if(values.size()!=count)FatalErrorInFunction<<"Surface link packet count mismatch"<<abort(FatalError);
            for(Label k=0;k<count;k++)receive(partner,colour<0&&!toSources?k:indices[k],values[k]);
        }
    }

public:
    bool matches(const std::vector<SurfaceLink>& candidate)const {
        if(candidate.size()!=links_.size())return false;
        for(std::size_t i=0;i<candidate.size();i++){
            const auto& a=candidate[i];const auto& b=links_[i];
            if(a.source!=b.source||a.receiverRank!=b.receiverRank||a.receiverCell!=b.receiverCell||
               a.gx!=b.gx||a.gy!=b.gy||a.rx!=b.rx||a.ry!=b.ry)return false;
        }
        return true;
    }
    SurfaceLinkExchange(const Foam::fvMesh& mesh,const std::vector<SurfaceLink>& links,double dx,double dy,
                        std::array<bool,2> periodicFlags={{false,false}}):
        mesh_(mesh),periodic_(periodicFlags),links_(links),colours_(links.size()),remoteSlot_(links.size(),-1) {
        using namespace Foam;
        PhaseStencil geometryCheck(mesh,dx,dy);
        const std::array<double,2> spacing{{dx,dy}};std::array<double,2> origin{};auto& extent=extent_;
        std::vector<std::array<Label,2>> coordinates(mesh.nCells());int invalid=0;
        for(int axis=0;axis<2;axis++){
            const scalarField coordinate(mesh.C().primitiveField().component(axis));
            origin[axis]=gMin(coordinate);extent[axis]=Label(std::llround((gMax(coordinate)-origin[axis])/spacing[axis]))+1;
        }
        forAll(mesh.C(),cell)for(int axis=0;axis<2;axis++){
            const double index=(mesh.C()[cell][axis]-origin[axis])/spacing[axis];coordinates[cell][axis]=Label(std::llround(index));invalid|=std::abs(index-coordinates[cell][axis])>1e-7;
        }
        labelList sendCounts(Pstream::nProcs(),0),receiveCounts(Pstream::nProcs(),0);
        for(Label i=0;i<Label(links_.size());i++){
            const auto& link=links_[i];
            const bool valid=link.source>=0&&link.source<mesh.nCells()&&link.receiverRank>=0&&link.receiverRank<Pstream::nProcs()&&link.receiverCell>=0&&
                link.gx>=0&&link.gx<extent[0]&&link.rx>=0&&link.rx<extent[0]&&link.gy>=0&&link.gy<extent[1]&&link.ry>=0&&link.ry<extent[1];
            const int colour=valid?colourOf(link.gx,link.gy,link.rx,link.ry):-1;
            if(!valid||colour<0){invalid=1;continue;}
            invalid|=coordinates[link.source][0]!=link.gx||coordinates[link.source][1]!=link.gy;
            colours_[i]=colour;groups_[colour].push_back(i);
            if(link.receiverRank==Pstream::myProcNo()){
                if(link.receiverCell>=mesh.nCells())invalid=1;
                else invalid|=coordinates[link.receiverCell][0]!=link.rx||coordinates[link.receiverCell][1]!=link.ry;
            }else{sendCounts[link.receiverRank]++;remoteSlot_[i]=remoteCount_++;}
        }
        requireCollectively(invalid,"Surface links require valid local sources and in-domain Cartesian receivers one or two grid steps away under the configured periodic axes");
        UPstream::allToAll(sendCounts,receiveCounts);
        labelList partnerIndex(Pstream::nProcs(),-1);
        for(Label rank=0;rank<Pstream::nProcs();rank++)if(sendCounts[rank]||receiveCounts[rank]){
            partnerIndex[rank]=partners_.size();partners_.push_back(Partner{});partners_.back().rank=rank;
        }
        neighbours_.setSize(partners_.size());for(Label p=0;p<Label(partners_.size());p++)neighbours_[p]=partners_[p].rank;
        for(Label i=0;i<Label(links_.size());i++)if(remoteSlot_[i]>=0){
            auto& partner=partners_[partnerIndex[links_[i].receiverRank]];partner.outgoing.push_back(i);partner.outColour[colours_[i]].push_back(i);
        }
        if(Pstream::parRun()){
            PstreamBuffers buffers(UPstream::commsTypes::nonBlocking);
            for(const auto& partner:partners_)if(!partner.outgoing.empty()){
                List<Request> requests(partner.outgoing.size());
                forAll(requests,k){const auto& link=links_[partner.outgoing[k]];requests[k]=Request({link.receiverCell,link.gx,link.gy,link.rx,link.ry,Label(colours_[partner.outgoing[k]])});}
                UOPstream stream(partner.rank,buffers);stream<<requests;
            }
            buffers.finishedNeighbourSends(neighbours_);
            for(auto& partner:partners_)if(receiveCounts[partner.rank]){
                List<Request> requests;UIPstream stream(partner.rank,buffers);stream>>requests;
                if(requests.size()!=receiveCounts[partner.rank])FatalErrorInFunction<<"Surface link request count mismatch"<<abort(FatalError);
                partner.incoming.assign(requests.begin(),requests.end());
                forAll(requests,k){
                    const auto& r=requests[k];const int colour=colourOf(r[1],r[2],r[3],r[4]);
                    if(r[0]<0||r[0]>=mesh.nCells()||colour<0||colour!=r[5]){invalid=1;continue;}
                    invalid|=coordinates[r[0]][0]!=r[3]||coordinates[r[0]][1]!=r[4];partner.inColour[colour].push_back(k);
                }
            }
        }
        requireCollectively(invalid,"Surface link receiver cell or Cartesian coordinates do not match the receiver partition");
        for(Label i=0;i<Label(links_.size());i++)if(remoteSlot_[i]<0){
            const auto& link=links_[i];scatterOrder_.push_back({link.receiverCell,link.gx,link.gy,i,true});
        }
        for(auto& partner:partners_){
            partner.incomingOffset=incomingCount_;
            for(const auto& request:partner.incoming)scatterOrder_.push_back({request[0],request[1],request[2],incomingCount_++,false});
        }
        // Geometric support feeds stochastic nucleation. Sum in the same
        // global source order on every decomposition, including when local
        // and remote sources exchange roles after repartitioning.
        std::sort(scatterOrder_.begin(),scatterOrder_.end(),[](const Contribution& a,const Contribution& b){
            return std::tie(a.cell,a.gx,a.gy)<std::tie(b.cell,b.gx,b.gy);
        });
        for(int colour=0;colour<colourCount;colour++){
            std::vector<unsigned char> used(mesh.nCells(),0);
            for(Label i:groups_[colour]){
                const auto& link=links_[i];invalid|=used[link.source]++!=0;
                if(link.receiverRank==Pstream::myProcNo())invalid|=used[link.receiverCell]++!=0;
            }
            for(const auto& partner:partners_)for(Label k:partner.inColour[colour])invalid|=used[partner.incoming[k][0]]++!=0;
            int present=!groups_[colour].empty();reduce(present,maxOp<int>());globalGroups_[colour]=present;
        }
        requireCollectively(invalid,"Surface link colouring has duplicate/shared endpoints; deduplicate source-to-receiver links before construction");
    }

    const std::vector<SurfaceLink>& links()const {return links_;}
    unsigned char colour(Label link)const {return colours_.at(link);}

    template<unsigned N>
    std::vector<Foam::FixedList<Foam::scalar,N>> gatherTargets(const std::vector<Foam::FixedList<Foam::scalar,N>>& cells)const {
        using Packet=Foam::FixedList<Foam::scalar,N>;
        requireCollectively(cells.size()!=std::size_t(mesh_.nCells()),"Surface target gather requires one packet per local cell");
        std::vector<Packet> result(links_.size());
        for(Label i=0;i<Label(links_.size());i++)if(remoteSlot_[i]<0)result[i]=cells[links_[i].receiverCell];
        communicate<N>(-1,true,[&](const Partner& p,Label k,Packet& value){value=cells[p.incoming[k][0]];},
                       [&](const Partner&,Label i,const Packet& value){result[i]=value;});
        return result;
    }

    template<unsigned N>
    void scatterAddTargets(const std::vector<Foam::FixedList<Foam::scalar,N>>& amounts,
                           std::vector<Foam::FixedList<Foam::scalar,N>>& cells)const {
        using Packet=Foam::FixedList<Foam::scalar,N>;
        int invalid=cells.size()!=std::size_t(mesh_.nCells())||amounts.size()!=links_.size();
        for(const auto& packet:amounts)for(double value:packet)invalid|=!std::isfinite(value);
        requireCollectively(invalid,"Surface target scatter requires finite per-link amounts and one packet per local cell");
        auto add=[&](Label cell,const Packet& value){for(unsigned k=0;k<N;k++){const double sum=cells[cell][k]+value[k];if(!std::isfinite(sum))invalid=1;else cells[cell][k]=sum;}};
        std::vector<Packet> received(incomingCount_);
        communicate<N>(-1,false,[&](const Partner&,Label i,Packet& value){value=amounts[i];},
                       [&](const Partner& p,Label k,const Packet& value){received[p.incomingOffset+k]=value;});
        for(const auto& entry:scatterOrder_)add(entry.cell,entry.local?amounts[entry.index]:received[entry.index]);
        requireCollectively(invalid,"Surface target scatter produced a nonfinite accumulated packet");
    }

    template<unsigned N,class Callback>
    void execute(std::vector<Foam::FixedList<Foam::scalar,N>>& cells,const Foam::boolList& activeLinks,Callback&& callback)const {
        using namespace Foam;using Packet=FixedList<scalar,N>;
        requireCollectively(cells.size()!=std::size_t(mesh_.nCells())||activeLinks.size()!=Label(links_.size()),
                            "Surface link execution requires one packet per local cell and one flag per local link");
        std::vector<Packet> targets(remoteCount_);
        for(int colour=0;colour<colourCount;colour++){
            if(!globalGroups_[colour])continue;
            int active=0;for(Label i:groups_[colour])active+=activeLinks[i];reduce(active,sumOp<int>());if(!active)continue;
            communicate<N>(colour,true,[&](const Partner& p,Label k,Packet& value){value=cells[p.incoming[k][0]];},
                           [&](const Partner&,Label i,const Packet& value){targets[remoteSlot_[i]]=value;});
            int failed=0;std::string firstError;
            for(Label i:groups_[colour])if(activeLinks[i]){
                const auto& link=links_[i];auto& source=cells[link.source];
                auto& target=remoteSlot_[i]<0?cells[link.receiverCell]:targets[remoteSlot_[i]];
                Packet a(source),b(target);
                try{
                    callback(i,a,b);
                    for(unsigned k=0;k<N;k++)if(!std::isfinite(a[k])||!std::isfinite(b[k]))throw std::runtime_error("Nonfinite surface link packet");
                    source=a;target=b;
                }catch(const std::exception& error){if(!failed)firstError=error.what();failed++;}
                catch(...){if(!failed)firstError="Unknown surface link callback failure";failed++;}
            }
            communicate<N>(colour,false,[&](const Partner&,Label i,Packet& value){value=targets[remoteSlot_[i]];},
                           [&](const Partner& p,Label k,const Packet& value){cells[p.incoming[k][0]]=value;});
            if(failed)Pout<<"SURFACE_LINK colour="<<colour<<" failures="<<failed<<" "<<firstError<<endl;
            reduce(failed,sumOp<int>());if(failed)throw std::runtime_error("Surface link callback failed in "+std::to_string(failed)+" links after collective exchange");
        }
    }
};
}
