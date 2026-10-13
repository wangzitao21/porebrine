// SPDX-License-Identifier: GPL-3.0-or-later
#include "SurfaceLinkExchange.H"
#include "FaceReactionExchange.H"
#include <fstream>
#include <iomanip>
#include <set>

int main(int argc,char** argv){
    using namespace Foam;using saltflow::SurfaceLink;using saltflow::SurfaceLinkExchange;
    argList args(argc,argv);if(!args.checkRootCase())return 1;
    Time runTime(Time::controlDictName,args);fvMesh mesh(IOobject(fvMesh::defaultRegion,runTime.timeName(),runTime,IOobject::MUST_READ));
    IOdictionary config(IOobject("linkTest",runTime.constant(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE));
    const label nx=readLabel(config.lookup("nx")),ny=readLabel(config.lookup("ny"));const double spacing=readScalar(config.lookup("spacing"));
    const std::array<bool,2> requestedPeriodic{{config.lookupOrDefault<bool>("periodicX",false),config.lookupOrDefault<bool>("periodicY",false)}};
    PhaseStencil stencil(mesh,spacing,spacing);saltflow::FaceReactionExchange faceExchange(mesh,stencil,spacing,spacing);
    const auto periodic=faceExchange.periodicAxes();if(periodic!=requestedPeriodic||faceExchange.gridExtent()!=std::array<label,2>{{nx,ny}})return 13;
    // Small independent fixture metadata; the production executor itself
    // discovers only sparse receiver requests, never global cell packets.
    labelList ranks(nx*ny,-1),indices(nx*ny,-1),ids(mesh.nCells());
    forAll(mesh.C(),i){const label gx=std::llround(mesh.C()[i].x()/spacing-.5),gy=std::llround(mesh.C()[i].y()/spacing-.5);ids[i]=gx+nx*gy;ranks[ids[i]]=Pstream::myProcNo();indices[ids[i]]=i;}
    for(label id=0;id<nx*ny;id++){reduce(ranks[id],maxOp<label>());reduce(indices[id],maxOp<label>());}
    const std::array<std::array<label,2>,12> displacements{{{{1,0}},{{-1,0}},{{0,1}},{{0,-1}},{{2,0}},{{-2,0}},{{0,2}},{{0,-2}},{{1,1}},{{1,-1}},{{-1,1}},{{-1,-1}}}};
    std::vector<SurfaceLink> links;
    forAll(ids,i){std::set<label> targets;for(const auto& d:displacements){
        const label gx=ids[i]%nx,gy=ids[i]/nx;label rx=gx+d[0],ry=gy+d[1];
        if(periodic[0])rx=(rx+nx)%nx;if(periodic[1])ry=(ry+ny)%ny;
        if(rx>=0&&rx<nx&&ry>=0&&ry<ny){const label target=rx+nx*ry;if(target!=ids[i]&&targets.insert(target).second)links.push_back({i,ranks[target],indices[target],gx,gy,rx,ry});}
    }
    }
    SurfaceLinkExchange exchange(mesh,links,spacing,spacing,periodic);
    auto distance=[&](label gx,label gy,label rx,label ry){label x=std::abs(gx-rx),y=std::abs(gy-ry);if(periodic[0])x=min(x,nx-x);if(periodic[1])y=min(y,ny-y);return x+y;};
    using Packet=FixedList<scalar,31>;
    std::vector<Packet> cells(mesh.nCells(),Packet(Zero));
    forAll(ids,i){const label id=ids[i];cells[i][30]=id%4==0;cells[i][23]=id;cells[i][22]=1234.56789;cells[i][15]=.9;if(cells[i][30]){for(int j=0;j<7;j++)cells[i][j]=(j+1)*(.4+.1*(id%5));cells[i][14]=cells[i][0];}}
    const auto initial=cells;
    const auto gathered=exchange.gatherTargets(cells);
    for(std::size_t i=0;i<links.size();i++)if(gathered[i][23]!=links[i].rx+nx*links[i].ry)return 2;
    std::vector<FixedList<scalar,7>> areas(links.size(),FixedList<scalar,7>(Zero)),cellAreas(mesh.nCells(),FixedList<scalar,7>(Zero));
    for(auto& area:areas)for(unsigned j=0;j<7;j++)area[j]=j+1;
    exchange.scatterAddTargets(areas,cellAreas);
    forAll(ids,i){
        label expected=0;for(label x=0;x<nx;x++)for(label y=0;y<ny;y++){const label d=distance(x,y,ids[i]%nx,ids[i]/nx);expected+=d>=1&&d<=2;}
        for(unsigned j=0;j<7;j++)if(cellAreas[i][j]!=expected*(j+1))return 3;
    }
    auto contribution=[](label id){return id%3==0?1e16:id%3==1?.1*(id+1):1e-9*(id+1);};
    auto moments=[&](label gx,label gy,label rx,label ry){
        const double weight=contribution(gx+nx*gy);label dx=rx-gx,dy=ry-gy;
        if(periodic[0]&&std::abs(dx)>nx/2)dx+=dx>0?-nx:nx;
        if(periodic[1]&&std::abs(dy)>ny/2)dy+=dy>0?-ny:ny;
        FixedList<scalar,7> result;for(unsigned j=0;j<4;j++)result[j]=(j+1)*weight;
        result[4]=dx*weight;result[5]=dy*weight;result[6]=-weight;return result;
    };
    for(std::size_t i=0;i<links.size();i++)areas[i]=moments(links[i].gx,links[i].gy,links[i].rx,links[i].ry);
    std::fill(cellAreas.begin(),cellAreas.end(),FixedList<scalar,7>(Zero));exchange.scatterAddTargets(areas,cellAreas);
    forAll(ids,i){
        FixedList<scalar,7> expected(Zero);const label rx=ids[i]%nx,ry=ids[i]/nx;
        for(label gx=0;gx<nx;gx++)for(label gy=0;gy<ny;gy++){
            const label d=distance(gx,gy,rx,ry);if(d<1||d>2)continue;
            const auto moment=moments(gx,gy,rx,ry);for(unsigned j=0;j<7;j++)expected[j]+=moment[j];
        }
        if(cellAreas[i]!=expected)return 12;
    }
    std::ofstream areaOutput((runTime.path()/"link_areas.csv").c_str());areaOutput<<std::setprecision(17)<<"id,a0,a1,a2,a3,a4,a5,a6\n";
    forAll(ids,i){areaOutput<<ids[i];for(double value:cellAreas[i])areaOutput<<','<<value;areaOutput<<'\n';}
    auto transfer=[](label,Packet& source,Packet& target){
        if(!source[30]||target[30])return;
        const double amount=std::min({.6,source[0],target[15]});
        for(int j=0;j<7;j++){source[j]-=(j+1)*amount;target[7+j]+=(j+1)*amount;}
        source[14]-=amount;target[14]+=amount;target[15]-=amount;source[24]+=amount;
    };
    boolList active(links.size(),true);for(int pass=0;pass<3;pass++)exchange.execute(cells,active,transfer);
    double conservation=0,minimum=GREAT;
    for(unsigned j=0;j<8;j++){
        double before=0,after=0;
        for(label i=0;i<label(cells.size());i++){before+=j<7?initial[i][j]+initial[i][7+j]:initial[i][14];after+=j<7?cells[i][j]+cells[i][7+j]:cells[i][14];minimum=min(minimum,j<7?min(cells[i][j],cells[i][7+j]):cells[i][14]);}
        reduce(before,sumOp<double>());reduce(after,sumOp<double>());conservation=max(conservation,std::abs(after-before)/max(std::abs(before),1.));
    }
    reduce(minimum,minOp<double>());if(conservation>5e-15||minimum<-1e-14)return 4;
    std::ofstream output((runTime.path()/"link_cells.csv").c_str());output<<std::setprecision(17)<<"id";for(int j=0;j<31;j++)output<<",p"<<j;output<<'\n';
    for(const auto& packet:cells){output<<packet[23];for(double value:packet)output<<','<<value;output<<'\n';}
    // A two-hop edge reaches rank 2 directly on a three-column/three-rank
    // grid. A failure is injected only at its source process.
    cells=initial;bool failed=false;
    try{exchange.execute(cells,active,[&](label link,Packet& a,Packet& b){
        if(links[link].gx==0&&links[link].gy==0&&links[link].rx==min(nx-1,label(2))&&links[link].ry==0){a[0]=-123;b[0]=-456;throw std::runtime_error("intentional two-hop callback failure");}
        transfer(link,a,b);
    });}catch(const std::runtime_error&){failed=true;}
    int allFailed=failed;reduce(allFailed,minOp<int>());if(!allFailed)return 5;
    for(const auto& packet:cells)if(packet[0]<-1e-14)return 6;
    // Sparse and empty local link lists still participate in requests from
    // other ranks. Exercise these with no local sources on most partitions.
    std::vector<SurfaceLink> sparse;for(const auto& link:links)if(ids[link.source]==0)sparse.push_back(link);
    SurfaceLinkExchange sparseExchange(mesh,sparse,spacing,spacing,periodic);cells=initial;
    boolList sparseActive(sparse.size(),true);sparseExchange.execute(cells,sparseActive,transfer);
    const auto sparseTargets=sparseExchange.gatherTargets(cells);for(std::size_t i=0;i<sparse.size();i++)if(sparseTargets[i][23]!=sparse[i].rx+nx*sparse[i].ry)return 7;
    SurfaceLinkExchange empty(mesh,{},spacing,spacing,periodic);const auto unchanged=cells;
    empty.execute(cells,boolList(),transfer);if(cells!=unchanged||!empty.gatherTargets(cells).empty())return 8;
    for(int mode=0;mode<3;mode++){
        auto invalid=links;
        if(!invalid.empty()&&Pstream::master()){
            if(mode==0)invalid.push_back(invalid.front());
            if(mode==1)invalid.front().receiverCell=labelMax;
            if(mode==2)invalid.front().rx=-1;
        }
        bool rejected=false;try{SurfaceLinkExchange bad(mesh,invalid,spacing,spacing,periodic);}catch(const std::invalid_argument&){rejected=true;}
        int allRejected=rejected;reduce(allRejected,minOp<int>());if(!allRejected)return 9+mode;
    }
    if(Pstream::master())Info<<"SURFACE_LINK_CHECK conservation="<<conservation<<" minimum="<<minimum<<" collective_failure="<<allFailed<<endl;
    return 0;
}
