// SPDX-License-Identifier: GPL-3.0-or-later
#include "FaceReactionExchange.H"
#include <fstream>
#include <iomanip>
#include <memory>

int main(int argc,char** argv){
    using namespace Foam;
    argList args(argc,argv);if(!args.checkRootCase())return 1;
    Time runTime(Time::controlDictName,args);
    fvMesh mesh(IOobject(fvMesh::defaultRegion,runTime.timeName(),runTime,IOobject::MUST_READ));
    IOdictionary config(IOobject("exchangeTest",runTime.constant(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE));
    const double spacing=readScalar(config.lookup("spacing"));
    const label nx=readLabel(config.lookup("nx")),ny=readLabel(config.lookup("ny"));
    const bool expectedUnsupported=config.lookupOrDefault<bool>("expectedUnsupported",false);
    PhaseStencil stencil(mesh,spacing,spacing);
    std::unique_ptr<saltflow::FaceReactionExchange> exchange;
    bool unsupported=false;
    try{exchange=std::make_unique<saltflow::FaceReactionExchange>(mesh,stencil,spacing,spacing);}
    catch(const std::invalid_argument&){unsupported=true;}
    if(expectedUnsupported)return unsupported?0:2;
    if(unsupported)return 3;
    using Packet=FixedList<scalar,31>;
    std::vector<Packet> cells(mesh.nCells(),Packet(Zero));
    forAll(mesh.C(),i){
        const label x=std::llround(mesh.C()[i].x()/spacing-.5),y=std::llround(mesh.C()[i].y()/spacing-.5),id=x+nx*y;
        // Every fourth cell is a finite donor; multiple adjacent sinks must
        // share its one inventory, including across processor/cyclic faces.
        const bool donor=id%4==0;
        cells[i][30]=donor;cells[i][23]=id;cells[i][22]=1234.56789;
        cells[i][15]=.9;
        if(donor){for(int j=0;j<7;j++)cells[i][j]=(j+1)*(.4+.1*(id%5));cells[i][14]=cells[i][0];}
    }
    boolList active(mesh.nFaces(),false);
    forAll(mesh.neighbour(),f)active[f]=true;
    forAll(mesh.boundaryMesh(),patch)if(mesh.boundaryMesh()[patch].coupled()){
        const auto& p=mesh.boundaryMesh()[patch];forAll(p,f)active[p.start()+f]=!(exchange->direction(p.start()+f)&1);
    }
    std::vector<FixedList<scalar,2>> faceData(exchange->coupledFaces().size());
    for(std::size_t k=0;k<faceData.size();k++){const label face=exchange->coupledFaces()[k];faceData[k][0]=cells[exchange->cell(face)][23];faceData[k][1]=faceData[k][0]+.25;}
    exchange->swapCoupled(faceData);
    label failureKey=labelMax;
    for(std::size_t k=0;k<faceData.size();k++){
        const label face=exchange->coupledFaces()[k],id=cells[exchange->cell(face)][23];
        const int direction=exchange->direction(face),sign=(direction&1)?1:-1;
        const label other=direction<2?(id/nx)*nx+(id%nx+sign+nx)%nx:((id/nx+sign+ny)%ny)*nx+id%nx;
        if(faceData[k][0]!=other||faceData[k][1]!=other+.25)return 8;
        const label patch=mesh.boundaryMesh().whichPatch(face);
        if(isA<processorPolyPatch>(mesh.boundaryMesh()[patch]))failureKey=min(failureKey,min(id,other));
    }
    reduce(failureKey,minOp<label>());if(!Pstream::parRun())failureKey=0;
    auto transfer=[](label,bool,Packet& a,Packet& b){
        if(a[30]==b[30])return;
        auto& donor=a[30]?a:b;auto& sink=a[30]?b:a;
        const double amount=std::min({.6,donor[0],sink[15]});
        for(int j=0;j<7;j++){donor[j]-=(j+1)*amount;sink[7+j]+=(j+1)*amount;}
        donor[14]-=amount;sink[14]+=amount;sink[15]-=amount;donor[24]+=amount;
    };
    const auto initial=cells;
    for(int pass=0;pass<3;pass++)exchange->execute(cells,active,transfer);
    std::array<double,8> before{},after{};double minimum=GREAT;
    for(label i=0;i<label(cells.size());i++){
        for(int j=0;j<7;j++){before[j]+=initial[i][j]+initial[i][7+j];after[j]+=cells[i][j]+cells[i][7+j];minimum=min(minimum,min(cells[i][j],cells[i][7+j]));}
        before[7]+=initial[i][14];after[7]+=cells[i][14];minimum=min(minimum,min(cells[i][14],cells[i][15]));
        if(cells[i][22]!=initial[i][22]||cells[i][23]!=initial[i][23]||cells[i][30]!=initial[i][30])return 4;
    }
    double error=0;for(int j=0;j<8;j++){reduce(before[j],sumOp<double>());reduce(after[j],sumOp<double>());error=max(error,std::abs(before[j]-after[j])/max(std::abs(before[j]),1.));}reduce(minimum,minOp<double>());
    if(error>5e-15||minimum<-1e-14)return 5;
    std::ofstream output((runTime.path()/"exchange_cells.csv").c_str());output<<std::setprecision(17)<<"id";
    for(int j=0;j<31;j++)output<<",p"<<j;output<<'\n';
    for(const auto& packet:cells){output<<packet[23];for(double value:packet)output<<','<<value;output<<'\n';}
    // Empty and spatially sparse work still follow one collective schedule.
    // In the sparse sweep only the rank containing cell zero enables faces.
    cells=initial;boolList sparse(mesh.nFaces(),false);label calls=0;
    exchange->execute(cells,sparse,[&](label,bool,Packet&,Packet&){calls++;});
    if(calls||cells!=initial)return 9;
    forAll(mesh.neighbour(),f)sparse[f]=cells[mesh.owner()[f]][23]==0||cells[mesh.neighbour()[f]][23]==0;
    for(label f:exchange->coupledFaces())sparse[f]=cells[exchange->cell(f)][23]==0;
    exchange->execute(cells,sparse,[&](label f,bool boundary,Packet& a,Packet& b){calls++;transfer(f,boundary,a,b);});
    reduce(calls,sumOp<label>());if(!calls)return 10;
    double sparseTotal=0,initialTotal=0;
    for(label i=0;i<label(cells.size());i++){sparseTotal+=cells[i][0]+cells[i][7];initialTotal+=initial[i][0]+initial[i][7];}
    reduce(sparseTotal,sumOp<double>());reduce(initialTotal,sumOp<double>());
    if(std::abs(sparseTotal-initialTotal)>1e-14*max(initialTotal,1.))return 11;
    // Throw only on the unique computing side of one cross-rank face. All
    // peers must still finish the halo return before receiving the failure.
    cells=initial;bool failed=false;
    try{
        exchange->execute(cells,active,[&](label,bool boundary,Packet& a,Packet& b){
            if(min(label(a[23]),label(b[23]))==failureKey&&(boundary||!Pstream::parRun())){
                a[0]=-123.;b[0]=-456.;throw std::runtime_error("intentional transactional callback failure");
            }
            transfer(0,false,a,b);
        });
    }catch(const std::runtime_error&){failed=true;}
    int allFailed=failed;reduce(allFailed,minOp<int>());if(!allFailed)return 6;
    for(const auto& packet:cells)if(packet[0]<-1e-14)return 7;
    if(Pstream::master())Info<<"FACE_EXCHANGE_CHECK cells="<<nx*ny<<" conservation="<<error<<" minimum="<<minimum<<" collective_failure="<<allFailed<<endl;
    return 0;
}
