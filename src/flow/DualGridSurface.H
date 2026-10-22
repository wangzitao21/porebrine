// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "LocalWaterReaction.hpp"
#include "LocalChemistry.H"
#include "SurfaceReaction.hpp"
#include <chrono>

namespace saltflow {
class DualGridSurface {
    const Foam::fvMesh& mesh;
    std::vector<int> index,foamIndex;
    bool unblocked=false;
    double geometryInterval=.01,minimumLiquid=.5,diffusivity=0;
    mutable bool pressureAnchorsReady=false;
    mutable Foam::labelList pressureAnchors;
    std::uint64_t seed=2024;
    std::vector<saltgeom::NucleationLaw> laws;
    using Areas=std::array<double,3>;
    saltdual::ReactionMap map;
    std::vector<saltdual::WaterTask> tasks;
    std::vector<std::array<unsigned char,4>> locallyTransported;
    void allocate(){map.rebuild(geometry,unblocked,minimumLiquid);tasks=saltdual::waterTasks(geometry,map);pressureAnchorsReady=false;}
    void topology()const{
        if(!geometry.waterInitialized)return;
        for(int p=0;p<geometry.nx*geometry.ny;p++)if(std::abs(geometry.closedFraction[p]-geometry.storedClosedFraction[p])>1e-9)
            throw std::runtime_error("Closed pore geometry and retained liquid inventory disagree");
    }
public:
    bool enabled=false;
    saltdual::Geometry geometry;
    std::vector<std::array<double,2>> reactionVectors;
    double calls=0,accepted=0,retries=0,geometrySeconds=0,maxVolumeError=0,maxLocalLiquidRegions=0,maxIsolatedLiquidFraction=0;
    double maximumReactors=0,minimumReceiverPore=1;
    double waterExchangeSeconds=0,waterRemapSeconds=0,waterExchangeIterations=0,maxWaterResidual=0,maxWaterPools=0,contactEvents=0,maxWaterCluster=0,maxLocalReactionResidual=0,maxLocalStateError=0;
    DualGridSurface(const Foam::dictionary& cfg,const Foam::fvMesh& m,double dx,double dy,double epsilon):mesh(m){
        enabled=cfg.lookupOrDefault<bool>("dualMineralGeometry",false);if(!enabled)return;
        if(Foam::Pstream::parRun())throw std::invalid_argument("Dual polygon geometry currently requires one MPI rank; chemistry can use OpenMP threads");
        if(cfg.lookupOrDefault<Foam::word>("mineralInterfaceModel","phaseVector")!="phaseVector"||cfg.found("initialSurfaceCoating"))throw std::invalid_argument("Dual geometry cannot be combined with the population coating backend");
        if(cfg.lookupOrDefault<bool>("adaptiveTimeStep",false))throw std::invalid_argument("Dual geometry uses event rejection, not step-doubling error estimation");
        unblocked=cfg.lookupOrDefault<bool>("dualIgnoreShielding",false);
        geometryInterval=cfg.lookupOrDefault<Foam::scalar>("dualGeometryInterval",.01);
        if(!(geometryInterval>0))throw std::invalid_argument("Invalid dual geometry interval");
        const auto path=m.time().value()==0?m.time().constant()/"dualGeometry.bin":m.time().timePath()/"dualGeometry.bin";
        // Covered-area totals always remain available. Individual buried
        // faces are used only by the shielding-disabled comparison mode.
        geometry=saltdual::Geometry::read(path.c_str(),unblocked);
        diffusivity=Foam::dimensionedScalar("diffusivity",Foam::dimViscosity,cfg).value();
        const double neck=cfg.lookupOrDefault<Foam::scalar>("dualWaterNeckFraction",geometry.waterNeckFraction);
        if(!(neck>=0&&neck<=1))throw std::invalid_argument("dualWaterNeckFraction must lie in [0,1]");
        if(neck!=geometry.waterNeckFraction){if(geometry.waterInitialized)throw std::invalid_argument("Restart requires the saved local water neck threshold");geometry.waterNeckFraction=neck;geometry.rebuild();}
        minimumLiquid=cfg.lookupOrDefault<Foam::scalar>("dualMinimumLiquidFraction",.5);
        if(!(minimumLiquid>0&&minimumLiquid<=1))throw std::invalid_argument("dualMinimumLiquidFraction must lie in (0,1]");
        if(minimumLiquid!=geometry.waterCouplingFraction){if(geometry.waterInitialized)throw std::invalid_argument("Restart requires the saved water coupling threshold");geometry.waterCouplingFraction=minimumLiquid;geometry.rebuild();}
        locallyTransported.resize(geometry.nx*geometry.ny);
        seed=cfg.lookupOrDefault<Foam::label>("nucleationSeed",2024);
        if(cfg.found("heterogeneousNucleation"))for(const auto& entry:cfg.subDict("heterogeneousNucleation")){
            const auto& d=entry.dict();saltgeom::NucleationLaw law;
            law.product=saltchem::mineralIndex(Foam::word(d.lookup("product")).c_str());law.substrate=saltchem::mineralIndex(Foam::word(d.lookup("substrate")).c_str());
            if(d.lookupOrDefault<Foam::word>("nucleationModel","effectiveBarrier")!="effectiveBarrier"||d.found("interfacialEnergyJm2"))throw std::invalid_argument("Dual planar nuclei currently require the effectiveBarrier law, not a 3D Kelvin/CNT population");
            law.prefactor=Foam::readScalar(d.lookup("prefactor"));law.barrier=Foam::readScalar(d.lookup("barrier"));
            law.contactAngle=d.lookupOrDefault<Foam::scalar>("contactAngle",60);law.radius=d.lookupOrDefault<Foam::scalar>("embryoRadius",1e-7);law.validate();
            if(law.product>=3||law.substrate>=3||law.radius>.2*geometry.h())throw std::invalid_argument("Dual nucleation requires the three supported minerals and embryo radius <=0.2h");
            for(const auto& old:laws)if(old.product==law.product&&old.substrate==law.substrate)throw std::invalid_argument("Duplicate dual nucleation pair");
            laws.push_back(law);
        }
        if(std::abs(dx/geometry.spacing-1)>1e-12||std::abs(dy/dx-1)>1e-12||epsilon!=geometry.epsilon||m.nCells()!=geometry.nx*geometry.ny)throw std::invalid_argument("Main mesh and fine geometry do not match");
        index.resize(m.nCells());foamIndex.resize(m.nCells(),-1);reactionVectors.resize(m.nCells());
        forAll(m.C(),i){int x=int(std::llround(m.C()[i].x()/dx-.5)),y=int(std::llround(m.C()[i].y()/dy-.5)),p=y*geometry.nx+x;
            if(x<0||y<0||x>=geometry.nx||y>=geometry.ny||foamIndex[p]>=0||std::abs(m.V()[i]/geometry.coarseVolume()-1)>1e-10)throw std::invalid_argument("Dual geometry requires a regular, origin-aligned 2D main mesh");
            index[i]=p;foamIndex[p]=i;
        }
        topology();allocate();
        maxLocalLiquidRegions=*std::max_element(geometry.liquidComponents.begin(),geometry.liquidComponents.end());
        Foam::Info<<"DUAL_GEOMETRY main_cells="<<m.nCells()<<" fine_cells="<<geometry.cells.size()<<" H_over_h="<<geometry.ratio<<" ignore_shielding="<<unblocked<<Foam::endl;
    }
    void check(const Foam::PtrList<Foam::volScalarField>& solid,const Foam::volScalarField& inert){
        if(!enabled)return;
        forAll(inert,i){if(inert[i]!=0)throw std::invalid_argument("Dual geometry prototype supports Halite, Sylvite and Carnallite only");
            for(int j=0;j<7;j++){double expected=j<3?geometry.fractions[index[i]][j]:0;double error=std::abs(solid[j][i]-expected);maxVolumeError=std::max(maxVolumeError,error);if(error>1e-9)throw std::runtime_error("Main/fine geometry inventory mismatch");}
        }
    }
    void configure(int i,const saltchem::Inventory& cell,saltchem::RateParameters& rates)const{
        if(!enabled)return;int p=index[i];rates.geometry.separateMineralAreas=true;rates.area={};rates.geometry.growthArea={};rates.geometry.surfaceExponent={};rates.geometry.lockedMoles={};
        rates.geometry.inertFraction=geometry.storedClosedFraction[p];
        for(int m=0;m<3;m++){rates.area[m]=map.dissolution[p][m]/mesh.V()[i];rates.geometry.growthArea[m]=map.growth[p][m]/mesh.V()[i];rates.referenceMoles[m]=std::max(cell.mineral[m],1e-30);}
    }
    void initialize(std::vector<saltchem::Inventory>& cells,const saltchem::MineralVector& volumes,Foam::volScalarField& epsilon){
        if(!enabled||geometry.waterInitialized)return;
        std::vector<saltchem::Inventory> ordered(cells.size());for(int p=0;p<int(cells.size());p++)ordered[p]=cells[foamIndex[p]];
        saltdual::initializeWater(geometry,ordered,volumes);
        for(int p=0;p<int(cells.size());p++){int i=foamIndex[p];epsilon[i]*=saltchem::aqueousMass(ordered[p])/saltchem::aqueousMass(cells[i]);cells[i]=ordered[p];}
        epsilon.correctBoundaryConditions();
    }
    void constrainDiffusion(Foam::surfaceScalarField& field,bool flow=false)const{
        if(!enabled)return;
        // Geometry supplies the open water cross section. It replaces the
        // volume-averaged epsilon^2 penalty; a cut cell's tiny stored water
        // fraction is not an additional physical diffusion resistance.
        forAll(mesh.neighbour(),f){int p=index[mesh.owner()[f]],q=index[mesh.neighbour()[f]];int direction=q-p==1?1:q-p==-1?0:q>p?3:2;
            field[f]*=(!flow&&locallyTransported[p][direction])?0:std::clamp(geometry.openings[p][direction],0.,1.);
        }
        forAll(mesh.boundary(),patch)if(mesh.boundary()[patch].type()!="empty")forAll(mesh.boundary()[patch],f){
            const auto sf=mesh.Sf().boundaryField()[patch][f];int direction=std::abs(sf.x())>std::abs(sf.y())?(sf.x()>0?1:0):(sf.y()>0?3:2);
            const int p=index[mesh.boundary()[patch].faceCells()[f]];field.boundaryFieldRef()[patch][f]*=std::clamp(geometry.openings[p][direction],0.,1.);
        }
    }
    void anchorPressure(Foam::fvScalarMatrix& equation,const Foam::volScalarField& pressure)const{
        if(!enabled)return;
        if(!pressureAnchorsReady){
            std::vector<unsigned char> seen(geometry.nx*geometry.ny),fixed(seen.size());std::vector<Foam::label> anchors;
            forAll(mesh.boundary(),patch)if(pressure.boundaryField()[patch].fixesValue())forAll(mesh.boundary()[patch],f){
                int p=index[mesh.boundary()[patch].faceCells()[f]];const auto normal=mesh.Sf().boundaryField()[patch][f];int d=std::abs(normal.x())>std::abs(normal.y())?(normal.x()>0?1:0):(normal.y()>0?3:2);if(geometry.openings[p][d]>0)fixed[p]=1;
            }
            for(int first=0;first<int(seen.size());first++)if(!seen[first]){
                std::vector<int> queue{first};seen[first]=1;bool anchored=false;
                for(std::size_t k=0;k<queue.size();k++){int p=queue[k];anchored|=fixed[p];for(int d=0;d<4;d++){
                    if(geometry.openings[p][d]<=1e-12||(d==0&&p%geometry.nx==0)||(d==1&&p%geometry.nx==geometry.nx-1)||(d==2&&p<geometry.nx)||(d==3&&p>=(geometry.ny-1)*geometry.nx))continue;
                    int q=p+(d==0?-1:d==1?1:d==2?-geometry.nx:geometry.nx);if(!seen[q]){seen[q]=1;queue.push_back(q);}
                }}
                if(!anchored)anchors.push_back(foamIndex[first]);
            }
            pressureAnchors.setSize(anchors.size());for(std::size_t k=0;k<anchors.size();k++)pressureAnchors[k]=anchors[k];pressureAnchorsReady=true;
        }
        double scale=0;int count=0;for(double d:equation.diag())if(std::abs(d)>1e-30){scale+=std::abs(d);count++;}scale=count?scale/count:1;
        Foam::scalarField values(pressureAnchors.size());forAll(pressureAnchors,k)values[k]=pressure[pressureAnchors[k]];
        // Preserve the pressure matrix's physical scale. Unit diagonals on
        // solid rows would dominate the solver's normalised residual.
        forAll(pressureAnchors,k)if(std::abs(equation.diag()[pressureAnchors[k]])<1e-30)equation.diag()[pressureAnchors[k]]=-scale;
        equation.setValues(pressureAnchors,values);
    }
    struct Stats {double accepted=0,evaluations=0,restarts=0,fallbacks=0,jacobians=0,jacobianFallbacks=0;};
    template<class Pressure>
    Stats advance(std::vector<saltchem::Inventory>& cells,const saltchem::RateParameters& base,const saltchem::MineralVector& volume,
                  LocalChemistry& chemistry,double temperature,Pressure pressure,double interval,const saltchem::IntegrationOptions& options){
        Stats total;double remaining=interval,step=std::min(interval,geometryInterval);int rejection=0;
        std::vector<saltchem::Inventory> ordered(cells.size());for(int p=0;p<int(cells.size());p++)ordered[p]=cells[foamIndex[p]];
        saltdual::initializeWater(geometry,ordered,volume);
        std::fill(locallyTransported.begin(),locallyTransported.end(),std::array<unsigned char,4>{});
        std::fill(reactionVectors.begin(),reactionVectors.end(),std::array<double,2>{});
        int attempts=0;
        bool capturedChemistry=false;
        while(remaining>interval*1e-12){
            if(++attempts>2048)throw std::runtime_error("Dual geometry exceeded its per-step event budget");
            step=std::min(step,remaining);const auto start=ordered;const auto original=geometry;
            std::vector<Areas> target;
            std::vector<saltdual::SurfaceBudget> budgets;
            auto saveFailure=[&](){if(target.empty())return;const auto output=mesh.time().globalPath()/"output";original.write((output/"dual_failure_start.bin").c_str());
                std::ofstream f((output/"dual_failure_target.f64").c_str(),std::ios::binary);for(const auto& a:target)f.write(reinterpret_cast<const char*>(a.data()),3*sizeof(double));
                std::ofstream detail((output/"dual_failure_budget.bin").c_str(),std::ios::binary);for(const auto& a:target)detail.write(reinterpret_cast<const char*>(a.data()),sizeof(a));std::uint32_t count=budgets.size();detail.write(reinterpret_cast<const char*>(&count),4);
                for(const auto& b:budgets){std::int32_t values[]{b.parent,b.mineral,b.liquidRegion,int(b.surfaces.size())};detail.write(reinterpret_cast<const char*>(values),sizeof(values));detail.write(reinterpret_cast<const char*>(&b.change),8);detail.write(reinterpret_cast<const char*>(b.surfaces.data()),4*b.surfaces.size());}
            };
            auto record=[&](const saltchem::IntegrationStats& s){total.accepted+=s.accepted;total.evaluations+=s.rateEvaluations;total.restarts+=s.historyRestarts;total.fallbacks+=s.integratorFallbacks;total.jacobians+=s.analyticJacobians;total.jacobianFallbacks+=s.jacobianFallbacks;calls++;};
            try{
                saltdual::applyMineralDebt(geometry,ordered,volume);
                std::vector<std::vector<saltdual::RegionReaction>> batch(tasks.size());std::vector<std::string> errors(tasks.size());std::vector<unsigned char> retryable(tasks.size());
                auto reactionTick=std::chrono::steady_clock::now();
                #pragma omp parallel for schedule(dynamic,1)
                for(int task=0;task<int(tasks.size());task++){
                    try{batch[task]=saltdual::reactWaterTask(tasks[task],geometry,map,ordered,base,volume,diffusivity,step,options,[&](int p)->ChemicalWorker&{return chemistry.at(temperature,pressure(foamIndex[p]));});}
                    catch(const saltdual::StepFailure& e){errors[task]=e.what();retryable[task]=1;}
                    catch(const std::exception& e){errors[task]=e.what();}
                }
                waterExchangeSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-reactionTick).count();
                if(std::getenv("SALT_DUAL_TRACE"))Foam::Info<<"DUAL_REACTION interval="<<step<<" attempt="<<attempts<<" tasks="<<tasks.size()<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-reactionTick).count()<<Foam::endl;
                std::vector<saltdual::RegionReaction> results;
                for(int task=0;task<int(tasks.size());task++){
                    if(!errors[task].empty()){
                        if(std::getenv("SALT_DUAL_CAPTURE")&&(!capturedChemistry||rejection>=39)){
                            const auto output=mesh.time().globalPath()/"output";geometry.write((output/"dual_chemistry_start.bin").c_str());std::ofstream f((output/"dual_chemistry_state.bin").c_str(),std::ios::binary);
                            std::uint32_t header[]{0x57415443,1,std::uint32_t(ordered.size()),std::uint32_t(task),std::uint32_t(options.maximumAttempts),std::uint32_t(options.method),std::uint32_t(options.analyticJacobian),0};f.write(reinterpret_cast<const char*>(header),sizeof(header));
                            double values[]{step,temperature,diffusivity,base.geometry.epsilonMin,base.geometry.closureWidth,base.geometry.directionWidth,options.relativeTolerance,options.absoluteTolerance};f.write(reinterpret_cast<const char*>(values),sizeof(values));
                            for(const auto* vector:{&base.neutral,&base.acid,&base.co2,&base.acidOrder,&base.dissolutionExponent,&base.affinityExponent,&base.geometry.weights})f.write(reinterpret_cast<const char*>(vector->data()),sizeof(*vector));
                            for(int p=0;p<int(ordered.size());p++){double pa=pressure(foamIndex[p]);f.write(reinterpret_cast<const char*>(&pa),8);auto stock=saltdual::waterStock(ordered[p]);f.write(reinterpret_cast<const char*>(stock.data()),sizeof(stock));f.write(reinterpret_cast<const char*>(ordered[p].mineral.data()),sizeof(ordered[p].mineral));}capturedChemistry=true;
                        }
                        if(retryable[task])throw saltdual::StepFailure(errors[task]);throw std::runtime_error(errors[task]);
                    }
                    maxWaterCluster=std::max(maxWaterCluster,double(tasks[task].regions.size()));
                    for(auto& result:batch[task]){if(result.stats.rateEvaluations>0)record(result.stats);waterExchangeIterations+=result.stats.newtonIterations;maxWaterResidual=std::max(maxWaterResidual,result.exchangeResidual);maxLocalReactionResidual=std::max(maxLocalReactionResidual,result.nonlinearResidualRatio);maxLocalStateError=std::max(maxLocalStateError,result.nonlinearStateErrorRatio);results.push_back(std::move(result));}
                }
                budgets=saltdual::commitRegions(geometry,ordered,results,volume);
                target.resize(geometry.nx*geometry.ny);for(int p=0;p<int(ordered.size());p++)for(int m=0;m<3;m++)target[p][m]=ordered[p].mineral[m]*volume[m];
                std::vector<saltdual::Point> vectors(ordered.size());
                auto tick=std::chrono::steady_clock::now();
                try{geometry.advance(target,unblocked,.2,&budgets);}catch(...){geometrySeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-tick).count();throw;}
                geometrySeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-tick).count();
                vectors=geometry.sweptVolumeVectors;
                maxVolumeError=std::max(maxVolumeError,geometry.chemicalVolumeError);
                for(int p=0;p<int(ordered.size());p++)for(int m=0;m<3;m++)ordered[p].mineral[m]=geometry.fractions[p][m]/volume[m];
                tick=std::chrono::steady_clock::now();if(geometry.lastMovementChanged)saltdual::remapWater(original,geometry,ordered,&geometry.pockets);waterRemapSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-tick).count();
                auto births=saltdual::nucleateSurface(geometry,map,ordered,laws,volume,step,temperature,seed,[&](int p)->ChemicalWorker&{return chemistry.at(temperature,pressure(foamIndex[p]));});

                for(int p=0;p<int(target.size());p++)for(int m=0;m<3;m++){
                    double error=std::abs(ordered[p].mineral[m]*volume[m]-geometry.fractions[p][m]);maxVolumeError=std::max(maxVolumeError,error);if(error>1e-9)throw saltdual::StepFailure("Nucleation/coarse inventory mismatch");
                }
                topology();
                for(int p=0;p<int(ordered.size());p++){auto v=vectors[p]+births.volumeVectors[p];reactionVectors[foamIndex[p]][0]+=v.x/interval;reactionVectors[foamIndex[p]][1]+=v.y/interval;}
                if(geometry.lastMovementChanged||std::accumulate(births.count.begin(),births.count.end(),std::uint64_t(0))>0)allocate();maxLocalLiquidRegions=std::max(maxLocalLiquidRegions,double(*std::max_element(geometry.liquidComponents.begin(),geometry.liquidComponents.end())));
                maxIsolatedLiquidFraction=std::max(maxIsolatedLiquidFraction,geometry.isolatedLiquidArea/geometry.cells.size());maxWaterPools=std::max(maxWaterPools,double(geometry.pockets.size()));
                maximumReactors=std::max(maximumReactors,double(map.active.size()));for(int r:map.active)minimumReceiverPore=std::min(minimumReceiverPore,geometry.waterRegions[r].fraction);
                for(int p=0;p<geometry.nx*geometry.ny;p++)for(int d=0;d<4;d++)if(original.openings[p][d]>original.diffusiveOpenings[p][d]+1e-12)locallyTransported[p][d]=1;
                remaining-=step;accepted++;rejection=0;step=std::min(geometryInterval,step*1.5);
            }catch(const saltdual::StepFailure& error){
                if(std::getenv("SALT_DUAL_TRACE"))Foam::Info<<"DUAL_RETRY interval="<<step<<" reason="<<error.what()<<Foam::endl;
                if(std::getenv("SALT_DUAL_CAPTURE")&&rejection==0)saveFailure();
                double factor=std::clamp(error.scale,1e-6,.999999);if(factor!=.5)contactEvents++;
                // The first geometric estimate targets the contact event.
                // A stiff local reaction is not proportional to interval
                // length: repeated near-unit estimates can otherwise use
                // the entire retry budget without approaching admissibility.
                if(rejection>0)factor=std::min(factor,.5);
                if(rejection>=40||step*factor<1e-12)saveFailure();geometry=original;ordered=start;retries++;step*=factor;
                if(++rejection>40||step<1e-12)throw;
            }catch(...){saveFailure();throw;}
        }
        for(int p=0;p<int(ordered.size());p++)cells[foamIndex[p]]=ordered[p];return total;
    }
    void restore(const saltdual::Geometry& g){if(enabled){geometry=g;allocate();}}
    void checkpoint()const{if(enabled)geometry.write((mesh.time().timePath()/"dualGeometry.bin").c_str());}
    void save(const Foam::fileName& output)const{
        if(!enabled)return;
        Foam::mkDir(output/"dual_geometry");geometry.write((output/"dual_geometry"/(mesh.time().timeName()+".bin")).c_str());
        std::vector<double> isolated(geometry.nx*geometry.ny);for(const auto& r:geometry.waterRegions)if(!r.connected)isolated[r.parent]+=r.fraction;
        std::ofstream f((output/"dual_geometry"/(mesh.time().timeName()+".csv")).c_str());f<<std::setprecision(17)<<"cell,liquid_regions,liquid_receiver,nuclei_Halite,nuclei_Sylvite,nuclei_Carnallite,closed_water_fraction,retained_water_fraction,mobile_water_fraction,isolated_water_fraction";
        for(int m=0;m<3;m++)for(auto suffix:{"_fraction","_area_m2","_covered_m2","_reaction_area_m2"})f<<','<<saltchem::mineralNames[m]<<suffix;f<<'\n';
        for(int p=0;p<geometry.nx*geometry.ny;p++){std::uint64_t h=0,s=0,c=0;for(int k=0;k<3;k++){h+=geometry.nucleationEvents[p][k];s+=geometry.nucleationEvents[p][3+k];c+=geometry.nucleationEvents[p][6+k];}
            f<<p<<','<<geometry.liquidComponents[p]<<','<<map.receiver[p]<<','<<h<<','<<s<<','<<c<<','<<geometry.storedClosedFraction[p]<<','<<geometry.storedClosedFraction[p]<<','<<map.pore[p]<<','<<isolated[p];for(int m=0;m<3;m++)f<<','<<geometry.fractions[p][m]<<','<<geometry.dissolution[p][m]<<','<<geometry.covered[p][m]<<','<<map.dissolution[p][m];f<<'\n';}
        std::ofstream water((output/"dual_geometry"/(mesh.time().timeName()+".water.csv")).c_str());water<<std::setprecision(17)<<"pool,cell,liquid_fraction,Na,K,Ca,Mg,Cl,S,C,water_equivalent\n";
        for(std::size_t k=0;k<geometry.pockets.size();k++){const auto& p=geometry.pockets[k];water<<k<<','<<p.parent<<','<<p.fraction;for(double a:p.inventory)water<<','<<a;water<<'\n';}
    }
    template<class Pressure>
    void saveReactiveWater(const Foam::fileName& output,const std::vector<saltchem::Inventory>& cells,LocalChemistry& chemistry,double temperature,Pressure pressure)const{
        if(!enabled)return;
        std::ofstream f((output/"dual_geometry"/(mesh.time().timeName()+".reactors.csv")).c_str());
        f<<std::setprecision(17)<<"region,cell,retained,connected,liquid_fraction,A_Halite,A_Sylvite,A_Carnallite,mNa,mK,mMg,mCl,SIHalite,SISylvite,SICarnallite,waterActivity,log_aK,log_aMg,log_aCl\n";
        for(int r:map.active){const auto& region=geometry.waterRegions[r];int i=foamIndex[region.parent];auto state=cells[i];
            if(geometry.waterInitialized&&!region.bulk)saltdual::setWater(state,geometry.pockets.at(geometry.pocketRegion.at(r)).inventory);
            if(state.waterEquivalent<=0)continue;
            double water;auto chemical=chemistry.at(temperature,pressure(i)).speciation.conserved(state,water,nullptr,false);auto totals=saltchem::componentTotals(chemical.molality);
            std::array<double,3> areas{};for(const auto& patch:map.patches[r])areas[patch.mineral]+=patch.dissolution;
            f<<r<<','<<region.parent<<','<<!region.bulk<<','<<region.connected<<','<<region.fraction;
            for(double a:areas)f<<','<<a;for(int j:{0,1,3,4})f<<','<<totals[j];for(int m=0;m<3;m++)f<<','<<chemical.saturationIndex[m];
            f<<','<<chemical.waterActivity<<','<<chemical.logActivity[saltchem::KIon]<<','<<chemical.logActivity[saltchem::MgIon]<<','<<chemical.logActivity[saltchem::ClIon]<<'\n';
        }
    }
    std::uint64_t nuclei(int mineral)const{std::uint64_t n=0;for(const auto& e:geometry.nucleationEvents)for(int s=0;s<3;s++)n+=e[3*mineral+s];return n;}
    double closedFraction(int i)const{return enabled?geometry.storedClosedFraction[index[i]]:0.;}
    void addClosedInventory(std::array<double,8>& total,const saltchem::MineralVector& volume)const{if(enabled){for(const auto& a:geometry.closedInventory)for(int j=0;j<8;j++)total[j]+=a[j]*geometry.coarseVolume();auto debt=saltdual::mineralDebtComponents(geometry,volume);for(int j=0;j<8;j++)total[j]+=debt[j]*geometry.coarseVolume();}}
    double mineralDebt(int mineral,double molarVolume)const{return enabled&&mineral<3?geometry.chemicalVolumeDebt[mineral]*geometry.coarseVolume()/molarVolume:0.;}
};
}
