// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <sstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include "GeometryCompetition.hpp"
#include "DissolutionSweep.hpp"

namespace saltdual {
inline void Geometry::advance(const std::vector<std::array<double,3>>& target,bool unblocked,double maxMove,const std::vector<SurfaceBudget>* supplied){
    int phases=0;if(supplied)for(const auto& budget:*supplied){if(budget.mineral<0||budget.mineral>=3)throw std::invalid_argument("Invalid interface mineral");if(budget.change>0)phases|=1<<budget.mineral;}
    if(!supplied||!phases||(phases&(phases-1))==0){advanceImpl(target,unblocked,maxMove,supplied,true);return;}
    const auto original=*this;
    try{advanceImpl(target,unblocked,maxMove,supplied,false);}
    catch(const StepFailure& error){*this=original;if(error.scale!=.5)throw;advanceImpl(target,unblocked,maxMove,supplied,true);}
}
inline void Geometry::advanceImpl(const std::vector<std::array<double,3>>& target,bool unblocked,double maxMove,const std::vector<SurfaceBudget>* supplied,bool fullGrowthBand){
    if(target.size()!=fractions.size())throw std::invalid_argument("Dual geometry target size mismatch");
    const bool connected=supplied!=nullptr;
    std::map<int,Cell> retainedDomains;
    if(connected)for(int r=0;r<int(waterRegions.size());r++)if(!waterRegions[r].bulk){
        auto& domain=retainedDomains[r];for(auto [cell,k]:waterRegions[r].pieces){auto piece=cells[cell][k];auto o=origin(cell);for(auto& point:piece.polygon)point=point+o;for(auto& ring:piece.holes)for(auto& point:ring)point=point+o;domain.push_back(std::move(piece));}
    }
    std::array<long double,3> targetTotal{};for(const auto& a:target)for(int m=0;m<3;m++)targetTotal[m]+=a[m];
    constexpr double volumeQuantum=1e-10;
    std::array<double,3> unrepresented{};
    std::vector<SurfaceBudget> automatic;
    std::vector<std::vector<const SurfaceBudget*>> work;
    std::vector<std::vector<int>> taskCells;
    std::map<const SurfaceBudget*,std::vector<int>> reach;
    std::vector<int> active;
    if(!connected){
        std::map<std::pair<int,int>,std::vector<int>> dissolutionFaces,growthFaces;
        for(int k=0;k<int(wet.size());k++){const auto& s=wet[k];dissolutionFaces[{parent(s.solid),s.mineral}].push_back(k);growthFaces[{parent(s.water),s.mineral}].push_back(k);}
        if(unblocked)for(int k=0;k<int(buried.size());k++){const auto& s=buried[k];dissolutionFaces[{parent(s.solid),s.mineral}].push_back(-k-1);}
        for(int p=0;p<nx*ny;p++)for(int m=0;m<3;m++){double change=target[p][m]-fractions[p][m];if(change!=0)automatic.push_back({p,m,change,(change>0?growthFaces:dissolutionFaces)[{p,m}]});}
        work.resize(nx*ny);taskCells.resize(nx*ny);
        for(const auto& b:automatic){if(work[b.parent].empty()){active.push_back(b.parent);for(int y=b.parent/nx*ratio;y<(b.parent/nx+1)*ratio;y++)for(int x=b.parent%nx*ratio;x<(b.parent%nx+1)*ratio;x++)taskCells[b.parent].push_back(y*fx()+x);}work[b.parent].push_back(&b);reach[&b]=taskCells[b.parent];}
    }else{
        // One conservative budget per liquid region and mineral. The moving
        // interface can cross a main-cell boundary; geometry owns the final
        // solid locations, while the reacting liquid owns its aqueous source.
        std::map<std::tuple<int,int,int>,SurfaceBudget> merged;int separate=-1;
        for(const auto& b:*supplied){if(b.change==0)continue;if(b.parent<0||b.parent>=nx*ny||b.mineral<0||b.mineral>=3||!std::isfinite(b.change))throw std::invalid_argument("Invalid interface reaction budget");
            auto key=std::make_tuple(b.liquidRegion>=0?b.liquidRegion:separate--,b.mineral,b.change>0?1:-1);auto entry=merged.emplace(key,SurfaceBudget{b.parent,b.mineral,0,{},b.liquidRegion});auto& destination=entry.first->second;
            destination.change+=b.change;destination.surfaces.insert(destination.surfaces.end(),b.surfaces.begin(),b.surfaces.end());}
        for(auto& entry:merged){auto& b=entry.second;std::sort(b.surfaces.begin(),b.surfaces.end());b.surfaces.erase(std::unique(b.surfaces.begin(),b.surfaces.end()),b.surfaces.end());automatic.push_back(std::move(b));}
        std::array<long double,3> difference{};for(int p=0;p<nx*ny;p++)for(int m=0;m<3;m++)difference[m]+=static_cast<long double>(target[p][m])-fractions[p][m];for(const auto& b:automatic)difference[b.mineral]-=b.change;
        for(int m=0;m<3;m++)if(std::abs(difference[m])>1e-13){auto best=automatic.end();for(auto it=automatic.begin();it!=automatic.end();++it)if(it->mineral==m&&(best==automatic.end()||std::abs(it->change)>std::abs(best->change)))best=it;if(best!=automatic.end())best->change+=double(difference[m]);
            else{SurfaceBudget pending{0,m,double(difference[m]),{},-1};for(int k=0;k<int(wet.size());k++)if(wet[k].mineral==m)pending.surfaces.push_back(k);if(!pending.surfaces.empty())automatic.push_back(std::move(pending));else unrepresented[m]=double(difference[m]);}}
        std::vector<int> roots(automatic.size());std::iota(roots.begin(),roots.end(),0);auto root=[&](int a){while(roots[a]!=a){roots[a]=roots[roots[a]];a=roots[a];}return a;};std::map<int,int> first;
        for(int i=0;i<int(automatic.size());i++){
            auto& b=automatic[i];std::set<int> affected;
            for(int id:b.surfaces){const auto& face=id>=0?wet.at(id):buried.at(-id-1);int x0=std::max(0,int(std::floor(std::min(face.a.x,face.b.x)-maxMove-1e-9))),x1=std::min(fx()-1,int(std::floor(std::max(face.a.x,face.b.x)+maxMove+1e-9)));int y0=std::max(0,int(std::floor(std::min(face.a.y,face.b.y)-maxMove-1e-9))),y1=std::min(fy()-1,int(std::floor(std::max(face.a.y,face.b.y)+maxMove+1e-9)));for(int y=y0;y<=y1;y++)for(int x=x0;x<=x1;x++)affected.insert(y*fx()+x);}
            if(affected.empty()&&std::abs(b.change)>1e-10)throw StepFailure("Reaction budget has no physical interface");
            reach[&b]={affected.begin(),affected.end()};for(int cell:affected){auto entry=first.emplace(cell,i);if(!entry.second)roots[root(i)]=root(entry.first->second);}
        }
        std::map<int,int> groups;
        for(int i=0;i<int(automatic.size());i++){if(reach[&automatic[i]].empty())continue;auto entry=groups.emplace(root(i),work.size());if(entry.second){active.push_back(work.size());work.emplace_back();taskCells.emplace_back();}int group=entry.first->second;work[group].push_back(&automatic[i]);const auto& ids=reach[&automatic[i]];taskCells[group].insert(taskCells[group].end(),ids.begin(),ids.end());}
        for(auto& ids:taskCells){std::sort(ids.begin(),ids.end());ids.erase(std::unique(ids.begin(),ids.end()),ids.end());}
    }
    std::vector<std::string> failures(work.size());std::vector<double> retryScale(work.size(),.5);std::vector<unsigned char> retryable(work.size());std::uint64_t repaired=0;
    const bool profiling=std::getenv("SALT_DUAL_GEOMETRY_PROFILE");std::vector<std::array<double,6>> profile(work.size());
    lastMovementChanged=false;std::vector<unsigned char> changed(work.size());
    createdWater.clear();std::vector<std::vector<CreatedWater>> creation(work.size());
    std::vector<std::array<double,3>> deferred(work.size());
    sweptVolumeVectors.assign(nx*ny,{});std::vector<std::map<int,Point>> movement(work.size());
    // Each owner prepares every mineral from one common material partition.
    // Intersecting swept regions trigger a first-contact event before commit.
    // Consequently neither mineral numbering nor loop order claims pore space.
    #pragma omp parallel for schedule(dynamic,1) reduction(+:repaired) if(active.size()>1)
    for(int task=0;task<int(active.size());task++){
        int p=active[task];const auto repairsBefore=booleanContext().fallbacks;
        auto profileStart=std::chrono::steady_clock::now();auto mark=[&](int k){if(profiling){auto now=std::chrono::steady_clock::now();profile[p][k]+=std::chrono::duration<double>(now-profileStart).count();profileStart=now;}};
        try{
            BooleanEnvironment environment;auto c=booleanContext().handle;
            const auto& members=taskCells[p];
            Cell protectedPieces;
            if(connected)for(const auto& entry:retainedDomains){bool nearby=false;for(auto [cell,k]:waterRegions[entry.first].pieces)if(std::binary_search(members.begin(),members.end(),cell)){nearby=true;break;}if(nearby)protectedPieces.insert(protectedPieces.end(),entry.second.begin(),entry.second.end());}
            auto protectedWater=unionPieces(protectedPieces);
            auto parameters=GEOSBufferParams_create_r(c);
            struct BufferGuard{GEOSContextHandle_t c;GEOSBufferParams* p;~BufferGuard(){GEOSBufferParams_destroy_r(c,p);}} guard{c,parameters};
            GEOSBufferParams_setEndCapStyle_r(c,parameters,GEOSBUF_CAP_ROUND);GEOSBufferParams_setJoinStyle_r(c,parameters,GEOSBUF_JOIN_MITRE);
            GEOSBufferParams_setMitreLimit_r(c,parameters,2);GEOSBufferParams_setQuadrantSegments_r(c,parameters,4);
            for(int stage=0;stage<2;stage++){
                if(stage==0&&connected){changed[p]|=sweepInterfaces(*this,members,work[p],maxMove,targetTotal,creation[p],movement[p],deferred[p]);mark(0);continue;}
                if(stage==1&&connected){int phases=0;for(const auto* budget:work[p])if(budget->change>0)phases|=1<<budget->mineral;
                    if(phases&&(phases&(phases-1))==0){changed[p]|=sweepInterfaces(*this,members,work[p],maxMove,targetTotal,creation[p],movement[p],deferred[p],true,&retainedDomains);mark(1);continue;}}
                std::map<int,Geos> sources;std::map<const GEOSGeometry*,Geos> restrictedSources;
                auto source=[&](int phase)->const GEOSGeometry*{
                    auto found=sources.find(phase);if(found!=sources.end())return found->second.get();
                    Cell pieces;for(int g:members)for(const auto& piece:cells[g])if(piece.phase==phase){auto q=piece;for(auto& v:q.polygon)v=v+origin(g);for(auto& h:q.holes)for(auto& v:h)v=v+origin(g);pieces.push_back(std::move(q));}
                    return sources.emplace(phase,unionPieces(pieces)).first->second.get();
                };
                struct Candidate{int mineral,from,liquid,reactor;double distance,wanted,initialWanted;Geos boundary,cut;std::vector<std::pair<Point,Point>> segments;std::vector<int> cells;Point normal;};std::vector<Candidate> candidates;
                auto swept=[&](const GEOSGeometry* boundary,int from,int liquid,double distance){auto band=owned(GEOSBufferWithParams_r(c,boundary,parameters,distance));auto found=restrictedSources.find(boundary);return intersection(found!=restrictedSources.end()?found->second.get():source(from),band.get());};
                for(const auto* budget:work[p]){
                    const auto& b=*budget;bool growing=b.change>0;if(growing!=(stage==1))continue;
                    if(std::abs(b.change)<volumeQuantum&&targetTotal[b.mineral]>0){deferred[p][b.mineral]+=b.change;continue;}
                    if(b.surfaces.empty()){if(std::abs(b.change)<1e-10)continue;throw StepFailure("Reaction budget has no liquid contact");}
                    std::vector<GEOSGeometry*> lines;std::vector<std::pair<Point,Point>> segments;Point normal{};double frontLength=0;
                    for(int id:b.surfaces){const auto& s=id>=0?wet.at(id):buried.at(-id-1);
                        segments.push_back({s.a,s.b});
                        double l=length(s.b-s.a);normal=normal+s.normal*l;frontLength+=l;
                        auto* seq=GEOSCoordSeq_create_r(c,2,2);GEOSCoordSeq_setX_r(c,seq,0,s.a.x);GEOSCoordSeq_setY_r(c,seq,0,s.a.y);GEOSCoordSeq_setX_r(c,seq,1,s.b.x);GEOSCoordSeq_setY_r(c,seq,1,s.b.y);lines.push_back(GEOSGeom_createLineString_r(c,seq));}
                    auto collection=owned(GEOSGeom_createCollection_r(c,GEOS_MULTILINESTRING,lines.data(),lines.size()));
                    auto united=owned(GEOSUnaryUnionPrec_r(c,collection.get(),overlayPrecision));auto boundary=owned(GEOSLineMerge_r(c,united.get()));
                    int from=growing?-1:b.mineral,liquid=-1;
                    if(connected){
                        Cell pieces;for(int cell:reach.at(budget))for(const auto& piece:cells[cell])if(piece.phase==from){auto q=piece;for(auto& v:q.polygon)v=v+origin(cell);for(auto& h:q.holes)for(auto& v:h)v=v+origin(cell);pieces.push_back(std::move(q));}
                        auto available=unionPieces(pieces);
                        if(growing&&retainedDomains.count(b.liquidRegion)){
                            // A local-water reaction consumes that reservoir's
                            // pore space. Crossing a narrow throat must not
                            // assign its tiny volume budget to neighbouring
                            // bulk water and leave an artificial empty film.
                            auto domain=retainedDomains.at(b.liquidRegion);
                            for(const auto& created:creation[p])if(created.region==b.liquidRegion)for(auto piece:created.shape){auto o=origin(created.cell);for(auto& point:piece.polygon)point=point+o;for(auto& ring:piece.holes)for(auto& point:ring)point=point+o;domain.push_back(std::move(piece));}
                            auto localWater=unionPieces(domain);available=intersection(available.get(),localWater.get());
                        }else if(growing&&!protectedPieces.empty())available=difference(available.get(),protectedWater.get());
                        if(GEOSGeomTypeId_r(c,available.get())==GEOS_MULTIPOLYGON){
                            std::vector<GEOSGeometry*> touching;
                            for(int k=0;k<GEOSGetNumGeometries_r(c,available.get());k++){auto part=GEOSGetGeometryN_r(c,available.get(),k);double distance=0;GEOSDistance_r(c,part,boundary.get(),&distance);if(distance<=2e-9)touching.push_back(GEOSGeom_clone_r(c,part));}
                            available=owned(GEOSGeom_createCollection_r(c,GEOS_MULTIPOLYGON,touching.data(),touching.size()));
                        }
                        restrictedSources.emplace(boundary.get(),std::move(available));
                    }
                    double wanted=std::abs(b.change)*ratio*ratio/(1-epsilon);
                    double bound=connected&&growing&&!fullGrowthBand&&frontLength>0&&!retainedDomains.count(b.liquidRegion)?std::min(maxMove,std::max(4*wanted/frontLength,1e-7)):maxMove;
                    auto upper=swept(boundary.get(),from,liquid,bound);
                    candidates.push_back({b.mineral,from,liquid,b.liquidRegion,bound,wanted,wanted,std::move(boundary),std::move(upper),std::move(segments),reach.at(budget),frontLength>0?normal*(1/frontLength):Point{}});
                }
                mark(1);
                // Partition shared normal columns once, before correcting
                // volume. Same-mineral local-water patches then have monotone,
                // independent scalar solves instead of a coupled corner solve.
                std::vector<Geos> columns;
                for(std::size_t partitionPass=0;partitionPass<=candidates.size();partitionPass++){
                columns.clear();for(const auto& candidate:candidates)columns.push_back(owned(GEOSGeom_clone_r(c,candidate.cut.get())));
                std::vector<std::array<double,4>> candidateBounds;for(const auto& candidate:candidates)candidateBounds.push_back(geometryBounds(candidate.cut.get()));
                for(std::size_t i=0;i<candidates.size();i++)for(std::size_t j=i+1;j<candidates.size();j++){
                    if(candidates[i].mineral!=candidates[j].mineral||candidates[i].from!=candidates[j].from||!boundsOverlap(candidateBounds[i],candidateBounds[j]))continue;
                    auto overlap=intersection(candidates[i].cut.get(),candidates[j].cut.get());double amount=0;GEOSArea_r(c,overlap.get(),&amount);if(amount<=1e-10)continue;
                    auto divided=divideArrival(overlap.get(),candidates[i],candidates[j]);columns[i]=difference(columns[i].get(),divided.second.get());columns[j]=difference(columns[j].get(),divided.first.get());
                }
                bool expanded=false;for(std::size_t k=0;k<candidates.size();k++){double capacity=0;GEOSArea_r(c,columns[k].get(),&capacity);auto& candidate=candidates[k];if(candidate.wanted>capacity+1e-10&&candidate.distance<maxMove){candidate.distance=maxMove;candidate.cut=swept(candidate.boundary.get(),candidate.from,candidate.liquid,maxMove);expanded=true;}}
                if(!expanded)break;
                }
                mark(2);
                for(std::size_t k=0;k<candidates.size();k++){
                    auto& candidate=candidates[k];double lo=0,hi=candidate.distance,vlo=0,vhi=0;GEOSArea_r(c,columns[k].get(),&vhi);
                    double frontLength=0;for(auto [a,b]:candidate.segments)frontLength+=length(b-a);
                    double perimeter=0;GEOSLength_r(c,columns[k].get(),&perimeter);
                    const double volumeTolerance=std::max(1e-11,overlayPrecision*(frontLength+perimeter));
                    restrictedSources.insert_or_assign(candidate.boundary.get(),std::move(columns[k]));
                    auto defer=[&](){deferred[p][candidate.mineral]+=(stage==1?1:-1)*candidate.initialWanted*(1-epsilon)/(ratio*ratio);candidate.initialWanted=0;candidate.wanted=0;candidate.distance=0;candidate.cut=owned(GEOSGeom_createEmptyPolygon_r(c));};
                    if(candidate.wanted>vhi+std::max(1e-10,2*volumeTolerance)){
                        if(connected&&targetTotal[candidate.mineral]>0&&(candidate.wanted*(1-epsilon)/(ratio*ratio)<=1e-9||(candidate.wanted-vhi)*(1-epsilon)/(ratio*ratio)<=volumeQuantum)){
                            if(vhi>1e-15){candidate.distance=maxMove;candidate.cut=swept(candidate.boundary.get(),candidate.from,candidate.liquid,maxMove);candidate.wanted=vhi;}else defer();continue;}
                        std::ostringstream message;message<<std::setprecision(12)<<"Interface displacement/stock limit: phase="<<candidate.mineral<<" region="<<candidate.reactor<<" wanted="<<candidate.wanted<<" capacity="<<vhi;throw StepFailure(message.str(),std::clamp(vhi/candidate.wanted,.02,.999999));}
                    double bestError=INFINITY,bestDistance=0;
                    for(int it=0;it<64;it++){
                        candidate.distance=it<8&&vhi>vlo?lo+(std::min(candidate.wanted,vhi)-vlo)*(hi-lo)/(vhi-vlo):(lo+hi)*.5;
                        candidate.cut=swept(candidate.boundary.get(),candidate.from,candidate.liquid,candidate.distance);double amount=0;GEOSArea_r(c,candidate.cut.get(),&amount);
                        double error=std::abs(amount-candidate.wanted);if(error<bestError){bestError=error;bestDistance=candidate.distance;}
                        if(error<volumeTolerance)break;if(amount<candidate.wanted){lo=candidate.distance;vlo=amount;}else{hi=candidate.distance;vhi=amount;}
                    }
                    if(bestError>std::max(1e-10,2*volumeTolerance)){
                        if(connected&&candidate.wanted*(1-epsilon)/(ratio*ratio)<=1e-9&&targetTotal[candidate.mineral]>0){
                            if(bestError*(1-epsilon)/(ratio*ratio)>1e-9){defer();continue;}
                        }else{
                        std::ostringstream message;message<<std::setprecision(12)<<"Interface volume residual: phase="<<candidate.mineral<<" wanted="<<candidate.wanted<<" error="<<bestError;throw StepFailure(message.str());}
                    }
                    if(candidate.distance!=bestDistance){candidate.distance=bestDistance;candidate.cut=swept(candidate.boundary.get(),candidate.from,candidate.liquid,candidate.distance);}
                    GEOSArea_r(c,candidate.cut.get(),&candidate.wanted);
                    if(candidate.wanted<=1e-15)defer();
                }
                mark(3);
                candidates.erase(std::remove_if(candidates.begin(),candidates.end(),[](const auto& candidate){return candidate.wanted==0;}),candidates.end());
                double eventScale=1,eventGap=0;bool junction=false;
                std::vector<std::array<double,4>> eventBounds;eventBounds.reserve(candidates.size());
                for(const auto& candidate:candidates)eventBounds.push_back(geometryBounds(candidate.cut.get()));
                for(std::size_t a=0;a<candidates.size();a++)for(std::size_t b=a+1;b<candidates.size();b++){
                    auto& x=candidates[a];auto& y=candidates[b];if(x.from!=y.from||x.wanted==0||y.wanted==0||!boundsOverlap(eventBounds[a],eventBounds[b]))continue;
                    auto common=intersection(x.cut.get(),y.cut.get());double overlap=0;GEOSArea_r(c,common.get(),&overlap);if(overlap<=1e-10)continue;
                    if(x.mineral==y.mineral){junction=true;continue;}
                    double separation=0;GEOSDistance_r(c,x.boundary.get(),y.boundary.get(),&separation);if(separation<=interfaceContactTolerance){junction=true;continue;}
                    auto& small=x.wanted<y.wanted?x:y;
                    if(connected&&small.wanted*(1-epsilon)/(ratio*ratio)<=1e-9&&targetTotal[small.mineral]>0){deferred[p][small.mineral]+=(stage==1?1:-1)*small.initialWanted*(1-epsilon)/(ratio*ratio);small.initialWanted=0;small.wanted=0;small.distance=0;small.cut=owned(GEOSGeom_createEmptyPolygon_r(c));continue;}
                    double low=0,high=1;
                    for(int it=0;it<40;it++){double mid=(low+high)*.5;auto u=swept(x.boundary.get(),x.from,x.liquid,x.distance*mid),v=swept(y.boundary.get(),y.from,y.liquid,y.distance*mid),both=intersection(u.get(),v.get());double amount=0;GEOSArea_r(c,both.get(),&amount);if(amount<=5e-11)low=mid;else high=mid;}
                    if(low<eventScale){eventScale=low;eventGap=separation;}
                }
                if(eventScale<1){std::ostringstream message;message<<"Competing mineral interfaces reached their first contact: fraction="<<eventScale<<" initial_gap_h="<<eventGap;throw StepFailure(message.str(),std::clamp(eventScale,1e-6,.999999));}
                candidates.erase(std::remove_if(candidates.begin(),candidates.end(),[](const auto& candidate){return candidate.wanted==0;}),candidates.end());
                if(candidates.empty())continue;
                if(junction){
                    // Existing triple junctions advance together. Partition
                    // their overlaps by arrival time, then correct all normal
                    // displacements simultaneously to match chemical volumes.
                    auto correctJunction=[&](std::vector<Candidate>& candidates){
                    const int n=candidates.size();std::vector<double> jac(n*n),residual(n),direction;
                    const double junctionTolerance=volumeQuantum*ratio*ratio/(1-epsilon);
                    for(auto& candidate:candidates)candidate.cut=swept(candidate.boundary.get(),candidate.from,candidate.liquid,candidate.distance);
                    auto evaluate=[&](){auto cuts=competingCuts(candidates);double error=0;for(int i=0;i<n;i++){double amount=0;GEOSArea_r(c,cuts[i].get(),&amount);residual[i]=amount-candidates[i].wanted;error=std::max(error,std::abs(residual[i]));}return std::make_pair(std::move(cuts),error);};
                    for(int i=0;i<n;i++)jac[i*n+i]=std::max(1e-12,candidates[i].wanted/std::max(candidates[i].distance,1e-12));
                    auto current=evaluate();bool converged=current.second<=junctionTolerance,needJacobian=false,freshJacobian=false;
                    for(int iteration=0;iteration<40&&!converged;iteration++){
                        if(needJacobian){
                            const auto base=residual;
                            for(int column=0;column<n;column++){
                                auto& candidate=candidates[column];double priorDistance=candidate.distance;
                                double delta=std::max(128*overlayPrecision,1e-3*priorDistance);if(priorDistance+delta>maxMove)delta=-std::min(delta,.5*priorDistance);
                                auto priorCut=std::move(candidate.cut);candidate.distance=priorDistance+delta;candidate.cut=swept(candidate.boundary.get(),candidate.from,candidate.liquid,candidate.distance);
                                auto perturbed=evaluate();for(int row=0;row<n;row++)jac[row*n+column]=(residual[row]-base[row])/delta;
                                candidate.distance=priorDistance;candidate.cut=std::move(priorCut);
                            }
                            residual=base;needJacobian=false;freshJacobian=true;
                        }
                        auto a=jac,rhs=residual;for(double& v:rhs)v=-v;direction.assign(n,0);
                        bool valid=true;
                        for(int k=0;k<n;k++){int pivot=k;for(int i=k+1;i<n;i++)if(std::abs(a[i*n+k])>std::abs(a[pivot*n+k]))pivot=i;if(std::abs(a[pivot*n+k])<1e-20){valid=false;break;}if(pivot!=k){for(int j=k;j<n;j++)std::swap(a[k*n+j],a[pivot*n+j]);std::swap(rhs[k],rhs[pivot]);}for(int i=k+1;i<n;i++){double factor=a[i*n+k]/a[k*n+k];for(int j=k+1;j<n;j++)a[i*n+j]-=factor*a[k*n+j];rhs[i]-=factor*rhs[k];}}
                        if(valid)for(int i=n-1;i>=0;i--){double v=rhs[i];for(int j=i+1;j<n;j++)v-=a[i*n+j]*direction[j];direction[i]=v/a[i*n+i];}
                        else for(int i=0;i<n;i++)direction[i]=-residual[i]/std::max(1e-12,candidates[i].wanted/std::max(candidates[i].distance,1e-12));
                        std::vector<double> before(n);for(int i=0;i<n;i++)before[i]=candidates[i].distance;const auto prior=residual;bool accepted=false;
                        for(double alpha=1;alpha>=1./4096;alpha*=.5){
                            bool feasible=true;for(int i=0;i<n;i++){double d=before[i]+alpha*direction[i];if(!(d>0&&d<=maxMove)){feasible=false;break;}}if(!feasible)continue;
                            for(int i=0;i<n;i++){auto& candidate=candidates[i];candidate.distance=before[i]+alpha*direction[i];candidate.cut=swept(candidate.boundary.get(),candidate.from,candidate.liquid,candidate.distance);}
                            auto trial=evaluate();if(trial.second<current.second||trial.second<=junctionTolerance){
                                double denominator=0;for(int i=0;i<n;i++)denominator+=std::pow(alpha*direction[i],2);
                                if(denominator>1e-60)for(int i=0;i<n;i++){double update=residual[i]-prior[i];for(int j=0;j<n;j++)update-=jac[i*n+j]*alpha*direction[j];for(int j=0;j<n;j++)jac[i*n+j]+=update*alpha*direction[j]/denominator;}
                                current=std::move(trial);accepted=true;freshJacobian=false;converged=current.second<=junctionTolerance;break;
                            }
                        }
                        if(!accepted){
                            for(int i=0;i<n;i++){auto& candidate=candidates[i];candidate.distance=before[i];candidate.cut=swept(candidate.boundary.get(),candidate.from,candidate.liquid,candidate.distance);}residual=prior;
                            if(!freshJacobian){needJacobian=true;continue;}
                            std::ostringstream message;message<<std::setprecision(12)<<"Competing interface junction group="<<p<<" stage="<<stage<<" error="<<current.second<<" candidates=";for(int i=0;i<n;i++)message<<candidates[i].mineral<<':'<<candidates[i].reactor<<':'<<candidates[i].wanted<<':'<<residual[i]<<',';throw StepFailure(message.str());}
                    }
                    if(!converged)throw StepFailure("Competing interface volume correction exceeded tolerance");
                    for(int i=0;i<n;i++)candidates[i].cut=std::move(current.first[i]);
                    };
                    // Correct only mutually competing fronts. Independent
                    // nuclei keep their already converged scalar volumes.
                    // If a correction reaches another cluster, merge those
                    // clusters and resolve their conservative contact together.
                    std::vector<int> roots(candidates.size());std::iota(roots.begin(),roots.end(),0);auto root=[&](int a){while(roots[a]!=a){roots[a]=roots[roots[a]];a=roots[a];}return a;};
                    auto mergeContacts=[&](){bool joined=false;std::vector<std::array<double,4>> bounds;for(const auto& candidate:candidates)bounds.push_back(geometryBounds(candidate.cut.get()));
                        for(int i=0;i<int(candidates.size());i++)for(int j=i+1;j<int(candidates.size());j++)if(root(i)!=root(j)&&candidates[i].from==candidates[j].from&&boundsOverlap(bounds[i],bounds[j])){auto common=intersection(candidates[i].cut.get(),candidates[j].cut.get());double amount=0;GEOSArea_r(c,common.get(),&amount);if(amount>1e-10){roots[root(i)]=root(j);joined=true;}}return joined;};
                    mergeContacts();
                    for(std::size_t pass=0;pass<candidates.size();pass++){
                        std::map<int,std::vector<int>> clusters;for(int i=0;i<int(candidates.size());i++)clusters[root(i)].push_back(i);
                        for(const auto& cluster:clusters)if(cluster.second.size()>1){std::vector<Candidate> local;for(int i:cluster.second)local.push_back(std::move(candidates[i]));correctJunction(local);for(std::size_t k=0;k<local.size();k++)candidates[cluster.second[k]]=std::move(local[k]);}
                        if(!mergeContacts())break;
                        if(pass+1==candidates.size())throw StepFailure("Competing interface clusters did not stabilise");
                    }
                }
                for(const auto& candidate:candidates){double amount=0;GEOSArea_r(c,candidate.cut.get(),&amount);deferred[p][candidate.mineral]+=(stage==1?1:-1)*(candidate.initialWanted-amount)*(1-epsilon)/(ratio*ratio);}
                mark(4);
                std::vector<CellEdit> edits;
                for(int g:members){
                    auto offset=origin(g);auto box=square();for(auto& v:box)v=v+offset;auto cellBox=polygon(box);
                    std::array<std::vector<Cell>,3> changes;
                    for(const auto& candidate:candidates){if(!std::binary_search(candidate.cells.begin(),candidate.cells.end(),g))continue;auto part=intersection(candidate.cut.get(),cellBox.get());double amount=0;GEOSArea_r(c,part.get(),&amount);if(amount<1e-15)continue;
                        auto local=owned(GEOSGeom_transformXY_r(c,part.get(),[](double* x,double* y,void* ptr){auto o=static_cast<Point*>(ptr);*x-=o->x;*y-=o->y;return 1;},&offset));
                        if(stage==0&&candidate.reactor>=0){CreatedWater water{candidate.reactor,g,{}};unpackPolygon(local.get(),-1,water.shape);creation[p].push_back(std::move(water));}
                        movement[p][parent(g)]=movement[p][parent(g)]+candidate.normal*((stage==0?1:-1)*amount*(1-epsilon)/(ratio*ratio));
                        Cell coordinates;unpackPolygon(local.get(),candidate.mineral,coordinates,false);changes[candidate.mineral].push_back(std::move(coordinates));}
                    if(std::all_of(changes.begin(),changes.end(),[](const auto& a){return a.empty();}))continue;
                    edits.push_back({g,std::move(changes)});
                }
                auto removed=commitCellEdits(cells,edits,stage==1);for(const auto& values:removed)for(int m=0;m<3;m++)deferred[p][m]+=values[m]*(1-epsilon)/(ratio*ratio);
                changed[p]|=!edits.empty();
                mark(5);
            }
        }catch(const StepFailure& error){failures[p]=error.what();retryable[p]=1;retryScale[p]=error.scale;}
        catch(const std::exception& error){failures[p]=error.what();}
        repaired+=booleanContext().fallbacks-repairsBefore;
    }
    booleanFallbacks+=repaired;
    if(profiling){std::array<double,6> sums{};for(const auto& times:profile)for(int k=0;k<6;k++)sums[k]+=times[k];std::cerr<<"GEOMETRY_PROFILE wet="<<wet.size()<<" groups="<<active.size();for(double value:sums)std::cerr<<' '<<value;std::cerr<<'\n';}
    double scale=1;std::string failure;
    for(int p:active)if(!failures[p].empty()){if(!retryable[p])throw std::runtime_error(failures[p]);failure=failures[p];scale=std::min(scale,retryScale[p]);}
    if(!failure.empty())throw StepFailure(failure,scale);
    for(int p:active)for(auto& water:creation[p])createdWater.push_back(std::move(water));
    for(int p:active)for(const auto& entry:movement[p])sweptVolumeVectors[entry.first]=sweptVolumeVectors[entry.first]+entry.second;
    lastMovementChanged=std::any_of(changed.begin(),changed.end(),[](unsigned char value){return value!=0;});if(lastMovementChanged)rebuild();
    chemicalVolumeError=0;
    if(connected){std::array<long double,3> error{};for(int p=0;p<nx*ny;p++)for(int m=0;m<3;m++)error[m]+=static_cast<long double>(fractions[p][m])-target[p][m];for(int m=0;m<3;m++){chemicalVolumeDebt[m]=-double(error[m]);error[m]+=unrepresented[m];for(const auto& a:deferred)error[m]+=a[m];chemicalVolumeError=std::max(chemicalVolumeError,double(std::abs(error[m])));}}
    else for(int p:active)for(int m=0;m<3;m++)chemicalVolumeError=std::max(chemicalVolumeError,std::abs(fractions[p][m]-target[p][m]));
    if(chemicalVolumeError>1e-9){std::ostringstream message;message<<std::setprecision(17)<<"Fine geometry and conserved mineral inventory disagree: "<<chemicalVolumeError;throw StepFailure(message.str());}
}
}
