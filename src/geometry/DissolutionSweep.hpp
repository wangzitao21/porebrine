// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
namespace saltdual {
inline Geos intersectLines(const GEOSGeometry* line,const GEOSGeometry* material){
    auto c=booleanContext().handle;auto value=GEOSIntersectionPrec_r(c,line,material,overlayPrecision);
    if(!value){booleanContext().fallbacks++;value=GEOSIntersection_r(c,line,material);}return owned(value);
}
inline void lineCells(const GEOSGeometry* line,const Geometry& g,double distance,std::set<int>& cells){
    auto c=booleanContext().handle;int type=GEOSGeomTypeId_r(c,line);
    if(type==GEOS_LINESTRING||type==GEOS_LINEARRING){
        auto coordinates=GEOSGeom_getCoordSeq_r(c,line);unsigned size=0;GEOSCoordSeq_getSize_r(c,coordinates,&size);
        Point a;for(unsigned k=0;k<size;k++){Point b;GEOSCoordSeq_getX_r(c,coordinates,k,&b.x);GEOSCoordSeq_getY_r(c,coordinates,k,&b.y);if(k){
            int x0=std::max(0,int(std::floor(std::min(a.x,b.x)-distance-1e-9))),x1=std::min(g.fx()-1,int(std::floor(std::max(a.x,b.x)+distance+1e-9)));
            int y0=std::max(0,int(std::floor(std::min(a.y,b.y)-distance-1e-9))),y1=std::min(g.fy()-1,int(std::floor(std::max(a.y,b.y)+distance+1e-9)));
            for(int y=y0;y<=y1;y++)for(int x=x0;x<=x1;x++)cells.insert(y*g.fx()+x);
        }a=b;}
    }else if(type==GEOS_MULTILINESTRING||type==GEOS_GEOMETRYCOLLECTION)for(int k=0;k<GEOSGetNumGeometries_r(c,line);k++)lineCells(GEOSGetGeometryN_r(c,line,k),g,distance,cells);
}

// A continuous mineral boundary has one moving envelope, even when several
// local liquids supply different fluxes along it. Union their weighted normal
// sweeps and correct the total volume on each connected mineral support.
// This avoids treating every chemical control-volume edge as a new collision.
inline bool sweepInterfaces(Geometry& g,const std::vector<int>& members,const std::vector<const SurfaceBudget*>& budgets,
    double maxMove,const std::array<long double,3>& targetTotal,std::vector<CreatedWater>& creation,
    std::map<int,Point>& movement,std::array<double,3>& deferred,bool growing=false,const std::map<int,Cell>* waterDomains=nullptr){
    auto context=booleanContext().handle;
    auto parameters=GEOSBufferParams_create_r(context);
    struct Guard{GEOSContextHandle_t c;GEOSBufferParams* p;~Guard(){GEOSBufferParams_destroy_r(c,p);}} guard{context,parameters};
    GEOSBufferParams_setEndCapStyle_r(context,parameters,GEOSBUF_CAP_ROUND);GEOSBufferParams_setJoinStyle_r(context,parameters,GEOSBUF_JOIN_MITRE);
    GEOSBufferParams_setMitreLimit_r(context,parameters,2);GEOSBufferParams_setQuadrantSegments_r(context,parameters,4);
    struct Flux {int mineral,region;double demand,length=0,assigned=0;Geos front;Point normal;};std::vector<Flux> fluxes;
    for(const auto* pointer:budgets){const auto& b=*pointer;if(b.change==0||(b.change>0)!=growing)continue;
        if(std::abs(b.change)<1e-10&&targetTotal[b.mineral]>0){deferred[b.mineral]+=b.change;continue;}
        std::map<int,std::vector<int>> regions;double totalLength=0;
        for(int id:b.surfaces){const auto& s=id>=0?g.wet.at(id):g.buried.at(-id-1);int r=id>=0?s.liquidRegion:b.liquidRegion;regions[r].push_back(id);totalLength+=length(s.b-s.a);}
        if(totalLength==0){if(std::abs(b.change)<=1e-9){deferred[b.mineral]+=b.change;continue;}throw StepFailure("Dissolution has no exposed mineral interface");}
        for(const auto& entry:regions){std::vector<GEOSGeometry*> lines;double l=0;Point normal{};
            for(int id:entry.second){const auto& s=id>=0?g.wet[id]:g.buried[-id-1];double weight=length(s.b-s.a);l+=weight;normal=normal+s.normal*weight;
                auto* coordinates=GEOSCoordSeq_create_r(context,2,2);GEOSCoordSeq_setX_r(context,coordinates,0,s.a.x);GEOSCoordSeq_setY_r(context,coordinates,0,s.a.y);GEOSCoordSeq_setX_r(context,coordinates,1,s.b.x);GEOSCoordSeq_setY_r(context,coordinates,1,s.b.y);lines.push_back(GEOSGeom_createLineString_r(context,coordinates));}
            auto collection=owned(GEOSGeom_createCollection_r(context,GEOS_MULTILINESTRING,lines.data(),lines.size()));auto united=owned(GEOSUnaryUnionPrec_r(context,collection.get(),overlayPrecision));auto front=owned(GEOSLineMerge_r(context,united.get()));
            fluxes.push_back({b.mineral,entry.first,std::abs(b.change)*g.ratio*g.ratio/(1-g.epsilon)*l/totalLength,l,0,std::move(front),normal*(1/l)});
        }
    }
    if(fluxes.empty())return false;
    struct Part{int flux;double length,step=0;Geos front,material,cut;double demand=0;};
    struct Body{int mineral;Geos material,cut;std::vector<Part> parts;std::vector<int> cells;double demand=0,capacity=0;};std::vector<Body> bodies;
    for(int m=0;m<3;m++){
        if(std::none_of(fluxes.begin(),fluxes.end(),[&](const Flux& f){return f.mineral==m;}))continue;
        Cell material;for(int cell:members)for(const auto& piece:g.cells[cell])if(growing?piece.phase<0:piece.phase==m){auto p=piece;for(auto& v:p.polygon)v=v+g.origin(cell);for(auto& h:p.holes)for(auto& v:h)v=v+g.origin(cell);material.push_back(std::move(p));}
        auto source=unionPieces(material);std::vector<std::pair<int,Geos>> domains;
        if(growing&&waterDomains){Cell protectedPieces;
            for(const auto& entry:*waterDomains){bool nearby=false;for(auto [cell,k]:g.waterRegions[entry.first].pieces)if(std::binary_search(members.begin(),members.end(),cell)){nearby=true;break;}if(!nearby)continue;
                Cell pieces=entry.second;for(const auto& created:creation)if(created.region==entry.first)for(auto piece:created.shape){auto o=g.origin(created.cell);for(auto& v:piece.polygon)v=v+o;for(auto& ring:piece.holes)for(auto& v:ring)v=v+o;pieces.push_back(std::move(piece));}
                protectedPieces.insert(protectedPieces.end(),pieces.begin(),pieces.end());auto shape=unionPieces(pieces);domains.emplace_back(entry.first,intersection(source.get(),shape.get()));
            }
            if(!protectedPieces.empty()){auto shape=unionPieces(protectedPieces);source=difference(source.get(),shape.get());}
        }
        domains.emplace_back(-1,std::move(source));
        for(auto& domain:domains){std::vector<GEOSGeometry*> components;collectPolygons(domain.second.get(),components);
        const bool single=components.size()==1;
        for(auto* component:components){Body body{m,owned(component),Geos(nullptr,{context}),{},{},0};
            auto support=owned(GEOSBufferWithParams_r(context,body.material.get(),parameters,4e-10));std::set<int> affected;
            for(int k=0;k<int(fluxes.size());k++)if(fluxes[k].mineral==m&&(!growing||(waterDomains&&waterDomains->count(fluxes[k].region)?fluxes[k].region:-1)==domain.first)){auto front=single?owned(GEOSGeom_clone_r(context,fluxes[k].front.get())):intersectLines(fluxes[k].front.get(),support.get());double l=0;GEOSLength_r(context,front.get(),&l);if(l<=tiny)continue;
                fluxes[k].assigned+=l;std::set<int> local;lineCells(front.get(),g,maxMove,local);affected.insert(local.begin(),local.end());
                int xmin=g.fx(),xmax=0,ymin=g.fy(),ymax=0;for(int cell:local){xmin=std::min(xmin,cell%g.fx());xmax=std::max(xmax,cell%g.fx()+1);ymin=std::min(ymin,cell/g.fx());ymax=std::max(ymax,cell/g.fx()+1);}
                auto box=polygon({{double(xmin),double(ymin)},{double(xmax),double(ymin)},{double(xmax),double(ymax)},{double(xmin),double(ymax)}});auto nearby=intersection(body.material.get(),box.get());
                body.parts.push_back({k,l,0,std::move(front),std::move(nearby),Geos(nullptr,{context})});}
            if(!body.parts.empty()){body.cells.assign(affected.begin(),affected.end());bodies.push_back(std::move(body));}
        }
        }
    }
    const double sign=growing?1:-1;
    for(const auto& flux:fluxes)if(flux.assigned==0){if(flux.demand*(1-g.epsilon)/(g.ratio*g.ratio)>1e-9)throw StepFailure("Exposed interface has no available material");deferred[flux.mineral]+=sign*flux.demand*(1-g.epsilon)/(g.ratio*g.ratio);}
    std::vector<double> excess(fluxes.size());
    for(auto& body:bodies){GEOSArea_r(context,body.material.get(),&body.capacity);
        for(auto& part:body.parts){part.demand=fluxes[part.flux].demand*part.length/fluxes[part.flux].assigned;body.demand+=part.demand;}
        if(body.demand>body.capacity){double factor=body.capacity/body.demand;for(auto& part:body.parts){double kept=part.demand*factor;excess[part.flux]+=part.demand-kept;part.demand=kept;}body.demand=body.capacity;}
    }
    std::vector<std::vector<std::pair<int,int>>> supports(fluxes.size());
    for(int b=0;b<int(bodies.size());b++)for(int p=0;p<int(bodies[b].parts.size());p++)supports[bodies[b].parts[p].flux].push_back({b,p});
    auto reassign=[&](int missing){
        // A water region with only one accessible crystal must not lose its
        // stock to a region that can also react with another crystal. Follow
        // a residual-capacity path and move those other allocations to their
        // own exposed supports, keeping every region's reaction amount fixed.
        std::vector<int> viaBody(fluxes.size(),-2),viaPart(fluxes.size(),-1),fromFlux(bodies.size(),-1),fromPart(bodies.size(),-1);
        std::queue<int> queue;queue.push(missing);viaBody[missing]=-1;int destination=-1;
        while(!queue.empty()&&destination<0){int f=queue.front();queue.pop();
            for(auto [b,p]:supports[f]){if(fromFlux[b]>=0)continue;fromFlux[b]=f;fromPart[b]=p;
                if(bodies[b].capacity>bodies[b].demand+1e-14){destination=b;break;}
                for(int k=0;k<int(bodies[b].parts.size());k++){const auto& part=bodies[b].parts[k];if(part.demand>1e-14&&viaBody[part.flux]==-2){viaBody[part.flux]=b;viaPart[part.flux]=k;queue.push(part.flux);}}
            }
        }
        if(destination<0)return false;
        double amount=std::min(excess[missing],bodies[destination].capacity-bodies[destination].demand);
        for(int b=destination;;){int f=fromFlux[b];if(f==missing)break;int prior=viaBody[f];amount=std::min(amount,bodies[prior].parts[viaPart[f]].demand);b=prior;}
        if(!(amount>0))return false;
        for(int b=destination;;){int f=fromFlux[b];bodies[b].parts[fromPart[b]].demand+=amount;bodies[b].demand+=amount;if(f==missing)break;
            int prior=viaBody[f];bodies[prior].parts[viaPart[f]].demand-=amount;bodies[prior].demand-=amount;b=prior;
        }
        excess[missing]=std::max(0.,excess[missing]-amount);return true;
    };
    // The chemical control volume supplies a phase amount, while individual
    // exposed crystals have finite stocks. Exhaust small crystals exactly and
    // distribute the remaining amount only over other exposed supports of the
    // same mineral in contact with that same reacting water.
    for(int flux=0;flux<int(excess.size());flux++)while(excess[flux]>1e-13){
        double weight=0;for(const auto& body:bodies)if(body.capacity>body.demand+1e-14)for(const auto& part:body.parts)if(part.flux==flux)weight+=part.length;
        if(weight==0){if(reassign(flux))continue;if(excess[flux]*(1-g.epsilon)/(g.ratio*g.ratio)<=1e-10){deferred[fluxes[flux].mineral]+=sign*excess[flux]*(1-g.epsilon)/(g.ratio*g.ratio);excess[flux]=0;break;}
            std::ostringstream message;message<<std::setprecision(17)<<"Local interface amount exceeds its connected mineral/pore stock: phase="<<fluxes[flux].mineral<<" region="<<fluxes[flux].region<<" growth="<<growing<<" demand_h2="<<fluxes[flux].demand<<" excess_h2="<<excess[flux]<<" exposed_h="<<fluxes[flux].assigned;
            throw StepFailure(message.str(),std::clamp(1-excess[flux]/fluxes[flux].demand,.02,.999999));}
        double available=excess[flux];for(auto& body:bodies)if(body.capacity>body.demand+1e-14)for(auto& part:body.parts)if(part.flux==flux){double amount=std::min(body.capacity-body.demand,available*part.length/weight);part.demand+=amount;body.demand+=amount;excess[flux]-=amount;}
    }
    for(auto& body:bodies){
        if(body.demand==0){body.cut=owned(GEOSGeom_createEmptyPolygon_r(context));continue;}
        double maximumStep=0,totalLength=0;
        for(auto& part:body.parts){part.step=part.demand/part.length;maximumStep=std::max(maximumStep,part.step);totalLength+=part.length;}
        if(body.demand*(1-g.epsilon)/(g.ratio*g.ratio)<1e-10&&targetTotal[body.mineral]>0){deferred[body.mineral]+=sign*body.demand*(1-g.epsilon)/(g.ratio*g.ratio);body.cut=owned(GEOSGeom_createEmptyPolygon_r(context));continue;}
        auto swept=[&](double scale){std::vector<GEOSGeometry*> parts;
            for(auto& part:body.parts){auto band=owned(GEOSBufferWithParams_r(context,part.front.get(),parameters,scale*part.step));part.cut=intersection(part.material.get(),band.get());parts.push_back(GEOSGeom_clone_r(context,part.cut.get()));}
            auto collection=owned(GEOSGeom_createCollection_r(context,GEOS_GEOMETRYCOLLECTION,parts.data(),parts.size()));auto result=GEOSUnaryUnionPrec_r(context,collection.get(),overlayPrecision);if(!result){booleanContext().fallbacks++;result=GEOSUnaryUnion_r(context,collection.get());}auto swept=polygonOnly(owned(result));
            if(GEOSGeomTypeId_r(context,swept.get())==GEOS_MULTIPOLYGON){std::vector<GEOSGeometry*> connected;
                for(int i=0;i<GEOSGetNumGeometries_r(context,swept.get());i++){auto piece=GEOSGetGeometryN_r(context,swept.get(),i);bool wet=false;for(const auto& part:body.parts){double d=0;GEOSDistance_r(context,piece,part.front.get(),&d);if(d<=4e-10){wet=true;break;}}if(wet)connected.push_back(GEOSGeom_clone_r(context,piece));}
                swept=owned(GEOSGeom_createCollection_r(context,GEOS_MULTIPOLYGON,connected.data(),connected.size()));
            }
            return swept;
        };
        double maximum=maxMove/maximumStep,lo=0,hi=std::min(1.,maximum),vlo=0,vhi=0,sourceArea=0;GEOSArea_r(context,body.material.get(),&sourceArea);
        const double tolerance=std::max(1e-11,overlayPrecision*(totalLength+1));
        body.cut=swept(hi);GEOSArea_r(context,body.cut.get(),&vhi);
        while(vhi<body.demand-tolerance&&hi<maximum&&vhi<sourceArea-tolerance){lo=hi;vlo=vhi;hi=std::min(maximum,hi*2);body.cut=swept(hi);GEOSArea_r(context,body.cut.get(),&vhi);}
        if(body.demand>vhi+std::max(1e-10,2*tolerance)){
            if((body.demand*(1-g.epsilon)/(g.ratio*g.ratio)<=1e-9||(body.demand-vhi)*(1-g.epsilon)/(g.ratio*g.ratio)<=1e-10)&&targetTotal[body.mineral]>0){deferred[body.mineral]+=sign*(body.demand-vhi)*(1-g.epsilon)/(g.ratio*g.ratio);continue;}
            std::ostringstream message;message<<std::setprecision(14)<<"Moving mineral reached its stock/displacement event: phase="<<body.mineral<<" growth="<<growing<<" wanted="<<body.demand<<" capacity="<<vhi<<" source="<<sourceArea;throw StepFailure(message.str(),std::clamp(vhi/body.demand,.02,.999999));
        }
        double bestError=std::abs(vhi-body.demand),best=hi;
        for(int iteration=0;iteration<64&&bestError>tolerance;iteration++){
            double scale=iteration<8&&vhi>vlo?lo+(std::min(body.demand,vhi)-vlo)*(hi-lo)/(vhi-vlo):(lo+hi)*.5;
            body.cut=swept(scale);double area=0;GEOSArea_r(context,body.cut.get(),&area);double error=std::abs(area-body.demand);if(error<bestError){bestError=error;best=scale;}
            if(area<body.demand){lo=scale;vlo=area;}else{hi=scale;vhi=area;}
        }
        if(bestError>std::max(1e-10,2*tolerance))throw StepFailure("Continuous dissolution sweep exceeded its volume tolerance");
        body.cut=swept(best);double area=0;GEOSArea_r(context,body.cut.get(),&area);deferred[body.mineral]+=sign*(body.demand-area)*(1-g.epsilon)/(g.ratio*g.ratio);
    }
    std::map<int,std::vector<int>> owners;for(int k=0;k<int(bodies.size());k++)if(bodies[k].cut)for(int cell:bodies[k].cells)owners[cell].push_back(k);
    std::vector<CellEdit> edits;
    for(const auto& entry:owners){int cell=entry.first;auto offset=g.origin(cell);auto box=square();for(auto& v:box)v=v+offset;auto rectangle=polygon(box);std::array<std::vector<Cell>,3> removed;
        for(int k:entry.second){auto& body=bodies[k];auto part=intersection(body.cut.get(),rectangle.get());double amount=0;GEOSArea_r(context,part.get(),&amount);if(amount<1e-15)continue;
            auto local=owned(GEOSGeom_transformXY_r(context,part.get(),[](double* x,double* y,void* pointer){auto o=static_cast<Point*>(pointer);*x-=o->x;*y-=o->y;return 1;},&offset));
            std::vector<double> weights(body.parts.size());std::vector<Geos> contributions;double sum=0;for(std::size_t i=0;i<body.parts.size();i++){auto contribution=intersection(part.get(),body.parts[i].cut.get());GEOSArea_r(context,contribution.get(),&weights[i]);sum+=weights[i];contributions.push_back(std::move(contribution));}
            Point normal{};for(std::size_t i=0;i<weights.size();i++)if(weights[i]>0&&sum>0){const auto& flux=fluxes[body.parts[i].flux];double share=weights[i]/sum;normal=normal+flux.normal*share;
                if(!growing&&flux.region>=0){
                    auto supplied=owned(GEOSGeom_transformXY_r(context,contributions[i].get(),[](double* x,double* y,void* pointer){auto o=static_cast<Point*>(pointer);*x-=o->x;*y-=o->y;return 1;},&offset));
                    CreatedWater water{flux.region,cell,{},1};unpackPolygon(supplied.get(),-1,water.shape);creation.push_back(std::move(water));
                }}
            movement[g.parent(cell)]=movement[g.parent(cell)]+normal*(-sign*amount*(1-g.epsilon)/(g.ratio*g.ratio));Cell coordinates;unpackPolygon(local.get(),body.mineral,coordinates,false);removed[body.mineral].push_back(std::move(coordinates));
        }
        if(std::all_of(removed.begin(),removed.end(),[](const auto& a){return a.empty();}))continue;
        edits.push_back({cell,std::move(removed)});
    }
    auto pruned=commitCellEdits(g.cells,edits,growing);for(const auto& values:pruned)for(int m=0;m<3;m++)deferred[m]+=values[m]*(1-g.epsilon)/(g.ratio*g.ratio);
    return !edits.empty();
}
}
