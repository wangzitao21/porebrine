// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
namespace saltdual {
// Below one thousandth of a fine cell, use the simultaneous contact solve.
// This coalesces geometrically indistinguishable events without moving any
// boundary or relaxing either mineral's conservative volume constraint.
constexpr double interfaceContactTolerance=1e-3;
inline std::array<double,4> geometryBounds(const GEOSGeometry* shape){std::array<double,4> box{};if(!GEOSGeom_getExtent_r(booleanContext().handle,shape,&box[0],&box[1],&box[2],&box[3]))return {1,1,0,0};return box;}
inline bool boundsOverlap(const std::array<double,4>& a,const std::array<double,4>& b){return a[0]<b[2]&&b[0]<a[2]&&a[1]<b[3]&&b[1]<a[3];}
inline double frontDistance(Point p,const std::vector<std::pair<Point,Point>>& segments){
    double distance=INFINITY;
    for(auto [a,b]:segments){auto v=b-a;double vv=dot(v,v);double t=vv>0?std::clamp(dot(p-a,v)/vv,0.,1.):0;distance=std::min(distance,length(p-(a+v*t)));}
    return distance;
}
inline Poly clipScores(const Poly& polygon,const std::vector<double>& score){
    Poly result;
    for(std::size_t i=0;i<polygon.size();i++){std::size_t j=(i+1)%polygon.size();auto a=polygon[i],b=polygon[j];double fa=score[i],fb=score[j];if(fa<=0)result.push_back(a);if((fa<0&&fb>0)||(fa>0&&fb<0))result.push_back(a+(b-a)*(fa/(fa-fb)));}
    return clean(std::move(result));
}
inline long double segmentDistanceSquared(Point p,Point a,Point b){
    long double vx=static_cast<long double>(b.x)-a.x,vy=static_cast<long double>(b.y)-a.y,px=static_cast<long double>(p.x)-a.x,py=static_cast<long double>(p.y)-a.y;long double vv=vx*vx+vy*vy,t=vv>0?std::clamp((px*vx+py*vy)/vv,0.L,1.L):0;long double x=px-t*vx,y=py-t*vy;return x*x+y*y;
}
inline long double frontDistanceSquared(Point p,const std::vector<std::pair<Point,Point>>& segments){
    long double result=std::numeric_limits<long double>::infinity();for(auto [a,b]:segments)result=std::min(result,segmentDistanceSquared(p,a,b));return result;
}
template<class Candidate>
std::pair<Geos,Geos> divideArrival(const GEOSGeometry* overlap,const Candidate& originalA,const Candidate& originalB){
    auto c=booleanContext().handle;auto triangles=owned(GEOSConstrainedDelaunayTriangulation_r(c,overlap));Cell aPieces,bPieces;
    double xmin,ymin,xmax,ymax;GEOSGeom_getExtent_r(c,overlap,&xmin,&ymin,&xmax,&ymax);
    struct Front {int mineral;double distance;std::vector<std::pair<Point,Point>> segments;};
    auto nearby=[&](const Candidate& input){Front front{input.mineral,input.distance,{}};double reach=input.distance+2e-10;
        for(auto segment:input.segments){auto [u,v]=segment;if(std::max(u.x,v.x)<xmin-reach||std::min(u.x,v.x)>xmax+reach||std::max(u.y,v.y)<ymin-reach||std::min(u.y,v.y)>ymax+reach)continue;front.segments.push_back(segment);}return front;};
    // Only fronts that can reach this overlap enter its distance and endpoint
    // queries. Distant crystal edges cannot change the arrival partition.
    const auto a=nearby(originalA),b=nearby(originalB);
    const BoundaryEdgeIndex aIndex(a.segments),bIndex(b.segments);
    auto score=[&](Point p){long double x=aIndex.nearestSquared(p,[&](std::size_t k){return segmentDistanceSquared(p,a.segments[k].first,a.segments[k].second);}),y=bIndex.nearestSquared(p,[&](std::size_t k){return segmentDistanceSquared(p,b.segments[k].first,b.segments[k].second);});if(a.mineral!=b.mineral){x/=static_cast<long double>(a.distance)*a.distance;y/=static_cast<long double>(b.distance)*b.distance;}return double(x-y);};
    std::vector<std::pair<Point,double>> projectionLines;
    auto projections=[&](const auto& segments,const auto& index){std::vector<std::size_t> nearby;for(auto [u,v]:segments){auto tangent=v-u;double l=length(tangent);if(l<=tiny)continue;tangent=tangent*(1/l);if(tangent.x<0||(tangent.x==0&&tangent.y<0))tangent=tangent*(-1);
        for(auto p:{u,v}){int degree=0;index.candidates(p,p,nearby);for(auto k:nearby){auto [a,b]=segments[k];degree+=length(p-a)<=2e-10;degree+=length(p-b)<=2e-10;if(degree>1)break;}if(degree>1)continue;
            double d=dot(tangent,p);bool duplicate=false;for(const auto& line:projectionLines)if(length(tangent-line.first)<1e-12&&std::abs(d-line.second)<tiny){duplicate=true;break;}if(!duplicate)projectionLines.push_back({tangent,d});}
    }};projections(a.segments,aIndex);projections(b.segments,bIndex);
    auto assign=[&](auto&& self,Poly triangle,int level)->void{
        std::vector<double> values;double scale=0;for(auto p:triangle){values.push_back(score(p));scale=std::max(scale,std::abs(values.back()));}
        Point centre=(triangle[0]+triangle[1]+triangle[2])*(1./3);double middle=score(centre),predicted=(values[0]+values[1]+values[2])/3;scale=std::max(scale,std::abs(middle));
        bool negative=middle<0,positive=middle>0;for(double v:values){negative|=v<0;positive|=v>0;}
        for(int i=0;i<3;i++){double v=score((triangle[i]+triangle[(i+1)%3])*.5);negative|=v<0;positive|=v>0;scale=std::max(scale,std::abs(v));}
        if(scale==0){
            Point origin{},normal{};double closest=INFINITY;
            for(auto [a0,a1]:a.segments)for(auto [b0,b1]:b.segments)for(int i=0;i<2;i++)for(int j=0;j<2;j++){
                Point u=i?a1:a0,v=j?b1:b0;if(length(u-v)>2e-10)continue;
                auto da=(i?a0:a1)-u,db=(j?b0:b1)-v;double la=length(da),lb=length(db);if(la<=tiny||lb<=tiny)continue;
                Point n=da*(1/la)-db*(1/lb);double d=length((u+v)*.5-centre);if(length(n)>tiny&&d<closest){closest=d;origin=(u+v)*.5;normal=n;}
            }
            if(!std::isfinite(closest))throw StepFailure("Coincident interface fronts require distinct geometric support");
            for(std::size_t k=0;k<triangle.size();k++)values[k]=-dot(normal,triangle[k]-origin);
            auto ap=clipScores(triangle,values);for(double& v:values)v=-v;auto bp=clipScores(triangle,values);
            if(!ap.empty())aPieces.push_back({0,std::move(ap),{}});if(!bp.empty())bPieces.push_back({0,std::move(bp),{}});return;
        }
        if(!positive){aPieces.push_back({0,std::move(triangle),{}});return;}
        if(!negative){bPieces.push_back({0,std::move(triangle),{}});return;}
        if(level<4&&std::abs(middle-predicted)>1e-4*std::max(scale,1e-24)){
            auto ab=(triangle[0]+triangle[1])*.5,bc=(triangle[1]+triangle[2])*.5,ca=(triangle[2]+triangle[0])*.5;
            self(self,{triangle[0],ab,ca},level+1);self(self,{ab,triangle[1],bc},level+1);self(self,{ca,bc,triangle[2]},level+1);self(self,{ab,bc,ca},level+1);return;
        }
        auto ap=clipScores(triangle,values);for(double& v:values)v=-v;auto bp=clipScores(triangle,values);
        if(!ap.empty())aPieces.push_back({0,std::move(ap),{}});if(!bp.empty())bPieces.push_back({0,std::move(bp),{}});
    };
    for(int i=0;i<GEOSGetNumGeometries_r(c,triangles.get());i++){
        auto piece=GEOSGetGeometryN_r(c,triangles.get(),i);auto ring=GEOSGetExteriorRing_r(c,piece);if(!ring)continue;auto triangle=readRing(ring,false);if(triangle.size()!=3)continue;
        // Insert each open patch endpoint's projection boundary before triangulating.
        // A short patch between two long patches then retains its narrow
        // normal column even when all original polygon vertices lie outside it.
        std::vector<Poly> pieces{triangle};
        for(auto [normal,d]:projectionLines){std::vector<Poly> split;
            for(auto& polygon:pieces){double lo=INFINITY,hi=-INFINITY;for(auto p:polygon){double v=dot(normal,p)-d;lo=std::min(lo,v);hi=std::max(hi,v);}
                if(lo< -tiny&&hi>tiny){auto left=clip(polygon,normal,d),right=clip(polygon,normal*(-1),-d);if(!left.empty())split.push_back(std::move(left));if(!right.empty())split.push_back(std::move(right));}
                else split.push_back(std::move(polygon));
            }pieces=std::move(split);
        }
        for(auto& polygon:pieces){Point centre{};for(auto p:polygon)centre=centre+p;centre=centre*(1./polygon.size());
            for(std::size_t k=0;k<polygon.size();k++){Poly triangle{centre,polygon[k],polygon[(k+1)%polygon.size()]};if(area(triangle)>1e-16)assign(assign,std::move(triangle),0);}
        }
    }
    return {unionPieces(aPieces),unionPieces(bPieces)};
}
template<class Candidates>
std::vector<Geos> competingCuts(const Candidates& candidates){
    auto c=booleanContext().handle;std::vector<Geos> cuts;for(const auto& candidate:candidates)cuts.push_back(owned(GEOSGeom_clone_r(c,candidate.cut.get())));
    std::vector<std::array<double,4>> boxes;for(const auto& candidate:candidates)boxes.push_back(geometryBounds(candidate.cut.get()));
    for(std::size_t i=0;i<candidates.size();i++)for(std::size_t j=i+1;j<candidates.size();j++){
        if(candidates[i].from!=candidates[j].from||!boundsOverlap(boxes[i],boxes[j]))continue;auto overlap=intersection(candidates[i].cut.get(),candidates[j].cut.get());double amount=0;GEOSArea_r(c,overlap.get(),&amount);if(amount<=1e-10)continue;
        auto parts=divideArrival(overlap.get(),candidates[i],candidates[j]);cuts[i]=difference(cuts[i].get(),parts.second.get());cuts[j]=difference(cuts[j].get(),parts.first.get());
    }return cuts;
}
}
