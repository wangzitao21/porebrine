// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "BrineSystem.hpp"
#include <algorithm>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

namespace saltgeom {
constexpr double pi=3.14159265358979323846;
constexpr double boltzmannConstant=1.380649e-23,avogadroConstant=6.02214076e23;
constexpr double molarGasConstant=boltzmannConstant*avogadroConstant;
inline double capHeightFactor(double contactAngle){
    if(contactAngle==90)return 1;
    if(contactAngle==180)return 2;
    const double sine=std::sin(contactAngle*pi/360);
    return 2*sine*sine; // h/R = 1-cos(theta), stable near complete wetting.
}
inline double capVolumeFactor(double heightFactor){return pi/3*heightFactor*heightFactor*(3-heightFactor);}
inline double capOuterAreaFactor(double heightFactor){return 2*pi*heightFactor;}
inline double capFootprintFactor(double heightFactor){return pi*heightFactor*(2-heightFactor);}
using saltchem::MineralVector;
struct Point {double x=0,y=0;};
inline Point operator+(Point a,Point b){return {a.x+b.x,a.y+b.y};}
inline Point operator-(Point a,Point b){return {a.x-b.x,a.y-b.y};}
inline Point operator*(Point a,double b){return {a.x*b,a.y*b};}
inline double dot(Point a,Point b){return a.x*b.x+a.y*b.y;}
inline double cross(Point a,Point b){return a.x*b.y-a.y*b.x;}
inline double length(Point a){return std::hypot(a.x,a.y);}
// Convex polygons need at most 4 + 8 vertices after the material cuts.
struct Polygon {
    std::array<Point,16> vertex{};int n=0;
    void add(Point p){if(n==16)throw std::runtime_error("PLIC polygon capacity exceeded");vertex[n++]=p;}
    double area()const {double a=0;for(int i=0;i<n;i++)a+=cross(vertex[i],vertex[(i+1)%n]);return std::abs(a)*.5;}
};
inline Polygon box(double dx,double dy){Polygon p;for(Point v: {Point{0,0},{dx,0},{dx,dy},{0,dy}})p.add(v);return p;}
inline Polygon clip(const Polygon& p,Point normal,double offset){
    Polygon out;
    for(int i=0;i<p.n;i++){
        const Point a=p.vertex[i],b=p.vertex[(i+1)%p.n];
        const double da=dot(normal,a)-offset,db=dot(normal,b)-offset;
        if(da<=0)out.add(a);
        if((da<0&&db>0)||(da>0&&db<0))out.add(a+(b-a)*(da/(da-db)));
    }return out;
}
inline double overlap(Point a,Point b,Point c,Point d,double tolerance){
    const double l=length(b-a);if(l<=tolerance)return 0;
    const Point t=(b-a)*(1/l);
    if(std::abs(cross(t,c-a))>tolerance||std::abs(cross(t,d-a))>tolerance)return 0;
    const double u=dot(c-a,t),v=dot(d-a,t);
    return std::max(0.,std::min(l,std::max(u,v))-std::max(0.,std::min(u,v)));
}
inline double contact(const Polygon& a,const Polygon& b,double tolerance){
    double total=0;for(int i=0;i<a.n;i++)for(int j=0;j<b.n;j++)
        total+=overlap(a.vertex[i],a.vertex[(i+1)%a.n],b.vertex[j],b.vertex[(j+1)%b.n],tolerance);
    return total;
}
struct Interval {double lo=0,hi=0;};
inline Interval aperture(const Polygon& p,int side,double dx,double dy){
    const bool x=side<2;const double wall=(side&1)?(x?dx:dy):0,tol=1e-9*std::min(dx,dy);
    Interval out{ x?dy:dx,0 };
    for(int i=0;i<p.n;i++){
        const auto a=p.vertex[i],b=p.vertex[(i+1)%p.n];
        if(std::abs((x?a.x:a.y)-wall)<tol&&std::abs((x?b.x:b.y)-wall)<tol){
            out.lo=std::min({out.lo,x?a.y:a.x,x?b.y:b.x});out.hi=std::max({out.hi,x?a.y:a.x,x?b.y:b.x});
        }
    }if(out.hi<out.lo)return {};return out;
}
inline double overlap(Interval a,Interval b){return std::max(0.,std::min(a.hi,b.hi)-std::max(a.lo,b.lo));}
// Local height candidates fitted to the nine volume fractions. Unlike a
// two-neighbour gradient, a planar cut remains a planar cut at every angle.
// Only cut cells pay for the small geometric fit; no global interface solve.
inline Point heightFitNormal(const std::array<double,9>& fraction,double dx,double dy){
    std::array<double,3> columns{},rows{};
    for(int y=0;y<3;y++)for(int x=0;x<3;x++){columns[x]+=fraction[3*y+x];rows[y]+=fraction[3*y+x];}
    const double gx=columns[0]-columns[2],gy=rows[0]-rows[2];
    Point best{gx/dx,gy/dy};double bestError=std::numeric_limits<double>::infinity();
    if(fraction[4]<=1e-12||fraction[4]>=1-1e-12)return length(best)>0?best*(1/length(best)):Point{};
    const auto square=box(dx,dy);
    auto fit=[&](Point normal){
        if(length(normal)==0)return;
        normal=normal*(1/length(normal));
        double lo=std::min(0.,normal.x*dx)+std::min(0.,normal.y*dy);
        double hi=std::max(0.,normal.x*dx)+std::max(0.,normal.y*dy);
        for(int iter=0;iter<44;iter++){const double mid=(lo+hi)/2;
            if(clip(square,normal,mid).area()<fraction[4]*dx*dy)lo=mid;else hi=mid;}
        const double offset=(lo+hi)/2;double error=0;
        for(int y=0;y<3;y++)for(int x=0;x<3;x++){
            const double predicted=clip(square,normal,offset-normal.x*(x-1)*dx-normal.y*(y-1)*dy).area()/(dx*dy);
            const double diff=predicted-fraction[3*y+x];error+=diff*diff;
        }
        if(error<bestError){bestError=error;best=normal;}
    };
    fit(best);
    for(double slope:{columns[1]-columns[0],columns[2]-columns[1],(columns[2]-columns[0])/2})
        fit({-slope*dy/dx,gy>=0?1.:-1.});
    for(double slope:{rows[1]-rows[0],rows[2]-rows[1],(rows[2]-rows[0])/2})
        fit({gx>=0?1.:-1.,-slope*dx/dy});
    return best;
}
struct Reconstruction {
    std::array<Polygon,8> solid;Polygon liquid;
    std::array<Point,7> normal{};MineralVector internalArea{};
};
// Sequential multi-material PLIC. Order is explicit because volume fractions
// alone do not determine subcell topology. Carnallite is the primary interior;
// secondary phases occupy the remaining polygon. Subgrid caps are excluded.
inline Reconstruction reconstruct(const MineralVector& fraction,double inert,
    const std::array<std::array<double,4>,7>& neighbour,const std::array<Point,7>& inherited,
    double dx,double dy,double depth,double epsilonMin){
    Reconstruction r;r.liquid=box(dx,dy);const double cellArea=dx*dy;
    const int order[8]={7,2,0,1,3,4,5,6};
    for(int j:order){
        const double f=j==7?inert:fraction[j];if(f<=0)continue;
        const double target=std::min(r.liquid.area(),f*cellArea/(1-epsilonMin));
        Point normal=j==7?Point{1,0}:Point{(neighbour[j][0]-neighbour[j][1])/(2*dx),(neighbour[j][2]-neighbour[j][3])/(2*dy)};
        if(j!=7&&length(inherited[j])>0)normal=inherited[j];
        if(length(normal)<1e-12)normal={1,0};else normal=normal*(1/length(normal));
        if(j!=7)r.normal[j]=normal;
        if(target>=r.liquid.area()-cellArea*1e-13){r.solid[j]=r.liquid;r.liquid={};continue;}
        double lo=1e300,hi=-1e300;for(int k=0;k<r.liquid.n;k++){double d=dot(normal,r.liquid.vertex[k]);lo=std::min(lo,d);hi=std::max(hi,d);}
        for(int k=0;k<48;k++){const double mid=.5*(lo+hi);if(clip(r.liquid,normal,mid).area()<target)lo=mid;else hi=mid;}
        const double d=.5*(lo+hi);r.solid[j]=clip(r.liquid,normal,d);r.liquid=clip(r.liquid,normal*(-1),-d);
    }
    for(int j=0;j<7;j++)r.internalArea[j]=contact(r.solid[j],r.liquid,1e-9*std::min(dx,dy))*depth;
    return r;
}

// A marked Poisson process for each (product, substrate). Hazard and event
// counters are checkpointed, making retry and MPI scheduling reproducible.
inline std::uint64_t mix(std::uint64_t x){x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);}
inline double exponential(std::uint64_t seed,std::uint64_t cell,int pair,std::uint64_t event){
    const auto bits=mix(seed^mix(cell)^mix(std::uint64_t(pair))^mix(event+0xa174));
    return -std::log((double(bits>>11)+.5)/9007199254740992.);
}
// Stable pre-film layout, also used by the version 6/7 checkpoint reader.
struct PopulationMoments {
    double fraction=0,count=0,hazard=0,threshold=0;
    std::uint64_t events=0;
    // Population-mean birth support and subsequent normal recession. These
    // are ensemble moments, not coordinates of individually tracked grains.
    double supportArea=0,retreat=0;
    // Fixed contact-angle geometry for this (product, substrate) population.
    // The default h/R=1 is a hemisphere for explicitly constructed test data.
    double heightFactor=1;
    double interfacialEnergy=0; // J/m2; 0 means the empirical closure supplies no Kelvin energy.
    // Partition of material between an unresolved interfacial film and the
    // spatially transported residue. This is not a mechanical contact law.
    // -1 marks a fresh population whose support has not yet been recorded.
    double localShare=1,partitionSupport=-1;
    Point birthNormal{};
};
struct Population : PopulationMoments {
    // Positive for a planar, fixed-footprint patch; zero for spherical caps.
    // Coplanar film patches on one substrate are disjoint (a mineral mosaic).
    // Normal growth/dissolution changes thickness, not this reference area.
    double filmArea=0;
};
inline void clearPopulationMaterial(Population& p){
    const double hazard=p.hazard,threshold=p.threshold;const auto events=p.events;
    p=Population{};p.hazard=hazard;p.threshold=threshold;p.events=events;
}
// Missing pairs have the exact default state. Read-only access never allocates;
// mutable access retains waiting-time history even after the mineral vanishes.
// Pointees stay valid when another pair is inserted, as with the former array.
class PopulationTable {
    std::uint64_t present_=0;
    std::vector<std::unique_ptr<Population>> values_;
    static std::uint64_t bit(int pair){
        if(pair<0||pair>=49)throw std::out_of_range("Invalid mineral population pair");
        return std::uint64_t(1)<<pair;
    }
    std::size_t position(std::uint64_t b)const {return __builtin_popcountll(present_&(b-1));}
public:
    PopulationTable()=default;
    PopulationTable(const PopulationTable& other):present_(other.present_){
        values_.reserve(other.values_.size());
        for(const auto& p:other.values_)values_.emplace_back(new Population(*p));
    }
    PopulationTable(PopulationTable&& other)noexcept:
        present_(std::exchange(other.present_,0)),values_(std::move(other.values_)){}
    PopulationTable& operator=(PopulationTable other)noexcept{
        std::swap(present_,other.present_);values_.swap(other.values_);return *this;
    }
    bool empty()const {return present_==0;}
    std::size_t size()const {return values_.size();}
    std::uint64_t mask()const {return present_;}
    Population* find(int pair){
        const auto b=bit(pair);return present_&b?values_[position(b)].get():nullptr;
    }
    const Population& get(int pair)const{
        const auto b=bit(pair);static const Population absent{};
        return present_&b?*values_[position(b)]:absent;
    }
    const Population& operator[](int pair)const {return get(pair);}
    Population& operator[](int pair){
        const auto b=bit(pair);const auto at=position(b);
        if(!(present_&b)){
            values_.insert(values_.begin()+at,std::make_unique<Population>());present_|=b;
        }
        return *values_[at];
    }
};
struct SurfaceState {
    PopulationTable population;
    std::array<Point,7> inheritedNormal{};
    MineralVector births{},promotions{};
};
struct NucleationLaw {
    int product=0,substrate=2;
    double prefactor=0,barrier=0,contactAngle=90,radius=2e-8;
    bool physicalCNT=false;
    double interfacialEnergy=std::numeric_limits<double>::quiet_NaN(); // J/m2; no universal salt value.
    double shapeFactor()const {const double h=capHeightFactor(contactAngle);return h*h*(3-h)/4;}
    void validatePhysical(double molarVolume,double temperatureK)const {
        if(!std::isfinite(interfacialEnergy)||interfacialEnergy<=0||!std::isfinite(molarVolume)||molarVolume<=0||
           !std::isfinite(temperatureK)||temperatureK<=0)
            throw std::invalid_argument("Physical CNT requires positive interfacial energy, molar volume and absolute temperature");
    }
    double criticalRadius(double si,double molarVolume,double temperatureK)const {
        if(!physicalCNT)return std::numeric_limits<double>::quiet_NaN();
        validatePhysical(molarVolume,temperatureK);
        if(!(si>0))return std::numeric_limits<double>::infinity();
        const double value=std::log(2.)+std::log(interfacialEnergy)+std::log(molarVolume)
            -std::log(molarGasConstant)-std::log(temperatureK)-std::log(saltchem::ln10*si);
        return value>std::log(std::numeric_limits<double>::max())?std::numeric_limits<double>::infinity():std::exp(value);
    }
    double barrierExponent(double si,double molarVolume=std::numeric_limits<double>::quiet_NaN(),
                           double temperatureK=std::numeric_limits<double>::quiet_NaN())const {
        if(!(si>0))return std::numeric_limits<double>::infinity();
        const double logOmega=saltchem::ln10*si;
        if(!physicalCNT)return barrier*shapeFactor()/(logOmega*logOmega);
        validatePhysical(molarVolume,temperatureK);
        // SI=log10(IAP/K) is per mineral formula unit, including hydrate
        // water activity. Delta mu=kBT ln(Omega), not kBT ln(Omega^(1/nu)).
        // This is Delta G*/kBT; Delta G* itself has only (kBT)^2 below.
        const double logarithm=std::log(16*pi/3)+3*std::log(interfacialEnergy)
            +2*(std::log(molarVolume)-std::log(avogadroConstant))+std::log(shapeFactor())
            -3*(std::log(boltzmannConstant)+std::log(temperatureK))-2*std::log(logOmega);
        return logarithm>std::log(std::numeric_limits<double>::max())?std::numeric_limits<double>::infinity():std::exp(logarithm);
    }
    double embryoRadius(double si,double molarVolume,double temperatureK)const {
        if(!physicalCNT)return radius;
        if(!std::isfinite(si)||si<=0)throw std::invalid_argument("Physical CNT embryo requires a positive finite saturation index");
        const double critical=criticalRadius(si,molarVolume,temperatureK);
        return std::max(radius,std::nextafter(critical,std::numeric_limits<double>::infinity()));
    }
    double frequency(double si,double molarVolume=std::numeric_limits<double>::quiet_NaN(),
                     double temperatureK=std::numeric_limits<double>::quiet_NaN())const {
        if(!(si>0)||prefactor==0)return 0;
        return prefactor*std::exp(-barrierExponent(si,molarVolume,temperatureK));
    }
    void validate()const {
        if(product<0||product>=7||substrate<0||substrate>=7||product==substrate||!std::isfinite(prefactor)||prefactor<0||!std::isfinite(barrier)||barrier<0||!std::isfinite(contactAngle)||contactAngle<=0||contactAngle>180||!std::isfinite(radius)||radius<=0)
            throw std::invalid_argument("Invalid heterogeneous nucleation law");
        if(!(capVolumeFactor(capHeightFactor(contactAngle))>0))throw std::invalid_argument("Contact angle is below representable spherical-cap geometry");
        if(physicalCNT&&(!std::isfinite(interfacialEnergy)||interfacialEnergy<=0))
            throw std::invalid_argument("Physical CNT requires an explicitly supplied positive interfacialEnergy in J/m2");
    }
};
inline MineralVector subgrid(const SurfaceState& s){MineralVector f{};if(s.population.empty())return f;for(int j=0;j<7;j++)for(int k=0;k<7;k++)f[j]+=s.population[7*j+k].fraction;return f;}
inline bool whollySubgrid(double total,double caps){
    return total==0?caps==0:caps>0&&std::abs(total-caps)<=64*std::numeric_limits<double>::epsilon()*std::max(total,caps);
}
// Align only a phase already known to have no resolved stock. The conserved
// mineral inventory is authoritative; this removes geometry-conversion ulps
// after each birth/update, before they can accumulate across many steps.
inline void alignPureSubgridInventory(SurfaceState& surface,int mineral,double totalFraction){
    int largest=-1;double sum=0;
    for(int k=0;k<7;k++){
        const auto& p=surface.population.get(7*mineral+k);sum+=p.fraction;
        if(p.fraction>0&&(largest<0||p.fraction>surface.population.get(7*mineral+largest).fraction))largest=k;
    }
    const double budget=64*std::numeric_limits<double>::epsilon()*std::max(totalFraction,sum);
    if(!std::isfinite(totalFraction)||totalFraction<0||std::abs(totalFraction-sum)>budget)
        throw std::runtime_error("Pure subgrid geometry exceeds the single-update roundoff budget");
    if(largest<0){if(totalFraction!=0)throw std::runtime_error("No cap remains to represent pure subgrid inventory");return;}
    long double other=0;for(int k=0;k<7;k++)if(k!=largest)other+=surface.population.get(7*mineral+k).fraction;
    const double corrected=double(static_cast<long double>(totalFraction)-other);
    if(corrected<0)throw std::runtime_error("Pure subgrid alignment would consume another population");
    surface.population[7*mineral+largest].fraction=corrected;
}
inline double radius(const Population& p,double cellVolume){return p.count>0&&p.fraction>0?std::cbrt(p.fraction*cellVolume/(capVolumeFactor(p.heightFactor)*p.count)):0;}
inline double capHeight(const Population& p,double cellVolume){return p.filmArea>0?p.fraction*cellVolume/p.filmArea:p.heightFactor*radius(p,cellVolume);}
inline double outerArea(const Population& p,double cellVolume){if(p.fraction<=0)return 0;if(p.filmArea>0)return p.filmArea;const double r=radius(p,cellVolume);return capOuterAreaFactor(p.heightFactor)*r*r*p.count;}
inline double footprintArea(const Population& p,double cellVolume){if(p.fraction<=0)return 0;if(p.filmArea>0)return p.filmArea;const double r=radius(p,cellVolume);return capFootprintFactor(p.heightFactor)*r*r*p.count;}
inline double capExtent(const Population& p,double cellVolume){return radius(p,cellVolume)*std::max(p.heightFactor,std::sqrt(p.heightFactor*(2-p.heightFactor)));}
struct ProductLayer {
    bool enabled=false,crossGrid=false;
    double diffusivity=1e-9,porosity=.35,tortuosity=2,affinityConcentration=1000;
    MineralVector kinetic{};
    void validate()const {
        if(!std::isfinite(diffusivity)||diffusivity<0||!std::isfinite(porosity)||porosity<0||porosity>=1||!std::isfinite(tortuosity)||tortuosity<1||!std::isfinite(affinityConcentration)||affinityConcentration<=0)
            throw std::invalid_argument("Invalid porous product-layer properties");
        for(double k:kinetic)if(!std::isfinite(k)||k<0)throw std::invalid_argument("Invalid product-layer kinetic coefficient");
    }
    // Linear affinity-film closure: C* converts the dimensionless TST
    // affinity into a mol/m3 concentration driving force. Only the excess
    // subgrid resistance relative to unobstructed water is added; resolved
    // mineral volume is excluded from this film calculation.
    double transmission(double thickness,int mineral)const {
        if(thickness<=0||kinetic[mineral]==0)return 1;
        if(diffusivity==0||porosity==0)return 0;
        const double excess=thickness/diffusivity*(tortuosity/porosity-1);
        return 1/(1+kinetic[mineral]*excess/affinityConcentration);
    }
};
// Additional resistivity of a retained laminar product/water region, in 1/D0
// units per unit cell length. Existing eps^2 hindrance is subtracted, not added
// for a second time. Primary/resolved solid is excluded from this region.
// Porosity is reduced if the nominal porous layer cannot fit the available
// water volume. The nonnegative excess preserves the baseline transport limit.
inline double retainedLayerResistance(double water,double product,const ProductLayer& layer){
    if(!(product>0))return 0;
    if(!(water>0))return 1e100;
    const double region=water+product,pore=water/region;
    const double occupied=std::min(region,product/(1-layer.porosity));
    const double porosity=std::max(0.,1-product/occupied),fraction=occupied/region;
    if(porosity<=0||layer.diffusivity<=0)return 1e100;
    const double series=1-fraction+fraction*layer.tortuosity/porosity;
    return region*std::max(0.,series-1/(pore*pore));
}
inline double faceDiffusion(double baseline,double excessOverDLeft,double excessOverDRight,double weight=.5){
    if(baseline<=0)return 0;
    return baseline/(1+baseline*(weight*excessOverDLeft+(1-weight)*excessOverDRight));
}
struct RetainedConductivity {
    double relative=1,relativeNormal=1,relativeParallel=1,coverage=0,porosity=0,thicknessFraction=0;
};
// Parallel uncovered water and covered columns; water and porous product
// are in series within a covered column. The occupied layer is shortened
// or its porosity reduced to fit all product, never by discarding inventory.
// coverage is supplied from the caller's directional footprint. This is a
// laminate closure, not a universal connectivity law for arbitrary grains.
inline RetainedConductivity retainedConductivity(double water,double product,double coverage,
                                                  double waterDiffusivity,const ProductLayer& layer){
    if(!std::isfinite(water)||water<0||!std::isfinite(product)||product<0||
       !std::isfinite(coverage)||coverage<0||coverage>1||!std::isfinite(waterDiffusivity)||waterDiffusivity<0)
        throw std::invalid_argument("Invalid retained-layer volume, coverage or water diffusivity");
    layer.validate();RetainedConductivity out;
    if(product==0)return out;
    if(water==0||waterDiffusivity==0){out.relative=out.relativeNormal=out.relativeParallel=0;out.coverage=1;out.thicknessFraction=1;return out;}
    const double fraction=product/(water+product);
    out.coverage=std::max(coverage,fraction);
    out.porosity=layer.enabled?layer.porosity:0;
    out.thicknessFraction=fraction/(out.coverage*(1-out.porosity));
    if(out.thicknessFraction>1){out.thicknessFraction=1;out.porosity=std::max(0.,1-fraction/out.coverage);}
    const double effective=layer.enabled?layer.diffusivity*out.porosity/layer.tortuosity:0;
    double covered=0;
    if(effective>0){const double scale=std::max(effective,waterDiffusivity),a=effective/scale,b=waterDiffusivity/scale;
        covered=a/((1-out.thicknessFraction)*a+out.thicknessFraction*b);}
    out.relative=out.relativeNormal=(1-out.coverage)+out.coverage*covered;
    const double occupied=out.coverage*out.thicknessFraction;
    const long double parallel=1-occupied+static_cast<long double>(occupied)*effective/waterDiffusivity;
    if(parallel>std::numeric_limits<double>::max())throw std::overflow_error("Retained-layer relative conductivity exceeds numerical range");
    out.relativeParallel=double(parallel);
    return out;
}
inline double partitionReference(const Population& p,double bare){
    return p.partitionSupport>=0?p.partitionSupport:(p.supportArea>0?p.supportArea:bare);
}
inline double attachedShare(const Population& p,double bare,double cellVolume){
    (void)cellVolume;
    if(!(p.fraction>0&&p.count>0&&bare>0&&p.localShare>0))return 0;
    if(p.heightFactor==2)return 0; // Tangent sphere: zero solid/substrate contact area.
    const double reference=partitionReference(p,bare);
    return reference>0?p.localShare*std::min(1.,bare/reference):0;
}
// Call after reconstructing current support and normals, before measure().
// Support contraction and a reversal relative to the birth normal retire
// material irreversibly. Increasing area cannot reattach existing residue.
// A single population cannot resolve distinct same-mineral patches; the
// normal guard rejects an obvious reversal, not every possible identity change.
inline void updatePartition(SurfaceState& s,const MineralVector& support,const std::array<Point,7>& normal){
    if(s.population.empty())return;
    for(int j=0;j<7;j++)for(int k=0;k<7;k++){
        if(!s.population.find(7*j+k))continue;
        auto& p=s.population[7*j+k];if(!(p.fraction>0&&p.count>0))continue;
        if(!std::isfinite(support[k])||support[k]<0||!std::isfinite(p.localShare)||p.localShare<0||p.localShare>1)
            throw std::invalid_argument("Invalid interfacial-film partition");
        const double currentNorm=length(normal[k]),bornNorm=length(p.birthNormal);
        if(p.partitionSupport<0)p.partitionSupport=p.supportArea>0?p.supportArea:support[k];
        if(!(support[k]>0)||p.heightFactor==2||
           (bornNorm>0&&currentNorm>0&&dot(p.birthNormal,normal[k])< -64*std::numeric_limits<double>::epsilon()*bornNorm*currentNorm)){
            p.localShare=0;p.partitionSupport=0;continue;
        }
        if(p.partitionSupport>0&&support[k]<p.partitionSupport){p.localShare*=support[k]/p.partitionSupport;p.partitionSupport=support[k];}
        if(bornNorm==0&&currentNorm>0&&p.localShare>0)p.birthNormal=normal[k]*(1/currentNorm);
    }
}
// A new embryo contributes only its new volume to the local-film account.
// Already retained material never enters that account again. The population
// retains its existing common-size/moment approximation after merging births.
inline void addPartitionVolume(Population& p,double fraction,double support){
    const double oldShare=attachedShare(p,support,1);
    const double freshShare=support>0&&p.heightFactor<2?1:0;
    p.localShare=(p.fraction*oldShare+fraction*freshShare)/(p.fraction+fraction);
    p.partitionSupport=std::max(0.,support);
    if(oldShare==0)p.birthNormal={};
}
struct SurfaceMeasure {
    MineralVector exposed{},capArea{},coverage{},footprint{},accessibility{};
    MineralVector layerCoverage{},layerThickness{},layerTransmission{},attachedCount{},attachedFraction{},residueFraction{};
    std::array<double,49> populationArea{};
};
// Random fixed-contact-angle spherical caps on an extruded 2D interface. Impingement is a
// population closure (Boolean coverage), not resolved individual crystals.
// It preserves solid volume even when footprints overlap or populations merge.
inline SurfaceMeasure measure(const SurfaceState& s,const MineralVector& bare,double cellVolume,bool blocking=true,const ProductLayer& layer={}){
    SurfaceMeasure r;
    r.accessibility.fill(1);r.layerTransmission.fill(1);
    if(s.population.empty()){r.exposed=bare;return r;}
    std::array<double,49> projected{},outer{},share{},domain{},logUncovered{},filmCovered{};
    MineralVector reference{},localVolume{};
    for(int j=0;j<7;j++)for(int k=0;k<7;k++){
        const int id=7*j+k;const auto& p=s.population[id];
        if(!(p.fraction>0&&p.count>0))continue;
        projected[id]=footprintArea(p,cellVolume);outer[id]=outerArea(p,cellVolume);r.footprint[k]+=projected[id];
        reference[k]=std::max(reference[k],p.supportArea>0?p.supportArea:bare[k]);
        share[id]=attachedShare(p,bare[k],cellVolume);
        r.attachedCount[j]+=p.count*share[id];r.attachedFraction[j]+=p.fraction*share[id];r.residueFraction[j]+=p.fraction*(1-share[id]);
        localVolume[k]+=p.fraction*share[id]*cellVolume;
        if(share[id]>0&&projected[id]>0){
            domain[id]=std::min(bare[k],partitionReference(p,bare[k]));
            const double optical=projected[id]*share[id]/domain[id];
            if(p.filmArea>0)filmCovered[id]=std::min(1.,optical);
            else logUncovered[id]=layer.enabled?-optical:
                (optical>=p.count?-1e300:p.count*std::log1p(-optical/p.count));
        }
    }
    for(int k=0;k<7;k++){
        // Nested representative support patches avoid redistributing old
        // footprints over newly exposed area. They are a low-order overlap
        // closure, not stored coordinates of individual crystal contacts.
        std::array<double,7> ends{};for(int j=0;j<7;j++)ends[j]=domain[7*j+k];std::sort(ends.begin(),ends.end());
        double start=0,coveredArea=0;
        for(double end:ends)if(end>start){double logEmpty=0,filmFraction=0;
            for(int j=0;j<7;j++)if(domain[7*j+k]>=end){logEmpty+=logUncovered[7*j+k];filmFraction+=filmCovered[7*j+k];}
            // Disjoint film patches, with random caps on the remaining area.
            const double sheet=std::min(1.,filmFraction);
            coveredArea+=(end-start)*(sheet+(1-sheet)*(-std::expm1(logEmpty)));start=end;}
        const double covered=bare[k]>0?std::clamp(coveredArea/bare[k],0.,1.):0;
        r.coverage[k]=r.layerCoverage[k]=covered;
        const double porosity=layer.enabled?layer.porosity:0;
        r.layerThickness[k]=coveredArea>0?localVolume[k]/(coveredArea*(1-porosity)):0;
        r.layerTransmission[k]=covered>0?(layer.enabled?layer.transmission(r.layerThickness[k],k):0):1;
        r.accessibility[k]=blocking?1-covered+covered*r.layerTransmission[k]:1;r.exposed[k]=bare[k]*r.accessibility[k];
        const double f=r.footprint[k],ref=reference[k];double birthCovered=0;
        if(ref>0){
            if(layer.enabled)birthCovered=-std::expm1(-f/ref);
            else {double logEmpty=0;for(int j=0;j<7;j++){const auto& p=s.population[7*j+k];if(p.count>0){const double q=projected[7*j+k]/(p.count*ref);logEmpty+=q>=1?-1e300:p.count*std::log1p(-q);}}
                birthCovered=-std::expm1(logEmpty);}
        }
        for(int j=0;j<7;j++)if(outer[7*j+k]>0){
            // Crowding of the product's outer surface is independent of the
            // bookkeeping split and never vanishes with its primary support.
            const double a=outer[7*j+k]*(s.population[7*j+k].filmArea==0&&ref>0&&f>0&&projected[7*j+k]>0&&(layer.enabled||share[7*j+k]>0)?ref*birthCovered/f:1);
            r.populationArea[7*j+k]=a;r.capArea[j]+=a;
        }
    }
    for(int j=0;j<7;j++)r.exposed[j]+=r.capArea[j];return r;
}
inline bool hasResolvedSubstrate(const saltchem::Inventory& cell,const SurfaceState& surface,
                                 int substrate,const MineralVector& molarVolume,const MineralVector* resolvedSubstrate=nullptr){
    // A fluid-side interface can be supported by a neighbouring solid owner.
    // The supplied current resolved inventory is a presence guard only; an
    // embryo never consumes the substrate mineral itself.
    if(resolvedSubstrate)return (*resolvedSubstrate)[substrate]>0;
    const double total=cell.mineral[substrate]*molarVolume[substrate];
    double caps=0;for(int k=0;k<7;k++)caps+=surface.population[7*substrate+k].fraction;
    // bare refers to resolved PLIC material, not to the separate cap population.
    // Only discard cancellation at machine precision; no physical size cutoff.
    return total-caps>64*std::numeric_limits<double>::epsilon()*std::max(total,caps);
}
inline bool addEmbryo(saltchem::Inventory& cell,SurfaceState& surface,const NucleationLaw& law,
                     double cellVolume,const MineralVector& molarVolume,double epsilonMin,double inert,double supportArea=0,
                     const MineralVector* resolvedSubstrate=nullptr,double si=std::numeric_limits<double>::quiet_NaN(),
                     double temperatureK=std::numeric_limits<double>::quiet_NaN()){
    if(!hasResolvedSubstrate(cell,surface,law.substrate,molarVolume,resolvedSubstrate))return false;
    const double heightFactor=capHeightFactor(law.contactAngle),r=law.embryoRadius(si,molarVolume[law.product],temperatureK);
    const double logVolume=std::log(capVolumeFactor(heightFactor))+3*std::log(r);
    const double logMoles=logVolume-std::log(cellVolume)-std::log(molarVolume[law.product]);
    if(logVolume-std::log(cellVolume)>=0||logMoles>std::log(std::numeric_limits<double>::max()))return false;
    const double moles=std::exp(logMoles),fraction=moles*molarVolume[law.product];
    if(!std::isfinite(fraction)||!std::isfinite(moles)||!(fraction>0&&moles>0))return false;
    auto& p=surface.population[7*law.product+law.substrate];
    if(p.filmArea>0)throw std::invalid_argument("A fixed-footprint film cannot merge new spherical-cap nuclei into the same population");
    double oldCaps=0;for(int k=0;k<7;k++)oldCaps+=surface.population.get(7*law.product+k).fraction;
    const bool purelySubgrid=whollySubgrid(cell.mineral[law.product]*molarVolume[law.product],oldCaps);
    if(p.count>0&&p.heightFactor!=heightFactor)
        throw std::invalid_argument("Cannot change contact angle of an existing mineral population");
    const double interfacialEnergy=law.physicalCNT?law.interfacialEnergy:0;
    if(p.count>0&&p.interfacialEnergy!=interfacialEnergy)
        throw std::invalid_argument("Cannot change interfacial energy of an existing mineral population");
    for(int i=0;i<7;i++)if(cell.aqueous[i]<saltchem::stoichiometry[law.product][i]*moles)return false;
    if(cell.waterEquivalent<=saltchem::hydration[law.product]*saltchem::waterMolarMass*moles)return false;
    double solid=inert;for(int j=0;j<7;j++)solid+=cell.mineral[j]*molarVolume[j];
    if(1-solid-epsilonMin<fraction)return false;
    // Same reaction extent as kinetic growth; no massless seed or hidden mass.
    for(int i=0;i<7;i++)cell.aqueous[i]-=saltchem::stoichiometry[law.product][i]*moles;
    cell.waterEquivalent-=saltchem::hydration[law.product]*saltchem::waterMolarMass*moles;
    cell.mineral[law.product]+=moles;
    p.heightFactor=heightFactor;
    p.interfacialEnergy=interfacialEnergy;
    addPartitionVolume(p,fraction,supportArea);
    p.supportArea=(p.supportArea*p.count+supportArea)/(p.count+1);p.retreat*=p.count/(p.count+1);
    p.fraction+=fraction;p.count+=1;surface.births[law.product]+=1;
    if(purelySubgrid)alignPureSubgridInventory(surface,law.product,cell.mineral[law.product]*molarVolume[law.product]);
    return true;
}
inline void nucleate(saltchem::Inventory& cell,SurfaceState& surface,const std::vector<NucleationLaw>& laws,
    const MineralVector& si,const MineralVector& bare,double dt,double cellVolume,const MineralVector& molarVolume,
    double epsilonMin,double inert,std::uint64_t seed,std::uint64_t cellKey,const ProductLayer& layer={},
    const MineralVector* resolvedSubstrate=nullptr,double temperatureK=std::numeric_limits<double>::quiet_NaN()){
    if(!std::isfinite(dt)||dt<0)throw std::invalid_argument("Invalid nucleation interval");
    if(laws.empty()||dt==0)return;
    if(laws.size()>49)throw std::invalid_argument("Too many independent nucleation laws");
    std::array<double,49> frequencies{},rates{};std::array<bool,49> blocked{};
    for(std::size_t l=0;l<laws.size();l++)
        frequencies[l]=laws[l].frequency(si[laws[l].product],molarVolume[laws[l].product],temperatureK);
    if(std::none_of(frequencies.begin(),frequencies.end(),[](double frequency){return frequency>0;}))return;
    double remaining=dt;int attempts=0;
    while(remaining>0){
        // Competing laws advance on one physical clock. Every accepted birth
        // immediately changes the uncovered area seen by all products, not
        // just the product whose event happened first.
        const auto surfaceMeasure=measure(surface,bare,cellVolume,true,layer);
        double advance=remaining;int event=-1;rates.fill(0);
        for(std::size_t l=0;l<laws.size();l++){
            const auto& law=laws[l];
            if(blocked[l]||frequencies[l]==0||!hasResolvedSubstrate(cell,surface,law.substrate,molarVolume,resolvedSubstrate))continue;
            const double available=bare[law.substrate]*(1-surfaceMeasure.coverage[law.substrate]);
            if(!(available>0))continue;
            if(available>1&&frequencies[l]>std::numeric_limits<double>::max()/available)
                throw std::runtime_error("Nucleation rate exceeds numerical range");
            const double rate=rates[l]=frequencies[l]*available;if(rate==0)continue;
            auto& p=surface.population[7*law.product+law.substrate];
            if(p.threshold==0)p.threshold=exponential(seed,cellKey,7*law.product+law.substrate,p.events);
            const double needed=std::max(0.,p.threshold-p.hazard);
            const double budget=rate>1&&remaining>std::numeric_limits<double>::max()/rate?
                std::numeric_limits<double>::infinity():rate*remaining;
            if(budget>=needed){const double wait=needed/rate;if(wait<advance||(wait==advance&&event<0)){advance=wait;event=int(l);}}
        }
        for(std::size_t l=0;l<laws.size();l++)if(rates[l]>0)
            surface.population[7*laws[l].product+laws[l].substrate].hazard+=rates[l]*advance;
        remaining=std::max(0.,remaining-advance);if(event<0)break;
        const auto& law=laws[event];auto& p=surface.population[7*law.product+law.substrate];
        p.hazard=0;p.events++;p.threshold=exponential(seed,cellKey,7*law.product+law.substrate,p.events);
        if(++attempts>10000)throw std::runtime_error("Too many nucleation events; reduce time step or prefactor");
        // No reagent or pore volume can be regained during this frozen-
        // chemistry substep. A failed birth therefore pauses its own clock.
        blocked[event]=!addEmbryo(cell,surface,law,cellVolume,molarVolume,epsilonMin,inert,bare[law.substrate],resolvedSubstrate,si[law.product],temperatureK);
    }
}
inline void promoteResolvedPopulations(SurfaceState& surface,double cellVolume,double resolution,
                                      const std::array<Point,7>& normal,bool resolve=true){
    if(!std::isfinite(cellVolume)||cellVolume<=0||!std::isfinite(resolution)||resolution<=0)
        throw std::invalid_argument("Invalid resolved cap transition scale");
    // A population has volume and area moments, but no footprint coordinates.
    // Its radius alone cannot justify replacing it by a cell-spanning plane.
    // The population closure retains these moments, including after support
    // loss; it does not claim to resolve individual crystal fronts.
    if(!resolve)return;
    for(int j=0;j<7;j++)for(int k=0;k<7;k++){
        if(!surface.population.find(7*j+k))continue;
        auto& p=surface.population[7*j+k];
        if(p.filmArea==0&&p.fraction>0&&capExtent(p,cellVolume)>=resolution){
            if(length(normal[k])>0)surface.inheritedNormal[j]=normal[k];
            surface.promotions[j]+=p.count;clearPopulationMaterial(p);
        }
    }
}
inline MineralVector inaccessibleFractions(const SurfaceState& s,const MineralVector& total,const SurfaceMeasure& measured){
    const auto caps=subgrid(s);MineralVector locked{};
    for(int j=0;j<7;j++){
        if(measured.exposed[j]<=measured.capArea[j])locked[j]=std::max(0.,total[j]-caps[j]);
        for(int k=0;k<7;k++)if(measured.populationArea[7*j+k]<=0)locked[j]+=s.population[7*j+k].fraction;
        locked[j]=std::min(total[j],locked[j]);
    }return locked;
}
inline void grow(SurfaceState& s,const MineralVector& before,const MineralVector& after,
    const SurfaceMeasure& measured,const MineralVector& totalGrowthArea,double cellVolume,double resolution,
    const std::array<Point,7>& normal,const MineralVector& bare={},const ProductLayer& layer={},double epsilonMin=0,
    const std::array<bool,7>* externalSupport=nullptr,bool updateSubstrateRetreat=true,bool resolve=true){
    const auto old=subgrid(s);MineralVector resolvedBefore{},resolvedAfter{};
    for(int j=0;j<7;j++){
        const double delta=after[j]-before[j];
        const double roundoff=64*std::numeric_limits<double>::epsilon()*std::max(before[j],after[j]);
        if(old[j]>before[j]+roundoff)throw std::runtime_error("Subgrid mineral exceeds total inventory");
        resolvedBefore[j]=std::max(0.,before[j]-old[j]);
        if(resolvedBefore[j]<=roundoff&&old[j]>0)resolvedBefore[j]=0;
        resolvedAfter[j]=resolvedBefore[j];
        if(delta==0)continue;
        const double resolvedArea=std::max(0.,(delta>0?totalGrowthArea[j]:measured.exposed[j])-measured.capArea[j]);
        std::array<double,7> inverseLength{},fraction{},area{};
        double accessible=resolvedArea>0?resolvedBefore[j]:0,slope=resolvedArea/cellVolume,maxRetreat=0;
        if(resolvedArea>0)maxRetreat=resolvedBefore[j]*cellVolume/resolvedArea;
        for(int k=0;k<7;k++){
            const auto& p=s.population.get(7*j+k);fraction[k]=p.fraction;area[k]=measured.populationArea[7*j+k];
            if(p.fraction>0&&p.count>0&&area[k]>0){
                // dV=A_eff*dl and dV/dR=3V/R for a fixed-contact-angle cap.
                // At equal exposure, equal-angle caps have the same dR.
                // Across angles dR/dl=A_outer/(3V/R), preserving the actual
                // surface-integrated reaction volume instead of assuming a
                // hemispherical radial velocity for every contact angle.
                inverseLength[k]=area[k]/(3*p.fraction*cellVolume);
                accessible+=p.fraction;slope+=area[k]/cellVolume;
                maxRetreat=std::max(maxRetreat,1/inverseLength[k]);
            }
        }
        if(delta< -accessible-roundoff)throw std::runtime_error("Reaction consumed inaccessible mineral; reduce geometry step");
        if(!(slope>0))throw std::runtime_error("Mineral changed without an accessible surface");
        auto changeAt=[&](double distance,double* derivative=nullptr){
            const double resolvedChange=resolvedArea*distance/cellVolume;
            double change=std::max(-resolvedBefore[j],resolvedChange);
            double d=resolvedChange> -resolvedBefore[j]?resolvedArea/cellVolume:0;
            for(int k=0;k<7;k++)if(inverseLength[k]>0){
                const double x=distance*inverseLength[k];
                if(x<=-1)change-=fraction[k];
                else {change+=fraction[k]*x*(3+x*(3+x));d+=area[k]/cellVolume*(1+x)*(1+x);}
            }
            if(derivative)*derivative=d;return change;
        };
        double distance=delta/slope;
        if(delta<0&&-delta>=accessible-roundoff)distance=-maxRetreat;
        else{
            double lo=delta>0?0:-maxRetreat,hi=delta>0?distance:0;
            distance=std::clamp(distance,lo,hi);
            for(int iteration=0;iteration<64;iteration++){
                double derivative=0;const double error=changeAt(distance,&derivative)-delta;
                if(std::abs(error)<=8*std::numeric_limits<double>::epsilon()*std::abs(delta))break;
                if(error>0)hi=distance;else lo=distance;
                const double candidate=derivative>0?distance-error/derivative:.5*(lo+hi);
                distance=candidate>lo&&candidate<hi?candidate:.5*(lo+hi);
            }
        }
        resolvedAfter[j]=std::max(0.,resolvedBefore[j]+resolvedArea*distance/cellVolume);
        int correction=-1;double sum=0;
        for(int k=0;k<7;k++){
            if(!s.population.find(7*j+k))continue;
            auto& p=s.population[7*j+k];
            if(inverseLength[k]>0){const double ratio=std::max(0.,1+distance*inverseLength[k]);p.fraction=fraction[k]*ratio*ratio*ratio;
                if(p.fraction>0&&(correction<0||p.fraction>s.population[7*j+correction].fraction))correction=k;}
            sum+=p.fraction;
        }
        // Close only floating-point residue on an exposed population. This
        // keeps a wholly subgrid mineral from creating a spurious PLIC film,
        // and never transfers reaction volume to buried material.
        const double residual=after[j]-resolvedAfter[j]-sum;
        if(correction>=0){
            if(std::abs(residual)>roundoff||s.population[7*j+correction].fraction+residual<0)
                throw std::runtime_error("Mineral population volume solve did not converge");
            s.population[7*j+correction].fraction+=residual;
        }
    }
    for(int j=0;j<7;j++)for(int k=0;k<7;k++){
            if(!s.population.find(7*j+k))continue;
            auto& p=s.population[7*j+k];
            if(p.fraction==0){clearPopulationMaterial(p);continue;}
            if(layer.enabled){
                if(updateSubstrateRetreat){
                    if(bare[k]>0)p.retreat+=std::max(0.,resolvedBefore[k]-resolvedAfter[k])*cellVolume/(bare[k]*(1-epsilonMin));
                    if(resolvedAfter[k]==0&&!(externalSupport&&(*externalSupport)[k]))p.retreat=std::max(p.retreat,2*capHeight(p,cellVolume));
                }
                if(length(normal[k])>0)s.inheritedNormal[j]=normal[k];
            }
            // Legacy mode resolves a population when its support disappears.
            // Porous mode retains its subgrid reacting area until resolution.
            if(resolve&&(capExtent(p,cellVolume)>=resolution||(!layer.enabled&&updateSubstrateRetreat&&resolvedAfter[k]==0&&!(externalSupport&&(*externalSupport)[k])))){
                if(length(normal[k])>0)s.inheritedNormal[j]=normal[k];
                s.promotions[j]+=p.count;clearPopulationMaterial(p);
            }
    }
}
}
