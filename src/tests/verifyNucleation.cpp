// SPDX-License-Identifier: GPL-3.0-or-later
#include "MineralInterface.hpp"
#include <iostream>
#include <iomanip>
using namespace saltgeom;
using namespace saltchem;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void close(double a,double b,double relative,const char* message){
    require(std::abs(a-b)<=relative*std::max({std::abs(a),std::abs(b),1e-100}),message);
}
std::array<double,8> conserved(const Inventory& state){
    std::array<double,8> value{};for(int j=0;j<7;j++)value[j]=state.aqueous[j];value[7]=state.waterEquivalent;
    for(int j=0;j<7;j++){
        for(int k=0;k<7;k++)value[k]+=state.mineral[j]*stoichiometry[j][k];
        value[7]+=state.mineral[j]*hydration[j]*waterMolarMass;
    }return value;
}
int main(){try{
    constexpr double cellVolume=1e-15,R=1e-7,supportArea=1e-10;
    std::cout<<std::setprecision(16)<<"test,value\n";
    for(double theta:{30.,60.,90.,120.,150.,180.}){
        const double cosine=std::cos(theta*pi/180),h=1-cosine;
        const double V=pi/3*R*R*R*(2-3*cosine+cosine*cosine*cosine);
        Population population;population.fraction=V/cellVolume;population.count=1;
        population.heightFactor=capHeightFactor(theta);population.supportArea=supportArea;
        close(radius(population,cellVolume),R,2e-14,"Spherical-cap curvature radius/volume mismatch");
        close(outerArea(population,cellVolume),2*pi*R*R*h,3e-14,"Spherical-cap outer area mismatch");
        close(footprintArea(population,cellVolume),pi*R*R*(1-cosine*cosine),3e-14,"Spherical-cap footprint mismatch");
        close(capHeight(population,cellVolume),R*h,3e-14,"Spherical-cap height mismatch");
        SurfaceState surface;surface.population[2]=population;MineralVector bare{};bare[2]=supportArea;
        const auto measured=measure(surface,bare,cellVolume);
        close(measured.capArea[0],2*pi*R*R*h,3e-14,"Isolated spherical-cap reacting area mismatch");
        close(measured.coverage[2],pi*R*R*(1-cosine*cosine)/supportArea,3e-14,"Isolated spherical-cap coverage mismatch");
    }
    std::cout<<"contact_angle_volume_footprint_and_area,1\n";
    {
        MineralVector volumes{};volumes.fill(2.7e-5);
        Inventory cell;cell.waterEquivalent=1000;cell.aqueous.fill(1000);cell.mineral[2]=cell.mineral[3]=.1/volumes[2];
        NucleationLaw first{0,2,1e12,0,45,1e-8},second{0,3,1e12,0,120,7e-9};
        auto expected=cell;SurfaceState surface;
        double birthMoles[2]{};
        for(int i=0;i<2;i++){auto copy=cell;SurfaceState single;const auto& law=i?second:first;
            require(addEmbryo(copy,single,law,cellVolume,volumes,1e-5,0,supportArea),"Birth fixture failed");birthMoles[i]=copy.mineral[0];}
        for(int i=0;i<20000;i++){
            const int kind=i%2;const auto& law=kind?second:first;
            require(addEmbryo(cell,surface,law,cellVolume,volumes,1e-5,0,supportArea),"Repeated small births exhausted an unexpected resource");
            expected.mineral[0]+=birthMoles[kind];expected.aqueous[Na]-=birthMoles[kind];expected.aqueous[Cl]-=birthMoles[kind];
            close(subgrid(surface)[0],cell.mineral[0]*volumes[0],4*std::numeric_limits<double>::epsilon(),"Pure cap fraction drift accumulated across births");
        }
        require(cell.aqueous==expected.aqueous&&cell.mineral==expected.mineral&&cell.waterEquivalent==expected.waterEquivalent,
            "Pure-cap geometric alignment modified conserved chemistry");
        std::cout<<"twenty_thousand_births_geometry_alignment,1\n";
    }
    {
        constexpr double footprint=0.25,dt=.2;
        NucleationLaw law{0,2,1e12,0,90,std::sqrt(footprint*supportArea/pi)};
        MineralVector volumes{},si{},bare{};volumes.fill(2.7e-5);si.fill(1);bare[2]=supportArea;
        Inventory initial;initial.waterEquivalent=1e6;initial.aqueous.fill(1e6);initial.mineral[2]=.1/volumes[2];
        SurfaceState whole,split;auto a=initial,b=initial;
        constexpr std::uint64_t seed=73,key=91;
        int expected=0;double elapsed=0;
        while(true){
            const double wait=exponential(seed,key,2,expected)/(law.prefactor*supportArea*std::pow(1-footprint,expected));
            if(elapsed+wait>dt)break;elapsed+=wait;expected++;
        }
        nucleate(a,whole,{law},si,bare,dt,cellVolume,volumes,1e-5,0,seed,key);
        nucleate(b,split,{law},si,bare,.07,cellVolume,volumes,1e-5,0,seed,key);
        nucleate(b,split,{law},si,bare,dt-.07,cellVolume,volumes,1e-5,0,seed,key);
        require(whole.births[0]==expected&&whole.births==split.births,"Nucleation clock ignored coverage changes during the same interval");
        close(whole.population[2].hazard,split.population[2].hazard,3e-13,"Event waiting-time clock depends on interval splitting");
        NucleationLaw second=law;second.product=1;second.radius*=.7;
        SurfaceState forward,reverse;auto c=initial,d=initial;
        nucleate(c,forward,{law,second},si,bare,dt,cellVolume,volumes,1e-5,0,seed,key);
        nucleate(d,reverse,{second,law},si,bare,dt,cellVolume,volumes,1e-5,0,seed,key);
        require(forward.births==reverse.births&&c.aqueous==d.aqueous&&c.mineral==d.mineral,"Competing nucleation laws depend on dictionary order");
        const auto before=conserved(initial),after=conserved(c);for(int j=0;j<8;j++)close(before[j],after[j],3e-14,"Competing nucleation lost stoichiometric mass");
        std::cout<<"coverage_dependent_event_clock_births,"<<expected<<"\ncompeting_law_order_independence,1\n";
    }
    {
        SurfaceState population;MineralVector before{},after{},bare{};
        const double angles[2]={30,120},radii[2]={1e-7,2e-7},increment=1e-9;
        for(int i=0;i<2;i++){
            auto& cap=population.population[2+i];cap.heightFactor=capHeightFactor(angles[i]);cap.count=1;cap.supportArea=supportArea;
            cap.fraction=capVolumeFactor(cap.heightFactor)*std::pow(radii[i],3)/cellVolume;
            before[0]+=cap.fraction;before[2+i]=.2;bare[2+i]=supportArea;
            const double finalRadius=radii[i]+2*increment/(cap.heightFactor*(3-cap.heightFactor));
            after[0]+=capVolumeFactor(cap.heightFactor)*std::pow(finalRadius,3)/cellVolume;
        }
        after[2]=after[3]=.2;const auto measured=measure(population,bare,cellVolume);
        grow(population,before,after,measured,measured.exposed,cellVolume,1.,{},bare);
        for(int i=0;i<2;i++){
            const double h=capHeightFactor(angles[i]);
            close(radius(population.population[2+i],cellVolume),radii[i]+2*increment/(h*(3-h)),5e-13,"Mixed-angle normal growth violates dV=A*dl");
        }
        close(subgrid(population)[0],after[0],2e-15,"Mixed-angle growth lost volume or created a resolved residue");
        std::cout<<"mixed_angle_conservative_normal_growth,1\n";
    }
    // Numerical fixtures below are deliberately synthetic, not fitted salt
    // interfacial energies. They verify units and formula-unit conventions.
    NucleationLaw cnt{0,2,1e13,0,60,1e-12};cnt.physicalCNT=true;cnt.interfacialEnergy=.02;
    cnt.validate();constexpr double Vm=2.7e-5,T=298.15;
    const double omega=2,si=std::log10(omega),logOmega=std::log(omega);
    const double thermal=boltzmannConstant*T,moleculeVolume=Vm/avogadroConstant;
    const double shape=(2+.5)*std::pow(1-.5,2)/4;
    const double expectedRadius=2*cnt.interfacialEnergy*Vm/(molarGasConstant*T*logOmega);
    const double expectedExponent=16*pi*std::pow(cnt.interfacialEnergy,3)*moleculeVolume*moleculeVolume*shape/(3*std::pow(thermal,3)*logOmega*logOmega);
    close(cnt.criticalRadius(si,Vm,T),expectedRadius,3e-14,"CNT critical radius has incorrect thermodynamic units");
    close(cnt.barrierExponent(si,Vm,T),expectedExponent,8e-14,"CNT barrier divided by kBT has incorrect thermodynamic units");
    close(cnt.frequency(si,Vm,T),cnt.prefactor*std::exp(-expectedExponent),1e-13,"CNT rate exponent mismatch");
    const double meanIonicS=std::sqrt(omega),nu=2;
    close(cnt.criticalRadius(si,Vm,T),2*cnt.interfacialEnergy*Vm/(molarGasConstant*T*nu*std::log(meanIonicS)),3e-14,"Mean ionic supersaturation lost its stoichiometric factor");
    close(cnt.criticalRadius(si,2*Vm,T),2*expectedRadius,3e-14,"CNT radius molar-volume scaling");
    close(cnt.barrierExponent(si,Vm,2*T),expectedExponent/8,8e-14,"CNT barrier temperature scaling");
    require(cnt.frequency(0,Vm,T)==0&&cnt.frequency(-1,Vm,T)==0,"CNT nucleation below saturation");
    std::cout<<"cnt_critical_radius_m,"<<expectedRadius<<"\ncnt_barrier_over_kBT,"<<expectedExponent<<'\n';
    {
        Inventory fluid;fluid.waterEquivalent=1000;fluid.aqueous.fill(.01);fluid.mineral[2]=.1/Vm;
        const auto old=conserved(fluid);MineralVector volumes{};volumes.fill(Vm);SurfaceState state;
        require(addEmbryo(fluid,state,cnt,cellVolume,volumes,1e-5,0,supportArea,nullptr,si,T),"Physical CNT embryo unexpectedly rejected");
        close(radius(state.population[2],cellVolume),expectedRadius,8e-14,"Physical CNT created an embryo below the critical-size construction");
        require(radius(state.population[2],cellVolume)>=expectedRadius*(1-8e-14),"Physical CNT embryo is subcritical");
        const auto now=conserved(fluid);for(int j=0;j<8;j++)close(old[j],now[j],3e-14,"Spherical-cap embryo lost stoichiometric mass");
        const auto saved=fluid;auto changedAngle=cnt;changedAngle.contactAngle=90;bool rejected=false;
        try{addEmbryo(fluid,state,changedAngle,cellVolume,volumes,1e-5,0,supportArea,nullptr,si,T);}catch(const std::invalid_argument&){rejected=true;}
        require(rejected&&fluid.aqueous==saved.aqueous&&fluid.mineral==saved.mineral,"Existing cap population was silently reinterpreted at another angle");
        MineralVector saturation{},bare{};saturation[0]=si;bare[2]=supportArea;
        nucleate(fluid,state,{cnt},saturation,bare,.1,cellVolume,volumes,1e-5,0,2024,0,{},nullptr,T);
        require(state.births[0]>1,"Physical CNT law failed to generate a stochastic population");
        const auto final=conserved(fluid);for(int j=0;j<8;j++)close(old[j],final[j],3e-14,"Physical CNT population lost stoichiometric mass");
        std::cout<<"cnt_conservative_births,"<<state.births[0]<<'\n';
    }
    {
        NucleationLaw legacy{0,2,1e12,1,45,1e-8};legacy.validate();
        const double c=std::cos(45*pi/180),expected=1e12*std::exp(-(2+c)*(1-c)*(1-c)/(4*logOmega*logOmega));
        close(legacy.frequency(si),expected,2e-15,"Legacy dimensionless-barrier mode changed its rate convention");
        require(std::isnan(legacy.criticalRadius(si,Vm,T)),"Legacy mode invented a critical-radius calibration");
        auto missing=legacy;missing.physicalCNT=true;bool rejected=false;
        try{missing.validate();}catch(const std::invalid_argument&){rejected=true;}
        require(rejected,"Physical CNT accepted an unprovided interfacial energy");
        rejected=false;try{cnt.frequency(si);}catch(const std::invalid_argument&){rejected=true;}
        require(rejected,"Physical CNT silently invented temperature or molar volume");
        std::cout<<"explicit_legacy_and_uncalibrated_parameter_guards,1\nall_passed,1\n";
    }
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
