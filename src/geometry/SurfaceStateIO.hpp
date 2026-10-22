// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "MineralInterface.hpp"
#include <istream>
#include <ostream>
#include <type_traits>

namespace saltgeom {
// Version 6 used this fixed, trivially-copyable layout. Keep its reader so
// existing runs remain restartable after compacting in-memory populations.
struct LegacySurfaceState {
    std::array<PopulationMoments,49> population{};
    std::array<Point,7> inheritedNormal{};
    MineralVector births{},promotions{};
};
template<class T> inline void writeSurfaceValue(std::ostream& out,const T& value){
    static_assert(std::is_trivially_copyable<T>::value,"Checkpoint values must not contain pointers");
    out.write(reinterpret_cast<const char*>(&value),sizeof(value));
}
template<class T> inline void readSurfaceValue(std::istream& in,T& value){
    static_assert(std::is_trivially_copyable<T>::value,"Checkpoint values must not contain pointers");
    in.read(reinterpret_cast<char*>(&value),sizeof(value));
    if(!in)throw std::runtime_error("Truncated mineral surface checkpoint");
}
inline std::array<double,13> populationValues(const Population& p){
    return {{p.fraction,p.count,p.hazard,p.threshold,p.supportArea,p.retreat,p.heightFactor,
        p.interfacialEnergy,p.localShare,p.partitionSupport,p.birthNormal.x,p.birthNormal.y,p.filmArea}};
}
inline void validatePopulation(const Population& p){
    for(double value:populationValues(p))if(!std::isfinite(value))
        throw std::runtime_error("Nonfinite mineral population checkpoint");
    if(p.fraction<0||p.count<0||p.hazard<0||p.threshold<0||p.supportArea<0||p.retreat<0||
       p.heightFactor<=0||p.heightFactor>2||p.interfacialEnergy<0||p.localShare<0||p.localShare>1||
       p.partitionSupport< -1||p.filmArea<0||(p.filmArea>0&&p.interfacialEnergy>0)||(p.fraction>0&&p.count==0))
        throw std::runtime_error("Invalid mineral population checkpoint");
}
inline void writeSurfaceState(std::ostream& out,const SurfaceState& state){
    writeSurfaceValue(out,state.population.mask());
    writeSurfaceValue(out,state.inheritedNormal);
    writeSurfaceValue(out,state.births);writeSurfaceValue(out,state.promotions);
    for(int pair=0;pair<49;pair++)if(state.population.mask()&(std::uint64_t(1)<<pair))
        writeSurfaceValue(out,state.population[pair]);
    if(!out)throw std::runtime_error("Failed to write mineral surface checkpoint");
}
inline SurfaceState readSurfaceState(std::istream& in,bool legacy=false,bool oldPopulation=false){
    SurfaceState state;
    if(legacy){
        LegacySurfaceState previous;readSurfaceValue(in,previous);
        state.inheritedNormal=previous.inheritedNormal;state.births=previous.births;state.promotions=previous.promotions;
        const auto defaults=populationValues(Population{});
        for(int pair=0;pair<49;pair++){
            Population p;static_cast<PopulationMoments&>(p)=previous.population[pair];validatePopulation(p);
            if(p.events||populationValues(p)!=defaults)state.population[pair]=p;
        }
    }else{
        std::uint64_t mask=0;readSurfaceValue(in,mask);
        if(mask>>49)throw std::runtime_error("Invalid mineral population checkpoint mask");
        readSurfaceValue(in,state.inheritedNormal);
        readSurfaceValue(in,state.births);readSurfaceValue(in,state.promotions);
        for(int pair=0;pair<49;pair++)if(mask&(std::uint64_t(1)<<pair)){
            Population p;
            if(oldPopulation)readSurfaceValue(in,static_cast<PopulationMoments&>(p));else readSurfaceValue(in,p);
            validatePopulation(p);state.population[pair]=p;
        }
    }
    for(const auto& normal:state.inheritedNormal)if(!std::isfinite(normal.x)||!std::isfinite(normal.y))
        throw std::runtime_error("Nonfinite mineral normal checkpoint");
    for(const auto* values:{&state.births,&state.promotions})for(double value:*values)
        if(!std::isfinite(value)||value<0)throw std::runtime_error("Invalid mineral event checkpoint");
    return state;
}
}
