// SPDX-License-Identifier: GPL-3.0-or-later
#include "SurfaceStateIO.hpp"
#include <iostream>
#include <sstream>

using namespace saltgeom;
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
void equal(const SurfaceState& a,const SurfaceState& b){
    require(a.births==b.births&&a.promotions==b.promotions,"Lost event totals");
    for(int j=0;j<7;j++)require(a.inheritedNormal[j].x==b.inheritedNormal[j].x&&
        a.inheritedNormal[j].y==b.inheritedNormal[j].y,"Lost interface normal");
    for(int pair=0;pair<49;pair++)require(populationValues(a.population[pair])==populationValues(b.population[pair])&&
        a.population[pair].events==b.population[pair].events,"Lost population geometry or stochastic history");
}
template<class Action> void rejects(Action action){
    bool failed=false;try{action();}catch(const std::runtime_error&){failed=true;}
    require(failed,"Corrupt checkpoint was accepted");
}
int main(){try{
    require(sizeof(SurfaceState)<512,"Empty interface state is no longer compact");
    SurfaceState empty;MineralVector zero{};std::array<Point,7> normals{};
    updatePartition(empty,zero,normals);
    grow(empty,zero,zero,measure(empty,zero,1),zero,1,1,normals);
    promoteResolvedPopulations(empty,1,1,normals);
    require(subgrid(empty)==zero&&empty.population.empty(),"Absent populations were allocated by geometry updates");
    const auto bare=MineralVector{{1,2,3,4,5,6,7}};
    const auto measured=measure(empty,bare,1);
    require(measured.exposed==bare&&measured.coverage==zero,"Empty surface changed bare area");

    SurfaceState state;
    auto& cap=state.population[2];cap.fraction=.01;cap.count=4;cap.supportArea=.5;
    cap.hazard=.25;cap.threshold=.7;cap.events=3;cap.birthNormal={1,0};
    state.population[9]=state.population[2];
    require(&cap==state.population.find(2),"Inserting another pair invalidated a population reference");
    state.population[9].fraction=.02;
    auto& waiting=state.population[48];waiting.hazard=.3;waiting.threshold=.9;waiting.events=42;
    state.inheritedNormal[0]={1,0};state.births[0]=8;state.promotions[1]=2;
    require(state.population.size()==3,"Sparse population count mismatch");
    const auto accepted=state;
    auto trial=state;
    trial.population[2].fraction*=2;clearPopulationMaterial(trial.population[48]);
    require(state.population[2].fraction==.01,"Rejected trial modified its accepted state");
    trial=accepted;equal(trial,accepted);
    auto moved=std::move(trial);equal(moved,accepted);
    require(trial.population.empty(),"Moved-from population table is inconsistent");

    std::stringstream checkpoint(std::ios::in|std::ios::out|std::ios::binary);
    writeSurfaceState(checkpoint,state);auto restored=readSurfaceState(checkpoint);equal(restored,state);
    require(restored.population.size()==3,"Checkpoint allocated absent pairs");
    clearPopulationMaterial(restored.population[48]);
    require(restored.population[48].events==42&&restored.population[48].hazard==.3,
        "Exhaustion reset a waiting-time clock");
    std::stringstream legacy(std::ios::in|std::ios::out|std::ios::binary);
    LegacySurfaceState old;
    for(int pair=0;pair<49;pair++)old.population[pair]=state.population.get(pair);
    old.inheritedNormal=state.inheritedNormal;old.births=state.births;old.promotions=state.promotions;
    writeSurfaceValue(legacy,old);
    auto imported=readSurfaceState(legacy,true);equal(imported,state);
    require(imported.population.size()==3,"Legacy import materialized unused pairs");
    std::stringstream previous(std::ios::in|std::ios::out|std::ios::binary);
    writeSurfaceValue(previous,state.population.mask());writeSurfaceValue(previous,state.inheritedNormal);
    writeSurfaceValue(previous,state.births);writeSurfaceValue(previous,state.promotions);
    for(int pair=0;pair<49;pair++)if(state.population.mask()&(std::uint64_t(1)<<pair))
        writeSurfaceValue(previous,static_cast<const PopulationMoments&>(state.population[pair]));
    equal(readSurfaceState(previous,false,true),state);
    auto film=state;film.population[2].filmArea=.25;
    std::stringstream sheets(std::ios::in|std::ios::out|std::ios::binary);
    writeSurfaceState(sheets,film);equal(readSurfaceState(sheets),film);

    const auto bytes=checkpoint.str();
    for(std::size_t length:{std::size_t(0),std::size_t(7),bytes.size()-1}){
        std::istringstream broken(bytes.substr(0,length),std::ios::binary);
        rejects([&](){readSurfaceState(broken);});
    }
    std::stringstream mask(std::ios::in|std::ios::out|std::ios::binary);
    writeSurfaceValue(mask,std::uint64_t(1)<<49);rejects([&](){readSurfaceState(mask);});
    auto invalid=state;invalid.population[2].hazard=std::numeric_limits<double>::quiet_NaN();
    std::stringstream nonfinite(std::ios::in|std::ios::out|std::ios::binary);
    writeSurfaceState(nonfinite,invalid);rejects([&](){readSurfaceState(nonfinite);});
    std::cout<<"SURFACE_STATE PASS empty_bytes="<<sizeof(SurfaceState)
             <<" legacy_bytes="<<sizeof(LegacySurfaceState)
             <<" populated_checkpoint_bytes="<<bytes.size()<<'\n';
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
