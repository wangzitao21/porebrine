// SPDX-License-Identifier: GPL-3.0-or-later
#include "SurfaceReaction.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

int main(int argc,char** argv){try{
    using namespace saltchem;
    if(argc<2||argc>3)throw std::invalid_argument("Supply pitzer.dat and optional iteration count");
    const int repetitions=argc==3?std::stoi(argv[2]):200;if(repetitions<1)throw std::invalid_argument("Positive iteration count required");
    PitzerModel model(argv[1]);Speciation spec(model);Reactions reactions(model);
    constexpr double volume=1e-6,area=1e-4,interval=.01;
    Inventory solid;solid.mineral[Halite]=.01;solid.waterEquivalent=.001;
    const auto dilute=spec.inventory({{.5,0,0,0,.5,0,0}},.001,{}),supersaturated=spec.inventory({{7,0,0,0,7,0,0}},.001,{});
    RateParameters rates;rates.geometry.enabled=true;rates.geometry.bulkVolume=volume;rates.neutral[Halite]=.02;
    MineralVector exposed{},locked{};exposed[Halite]=area;
    IntegrationOptions options;options.relativeTolerance=2e-10;options.absoluteTolerance=1e-14;
    struct Scenario{std::string name;SurfaceFilmParameters film;bool useFilm=true,growth=false;};
    std::vector<Scenario> scenarios(9);
    scenarios[0].name="bare_no_film";scenarios[0].useFilm=false;
    scenarios[1].name="bare_default_film";
    scenarios[2].name="zero_coverage_unused_film";scenarios[2].film.coatedVelocity[0]=1e-5;
    scenarios[3].name="fully_sealed";scenarios[3].film.coverage[0]=1;scenarios[3].film.coatedVelocity[0]=0;
    scenarios[4].name="partial_sealed";scenarios[4].film.coverage[0]=.7;scenarios[4].film.coatedVelocity[0]=0;
    scenarios[5].name="fully_covered_conducting";scenarios[5].film.coverage[0]=1;scenarios[5].film.coatedVelocity[0]=1e-5;
    scenarios[6].name="two_conducting_branches";scenarios[6].film.coverage[0]=.7;scenarios[6].film.bareVelocity[0]=5e-5;scenarios[6].film.coatedVelocity[0]=1e-5;
    scenarios[7].name="growth_no_film";scenarios[7].growth=true;scenarios[7].useFilm=false;
    scenarios[8].name="growth_conducting_film";scenarios[8].growth=true;scenarios[8].film.coverage[0]=1;scenarios[8].film.coatedVelocity[0]=1e-5;
    std::cout<<std::setprecision(12)<<"scenario,repetitions,microseconds_per_cell,film_speciations_per_cell,native_rate_evaluations_per_cell,mean_extent_mol\n";
    double checksum=0;
    for(const auto& scenario:scenarios){
        double filmCalls=0,rateCalls=0,extent=0;
        auto once=[&](){
            if(scenario.growth){
                auto cell=supersaturated;auto parameters=rates;parameters.geometry.separateMineralAreas=true;parameters.geometry.growthArea[0]=area;
                SurfaceFilmRateResult film;if(scenario.useFilm)film=applySurfaceFilmRates(reactions,cell,parameters,scenario.film);
                const auto stats=reactions.advance(cell,interval,parameters,options);filmCalls+=film.filmEvaluations;rateCalls+=stats.rateEvaluations;extent+=cell.mineral[0];
            }else{
                auto primary=solid,liquid=dilute;
                const auto result=advanceSurfaceDissolution(reactions,primary,liquid,rates,exposed,locked,model.molarVolume,volume,interval,options,scenario.useFilm?&scenario.film:nullptr);
                filmCalls+=result.filmEvaluations;rateCalls+=result.stats.rateEvaluations;extent+=result.dissolved[0];
            }
        };
        for(int i=0;i<5;i++)once();filmCalls=rateCalls=extent=0;
        const auto start=std::chrono::steady_clock::now();for(int i=0;i<repetitions;i++)once();
        const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();checksum+=extent;
        std::cout<<scenario.name<<','<<repetitions<<','<<elapsed*1e6/repetitions<<','<<filmCalls/repetitions<<','<<rateCalls/repetitions<<','<<extent/repetitions<<'\n';
    }
    if(!std::isfinite(checksum))throw std::runtime_error("Nonfinite benchmark result");
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
