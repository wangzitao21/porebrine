// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "BrineReactions.hpp"
#include <memory>
#include <vector>
#include <omp.h>

namespace saltflow {
// Each worker owns the temperature/pressure-dependent coefficients. Neither
// the database nor a mutable Pitzer model is shared between OpenMP workers.
struct ChemicalWorker {
    saltchem::PitzerModel model;
    saltchem::Speciation speciation;
    saltchem::Reactions reactions;
    explicit ChemicalWorker(const saltchem::PitzerModel& master):
        model(master),speciation(model),reactions(model) {}
};
class LocalChemistry {
    std::vector<std::unique_ptr<ChemicalWorker>> workers;
public:
    LocalChemistry(const saltchem::PitzerModel& master,int threads) {
        if(threads<1)throw std::invalid_argument("chemistryThreads must be positive");
        for(int i=0;i<threads;i++)workers.emplace_back(new ChemicalWorker(master));
    }
    ChemicalWorker& at(double temperatureK,double pressureAtm) {
        auto& worker=*workers.at(omp_get_thread_num());
        worker.model.setConditions(temperatureK,pressureAtm);
        return worker;
    }
};
}
