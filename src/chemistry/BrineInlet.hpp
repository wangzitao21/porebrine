// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "BrineSpeciation.hpp"
namespace saltchem {
inline Inventory haliteSaturatedInlet(const Speciation& spec){
    auto state=[&](double m){return spec.solve(Totals{{m,0,0,0,m,0,0}});};
    double lo=.01,hi=15;
    if(state(lo).saturationIndex[Halite]>=0||state(hi).saturationIndex[Halite]<=0)
        throw std::runtime_error("Could not bracket the halite-saturated inlet between 0.01 and 15 mol/kg");
    for(int i=0;i<48;i++){const double mid=.5*(lo+hi);if(state(mid).saturationIndex[Halite]>0)hi=mid;else lo=mid;}
    const double m=.5*(lo+hi);return spec.inventory(Totals{{m,0,0,0,m,0,0}},1);
}
}
