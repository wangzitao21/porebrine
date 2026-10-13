// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "SurfaceFilm.hpp"

namespace saltchem {
struct SurfaceReactionResult {
    IntegrationStats stats;
    MineralVector dissolved{},retreat{};
    int filmEvaluations=0;
    MineralVector filmTransmission{{1,1,1,1,1,1,1}};
};

// Inventories on both sides use the same representative bulk-volume basis.
// Face area and cellVolume describe the physical mesh cell. Thus a density
// inventory uses geometry.bulkVolume=1, whereas an extensive inventory uses
// geometry.bulkVolume=cellVolume. Each caller must serialize faces sharing
// either endpoint; this function consumes the current liquid inventory once.
inline SurfaceReactionResult advanceSurfaceDissolution(
    const Reactions& reactions,Inventory& solidOwner,Inventory& liquidOwner,
    const RateParameters& solidRates,const MineralVector& faceArea,
    const MineralVector& lockedCaps,const MineralVector& molarVolume,
    double cellVolume,double interval,const IntegrationOptions& options={},const SurfaceFilmParameters* film=nullptr){
    if(&solidOwner==&liquidOwner)throw std::invalid_argument("Surface reaction requires distinct solid and liquid owners");
    if(!std::isfinite(cellVolume)||cellVolume<=0||!solidRates.geometry.enabled)
        throw std::invalid_argument("Surface reaction requires a positive cell volume and solid-owner geometry");
    if(!std::isfinite(interval)||interval<0)throw std::invalid_argument("Invalid surface-reaction interval");
    RateParameters parameters=solidRates;
    parameters.geometry.separateMineralAreas=true;
    const double bulkVolume=parameters.geometry.bulkVolume;
    if(!std::isfinite(bulkVolume)||bulkVolume<=0)throw std::invalid_argument("Invalid surface-reaction inventory volume basis");
    Inventory initial=liquidOwner;initial.mineral=solidOwner.mineral;
    for(int j=0;j<nMinerals;j++){
        if(!std::isfinite(faceArea[j])||faceArea[j]<0||!std::isfinite(lockedCaps[j])||lockedCaps[j]<0||
           lockedCaps[j]>initial.mineral[j]||!std::isfinite(molarVolume[j])||molarVolume[j]<=0)
            throw std::invalid_argument("Invalid exposed surface, locked cap inventory or mineral volume");
        parameters.area[j]=faceArea[j]*bulkVolume/cellVolume;
        parameters.geometry.growthArea[j]=0;
        parameters.geometry.surfaceExponent[j]=0;
        parameters.geometry.curvatureLogOmega[j]=0;
        parameters.geometry.lockedMoles[j]=lockedCaps[j];
        parameters.referenceMoles[j]=std::max(initial.mineral[j]-lockedCaps[j],1e-30);
    }
    // The liquid provides composition and solvent, while the solid owner
    // supplies the pore-volume constraint. A previously cached porosity from
    // another face/chemical stage cannot replace the current solid inventory.
    parameters.geometry.initialPorosity=poreFraction(solidOwner,parameters.geometry,molarVolume);
    Inventory reacted=initial;SurfaceReactionResult result;
    if(film){
        // Stock availability belongs to this conservative two-owner step;
        // the reusable film helper also supports stock-free growth/probes.
        for(int j=0;j<nMinerals;j++)if(interval==0||initial.mineral[j]==lockedCaps[j])parameters.area[j]=0;
        const auto applied=applySurfaceFilmRates(reactions,initial,parameters,*film);
        result.filmEvaluations=applied.filmEvaluations;result.filmTransmission=applied.filmTransmission;
    }
    // The film multiplier is frozen over this interval. The seven-mineral
    // bulk ODE remains coupled, but this is first-order film/bulk splitting,
    // not a fully implicit evolving surface-composition solve.
    result.stats=reactions.advance(reacted,interval,parameters,options);
    for(int j=0;j<nMinerals;j++){
        const double extent=initial.mineral[j]-reacted.mineral[j];
        const double available=initial.mineral[j]-lockedCaps[j];
        const double roundoff=64*std::numeric_limits<double>::epsilon()*std::max(initial.mineral[j],reacted.mineral[j]);
        if(extent< -roundoff||extent>available+roundoff)
            throw std::runtime_error("Dissolution-only surface step exceeded accessible mineral bounds");
        result.dissolved[j]=std::clamp(extent,0.,available);
        result.retreat[j]=faceArea[j]>0?result.dissolved[j]*molarVolume[j]*cellVolume/(bulkVolume*faceArea[j]):0;
    }
    // Reconstruct both sides from one set of reaction extents. This retains
    // exact stoichiometric pairing, including hydrate water, rather than
    // separately copying two independently rounded inventory differences.
    const auto transferred=Reactions::changed(initial,result.dissolved);
    result.stats.dissolved=result.dissolved;result.stats.precipitated.fill(0);
    solidOwner.mineral=transferred.mineral;
    liquidOwner.aqueous=transferred.aqueous;
    liquidOwner.waterEquivalent=transferred.waterEquivalent;
    return result;
}
}
