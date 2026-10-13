// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <array>
#include <cmath>
#include <string>
#include <stdexcept>
#include <limits>
#include "BrineDual.hpp"
namespace saltchem {
constexpr int nComponents=7,nSpecies=14,nMinerals=7;
enum Component {Na,K,Ca,Mg,Cl,S,C};
enum Species {NaIon,KIon,CaIon,MgIon,ClIon,Sulfate,Carbonate,H,OH,Bicarbonate,CO2,Bisulfate,MgOH,MgCO3};
enum Mineral {Halite,Sylvite,Carnallite,Polyhalite,Gypsum,Bischofite,Calcite};
using Totals=std::array<double,nComponents>;
using SpeciesVector=std::array<double,nSpecies>;
using MineralVector=std::array<double,nMinerals>;
constexpr double waterMolarMass=.018016,ln10=2.302585092994045684;
inline constexpr std::array<const char*,nSpecies> speciesNames{{"Na+","K+","Ca+2","Mg+2","Cl-","SO4-2","CO3-2","H+","OH-","HCO3-","CO2","HSO4-","MgOH+","MgCO3"}};
inline constexpr std::array<const char*,nComponents> componentNames{{"Na","K","Ca","Mg","Cl","S(6)","C(4)"}};
inline constexpr std::array<const char*,nMinerals> mineralNames{{"Halite","Sylvite","Carnallite","Polyhalite","Gypsum","Bischofite","Calcite"}};
inline constexpr std::array<int,nSpecies> charge{{1,1,2,2,-1,-2,-2,1,-1,-1,0,-1,1,0}};
inline constexpr SpeciesVector speciesMass{{.0229898,.0390983,.04008,.024305,.035453,.096064,.0600111,.001008,.017008,.0610191,.0440111,.097072,.041313,.0843161}};
// Products for one mole of mineral dissolution. S is sulfate-S, C is inorganic C.
inline constexpr std::array<Totals,nMinerals> stoichiometry{{
    {{1,0,0,0,1,0,0}},{{0,1,0,0,1,0,0}},{{0,1,0,1,3,0,0}},
    {{0,2,2,1,0,4,0}},{{0,0,1,0,0,1,0}},{{0,0,0,1,2,0,0}},{{0,0,1,0,0,0,1}}
}};
inline constexpr MineralVector hydration{{0,0,6,2,2,6,0}};
// Hydrogen-balanced component mass basis. Together with waterEquivalent these
// are an exact linear representation of aqueous mass, including acid/base and
// carbonate interconversion. Each value is kg per mol of conserved component.
inline constexpr Totals componentMass{{.0219818,.0380903,.038064,.022289,.036461,.09808,.0620271}};
inline constexpr MineralVector mineralMass{{.0584428,.0745513,.2778583,.6029496,.172176,.203307,.1000911}};
inline int speciesIndex(const std::string& s){for(int i=0;i<nSpecies;i++)if(s==speciesNames[i])return i;return -1;}
inline int mineralIndex(const std::string& s){for(int i=0;i<nMinerals;i++)if(s==mineralNames[i])return i;return -1;}
inline Totals componentTotals(const SpeciesVector& m){return {{m[NaIon],m[KIon],m[CaIon],m[MgIon]+m[MgOH]+m[MgCO3],m[ClIon],m[Sulfate]+m[Bisulfate],m[Carbonate]+m[Bicarbonate]+m[CO2]+m[MgCO3]}};}
inline double chargeBalance(const SpeciesVector& m){double v=0;for(int i=0;i<nSpecies;i++)v+=charge[i]*m[i];return v;}
inline double checkedExp(double x){
    // Reject inadmissible Newton/ODE trial states before OpenFOAM's enabled
    // floating-point traps can abort the process. Callers can then backtrack.
    if(std::isnan(x)||x>std::log(std::numeric_limits<double>::max())-1)throw std::domain_error("Exponential exceeds floating-point range");
    return std::exp(x);
}
inline double pow10(double x){return checkedExp(ln10*x);}
template<int N> Dual<N> checkedExp(const Dual<N>& x){const double v=checkedExp(x.value);return chain(x,v,v);}
template<int N> Dual<N> pow10(const Dual<N>& x){return checkedExp(ln10*x);}
struct AqueousState {
    SpeciesVector molality{},logGamma{},logActivity{};
    MineralVector saturationIndex{};
    double waterActivity=1,density=997,ionicStrength=0,pH=7;
    double chargeResidual=0,massResidual=0,equilibriumResidual=0;
    int iterations=0;
};
// Conservative inventory: mol of each component and kg of water-equivalent.
// W_eq = W_free + M_H2O * (n_OH + n_MgOH - n_CO2).
// This conserves O (and H through charge balance) during carbonate speciation.
struct Inventory {Totals aqueous{};MineralVector mineral{};double waterEquivalent=1;};
inline double aqueousMass(const Inventory& s){double mass=s.waterEquivalent;for(int j=0;j<nComponents;j++)mass+=componentMass[j]*s.aqueous[j];return mass;}
inline double strongCharge(const Totals& n){return n[Na]+n[K]+2*n[Ca]+2*n[Mg]-n[Cl]-2*n[S]-2*n[C];}
}
