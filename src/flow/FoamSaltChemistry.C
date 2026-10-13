// SPDX-License-Identifier: GPL-3.0-or-later
#include "FoamSaltChemistry.H"
#include <vector>
namespace Foam {
void advanceSaltChemistry(const fvMesh& mesh,const SaltChemistryHandle* handle,
    PtrList<volScalarField>& aqueousMoles,PtrList<volScalarField>& mineralMoles,
    volScalarField& waterEquivalent,const PtrList<volScalarField>& reactiveArea,
    const PtrList<volScalarField>& referenceMoles,const SaltRates& coefficients,
    double dt,const SaltOptions& options,volScalarField& rho,
    volScalarField& aqueousMassSource,PtrList<volScalarField>& componentSources){
    if(aqueousMoles.size()!=7||mineralMoles.size()!=7||reactiveArea.size()!=7||referenceMoles.size()!=7||componentSources.size()!=7)FatalErrorInFunction<<"Expected seven component and mineral fields"<<exit(FatalError);
    const dimensionSet concentration=dimMoles/dimVolume;
    if(waterEquivalent.dimensions()!=dimDensity||rho.dimensions()!=dimDensity||aqueousMassSource.dimensions()!=dimDensity/dimTime)FatalErrorInFunction<<"Incorrect density/source field dimensions"<<exit(FatalError);
    for(int j=0;j<7;j++)if(aqueousMoles[j].dimensions()!=concentration||mineralMoles[j].dimensions()!=concentration||referenceMoles[j].dimensions()!=concentration||reactiveArea[j].dimensions()!=dimless/dimLength||componentSources[j].dimensions()!=concentration/dimTime)FatalErrorInFunction<<"Incorrect chemical field dimensions"<<exit(FatalError);
    std::vector<label> wet;label failed=0;std::string message="Failure on another MPI rank";
    forAll(waterEquivalent,i){
        if(waterEquivalent[i]==0){for(int j=0;j<7;j++)if(aqueousMoles[j][i]!=0){failed=1;message="Dry cell has nonzero aqueous inventory";}}
        else wet.push_back(i);
    }
    const size_t n=wet.size();std::vector<SaltCell> cells(n);std::vector<SaltRates> rates(n,coefficients);std::vector<SaltResult> result(n);std::vector<int> status(n);
    for(size_t k=0;k<n;k++){label i=wet[k];double volume=mesh.V()[i];cells[k].water_equivalent_kg=waterEquivalent[i]*volume;for(int j=0;j<7;j++){cells[k].aqueous_mol[j]=aqueousMoles[j][i]*volume;cells[k].mineral_mol[j]=mineralMoles[j][i]*volume;rates[k].area_m2[j]=reactiveArea[j][i]*volume;rates[k].reference_mol[j]=referenceMoles[j][i]*volume;}}
    if(!failed){failed=salt_chemistry_advance(handle,n,cells.data(),rates.data(),dt,&options,result.data(),status.data());if(failed)message=salt_chemistry_last_error();}
    reduce(failed,maxOp<label>());
    if(failed)FatalErrorInFunction<<message<<exit(FatalError);
    // Commit after all MPI ranks succeeded. The source and rho can then be used
    // by the host mass-continuity and momentum equations.
    forAll(waterEquivalent,i){aqueousMassSource[i]=0;for(int j=0;j<7;j++)componentSources[j][i]=0;}
    for(size_t k=0;k<n;k++){label i=wet[k];double volume=mesh.V()[i];waterEquivalent[i]=cells[k].water_equivalent_kg/volume;rho[i]=result[k].density_kg_m3;aqueousMassSource[i]=result[k].aqueous_mass_source_kg_s/volume;
        for(int j=0;j<7;j++){aqueousMoles[j][i]=cells[k].aqueous_mol[j]/volume;mineralMoles[j][i]=cells[k].mineral_mol[j]/volume;componentSources[j][i]=result[k].component_source_mol_s[j]/volume;}}
    waterEquivalent.correctBoundaryConditions();rho.correctBoundaryConditions();aqueousMassSource.correctBoundaryConditions();
    for(int j=0;j<7;j++){aqueousMoles[j].correctBoundaryConditions();mineralMoles[j].correctBoundaryConditions();componentSources[j].correctBoundaryConditions();}
}
}
