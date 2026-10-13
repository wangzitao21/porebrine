// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "fvCFD.H"
#include "SaltChemistryC.h"
namespace Foam {
// Chemical fields use mol/m3 bulk for components, minerals and reference moles,
// kg/m3 bulk for waterEquivalent, and m2/m3 bulk for reactiveArea.
// rho is kg/m3 liquid; aqueousMassSource is kg/(m3 bulk s).
// No transport or pressure equation is advanced by this adapter.
// Chemical increments are already committed to the inventory fields. Do not
// apply componentSources a second time to those updated fields. rho and the
// aqueous mass source are supplied for the host's consistent flow coupling.
void advanceSaltChemistry(
    const fvMesh& mesh,const SaltChemistryHandle* handle,
    PtrList<volScalarField>& aqueousMoles,PtrList<volScalarField>& mineralMoles,
    volScalarField& waterEquivalent,const PtrList<volScalarField>& reactiveArea,
    const PtrList<volScalarField>& referenceMoles,const SaltRates& coefficients,
    double dt,const SaltOptions& options,volScalarField& rho,
    volScalarField& aqueousMassSource,PtrList<volScalarField>& componentSources);
}
