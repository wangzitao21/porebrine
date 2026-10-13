// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
#define SALT_CHEMISTRY_ABI_VERSION 1
#define SALT_COMPONENT_COUNT 7
#define SALT_MINERAL_COUNT 7
// Component order: Na, K, Ca, Mg, Cl, sulfate-S, inorganic-C.
// Mineral order: Halite, Sylvite, Carnallite, Polyhalite, Gypsum, Bischofite, Calcite.
typedef struct SaltChemistryHandle SaltChemistryHandle;
typedef struct {
    double aqueous_mol[7],mineral_mol[7];
    // kg water-equivalent; conserved O basis, NOT kg free liquid water.
    double water_equivalent_kg;
} SaltCell;
typedef struct {
    double k_neutral[7],k_acid[7],k_co2[7],acid_order[7];
    double area_m2[7],reference_mol[7],dissolution_exponent[7],affinity_exponent[7];
} SaltRates;
typedef struct {double relative_tolerance,absolute_tolerance_mol_per_kg;int maximum_attempts,threads;} SaltOptions;
typedef struct {
    double density_kg_m3,free_water_kg,aqueous_mass_kg,aqueous_volume_m3,pH,water_activity;
    double saturation_index[7],component_source_mol_s[7],mineral_change_mol[7];
    double aqueous_mass_source_kg_s;
    int accepted_steps,rejected_steps,rate_evaluations;
} SaltResult;
// Status-returning functions return 0 on success. No exceptions cross the C boundary.
int salt_chemistry_create(const char* database,double temperature_K,double pressure_atm,SaltChemistryHandle** handle);
// Change uniform thermodynamic conditions between calls; retains the parsed
// database. Do not call concurrently with another operation on this handle.
int salt_chemistry_set_conditions(SaltChemistryHandle* handle,double temperature_K,double pressure_atm);
// strict=1 rejects water-property temperature capping / pressure elevation.
// The compatible default (0) follows PHREEQC. Diagnostics use converged aw.
typedef struct {double requested_temperature_K,water_temperature_K,requested_pressure_atm,effective_pressure_atm;int temperature_capped,pressure_raised;} SaltWaterState;
int salt_chemistry_set_water_policy(SaltChemistryHandle* handle,int strict);
int salt_chemistry_water_state(const SaltChemistryHandle* handle,double water_activity,SaltWaterState* result);
void salt_chemistry_destroy(SaltChemistryHandle* handle);
const char* salt_chemistry_last_error(void);
// Initializes from mol/kg free water. Mineral inventory is in mol per cell.
int salt_chemistry_initialize(const SaltChemistryHandle* handle,const double totals_mol_kg[7],double free_water_kg,const double mineral_mol[7],SaltCell* cell,SaltResult* result);
// Whole-batch transaction: cells and results are committed only if every cell succeeds.
// Each rates[i] specifies the surface areas and kinetic parameters of cell i.
// dt=0 computes properties without advancing chemistry. Density is the value to
// return to the host flow solver; mass_source is the conservative aqueous source.
// count=0 is a valid no-op with null arrays. Wet cells require positive water;
// dry solid cells are excluded by the OpenFOAM adapter.
int salt_chemistry_advance(const SaltChemistryHandle* handle,size_t count,SaltCell* cells,const SaltRates* rates,double dt,const SaltOptions* options,SaltResult* results,int* cell_status);
#ifdef __cplusplus
}
#endif
