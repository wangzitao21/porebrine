#!/usr/bin/env bash
# OpenFOAM initialization uses unset variables and commands with nonzero status.
SALT_PROJECT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export SALT_PROJECT
foam_environment=${OPENFOAM_BASHRC:-${WM_PROJECT_DIR:-"$HOME/OpenFOAM/OpenFOAM-v2606"}/etc/bashrc}
if [[ ! -f "$foam_environment" ]]; then
    printf 'Set OPENFOAM_BASHRC to the OpenFOAM etc/bashrc file.\n' >&2
    return 1
fi
porebrine_shell_options=$-
set +eu
source "$foam_environment"
porebrine_environment_status=$?
[[ "$porebrine_shell_options" != *e* ]] || set -e
[[ "$porebrine_shell_options" != *u* ]] || set -u
if (( porebrine_environment_status )); then
    unset foam_environment porebrine_shell_options porebrine_environment_status
    return 1
fi
unset foam_environment porebrine_shell_options porebrine_environment_status
RM_PREFIX=${RM_PREFIX:-"$HOME/OpenFOAM/dependencies/phreeqcrm"}
SUNDIALS_PREFIX=${SUNDIALS_PREFIX:-"$HOME/OpenFOAM/dependencies/sundials"}
export RM_PREFIX SUNDIALS_PREFIX
