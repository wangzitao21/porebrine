#!/usr/bin/env bash
# Source before enabling nounset: OpenFOAM's environment uses unset variables.
SALT_PROJECT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export SALT_PROJECT
source "${OPENFOAM_BASHRC:-/home/wang/OpenFOAM/OpenFOAM-v2606/etc/bashrc}"
RM_PREFIX=${RM_PREFIX:-/home/wang/OpenFOAM/dependencies/phreeqcrm}
SUNDIALS_PREFIX=${SUNDIALS_PREFIX:-/home/wang/OpenFOAM/dependencies/sundials}
export RM_PREFIX SUNDIALS_PREFIX
