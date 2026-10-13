#!/usr/bin/env bash
set +eu
source "$(dirname "$0")/environment.sh"
set -euo pipefail
cd "$SALT_PROJECT"
reference=0;checks=0;flow_only=0;tests_only=0;native_cpu=0
for arg in "$@"; do
 case "$arg" in
  --reference) reference=1;;
  --checks) checks=1;;
  --tests-only) tests_only=1;checks=1;;
  --flow-only) flow_only=1;;
  --native) native_cpu=1;;
  --help) echo 'Usage: bash scripts/build.sh [--reference] [--checks | --tests-only | --flow-only] [--native]'; exit 0;;
  *) echo "Unknown build option: $arg" >&2; exit 2;;
 esac
done
if (( flow_only && checks )); then echo '--flow-only cannot be combined with --checks or --tests-only' >&2; exit 2; fi
[[ -f "$SUNDIALS_PREFIX/include/cvode/cvode.h" ]] || { echo 'Set SUNDIALS_PREFIX to a SUNDIALS installation with CVODE.' >&2; exit 1; }
mkdir -p bin
stage=$(mktemp -d "$SALT_PROJECT/bin/.build.XXXXXX")
trap 'rm -rf -- "$stage"' EXIT
common=(-std=c++17 -O3 -g -fopenmp -Wall -Wextra -Wno-misleading-indentation
 -Isrc/chemistry -Isrc/geometry -Isrc/flow -Isrc/reference)
if (( native_cpu )); then common+=(-march=native -mtune=native); fi
chemistry=(-DSALT_HAS_CVODE -I"$SUNDIALS_PREFIX/include")
kinetic_libs=(-L"$SUNDIALS_PREFIX/lib" -Wl,--exclude-libs,ALL -lsundials_cvode
 -lsundials_nvecserial -lsundials_sunlinsoldense -lsundials_sunmatrixdense -lsundials_core)
foam=(-DWM_DP -DWM_LABEL_SIZE=32 -DNoRepository -I"$FOAM_SRC/OpenFOAM/lnInclude"
 -I"$FOAM_SRC/OSspecific/POSIX/lnInclude" -I"$FOAM_SRC/finiteVolume/lnInclude" -I"$FOAM_SRC/meshTools/lnInclude")
foam_libs=(-L"$FOAM_LIBBIN" -L"$FOAM_LIBBIN/$FOAM_MPI" -Wl,-rpath,"$FOAM_LIBBIN"
 -Wl,-rpath,"$FOAM_LIBBIN/$FOAM_MPI" -lfiniteVolume -lmeshTools -lOpenFOAM -lPstream)
if (( ! flow_only && ! tests_only )); then
 g++ "${common[@]}" "${chemistry[@]}" -fPIC -shared src/chemistry/SaltChemistryC.cpp "${kinetic_libs[@]}" -o "$stage/libSaltChemistry.so"
 g++ "${common[@]}" src/apps/saltChemistry.cpp -L"$stage" -Wl,-rpath,'$ORIGIN' -lSaltChemistry -o "$stage/saltChemistry"
fi
if (( ! tests_only )); then
 g++ "${common[@]}" "${chemistry[@]}" "${foam[@]}" src/flow/saltLakeFoam.C "${foam_libs[@]}" "${kinetic_libs[@]}" -o "$stage/saltLakeFoam"
fi
if (( ! flow_only && ! tests_only )); then
 g++ "${common[@]}" "${foam[@]}" -fPIC -shared src/flow/FoamSaltChemistry.C -L"$stage" -Wl,-rpath,'$ORIGIN' -lSaltChemistry "${foam_libs[@]}" -o "$stage/libFoamSaltChemistry.so"
fi
if (( checks )); then
 test_library="$stage"
 if (( tests_only )); then
  test_library="$SALT_PROJECT/bin"
  for library in libSaltChemistry.so libFoamSaltChemistry.so; do
   [[ -f "$test_library/$library" ]] || { echo "Build $library before --tests-only." >&2; exit 1; }
  done
 fi
 mkdir -p "$stage/tests"
 for name in verifyContracts verifyDerivatives verifyMineralInterface verifyProductTransport verifySurfaceReaction verifyNucleation verifyCapReactions verifySurfaceFilm verifySurfacePartition benchmarkSurfaceFilm; do
  g++ "${common[@]}" "${chemistry[@]}" "src/tests/$name.cpp" "${kinetic_libs[@]}" -o "$stage/tests/$name"
 done
 g++ "${common[@]}" src/tests/verifyChemistryApi.cpp -L"$test_library" -Wl,-rpath,'$ORIGIN/..' -lSaltChemistry -o "$stage/tests/verifyChemistryApi"
 g++ "${common[@]}" "${foam[@]}" src/tests/verifyFoamChemistry.C -L"$test_library" -Wl,-rpath,'$ORIGIN/..' -lFoamSaltChemistry -lSaltChemistry "${foam_libs[@]}" -o "$stage/tests/verifyFoamChemistry"
 g++ "${common[@]}" "${foam[@]}" src/tests/verifyTransportOperator.C "${foam_libs[@]}" -o "$stage/tests/verifyTransportOperator"
 g++ "${common[@]}" "${foam[@]}" src/tests/verifyFaceReactionExchange.C "${foam_libs[@]}" -o "$stage/tests/verifyFaceReactionExchange"
 g++ "${common[@]}" "${foam[@]}" src/tests/verifySurfaceLinkExchange.C "${foam_libs[@]}" -o "$stage/tests/verifySurfaceLinkExchange"
 g++ "${common[@]}" "${foam[@]}" src/tests/verifySurfaceRouting.C "${foam_libs[@]}" -o "$stage/tests/verifySurfaceRouting"
fi
if (( reference )); then
 [[ -f "$RM_PREFIX/include/PhreeqcRM.h" ]] || { echo 'Set RM_PREFIX to a PhreeqcRM installation.' >&2; exit 1; }
 mkdir -p "$stage/reference"
 rm_libs=(-L"$RM_PREFIX/lib" -Wl,-rpath,"$RM_PREFIX/lib" -lPhreeqcRM)
 if (( ! tests_only )); then
  g++ "${common[@]}" "${chemistry[@]}" "${foam[@]}" -I"$RM_PREFIX/include" -DSALT_REFERENCE_BACKEND src/flow/saltLakeFoam.C "${foam_libs[@]}" "${kinetic_libs[@]}" "${rm_libs[@]}" -o "$stage/reference/saltLakePhreeqcFoam"
 fi
 if (( ! flow_only )); then
  g++ "${common[@]}" "${chemistry[@]}" -I"$RM_PREFIX/include" src/tests/verifyTemperaturePressure.cpp "${kinetic_libs[@]}" "${rm_libs[@]}" -o "$stage/reference/verifyTemperaturePressure"
 fi
fi
for file in "$stage"/*; do
 if [[ -d "$file" ]]; then mkdir -p "bin/$(basename "$file")"; mv "$file"/* "bin/$(basename "$file")/";
 else mv "$file" bin/; fi
done
if (( tests_only )); then echo 'Built verification executables.'
else echo 'Built saltLakeFoam, saltChemistry and reusable C/OpenFOAM libraries.'; fi
