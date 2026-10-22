#!/usr/bin/env bash
set -euo pipefail
if [[ "${1:-}" == --help && $# == 1 ]]; then
 printf 'Usage: bash scripts/build_geometry.sh\nBuild pinned GEOS into bin/geos. Set GEOS_BUILD_JOBS to control parallel compilation.\n'
 exit 0
fi
if (( $# )); then printf 'Unknown geometry build option: %s\n' "$1" >&2; exit 2; fi
geometry_root=$(cd "$(dirname "$0")/.." && pwd)
geometry_work=$(mktemp -d /tmp/porebrine-geos.XXXXXX)
trap 'rm -rf -- "$geometry_work"' EXIT
geometry_archive="$geometry_work/geos.tar.bz2"
if [[ -f "$geometry_root/bin/.geos-source/geos.tar.bz2" ]]; then
 cp "$geometry_root/bin/.geos-source/geos.tar.bz2" "$geometry_archive"
else
 curl -L -f --retry 2 https://download.osgeo.org/geos/geos-3.13.1.tar.bz2 -o "$geometry_archive"
fi
printf '%s  %s\n' df2c50503295f325e7c8d7b783aca8ba4773919cde984193850cf9e361dfd28c "$geometry_archive" | sha256sum -c -
tar -xjf "$geometry_archive" -C "$geometry_work"
# GEOS only uses this RNG to balance its hot-pixel spatial index. A fixed
# permutation seed is independent of chemistry/nucleation and reproducible.
# This also avoids repeated hardware-entropy setup in hot geometry operations.
python3 - "$geometry_work/geos-3.13.1/src/noding/snapround/HotPixelIndex.cpp" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1]);text=p.read_text()
old='std::random_device rd;\n    std::mt19937 g(rd());'
if text.count(old)!=2:
    raise RuntimeError('GEOS source version does not match the pinned index patch')
p.write_text(text.replace(old,'std::mt19937 g(0x5A17u);'))
PY
cmake -S "$geometry_work/geos-3.13.1" -B "$geometry_work/build" \
 -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF -DBUILD_GEOSOP=OFF \
 -DBUILD_BENCHMARKS=OFF -DBUILD_DOCUMENTATION=OFF -DBUILD_SHARED_LIBS=ON \
 -DCMAKE_INSTALL_PREFIX="$geometry_root/bin/geos" -DCMAKE_INSTALL_LIBDIR=lib \
 -DCMAKE_INSTALL_RPATH='$ORIGIN'
cmake --build "$geometry_work/build" --parallel "${GEOS_BUILD_JOBS:-8}"
cmake --install "$geometry_work/build"
printf 'Built local GEOS 3.13.1 with a deterministic spatial-index shuffle.\n'
