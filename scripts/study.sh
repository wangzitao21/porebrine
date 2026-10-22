#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/environment.sh" || exit 1
cd "$SALT_PROJECT"
exec python3 -B scripts/experiments/dual_resolution.py "$@"
