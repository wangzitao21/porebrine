#!/usr/bin/env bash
set +eu
source "$(dirname "$0")/environment.sh" || exit 1
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
cd "$SALT_PROJECT"
exec python3 -B scripts/run.py "$@"
