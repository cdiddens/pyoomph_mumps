#!/usr/bin/env bash

# Editable install of pyoomph_mumps, modelled on pyoomph's own build_for_develop.sh - including the
# trap, because pip signals a failed compile only through its exit status and the last lines of its
# output are CMake/ninja noise by which point the first "error:" has long scrolled past.
set -euo pipefail

verdict() {
  local rc=$?
  if [ "$rc" -eq 0 ]; then
    echo "build_for_develop: BUILD OK"
  else
    echo "build_for_develop: BUILD FAILED (exit $rc) - search the log above for 'error:'" >&2
  fi
  exit "$rc"
}
trap verdict EXIT

cd "$(dirname "$0")"

PYTHON=${PYTHON:-python3}

# MUST match how pyoomph itself was built. pyoomph records this as the presence or absence of a
# pyoomph/NO_MPI marker file, so it can be read without importing pyoomph - which matters here,
# because importing pyoomph initialises MPI and this script may run before pyoomph_mumps exists.
if [ -z "${PYOOMPH_MUMPS_USE_MPI:-}" ]; then
  PYOOMPH_MUMPS_USE_MPI=$($PYTHON - <<'PY'
import importlib.util, pathlib
spec = importlib.util.find_spec("pyoomph")
if spec is None or not spec.submodule_search_locations:
    print("OFF")
else:
    root = pathlib.Path(list(spec.submodule_search_locations)[0])
    print("OFF" if (root / "NO_MPI").exists() else "ON")
PY
)
  echo "build_for_develop: pyoomph was built with MPI=$PYOOMPH_MUMPS_USE_MPI, matching it"
fi

# system   - link the MUMPS already installed on this machine. Much the faster loop, and what to use
#            while developing; needs libmumps-dev (Debian/Ubuntu) or equivalent, WITH the complex
#            arithmetic (libzmumps).
# download - the package default: fetch and build MUMPS (and, if there is no BLAS, a reference
#            LAPACK) from source. Slow the first time, then cached in build/.
MODE=${PYOOMPH_MUMPS_BUILD_MODE:-system}
case "$MODE" in
  system)   MUMPS_ARGS=(--config-settings=cmake.define.PYOOMPH_MUMPS_DOWNLOAD=OFF) ;;
  download) MUMPS_ARGS=(--config-settings=cmake.define.PYOOMPH_MUMPS_DOWNLOAD=ON) ;;
  *) echo "PYOOMPH_MUMPS_BUILD_MODE must be 'system' or 'download', got '$MODE'" >&2; exit 2 ;;
esac
echo "build_for_develop: MUMPS acquisition mode: $MODE"

# See pyoomph's script: Ubuntu's PEP 668 marker makes pip refuse the editable install outright, and
# the install lands in the user site-packages either way.
PEP668=()
if $PYTHON -c 'import os, sys, sysconfig; sys.exit(0 if sys.prefix == sys.base_prefix and os.path.exists(os.path.join(sysconfig.get_path("stdlib"), "EXTERNALLY-MANAGED")) else 1)'; then
  PEP668=(--break-system-packages)
fi

$PYTHON -m pip install ${PEP668[@]+"${PEP668[@]}"} --no-build-isolation -e . -v \
    --config-settings=editable.mode=redirect \
    --config-settings=build-dir=build \
    --config-settings=build.verbose=true \
    --config-settings=build.tool-args=-j4 \
    --config-settings=cmake.build-type=RelWithDebInfo \
    --config-settings=cmake.define.PYOOMPH_MUMPS_USE_MPI="$PYOOMPH_MUMPS_USE_MPI" \
    "${MUMPS_ARGS[@]}"
