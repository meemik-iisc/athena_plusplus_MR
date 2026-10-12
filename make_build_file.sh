#!/usr/bin/env bash
# Configure Athena++ for a given problem generator with this machine's flags.
# Only runs configure.py; run `make clean && make -j` yourself afterwards.
#
# Usage:
#   ./make_build_file.sh wind_rad_bondi                    # cylindrical, MPI + parallel HDF5
#   ./make_build_file.sh blast --coord cartesian           # different coordinates
#   ./make_build_file.sh wind_rad_bondi --no-mpi           # serial build, serial HDF5
#   ./make_build_file.sh wind_rad_bondi --dry-run          # print the commands only
#   ./make_build_file.sh wind_rad_bondi -- -b --flux=hlld  # extra configure.py flags after --

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# System HDF5 builds: the openmpi one is needed for MPI runs (parallel HDF5 output)
HDF5_MPI_INC=/usr/include/hdf5/openmpi
HDF5_MPI_LIB=/usr/lib/x86_64-linux-gnu/hdf5/openmpi
HDF5_SERIAL_INC=/usr/include/hdf5/serial
HDF5_SERIAL_LIB=/usr/lib/x86_64-linux-gnu/hdf5/serial

# g++ >= 12 on x86 defines __FLT16_MAX__ but has no __fp16 type, which breaks
# src/athena.hpp when HDF5 output is enabled; undefining the macro disables fp16.
FP16_WORKAROUND=-U__FLT16_MAX__

usage() {
  sed -n '2,10p' "$0" | sed 's/^# \{0,1\}//'
  echo
  echo "Options: --coord <system> (default cylindrical), --no-mpi, --no-hdf5,"
  echo "         --dry-run, -h/--help"
  exit "${1:-0}"
}

[[ $# -lt 1 ]] && usage 1

prob=""
coord=cylindrical
mpi=1
hdf5=1
dry_run=0
extra=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage 0 ;;
    --coord)   coord="$2"; shift 2 ;;
    --no-mpi)  mpi=0; shift ;;
    --no-hdf5) hdf5=0; shift ;;
    --dry-run) dry_run=1; shift ;;
    --)        shift; extra=("$@"); break ;;
    -*)        echo "Unknown option: $1" >&2; usage 1 ;;
    *)         if [[ -n "$prob" ]]; then echo "Only one problem generator allowed" >&2; usage 1; fi
               prob="${1%.cpp}"; shift ;;
  esac
done

[[ -z "$prob" ]] && usage 1
if [[ ! -f "$ROOT/src/pgen/$prob.cpp" ]]; then
  echo "Problem generator not found: $ROOT/src/pgen/$prob.cpp" >&2
  exit 1
fi

configure=(python3 configure.py "--prob=$prob" "--coord=$coord")
if [[ $mpi -eq 1 ]]; then
  configure+=(-mpi)
fi
if [[ $hdf5 -eq 1 ]]; then
  if [[ $mpi -eq 1 ]]; then
    configure+=(-hdf5 "--include=$HDF5_MPI_INC" "--lib_path=$HDF5_MPI_LIB")
  else
    configure+=(-hdf5 "--include=$HDF5_SERIAL_INC" "--lib_path=$HDF5_SERIAL_LIB")
  fi
  configure+=("--cflag=$FP16_WORKAROUND")
fi
configure+=("${extra[@]}")

run() {
  echo ">> $*"
  if [[ $dry_run -eq 0 ]]; then
    (cd "$ROOT" && "$@")
  fi
}

run "${configure[@]}"

if [[ $dry_run -eq 0 ]]; then
  echo
  echo "Configured $prob. To build:"
  echo "  cd $ROOT && make clean && make -j"
fi
