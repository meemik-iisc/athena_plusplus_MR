#!/usr/bin/env bash
# Launch an Athena++ run in the background with MPI.
#
# Put this script next to a copy of the athena executable (./athena), edit the
# settings below, then:
#   ./run_file.sh            # fresh run
#   ./run_file.sh restart    # restart from RST_FILENAME
set -euo pipefail

# ---------------------------- settings ---------------------------------------
JOB_NAME="wind_rad_bondi"
BUILD_NAME="build_wind_rad_bondi"
INPUT_FILENAME="athinput.wind_rad_bondi"
OUTPUT_FILENAME="wind_rad_bondi_10kyr"
RST_FILENAME="RunawayBH.00001.rst"     # <problem_id>.NNNNN.rst, for restarts
NPROC=1                                # MPI ranks; must be <= number of MeshBlocks

# Extra command-line parameter overrides, e.g. ("time/tlim=30" "problem/mu=0.6")
OVERRIDES=()

ATHENA_ROOT="/home/meemik/Meemik/athena_plusplus_MR"
# -----------------------------------------------------------------------------

MODE="${1:-fresh}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


INPUT_PATH="${SCRIPT_DIR}/${INPUT_FILENAME}"
OUTPUT_PATH="${SCRIPT_DIR}/${OUTPUT_FILENAME}"
RESTART_PATH="${OUTPUT_PATH}/${RST_FILENAME}"   # Athena++ writes .rst into -d dir
EXE_PATH="${SCRIPT_DIR}/athena"

echo "ATHENA_ROOT       = ${ATHENA_ROOT}"
echo "SCRIPT_DIR        = ${SCRIPT_DIR}"
echo "EXECUTABLE_PATH   = ${EXE_PATH}"
echo "INPUT_PATH        = ${INPUT_PATH}"
echo "OUTPUT_PATH       = ${OUTPUT_PATH}"
echo "MODE              = ${MODE}"
echo "NPROC             = ${NPROC}"
echo "OVERRIDES         = ${OVERRIDES[*]:-(none)}"

# Sanity checks
if [ ! -x "${EXE_PATH}" ]; then
    echo "Error: athena binary not found (or not executable) at ${EXE_PATH}"
    exit 1
fi
if ! command -v mpirun > /dev/null; then
    echo "Error: mpirun not found in PATH"
    exit 1
fi

mkdir -p "${OUTPUT_PATH}"
cd "${SCRIPT_DIR}"

case "${MODE}" in
  fresh)
    if [ ! -f "${INPUT_PATH}" ]; then
        echo "Error: input file not found at ${INPUT_PATH}"
        exit 1
    fi
    # Keep a copy of the input that produced this run
    cp "${INPUT_PATH}" "${OUTPUT_PATH}/"
    START_ARGS=(-i "${INPUT_PATH}")
    LOG_BASE="${SCRIPT_DIR}/${JOB_NAME}"
    ;;
  restart)
    if [ ! -f "${RESTART_PATH}" ]; then
        echo "Error: restart file not found at ${RESTART_PATH}"
        exit 1
    fi
    START_ARGS=(-r "${RESTART_PATH}")
    LOG_BASE="${SCRIPT_DIR}/${JOB_NAME}_restart"
    ;;
  *)
    echo "Usage: $0 [fresh|restart]"
    exit 1
    ;;
esac

nohup mpirun -np "${NPROC}" "${EXE_PATH}" \
    "${START_ARGS[@]}" \
    -d "${OUTPUT_PATH}" \
    "${OVERRIDES[@]}" \
    > "${LOG_BASE}.out" \
    2> "${LOG_BASE}.err" \
    < /dev/null &

echo "Started ${MODE} run with PID $!"
echo "Follow log:   tail -f ${LOG_BASE}.out"
echo "Follow err:   tail -f ${LOG_BASE}.err"
echo "Stop run:     kill $!"
