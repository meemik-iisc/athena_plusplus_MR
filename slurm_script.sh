#!/usr/bin/env bash
# Submit an Athena++ run to Slurm (compute nodes cn01-cn08, 24 cores each).
#
# Put this script next to a copy of the athena executable (./athena), edit the
# settings below, then:
#   ./slurm_script.sh                      # fresh run
#   ./slurm_script.sh restart              # restart from RST_FILENAME
#   ./slurm_script.sh [restart] --dry-run  # print the job script, don't submit
set -euo pipefail

# ---------------------------- settings ---------------------------------------
JOB_NAME="wind_rad_bondi"
INPUT_FILENAME="athinput.wind_rad_bondi"
OUTPUT_FILENAME="wind_rad_bondi_10kyr"
RST_FILENAME="RunawayBH.00001.rst"     # <problem_id>.NNNNN.rst, for restarts
NPROC=8                                # MPI ranks; must be <= number of MeshBlocks

# Extra command-line parameter overrides, e.g. ("time/tlim=30" "problem/mu=0.6")
OVERRIDES=()

# Slurm: nodes are allocated whole (24 cores, 128 GB each), max time 15 days
PARTITION="normal"
NODES=1                                # need NODES*24 >= NPROC
TIME_LIMIT="2-00:00:00"                # D-HH:MM:SS
# -----------------------------------------------------------------------------

MODE=fresh
DRY_RUN=0
for arg in "$@"; do
  case "${arg}" in
    fresh|restart) MODE="${arg}" ;;
    --dry-run)     DRY_RUN=1 ;;
    *) echo "Usage: $0 [fresh|restart] [--dry-run]"; exit 1 ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

INPUT_PATH="${SCRIPT_DIR}/${INPUT_FILENAME}"
OUTPUT_PATH="${SCRIPT_DIR}/${OUTPUT_FILENAME}"
RESTART_PATH="${OUTPUT_PATH}/${RST_FILENAME}"   # Athena++ writes .rst into -d dir
EXE_PATH="${SCRIPT_DIR}/athena"

echo "SCRIPT_DIR        = ${SCRIPT_DIR}"
echo "EXECUTABLE_PATH   = ${EXE_PATH}"
echo "INPUT_PATH        = ${INPUT_PATH}"
echo "OUTPUT_PATH       = ${OUTPUT_PATH}"
echo "MODE              = ${MODE}"
echo "NPROC             = ${NPROC} on ${NODES} node(s), partition ${PARTITION}, time ${TIME_LIMIT}"
echo "OVERRIDES         = ${OVERRIDES[*]:-(none)}"

# Sanity checks
if [ ! -x "${EXE_PATH}" ]; then
    echo "Error: athena binary not found (or not executable) at ${EXE_PATH}"
    exit 1
fi
if ! command -v sbatch > /dev/null; then
    echo "Error: sbatch not found in PATH"
    exit 1
fi

case "${MODE}" in
  fresh)
    if [ ! -f "${INPUT_PATH}" ]; then
        echo "Error: input file not found at ${INPUT_PATH}"
        exit 1
    fi
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
esac

# Shell-quoted athena command line for the job script
ATHENA_CMD=$(printf '%q ' "${EXE_PATH}" "${START_ARGS[@]}" -d "${OUTPUT_PATH}" "${OVERRIDES[@]}")

JOB_SCRIPT=$(cat <<EOF
#!/bin/bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --partition=${PARTITION}
#SBATCH --nodes=${NODES}
#SBATCH --ntasks=${NPROC}
#SBATCH --time=${TIME_LIMIT}
#SBATCH --output=${LOG_BASE}.%j.out
#SBATCH --error=${LOG_BASE}.%j.err

# Compute nodes only see /opt/ohpc/pub/modulefiles; athena has rpaths to
# libmpi and libhdf5, so only mpirun is needed from the module.
command -v module > /dev/null || source /etc/profile.d/lmod.sh
module purge
module load openmpi/4.1.1

cd $(printf '%q' "${SCRIPT_DIR}")
echo "Job \${SLURM_JOB_ID} on \${SLURM_JOB_NODELIST}, started \$(date)"
mpirun -np ${NPROC} ${ATHENA_CMD}
echo "Finished \$(date)"
EOF
)

if [ "${DRY_RUN}" -eq 1 ]; then
    echo
    echo "${JOB_SCRIPT}"
    exit 0
fi

mkdir -p "${OUTPUT_PATH}"
if [ "${MODE}" = fresh ]; then
    # Keep a copy of the input that produced this run
    cp "${INPUT_PATH}" "${OUTPUT_PATH}/"
fi

JOB_ID=$(echo "${JOB_SCRIPT}" | sbatch --parsable)

echo "Submitted ${MODE} run as Slurm job ${JOB_ID}"
echo "Job status:   squeue -j ${JOB_ID}"
echo "Follow log:   tail -f ${LOG_BASE}.${JOB_ID}.out"
echo "Follow err:   tail -f ${LOG_BASE}.${JOB_ID}.err"
echo "Stop run:     scancel ${JOB_ID}"
