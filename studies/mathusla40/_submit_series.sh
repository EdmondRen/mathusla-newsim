## Handel input arguments

## Default parameters
SCRIPT_FILENAME=""

SIM_REPO_DIR=$(realpath ../..)
SIM_TEMP_DIR=""
MERGE_DIR=""
N_EVENTS_PER_RUN=1000
N_RUNS_PER_SERIES=3
N_JOB=1
JOB_TIME_HOURS=""
JOB_NAME="cosmic"

SUBMIT="False" # {True, False, Run}
CLEAN_FILES="y"

usage() {
    echo "Usage: $0 [-p SIM_REPO_DIR] [-o SIM_TEMP_DIR] [-m MERGE_DIR] [-e N_EVENTS_PER_RUN] [-r N_RUNS_PER_SERIES] [-j N_JOBs] [-t JOB_TIME_HOURs] [-c y or n]" 1>&2
    echo " [-c] cleans up simulation and only keep the merged reconstruction file"
    exit 1
}
while getopts f:p:o:m:e:r:j:t:s:c:n: flag; do
    case "${flag}" in
    f) SCRIPT_FILENAME=${OPTARG} ;;
    p) SIM_REPO_DIR=${OPTARG} ;;
    o) SIM_TEMP_DIR=${OPTARG} ;;
    m) MERGE_DIR=${OPTARG} ;;
    e) N_EVENTS_PER_RUN=${OPTARG} ;;
    r) N_RUNS_PER_SERIES=${OPTARG} ;;
    j) N_JOB=${OPTARG} ;;
    t) JOB_TIME_HOURS=${OPTARG} ;;
    s) SUBMIT=${OPTARG} ;;
    c) CLEAN_FILES=${OPTARG} ;;
    n) JOB_NAME=${OPTARG} ;;
    h) usage ;;
    esac
done

# Shift positional arguments after options
shift $((OPTIND - 1))

# Make a slurm script

if [[ ! -z $SIM_TEMP_DIR ]]; then
    cmd_start_series="${SCRIPT_FILENAME} \\
        -p ${SIM_REPO_DIR} \\
        -o ${SIM_TEMP_DIR} \\ 
        -m ${MERGE_DIR} \\
        -e ${N_EVENTS_PER_RUN} \\
        -r ${N_RUNS_PER_SERIES} \\
        -c ${CLEAN_FILES}"
else
    cmd_start_series="${SCRIPT_FILENAME} \\
        -p ${SIM_REPO_DIR} \\
        -m ${MERGE_DIR} \\
        -e ${N_EVENTS_PER_RUN} \\
        -r ${N_RUNS_PER_SERIES} \\
        -c ${CLEAN_FILES}"
fi

# Make directory for logs
mkdir -p ${MERGE_DIR}/log
# Make a separate log file to record finished jobs
echo "Finished Jobs" > ${MERGE_DIR}/finshed_list.log

slurm_text="#!/bin/bash
#SBATCH --time=${JOB_TIME_HOURS}:00:00
#SBATCH --account=rrg-mdiamond
#SBATCH --array=1-${N_JOB}
#SBATCH --mem=4G
#SBATCH --job-name=musim-${JOB_NAME}
#SBATCH --constraint=[skylake|cascade]
#SBATCH --output=${MERGE_DIR}/log/${JOB_NAME}_%a.out

echo current job id is: \$SLURM_ARRAY_TASK_ID

bash ${cmd_start_series} -s \$SLURM_ARRAY_TASK_ID"

# Run or print the slurm commands
if [[ "$SUBMIT" == "True" ]]; then
    echo "$slurm_text" >slurm_script_temp.txt
    sbatch slurm_script_temp.txt
    rm -f slurm_script_temp.txt
elif [[ "$SUBMIT" == "Run" ]]; then
    # Run a single series with seriese number of 0
    echo "Running command: " "$cmd_start_series"
    echo "$cmd_start_series" -s 0 > slurm_script_temp.txt
    bash slurm_script_temp.txt
    rm -f slurm_script_temp.txt
else
    echo Debug: slurm script to be sumbitted
    echo "-------------------------------------------------------"
    echo "$slurm_text"
    echo "-------------------------------------------------------"
fi
