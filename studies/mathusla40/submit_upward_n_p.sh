basedir="/project/6049244/data/MATHUSLA/simulation_v2"

## Special settings for the Dry run
# 1. -c n : Do no delete simulation and digi files,
# 2. -o n : Do no delete simulation and digi files,

SUBMIT="False" # {True, False, Run}

events_per_job=400000
N_JOB=20
JOB_TIME_HOURS=3
RUN_NAMES=(proton neutron)
RUN_ENERGIES=(2000 5000 10000)
SIM_REPO_DIR=$(realpath ../..)

for run_name in "${RUN_NAMES[@]}"; do

    for run_energy in "${RUN_ENERGIES[@]}"; do

        data_directory=$basedir/upward/${run_name}_$run_energy
        mkdir -p $data_directory

        slurm_text="#!/bin/bash
#SBATCH --time=${JOB_TIME_HOURS}:00:00
#SBATCH --account=rrg-mdiamond
#SBATCH --array=1-${N_JOB}
#SBATCH --mem=4G
#SBATCH --job-name=musim-${run_name}
#SBATCH --constraint=[skylake|cascade]
#SBATCH --output=$basedir/upward/log/${run_name}_%a.out
echo current job id is: \$SLURM_ARRAY_TASK_ID
./upward_n_and_p/_start_single_run.sh $SIM_REPO_DIR $data_directory $events_per_job \$SLURM_ARRAY_TASK_ID $run_name $run_energy
"

        # Run or print the slurm commands
        if [[ "$SUBMIT" == "True" ]]; then
            echo "$slurm_text" >slurm_script_temp.txt
            sbatch slurm_script_temp.txt
            rm -f slurm_script_temp.txt
        elif [[ "$SUBMIT" == "Run" ]]; then
            # Run a single series with seriese number of 0
            echo "Running command: " "$cmd_start_series"
            echo "$cmd_start_series" -s 0 >slurm_script_temp.txt
            bash slurm_script_temp.txt
            rm -f slurm_script_temp.txt
        else
            echo Debug: slurm script to be sumbitted
            echo "-------------------------------------------------------"
            echo "$slurm_text"
            echo "-------------------------------------------------------"
        fi
    done
done
