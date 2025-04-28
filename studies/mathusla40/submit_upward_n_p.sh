basedir_data="/project/6049244/data/MATHUSLA/simulation_v2"

## Special settings for the Dry run
# 1. -c n : Do no delete simulation and digi files,
# 2. -o n : Do no delete simulation and digi files,

SUBMIT="Run" # {True, False, Run}

events_per_job=1000000
N_JOB=1
JOB_TIME_HOURS=3
RUN_NAMES=(proton neutron)
RUN_ENERGIES=(1000)
# RUN_ENERGIES=(2000 5000 10000)
SIM_REPO_DIR=$(realpath ../..)

for run_name in "${RUN_NAMES[@]}"; do

    for run_energy in "${RUN_ENERGIES[@]}"; do

        data_directory=$basedir_data/upward/${run_name}_$run_energy
        mkdir -p $data_directory

        slurm_text="#!/bin/bash
#SBATCH --time=${JOB_TIME_HOURS}:00:00
#SBATCH --account=rrg-mdiamond
#SBATCH --array=1-${N_JOB}
#SBATCH --mem=4G
#SBATCH --job-name=musim-${run_name}
#SBATCH --constraint=[skylake|cascade]
#SBATCH --output=$basedir_data/upward/log/${run_name}_%a.out
echo current job id is: \$SLURM_ARRAY_TASK_ID
./upward_n_and_p/_start_single_run.sh $SIM_REPO_DIR $data_directory $events_per_job \$SLURM_ARRAY_TASK_ID $run_name $run_energy
"

        # Run or print the slurm commands
        if [[ "$SUBMIT" == "True" ]]; then
            echo "$slurm_text" >slurm_script_temp.txt
            sbatch slurm_script_temp.txt
            rm -f slurm_script_temp.txt
        elif [[ "$SUBMIT" == "Run" ]]; then
            echo "Running command: " "$slurm_text"
            echo "$slurm_text" -s 0 >slurm_script_temp.txt        
            for ((i = 0; i < N_JOB; i++)); do
                # Manually set SLURM_ARRAY_TASK_ID
                export SLURM_ARRAY_TASK_ID=$i
                bash slurm_script_temp.txt
            done
            rm -f slurm_script_temp.txt
        else
            echo Debug: slurm script to be sumbitted
            echo "-------------------------------------------------------"
            echo "$slurm_text"
            echo "-------------------------------------------------------"
        fi
    done
done
