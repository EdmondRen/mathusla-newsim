#!/bin/bash
#SBATCH --time=0:15:00
#SBATCH --account=rrg-mdiamond
#SBATCH --array=1
#SBATCH --mem=4G
#SBATCH --job-name=muLHCv
#SBATCH --output=/project/6075887/MATHUSLA/simulation/llp/2_sim_output/log/hxx_%a.out

DATA_DIR="/project/6075887/MATHUSLA/simulation/"
# SIM_REPO_DIR="/home/tomren/geant_projects/mathusla-newsim/"
SIM_REPO_DIR="/project/6035200/tomren/jupyter/mathusla-newsim"

NOSIE_RATE_PER_BAR_HZ="27.4"
N_COSMIC_PER_EVENT=0.86

SEED=1

FILENAME_LIST=()
# Set a default value for SLURM_ARRAY_TASK_ID if not running as a job array
if [ -z ${SLURM_ARRAY_TASK_ID+x} ]; then
    SLURM_ARRAY_TASK_ID=0
fi
echo current job id is: $SLURM_ARRAY_TASK_ID


# Loop through each element of source_array and append to destination_array
for p in "${DATA_DIR}/genie/genie_output_2025_10/g4_output_atmos/run_*/run_0_digionly_recon_skim.root"; do
    FILENAME_LIST+=($p)
done
merged_file_name="$DATA_DIR/genie/genie_output_2025_10/g4_output_atmos/recon_skim_merged.root"


# Run Sim+Digi+Recon on each filename in the list
# i=0
# for filename in "${FILENAME_LIST[@]}"; do
pushd $SIM_REPO_DIR/build
for i in {1..100}; do
    filename=${FILENAME_LIST[i]}
    echo "Processing file: $filename"

    current_file_name_recon=${FILENAME_LIST[i]}

    # Merge the reconstruction result
    if ((i > 1)); then
        ./merge $current_file_name_recon $merged_file_name "data"
        echo " Results merged to $merged_file_name"
    else
        echo $current_file_name_recon
        echo $merged_file_name
        cp  $current_file_name_recon $merged_file_name
    fi


    ((i++))
done
popd
