#!/bin/bash
#SBATCH --time=0:15:00
#SBATCH --account=rrg-mdiamond
#SBATCH --array=1-120
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


############################################################################################
# Simulation of LHC neutrino


############################################################################################
# 2. Signal
# Loop through each element of source_array and append to destination_array
for p in "${DATA_DIR}/genie/genie_output_2025_10/merged/*.root"; do
    FILENAME_LIST+=($p)
done

# Run Sim+Digi+Recon on each filename in the list
# i=0
# for filename in "${FILENAME_LIST[@]}"; do
for i in $SLURM_ARRAY_TASK_ID; do
    filename=${FILENAME_LIST[i]}
    echo "Processing file: $filename"

    pushd $SIM_REPO_DIR/build
    filename_base=$(basename -- "$filename")
    filename_noext="${filename_base%.*}"
    RUN_NUMBER=0

    data_dir_this="${DATA_DIR}/genie/genie_output_2025_10/g4_output/${filename_noext:0:-7}/"
    mkdir -p $data_dir_this

    ## Simulation
    env G4RUN_MANAGER_TYPE=Serial ./simulation -r $RUN_NUMBER \
    ./simulation -r 0 \
        -s $SEED \
        -m $SIM_REPO_DIR/studies/mathusla40/llp_hxx/g4_llp_hxx_arg.mac,sourcefile,${filename} \
        -o $data_dir_this

    if [ -f $data_dir_this/run_${RUN_NUMBER}_digi.root ]; then
        echo "Already processed, skip simulation and digitization."
        exit 1
    fi    

    ## Digitizer
    ./digitizer $data_dir_this/run_${RUN_NUMBER}.root \
        -s $SEED \
        -p 100 \
        -n $NOSIE_RATE_PER_BAR_HZ


    ## Reconstruction
    # -k 0: save all events
    # -p 100: print every 100 events
    ./tracker $data_dir_this/run_${RUN_NUMBER}_digi.root \
        -r $data_dir_this/run_${RUN_NUMBER}.root \
        -k 0 \
        -p 100


    ## Reconstruction
    # -k 0: save all events
    # -R: discard simulation truth except for the generator status to reproduce the event later
    # -p 100: print every 100 events
    ./tracker $data_dir_this/run_${RUN_NUMBER}_digi.root \
        -o $data_dir_this/run_${RUN_NUMBER}_digionly_recon.root \
        -r $data_dir_this/run_${RUN_NUMBER}.root \
        -k 3 -R \
        -p 100


    ## Select events with at least one vertex in the box
    ./skim $data_dir_this/run_${RUN_NUMBER}_digionly_recon.root -p 1000

    popd

    ((i++))
done
