#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --account=rrg-mdiamond
#SBATCH --array=0-39
#SBATCH --mem=4G
#SBATCH --job-name=muLLP
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
# Simulation of pp->XX, X->bb signal
# 1. Generate cosmic ray simulations+digi that matches the number of signals
#   - Cosmic ray time is sampled uniformly +/- 1000 ns around 0.
#   - It can be set via /gen/cry/offset_t_low and /gen/cry/offset_t_high.
#   - [+/- 1000 ns] matches the default used in digitizer.
#   - 30k events will results in about 13k digitized events that have any hits in detector.
#   - CRY recorded 0.0302s real world time for 30k, so there about 13e3*2e-6/0.03 = **0.86** cosmic events per signal event.
#   - This number is used when combining cosmic events with signal
# 3. Run simulation + digitization on signal
# 4. Combine the digitzation results from cosmic to signal with given probability
#   - The average cosmic event per signal sample is **0.86**
# 5. Run reconsturction on combined files


############################################################################################
# 1. Generate cosmic simulation
## Refer to ../trigger/run_cosmic_all.sh for details of cosmic simulation

############################################################################################
# 2. Signal
# Loop through each element of source_array and append to destination_array
for p in "${DATA_DIR}/llp/1_decay_files_transformed/MATHUSLA_LLPfiles_HXX/*.root"; do
    FILENAME_LIST+=($p)
done

EXAMPLE_MACRO="g4_llp_hxx_example.mac"
# Check if the input file exists
if [ ! -f "$EXAMPLE_MACRO" ]; then
    echo "Error: $EXAMPLE_MACRO not found!"
    exit 1
fi

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

    data_dir_this="${DATA_DIR}/llp/2_sim_output/hxx/${filename_noext}/"
    mkdir -p $data_dir_this
    data_dir_cosmic="${DATA_DIR}/cosmic/cosmic_all/"

    ## Simulation
    # env G4RUN_MANAGER_TYPE=Serial ./simulation -r $RUN_NUMBER \
    # ./simulation -r 0 \
    #     -s $SEED \
    #     -m $SIM_REPO_DIR/studies/mathusla40/llp_hxx/g4_llp_hxx_arg.mac,sourcefile,${filename} \
    #     -o $data_dir_this

    # if [ -f $data_dir_this/run_${RUN_NUMBER}_digi.root ]; then
    #     echo "Already processed, skip simulation and digitization."
    #     exit 1
    # fi    

    # ## Digitizer
    # ./digitizer $data_dir_this/run_${RUN_NUMBER}.root \
    #     -s $SEED \
    #     -p 100 \
    #     -n $NOSIE_RATE_PER_BAR_HZ

    # ## Add cosmic ray
    # ./attach_cosmic $data_dir_this/run_${RUN_NUMBER}_digi.root \
    #     $data_dir_cosmic/run_${i}_digi.root \
    #     $data_dir_this/run_${RUN_NUMBER}_digi_cosmic.root \
    #     $N_COSMIC_PER_EVENT


    # ## Reconstruction on vanila digitzation (no cosmic)
    # # -k 0: save all events
    # # -R: discard simulation truth except for the generator status to reproduce the event later
    # # -p 100: print every 100 events
    # ./tracker $data_dir_this/run_${RUN_NUMBER}_digi.root \
    #     -r $data_dir_this/run_${RUN_NUMBER}.root \
    #     -k 0 \
    #     -p 100

    ## Reconstruction on cosmic-added digitzation
    # -k 0: save all events
    # -R: discard simulation truth except for the generator status to reproduce the event later
    # -p 100: print every 100 events
    ./tracker $data_dir_this/run_${RUN_NUMBER}_digi_cosmic.root \
        -r $data_dir_this/run_${RUN_NUMBER}.root \
        -k 0 \
        -p 100

    ## Select events with at least one vertex in the box
    ./skim $data_dir_this/run_${RUN_NUMBER}_digi_cosmic_recon.root -p 1000

    popd

    ((i++))
done
