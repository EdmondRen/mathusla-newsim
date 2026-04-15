#!/bin/bash
#SBATCH --time=2:30:00
#SBATCH --account=rrg-mdiamond
#SBATCH --array=10-99
#SBATCH --mem=4G
#SBATCH --job-name=mucosmic_all
#SBATCH --output=/project/6075887/MATHUSLA/simulation/cosmic/cosmic_all/log/%a.out

DATA_DIR="/project/6075887/MATHUSLA/simulation/cosmic/cosmic_all_sidewall"
SIM_REPO_DIR="/project/6035200/tomren/jupyter/mathusla-newsim"

NOSIE_RATE_PER_BAR_HZ="27.4"
SEED=1
NITER=1

# Set a default value for SLURM_ARRAY_TASK_ID if not running as a job array
if [ -z ${SLURM_ARRAY_TASK_ID+x} ]; then
    SLURM_ARRAY_TASK_ID=0
fi
echo current job id is: $SLURM_ARRAY_TASK_ID

############################################################################################
# Simulation of all cosmic ray
# - Generator: cry
#   - Sample all supported particles
############################################################################################


pushd $SIM_REPO_DIR/build

NEVENT_COSMIC=500000
for ((i = SLURM_ARRAY_TASK_ID*NITER; i < (SLURM_ARRAY_TASK_ID+1)*NITER; i++)); do

    RUN_NUMBER=$i
    SEED=$RUN_NUMBER
    ## Simulation
    #env G4RUN_MANAGER_TYPE=Serial ./simulation -r $RUN_NUMBER \
    ./simulation -r $RUN_NUMBER \
        -s $SEED \
        -m $SIM_REPO_DIR/studies/mathusla40/trigger/g4config_cry_all_mathusla40_abstime_sidewall.mac,events,$NEVENT_COSMIC \
        -o $DATA_DIR


    ## Digitizer
    ./digitizer $DATA_DIR/run_${RUN_NUMBER}.root \
        -s $SEED \
        -p 5000 \
        -n $NOSIE_RATE_PER_BAR_HZ

    ## Reconstruction on digitzation
    # -k 0: save all events
    # -p 1000: print every 1000 events
    ./tracker $DATA_DIR/run_${RUN_NUMBER}_digi.root \
        -r $DATA_DIR/run_${RUN_NUMBER}.root \
        -k 0 \
        -p 5000

done
popd
