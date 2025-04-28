#!/bin/bash

DATA_DIR="/home/tomren/geant_projects/musim_test/cosmic_all/"
SIM_REPO_DIR="/home/tomren/geant_projects/mathusla-newsim/"

NOSIE_RATE_PER_BAR_HZ="27.4"
SEED=1
NITER=10


############################################################################################
# Simulation of all cosmic ray
# - Generator: cry
#   - Sample all supported particles
############################################################################################


pushd $SIM_REPO_DIR/build

NEVENT_COSMIC=500000
for ((i = 0; i < NITER; i++)); do

    RUN_NUMBER=$i
    SEED=$RUN_NUMBER
    ## Simulation
    env G4RUN_MANAGER_TYPE=Serial ./simulation -r $RUN_NUMBER \
        -s $SEED \
        -m $SIM_REPO_DIR/studies/mathusla40/trigger/g4config_cry_all_mathusla40_abstime.mac,events,$NEVENT_COSMIC \
        -o $DATA_DIR


    ## Digitizer
    ./digitizer $DATA_DIR/run_${RUN_NUMBER}.root \
        -s $SEED \
        -p 100 \
        -n $NOSIE_RATE_PER_BAR_HZ


    ## Reconstruction on digitzation
    # -k 0: save all events
    # -p 1000: print every 1000 events
    ./tracker $DATA_DIR/run_${RUN_NUMBER}_digi.root \
        -r $DATA_DIR/run_${RUN_NUMBER}.root \
        -k 0 \
        -p 1000

done
popd
