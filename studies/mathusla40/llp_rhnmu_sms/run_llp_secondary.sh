#!/bin/bash

DATA_DIR="/home/tomren/geant_projects/musim_test/llp_sms/"
SIM_REPO_DIR="/home/tomren/geant_projects/mathusla-newsim/"

NOSIE_RATE_PER_BAR_HZ="27.4"
SEED=1

momentum_list=(deweighted_LLP_0.0774597.root deweighted_LLP_0.274545.root deweighted_LLP_0.678553.root deweighted_LLP_0.878742.root deweighted_LLP_1.9.root deweighted_LLP_3.47544.root deweighted_LLP_3.74996.root deweighted_LLP_4.7.root deweighted_LLP_RHN_Umu_0.1.root deweighted_LLP_RHN_Umu_0.273842.root deweighted_LLP_RHN_Umu_0.649382.root deweighted_LLP_RHN_Umu_0.865964.root deweighted_LLP_RHN_Umu_1.90365.root deweighted_LLP_RHN_Umu_3.08517.root deweighted_LLP_RHN_Umu_3.85.root deweighted_LLP_RHN_Umu_4.8028.root )
FILENAME_LIST=()



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

N_COSMIC_PER_EVENT=0.86

############################################################################################
# 1. Generate cosmic simulation for each signal samples

# NITER=${#momentum_list[@]}
# pushd $SIM_REPO_DIR/build

# NEVENT_COSMIC=30000
# for ((i = 0; i < NITER; i++)); do

#     RUN_NUMBER=$i
#     SEED=$RUN_NUMBER
#     ## Simulation
#     env G4RUN_MANAGER_TYPE=Serial ./simulation -r $RUN_NUMBER \
#         -s $SEED \
#         -m $SIM_REPO_DIR/studies/mathusla40/llp_hxx/run_cry_all_mathusla40.mac,events,$NEVENT_COSMIC \
#         -o $DATA_DIR/cosmic

#     ## Digitizer
#     ./digitizer $DATA_DIR/cosmic/run_${RUN_NUMBER}.root \
#         -s $SEED \
#         -p 3000 \
#         -n 0
# done
# popd

############################################################################################
# 2. Signal
# Loop through each element of source_array and append to destination_array
for p in "${momentum_list[@]}"; do
    FILENAME_LIST+=("$DATA_DIR/scripts/$p")
done


# Run Sim+Digi+Recon on each filename in the list
i=0
for filename in "${FILENAME_LIST[@]}"; do

    pushd $SIM_REPO_DIR/build
    RUN_NUMBER=$i

    ## Simulation
    echo
    echo $SIM_REPO_DIR/studies/mathusla40/llp_rhnmu_sms/g4_filereader.mac,filename,$filename
    env G4RUN_MANAGER_TYPE=Serial ./simulation -r $RUN_NUMBER \
        -s $SEED \
        -m $SIM_REPO_DIR/studies/mathusla40/llp_rhnmu_sms/g4_filereader.mac,filename,$filename \
        -o $DATA_DIR

    ## Digitizer
    ./digitizer $DATA_DIR/run_${RUN_NUMBER}.root \
        -s $SEED \
        -p 1000 \
        -n $NOSIE_RATE_PER_BAR_HZ

    ## Add cosmic ray
    ./attach_cosmic $DATA_DIR/run_${RUN_NUMBER}_digi.root \
        $DATA_DIR/cosmic/run_${RUN_NUMBER}_digi.root \
        $DATA_DIR/run_${RUN_NUMBER}_digi_cosmic.root \
        $N_COSMIC_PER_EVENT


    ## Reconstruction on vanila digitzation (no cosmic)
    # -k 0: save all events
    # -R: discard simulation truth except for the generator status to reproduce the event later
    # -p 100: print every 100 events
    ./tracker $DATA_DIR/run_${RUN_NUMBER}_digi.root \
        -r $DATA_DIR/run_${RUN_NUMBER}.root \
        -k 0 \
        -p 1000

    ## Reconstruction on cosmic-added digitzation
    # -k 0: save all events
    # -R: discard simulation truth except for the generator status to reproduce the event later
    # -p 100: print every 100 events
    ./tracker $DATA_DIR/run_${RUN_NUMBER}_digi_cosmic.root \
        -r $DATA_DIR/run_${RUN_NUMBER}.root \
        -k 0 \
        -p 1000 -R

    ## Select events with at least one vertex in the box
    ./skim $DATA_DIR/run_${RUN_NUMBER}_digi_cosmic_recon.root -p 1000

    popd

    ((i++))
done
