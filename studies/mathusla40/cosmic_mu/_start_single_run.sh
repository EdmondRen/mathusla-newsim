# Usage: ./_start_single_run.sh  SIM_REPO_DIR OUTPUT_DIR EVENTS_PER_RUN RUN_NUMBER 
# example: ./_start_single_run.sh ~/geant_projects/mathusla-newsim $data_dir 10000 1

SIM_REPO_DIR=$1
OUTPUT_DIR=$2
EVENTS_PER_RUN=$3
RUN_NUMBER=$4

SEED=$RUN_NUMBER
NOSIE_RATE_PER_BAR_HZ="27.4"

cd $SIM_REPO_DIR/build

## Simulation
env G4RUN_MANAGER_TYPE=Serial ./simulation -r $RUN_NUMBER \
    -s $SEED \
    -m $SIM_REPO_DIR/studies/mathusla40/cosmic_mu/scripts/run_cry_muon_mathusla40.mac,events,$EVENTS_PER_RUN \
    -o $OUTPUT_DIR

## Digitizer
./digitizer $OUTPUT_DIR/run_${RUN_NUMBER}.root \
    -s $SEED \
    -p 10000 \
    -n $NOSIE_RATE_PER_BAR_HZ


## Reconstruction
# -R: discard simulation truth except for the generator status to reproduce the event later
# -p 10000: print every 100000 events
./tracker $OUTPUT_DIR/run_${RUN_NUMBER}_digi.root \
    -o $OUTPUT_DIR/run_${RUN_NUMBER}_digi_recon_all.root \
    -r $OUTPUT_DIR/run_${RUN_NUMBER}.root \
    -k 2 \
    -p 10000    

## Reconstruction, only save events with upward vertices
# -k 3: save only events with UPWARD vertices
# -R: discard simulation truth except for the generator status to reproduce the event later
# -p 10000: print every 100000 events
./tracker $OUTPUT_DIR/run_${RUN_NUMBER}_digi.root \
    -o $OUTPUT_DIR/run_${RUN_NUMBER}_digi_recon_upward.root \
    -r $OUTPUT_DIR/run_${RUN_NUMBER}.root \
    -k 3 \
    -p 10000


