# Usage: ./_start_single_run.sh  SIM_REPO_DIR OUTPUT_DIR EVENTS_PER_RUN RUN_NUMBER  PARTICLE ENERGY_MeV
# example: ./_start_single_run.sh ~/geant_projects/mathusla-newsim $data_dir 1000000 1 proton 10000

SIM_REPO_DIR=$1
OUTPUT_DIR=$2
EVENTS_PER_RUN=$3
RUN_NUMBER=$4
PARTICLE=$5
ENERGY=$6

SEED=$RUN_NUMBER
NOSIE_RATE_PER_BAR_HZ="27.4"

cd $SIM_REPO_DIR/build

## Simulation
echo 
echo $SIM_REPO_DIR/studies/mathusla40/upward_n_and_p/scripts/run_gun.mac,particle,$PARTICLE,energy,$ENERGY,nevents,$EVENTS_PER_RUN
echo 
env G4RUN_MANAGER_TYPE=Serial ./simulation -r $RUN_NUMBER \
    -s $SEED \
    -m $SIM_REPO_DIR/studies/mathusla40/upward_n_and_p/scripts/run_gun.mac,particle,$PARTICLE,energy,$ENERGY,nevents,$EVENTS_PER_RUN \
    -o $OUTPUT_DIR

## Digitizer
./digitizer $OUTPUT_DIR/run_${RUN_NUMBER}.root \
    -s $SEED \
    -p 30000 \
    -n $NOSIE_RATE_PER_BAR_HZ

## Reconstruction
# -k 3: save only events with UPWARD vertices
# -R: discard simulation truth except for the generator status to reproduce the event later
# -p 10000: print every 100000 events
./tracker $OUTPUT_DIR/run_${RUN_NUMBER}_digi.root \
    -r $OUTPUT_DIR/run_${RUN_NUMBER}.root \
    -k 3 \
    -p 30000

## Data skimming
./skim $OUTPUT_DIR/run_${RUN_NUMBER}_digi_recon.root --output $OUTPUT_DIR/run_${RUN_NUMBER}_digi_recon_skim.root -p 2000

## Delete raw file, only keep the one after skimming
\rm $OUTPUT_DIR/run_${RUN_NUMBER}.root
\rm $OUTPUT_DIR/run_${RUN_NUMBER}_digi.root
\rm $OUTPUT_DIR/run_${RUN_NUMBER}_digi_recon.root