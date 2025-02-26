#!/bin/bash


DATA_DIR="/home/tomren/geant_projects/temp/"
SIM_REPO_DIR="/home/tomren/geant_projects/mathusla-newsim/"

NOSIE_RATE_PER_BAR_HZ="27.4"
SEED=1

# momentum_list=( 15 )
momentum_list=( 15 25 35 45 55)
FILENAME_LIST=()

# Loop through each element of source_array and append to destination_array
for p in "${momentum_list[@]}"; do
    FILENAME_LIST+=("$DATA_DIR/LLP_bb_${p}_ctau_1000.root")
done

EXAMPLE_MACRO="g4_llp_hxx_example.mac"
# Check if the input file exists
if [ ! -f "$EXAMPLE_MACRO" ]; then
    echo "Error: $EXAMPLE_MACRO not found!"
    exit 1
fi

# Process each filename in the list
i=0
for filename in "${FILENAME_LIST[@]}"; do

    # Make a macro for geant
    if [ -n "$filename" ]; then
        OUTPUT_FILE="g4_llp_hxx_${momentum_list[i]}_GeV.mac"
        sed -E "s|(/gen/recreate/pathname )[^ ]+|\1${filename}|g" "$EXAMPLE_MACRO" > "$OUTPUT_FILE"
        echo "Created: $OUTPUT_FILE"
    fi

    pushd $SIM_REPO_DIR/build
    RUN_NUMBER=${momentum_list[i]}
    ## Simulation
    echo $SIM_REPO_DIR/studies/mathusla40/llp_hxx/g4_llp_hxx_${momentum_list[i]}_GeV.mac 
    env G4RUN_MANAGER_TYPE=Serial ./simulation -r $RUN_NUMBER \
        -s $SEED \
        -m $SIM_REPO_DIR/studies/mathusla40/llp_hxx/g4_llp_hxx_${momentum_list[i]}_GeV.mac \
        -o $DATA_DIR
    
    ## Digitizer
    ./digitizer $DATA_DIR/run_${RUN_NUMBER}.root \
        -s $SEED \
        -p 100 \
        -n $NOSIE_RATE_PER_BAR_HZ
    
    ## Reconstruction
    # -k 3: save only events with UPWARD vertices
    # -R: discard simulation truth except for the generator status to reproduce the event later
    # -p 10000: print every 100000 events
    ./tracker $DATA_DIR/run_${RUN_NUMBER}_digi.root \
        -r $DATA_DIR/run_${RUN_NUMBER}.root \
        -k 0 \
        -p 50

    popd

    (( i++ ))    
done