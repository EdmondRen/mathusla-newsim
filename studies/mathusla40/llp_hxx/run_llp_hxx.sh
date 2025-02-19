#!/bin/bash


event_source_path="/home/tomren/geant_projects/temp/"
momentum_list=( 25 )
# momentum_list=( 15 25 35 45 55)
FILENAME_LIST=()

# Loop through each element of source_array and append to destination_array
for p in "${momentum_list[@]}"; do
    FILENAME_LIST+=("$event_source_path/LLP_bb_${p}_ctau_1000.root")
done


EXAMPLE_MACRO="g4_llp_hxx_example.mac"

# Check if the input file exists
if [ ! -f "$EXAMPLE_MACRO" ]; then
    echo "Error: $EXAMPLE_MACRO not found!"
    exit 1
fi

# Process each filename in the list
while IFS= read -r filename; do
    if [ -n "$filename" ]; then
        OUTPUT_FILE="g4_llp_hxx_${filename}_GeV.txt"
        sed -E "s|(/gen/recreate/pathname )[^ ]+|\1${filename}|g" "$EXAMPLE_MACRO" > "$OUTPUT_FILE"
        echo "Created: $OUTPUT_FILE"
    fi
done < "$FILENAME_LIST"
