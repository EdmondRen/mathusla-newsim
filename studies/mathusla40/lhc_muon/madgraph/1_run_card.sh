#!/bin/bash

# A script to prepare a MadGraph parcard and launch the process.
# It takes three arguments:
# 1. The path to a template parcard.
# 2. The desired name/path for the output directory.
# 3. The random seed for the run.

# --- 1. Argument Validation ---
if [ "$#" -ne 3 ]; then
    echo "ERROR: Invalid number of arguments."
    echo "Usage: $0 <path_to_parcard> <output_directory> <seed>"
    exit 1
fi

PARCARD_PATH=$1
OUTPUT_DIR=$2
SEED=$3

# Check if the parcard file exists
if [ ! -f "$PARCARD_PATH" ]; then
    echo "ERROR: Parcard file not found at '$PARCARD_PATH'"
    exit 1
fi

# --- 2. Setup Environment ---
echo "--- MadGraph Job Setup ---"
echo "Template Card: $PARCARD_PATH"
echo "Output Directory: $OUTPUT_DIR"
echo "Random Seed: $SEED"
echo "--------------------------"

# Create the output directory. The -p flag prevents errors if it already exists.
mkdir -p "$OUTPUT_DIR"

# Define the path for the new, modified card that will be used for the run.
MODIFIED_CARD="$OUTPUT_DIR/run_card.txt"

# Copy the template card to the output directory.
cp "$PARCARD_PATH" "$MODIFIED_CARD"
echo "Copied template card to $MODIFIED_CARD"


# --- 3. Modify the Parcard ---
# Use sed (stream editor) to find and replace the placeholder lines in-place.
# We use '|' as the delimiter instead of '/' to avoid conflicts with directory paths.
# s|pattern|replacement|g
echo "Modifying the run card with user-defined parameters..."

sed -i "s|^output.*|output $OUTPUT_DIR|" "$MODIFIED_CARD"
sed -i "s|^set iseed.*|set iseed = $SEED|" "$MODIFIED_CARD"

echo "Card modification complete."


# --- 4. Launch MadGraph ---
echo "Launching MadGraph (mg5_aMC)..."
echo "This may take a significant amount of time."

# Execute MadGraph with the newly created and modified parameter card.
# Make sure 'mg5_aMC' is in your system's PATH.
mg5_aMC "$MODIFIED_CARD"

echo "--- MadGraph process finished. ---"

rm py.py


