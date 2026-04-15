#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --account=rrg-mdiamond
#SBATCH --mem=16G
#SBATCH --job-name=musim-merge
#SBATCH --output=~/log/%x_%a.out

###
## muon: 
# SIM_DATA_DIR=/project/6049244/data/MATHUSLA/simulation_v2
# ./submit_skim_merge_files.sh "$SIM_DATA_DIR/cosmic/cosmic_mu/series_*_digi_recon.root" 1 2627 $SIM_DATA_DIR/cosmic/cosmic_mu/merged/
## Submitting n+p: 
# ./submit_skim_merge_files.sh "$SIM_DATA_DIR/cosmic/cosmic_np/series_*_digi_recon.root" 1 500 $SIM_DATA_DIR/cosmic/cosmic_np/merged/

SIM_REPO_DIR="/project/6035200/tomren/jupyter/mathusla-newsim"

## Usage
if [ $# -ne 4 ]; then
    echo "Merge files that match a wildcard, within a given start and end index."
    echo "Usage: $0 <FILE_PATH_WILDCARD> <START_INDEX> <END_INDEX> <MERGE_DIR>"
    exit 1
fi

FILE_PATH_WILDCARD=$1
START_INDEX=$2
END_INDEX=$3
MERGE_DIR=$4

## TEMPORARY DIRECTORIES 
if [[ -z "$SLURM_TMPDIR" ]]; then 
    echo "No SLURM_TMPDIR detected, using MERGE_DIR for temporary files."
    TEMP_DIR="$MERGE_DIR"
    SLURM_ARRAY_TASK_ID=0
else
    TEMP_DIR="$SLURM_TMPDIR"
fi

mkdir -p "$TEMP_DIR"
mkdir -p "$MERGE_DIR"

echo "Using temporary directory: $TEMP_DIR"
echo "Merging results will be stored in: $MERGE_DIR"

## Find all matching files
# Sort the files to ensure consistent order
FILES=($(ls $FILE_PATH_WILDCARD 2>/dev/null | sort))
NUM_FILES=${#FILES[@]}

if (( NUM_FILES == 0 )); then
    echo "No files found matching pattern: $FILE_PATH_WILDCARD"
    exit 1
fi

echo "Found $NUM_FILES matching files."

## Select files in the given range
SELECTED_FILES=("${FILES[@]:$START_INDEX:$((END_INDEX - START_INDEX + 1))}")

if (( ${#SELECTED_FILES[@]} == 0 )); then
    echo "No files in the specified range ($START_INDEX–$END_INDEX)"
    exit 1
fi

echo "Merging ${#SELECTED_FILES[@]} files (indices $START_INDEX–$END_INDEX)..."

## Merge files
OUT_FILE_final="$MERGE_DIR/merged_${START_INDEX}_${END_INDEX}.root"
OUT_FILE_tmp="$TEMP_DIR/merged_${START_INDEX}_${END_INDEX}.root"
OUT_FILE_skim_tmp="$TEMP_DIR/skim_tmp.root"

# Initialize merged file as the first one
cp "${SELECTED_FILES[0]}" "$OUT_FILE_tmp"

# Merge subsequent files
for ((i = 1; i < ${#SELECTED_FILES[@]}; i++)); do
    echo "Merging ${SELECTED_FILES[i]}..."
    $SIM_REPO_DIR/build/skim ${SELECTED_FILES[i]} --output $OUT_FILE_skim_tmp -p 1000
    $SIM_REPO_DIR/build/merge $OUT_FILE_skim_tmp "$OUT_FILE_tmp" "data"
done

# Move the merged result to the final directory
mv "$OUT_FILE_tmp" "$OUT_FILE_final"

echo "Merged file saved to: $OUT_FILE_final"
echo "Done."
