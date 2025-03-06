#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --account=ctb-ut-atlas
#SBATCH --array=0
#SBATCH --mem=8G
#SBATCH --job-name=merge_lhc
#SBATCH --output=/home/tomren/log/merge_%a.out

## It is fast to process (0.5s each), but IO takes most of the time. 70 files are done in 5 minutes

##===============================================================================
## Those are the only parameters you should change:
basedir="/project/6049244/data/MATHUSLA/simulation/"
# basedir="/home/tomren/geant_projects/musim_test/cedar"
##===============================================================================

SIM_REPO_DIR=$(realpath ../..)
cd $SIM_REPO_DIR/build

run_names=(run-2024-07-cosmic-neutrino)

if [ -z "$SLURM_TMPDIR" ]; then
    echo "No SLURM_TMPDIR, running locally"
    NITER=${#run_names[@]}
    echo "N iteration: $NITER"
    SLURM_TMPDIR=$(pwd)
else
    NITER=1
fi

for ((i = 0; i < NITER; i++)); do

    if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
        run_name=${run_names[i]}
    else
        run_name=${run_names[$SLURM_ARRAY_TASK_ID]}
    fi

    data_directory=$basedir/$run_name
    mkdir -p $data_directory/v2_recon_and_skim

    # Get a list of reconstruction files
    # echo "Searching for files: $data_directory/DigiOutput/*/0/stat0.root"
    # files=($(ls $data_directory/DigiOutput/*/0/stat0.root))
    # echo "  ${#files[@]} found"
    # processed=0

    # Loop all files, select events and merge the output
    # Save output to temp file
    filename_merged=$data_directory/v2_recon_and_skim/merged.root
    filename_merged_tmp=$SLURM_TMPDIR/merged.root
    \rm $filename_merged
    # for filename in "${files[@]}"; do
    #     ((processed += 1))
    nfiles=10
    for ((processed=0;processed<=nfiles;processed++)); do
        echo Reconstruction
        filename=$data_directory/DigiOutput/run_$processed/0/stat0.root

        [ ! -f $filename ] && continue

        echo processing $processed/$nfiles "$filename"

        ## Reonstruction
        new_filename_recon="${filename%.*}_recon.${filename##*.}"
        new_filename_recon=$SLURM_TMPDIR/$(basename $new_filename_recon)
        ./tracker_old $filename --output $new_filename_recon -p 1000 -k 3 

        ## Pre-selection
        echo
        echo
        new_filename_skim="${filename%.*}_recon_skim.${filename##*.}"
        new_filename_skim=$SLURM_TMPDIR/$(basename $new_filename_skim)

        ./skim_old $new_filename_recon --output $new_filename_skim -p 1000

        echo
        echo "Merge file"
        echo ./merge $new_filename_skim $filename_merged_tmp data
        ./merge $new_filename_skim $filename_merged_tmp data
    done

    # Move all files back to perm storage
    echo "Move file back from SLURM_TMPDIR to $filename_merged"
    \mv $filename_merged_tmp $filename_merged
    chgrp rrg-mdiamond $filename_merged
    echo "Finished copying merged file"
    # for filename in "${files[@]}"; do
    #     new_filename="${filename%.*}_skim.${filename##*.}"
    #     new_filename=$SLURM_TMPDIR/$(basename $new_filename)
    #     \mv $new_filename $data_directory/v2_recon_and_skim/
    #     echo $filename
    # done

done
