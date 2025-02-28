#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --account=ctb-ut-atlas
#SBATCH --array=0-1
#SBATCH --mem=8G
#SBATCH --job-name=merge
#SBATCH --output=/home/tomren/log/merge_%a.out

## It is fast to process (0.5s each), but IO takes most of the time. 70 files are done in 5 minutes

##===============================================================================
## Those are the only parameters you should change:
basedir="/project/6049244/data/MATHUSLA/simulation_v2"
# basedir="/home/tomren/geant_projects/musim_test/cedar"
##===============================================================================

SIM_REPO_DIR=$(realpath ../..)
cd $SIM_REPO_DIR/build

run_names=(backup_cosmic_p backup_cosmic_n_partial)
# run_names=(cosmic_p cosmic_n)
# run_names=(cosmic_p)

if [ -z "$SLURM_TMPDIR" ]; then
    echo "No SLURM_TMPDIR, running locally"
    NITER=${#run_names[@]}
    echo "N iteration: $NITER"
    SLURM_TMPDIR=`pwd`
else
    NITER=1
fi

for ((i = 0; i < NITER; i++)); do

    if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
        run_name=${run_names[i]}
    else
        run_name=${run_names[$SLURM_ARRAY_TASK_ID]}
    fi

    data_directory=$basedir/cosmic/$run_name
    mkdir -p $data_directory/skim

    # Get a list of reconstruction files
    echo "Searching for files: $data_directory/series_*_digi_recon.root"
    files=($(ls $data_directory/series_*_digi_recon.root))
    echo "  ${#files[@]} found"
    processed=0

    # Loop all files, select events and merge the output
    # Save output to temp file
    filename_merged=$data_directory/skim/merged.root
    filename_merged_tmp=$SLURM_TMPDIR/merged.root
    \rm $filename_merged
    for filename in "${files[@]}"; do
        ((processed += 1))
        echo
        echo
        echo processing $processed/"${#files[@]}" "$filename"
        new_filename="${filename%.*}_skim.${filename##*.}"
        new_filename=$SLURM_TMPDIR/$(basename $new_filename)

        ./skim $filename --output $new_filename -p 1000

        echo
        echo "Merge file"
        echo ./merge $new_filename $filename_merged_tmp data
        ./merge $new_filename $filename_merged_tmp data
    done

    # Move all files back to perm storage
    echo "Move file back from SLURM_TMPDIR"
    \mv $filename_merged_tmp $filename_merged
    echo "Finished copying merged file"
    for filename in "${files[@]}"; do
        new_filename="${filename%.*}_skim.${filename##*.}"
        new_filename=$SLURM_TMPDIR/$(basename $new_filename)
        \mv $new_filename $data_directory/skim/
        echo $filename
    done

done
