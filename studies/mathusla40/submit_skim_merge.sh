##===============================================================================
## Those are the only parameters you should change:
# basedir="/project/6049244/data/MATHUSLA/simulation_v2"
basedir="/home/tomren/geant_projects/musim_test/cedar"
##===============================================================================

SIM_REPO_DIR=$(realpath ../..)
cd $SIM_REPO_DIR/build

# run_names=(cosmic_p)
run_names=(cosmic_p cosmic_n)

for run_name in "${run_names[@]}"; do

    data_directory=$basedir/cosmic/$run_name
    mkdir -p $data_directory/skim

    # Get a list of reconstruction files
    files=($(ls $data_directory/series_*_digi_recon.root))
    processed=0

    # Loop all files, select events and merge the output
    filename_merged=$data_directory/skim/merged.root
    \rm $filename_merged
    for filename in "${files[@]}"; do
        echo
        echo
        echo processing $processed/"${#files[@]}" "$filename"
        new_filename="${filename%.*}_skim.${filename##*.}"
        new_filename=$data_directory/skim/$(basename $new_filename)

        ./skim $filename --output $new_filename

        echo
        echo "Merge file"
        echo ./merge $new_filename $filename_merged data
        ./merge $new_filename $filename_merged data
    done

done
