# Usage: ./_start_single_run.sh  SIM_REPO_DIR OUTPUT_DIR EVENTS_PER_RUN SERIES_NUMBER  NRUN_PER_SERIES
#   example: ./_start_run_series.sh ~/geant_projects/mathusla-newsim $data_dir 200000 0 10 > 1.log &
# [1] 425088

SIM_REPO_DIR=$1
OUTPUT_DIR=$2
EVENTS_PER_RUN=$3
SERIES_NUMBER=$4
NRUN_PER_SERIES=$5

cd $SIM_REPO_DIR/build
series_output_dir=$OUTPUT_DIR/series_$SERIES_NUMBER

# Loop each run number
merged_file_name=$series_output_dir/series_${SERIES_NUMBER}_digi_recon.root
for irun in $(seq 1 $NRUN_PER_SERIES); do
    RUN=$((irun - 1))
    # Print the progress
    echo "Running series $SERIES_NUMBER, run number $RUN/$NRUN_PER_SERIES"

    actual_run_number=$((NRUN_PER_SERIES * SERIES_NUMBER + RUN))
    mkdir -p $series_output_dir

    # Run simulation and recon
    bash $SIM_REPO_DIR/studies/mathusla40/cosmic_p/_start_single_run.sh $SIM_REPO_DIR $series_output_dir $EVENTS_PER_RUN $actual_run_number

    # filenames
    current_file_name_truth=$series_output_dir/run_${actual_run_number}.root
    current_file_name_digi=$series_output_dir/run_${actual_run_number}_digi.root
    current_file_name_recon=$series_output_dir/run_${actual_run_number}_digi_recon.root

    # Merge the reconstruction result
    if ((irun > 1)); then
        ./merge $current_file_name_recon $merged_file_name
        rm -f $current_file_name_recon
    else
        echo $current_file_name_recon
        echo $merged_file_name
        mv -f $current_file_name_recon $merged_file_name
    fi

    # Delete the sim truth and digi
    rm -f $current_file_name_truth
    rm -f $current_file_name_digi

    echo "Results merged. Sim files deleted"
done
