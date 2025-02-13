
SIM_REPO_DIR=""
SIM_TEMP_DIR=""
MERGE_DIR=""
N_EVENTS_PER_RUN=""
N_RUNS_PER_SERIES=""
SERIES_NUMBER=""
CLEAN_FILES="n"

usage() { 
    echo "Usage: $0 [-p SIM_REPO_DIR] [-o SIM_TEMP_DIR] [-m MERGE_DIR] [-e N_EVENTS_PER_RUN] [-r N_RUNS_PER_SERIES] [-s SERIES_NUMBER] [-c y or n]" 1>&2; 
    exit 1; }
while getopts p:o:m:e:r:s:c: flag
do
    case "${flag}" in
        p) SIM_REPO_DIR=${OPTARG};;
        o) SIM_TEMP_DIR=${OPTARG};;
        m) MERGE_DIR=${OPTARG};;
        e) N_EVENTS_PER_RUN=${OPTARG};;
        r) N_RUNS_PER_SERIES=${OPTARG};;
        s) SERIES_NUMBER=${OPTARG};;
        c) CLEAN_FILES=${OPTARG};;
        h) usage;;
    esac
done

# Shift positional arguments after options
shift $((OPTIND - 1))

# TEMPORARY DIRECTORIES 
if [[ -z "$SLURM_TMPDIR" && -z "$SIM_TEMP_DIR" ]]; then 
    echo "No SLURM_TMPDIR, and no output_dir, save output to MERGE_DIR"
    SIM_TEMP_DIR=$MERGE_DIR/series_$SERIES_NUMBER
    mkdir -p $SIM_TEMP_DIR
elif [[ -z "$SIM_TEMP_DIR" ]]; then
    echo "No output_dir, use SLURM_TMPDIR to temporarily save output"
    SIM_TEMP_DIR=$SLURM_TMPDIR/series_$SERIES_NUMBER
    mkdir -p $SIM_TEMP_DIR

fi


cd $SIM_REPO_DIR/build

# Loop each run number
merged_file_name=$MERGE_DIR/series_${SERIES_NUMBER}_digi_recon.root
for irun in $(seq 1 $N_RUNS_PER_SERIES); do
    RUN=$((irun - 1))
    # Print the progress
    echo "Running $((RUN+1))/$N_RUNS_PER_SERIES of series #$SERIES_NUMBER"

    actual_run_number=$((N_RUNS_PER_SERIES * SERIES_NUMBER + RUN))
    mkdir -p $SIM_TEMP_DIR

    # Run simulation and recon
    bash $SIM_REPO_DIR/studies/mathusla40/cosmic_n/_start_single_run.sh $SIM_REPO_DIR $SIM_TEMP_DIR $N_EVENTS_PER_RUN $actual_run_number

    # filenames
    current_file_name_truth=$SIM_TEMP_DIR/run_${actual_run_number}.root
    current_file_name_digi=$SIM_TEMP_DIR/run_${actual_run_number}_digi.root
    current_file_name_recon=$SIM_TEMP_DIR/run_${actual_run_number}_digi_recon.root

    # Merge the reconstruction result
    if ((irun > 1)); then
        ./merge $current_file_name_recon $merged_file_name "data"
        echo " Results merged to $merged_file_name"
    else
        echo $current_file_name_recon
        echo $merged_file_name
        cp  $current_file_name_recon $merged_file_name
    fi

    # Delete the sim truth and digi
    if [[ ! -z "$CLEAN_FILES" && "$CLEAN_FILES" == "y" ]]; then
        rm -f $current_file_name_truth
        rm -f $current_file_name_digi
        rm -f $current_file_name_recon
        echo " Sim files deleted"
    fi

done

echo ""
echo ""
echo "Series #$SERIES_NUMBER finished!"

