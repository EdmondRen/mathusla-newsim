basedir="/project/6049244/data/MATHUSLA/simulation_v2"

## Special settings for the Dry run
# 1. -c n : Do no delete simulation and digi files,
# 2. -o n : Do no delete simulation and digi files,

##===============================================================================
## Cosmic protons
proton_events_per_run=3000000
proton_runs_per_job=29
proton_n_jobs=1200
proton_job_time_hours=109
proton_run_name="cosmic_p"
proton_data_directory=$basedir/cosmic/$proton_run_name
mkdir -p $proton_data_directory

bash _submit_series.sh \
    -f $proton_run_name/_start_series_of_run.sh \
    -m $proton_data_directory \
    -e $proton_events_per_run \
    -r $proton_runs_per_job \
    -j $proton_n_jobs \
    -t $proton_job_time_hours \
    -n $proton_run_name \
    -c y \
    -s True

##===============================================================================
## Cosmic neutrons
neutron_events_per_run=8000000
neutron_runs_per_job=36
neutron_n_jobs=2800
neutron_job_time_hours=109
neutron_run_name="cosmic_n"
neutron_data_directory=$basedir/cosmic/$neutron_run_name
mkdir -p $neutron_data_directory

bash _submit_series.sh \
    -f $neutron_run_name/_start_series_of_run.sh \
    -m $neutron_data_directory \
    -e $neutron_events_per_run \
    -r $neutron_runs_per_job \
    -j $neutron_n_jobs \
    -t $neutron_job_time_hours \
    -n $neutron_run_name \
    -c y \
    -s True
