basedir="/project/6049244/data/MATHUSLA/simulation_v2"

##===============================================================================
## Cosmic protons
proton_events_per_run=3000000
proton_runs_per_job=4
proton_n_jobs=3
proton_job_time_hours=18
proton_run_name="cosmic_p"
proton_data_directory=$basedir/test/$proton_run_name
mkdir -p $proton_data_directory

bash _submit_series.sh \
    -f $proton_run_name/_start_series_of_run.sh \
    -m $proton_data_directory \
    -e $proton_events_per_run \
    -r $proton_runs_per_job \
    -j $proton_n_jobs \
    -t $proton_job_time_hours \
    -n $proton_run_name \
    -s True    



##===============================================================================
## Cosmic neutrons
neutron_events_per_run=8000000
neutron_runs_per_job=4
neutron_n_jobs=3
neutron_job_time_hours=12
neutron_run_name="cosmic_n"
neutron_data_directory=$basedir/test/$neutron_run_name
mkdir -p $neutron_data_directory

bash _submit_series.sh \
    -f $neutron_run_name/_start_series_of_run.sh \
    -m $neutron_data_directory \
    -e $neutron_events_per_run \
    -r $neutron_runs_per_job \
    -j $neutron_n_jobs \
    -t $neutron_job_time_hours \
    -n $neutron_run_name \
    -s True        