./_start_single_run.sh ~/geant_projects/mathusla-newsim data 10000 1 proton 10000

./upward_n_and_p/_start_single_run.sh /home/tomren/geant_projects/mathusla-newsim data 1000 3 neutron 10000

Total: 10k Events  
With vertex: 460
Upward+inside: 100
Time: 1 minutes on my PC; 3 minutes on cedar
Storage: 1.4MB (final)

For submitting:
Let's do 40x of the dry run:

400k event per job (3 hours)
10 jobs each
Loop (n,p)x(2GeV, 5GeV, 10GeV)


