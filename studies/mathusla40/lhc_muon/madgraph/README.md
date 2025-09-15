
This folder contains madgraph scripts for generating LHC backgrounds.

Card `card_wz_matched.dat` 

```
define lepall= mu+ mu-
define vorlepall = vl vl~ lepall
generate p p > lepall vorlepall @0
add process p p > lepall vorlepall j @1
output proc_sm_muprod_wz_matched
launch
1
0
set ickkw = 1
set xqcut = 20.0
set ptl   = 20.0
set etal  = 2.0
set drll  = 0
set iseed = 1
set use_syst = False
set pythia8_card HEPMCoutput:file hepmc
done
```