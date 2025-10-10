
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

# Make a single run of tt,bb,EW

```bash
Extractor=muon_extract.py
OUTPUT_DIR=~/tmp/
pTcut=0

# tt, cross section 597 pb
TMP_DIR=$OUTPUT_DIR/mg5_tt
./1_run_card.sh card_tt_nomatch.dat $TMP_DIR 1
python3 ${Extractor} "${TMP_DIR}/Events/run_01/hepmc" "${OUTPUT_DIR}/bkg_muon_tt.joblib" $pTcut

# bb 3.632e+06 +- 1.018e+04 pb
TMP_DIR=$OUTPUT_DIR/mg5_bb
./1_run_card.sh card_bb_nomatch.dat $TMP_DIR 1
python3 ${Extractor} "${TMP_DIR}/Events/run_01/hepmc" "${OUTPUT_DIR}/bkg_muon_bb.joblib" $pTcut

# EW 7504 pb
TMP_DIR=$OUTPUT_DIR/mg5_ew
./1_run_card.sh card_wz_nomatch.dat $TMP_DIR 1
python3 ${Extractor} "${TMP_DIR}/Events/run_01/hepmc" "${OUTPUT_DIR}/bkg_muon_ew.joblib" $pTcut

```