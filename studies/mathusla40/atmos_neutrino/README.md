GENIE

Modify gAtmoEvGen.cxx to use the histogram flux driver.

configure with --enable-atmo
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:\/opt/mathusla/release/lib

Let's first test the built-in atmospheric neutrino generator

``` bash
export GENIE_DATA_DIR=/home/tomren/geant_projects/genie/genie_data
export SIM_DIR=/home/tomren/geant_projects/mathusla-newsim
export FLUKA=$GENIE_DATA_DIR/flux_table_fluka

$GENIE/bin/gevgen_atmo -f FLUKA:$FLUKA/ok_nue02.dat[12],$FLUKA/ok_numu02.dat[14],$FLUKA/ok_anue02.dat[-12],$FLUKA/ok_anumu02.dat[-14] \
-g $GENIE_DATA_DIR/geometry/world.box.root \
-n 1 \
-r 3 -E 0.1,50  \
-R X:0.0,3.14159265,0.0 \
--flux-ray-generation-surface-distance 20 \
--flux-ray-generation-surface-radius 35  \
--cross-sections $GENIE_DATA_DIR/genie_xsec/v3_04_00/NULL/G1802a00000-k250-e1000/data/gxspl-NUsmall.xml \
 --message-thresholds $GENIE/config/Messenger_laconic.xml > log.txt
```


### Test new generator with bulk material

```bash
export GENIE_DATA_DIR=/home/tomren/geant_projects/genie/genie_data
export SIM_DIR=/home/tomren/geant_projects/mathusla-newsim

$GENIE/bin/gevgen_flux \
-f $SIM_DIR/studies/mathusla40/lhc_neutrino/flux_data/combined_flux.root,14[NHLLHCνμAtEachEnergyGeV],12[NHLLHCνeAtEachEnergyGeV] \
-g 1000080160[0.8879],1000010010[0.1121] \
-n 3000    -r 1    -E 0.1,50 \
--cross-sections $GENIE_DATA_DIR/genie_xsec/v3_04_00/NULL/G1802a00000-k250-e1000/data/gxspl-NUsmall.xml \
 --message-thresholds $GENIE/config/Messenger_laconic.xml
```


### Test with real geometry

**Generate 1 event** 

```bash
export GENIE_DATA_DIR=/home/tomren/geant_projects/genie/genie_data
export SIM_DIR=/home/tomren/geant_projects/mathusla-newsim

$GENIE/bin/gevgen_flux \
-f $SIM_DIR/studies/mathusla40/lhc_neutrino/flux_data/combined_flux.root,12[NHLLHCνeAtEachEnergyGeV],-12[NHLLHCνeAtEachEnergyGeV] \
-n 1    -r 1  \
-g $GENIE_DATA_DIR/geometry/world.box.root \
--cross-sections $GENIE_DATA_DIR/genie_xsec/v3_04_00/NULL/G1802a00000-k250-e1000/data/gxspl-NUsmall.xml \
 --message-thresholds $GENIE/config/Messenger_laconic.xml
```

**Generate 1 day of exposure**

Replace `-n 1` with `-T 86400`
It takes 3.5 minutes on my compuerter and took 2MB of disk space


```bash
$GENIE/bin/gevgen_flux \
-f $SIM_DIR/studies/mathusla40/lhc_neutrino/flux_data/combined_flux.root,12[NHLLHCνeAtEachEnergyGeV],-12[NHLLHCνeAtEachEnergyGeV] \
-T 86400    -r 1    \
-g $GENIE_DATA_DIR/geometry/world.box.root \
--cross-sections $GENIE_DATA_DIR/genie_xsec/v3_04_00/NULL/G1802a00000-k250-e1000/data/gxspl-NUsmall.xml \
 --message-thresholds $GENIE/config/Messenger_laconic.xml


 gntpc  -f rootracker -i gntp.1.ghep.root 
 --message-thresholds $GENIE/config/Messenger_laconic.xml
```


### Run everything together

run_genie.sh runs all four files. It takes 8 minutes and 6 MB for one day of exposure. 

100 days: 13 hours
1 year: 48 hours, 6GB
5 years: 240 hours, 30GB --> 100 day x 20










 
