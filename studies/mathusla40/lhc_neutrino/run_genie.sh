
# Usage
if [ $# -ne 1 ]; then
	echo "Usage: $0  <RUN_NUMBER>" # Print help message if number of arguments is not equal to 1
	exit 1
fi
RUN_NUMBER=${1}
RUN_TIME=86400 # seconds
SIM_DIR=/home/tomren/geant_projects/mathusla-newsim
SIM_GENIE_DIR=/home/tomren/tmp
GENIE_INSTALL_DIR=/home/tomren/geant_projects/genie/Generator-R-3_04_00
export GENIE=$GENIE_INSTALL_DIR
export GENIE_DATA_DIR=/home/tomren/geant_projects/genie/genie_data


# ------------------------------------------------------------------------------------------
# NHLLHCνe
DATA_DIR=$SIM_GENIE_DIR/genie_output_NHLLHCνe
mkdir $DATA_DIR

#  The output filename is built as:
# [prefix].[run_number].[event_tree_format].[file_format]
$GENIE/bin/gevgen_flux \
-f $SIM_DIR/studies/mathusla40/lhc_neutrino/flux_data/combined_flux.root,12[NHLLHCνeAtEachEnergyGeV],-12[NHLLHCνeAtEachEnergyGeV] \
-T $RUN_TIME    -r $RUN_NUMBER    -o $DATA_DIR/gntp \
-g $GENIE_DATA_DIR/geometry/world.box.root \
--cross-sections $GENIE_DATA_DIR/genie_xsec/v3_04_00/NULL/G1802a00000-k250-e1000/data/gxspl-NUsmall.xml \
 --message-thresholds $GENIE/config/Messenger_laconic.xml


$GENIE/bin/gntpc  -f rootracker -i $DATA_DIR/gntp.$RUN_NUMBER.ghep.root -o $DATA_DIR/gntp.$RUN_NUMBER.gtrac.root
 --message-thresholds $GENIE/config/Messenger_laconic.xml


# ------------------------------------------------------------------------------------------
# NHLLHCνmu
DATA_DIR=$SIM_GENIE_DIR/genie_output_NHLLHCνμ
mkdir $DATA_DIR

#  The output filename is built as:
# [prefix].[run_number].[event_tree_format].[file_format]
$GENIE/bin/gevgen_flux \
-f $SIM_DIR/studies/mathusla40/lhc_neutrino/flux_data/combined_flux.root,14[NHLLHCνμAtEachEnergyGeV],-14[NHLLHCνμAtEachEnergyGeV] \
-T $RUN_TIME    -r $RUN_NUMBER    -o $DATA_DIR/gntp \
-g $GENIE_DATA_DIR/geometry/world.box.root \
--cross-sections $GENIE_DATA_DIR/genie_xsec/v3_04_00/NULL/G1802a00000-k250-e1000/data/gxspl-NUsmall.xml \
 --message-thresholds $GENIE/config/Messenger_laconic.xml


$GENIE/bin/gntpc  -f rootracker -i $DATA_DIR/gntp.$RUN_NUMBER.ghep.root -o $DATA_DIR/gntp.$RUN_NUMBER.gtrac.root
 --message-thresholds $GENIE/config/Messenger_laconic.xml 



# ------------------------------------------------------------------------------------------
# NHLLHCHardνμ
DATA_DIR=$SIM_GENIE_DIR/genie_output_NHLLHCHardνμ
mkdir $DATA_DIR

#  The output filename is built as:
# [prefix].[run_number].[event_tree_format].[file_format]
$GENIE/bin/gevgen_flux \
-f $SIM_DIR/studies/mathusla40/lhc_neutrino/flux_data/combined_flux.root,14[NHLLHCHardνμAtEachEnergyGeV],-14[NHLLHCHardantiνμAtEachEnergyGeV] \
-T $RUN_TIME    -r $RUN_NUMBER    -o $DATA_DIR/gntp  -E 0.3,2000 \
-g $GENIE_DATA_DIR/geometry/world.box.root \
--cross-sections $GENIE_DATA_DIR/genie_xsec/v3_04_00/NULL/G1802a00000-k250-e1000/data/gxspl-NUsmall.xml \
 --message-thresholds $GENIE/config/Messenger_laconic.xml


$GENIE/bin/gntpc  -f rootracker -i $DATA_DIR/gntp.$RUN_NUMBER.ghep.root -o $DATA_DIR/gntp.$RUN_NUMBER.gtrac.root
 --message-thresholds $GENIE/config/Messenger_laconic.xml 




# ------------------------------------------------------------------------------------------
# NHLLHCHardνμ
DATA_DIR=$SIM_GENIE_DIR/genie_output_NHLLHCBtoνμ
mkdir $DATA_DIR

#  The output filename is built as:
# [prefix].[run_number].[event_tree_format].[file_format]
$GENIE/bin/gevgen_flux \
-f $SIM_DIR/studies/mathusla40/lhc_neutrino/flux_data/combined_flux.root,14[NHLLHCBtoνμAtEachEnergyGeV],-14[NHLLHCBtoνμAtEachEnergyGeV] \
-T $RUN_TIME    -r $RUN_NUMBER    -o $DATA_DIR/gntp \
-g $GENIE_DATA_DIR/geometry/world.box.root \
--cross-sections $GENIE_DATA_DIR/genie_xsec/v3_04_00/NULL/G1802a00000-k250-e1000/data/gxspl-NUsmall.xml \
 --message-thresholds $GENIE/config/Messenger_laconic.xml


$GENIE/bin/gntpc  -f rootracker -i $DATA_DIR/gntp.$RUN_NUMBER.ghep.root -o $DATA_DIR/gntp.$RUN_NUMBER.gtrac.root
 --message-thresholds $GENIE/config/Messenger_laconic.xml  
 



 