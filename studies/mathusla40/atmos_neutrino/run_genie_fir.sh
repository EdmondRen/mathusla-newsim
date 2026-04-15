
# Usage
if [ $# -ne 1 ]; then
	echo "Usage: $0  <RUN_NUMBER>" # Print help message if number of arguments is not equal to 1
	exit 1
fi


export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:\/opt/mathusla/release/lib

RUN_NUMBER=${1}
RUN_TIME=$((86400*365)) # seconds, 365 days of exposure. Should take 2.5*120/60 ~ 5 hours
SIM_CODE_DIR=/project/6035200/tomren/jupyter/mathusla-newsim
SIM_GENIE_DIR=/project/6075887/MATHUSLA/simulation/genie/genie_output_2025_10
GENIE_INSTALL_DIR=/project/6075887/MATHUSLA/simulation/software/Genie_R-3_04_00
export GENIE=$GENIE_INSTALL_DIR
export GENIE_DATA_DIR=/project/6075887/MATHUSLA/simulation/genie/genie_input
FLUKA=$GENIE_DATA_DIR/flux_table_fluka


# ------------------------------------------------------------------------------------------
# NHLLHCνe
DATA_DIR=$SIM_GENIE_DIR/genie_output_atmos
mkdir $DATA_DIR

#  The output filename is built as:
# [prefix].[run_number].[event_tree_format].[file_format]
$GENIE/bin/gevgen_atmo \
-f FLUKA:$FLUKA/ok_nue02.dat[12],$FLUKA/ok_numu02.dat[14],$FLUKA/ok_anue02.dat[-12],$FLUKA/ok_anumu02.dat[-14]  \
-T $RUN_TIME    -r $RUN_NUMBER  --seed  $RUN_NUMBER -o $DATA_DIR/gntp \
-g $GENIE_DATA_DIR/geometry/world.box.root \
-E 0.1,50  \
-R X:0.0,3.14159265,0.0 \
--flux-ray-generation-surface-distance 35 \
--flux-ray-generation-surface-radius 35  \
--cross-sections $GENIE_DATA_DIR/genie_xsec/v3_04_00/NULL/G1802a00000-k250-e1000/data/gxspl-NUsmall.xml \
 --message-thresholds $GENIE/config/Messenger_laconic.xml


$GENIE/bin/gntpc  -f rootracker -i $DATA_DIR/gntp.$RUN_NUMBER.ghep.root -o $DATA_DIR/gntp.$RUN_NUMBER.gtrac.root
 --message-thresholds $GENIE/config/Messenger_laconic.xml



 