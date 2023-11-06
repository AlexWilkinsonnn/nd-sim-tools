#!/bin/bash
################################################################################
# Script to produce edep-sim from a single particle muon macro ready for passing
# to larnd-sim.
# Purpose is for contrastive learning particle classification dataset
################################################################################
# Options

EDEP_OUTPATH="/pnfs/dune/scratch/users/awilkins/gps_edep"

SAVE_EDEP=true # edep-sim output
SAVE_EDEP_H5=true # dumped to hdf5 for larnd-sim

INPUTS_DIR="sim_inputs_singleparticle"

GEOMETRY="MPD_SPY_LAr.gdml"
TOPVOL="volArgonCubeActive"
EDEP_MAC="multi_muon_shower.mac"

FIRST=$1
NEVENTS=$2

################################################################################

PROCESS=0
RUNNO=$((${PROCESS}+${FIRST}))
RNDSEED=$(expr ${RUNNO} | bc)
RNDSEED=`echo "$RNDSEED" | cut -f 1 -d '.'`

# echo "Running on $(hostname) at ${GLIDEIN_Site}. GLIDEIN_DUNESite = ${GLIDEIN_DUNESite}"
# ls

# cp -r ${INPUT_TAR_DIR_LOCAL}/${INPUTS_DIR} .

# Don't try over and over again to copy a file when it isn't going to work
export IFDH_CP_UNLINK_ON_ERROR=1
export IFDH_CP_MAXRETRIES=1
export IFDH_DEBUG=0

# Dump the current clean environment to to env.sh so it can be restored when needed
echo "Saving environment env.sh"
declare -px > env.sh

# Setting up for genie stuff
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup ifdhc v2_6_6
setup geant4 v4_10_3_p01b -q e15:prof

# edep-sim needs to know where a certain GEANT .cmake file is...
G4_cmake_file=`find ${GEANT4_FQ_DIR}/lib64 -name 'Geant4Config.cmake'`
export Geant4_DIR=`dirname $G4_cmake_file`

# edep-sim needs to have the GEANT bin directory in the path
export PATH=$PATH:$GEANT4_FQ_DIR/bin

cp ${INPUTS_DIR}/* .

setup edepsim v3_2_0 -q e20:prof

echo "Running edepsim"
edep-sim -C \
         -g $GEOMETRY \
         -o edep.muon_showers.${RNDSEED}.root \
         -e ${NEVENTS} \
         $EDEP_MAC

echo "Resetting env with env.sh"
cat env.sh
unset $(comm -2 -3 <(\
        printenv | sed 's/=.*//' | sort) <(\
        sed -e 's/=.*//' -e 's/declare -x //' env.sh | sort))
source env.sh

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup edepsim v3_2_0 -q e20:prof
setup ifdhc v2_6_6

python -m venv .venv_dumpTree
source .venv_dumpTree/bin/activate
pip install fire h5py numpy

echo "Running larndsim dumpTree"
python dumpTree_larndsimv0_3_4_multiprimaryvtx.py --input_file edep.muon_showers.${RNDSEED}.root \
                                                  --output_file edep.muon_showers.${RNDSEED}.h5

if [ "$SAVE_EDEP" = true ] ; then
  ifdh cp edep.muon_showers.${RNDSEED}.root ${EDEP_OUTPATH}/edep.muon_showers.${RNDSEED}.root
fi
if [ "$SAVE_EDEP_H5" = true ] ; then
  ifdh cp edep.muon_showers.${RNDSEED}.h5 ${EDEP_OUTPATH}/edep.muon_showers.${RNDSEED}.h5
fi

