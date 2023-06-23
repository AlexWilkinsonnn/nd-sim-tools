#!/bin/bash
################################################################################
# Script to run genie + edepsim
################################################################################
# Options

INPUTS_DIR="sim_inputs"

GEOMETRY="MPD_SPY_LAr.gdml"
TOPVOL="volArgonCubeActive"
EDEP_MAC="dune-nd.mac"

MODE="neutrino"
HORN="FHC"
RHC=""
FLUX="dk2nu"
FLUXOPT="--dk2nu"
FLUXDIR="/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/Flux/g4lbne/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017"
OFFAXIS=0
OADIR="0m"

NPOT=$1

################################################################################

RNDSEED=$RANDOM
RNDSEED=`echo "$RNDSEED" | cut -f 1 -d '.'`

NEVENTS="-e ${NPOT}"

ls

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
setup dk2nugenie   v01_06_01f -q debug:e15
setup genie_xsec   v2_12_10   -q DefaultPlusValenciaMEC
setup genie_phyopt v2_12_10   -q dkcharmtau
setup geant4 v4_10_3_p01b -q e15:prof

# edep-sim needs to know where a certain GEANT .cmake file is...
G4_cmake_file=`find ${GEANT4_FQ_DIR}/lib64 -name 'Geant4Config.cmake'`
export Geant4_DIR=`dirname $G4_cmake_file`

# edep-sim needs to have the GEANT bin directory in the path
export PATH=$PATH:$GEANT4_FQ_DIR/bin

cp ${INPUTS_DIR}/* .

chmod +x copy_dune_flux
./copy_dune_flux --top ${FLUXDIR} --output flux_files --flavor ${MODE} --maxmb=300 ${FLUXOPT}

echo "flux_files:"
ls flux_files

# Modify GNuMIFlux.xml to the specified off-axis position
sed -i "s/<beampos> ( 0.0, 0.05387, 6.66 )/<beampos> ( ${OFFAXIS}, 0.05387, 6.66 )/g" GNuMIFlux.xml

export GXMLPATH=${PWD}:${GXMLPATH}
export GNUMIXML="GNuMIFlux.xml"

# Run GENIE
echo "Running gevgen"
TIME_GENIE=`date +%s`
gevgen_fnal \
    -f flux_files/${FLUX}*,DUNEND \
    -g ${GEOMETRY} \
    -t ${TOPVOL} \
    -L cm -D g_cm3 \
    ${NEVENTS} \
    --seed ${RNDSEED} \
    -r ${RNDSEED} \
    -o ${MODE} \
    --message-thresholds Messenger_production.xml \
    --cross-sections ${GENIEXSECPATH}/gxspl-FNALsmall.xml \
    --event-record-print-level 0 \
    --event-generator-list Default+CCMEC
    #-m ${GEOMETRY}.${TOPVOL}.maxpl.xml \

# Convert the genie output to rootracker
echo "Running gntpc"
TIME_ROOTRACKER=`date +%s`
cp ${MODE}.${RNDSEED}.ghep.root input_file.ghep.root
gntpc -i input_file.ghep.root -f rootracker \
      --event-record-print-level 0 \
      --message-thresholds Messenger_production.xml

# edep-sim wants number of events, but we are doing POT so the files will be slightly different
# get the number of events from the GENIE files to pass it to edep-sim
NPER=$(echo "std::cout << gtree->GetEntries() << std::endl;" | \
       genie -l -b input_file.ghep.root 2>/dev/null | \
       tail -1)

setup edepsim v3_2_0 -q e20:prof

echo "Running edepsim"
edep-sim -C \
         -g $GEOMETRY \
         -o edep.${RNDSEED}.root \
         -e ${NPER} \
         ${EDEP_MAC}

