#!/bin/bash
################################################################################
# Script to produce edep-sim in LArBath geometry to process through
# rotation + translation throws to make ECC corrected ND-FD pairs.
# Also run through normal ND geometry to produce parameterised reco.
################################################################################
# Options

# GENIE_OUTPATH="/pnfs/dune/scratch/users/awilkins/lep_contained_pairs/genie"
# EDEP_OUTPATH="/pnfs/dune/scratch/users/awilkins/lep_contained_pairs/edep"
# CAF_OUTPATH="/pnfs/dune/scratch/users/awilkins/lep_contained_pairs/caf"

SAVE_GENIE=false
SAVE_EDEP=false # edep-sim output
SAVE_EDEP_MAKECAF=false # summarised edep-sim for parameterised reco in mackeCAF
SAVE_CAF=false # currently parameterised reco caf

INPUTS_DIR="sim_inputs_larbath"
ND_CAFMAKER_DIR="ND_CAFMaker"

GEOMETRY_LARBATH="LArBath_ndtopvol.gdml"
GEOMETRY_ND="MPD_SPY_LAr.gdml"
TOPVOL_ND="volArgonCubeActive"
EDEP_MAC="dune-nd.mac"

MODE="neutrino"
HORN="FHC"
RHC=""
FLUX="dk2nu"
FLUXOPT="--dk2nu"
FLUXDIR="/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/Flux/g4lbne/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017"
OFFAXIS=0
OADIR="0m"

INTERACTIVE=true # not running on grid

FIRST=$1
NPOT=$2 # 1e14 ~ 20-30 events

################################################################################

if [ "$INTERACTIVE" = true ]; then
  PROCESS=0 # or 1?
else
  echo "Running on $(hostname) at ${GLIDEIN_Site}. GLIDEIN_DUNESite = ${GLIDEIN_DUNESite}"
fi

RUNNO=$((${PROCESS}+${FIRST}))
RNDSEED=$(expr 1000000*${OFFAXIS}+${RUNNO}+1000000 | bc)
RNDSEED=`echo "$RNDSEED" | cut -f 1 -d '.'`

NEVENTS="-e ${NPOT}"

cp ${INPUTS_DIR}/* .
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
gevgen_fnal \
    -f flux_files/${FLUX}*,DUNEND \
    -g ${GEOMETRY_ND} \
    -t ${TOPVOL_ND} \
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
cp ${MODE}.${RNDSEED}.ghep.root input_file.ghep.root

# Convert the genie output to rootracker
echo "Running gntpc"
gntpc -i input_file.ghep.root -f rootracker \
      --event-record-print-level 0 \
      --message-thresholds Messenger_production.xml

# edep-sim wants number of events, but we are doing POT so the files will be slightly different
# get the number of events from the GENIE files to pass it to edep-sim
NPER=$(echo "std::cout << gtree->GetEntries() << std::endl;" | \
       genie -l -b input_file.ghep.root 2>/dev/null | \
       tail -1)

echo "NPER=${NPER}"

setup edepsim v3_2_0 -q e20:prof

echo "Running edepsim"
edep-sim -C \
         -g ${GEOMETRY_LARBATH} \
         -o edep_larbath.${RNDSEED}.root \
         -e ${NPER} \
         $EDEP_MAC
edep-sim -C \
         -g ${GEOMETRY_ND} \
         -o edep_nd.${RNDSEED}.root \
         -e ${NPER} \
         $EDEP_MAC

# Want to rollback environment to use old ND_CAFMaker scripts
# Unset all new env vars and then source the old env - probably overkill but it works
echo "Resetting env with env.sh"
unset $(comm -2 -3 <(\
        printenv | sed 's/=.*//' | sort) <(\
        sed -e 's/=.*//' -e 's/declare -x //' env.sh | sort))
source env.sh

source ${ND_CAFMAKER_DIR}/ndcaf_setup.sh
export GXMLPATH=${PWD}:${GXMLPATH}
export GNUMIXML="GNuMIFlux.xml"

echo "Running makeCAF dumpTree"
python dumpTree_nogeoeff.py --infile edep_nd.${RNDSEED}.root \
                            --outfile edep_dump_nd.${RNDSEED}.root \

echo "Running makeCAF"
cd $ND_CAFMAKER_DIR
./makeCAF --infile ../edep_dump_nd.${RNDSEED}.root \
          --gfile ../${MODE}.${RNDSEED}.ghep.root \
          --outfile ../${HORN}.${RNDSEED}.nd.CAF.root \
          --fhicl ../fhicl.fcl \
          --seed ${RNDSEED} \
          ${RHC} \
          --oa ${OFFAXIS}
cd ..


if [ "$INTERACTIVE" = false ]; then
  if [ "$SAVE_GENIE" = true ]; then
    ifdh cp ${MODE}.${RNDSEED}.ghep.root ${GENIE_OUTPATH}/${HORN}.${RNDSEED}.ghep.root
  fi
  if [ "$SAVE_EDEP" = true ]; then
    ifdh cp edep_larbath.${RNDSEED}.root ${EDEP_OUTPATH}/${HORN}.${RNDSEED}.edep_larbath.root
    ifdh cp edep_nd.${RNDSEED}.root ${EDEP_OUTPATH}/${HORN}.${RNDSEED}.edep_nd.root
  fi
  if [ "$SAVE_EDEP_MAKECAF" = true ]; then
    ifdh cp edep_dump.${RNDSEED}.root ${EDEP_OUTPATH}/${HORN}.${RNDSEED}.edep_dump.root
  fi
  if [ "$SAVE_CAF" = true ]; then
    ifdh cp ${HORN}.${RNDSEED}.nd.CAF.root ${CAF_OUTPATH}/${HORN}.${RNDSEED}.nd.CAF.root
  fi
fi

