#!/bin/bash
################################################################################
# Script to produce edep-sim ready for loading into FD SimEnergyDeposits +
# larnd-sim
# Only want events with lepton ending inside the LAr
################################################################################

################################################################################
# Options

GENIE_OUTPATH="/pnfs/dune/scratch/users/awilkins/lep_contained_pairs/genie"
EDEP_OUTPATH="/pnfs/dune/scratch/users/awilkins/lep_contained_pairs/edep"
CAF_OUTPATH="/pnfs/dune/scratch/users/awilkins/lep_contained_pairs/caf"

SAVE_GENIE=true
SAVE_EDEP=true # edep-sim output
SAVE_EDEP_FLAT=true # dumped to flat tree for easier reading into and FD data product
SAVE_EDEP_H5=true # dumped to hdf5 for larnd-sim
SAVE_EDEP_MAKECAF=false # summarised edep-sim for parameterised reco in mackeCAF
SAVE_CAF=true # currently parameterised reco caf

INPUTS_DIR="sim_inputs"
ND_CAFMAKER_DIR="ND_CAFMaker"

GEOMETRY="MPD_SPY_LArFullSensDet.gdml"
TOPVOL="volArgonCubeActive"
EDEP_MAC="dune-nd_smallsteps.mac"

MODE="neutrino"
HORN="FHC"
RHC=""
FLUX="dk2nu"
FLUXOPT="--dk2nu"
FLUXDIR="/pnfs/dune/persistent/users/ljf26/fluxfiles/g4lbne/v3r5p4/QGSP_BERT"
OFFAXIS=0
OADIR="0m"

FIRST=$1
NPOT=$2 # 1e14 ~ 20-30 events

################################################################################

RUNNO=$((${PROCESS}+${FIRST}))
RNDSEED=$(expr 1000000*${OFFAXIS}+${RUNNO}+1000000 | bc)
RNDSEED=`echo "$RNDSEED" | cut -f 1 -d '.'`

NEVENTS="-e ${NPOT}"

echo "Running on $(hostname) at ${GLIDEIN_Site}. GLIDEIN_DUNESite = ${GLIDEIN_DUNESite}"
ls

mv ${INPUT_TAR_DIR_LOCAL}/${INPUTS_DIR} .
mv ${INPUT_TAR_DIR_LOCAL}/${ND_CAFMAKER_DIR} .

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

# Get flux files to local node
# dk2nu files: /pnfs/dune/persistent/users/ljf26/fluxfiles/g4lbne/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017/neutrino/flux/dk2nu
# gsimple files: /pnfs/dune/persistent/users/dbrailsf/flux/nd/gsimple/v2_8_6d/OptimizedEngineeredNov2017/neutrino/
# NOTE ifdh ls breaks all the time
echo "Cheking ifdh ls is working"
echo "ifdh ls $FLUXDIR:"
ifdh ls $FLUXDIR

chmod +x copy_dune_ndtf_flux
./copy_dune_ndtf_flux --top ${FLUXDIR} --output local_flux_files --flavor ${MODE} --base OptimizedEngineeredNov2017 --maxmb=300 ${FLUXOPT}

echo "local_flux_files:"
ls local_flux_files

if [ "${FLUX}" = "dk2nu" ]; then
cd local_flux_files
for f in *.dk2nu.root
do
    mv "$f" "dk2nu_$f"
done
cd ..
fi

# Modify GNuMIFlux.xml to the specified off-axis position
sed -i "s/<beampos> ( 0.0, 0.05387, 6.66 )/<beampos> ( ${OFFAXIS}, 0.05387, 6.66 )/g" GNuMIFlux.xml

export GXMLPATH=${PWD}:${GXMLPATH}
export GNUMIXML="GNuMIFlux.xml"

# Run GENIE
echo "Running gevgen"
TIME_GENIE=`date +%s`
gevgen_fnal \
    -f local_flux_files/${FLUX}*,DUNEND \
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
NPER=$(echo "std::cout << gtree->GetEntries() << std::endl;" | genie -l -b input_file.ghep.root 2>/dev/null  | tail -1)

setup edepsim v3_2_0 -q e20:prof

echo "Running edepsim"
edep-sim -C \
         -g $GEOMETRY \
         -o edep.${RNDSEED}.root \
         -e ${NPER} \
         $EDEP_MAC

# Want to rollback environment to use old ND_CAFMaker scripts
# Unset all new env vars and then source the old env - probably overkill but it works
echo "Resetting env with env.sh"
unset $(comm -2 -3 <(printenv | sed 's/=.*//' | sort) <(sed -e 's/=.*//' -e 's/declare -x //' env.sh | sort))
source env.sh

source ${ND_CAFMAKER_DIR}/ndcaf_setup.sh
export GXMLPATH=${PWD}:${GXMLPATH}
export GNUMIXML="GNuMIFlux.xml"

echo "Running makeCAF dumpTree"
python dumpTree_nogeoeff.py --infile edep.${RNDSEED}.root \
                            --outfile edep_dump.${RNDSEED}.root \
                            --muoncontained_ids_file muoncontained_eventids.${RNDSEED}.txt

echo "Running makeCAF"
cd $ND_CAFMAKER_DIR
./makeCAF --infile ../edep_dump.${RNDSEED}.root \
          --gfile ../${MODE}.${RNDSEED}.ghep.root \
          --outfile ../${HORN}.${RNDSEED}.CAF.root \
          --fhicl ../fhicl.fcl \
          --seed ${RNDSEED} \
          ${RHC} \
          --oa ${OFFAXIS}
cd ..

echo "Resetting env with env.sh"
cat env.sh
unset $(comm -2 -3 <(printenv | sed 's/=.*//' | sort) <(sed -e 's/=.*//' -e 's/declare -x //' env.sh | sort))
source env.sh

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup edepsim v3_2_0 -q e20:prof
setup ifdhc v2_6_6

python -m venv .venv_dumpTree
source .venv_dumpTree/bin/activate
pip install fire h5py numpy

echo "Running larndsim dumpTree"
python dumpTree_larndsimv0_3_4.py --input_file edep.${RNDSEED}.root \
                                  --output_file edep.${RNDSEED}.h5 \
                                  --muoncontained_ids_file muoncontained_eventids.${RNDSEED}.txt

echo "Running FD dumpTree"
python dump_edep_to_flattree.py edep.${RNDSEED}.root \
                                edep_flat.${RNDSEED}.root \
                                --muoncontained_ids_file muoncontained_eventids.${RNDSEED}.txt

if [ "$SAVE_GENIE" = true ] ; then
  ifdh cp ${MODE}.${RNDSEED}.ghep.root ${GENIE_OUTPATH}/${HORN}.${RNDSEED}.ghep.root
fi
if [ "$SAVE_EDEP" = true ] ; then
  ifdh cp edep.${RNDSEED}.root ${EDEP_OUTPATH}/${HORN}.${RNDSEED}.edep.root
fi
if [ "$SAVE_EDEP_FLAT" = true ] ; then
  ifdh cp edep_flat.${RNDSEED}.root ${EDEP_OUTPATH}/${HORN}.${RNDSEED}.edep_flat.root
fi
if [ "$SAVE_EDEP_H5" = true ] ; then
  ifdh cp edep.${RNDSEED}.h5 ${EDEP_OUTPATH}/${HORN}.${RNDSEED}.edep.h5
fi
if [ "$SAVE_MAKECAF" = true ] ; then
  ifdh cp ${HORN}.${RNDSEED}.CAF.root ${CAF_OUTPATH}/${HORN}.${RNDSEED}.CAF.root
fi

