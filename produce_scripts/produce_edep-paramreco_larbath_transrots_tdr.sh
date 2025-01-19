#!/bin/bash
################################################################################
# Script to produce edep-sim in LArBath geometry to process through
# rotation + translation throws to make ECC corrected ND-FD pairs.
# Using a configuration as close to the TDR production as I can.
# Also run through normal ND geometry to produce parameterised reco.
# NDCAFMaker dumper script heavily modified to work with LArBath
# events in the selection throws format.
################################################################################
# Options

GENIE_OUTPATH="/pnfs/dune/scratch/users/colweber/larbath_ndfd_pairs/tdr_sample/RHC/genie"
EDEP_OUTPATH="/pnfs/dune/scratch/users/colweber/larbath_ndfd_pairs/tdr_sample/RHC/edep"
CAF_OUTPATH="/pnfs/dune/scratch/users/colweber/larbath_ndfd_pairs/tdr_sample/RHC/caf"
NDFD_ROOT_OUTPUT="/pnfs/dune/scratch/users/colweber/larbath_ndfd_pairs/tdr_sample/RHC/pair_root"
PAIR_H5_OUTPUT="/pnfs/dune/scratch/users/colweber/larbath_ndfd_pairs/tdr_sample/FHC/pair_allinfo_h5"

SAVE_GENIE=false
SAVE_GTRAC=false
SAVE_EDEP=false # edep-sim output
SAVE_EDEP_MAKECAF=false # summarised edep-sim for parameterised reco in mackeCAF
SAVE_CAF=false # currently parameterised reco caf
SAVE_NDFD_ROOT=false # nd and fd depo data after translation + rotation throws to select
SAVE_PAIR_H5=true # nd-fd paired data file

# These dirs need to be in the job tarball
INPUTS_DIR="sim_inputs_larbath_selected_ndfd_pairs"
ND_CAFMAKER_DIR="ND_CAFMaker"
TRANSROTS_DIR="DUNE_ND_GeoEff"

GEOMETRY_ND="MPD_SPY_LAr.gdml"

TOPVOL_ND="volArgonCubeActive"
EDEP_MAC="dune-nd.mac"
EDEPSIM_ANA_CFG="UserConfig_tdr_nofdthrows.py"

MODE="antineutrino" # "neutrino" or "antineutrino" (controls horn current)
HORN="RHC" # "FHC" or "RHC" (controls naming of files)
RHC="--rhc" # "" or "--rhc" (controls ND parametrized reconstruction algorithm)
FLUX="dk2nu"
FLUXOPT="--dk2nu"
# FLUX="gsimple"
# FLUXOPT="--gsimple"
FLUXDIR="/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/Flux/g4lbne/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017"
OFFAXIS=0
OADIR="0m"

INTERACTIVE=false # not running on grid

FIRST=$1
NPOT=$2 # 1e14 ~ 20-30 events

################################################################################

if [ "$INTERACTIVE" = true ]; then
  PROCESS=0 # or 1?
else
  echo "Running on $(hostname) at ${GLIDEIN_Site}. GLIDEIN_DUNESite = ${GLIDEIN_DUNESite}"
  cp -r ${INPUT_TAR_DIR_LOCAL}/${INPUTS_DIR} .
  cp -r ${INPUT_TAR_DIR_LOCAL}/${ND_CAFMAKER_DIR} .
  cp -r ${INPUT_TAR_DIR_LOCAL}/${TRANSROTS_DIR} .
fi

RUNNO=$((${PROCESS}+${FIRST}))
RNDSEED=$(expr 1000000*${OFFAXIS}+${RUNNO}+1000000 | bc)
RNDSEED=`echo "$RNDSEED" | cut -f 1 -d '.'`

NEVENTS="-e ${NPOT}"

cp ${INPUTS_DIR}/* .
echo "LS-ing inputs after initial copy."
ls -rt

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

# I have no idea what changed to make the above cause gevgen to not find the GNuMIFlux.xml config.
# Note if you do export GXMLPATH=${PWD} instead it will find the correct xml but not run because
# it cannot find any splines
# NOTE Turns out I didnt need this, just having a fresh shell solved my problems (I didnt think
# this script could inherit anything important from my env but idk im dumb)
# export GXMLPATH=${PWD}:${GXMLPATH}
# export GNUMIXML="GNuMIFlux.xml"
export GNUMIFLUXXML="${PWD}/GNuMIFlux.xml"
export GDK2NUFLUXXML="${PWD}/GNuMIFlux.xml"
echo "LS-ing inputs pre-gevgen"
ls -rt
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
echo "LS-ing inputs post-gevgen, pre gntpc"
ls -rt
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

setup edepsim v3_0_0 -q e19:prof
echo "LS-ing inputs post gntpc, pre edep-sim on LAr world"
ls -rt
echo "Running edepsim"

# When running edep-sim in LAr world, we want the hits to automatically 
# break at the points where they would be crossing a material boundary in 
# the real world. This is so the individual hits don't cross boundaries when 
# examined in the real world, preventing us from having to interpolate the 
# material they passed through.
# To accomplish this, load up the ND geometry file and change all the 
# material references to "LAr". Use this file for edep-sim instead of the 
# infinite LAr bath.
# We also want the entire detective to be active, so we add the string 
# <auxiliary auxtype="SensDet" auxvalue="SimEnergyDeposit"/> where 
# appropriate
ALL_LAr_GEOMETRY_ND="ALL_LAr_$GEOMETRY_ND"
if [ ! -f $ALL_LAr_GEOMETRY_ND ];
then
	bash activate_all_lar_geometry_nd.sh -i $GEOMETRY_ND -o $ALL_LAr_GEOMETRY_ND
fi

edep-sim -C \
         -g ${ALL_LAr_GEOMETRY_ND} \
         -o edep_larbath.${RNDSEED}.root \
         -e ${NPER} \
         $EDEP_MAC

echo "LS-ing inputs post edep-sim on LAr world, pre edep-sim on ND"
ls -rt
# The ND hall is not best-represented by an infinite LAr bath, so the 
# simulation is improved by re-simulating just the muon in the ND hall. We 
# run edep-sim just on the muon in the ND hall so that the hadronic energy 
# deposits remain identical.
edep-sim -C \
         -g ${GEOMETRY_ND} \
         -o edep_ND.${RNDSEED}.root \
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
export GNUMIFLUXXML="${PWD}/GNuMIFlux.xml"
export GDK2NUFLUXXML="${PWD}/GNuMIFlux.xml"
echo "LS-ing inputs post edep-sim on ND, pre dumpTree on LAr"
ls -rt
echo "Running makeCAF dumpTree"
python dumpTree_tdr_nogeoeff_larbath.py --infile_edepsim edep_larbath.${RNDSEED}.root \
                                        --edepsim_geometry edep_ND.${RNDSEED}.root \
                                        --outfile edep_dump_larbath_nd.${RNDSEED}.root

echo "LS-ing inputs post dumpTree on LAr, pre dumpTree on ND"
ls -rt
python dumpTree_tdr_nogeoeff_larbath.py --infile_edepsim edep_ND.${RNDSEED}.root \
                                        --edepsim_geometry edep_ND.${RNDSEED}.root \
                                        --outfile edep_dump_ND_nd.${RNDSEED}.root

echo "LS-ing inputs post dumpTree on ND, pre makeCAF"
ls -rt
echo "Running makeCAF"
cd $ND_CAFMAKER_DIR
./makeCAF_resim-muon --infile ../edep_dump_larbath_nd.${RNDSEED}.root \
					 --infile_resim ../edep_dump_ND_nd.${RNDSEED}.root \
          			 --gfile ../${MODE}.${RNDSEED}.ghep.root \
          			 --outfile ../${HORN}.${RNDSEED}.nd.CAF.root \
          			 --fhicl ../fhicl.fcl \
          			 --seed ${RNDSEED} \
          			 ${RHC} \
          			 --oa ${OFFAXIS}
cd ..
echo "LS-ing inputs after makeCAF"
ls -rt
# Reset env again for GeoEff rotations+translation
echo "Resetting env with env.sh"
unset $(comm -2 -3 <(printenv | sed 's/=.*//' | sort) <(sed -e 's/=.*//' -e 's/declare -x //' env.sh | sort))
source env.sh

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup dunetpc v09_41_00_02 -q e20:prof
setup ifdhc v2_6_6
setup root v6_22_08d -q e20:p392:prof
setup duneutil v09_65_01d00 -q e20:prof
setup cmake v3_24_1
setup gcc v6_4_0
setup eigen v3_3_5
setup geant4 v4_10_6_p01e -q e20:prof
setup edepsim v3_2_0 -q e20:prof

python -m venv .venv_3.9.2_ndfd_pairs
source .venv_3.9.2_ndfd_pairs/bin/activate
# scipy 1.10 is the latest version compatible with edep-sim's numpy version 1.20.1 which the
# the venv cannot overwrite
pip install torch scipy==1.10 h5py fire
if [ "$INTERACTIVE" = true ]; then
  export PYTHONPATH=${PWD}/${TRANSROTS_DIR}/lib:${PYTHONPATH}
  export LD_LIBRARY_PATH=${PWD}/${TRANSROTS_DIR}/lib:${LD_LIBRARY_PATH}
else
  export PYTHONPATH=${_CONDOR_JOB_IWD}/${TRANSROTS_DIR}/lib:${PYTHONPATH}
  export LD_LIBRARY_PATH=${_CONDOR_JOB_IWD}/${TRANSROTS_DIR}/lib:${LD_LIBRARY_PATH}
fi

echo "Running translation + rotation throws to get selected nd-fd pairs"
mkdir n2fd_outputs
cd ${TRANSROTS_DIR}/app
python Edepsim_ana.py --config ../../${EDEPSIM_ANA_CFG} \
                      --out_dir ../../n2fd_outputs \
					  --caf_file ../../${HORN}.${RNDSEED}.nd.CAF.root \
                      ../../edep_larbath.${RNDSEED}.root # 1> /dev/null 2/ /dev/null
cd ../../

echo "Running nd-fd pair maker"
echo "LS-ing inputs post throws"
ls -lrth
echo "LS-ing n2fd_outputs/* post thtows"
ls -lrth n2fd_outputs/*
python dumpTree_larndsimv0_3_4_transrots-paramreco.py --param_reco_file ${HORN}.${RNDSEED}.nd.CAF.root \
                                                      n2fd_outputs/root_out/n2fd_paired_out.root \
                                                      ${HORN}.${RNDSEED}.ndfd_preco_pairs.h5

echo "Copying files to dCache..."
if [ "$SAVE_GENIE" = true ]; then
  ifdh cp ${MODE}.${RNDSEED}.ghep.root ${GENIE_OUTPATH}/${HORN}.${RNDSEED}.ghep.root
fi
if [ "$SAVE_GTRAC" = true ]; then
  ifdh cp ${MODE}.${RNDSEED}.gtrac.root ${GENIE_OUTPATH}/${HORN}.${RNDSEED}.gtrac.root
fi
if [ "$SAVE_EDEP" = true ]; then
  ifdh cp edep_larbath.${RNDSEED}.root ${EDEP_OUTPATH}/${HORN}.${RNDSEED}.edep_larbath.root
  ifdh cp edep_ND.${RNDSEED}.root ${EDEP_OUTPATH}/${HORN}.${RNDSEED}.edep_ND.root
fi
if [ "$SAVE_EDEP_MAKECAF" = true ]; then
  ifdh cp edep_dump_larbath_nd.${RNDSEED}.root ${EDEP_OUTPATH}/${HORN}.${RNDSEED}.edep_dump_larbath_nd.root
  ifdh cp edep_dump_ND_nd.${RNDSEED}.root ${EDEP_OUTPATH}/${HORN}.${RNDSEED}.edep_dump_ND_nd.root
fi
if [ "$SAVE_CAF" = true ]; then
  ifdh cp ${HORN}.${RNDSEED}.nd.CAF.root ${CAF_OUTPATH}/${HORN}.${RNDSEED}.nd.CAF.root
fi
if [ "$SAVE_NDFD_ROOT" = true ]; then
  ifdh cp n2fd_outputs/root_out/n2fd_paired_out.root ${NDFD_ROOT_OUTPUT}/${HORN}.${RNDSEED}.n2fd_paired_out.root
fi
if [ "$SAVE_PAIR_H5" = true ]; then
  ifdh cp ${HORN}.${RNDSEED}.ndfd_preco_pairs.h5 ${PAIR_H5_OUTPUT}/${HORN}.${RNDSEED}.ndfd_preco_pairs.h5
fi

