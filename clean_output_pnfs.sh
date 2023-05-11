#!/bin/bash
################################################################################
# Script to clean output of produce_edep-paramreco_muoncontained.sh.
# Deletes incomplete sets genie/edep/caf files
################################################################################
# Options

GENIE_PATH="/pnfs/dune/scratch/users/awilkins/lep_contained_pairs/genie"
EDEP_PATH="/pnfs/dune/scratch/users/awilkins/lep_contained_pairs/edep"
CAF_PATH="/pnfs/dune/scratch/users/awilkins/lep_contained_pairs/caf"

REQUIRE_GENIE=true # .ghep.root
REQUIRE_EDEP=true # .edep.root
REQUIRE_EDEP_FLAT=true # .edep_flat.root
REQUIRE_EDEP_H5=true # .edep.h5
REQUIRE_EDEP_MAKECAF=false # .edep_dump.root
REQUIRE_CAF=true # .CAF.root

DRY_RUN=false # don't actually delete files

################################################################################

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup ifdhc v2_6_8

num_arr=()
if [ "$REQUIRE_GENIE" = true ]; then
  for filepath in ${GENIE_PATH}/*; do
    file=${filepath##*/}
    num=$(echo $file | sed -e "s/FHC.\(.*\).ghep.root/\1/")
    num_arr+=($num)
  done
fi
if [ "$REQUIRE_EDEP" = true ] || \
   [ "$REQUIRE_EDEP_FLAT" = true ] || \
   [ "$REQUIRE_EDEP_H5" = true ] || 
   [ "$REQUIRE_EDEP_MAKECAF" = true ]; then
  for filepath in ${EDEP_PATH}/*; do
    file=${filepath##*/}
    if [ "$REQUIRE_EDEP" = true ] && [[ "$file" == *".edep.root" ]]; then
      num=$(echo $file | sed -e "s/FHC.\(.*\).edep.root/\1/")
      num_arr+=($num)
    elif [ "$REQUIRE_EDEP_FLAT" = true ] && [[ "$file" == *".edep_flat.root" ]]; then
      num=$(echo $file | sed -e "s/FHC.\(.*\).edep_flat.root/\1/")
      num_arr+=($num)
    elif [ "$REQUIRE_EDEP_H5" = true ] && [[ "$file" == *".edep.h5" ]]; then
      num=$(echo $file | sed -e "s/FHC.\(.*\).edep.h5/\1/")
      num_arr+=($num)
    elif [ "$REQUIRE_EDEP_MAKECAF" = true ] && [[ "$file" == *".edep_dump.root" ]]; then
      num=$(echo $file | sed -e "s/FHC.\(.*\).edep_dump.root/\1/")
      num_arr+=($num)
    fi
  done
fi
if [ "$REQUIRE_CAF" = true ]; then
  for filepath in ${CAF_PATH}/*; do
    file=${filepath##*/}
    num=$(echo $file | sed -e "s/FHC.\(.*\).CAF.root/\1/")
    num_arr+=($num)
  done
fi

uniq_num_arr=($(for num in "${num_arr[@]}"; do echo "${num}"; done | sort -u))

# for num in "${uniq_num_arr[@]}"; do
#   echo $num
# done

for num in "${uniq_num_arr[@]}"; do
  genie_filepath=${GENIE_PATH}/FHC.${num}.ghep.root
  edep_filepath=${EDEP_PATH}/FHC.${num}.edep.root
  edep_flat_filepath=${EDEP_PATH}/FHC.${num}.edep_flat.root
  edep_h5_filepath=${EDEP_PATH}/FHC.${num}.edep.h5
  edep_makecaf_filepath=${EDEP_PATH}/FHC.${num}.edep_dump.root
  caf_filepath=${CAF_PATH}/FHC.${num}.CAF.root

  if [ "$REQUIRE_GENIE" = true ] && [ ! -f $genie_filepath ]; then
    :;
  elif [ "$REQUIRE_EDEP" = true ] && [ ! -f $edep_filepath ]; then
    :;
  elif [ "$REQUIRE_EDEP_FLAT" = true ] && [ ! -f $edep_flat_filepath ]; then
    :;
  elif [ "$REQUIRE_EDEP_H5" = true ] && [ ! -f $edep_h5_filepath ]; then
    :;
  elif [ "$REQUIRE_EDEP_MAKECAF" = true ] && [ ! -f $edep_makecaf_filepath ]; then
    :;
  elif [ "$REQUIRE_CAF" = true ] && [ ! -f $caf_filepath ]; then
    :;
  else
    continue
  fi 

  if [ -f $genie_filepath ]; then
    echo "Deleting $genie_filepath"
    if [ "${DRY_RUN}" = false ]; then
      ifdh rm $genie_filepath
    fi
  fi
  if [ -f $edep_filepath ]; then
    echo "Deleting $edep_filepath"
    if [ "${DRY_RUN}" = false ]; then
      ifdh rm $edep_filepath
    fi
  fi
  if [ -f $edep_flat_filepath ]; then
    echo "Deleting $edep_flat_filepath"
    if [ "${DRY_RUN}" = false ]; then
      ifdh rm $edep_flat_filepath
    fi
  fi
  if [ -f $edep_h5_filepath ]; then
    echo "Deleting $edep_h5_filepath"
    if [ "${DRY_RUN}" = false ]; then
      ifdh rm $edep_h5_filepath
    fi
  fi
  if [ -f $edep_makecaf_filepath ]; then
    echo "Deleting $edep_makecaf_filepath"
    if [ "${DRY_RUN}" = false ]; then
      ifdh rm $edep_makecaf_filepath
    fi
  fi
  if [ -f $caf_filepath ]; then
    echo "Deleting $caf_filepath"
    if [ "${DRY_RUN}" = false ]; then
      ifdh rm $caf_filepath
    fi
  fi
done

