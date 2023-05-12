#!/bin/bash
################################################################################
# Script to dump the file numbers for running over with a grid job
# ... could use samweb but this is somehow easier
# Expects clean_output_pnfs.sh has been run so that all files have pairs across
# the genie/edep/caf directories as required. This way can just use one dir
# to see all the numbers
################################################################################
# Options

TARGET_DIR="/pnfs/dune/scratch/users/awilkins/lep_contained_pairs/genie"

################################################################################

ls ${TARGET_DIR} | sed "s/FHC\.\([0-9]*\)\..*/\1/" > file_nums.txt

