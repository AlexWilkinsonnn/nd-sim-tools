#!/bin/bash

VENV_DIR="${1:-.}"

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup edepsim v3_2_0 -q e20:prof

if [ -d "${VENV_DIR}/.venv_dumpTree" ]; then
  echo "Using venv at ${VENV_DIR}/.venv_dumpTree"
  source .venv_dumpTree/bin/activate
else
  echo "Creating venv at ${VENV_DIR}/.venv_dumpTree"
  python -m venv ${VENV_DIR}/.venv_dumpTree
  source .venv_dumpTree/bin/activate
  pip install fire h5py numpy
fi

