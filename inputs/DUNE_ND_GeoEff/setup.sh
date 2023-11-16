source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup dunetpc v09_41_00_02 -q e20:prof
setup root v6_22_08d -q e20:p392:prof
setup duneutil v09_65_01d00 -q e20:prof
#setup sam_web_client
setup cmake v3_24_1
setup gcc v6_4_0
setup eigen v3_3_5

# only need this line when doing interactively on dunegpvm
export PYTHONPATH=/dune/app/users/weishi/python3libs:$PYTHONPATH
setup geant4 v4_10_6_p01e -q e20:prof
setup edepsim v3_2_0 -q e20:prof

export PYTHONPATH=${PYTHONPATH}:${PWD}/lib/
