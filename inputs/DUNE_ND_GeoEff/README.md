# Set up

[First time only]
```
cd /dune/app/users/<your_username>
mkdir testn2fd
cd testn2fd

git clone --recurse-submodules -b N2FD https://github.com/weishi10141993/DUNE_ND_GeoEff.git      # Get geoEff library
# Note for git version (git --version) before 2.13, use: git clone --recursive -b N2FD https://github.com/weishi10141993/DUNE_ND_GeoEff.git
cd DUNE_ND_GeoEff
source setup.sh                                                                                  # Necessary setups for build
cmake -DPYTHON_EXECUTABLE:FILEPATH=`which python` .
make -j geoEff                                                                                   # Build geoEff
make -j pyGeoEff                                                                                 # Build pygeoEff
```

# Analyze edepsim data

```
source setup.sh

cd app
nohup python3 Edepsim_ana.py /dune/data/users/awilkins/extrapolation/edep.LArBath.NDGenieGen.root >& output.log &

# The first time you may need to install a few packages via pip install, depending on what it complains when you run, e.g.:
pip install --target=/dune/app/users/weishi/python3libs torch --upgrade
pip install --target=/dune/app/users/weishi/python3libs scipy --upgrade
```
