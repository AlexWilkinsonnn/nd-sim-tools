#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!                           USER Configure Below                         !
#!      Some parameters are used in geometric efficiency correction       !
#!      Some are used in near to far event translation                    !
#!      Do not delete even if it's not familiar to you                    !
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
from array import array

# Convert from EDepSim default units (mm, ns)
edep2cm = 0.1   # convert to cm
edep2us = 0.001 # convert to microseconds
# Convert GENIE to common units
gev2mev  = 1000 # convert GeV to MeV
# LArBath length each side
dimension = 40000 #cm
# Beamline point downward: clockwise rotate around ND det coordinate x axis
beamLineRotation = -0.101 #rad
# Coordinate transformation, units: meters
beamRefDetCoord = [0.0, 0.05387, 6.6] # beam center crossing ND, (0, 0, 0) is ND detector origin
detRefBeamCoord = [0.0, 0.0,    574.] # (0, 0, 0) is beam origin

verbose = False # default to False, true for debugging, a lot of printouts
myfileVerbose = False
throwfileVerbose = False
hadronhitVerbose = False
ntupleVerbose = False

N_throws = 64*64
N_THROWS_FD = 1
max_nd_throws = 20
max_fd_throws = 20
###
# Turning off selection throws by setting threshold s.t. the throw always passes
# This will do one throw at ND, ECC the event and do one throw at FD to get a vertex.
# For the tdr production the ND throw is not used for anything anyway.
###
# Cut on mu nn probability 
nd_contained_mu_cut = 0.0
nd_tracker_mu_cut = 0.0
# Hadronic veto selection crtieria
nd_veto_region_size = 1 # cm
nd_veto_threshold_energy = 1e9 # MeV
fd_veto_region_size = 1 # cm
fd_veto_threshold_energy = 1e9 # MeV

# Active volume for FD unit [cm]
# These correspond to the prism vetex cuts used for FD selection
FDActiveVol_min = [-370., -600.,    0.]
FDActiveVol_max = [ 370.,  600., 1400.]

# Active volume for ND unit [cm]
NDActiveVol_min = [-350., -150.,   0.]
NDActiveVol_max = [ 350.,  150., 500.]
# this offset [cm] is only for ND MC when On-Axis, use 0 for FD MC
# !!! however this offset is removed in newest duneana dumptree by jeremy!!!
# https://github.com/AlexWilkinsonnn/ND_CAFMaker/commit/b4b1efc547c59c479a2aefa32e8d6a1452c45028
NDLAr_OnAxis_offset = [0.,  5.5, 411.0]

# Fiducial volume for ND unit [cm]
ND_FV_min = [-300., -100., 50.]
ND_FV_max = [ 300.,  100., 350.]

# Fiducial volume for FD ( minus the same amount from FDActiveVol) unit [cm]
FD_FV_min = [-310., -550.,   50.]
FD_FV_max = [ 310.,  550., 1244.]

random_ND_off_axis_pos   = False # Set to true will only use a random ND off axis position per event in runGeoEffFDEvtSim
ND_off_axis_pos_stepsize = 2.5   # unit meters
# Average neutrino decay position in beam coordinates as a function of off axis x (from Luke)
# This is used to interpolate the decay position unit meters
meanPDPZ      = array('f', [57.7352, 59.1433, 60.8336, 62.4821, 65.005, 67.5822, 69.8348, 73.0865, 76.6664, 81.1443, 85.6266, 90.346, 93.362, 93.6072, 93.362,  90.346, 85.6266, 81.1443, 76.6664, 73.0865, 69.8348, 67.5822, 65.005, 62.4821, 60.8336, 59.1433, 57.7352]) # unit: meters
OffAxisPoints = array('f', [-30.5,   -28,     -25.5,   -23,     -20.5,  -18,     -15.5,   -13,     -10.5,   -8,      -5.5,    -3,     -0.5,   0,       0.5,     3,      5.5,     8,       10.5,    13,      15.5,    18,      20.5,   23,      25.5,    28,      30.5])
NDLarPos_new1 = [0., -28.] # unit meter
NDLarPos_new2 = [-1.75, -25.75]
ND_Lar_dtctr_pos_new_stepsize = 4

random_vtx_vx       = False # Set to true will only use a random vtx x per event in runGeoEffFDEvtSim
ND_local_x_stepsize = 48. # unit cm, must be a positive number below 200
ND_local_x_min      = -312. # unit cm,
ND_local_x_max      = 312.  # unit cm,

def IsFromPrimaryLep(trkid, parentid, primaryleptrkid):
    if ( trkid == primaryleptrkid ) or ( parentid[trkid] == primaryleptrkid ): return True
    elif ( parentid[trkid] == -1 ): return False
    else: return IsFromPrimaryLep(parentid[trkid], parentid, primaryleptrkid)

# Dont care where the ND vertex goes as we dont use it for making the tdr reco-reco pairs
def IsInNDFV(pos_x_cm, pos_y_cm, pos_z_cm):
    return True
