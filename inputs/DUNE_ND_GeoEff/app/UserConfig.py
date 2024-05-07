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
# Cut on mu nn probability 
nd_contained_mu_cut = 0.5
nd_tracker_mu_cut = 0.5
# Hadronic veto selection crtieria
nd_veto_region_size = 30 # cm
nd_veto_threshold_energy = 30 # MeV
fd_veto_region_size = 30 # cm
fd_veto_threshold_energy = 30 # MeV

# Active volume for FD unit [cm]
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
FD_FV_min = [-145., -375., 100.]
FD_FV_max = [ 145.,  375., 960.]

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

def IsInNDFV(pos_x_cm, pos_y_cm, pos_z_cm):
    # This ND FV cut function is copied from CAFAna:
    # https://github.com/DUNE/lblpwgtools/blob/master/CAFAna/Cuts/TruthCuts.h#L65-L98
    # To confirm: 1.6cm or 1.3cm dead region b/t modules in x?
    # Function for on-axis!!! geo Eff
    inDeadRegion = False
    for i in range(-3, 3):
        # 0.5cm cathode in the middle of each module, plus 0.5cm buffer
        cathode_center = i * 102.1
        if (pos_x_cm > cathode_center - 0.75 and pos_x_cm < cathode_center + 0.75):
            inDeadRegion = True

        # 1.6cm dead region between modules (0.5cm module wall and 0.3cm pixel
        # plane, x2) don't worry about outer boundary because events are only
        # generated in active Ar + insides
        module_boundary = i * 102.1 + 51.05
        if (i <= 2 and pos_x_cm > module_boundary - 1.3 and pos_x_cm < module_boundary + 1.3):
            inDeadRegion = True

    for i in range(1, 4):
        # module boundaries in z are 1.8cm (0.4cm ArCLight plane + 0.5cm module
        # wall, x2) module is 102.1cm wide, but only 101.8cm long due to cathode
        # (0.5cm) being absent in length but ArCLight is 0.1cm thicker than pixel
        # plane so it's 0.3cm difference positions are off-set by 0.6 because I
        # defined 0 to be the upstream edge based on the active volume by
        # inspecting a plot, and aparently missed by 3 mm, but whatever add 8mm =
        # 2 pad buffer due to worse position resolution in spatial dimension z
        # compared to timing direction x so total FV gap will be 1.8 + 2*0.8
        # = 3.4cm
        module_boundary = i * 101.8 - 0.6
        if (pos_z_cm > module_boundary - 1.7 and pos_z_cm < module_boundary + 1.7):
            inDeadRegion = True

    # return fiducial volume
    return (abs(pos_x_cm) < ND_FV_max[0] and abs(pos_y_cm) < ND_FV_max[1] and pos_z_cm > ND_FV_min[2] and pos_z_cm < ND_FV_max[2] and inDeadRegion == False)
