#! /usr/bin/env python
"""
Runs over edepsim files and produce paired events in near and far detector.
"""
import time
start_time = time.time()

import numpy as np
import os
import os.path
import sys
import torch
from scipy.spatial.transform import Rotation as R
from muonEffModel import muonEffModel
import ROOT
from ROOT import TG4Event, TFile, TTree, TGraph
from ROOT import gROOT # for creating the output file
from array import array
from math import cos, sin
import random

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("input")
parser.add_argument("--config", type=str, default="UserConfig.py")
parser.add_argument(
    "--out_dir", type=str, default="/dune/app/users/weishi/testn2fd/DUNE_ND_GeoEff/app/output"
)
args = parser.parse_args()
a, config_file, out_path = args.input, args.config, args.out_dir

with open(config_file) as infile:
    exec(infile.read())
import pyGeoEff

ROOT.gROOT.ProcessLine('#include<vector>') # for outputting data to ROOT file

# this works interactively
#edep_file = TFile("edep.*.root")
edep_file = TFile(a, "OPEN")
if not edep_file:
    print ("Error: could not open file", a)
    sys.exit()
inputTree = edep_file.Get("EDepSimEvents")
genieTree = edep_file.Get("DetSimPassThru/gRooTracker")

event = TG4Event()
inputTree.SetBranchAddress("Event", event)

entries = inputTree.GetEntriesFast()
genie_entries = genieTree.GetEntriesFast()
if entries != genie_entries:
    print("Edep-sim tree and GENIE tree number of entries do not match!")
    sys.exit()

gDecayZ = TGraph(27, OffAxisPoints, meanPDPZ)
# Load muon neural network
net = muonEffModel()
net.load_state_dict(torch.load("./muonEff30.nn", map_location=torch.device('cpu')))
net.eval()

#########
# OUTPUT
#########

# create directory for plots to be stored if it doesn't already exist
plotpath = out_path + "/plots"
rootpath = out_path + "/root_out" # for output ROOT files
if not os.path.exists( out_path):
    os.makedirs( out_path)
    print("out_path '" + out_path + "' did not exist. It has been created!")
if not os.path.exists( plotpath):
    os.makedirs( plotpath)
    print("plotpath '" + plotpath + "' did not exist. It has been created!")
if not os.path.exists( rootpath):
    os.makedirs( rootpath)
    print("rootpath '" + rootpath + "' did not exist. It has been created!")
print(" \n") # separate output

# Create the ROOT file that will hold the output of this script
f_out = TFile('{0}/n2fd_paired_out.root'.format(rootpath), 'RECREATE')
myEvents = TTree('myEvents', 'myEvents')
maxEdeps = 100000 # max number of edeps for initialize arrays
maxGenieParts = 1000

##################################
# Genie evt info
##################################
Genie_nParts = array('i', [0])
myEvents.Branch('Genie_nParts', Genie_nParts, 'Genie_nParts/I')
GenieParts_pdg = np.zeros((maxGenieParts,), dtype=np.int32)
myEvents.Branch('GenieParts_pdg', GenieParts_pdg, 'GenieParts_pdg[Genie_nParts]/I')
GenieParts_p0_MeV = np.zeros((maxGenieParts,), dtype=np.float32)
myEvents.Branch('GenieParts_p0_MeV', GenieParts_p0_MeV, 'GenieParts_p0_MeV[Genie_nParts]/F')
GenieParts_p1_MeV = np.zeros((maxGenieParts,), dtype=np.float32)
myEvents.Branch('GenieParts_p1_MeV', GenieParts_p1_MeV, 'GenieParts_p1_MeV[Genie_nParts]/F')
GenieParts_p2_MeV = np.zeros((maxGenieParts,), dtype=np.float32)
myEvents.Branch('GenieParts_p2_MeV', GenieParts_p2_MeV, 'GenieParts_p2_MeV[Genie_nParts]/F')
GenieParts_p3_MeV = np.zeros((maxGenieParts,), dtype=np.float32)
myEvents.Branch('GenieParts_p3_MeV', GenieParts_p3_MeV, 'GenieParts_p3_MeV[Genie_nParts]/F')
Genie_final_lep_pdg = array('i', [0])
myEvents.Branch('Genie_final_lep_pdg', Genie_final_lep_pdg, 'Genie_final_lep_pdg/I')
Genie_final_lep_p0_MeV = array('f', [0.0])
myEvents.Branch('Genie_final_lep_p0_MeV', Genie_final_lep_p0_MeV, 'Genie_final_lep_p0_MeV/F')
Genie_final_lep_p1_MeV = array('f', [0.0])
myEvents.Branch('Genie_final_lep_p1_MeV', Genie_final_lep_p1_MeV, 'Genie_final_lep_p1_MeV/F')
Genie_final_lep_p2_MeV = array('f', [0.0])
myEvents.Branch('Genie_final_lep_p2_MeV', Genie_final_lep_p2_MeV, 'Genie_final_lep_p2_MeV/F')
Genie_final_lep_p3_MeV = array('f', [0.0])
myEvents.Branch('Genie_final_lep_p3_MeV', Genie_final_lep_p3_MeV, 'Genie_final_lep_p3_MeV/F')

##################################
# LArBath evt info
##################################
nEdeps = array('i', [0])
myEvents.Branch('nEdeps', nEdeps, 'nEdeps/I')
deps_trkID = np.zeros((maxEdeps,), dtype=np.int32)
myEvents.Branch('deps_trkID', deps_trkID, 'deps_trkID[nEdeps]/I')
deps_parentID = np.zeros((maxEdeps,), dtype=np.int32)
myEvents.Branch('deps_parentID', deps_parentID, 'deps_parentID[nEdeps]/I')
deps_pdg = np.zeros((maxEdeps,), dtype=np.int32)
myEvents.Branch('deps_pdg', deps_pdg, 'deps_pdg[nEdeps]/I')
deps_E_MeV = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('deps_E_MeV', deps_E_MeV, 'deps_E_MeV[nEdeps]/F')
deps_start_t_us = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('deps_start_t_us', deps_start_t_us, 'deps_start_t_us[nEdeps]/F')
deps_stop_t_us = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('deps_stop_t_us', deps_stop_t_us, 'deps_stop_t_us[nEdeps]/F')
larbath_vtx_cm = array('f', 3*[0.0])
myEvents.Branch('larbath_vtx_cm', larbath_vtx_cm, 'larbath_vtx_cm[3]/F')
# edeps generated in LArBath (start point)
larbath_deps_start_x_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('larbath_deps_start_x_cm', larbath_deps_start_x_cm, 'larbath_deps_start_x_cm[nEdeps]/F') # larbath edeps x
larbath_deps_start_y_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('larbath_deps_start_y_cm', larbath_deps_start_y_cm, 'larbath_deps_start_y_cm[nEdeps]/F')
larbath_deps_start_z_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('larbath_deps_start_z_cm', larbath_deps_start_z_cm, 'larbath_deps_start_z_cm[nEdeps]/F')
# edeps generated in LArBath (stop point)
larbath_deps_stop_x_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('larbath_deps_stop_x_cm', larbath_deps_stop_x_cm, 'larbath_deps_stop_x_cm[nEdeps]/F')
larbath_deps_stop_y_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('larbath_deps_stop_y_cm', larbath_deps_stop_y_cm, 'larbath_deps_stop_y_cm[nEdeps]/F')
larbath_deps_stop_z_cm = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('larbath_deps_stop_z_cm', larbath_deps_stop_z_cm, 'larbath_deps_stop_z_cm[nEdeps]/F')

####################
# ND evt non-ecc
####################
# Random thrown vertex in ND (paired evt)
nd_vtx_cm_nonecc = array('f', 3*[0.0])
myEvents.Branch('nd_vtx_cm_nonecc', nd_vtx_cm_nonecc, 'nd_vtx_cm_nonecc[3]/F')
# Muon selection probability in ND LAr and tracker
nd_lep_contained_prob_nonecc = array('f', [0.0])
myEvents.Branch('nd_lep_contained_prob_nonecc', nd_lep_contained_prob_nonecc, 'nd_lep_contained_prob_nonecc/F')
nd_lep_tracker_prob_nonecc = array('f', [0.0])
myEvents.Branch('nd_lep_tracker_prob_nonecc', nd_lep_tracker_prob_nonecc, 'nd_lep_tracker_prob_nonecc/F')
nd_lep_ke_MeV_exit_ndlar_nonecc = array('f', [0.0])
myEvents.Branch('nd_lep_ke_MeV_exit_ndlar_nonecc', nd_lep_ke_MeV_exit_ndlar_nonecc, 'nd_lep_ke_MeV_exit_ndlar_nonecc/F')
# Edeps in ND random throw (start points)
nd_deps_start_x_cm_nonecc = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('nd_deps_start_x_cm_nonecc', nd_deps_start_x_cm_nonecc, 'nd_deps_start_x_cm_nonecc[nEdeps]/F')
nd_deps_start_y_cm_nonecc = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('nd_deps_start_y_cm_nonecc', nd_deps_start_y_cm_nonecc, 'nd_deps_start_y_cm_nonecc[nEdeps]/F')
nd_deps_start_z_cm_nonecc = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('nd_deps_start_z_cm_nonecc', nd_deps_start_z_cm_nonecc, 'nd_deps_start_z_cm_nonecc[nEdeps]/F')
nd_deps_start_InNDLAr_nonecc = np.zeros((maxEdeps,), dtype=np.int32)
myEvents.Branch('nd_deps_start_InNDLAr_nonecc', nd_deps_start_InNDLAr_nonecc, 'nd_deps_start_InNDLAr_nonecc[nEdeps]/I')
# Edeps in ND random throw (stop points)
nd_deps_stop_x_cm_nonecc = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('nd_deps_stop_x_cm_nonecc', nd_deps_stop_x_cm_nonecc, 'nd_deps_stop_x_cm_nonecc[nEdeps]/F')
nd_deps_stop_y_cm_nonecc = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('nd_deps_stop_y_cm_nonecc', nd_deps_stop_y_cm_nonecc, 'nd_deps_stop_y_cm_nonecc[nEdeps]/F')
nd_deps_stop_z_cm_nonecc = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('nd_deps_stop_z_cm_nonecc', nd_deps_stop_z_cm_nonecc, 'nd_deps_stop_z_cm_nonecc[nEdeps]/F')
nd_deps_stop_InNDLAr_nonecc = np.zeros((maxEdeps,), dtype=np.int32)
myEvents.Branch('nd_deps_stop_InNDLAr_nonecc', nd_deps_stop_InNDLAr_nonecc, 'nd_deps_stop_InNDLAr_nonecc[nEdeps]/I')

##################################
# FD paired evt with nd non-ecc
##################################
# Random thrown vertex in FD
fd_vtx_cm_pair_nd_nonecc = array('f', 3*[0.0])
myEvents.Branch('fd_vtx_cm_pair_nd_nonecc', fd_vtx_cm_pair_nd_nonecc, 'fd_vtx_cm_pair_nd_nonecc[3]/F')
# Edeps in FD random throw (start points)
fd_deps_start_x_cm_pair_nd_nonecc = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('fd_deps_start_x_cm_pair_nd_nonecc', fd_deps_start_x_cm_pair_nd_nonecc, 'fd_deps_start_x_cm_pair_nd_nonecc[nEdeps]/F')
fd_deps_start_y_cm_pair_nd_nonecc = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('fd_deps_start_y_cm_pair_nd_nonecc', fd_deps_start_y_cm_pair_nd_nonecc, 'fd_deps_start_y_cm_pair_nd_nonecc[nEdeps]/F')
fd_deps_start_z_cm_pair_nd_nonecc = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('fd_deps_start_z_cm_pair_nd_nonecc', fd_deps_start_z_cm_pair_nd_nonecc, 'fd_deps_start_z_cm_pair_nd_nonecc[nEdeps]/F')
# Edeps in FD random throw (stop points)
fd_deps_stop_x_cm_pair_nd_nonecc = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('fd_deps_stop_x_cm_pair_nd_nonecc', fd_deps_stop_x_cm_pair_nd_nonecc, 'fd_deps_stop_x_cm_pair_nd_nonecc[nEdeps]/F')
fd_deps_stop_y_cm_pair_nd_nonecc = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('fd_deps_stop_y_cm_pair_nd_nonecc', fd_deps_stop_y_cm_pair_nd_nonecc, 'fd_deps_stop_y_cm_pair_nd_nonecc[nEdeps]/F')
fd_deps_stop_z_cm_pair_nd_nonecc = np.zeros((maxEdeps,), dtype=np.float32)
myEvents.Branch('fd_deps_stop_z_cm_pair_nd_nonecc', fd_deps_stop_z_cm_pair_nd_nonecc, 'fd_deps_stop_z_cm_pair_nd_nonecc[nEdeps]/F')

nd_fd_throws_passed = array('i', [0])
myEvents.Branch('nd_fd_throws_passed', nd_fd_throws_passed, 'nd_fd_throws_passed/I')

###########################
# Loop over edepsim events
##########################

for jentry in range(entries):
    print("jentry = " + str(jentry))
    nb = inputTree.GetEntry(jentry)
    gb = genieTree.GetEntry(jentry)
    #print("nb =" + str(nb))
    #print("event number: ", event.EventId)
    #print("number of primaries: ", event.Primaries.size())
    #print("number of trajectories: ", event.Trajectories.size())
    #print("number of segments: ", event.SegmentDetectors.size())
    all_dep_startpos_list = list()
    all_dep_stoppos_list = list()
    all_dep_trkID_list  = list()
    all_dep_parentID_list = list()
    all_dep_pdg_list = list()
    all_edep_list   = list()
    all_dep_starttime_list = list()
    all_dep_stoptime_list = list()
    had_dep_pos_list = list()
    had_edep_list   = list()
    lep_dep_pos_list = list()
    lep_edep_list   = list()
    dep_startpos_in_ndlar_list = list()
    dep_stoppos_in_ndlar_list = list()

    # initialize
    # Define here but do not write out
    nLepEdeps = array('i', [0])
    lepdeps_E_MeV = np.zeros((maxEdeps,), dtype=np.float32)
    nd_lepdeps_x_cm_nonecc = np.zeros((maxEdeps,), dtype=np.float32)
    nd_lepdeps_y_cm_nonecc = np.zeros((maxEdeps,), dtype=np.float32)
    nd_lepdeps_z_cm_nonecc = np.zeros((maxEdeps,), dtype=np.float32)
    larbath_vtx_cm[0] = 0; larbath_vtx_cm[1] = 0; larbath_vtx_cm[2] = 0;
    nEdeps[0] = 0
    nLepEdeps[0] = 0
    nLepEdeps_insideNDLAr = 0
    nd_lep_contained_prob_nonecc[0] = -1; nd_lep_tracker_prob_nonecc[0] = -1;
    nd_lep_ke_MeV_exit_ndlar_nonecc[0] = 0
    nd_vtx_cm_nonecc[0] = 0; nd_vtx_cm_nonecc[1] = 0; nd_vtx_cm_nonecc[2] = 0;
    fd_vtx_cm_pair_nd_nonecc[0] = 0; fd_vtx_cm_pair_nd_nonecc[1] = 0; fd_vtx_cm_pair_nd_nonecc[2] = 0;
    nd_fd_throws_passed[0] = 0

    for primary in event.Primaries:
        #print("number of particles: ", primary.Particles.size())
        for ipart, particle in enumerate(primary.Particles):
            #print("ipart here = " +str(ipart))
            PDGCode = particle.GetPDGCode()
            #print("pdgcode: ", PDGCode)
            if abs(PDGCode) >= 11 and abs(PDGCode) <= 16:
                posx = primary.GetPosition().X() * edep2cm
                posy = primary.GetPosition().Y() * edep2cm
                posz = primary.GetPosition().Z() * edep2cm
                momx = particle.GetMomentum().X() # MeV
                momy = particle.GetMomentum().Y()
                momz = particle.GetMomentum().Z()
                # Generated lepton momentum in edepsim
                lep_p = [momx/gev2mev, momy/gev2mev, momz/gev2mev]
                PrimaryLepTrackID = particle.GetTrackId()
                print("edepsim lep pos x: ", posx, "  y: ", posy, "  z: ", posz, " [cm]")
                print("edepsim lep    px: ", lep_p[0], " py: ", lep_p[1], " pz: ", lep_p[2], " [GeV]")
                """
                print("pos x: ", posx, " y: ", posy, " z: ", posz, " [cm]")
                print("mom x: ", momx, " y: ", momy, " z: ", momz, " [MeV]")
                print("PrimaryLepTrackID: ", PrimaryLepTrackID)
                """

    trajectories_parentid = np.empty(len(event.Trajectories), dtype=np.int32)
    trajectories_pdg = np.empty(len(event.Trajectories), dtype=np.int32)
    for iTraj, trajectory in enumerate(event.Trajectories):
        # iTraj is same as TrackId, starts from #0
        trajectories_parentid[iTraj] = trajectory.GetParentId()
        trajectories_pdg[iTraj] = trajectory.GetPDGCode()
        #print("iTraj #", iTraj, " parentID: ", trajectories_parentid[iTraj], " pdg: ", trajectories_pdg[iTraj] )

    for containerName, hitSegments in event.SegmentDetectors:
        # iHit is the index, hitSEgment is the data stored at the index in the second item in event.SegementDectectors
        for iHit, hitSegment in enumerate(hitSegments):

            # Energy deposit from primary particles
            edep_start_x = hitSegment.GetStart().X() * edep2cm
            edep_start_y = hitSegment.GetStart().Y() * edep2cm
            edep_start_z = hitSegment.GetStart().Z() * edep2cm
            edep_start_t = hitSegment.GetStart().T() * edep2us
            edep_stop_x = hitSegment.GetStop().X() * edep2cm
            edep_stop_y = hitSegment.GetStop().Y() * edep2cm
            edep_stop_z = hitSegment.GetStop().Z() * edep2cm
            edep_stop_t = hitSegment.GetStop().T() * edep2us
            edep_x = (edep_start_x + edep_stop_x)/2 # use this for hadronic veto
            edep_y = (edep_start_y + edep_stop_y)/2
            edep_z = (edep_start_z + edep_stop_z)/2
            edep_trkID = hitSegment.Contrib[0] # what is different if use GetPrimaryId()?
            edep_parentID = trajectories_parentid[edep_trkID]
            edep_pdg = trajectories_pdg[edep_trkID]
            edep = hitSegment.GetEnergyDeposit()

            all_dep_startpos_list.append(edep_start_x)
            all_dep_startpos_list.append(edep_start_y)
            all_dep_startpos_list.append(edep_start_z)
            all_dep_starttime_list.append(edep_start_t)
            all_dep_stoppos_list.append(edep_stop_x)
            all_dep_stoppos_list.append(edep_stop_y)
            all_dep_stoppos_list.append(edep_stop_z)
            all_dep_stoptime_list.append(edep_stop_t)
            all_dep_trkID_list.append(edep_trkID)
            all_dep_parentID_list.append(edep_parentID)
            all_dep_pdg_list.append(edep_pdg)
            all_edep_list.append(edep)

            # Separate hadronic and leptonic part of edepsim
            # Later need these hadronic edeps to evaluate hadronic veto
            # Need to calculate remaining kinetic energy for muons exit ND LAr
            if IsFromPrimaryLep(edep_trkID, trajectories_parentid, PrimaryLepTrackID) == False:
                had_dep_pos_list.append(edep_x)
                had_dep_pos_list.append(edep_y)
                had_dep_pos_list.append(edep_z)
                had_edep_list.append(edep)
            else:
                lep_dep_pos_list.append(edep_x)
                lep_dep_pos_list.append(edep_y)
                lep_dep_pos_list.append(edep_z)
                lep_edep_list.append(edep)

    # for use in processing events before and after transformations
    nEdeps[0] = len(all_edep_list)
    nLepEdeps[0] = len(lep_edep_list)

    # Also save truth info from GENIE
    Genie_nParts[0] = 0
    Genie_final_lep_pdg[0] = 0
    Genie_final_lep_p0_MeV[0] = 0
    Genie_final_lep_p1_MeV[0] = 0
    Genie_final_lep_p2_MeV[0] = 0
    Genie_final_lep_p3_MeV[0] = 0
    all_genie_pdg_list = list()
    all_genie_p0_list = list()
    all_genie_p1_list = list()
    all_genie_p2_list = list()
    all_genie_p3_list = list()
    for p in range(genieTree.StdHepN):
        # Get only final state particles
        if genieTree.StdHepStatus[p] == 1:
            #print("GENIE pdg: ", genieTree.StdHepPdg[p], " status: ", genieTree.StdHepStatus[p], " p0: ", genieTree.StdHepP4[p*4 + 0]*gev2mev, " p1: ", genieTree.StdHepP4[p*4 + 1]*gev2mev, " p2: ", genieTree.StdHepP4[p*4 + 2]*gev2mev, " p3: ", genieTree.StdHepP4[p*4 + 3]*gev2mev)
            all_genie_pdg_list.append(genieTree.StdHepPdg[p])
            all_genie_p0_list.append(genieTree.StdHepP4[p*4 + 0]*gev2mev)
            all_genie_p1_list.append(genieTree.StdHepP4[p*4 + 1]*gev2mev)
            all_genie_p2_list.append(genieTree.StdHepP4[p*4 + 2]*gev2mev)
            all_genie_p3_list.append(genieTree.StdHepP4[p*4 + 3]*gev2mev)

            if abs(genieTree.StdHepPdg[p]) >= 11 and abs(genieTree.StdHepPdg[p]) <= 16:
                print("GENIE final lep pdg: ", genieTree.StdHepPdg[p], " status: ", genieTree.StdHepStatus[p], " p0 MeV: ", genieTree.StdHepP4[p*4 + 0]*gev2mev, " p1: ", genieTree.StdHepP4[p*4 + 1]*gev2mev, " p2: ", genieTree.StdHepP4[p*4 + 2]*gev2mev, " p3: ", genieTree.StdHepP4[p*4 + 3]*gev2mev)
                Genie_final_lep_pdg[0] = genieTree.StdHepPdg[p]
                Genie_final_lep_p0_MeV[0] = genieTree.StdHepP4[p*4 + 0]*gev2mev
                Genie_final_lep_p1_MeV[0] = genieTree.StdHepP4[p*4 + 1]*gev2mev
                Genie_final_lep_p2_MeV[0] = genieTree.StdHepP4[p*4 + 2]*gev2mev
                Genie_final_lep_p3_MeV[0] = genieTree.StdHepP4[p*4 + 3]*gev2mev
    # End particle loop in genie
    Genie_nParts[0] = len(all_genie_pdg_list)

    ##################################################################
    # Initialize geometric efficiency module to manipulate energy deps
    # Only consider ND LAr on axis
    ##################################################################
    seed = random.randrange(1024)
    print ("-- seed:", seed)

    geoEff = pyGeoEff.geoEff(seed, False)

    # Use neutrino decay position, rather than fixed neutrino direction as symmetry axis
    geoEff.setUseFixedBeamDir(False)

    # Near detector active dimensions for hadronic veto
    geoEff.setActiveX(NDActiveVol_min[0], NDActiveVol_max[0])
    geoEff.setActiveY(NDActiveVol_min[1], NDActiveVol_max[1])
    geoEff.setActiveZ(NDActiveVol_min[2], NDActiveVol_max[2])

    # Far detector active dimensions for hadronic veto
    geoEff.setFDActiveX(FDActiveVol_min[0], FDActiveVol_max[0])
    geoEff.setFDActiveY(FDActiveVol_min[1], FDActiveVol_max[1])
    geoEff.setFDActiveZ(FDActiveVol_min[2], FDActiveVol_max[2])

    # Range for random translation throws in ND fiducial volume
    geoEff.setRangeX(ND_FV_min[0], ND_FV_max[0])
    geoEff.setRangeY(ND_FV_min[1], ND_FV_max[1])
    geoEff.setRangeZ(ND_FV_min[2], ND_FV_max[2])

    # Range for random translation throws in FD fiducial volume
    geoEff.setRangeXFD(FD_FV_min[0], FD_FV_max[0])
    geoEff.setRangeYFD(FD_FV_min[1], FD_FV_max[1])
    geoEff.setRangeZFD(FD_FV_min[2], FD_FV_max[2])

    # Set offset between ND MC coordinate system and volumes defined above.
    # Is this still needed for Alex's Genie Gen??? To be validated/discussed
    geoEff.setOffsetX(NDLAr_OnAxis_offset[0])
    geoEff.setOffsetY(NDLAr_OnAxis_offset[1])
    geoEff.setOffsetZ(NDLAr_OnAxis_offset[2])

    # Original Genie Gen event vtx, not guaranteed to be at (0,0,0)!
    # Later we need the vertex to calculate the rotation axis when event is randomly put at a new position in ND
    geoEff.setVertex(posx, posy, posz) #cm

    # Interpolate event neutrino production point (beam coordinate)
    # Input needs to be in unit of meters
    decayZbeamCoord = gDecayZ.Eval( posx/100. - NDLAr_OnAxis_offset[0]/100. - detRefBeamCoord[0] )

    # Calculate neutrino production point in detector coordinate
    decayXdetCoord = beamRefDetCoord[0] - detRefBeamCoord[0]
    decayYdetCoord = beamRefDetCoord[1] - detRefBeamCoord[1]*cos(beamLineRotation) + ( decayZbeamCoord - detRefBeamCoord[2] )*sin(beamLineRotation)
    decayZdetCoord = beamRefDetCoord[2] + detRefBeamCoord[1]*sin(beamLineRotation) + ( decayZbeamCoord - detRefBeamCoord[2] )*cos(beamLineRotation)

    # Set production point in unit: cm
    # Later we need this decay position to calculate the rotation axis when event is randomly put at a new position in ND
    geoEff.setDecayPos(decayXdetCoord*100., decayYdetCoord*100., decayZdetCoord*100.)

    # Randomly throw the Genie Gen event in ND
    tot_nd_throw = 0

    # If after max_nd_throws still don't pass at nd, stop and move to next event (otherwise too much computing resources)
    while tot_nd_throw < max_nd_throws:
        print ("-- tot nd throw:", tot_nd_throw)
        # Configure veto region for hadronic veto for ND throws
        geoEff.setVetoSizes([nd_veto_region_size])
        # Configure E threshold for hadronic veto for ND throws
        geoEff.setVetoEnergyThresholds([nd_veto_threshold_energy])

        ####################################
        # Only do one throw in ND at a time
        ####################################
        geoEff.setNthrows(1)
        geoEff.throwTransforms() # this randomly generates new vtx position and a rotation angle w.r.t. the neutrino direction

        # Get the randomly generated vtx x, y, z, and the angle
        throwVtxX_nd = geoEff.getCurrentThrowTranslationsX() # cm
        throwVtxY_nd = geoEff.getCurrentThrowTranslationsY() # cm
        throwVtxZ_nd = geoEff.getCurrentThrowTranslationsZ() # cm
        throwAngle = geoEff.getCurrentThrowRotations()

        # Require random thrown vtx pos outside in ND dead regions
        if ( IsInNDFV( throwVtxX_nd[0] - NDLAr_OnAxis_offset[0], throwVtxY_nd[0] - NDLAr_OnAxis_offset[1], throwVtxZ_nd[0] - NDLAr_OnAxis_offset[2]) ):
            print ("-- nd throw", tot_nd_throw, "is in nd FV")

            # Interpolate neutrino production point (beam coordinate) for the random throw, unit meter
            RandomthrowdecayZbeamCoord = gDecayZ.Eval( throwVtxX_nd[0]/100. - NDLAr_OnAxis_offset[0]/100. - detRefBeamCoord[0] )

            # Calculate neutrino production point in detector coordinate, unit meter
            RandomthrowdecayXdetCoord = beamRefDetCoord[0] - detRefBeamCoord[0]
            RandomthrowdecayYdetCoord = beamRefDetCoord[1] - detRefBeamCoord[1]*cos(beamLineRotation) + ( RandomthrowdecayZbeamCoord - detRefBeamCoord[2] )*sin(beamLineRotation)
            RandomthrowdecayZdetCoord = beamRefDetCoord[2] + detRefBeamCoord[1]*sin(beamLineRotation) + ( RandomthrowdecayZbeamCoord - detRefBeamCoord[2] )*cos(beamLineRotation)

            # We have a new decay position because the above random throw changes x in ND
            # Need this to calculate the rotation axis for the random throw
            geoEff.setDecayPos4RandomThrowX(RandomthrowdecayXdetCoord*100., RandomthrowdecayYdetCoord*100., RandomthrowdecayZdetCoord*100.)

            # Set the hadronic part list of edeps to evaluate ND hadronic veto
            # NOTE: in the future may need to add lepton selection at ND as well (to be discussed)
            geoEff.setHitSegEdeps(had_edep_list)
            geoEff.setHitSegPoss(had_dep_pos_list)
            ndrandthrowresulthad = geoEff.getNDContainment4RandomThrowX() # returns a struct

            if (ndrandthrowresulthad.containresult[0][0][0] != 0):
                print ("-- nd throw", tot_nd_throw, "passed nd had veto")
                ###########################################################
                # Random throw passed ND hadronic veto !
                ###########################################################

                # Now change to the full list of edeps start points
                # the random thrown x/y/z/angle should remain the same because throw was done above already
                geoEff.setHitSegEdeps(all_edep_list)
                geoEff.setHitSegPoss(all_dep_startpos_list)
                # And call the function again to get new transformed positions for all edeps
                ndrandthrowresultall_start = geoEff.getNDContainment4RandomThrowX()

                # Repeat for edepsim stop points !!!
                geoEff.setHitSegEdeps(all_edep_list)
                geoEff.setHitSegPoss(all_dep_stoppos_list)
                ndrandthrowresultall_stop = geoEff.getNDContainment4RandomThrowX()

                # Repeat for lepton edeps(probably don't care start/stop points)
                # we need this to calculate muon ke after it exits NDLAr active vol
                geoEff.setHitSegEdeps(lep_edep_list)
                geoEff.setHitSegPoss(lep_dep_pos_list)
                ndrandthrowresultlep = geoEff.getNDContainment4RandomThrowX()

                ###########################################################
                #  (Oct 25, 2023)
                # Since we don't have a good TMS reco now, interested to know
                # probability of muons selected in tracker and its kinetic energy after exiting NDLAr
                ###########################################################

                # Perform same transformation on lepton momentum vec: needed for evaluate muon nn
                # this should go in a function similar to hadron manipulation
                # input would be decayToVertex and decayToTranslated

                # Vector from neutrino production point to original event vertex
                decayToVertex     = [posx - decayXdetCoord*100, posy - decayYdetCoord*100, posz - decayZdetCoord*100]
                # Vector from neutrino production point to randomly thrown vertex.
                decayToTranslated = [throwVtxX_nd[0] - RandomthrowdecayXdetCoord*100, throwVtxY_nd[0] - RandomthrowdecayYdetCoord*100, throwVtxZ_nd[0] - RandomthrowdecayZdetCoord*100]

                magDecayToVertex     = np.sqrt(np.sum(np.square(decayToVertex)))
                magDecayToTranslated = np.sqrt(np.sum(np.square(decayToTranslated)))

                translationAngle = np.dot(decayToTranslated, decayToVertex)
                translationAngle = np.divide(translationAngle, np.multiply(magDecayToVertex, magDecayToTranslated));
                #for angleval in translationAngle:if angleval<=-1 or angleval>=1: print(i_event, angleval)
                translationAngle = np.arccos(translationAngle);
                translationAxis = np.cross(decayToTranslated, decayToVertex)
                translationAxis = translationAxis/np.linalg.norm(translationAxis)
                translation_rot_vec = np.multiply(translationAxis, translationAngle)

                decayToTranslated = decayToTranslated/np.linalg.norm(decayToTranslated)
                phi_rot_vec = np.multiply(decayToTranslated, throwAngle)

                # Get rotation matrices due to:
                # Vertex translation (which "rotates" the average neutrino direction)
                translation_rot = R.from_rotvec(translation_rot_vec)
                # Random phi rotation around average neutrino direction
                phi_rot = R.from_rotvec(phi_rot_vec)
                # Rotate momentum
                lep_p = translation_rot.apply(lep_p)
                lep_p = phi_rot.apply(lep_p)

                print ("--- mu pos x cm: ", throwVtxX_nd[0] - NDLAr_OnAxis_offset[0], ", y: ", throwVtxY_nd[0] - NDLAr_OnAxis_offset[1], ", z: ", throwVtxZ_nd[0] - NDLAr_OnAxis_offset[2])
                print ("--- mu px GeV: ", lep_p[0], ", py: ", lep_p[1], ", pz: ", lep_p[2])

                # Input features for nn: momentum and vertex
                inputfeatures = np.column_stack((lep_p[0], lep_p[1], lep_p[2], throwVtxX_nd[0] - NDLAr_OnAxis_offset[0], throwVtxY_nd[0] - NDLAr_OnAxis_offset[1], throwVtxZ_nd[0] - NDLAr_OnAxis_offset[2]))
                # Convert to Pytorch tensor
                inputfeatures = torch.as_tensor(inputfeatures).type(torch.FloatTensor)
                # Evaluate neural network
                with torch.no_grad() :
                    netOut = net(inputfeatures)
                    netOut = torch.nn.functional.softmax(netOut, dim=1).detach().numpy()
                # Get contained probability for this throw
                nd_lep_contained_prob_nonecc[0] = np.array(netOut[:,0], dtype = float)
                # Get tracker probability for this throw
                nd_lep_tracker_prob_nonecc[0] = np.array(netOut[:,1], dtype = float)
                print ("--- nn prob. contained: ", nd_lep_contained_prob_nonecc[0], ", tracker: ", nd_lep_tracker_prob_nonecc[0])

                # Apply lepton selection for the throw
                if nd_lep_contained_prob_nonecc[0] >= nd_contained_mu_cut or nd_lep_tracker_prob_nonecc[0] >= nd_tracker_mu_cut:
                    # All info we need for transformed lep e deposits
                    lepdeps_E_MeV[:nLepEdeps[0]] = np.array(lep_edep_list, dtype=np.float32)
                    nd_lepdeps_x_cm_nonecc[:nLepEdeps[0]] = np.array(ndrandthrowresultlep.thrownEdepspos[0][0,:], dtype=np.float32)
                    nd_lepdeps_y_cm_nonecc[:nLepEdeps[0]] = np.array(ndrandthrowresultlep.thrownEdepspos[0][1,:], dtype=np.float32)
                    nd_lepdeps_z_cm_nonecc[:nLepEdeps[0]] = np.array(ndrandthrowresultlep.thrownEdepspos[0][2,:], dtype=np.float32)
                    # Now loop over nLepEdeps
                    print ("---- lep array KE tot (MeV):", np.sum(lepdeps_E_MeV))
                    print ("---- lep list KE tot (MeV):", np.sum(lep_edep_list))
                    print ("---- nLepEdeps:", nLepEdeps[0])
                    for iedep in range(0, nLepEdeps[0]):
                        # Deposit inside ND LAr, continue to next
                        if nd_lepdeps_x_cm_nonecc[iedep]-NDLAr_OnAxis_offset[0] > NDActiveVol_min[0] and nd_lepdeps_x_cm_nonecc[iedep]-NDLAr_OnAxis_offset[0] < NDActiveVol_max[0] and \
                           nd_lepdeps_y_cm_nonecc[iedep]-NDLAr_OnAxis_offset[1] > NDActiveVol_min[1] and nd_lepdeps_y_cm_nonecc[iedep]-NDLAr_OnAxis_offset[1] < NDActiveVol_max[1] and \
                           nd_lepdeps_z_cm_nonecc[iedep]-NDLAr_OnAxis_offset[2] > NDActiveVol_min[2] and nd_lepdeps_z_cm_nonecc[iedep]-NDLAr_OnAxis_offset[2] < NDActiveVol_max[2]:
                           nLepEdeps_insideNDLAr = nLepEdeps_insideNDLAr + 1
                           continue

                        # if not in ND LAr active vol, add up energy
                        nd_lep_ke_MeV_exit_ndlar_nonecc[0] = nd_lep_ke_MeV_exit_ndlar_nonecc[0] + lepdeps_E_MeV[iedep]

                    print ("---- lep number of edeps inside active NDLAr:", nLepEdeps_insideNDLAr)
                    print ("---- lep KE outside active NDLAr (MeV):", nd_lep_ke_MeV_exit_ndlar_nonecc[0])

                    # This set of info will be saved to ttree
                    nd_deps_start_x_cm_nonecc[:nEdeps[0]] = np.array(ndrandthrowresultall_start.thrownEdepspos[0][0,:], dtype=np.float32)
                    nd_deps_start_y_cm_nonecc[:nEdeps[0]] = np.array(ndrandthrowresultall_start.thrownEdepspos[0][1,:], dtype=np.float32)
                    nd_deps_start_z_cm_nonecc[:nEdeps[0]] = np.array(ndrandthrowresultall_start.thrownEdepspos[0][2,:], dtype=np.float32)
                    nd_deps_stop_x_cm_nonecc[:nEdeps[0]] = np.array(ndrandthrowresultall_stop.thrownEdepspos[0][0,:], dtype=np.float32)
                    nd_deps_stop_y_cm_nonecc[:nEdeps[0]] = np.array(ndrandthrowresultall_stop.thrownEdepspos[0][1,:], dtype=np.float32)
                    nd_deps_stop_z_cm_nonecc[:nEdeps[0]] = np.array(ndrandthrowresultall_stop.thrownEdepspos[0][2,:], dtype=np.float32)

                    # Add flags to inform downstream LAr-nd-sim which deposits are in active NDLAr volume
                    for iedep in range(0, nEdeps[0]):
                        # Deposit inside ND LAr
                        if nd_deps_start_x_cm_nonecc[iedep]-NDLAr_OnAxis_offset[0] > NDActiveVol_min[0] and nd_deps_start_x_cm_nonecc[iedep]-NDLAr_OnAxis_offset[0] < NDActiveVol_max[0] and \
                           nd_deps_start_y_cm_nonecc[iedep]-NDLAr_OnAxis_offset[1] > NDActiveVol_min[1] and nd_deps_start_y_cm_nonecc[iedep]-NDLAr_OnAxis_offset[1] < NDActiveVol_max[1] and \
                           nd_deps_start_z_cm_nonecc[iedep]-NDLAr_OnAxis_offset[2] > NDActiveVol_min[2] and nd_deps_start_z_cm_nonecc[iedep]-NDLAr_OnAxis_offset[2] < NDActiveVol_max[2]:

                           dep_startpos_in_ndlar_list.append(1)
                        else:
                           dep_startpos_in_ndlar_list.append(0)

                        # repeat for stopping point
                        if nd_deps_stop_x_cm_nonecc[iedep]-NDLAr_OnAxis_offset[0] > NDActiveVol_min[0] and nd_deps_stop_x_cm_nonecc[iedep]-NDLAr_OnAxis_offset[0] < NDActiveVol_max[0] and \
                           nd_deps_stop_y_cm_nonecc[iedep]-NDLAr_OnAxis_offset[1] > NDActiveVol_min[1] and nd_deps_stop_y_cm_nonecc[iedep]-NDLAr_OnAxis_offset[1] < NDActiveVol_max[1] and \
                           nd_deps_stop_z_cm_nonecc[iedep]-NDLAr_OnAxis_offset[2] > NDActiveVol_min[2] and nd_deps_stop_z_cm_nonecc[iedep]-NDLAr_OnAxis_offset[2] < NDActiveVol_max[2]:
                           dep_stoppos_in_ndlar_list.append(1)
                        else:
                           dep_stoppos_in_ndlar_list.append(0)

                    nd_deps_start_InNDLAr_nonecc[:nEdeps[0]] = np.array(dep_startpos_in_ndlar_list, dtype=np.int32)
                    nd_deps_stop_InNDLAr_nonecc[:nEdeps[0]] = np.array(dep_stoppos_in_ndlar_list, dtype=np.int32)
                    print ("---- nd_deps_start_InNDLAr_nonecc:", nd_deps_start_InNDLAr_nonecc[:nEdeps[0]], "size: ", nd_deps_start_InNDLAr_nonecc.size, "len of list: ", len(dep_startpos_in_ndlar_list))
                    print ("---- nd_deps_stop_InNDLAr_nonecc:", nd_deps_stop_InNDLAr_nonecc[:nEdeps[0]], "size: ", nd_deps_stop_InNDLAr_nonecc.size, "len of list: ", len(dep_stoppos_in_ndlar_list))

                    ###########################################################
                    # Now translate this event vertex to ND det coordinate (0,0,0) (and the edeps accordingly)
                    # Do it for all edeps and also had part edeps (because later we need had part only to evaluate FD had veto)
                    ###########################################################

                    # First tell the module where is the random thrown vertex
                    geoEff.setNDrandVertex(throwVtxX_nd[0], throwVtxY_nd[0], throwVtxZ_nd[0])
                    print ("-- nd throw x: ", throwVtxX_nd[0], "y: ", throwVtxY_nd[0], ", z: ", throwVtxZ_nd[0])
                    all_startposdep_ndorig_matrix = geoEff.move2ndorigin(ndrandthrowresultall_start.thrownEdepspos[0]) # returns Eigen::Matrix3Xf
                    all_stopposdep_ndorig_matrix = geoEff.move2ndorigin(ndrandthrowresultall_stop.thrownEdepspos[0]) # repeat for stop points
                    had_posdep_ndorig_matrix = geoEff.move2ndorigin(ndrandthrowresulthad.thrownEdepspos[0])

                    ####################################################################################################################
                    # Apply earth curvature correction to translate into FD coordinate system, vtx now at (0,0,0) in FD det coordinate,
                    # this info will be used by both leg 1 and leg 2 transformations below
                    ####################################################################################################################
                    had_posdep_fdorig_matrix = geoEff.getn2fEarthCurvatureCorr(had_posdep_ndorig_matrix, beamLineRotation) # returns Eigen::Matrix3Xf
                    all_startposdep_fdorig_matrix = geoEff.getn2fEarthCurvatureCorr(all_startposdep_ndorig_matrix, beamLineRotation)
                    all_stopposdep_fdorig_matrix = geoEff.getn2fEarthCurvatureCorr(all_stopposdep_ndorig_matrix, beamLineRotation) # repeat for stop points

                    ################
                    # Paired FD evt
                    ################
                    # non ecc nd evt + fd evt with same random rotation with non ecc nd
                    # ND event is without any earth curvature correction (already obtained above),
                    # after earth curvature correction, only random translate the nd event in fd to obtain paired FD event

                    # Tell the module where the vertex is in FD
                    geoEff.setVertexFD(0, 0, 0) # it's at FD origin because we moved it to origin and then just rotated at there

                    # Configure veto region for hadronic veto for FD throws
                    geoEff.setVetoSizes([fd_veto_region_size])
                    # Configure E threshold for hadronic veto for FD throws
                    geoEff.setVetoEnergyThresholds([fd_veto_threshold_energy])

                    tot_fd_throw_pair_nd_nonecc = 0

                    while tot_fd_throw_pair_nd_nonecc < max_fd_throws:
                        print ("---- tot fd throw to pair nd-non-ecc:", tot_fd_throw_pair_nd_nonecc)
                        ##########################################################################################
                        # Below do random throw (translate only) in FD similar to ND: only one throw in FD at a time
                        ##########################################################################################
                        geoEff.setNthrowsFD(1)
                        geoEff.throwTransformsFD() # this randomly generates new vtx position in FD FV

                        fd_vtx_x_cm_pair_nd_nonecc = geoEff.getCurrentFDThrowTranslationsX()
                        fd_vtx_y_cm_pair_nd_nonecc = geoEff.getCurrentFDThrowTranslationsY()
                        fd_vtx_z_cm_pair_nd_nonecc = geoEff.getCurrentFDThrowTranslationsZ()

                        # Check if it passes FD hadronic veto
                        geoEff.setHitSegEdeps(had_edep_list) # use the same had edep list
                        fdthrowresulthad_pair_nd_nonecc = geoEff.getFDContainment4RandomThrow(had_posdep_fdorig_matrix)

                        if (fdthrowresulthad_pair_nd_nonecc.containresult[0][0][0] != 0):
                            print ("---- tot fd throw to pair nd-non-ecc:", tot_fd_throw_pair_nd_nonecc, "passed fd had veto")
                            print ("---- throw x: ", fd_vtx_x_cm_pair_nd_nonecc[0], "y: ", fd_vtx_y_cm_pair_nd_nonecc[0], ", z: ", fd_vtx_z_cm_pair_nd_nonecc[0])
                            ###########################################################
                            # FD rand throw passes veto, write paired evt info
                            ###########################################################

                            # Now change to the full list of edeps
                            # the random thrown x/y/z should reamin the same because throw is done above already
                            geoEff.setHitSegEdeps(all_edep_list)
                            fdthrowresultall_start_pair_nd_nonecc = geoEff.getFDContainment4RandomThrow(all_startposdep_fdorig_matrix)
                            # Repeat for edepsim stop points !!!
                            geoEff.setHitSegEdeps(all_edep_list)
                            fdthrowresultall_stop_pair_nd_nonecc = geoEff.getFDContainment4RandomThrow(all_stopposdep_fdorig_matrix)

                            print ("Found paired fd-nd non ecc event")

                            nd_fd_throws_passed[0] = 1

                            #################################
                            # Unpack info and store to output
                            #################################
                            print ("Saving...")

                            # Genie info
                            GenieParts_pdg[:Genie_nParts[0]] = np.array(all_genie_pdg_list, dtype=np.int32)
                            GenieParts_p0_MeV[:Genie_nParts[0]] = np.array(all_genie_p0_list, dtype=np.float32)
                            GenieParts_p1_MeV[:Genie_nParts[0]] = np.array(all_genie_p1_list, dtype=np.float32)
                            GenieParts_p2_MeV[:Genie_nParts[0]] = np.array(all_genie_p2_list, dtype=np.float32)
                            GenieParts_p3_MeV[:Genie_nParts[0]] = np.array(all_genie_p3_list, dtype=np.float32)

                            # LArBath info
                            deps_trkID[:nEdeps[0]] = np.array(all_dep_trkID_list, dtype=np.int32)
                            deps_parentID[:nEdeps[0]] = np.array(all_dep_parentID_list, dtype=np.int32)
                            deps_pdg[:nEdeps[0]] = np.array(all_dep_pdg_list, dtype=np.int32)
                            deps_E_MeV[:nEdeps[0]] = np.array(all_edep_list, dtype=np.float32)
                            deps_start_t_us[:nEdeps[0]] = np.array(all_dep_starttime_list, dtype=np.float32)
                            deps_stop_t_us[:nEdeps[0]] = np.array(all_dep_stoptime_list, dtype=np.float32)
                            larbath_vtx_cm[0] = posx
                            larbath_vtx_cm[1] = posy
                            larbath_vtx_cm[2] = posz
                            larbath_deps_start_x_cm[:nEdeps[0]] = np.array(all_dep_startpos_list[::3], dtype=np.float32) # every 3 element: x list
                            larbath_deps_start_y_cm[:nEdeps[0]] = np.array(all_dep_startpos_list[1::3], dtype=np.float32) # y list
                            larbath_deps_start_z_cm[:nEdeps[0]] = np.array(all_dep_startpos_list[2::3], dtype=np.float32) # z list
                            larbath_deps_stop_x_cm[:nEdeps[0]] = np.array(all_dep_stoppos_list[::3], dtype=np.float32)
                            larbath_deps_stop_y_cm[:nEdeps[0]] = np.array(all_dep_stoppos_list[1::3], dtype=np.float32)
                            larbath_deps_stop_z_cm[:nEdeps[0]] = np.array(all_dep_stoppos_list[2::3], dtype=np.float32)

                            # Paired event in ND from random throw
                            nd_vtx_cm_nonecc[0] = throwVtxX_nd[0]
                            nd_vtx_cm_nonecc[1] = throwVtxY_nd[0]
                            nd_vtx_cm_nonecc[2] = throwVtxZ_nd[0]

                            # Paired fd event to nd non-ecc evt
                            fd_vtx_cm_pair_nd_nonecc[0] = fd_vtx_x_cm_pair_nd_nonecc[0]
                            fd_vtx_cm_pair_nd_nonecc[1] = fd_vtx_y_cm_pair_nd_nonecc[0]
                            fd_vtx_cm_pair_nd_nonecc[2] = fd_vtx_z_cm_pair_nd_nonecc[0]
                            fd_deps_start_x_cm_pair_nd_nonecc[:nEdeps[0]] = np.array(fdthrowresultall_start_pair_nd_nonecc.thrownEdepspos[0][0,:], dtype=np.float32)
                            fd_deps_start_y_cm_pair_nd_nonecc[:nEdeps[0]] = np.array(fdthrowresultall_start_pair_nd_nonecc.thrownEdepspos[0][1,:], dtype=np.float32)
                            fd_deps_start_z_cm_pair_nd_nonecc[:nEdeps[0]] = np.array(fdthrowresultall_start_pair_nd_nonecc.thrownEdepspos[0][2,:], dtype=np.float32)
                            fd_deps_stop_x_cm_pair_nd_nonecc[:nEdeps[0]] = np.array(fdthrowresultall_stop_pair_nd_nonecc.thrownEdepspos[0][0,:], dtype=np.float32)
                            fd_deps_stop_y_cm_pair_nd_nonecc[:nEdeps[0]] = np.array(fdthrowresultall_stop_pair_nd_nonecc.thrownEdepspos[0][1,:], dtype=np.float32)
                            fd_deps_stop_z_cm_pair_nd_nonecc[:nEdeps[0]] = np.array(fdthrowresultall_stop_pair_nd_nonecc.thrownEdepspos[0][2,:], dtype=np.float32)

                            # Break the while loop, move on to next evt
                            print ("Paired data saved, breaking fd throw loop")
                            break
                        else:
                            print ("---- tot fd throw to pair nd-non-ecc:", tot_fd_throw_pair_nd_nonecc, "failed fd had veto!")

                        # indentation is important!
                        # if don't, put it in another random FD pos...until it passes FD veto
                        tot_fd_throw_pair_nd_nonecc = tot_fd_throw_pair_nd_nonecc + 1

                    # if reached max fd throw and still didn't pass FD veto, try next nd throw
                    if tot_fd_throw_pair_nd_nonecc == max_fd_throws:
                        print ("Reached max fd throw to pair nd-non-ecc", max_fd_throws, ", continue to next nd throw")
                        tot_nd_throw = tot_nd_throw + 1
                        continue

                    # if found paired fd evts, break nd throw loop
                    print ("And breaking nd throw loop")
                    break

                else:
                    print ("-- nd throw", tot_nd_throw, "failed nd mu selection!")

            else:
                print ("-- nd throw", tot_nd_throw, "failed nd had veto!")

        else:
            print ("-- nd throw", tot_nd_throw, "outside nd FV!")

        # indentation is important!
        tot_nd_throw = tot_nd_throw + 1

    if not nd_fd_throws_passed[0]:
        print("-- no nd fd throw pairs found after max tries. Giving up!")
        nEdeps[0] = 0
        nd_lep_contained_prob_nonecc[0] = 0
        nd_lep_tracker_prob_nonecc[0] = 0
        nd_lep_ke_MeV_exit_ndlar_nonecc[0] = 0

    # event level
    myEvents.Fill()

f_out.cd()
myEvents.Write()

print("\n")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
