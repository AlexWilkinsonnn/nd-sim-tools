#!/usr/bin/env python
"""
dumpTree.py taken from last_unstructured_cafs tag of ND_CAFMaker.
Removed geometric efficiency stuff.
Uses a separate gdml file for checking which ND hall volumes segments are in. Doing this so the
reco can be applied to LArBath edep-sim events so long as the vertex was generated in the ND hall
geometry coordinates. The separate gdml file is passed through a dummy edep-sim files that has
an "EDepSimGeometry" branch. When passing the actual gdml I ran into units (cm and mm) issues.
NOTE: there are problems with doing this. The magnetic field of ND-GAr is not being simulated so
the variables that determine if the muon is tracker matched with GAr/ECAL
(distance muon travels in respective geometry) will be a bit wrong since the track will not be
bending.
"""
import sys
import os.path
import os
import ROOT
from optparse import OptionParser
import xml.etree.ElementTree as ET
from array import array
from math import cos, sin

# Just grepped the gdml for 'auxval="BarrelECal_vol"' and 'auxvalue="EndcapECal_vol"' to get these
SENSDET_BARRELECALS = set(
    "BarrelECal_stave%02d_module%02d_layer_%02d_slice2_vol" % (stave, module, layer)
    for stave in range(1, 9)
        for module in range(1,6)
            for layer in range(1, 61)
)
SENSDET_ENDCAPECALS = set(
    "EndcapECal_stave%02d_module%02d_layer_%02d_slice2_vol" % (stave, module, layer)
    for stave in range(1, 5)
        for module in [0, 6]
            for layer in range(1, 61)
)

def loop( evt, tgeo, tout ):

    offset = [ 0., 5.5, 411. ]
    collarLo = [ -320., -120., 30. ]
    collarHi = [ 320., 120., 470. ]

    # Average neutrino decay position in beam coordinates as a function of vertex x (from Luke): Will be used to set the decay position event-by-event.
    OffAxisPoints = array('f', [-2, 0.5, 3,    5.5, 8, 10.5, 13, 15.5, 18,  20.5, 23,  25.5, 28,   30.5])
    meanPDPZ = array('f', [ 93.6072, 93.362,  90.346, 85.6266, 81.1443, 76.6664, 73.0865, 69.8348, 67.5822, 65.005, 62.4821, 60.8336, 59.1433, 57.7352])
    gDecayZ = ROOT.TGraph(14, OffAxisPoints, meanPDPZ)
    
    # Get beam parameters
    try :
        GNuMIFluxTree = ET.parse(os.environ["GNUMIXML"])
        beamLineRotation = float(GNuMIFluxTree.getroot()[0].findall('beamdir')[0].findall('rotation')[0].text)
        posText = GNuMIFluxTree.getroot()[0].findall('beampos')[0].text
        posText = posText.replace("(", "")
        posText = posText.replace(")", "")
        posText = posText.replace("=", ",")
        posText = posText.split(",")
        beamRefDetCoord = [ float(posText[i]) for i in range(3)]
        detRefBeamCoord = [ float(posText[i]) for i in range(3,6)]
    except Exception as e:
        print(str(e))
        print("dumpTree.py WARNING!!! Error reading GNUMIXML file. Set GNUMIXML enviornment variable to point at GNuMI flux configuration xml file. Exitting!")
        exit(-1)
    
    event = ROOT.TG4Event()
    events.SetBranchAddress("Event",ROOT.AddressOf(event))

    N = events.GetEntries()

    print("Starting loop over %d entries" % N)
    ient = 0 # This is unnecessary
    iwritten = 0
    for ient in range(N):

        if ient % 100 == 0:
            print("Event %d of %d..." % (ient,N))
        events.GetEntry(ient)

        t_eventID[0] = -1;
        t_eventID[0] = event.EventId;

        for ivtx,vertex in enumerate(event.Primaries):

            ## initialize output variables
            t_ievt[0] = ient;
            t_vtx[0]=0.0; t_vtx[1]=0.0; t_vtx[2]=0.0;
            t_p3lep[0]=0.0; t_p3lep[1]=0.0; t_p3lep[2]=0.0;
            t_lepDeath[0]=0.0; t_lepDeath[1]=0.0; t_lepDeath[2]=0.0;
            t_lepPdg[0] = 0
            t_lepKE[0] = 0.
            t_muonExitPt[0] = 0.0; t_muonExitPt[1] = 0.0; t_muonExitPt[2] = 0.0; 
            t_muonExitMom[0] = 0.0; t_muonExitMom[1] = 0.0; t_muonExitMom[2] = 0.0; 
            t_muonReco[0] = -1;
            t_muon_endVolName.replace(0, ROOT.std.string.npos, "")
            t_muGArLen[0]=0.0;
            t_hadTot[0] = 0.
            t_hadP[0] = 0.
            t_hadN[0] = 0.
            t_hadPip[0] = 0.
            t_hadPim[0] = 0.
            t_hadPi0[0] = 0.
            t_hadOther[0] = 0.
            t_hadCollar[0] = 0.
            t_nFS[0] = 0
            ## done

            # now ID numucc
            reaction=vertex.Reaction

            # set the vertex location for output
            for i in range(3): 
                t_vtx[i] = vertex.Position[i] / 10. - offset[i] # cm

            # fiducial vertex pre-cut
            #if abs(t_vtx[0]) > 310. or abs(t_vtx[1]) > 110. or t_vtx[2] < 40. or t_vtx[2] > 360.:
            #    continue
            
            # Use vertex to determine mean decay point
            decayZbeamCoord = gDecayZ.Eval(vertex.Position[0] / 1000 - detRefBeamCoord[0])*100 # in cm
            decayZdetCoord = (-1*detRefBeamCoord[2]*100+decayZbeamCoord)*cos(beamLineRotation) - (-1*detRefBeamCoord[1]*100*sin(beamLineRotation)) + beamRefDetCoord[2]*100
            decayYdetCoord = (-1*detRefBeamCoord[1]*100*cos(beamLineRotation)) + (-1*detRefBeamCoord[2]*100+decayZbeamCoord)*sin(beamLineRotation) + beamRefDetCoord[1]*100
            decayXdetCoord = -1*detRefBeamCoord[0]*100 + beamRefDetCoord[0]*100

            ileptraj = -1
            nfsp = 0
            nHadrons = 0
            # get the lepton kinematics from the edepsim file
            fsParticleIdx = {}
            for ipart,particle in enumerate(vertex.Particles):
                e = particle.Momentum[3]
                p = (particle.Momentum[0]**2 + particle.Momentum[1]**2 + particle.Momentum[2]**2)**0.5
                m = (e**2 - p**2)**0.5
                t_fsPdg[nfsp] = particle.PDGCode
                t_fsPx[nfsp] = particle.Momentum[0]
                t_fsPy[nfsp] = particle.Momentum[1]
                t_fsPz[nfsp] = particle.Momentum[2]
                t_fsE[nfsp] = e
                fsParticleIdx[particle.TrackId] = nfsp
                nfsp += 1
                pdg = particle.PDGCode
                if abs(pdg) in [11,12,13,14]:
                    ileptraj = particle.TrackId
                    t_lepPdg[0] = pdg
                    # set the muon momentum for output
                    for i in range(3): t_p3lep[i] = particle.Momentum[i]
                    t_lepKE[0] = e - m

            assert ileptraj != -1, "There isn't a lepton??"
            t_nFS[0] = nfsp

            # If there is a muon, determine how to reconstruct its momentum and charge
            if abs(t_lepPdg[0]) == 13:
                leptraj = event.Trajectories[ileptraj]
                for p in leptraj.Points:
                    pt = p.Position
                    node = tgeo.FindNode( pt.X(), pt.Y(), pt.Z() )
                    volName = node.GetName()
                    active = False
                    if "LAr" in volName or "PixelPlane" in volName or "sPlane" in volName: # in active volume, update exit points
                        t_muonExitPt[0] = pt.X() / 10. - offset[0]
                        t_muonExitPt[1] = pt.Y() / 10. - offset[1]
                        t_muonExitPt[2] = pt.Z() / 10. - offset[2]
                        t_muonExitMom[0] = p.Momentum.x()
                        t_muonExitMom[1] = p.Momentum.y()
                        t_muonExitMom[2] = p.Momentum.z()
                    else:
                        t_muonExitPt[0] = pt.X() / 10. - offset[0]
                        t_muonExitPt[1] = pt.Y() / 10. - offset[1]
                        t_muonExitPt[2] = pt.Z() / 10. - offset[2]
                        break

                endpt = leptraj.Points[-1].Position

                node = tgeo.FindNode( endpt.X(), endpt.Y(), endpt.Z() )

                t_lepDeath[0] = endpt.X()/10. - offset[0]
                t_lepDeath[1] = endpt.Y()/10. - offset[1]
                t_lepDeath[2] = endpt.Z()/10. - offset[2]

                endVolName = node.GetName() if node else "larbath"
                t_muon_endVolName.replace(0, ROOT.std.string.npos, endVolName)
                if "LArActive" in endVolName: t_muonReco[0] = 1 # contained
                elif "ECal" in endVolName: t_muonReco[0] = 3 # ECAL stopper
                else: t_muonReco[0] = 0 # endpoint not in active material, but might still be reconstructed by curvature if GAr length > 0

                # look for muon hits in the gas TPC
                hits = []
                # cntr_volNames = {}
                for key in event.SegmentDetectors:
                    for hit in key.second:
                        hMid = ROOT.TVector3(
                            (hit.Start[0] + hit.Stop[0])/2,
                            (hit.Start[1] + hit.Stop[1])/2,
                            (hit.Start[2] + hit.Stop[2])/2,
                        )
                        node = tgeo.FindNode(hMid.X(), hMid.Y(), hMid.Z())
                        if not node:
                            # print hMid.X()/10., hMid.Y()/10., hMid.Z()/10.
                            # print ""
                            continue
                        volName = node.GetName()
                        # if volName not in cntr_volNames:
                        #     cntr_volNames[volName] = 1
                        # else:
                        #     cntr_volNames[volName] += 1
                        # vol names are like "name_PV_{0,1,2,...}" since volumes can be referenced
                        # multiple times in the gdml. Strip off the PV_{0,...} bit.
                        if ("_").join(volName.split("_")[:-2]) in ["TPC_Drift1", "TPC_Drift2"]:
                            hits.append(hit)

                # print cntr_volNames
                # print ""

                tot_length = 0.0
                for hit in hits:
                    if hit.PrimaryId == ileptraj: # hit is due to the muon
                        # TG4HitSegment::TrackLength includes all delta-rays, which spiral in gas TPC and give ridiculously long tracks
                        hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                        hStop = ROOT.TVector3( hit.Stop[0]/10.-offset[0], hit.Stop[1]/10.-offset[1], hit.Stop[2]/10.-offset[2] )
                        tot_length += (hStop-hStart).Mag()

                # look for muon hits in the ECAL
                ehits = []
                for key in event.SegmentDetectors:
                    for hit in key.second:
                        hMid = ROOT.TVector3(
                            (hit.Start[0] + hit.Stop[0])/2,
                            (hit.Start[1] + hit.Stop[1])/2,
                            (hit.Start[2] + hit.Stop[2])/2
                        )
                        node = tgeo.FindNode(hMid.X(), hMid.Y(), hMid.Z())
                        if not node:
                            continue
                        volName = node.GetName()
                        if (
                            ("_").join(volName.split("_")[:-2]) in
                            SENSDET_BARRELECALS.union(SENSDET_ENDCAPECALS)
                        ):
                            ehits.append(hit)

                etot_length = 0.0
                for hit in ehits:
                    if hit.PrimaryId == ileptraj: # hit is due to the muon
                        # TG4HitSegment::TrackLength includes all delta-rays, which spiral in gas TPC and give ridiculously long tracks
                        hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                        hStop = ROOT.TVector3( hit.Stop[0]/10.-offset[0], hit.Stop[1]/10.-offset[1], hit.Stop[2]/10.-offset[2] )
                        etot_length += (hStop-hStart).Mag()
                
                t_muGArLen[0] = tot_length
                t_muECalLen[0] = etot_length

                if tot_length > 50.: t_muonReco[0] = 2 # tracker

            # hadronic containment -- find hits in ArgonCube. 
            # Sort into active vs. inactive. 
            # Inactive hits will be used to correct the active hits to 
            # make the LArbath reco quantities closer to what would 
            # actually be seen in the ND
            hits = []
            inactive_hits = []
            rho_LAr = 8.73811350597e18 # Extracted from edep-sim output. 
                                       # Units of MeV ns^2 mm^-5
                                       
            # For each detector volume that has hits:
            for key in event.SegmentDetectors:
            	# For each hit in that volume:
                for hit in key.second:
                    # Find the midpoint and the volume it's in. 
                    # Check that it's a valid volume
                    hMid = ROOT.TVector3(
                        (hit.Start[0] + hit.Stop[0])/2,
                        (hit.Start[1] + hit.Stop[1])/2,
                        (hit.Start[2] + hit.Stop[2])/2
                    )
                    node = tgeo.FindNode(hMid.X(), hMid.Y(), hMid.Z())
                    if not node:
                        continue
                    
                    # We determine if the hit is in the active v. inactive 
                    # region using the volume name.
                    # The volume name for the active regions is 
                    # "volLarActive". Everything else is inactive.
                    volName = node.GetName()
                    if ("_").join(volName.split("_")[:-2]) == "volLArActive":
                        hits.append(hit)
                    # Correction should only include hits in ND-LAr.
                    elif (hMid.X() > -3573.5 and hMid.X() < 3573.5 and 
                        hMid.Y() > -1451.23 and hMid.Y() < 1558.97 and 
                        hMid.Z() > 4114.5 and hMid.Z() < 9205.5):
                        inactive_hits.append(hit)

            # Truth-matching energy
            # We want to associate all of the energy from simulated 
            # particles in Geant (trajectories) to the ultimate parent 
            # particle, that is, the primary particle from whence that 
            # energy came.
            
            # Make dictionaries of trajectory --> primary pdg, active 
            # energy, inactive energy
            traj_to_pdg = {}
            traj_to_Eact = {}
            traj_to_Einact = {}
            # For pi0s, stop at the photons to determine the visible 
            # energy due to gamma1/gamma2. To do this, make a dictionary 
            # where the keys are the track IDs of particles that came from 
            # a pi0 gamma, and the values are the trackIDs of the gammas 
            # that correspond to the trackID key. The list is to record 
            # the trackIDs of the gammas forming from a primary pi0.
            tid_to_gamma = {}
            gamma_tids = []
            # For each trajectory (aka particle in Geant):
            for traj in event.Trajectories:
                # Mom is the track ID for the trajectory that came right 
                # before the trajectory in question. If it's a primary 
                # particle, the value is -1. 
                mom = traj.ParentId
                # Each particle has a unique track ID. 
                tid = traj.TrackId
                # The first n trajectories are the initial trajectories 
                # for the n primary particles in an event. The next line 
                # therefore considers if the trajectories in question are 
                # gammas from a primary pi0.
                if event.Trajectories[mom].PDGCode == 111 and event.Trajectories[tid].PDGCode == 22 and event.Trajectories[mom].ParentId == -1:
                    # If yes, then append the trajectory's trackID to the 
                    # list of gamma trackIDs.
                    gamma_tids.append(tid)
                # If we are not looking at a primary, we want to know the 
                # primary PDG that led to that trajectory. This block 
                # walks back through the chain of trajectories to get that 
                # info. If the trajectory in question is not a primary:
                while mom != -1:
                    # The next trackId to consider is the mother of the 
                    # current one
                    tid = mom
                    # The next mother to consider is the parentID of the 
                    # current mother
                    mom = event.Trajectories[mom].ParentId
                    # If the next mother to consider is one of the gammas 
                    # from a primary pi0:
                    if mom in gamma_tids:
                        # This is a particle that is directly due to a 
                        # gamma from a pi0, so add it to the dictionary.
                        tid_to_gamma[tid] = mom
                # We now have a traj, and the trackID (tid) of its 
                # ultimate mother. Record the PDG in the dictionary and 
                # set the energy counter for that trajectory to 0.0
                traj_to_pdg[traj] = event.Trajectories[tid].PDGCode
                traj_to_Eact[traj] = 0.0
                traj_to_Einact[traj] = 0.0

            # Initialize collar energy, as it will be handled like a 
            # trajectory
            collar_energy = 0.

            # Each of these lists has one entry per primary particle
            track_length = [0. for i in range(nfsp)]
            dEdX = [[] for i in range(nfsp)]
            this_step = [[0.,0.] for i in range(nfsp)]
            end_point = [None for i in range(nfsp)]
            int_energy = [0. for i in range(nfsp)]
            trk_calo = [0. for i in range(nfsp)]
            
            # The dictionary gamma energy is for keeping the energy 
            # associated with each gamma from a primary pi0 separate from 
            # each other.
            gamma_energy = {}
            for g in gamma_tids:
                gamma_energy[g] = 0.
                
            # Loop over the active hits first and assign them to their 
            # appropriate places.
            for hit in hits:
                hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                hStop = ROOT.TVector3( hit.Stop[0]/10.-offset[0], hit.Stop[1]/10.-offset[1], hit.Stop[2]/10.-offset[2] )

                # Get the trajectory that made this hit.
                # Don't use edep-sim's PrimaryId, which thinks you want to 
                # associate absoltely everything with the primary. 
                # Instead, get the actual contributors (usually only one; 
                # may be several for certain run settings, but this is a 
                # very unusual situation) and take the biggest, which is 
                # always the zeroth entry.
                tid = hit.Contrib[0]
                traj = event.Trajectories[tid]
                
                # If it's a hit from a primary particle:
                if traj.ParentId == -1:
                    # Get the index of the primary particle that 
                    # ultimately led to this hit. fsParticleIdx is a 
                    # dictionary where the keys are the trackIDs of the 
                    # primary trajectories, and the values are indices 0, 
                    # 1, ...
                    idx = fsParticleIdx[hit.PrimaryId]
                    
                    # For this hit due to a primary, record its energy 
                    # deposited, and its length for dE/dx information
                    trk_calo[idx] += hit.EnergyDeposit
                    end_point[idx] = hStop
                    dx = (hStop-hStart).Mag()
                    track_length[idx] += dx
                    this_step[idx][1] += dx
                    this_step[idx][0] += hit.EnergyDeposit
                    # If the hit is longer than half a millimeter, then it 
                    # contributed significantly to the track:
                    if this_step[idx][1] > 0.5:
                        # Append the info we need for track reco
                        dEdX[idx].append( (this_step[idx][0], this_step[idx][1], hStart) ) # MeV/mm
                        # Reset this step
                        this_step[idx] = [0., 0.]
                # If this hit is not due to a primary particle, but its 
                # start is close to the current end of any primary track, 
                # then add its energy to the internal energy of that 
                # primary particle, as it is likely the decay product of 
                # the primary.
                else:
                    # For each primary particle's current track endpoint:
                    for k,ep in enumerate(end_point):
                        # If we haven't yet updated the endpoint, ignore 
                        # the hit
                        if ep is None: continue
                        # If it's within 10 mm of an endpoint, then 
                        # associate its energy with the internal energy of 
                        # the particle
                        if (hStart-ep).Mag() < 10.:
                            int_energy[k] += hit.EnergyDeposit

                # If the hit is due to a trajectory that ultimately came 
                # from a gamma from a pi0 decay:
                if tid in tid_to_gamma:
                    # Add this hit's energy to that gamma's index
                    gamma_energy[tid_to_gamma[tid]] += hit.EnergyDeposit

                # If the hit is not due to the lepton, but some other 
                # primary particle (here we do want to check lepton 
                # association and only consider those hits not due to the 
                # lepton)
                if hit.PrimaryId != ileptraj:
                    hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )

                    # check if hit is in collar region
                    if hStart.x() < collarLo[0] or hStart.x() > collarHi[0] or hStart.y() < collarLo[1]  or hStart.y() > collarHi[1] or hStart.z() < collarLo[2] or hStart.z() > collarHi[2]:
                        collar_energy += hit.EnergyDeposit
                    
                    # Add the energy to the trajectory's active energy 
                    # register
                    traj_to_Eact[traj] += hit.EnergyDeposit
            
            # Repeat but for the inactive hits
            # For the hit-based correction, we track the corrections for 
            # each trajectory separately. When all have been calculated, 
            # we add them to the uncorrected values to get the corrected 
            # ones. If the result is negative, the correction was too 
            # large, typically indicating that the majority of a 
            # particle's energy deposits were in an inactive region. In 
            # this case, we'd be better off not correcting at all. So, if 
            # the result is negative, we discard the correction. We sum 
            # over the corrected trajectories to get the reconstructed, 
            # corrected hadronic energy, proton energy, neutron energy, 
            # etc.
            # The collar is treated like its own trajectory
            collar_energy_corr = 0.0
                    
            # Account for hits in inactive regions. For each of these hits
            for hit in inactive_hits:
                Edep = hit.EnergyDeposit
                hMid = ROOT.TVector3(
                    (hit.Start[0] + hit.Stop[0])/2,
                    (hit.Start[1] + hit.Stop[1])/2,
                    (hit.Start[2] + hit.Stop[2])/2
                )
                # Get the GDML volume it was in
                node = tgeo.FindNode(hMid.X(), hMid.Y(), hMid.Z())
                # Get the density
                rho = node.GetMedium().GetMaterial().GetDensity()
                # Calculate the correction
                corr = 1 - (rho / rho_LAr)
                Edep_corr = Edep * corr
                
                # The hit starts and stops are translated into a 
                # coordinate system in which the center of the upstream 
                # detector wall is (0, 0, 0). This makes the veto region 
                # easier to calculate, as it's the outer 30cm of ND-LAr.
                hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                hStop = ROOT.TVector3( hit.Stop[0]/10.-offset[0], hit.Stop[1]/10.-offset[1], hit.Stop[2]/10.-offset[2] )

                # Get the trajectory that made this hit.
                # Don't use edep-sim's PrimaryId, which thinks you want to 
                # associate absoltely everything with the primary. 
                # Instead, get the actual contributors (usually only one; 
                # may be several for certain run settings, but this is a 
                # very unusual situation) and take the biggest, which is 
                # always the zeroth entry.
                tid = hit.Contrib[0]
                traj = event.Trajectories[tid]
                
                # If it's a hit from a primary particle:
                if traj.ParentId == -1:
                    # Get the index of the primary particle that 
                    # ultimately led to this hit. fsParticleIdx is a 
                    # dictionary where the keys are the trackIDs of the 
                    # primary trajectories, and the values are indices 0, 
                    # 1, ...
                    idx = fsParticleIdx[hit.PrimaryId]
                    
                    # For this hit due to a primary, record its energy 
                    # deposited, and its length for dE/dx information
                    trk_calo[idx] += Edep_corr
                    end_point[idx] = hStop
                    dx = (hStop-hStart).Mag()
                    track_length[idx] += dx
                    this_step[idx][1] += dx
                    this_step[idx][0] += Edep_corr
                    # If the hit is longer than half a millimeter, then it 
                    # contributed significantly to the track:
                    if this_step[idx][1] > 0.5:
                        # Append the info we need for track reco
                        dEdX[idx].append( (this_step[idx][0], this_step[idx][1], hStart) ) # MeV/mm
                        # Reset this step
                        this_step[idx] = [0., 0.]
                # If this hit is not due to a primary particle, but its 
                # start is close to the current end of any primary track, 
                # then add its energy to the internal energy of that 
                # primary particle, as it is likely the decay product of 
                # the primary.
                else:
                    # For each primary particle's current track endpoint:
                    for k,ep in enumerate(end_point):
                        # If we haven't yet updated the endpoint, ignore 
                        # the hit
                        if ep is None: continue
                        # If it's within 10 mm of an endpoint, then 
                        # associate its energy with the internal energy of 
                        # the particle
                        if (hStart-ep).Mag() < 10.:
                            int_energy[k] += Edep_corr

                # If the hit is due to a trajectory that ultimately came 
                # from a gamma from a pi0 decay:
                if tid in tid_to_gamma:
                    # Add this hit's energy to that gamma's index
                    gamma_energy[tid_to_gamma[tid]] += Edep_corr

                # If the hit is not due to the lepton, but some other 
                # primary particle (here we do want to check lepton 
                # association and only consider those hits not due to the 
                # lepton)
                if hit.PrimaryId != ileptraj:
                
                    # check if hit is in collar region
                    if hStart.x() < collarLo[0] or hStart.x() > collarHi[0] or hStart.y() < collarLo[1] or hStart.y() > collarHi[1] or hStart.z() < collarLo[2] or hStart.z() > collarHi[2]:
                        collar_energy_corr += Edep_corr
                    
                    # Determine primary particle
                    traj_to_Einact[traj] += Edep_corr

            # Perform the check to make sure none of the corrections 
            # result in negative energies for the trajectories.
            # This chunk takes the 'HitCorr' approach, where only the 
            # reported quantities (total had. energy of protons, total 
            # had energy of pions, etc.) are checked for negative 
            # values, and if one is found, no density correction is 
            # applied.
            p_act = 0.
            n_act = 0.
            pip_act = 0.
            pim_act = 0.
            pi0_act = 0.
            other_act = 0.
            p_inact = 0.
            n_inact = 0.
            pip_inact = 0.
            pim_inact = 0.
            pi0_inact = 0.
            other_inact = 0.
            for traj in traj_to_pdg.keys():
                pdg = traj_to_pdg[traj]
                if pdg == 2212:
                    p_act += traj_to_Eact[traj]
                    p_inact += traj_to_Einact[traj]
                elif pdg == 2112:
                    n_act += traj_to_Eact[traj]
                    n_inact += traj_to_Einact[traj]
                elif pdg == 211:
                    pip_act += traj_to_Eact[traj]
                    pip_inact += traj_to_Einact[traj]
                elif pdg == -211:
                    pim_act += traj_to_Eact[traj]
                    pim_inact += traj_to_Einact[traj]
                elif pdg == 111:
                    pi0_act += traj_to_Eact[traj]
                    pi0_inact += traj_to_Einact[traj]
                else:
                    other_act += traj_to_Eact[traj]
                    other_inact += traj_to_Einact[traj]
            if p_act + p_inact > 0: t_hadP[0] = p_act + p_inact
            else: t_hadP[0] = p_act
            if n_act + n_inact > 0: t_hadN[0] = n_act + n_inact
            else: t_hadN[0] = n_act
            if pip_act + pip_inact > 0: t_hadPip[0] = pip_act + pip_inact
            else: t_hadPip[0] = pip_act
            if pim_act + pim_inact > 0: t_hadPim[0] = pim_act + pim_inact
            else: t_hadPim[0] = pim_act
            if pi0_act + pi0_inact > 0: t_hadPi0[0] = pi0_act + pi0_inact
            else: t_hadPi0[0] = pi0_act
            if other_act + other_inact > 0: t_hadOther[0] = other_act + other_inact
            else: t_hadOther[0] = other_act
            
            t_hadCollar[0] = (collar_energy + collar_energy_corr if collar_energy + collar_energy_corr > 0 else collar_energy)
            
            # The following chunk implements the 'HybridCorr', where 
            # each trajectory is checked for a negative corrected energy,
            # and if it has one, the correction is not applied for 
            # that trajectory.
            ''' traj_to_Ereco = {}
            for traj in traj_to_pdg.keys():
                traj_E = traj_to_Eact[traj] + traj_to_Einact[traj]
                if traj_E < 0:
            	    traj_to_Ereco[traj] = traj_to_Eact[traj]
                else:
                    traj_to_Ereco[traj] = traj_E
                    t_totCorr[0] += traj_to_Einact[traj]
            t_hadCollar[0] = (collar_energy + collar_energy_corr if collar_energy + collar_energy_corr > 0 else collar_energy)
            
            # Assign the trajectories' energies to the appropriate places
            for traj in traj_to_pdg.keys():
                traj_E = traj_to_Ereco[traj]
                traj_pdg = traj_to_pdg[traj]
                t_hadTot[0] += traj_E
                if traj_pdg in [11, -11, 13, -13]: continue # lepton
                elif traj_pdg == 2212: t_hadP[0] += traj_E
                elif traj_pdg == 2112: t_hadN[0] += traj_E
                elif traj_pdg == 211: t_hadPip[0] += traj_E
                elif traj_pdg == -211: t_hadPim[0] += traj_E
                elif traj_pdg == 111: t_hadPi0[0] += traj_E
                else: t_hadOther[0] += traj_E '''

            for i in range(nfsp):
                t_fsTrkLen[i] = track_length[i]
                # Average of anything in the last 3cm of track
                t_fsTrkEnddEdX[i] = 0.
                t_fsTrkFrontdEdX[i] = 0.
                totlen = 0.
                j = len(dEdX[i])-1
                while j >= 0:
                    totlen += dEdX[i][j][1]
                    t_fsTrkEnddEdX[i] += dEdX[i][j][0]
                    j -= 1
                    if totlen > 3.: break
                if totlen > 0.: t_fsTrkEnddEdX[i] /= totlen
              
                totlen = 0.
                for k in range(j+1):
                    totlen += dEdX[i][k][1]
                    t_fsTrkFrontdEdX[i] += dEdX[i][k][0]
                if totlen > 0.: t_fsTrkFrontdEdX[i] /= totlen

                t_fsTrkEndpointBall[i] = int_energy[i]
                t_fsTrkCalo[i] = trk_calo[i]

                t_fsGamma1[i] = 0.
                t_fsGamma2[i] = 0.

            # Photons
            for t in gamma_tids:
                mom = event.Trajectories[t].ParentId
                if t_fsGamma1[mom] == 0.: t_fsGamma1[mom] = gamma_energy[t]
                elif t_fsGamma2[mom] == 0.: t_fsGamma2[mom] = gamma_energy[t]
                else:
                    print("Pi0 has more than two photons wtf")


            tout.Fill()
            iwritten += 1
        ient += 1 # This is unnecessary
        

if __name__ == "__main__":

    ROOT.gROOT.SetBatch(1)

    parser = OptionParser()
    parser.add_option('--infile_edepsim', help='Input edep-sim file name', default="edep.root")
    parser.add_option(
        '--edepsim_geometry',
        help=(
            'dummy edep-sim file that contrains on a branch' +
            'the gdml to use for checking which volumes energy depositions are in'
        ),
        default="geometry.gdml"
    )
    parser.add_option('--outfile', help='Output file name', default="out.root")
    parser.add_option('--seed', help='Seed for geometric efficiency throws', default=0, type = "int")

    (args, dummy) = parser.parse_args()

    # make an output ntuple
    fout = ROOT.TFile( args.outfile, "RECREATE" )
    tout = ROOT.TTree( "tree","tree" )
    # include edep-sim eventID for matching translation results
    t_eventID = array('i',[0])
    tout.Branch('eventID', t_eventID, 'eventID/I')
    t_ievt = array('i',[0])
    tout.Branch('ievt',t_ievt,'ievt/I')
    t_Ev = array('f', [0.])
    tout.Branch('Ev',t_Ev,'Ev/F')
    t_p3lep = array('f',3*[0.0])
    tout.Branch('p3lep',t_p3lep,'p3lep[3]/F')
    t_vtx = array('f',3*[0.0])
    tout.Branch('vtx',t_vtx,'vtx[3]/F')
    t_lepDeath = array('f',3*[0.0])
    tout.Branch('lepDeath',t_lepDeath,'lepDeath[3]/F')
    t_lepPdg = array('i',[0])
    tout.Branch('lepPdg',t_lepPdg,'lepPdg/I')
    t_lepKE = array('f',[0])
    tout.Branch('lepKE',t_lepKE,'lepKE/F')
    t_muonExitPt = array('f',3*[0.0])
    tout.Branch('muonExitPt',t_muonExitPt,'muonExitPt[3]/F')
    t_muonExitMom = array('f',3*[0.0])
    tout.Branch('muonExitMom',t_muonExitMom,'muonExitMom[3]/F')
    t_muonReco = array('i',[0])
    tout.Branch('muonReco',t_muonReco,'muonReco/I')
    t_muon_endVolName = ROOT.std.string()
    tout.Branch('muon_endVolName', t_muon_endVolName)
    t_muGArLen = array('f',[0])
    tout.Branch('muGArLen',t_muGArLen,'muGArLen/F')
    t_muECalLen = array('f',[0])
    tout.Branch('muECalLen',t_muECalLen,'muECalLen/F')
    t_hadTot = array('f', [0.] )
    tout.Branch('hadTot', t_hadTot, 'hadTot/F' )
    t_hadP = array('f', [0.] )
    tout.Branch('hadP', t_hadP, 'hadP/F' )
    t_hadN = array('f', [0.] )
    tout.Branch('hadN', t_hadN, 'hadN/F' )
    t_hadPip = array('f', [0.] )
    tout.Branch('hadPip', t_hadPip, 'hadPip/F' )
    t_hadPim = array('f', [0.] )
    tout.Branch('hadPim', t_hadPim, 'hadPim/F' )
    t_hadPi0 = array('f', [0.] )
    tout.Branch('hadPi0', t_hadPi0, 'hadPi0/F' )
    t_hadOther = array('f', [0.] )
    tout.Branch('hadOther', t_hadOther, 'hadOther/F' )
    t_hadCollar = array('f', [0.] )
    tout.Branch('hadCollar', t_hadCollar, 'hadCollar/F' )
    t_totCorr = array('f', [0.] )
    tout.Branch('totCorr', t_totCorr, 'totCorr/F' )
    t_nFS = array('i',[0])
    tout.Branch('nFS',t_nFS,'nFS/I')
    t_fsPdg = array('i',100*[0])
    tout.Branch('fsPdg',t_fsPdg,'fsPdg[nFS]/I')
    t_fsPx = array('f',100*[0.])
    tout.Branch('fsPx',t_fsPx,'fsPx[nFS]/F')
    t_fsPy = array('f',100*[0.])
    tout.Branch('fsPy',t_fsPy,'fsPy[nFS]/F')
    t_fsPz = array('f',100*[0.])
    tout.Branch('fsPz',t_fsPz,'fsPz[nFS]/F')
    t_fsE = array('f',100*[0.])
    tout.Branch('fsE',t_fsE,'fsE[nFS]/F')
    t_fsTrkLen = array('f',100*[0.])
    tout.Branch('fsTrkLen',t_fsTrkLen,'fsTrkLen[nFS]/F')
    t_fsTrkFrontdEdX = array('f',100*[0.])
    tout.Branch('fsTrkFrontdEdX',t_fsTrkFrontdEdX,'fsTrkFrontdEdX[nFS]/F')
    t_fsTrkEnddEdX = array('f',100*[0.])
    tout.Branch('fsTrkEnddEdX',t_fsTrkEnddEdX,'fsTrkEnddEdX[nFS]/F')
    t_fsTrkEndpointBall = array('f', 100*[0.])
    tout.Branch('fsTrkEndpointBall',t_fsTrkEndpointBall,'fsTrkEndpointBall[nFS]/F')
    t_fsGamma1 = array('f',100*[0.])
    tout.Branch('fsGamma1',t_fsGamma1,'fsGamma1[nFS]/F')
    t_fsGamma2 = array('f',100*[0.])
    tout.Branch('fsGamma2',t_fsGamma2,'fsGamma2[nFS]/F')
    t_fsTrkCalo = array('f',100*[0.])
    tout.Branch('fsTrkCalo',t_fsTrkCalo,'fsTrkCalo[nFS]/F')

    # Geometric efficiency stuff
    t_geoEffThrowResults = ROOT.std.vector('std::vector< std::vector < uint64_t > >')()
    tout.Branch('geoEffThrowResults', t_geoEffThrowResults)

    # Separate TTree to store translations and rotations of throws
    tGeoEfficiencyThrowsOut = ROOT.TTree( "geoEffThrows","geoEffThrows")

    events = ROOT.TChain( "EDepSimEvents", "main event tree" )
    #dspt = ROOT.TChain( "DetSimPassThru/gRooTracker", "other thing" )

    tf = ROOT.TFile( args.infile_edepsim )
    tf.MakeProject("EDepSimEvents","*","RECREATE++")

    events = tf.Get( "EDepSimEvents" )
    tf_dummy = ROOT.TFile(args.edepsim_geometry)
    tgeo = tf_dummy.Get("EDepSimGeometry")
    loop( events, tgeo, tout )

    fout.cd()
    tout.Write()
    tGeoEfficiencyThrowsOut.Write()
    



