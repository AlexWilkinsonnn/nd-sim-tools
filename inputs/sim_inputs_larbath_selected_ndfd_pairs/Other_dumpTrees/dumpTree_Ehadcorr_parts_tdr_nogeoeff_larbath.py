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
    except Exception, e:
        print(str(e))
        print("dumpTree.py WARNING!!! Error reading GNUMIXML file. Set GNUMIXML enviornment variable to point at GNuMI flux configuration xml file. Exitting!")
        exit(-1)
    
    event = ROOT.TG4Event()
    events.SetBranchAddress("Event",ROOT.AddressOf(event))

    N = events.GetEntries()

    print "Starting loop over %d entries" % N
    ient = 0 # This is unnecessary
    iwritten = 0
    for ient in range(N):

        if ient % 100 == 0:
            print "Event %d of %d..." % (ient,N)
        events.GetEntry(ient)

        t_eventID[0] = -1;
        t_eventID[0] = event.EventId;

        # For location information
        for itraj, trajectory in enumerate(event.Trajectories):
            momentum = trajectory.InitialMomentum.Vect().Mag()
            e = trajectory.InitialMomentum.E()
            PDG = trajectory.PDGCode
            t_TrueE[0] = e
            t_PDG[0] = PDG
            
            primaryID = trajectory.TrackId
            
            StopInActive = 0
            DistToCollar = 0.
            
            points = trajectory.Points
            
            start_point = points[0]
            start_point_loc = start_point.Position.Vect()
            last_point = points[-1]
            last_point_loc = last_point.Position.Vect()
        
            # Handle the last hit
            x = last_point_loc.X()
            y = last_point_loc.Y()
            z = last_point_loc.Z()
            node = tgeo.FindNode(x, y, z)
            volName = node.GetName()
            if ("_").join(volName.split("_")[:-2]) == "volLArActive":
                last_hit_loc = 0
            elif (x < 3573.5 and x > -3573.5 and y < 1558.97 and y > -1451.23 and z < 9205.5 and z > 4114.5):
                last_hit_loc = 1
            else:
                last_hit_loc = 2
                    
            # Handle the last hit before the collar
            i_last_point_pre_collar = 0
            for ipoint, point in enumerate(points):
                point_loc = point.Position.Vect()
                x_corr = point_loc.X()/10. - offset[0]
                y_corr = point_loc.Y()/10. - offset[1]
                z_corr = point_loc.Z()/10. - offset[2]
                x = point_loc.X()
                y = point_loc.Y()
                z = point_loc.Z()
                if (x_corr < collarLo[0] or x_corr > collarHi[0] or y_corr < collarLo[1] or y_corr > collarHi[1] or z_corr < collarLo[2] or z_corr > collarHi[2]) and (x < 3573.5 and x > -3573.5 and y < 1558.97 and y > -1451.23 and z < 9205.5 and z > 4114.5):
                    i_last_point_pre_collar = ipoint - 1
                    break
            last_point_pre_collar = points[i_last_point_pre_collar]
            last_point_pre_collar_loc = last_point_pre_collar.Position.Vect()
            DistToCollar = (last_point_pre_collar_loc - start_point_loc).Mag()
                    
            # For energy information
            RecoE = 0.
            CorrRecoE = 0.
            TrueEhadveto = 0.
            RecoEhadveto = 0.
            CorrRecoEhadveto = 0.
            Corr = 0.
            CorrCollar = 0.
            
            hits = []
            inactive_hits = []
            rho_LAr = 8.73811350597e18 # Extracted from edep-sim output
            
            for key in event.SegmentDetectors:
                for hit in key.second:
                    ID = hit.Contrib[0]
                    if ID != primaryID:
                        continue
                    hMid = ROOT.TVector3(
                        (hit.Start[0] + hit.Stop[0])/2,
                        (hit.Start[1] + hit.Stop[1])/2,
                        (hit.Start[2] + hit.Stop[2])/2
                    )
                    x = hMid.X()
                    y = hMid.Y()
                    z = hMid.Z()
                    node = tgeo.FindNode(hMid.X(), hMid.Y(), hMid.Z())
                    if not node:
                        continue
                    volName = node.GetName()
                    if ("_").join(volName.split("_")[:-2]) == "volLArActive":
                        hits.append(hit)
                    else:
                        inactive_hits.append(hit)
        
            for hit in hits:
                # We already know it's in the active region inside ND 
                # LAr, but need start to determine if it's in the 
                # collar
                hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                RecoE += hit.EnergyDeposit
                # Check if hit is in collar region
                if hStart.x() < collarLo[0] or hStart.x() > collarHi[0] or hStart.y() < collarLo[1] or hStart.y() > collarHi[1] or hStart.z() < collarLo[2] or hStart.z() > collarHi[2]:
                    TrueEhadveto += hit.EnergyDeposit
                    RecoEhadveto += hit.EnergyDeposit
                    
            CorrRecoE = RecoE
            CorrRecoEhadveto = RecoEhadveto
                    
            for hit in inactive_hits:
                Edep = hit.EnergyDeposit
                hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                hMid = ROOT.TVector3(
                    (hit.Start[0] + hit.Stop[0])/2,
                    (hit.Start[1] + hit.Stop[1])/2,
                    (hit.Start[2] + hit.Stop[2])/2
                )
                node = tgeo.FindNode(hMid.X(), hMid.Y(), hMid.Z())
                rho = node.GetMedium().GetMaterial().GetDensity()
                corr = 1 - (rho / rho_LAr)
                Edep_corr = Edep * corr
                Corr += Edep_corr
                CorrRecoE += Edep_corr
                
                if hStart.x() < collarLo[0] or hStart.x() > collarHi[0] or hStart.y() < collarLo[1] or hStart.y() > collarHi[1] or hStart.z() < collarLo[2] or hStart.z() > collarHi[2]:
                    TrueEhadveto += Edep
                    CorrRecoEhadveto += Edep_corr
                    CorrCollar += Edep_corr
            
            t_RecoE[0] = RecoE * 0.001
            t_CorrRecoE[0] = CorrRecoE * 0.001
            t_TrueEhadveto[0] = TrueEhadveto
            t_RecoEhadveto[0] = RecoEhadveto
            t_CorrRecoEhadveto[0] = CorrRecoEhadveto
            t_StopInActive[0] = last_hit_loc
            t_Corr[0] = Corr * 0.001
            t_CorrCollar[0] = CorrCollar
            t_DistToCollar[0] = DistToCollar

            tout.Fill()
        

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
    t_PDG = array('i', [0])
    tout.Branch('PDG', t_PDG, 'PDG/I')
    t_TrueE = array('f', [0])
    tout.Branch('TrueE', t_TrueE, 'TrueE/F')
    t_RecoE = array('f', [0])
    tout.Branch('RecoE', t_RecoE, 'RecoE/F')
    t_CorrRecoE = array('f', [0])
    tout.Branch('CorrRecoE', t_CorrRecoE, 'CorrRecoE/F')
    t_TrueEhadveto = array('f', [0])
    tout.Branch('TrueEhadveto', t_TrueEhadveto, 'TrueEhadveto/F')
    t_RecoEhadveto = array('f', [0])
    tout.Branch('RecoEhadveto', t_RecoEhadveto, 'RecoEhadveto/F')
    t_CorrRecoEhadveto = array('f', [0])
    tout.Branch('CorrRecoEhadveto', t_CorrRecoEhadveto, 'CorrRecoEhadvet0/F')
    t_StopInActive = array('i', [0])
    tout.Branch('StopInActive', t_StopInActive, 'StopInActive/I')
    t_Corr = array('f', [0])
    tout.Branch('Corr', t_Corr, 'Corr/F')
    t_CorrCollar = array('f', [0])
    tout.Branch('CorrCollar', t_CorrCollar, 'CorrCollar/F')
    t_DistToCollar = array('f', [0])
    tout.Branch('DistToCollar', t_DistToCollar, "DistToCollar/F")

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



