#!/usr/bin/env python
"""
dumpTree_hits.py adapted from dumpTree.py to save a ROOT file where each line is a hit from a hadron and only length, edep, IsActive, and IsLAr are kept. 
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
    ient = 0
    for ient in range(N):

        if ient % 100 == 0:
            print "Event %d of %d..." % (ient,N)
        events.GetEntry(ient)

        for ivtx,vertex in enumerate(event.Primaries):
            particles = vertex.Particles
            n_particles = len(particles)
            Etot = 0.0
            for ipart, particle in enumerate(particles):
                momentum = particle.Momentum
                e = momentum.E()
                Etot += e

            ## initialize output variables
            ihit = 0

            for key in event.SegmentDetectors:
                for hit in key.second:
                    E = hit.EnergyDeposit
                    hMid = ROOT.TVector3(
                        (hit.Start[0] + hit.Stop[0])/2,
                        (hit.Start[1] + hit.Stop[1])/2,
                        (hit.Start[2] + hit.Stop[2])/2
                    )
                    x = hMid.X()
                    y = hMid.Y()
                    z = hMid.Z()
                    '''Don't fill if the hit is outside ND LAr. These 
                    coordinates were found by digging through the GDML 
                    file by hand.'''
                    if (x < 3573.5 and x > -3573.5 and y < 1558.97 and y > -1451.23 and z < 9205.5 and z > 4114.5 and E >= 0):
                        t_eventID[0] = ient
                        t_n_part[0] = n_particles
                        t_TrueE[0] = Etot
                        t_hitID[0] = ihit
                	
                        # Length
                        leng = ((hit.Stop[0] - hit.Start[0]) ** 2 + 
                            (hit.Stop[1] - hit.Start[1]) ** 2 + 
                            (hit.Stop[2] - hit.Start[2]) ** 2) ** 0.5
                        if leng > 1000:
                            continue
                        t_len[0] = leng
                	
                        t_edep[0] = hit.EnergyDeposit
                        
                        # IsActive
                        node = tgeo.FindNode(x, y, z)
                        if not node:
                            continue
                        volName = node.GetName()
                        isactive = 0
                        if ("_").join(volName.split("_")[:-2]) == "volLArActive":
                            isactive = 1
                        t_IsActive[0] = isactive
                        
                        # IsLArInND
                        IsLArInND = 0
                        material = node.GetMedium().GetMaterial().GetName()
                        if material == "LAr":
                            IsLArInND = 1
                        t_IsLArInND[0] = IsLArInND

                        # Fill
                        tout.Fill()
                        ihit += 1
        ient += 1        


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
    tout.Branch('eventID',t_eventID,'eventID/I')
    t_n_part = array('i',[0])
    tout.Branch('n_part',t_n_part,'n_part/I')
    t_TrueE = array('f',[0])
    tout.Branch('TrueE', t_TrueE, 'TrueE/F')
    t_hitID = array('i', [0])
    tout.Branch('HitID', t_hitID, 'HitID/I')
    t_len = array('f',[0])
    tout.Branch('len',t_len,'len/F')
    t_edep = array('f',[0])
    tout.Branch('edep',t_edep,'edep/F')
    t_IsActive = array('i',[0])
    tout.Branch('IsActive',t_IsActive,'IsActive/I')
    t_IsLArInND = array('i',[0])
    tout.Branch('IsLArInND',t_IsLArInND,'IsLArInND/I')

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
    
    
