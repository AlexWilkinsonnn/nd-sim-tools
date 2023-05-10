"""
Converts ROOT file created by edep-sim into a flat tree of energy deposits ready for reading into
art-root and running FD detsim
"""

import argparse, math
from array import array

import ROOT

EDEP2CM = 0.1 # convert to cm
EDEP2US = 0.001 # convert to microseconds

# To compute number of ionisation elections from energy
BIRKS_AB = 0.800 # Birks Model constant
BIRKS_KB = 0.0486 # ""
E_FIELD = 0.5 # kV/cm
LAR_DENSITY = 1.38 # g/cm^3
W_ION = 23.6e-6 # Energy expended per ion pair in LAr


def main(args):
    ROOT.gROOT.ProcessLine("#include<vector>")
    out_f = ROOT.TFile.Open(args.output_file, "RECREATE")
    out_t = ROOT.TTree("nd_depos", "nd_depos")

    depos = ROOT.vector("std::vector<double>")()
    out_t.Branch("nd_depos", depos)
    vertex = ROOT.vector("double")(4)
    out_t.Branch("vertex", vertex)
    eventID = array("i", [0])
    out_t.Branch("eventID", eventID, "eventID/I")

    in_f = ROOT.TFile(args.input_file)
    in_t = in_f.Get("EDepSimEvents")

    event = ROOT.TG4Event()
    in_t.SetBranchAddress("Event", event)
    entries = in_t.GetEntries()

    if args.muoncontained_ids_file:
        with open(args.muoncontained_ids_file, "r") as f:
            valid_evids = { int(line.rstrip()) for line in f }

    print(valid_evids)
    for i_event in range(entries):
        depos.clear()
        eventID[0] = -999
        for i in range(4):
            vertex[i] = -999

        in_t.GetEntry(i_event)

        if args.muoncontained_ids_file and int(event.EventId) not in valid_evids:
            continue

        eventID[0] = event.EventId

        if len(event.Primaries) != 1:
            raise Exception("I don't understand something :)")

        for primary_vertex in event.Primaries:
            vertex[0] = primary_vertex.GetPosition().X() / 10
            vertex[1] = primary_vertex.GetPosition().Y() / 10
            vertex[2] = primary_vertex.GetPosition().Z() / 10
            vertex[3] = 0.0 # Assume t0 is perfect

        for _, segments in event.SegmentDetectors:
            for segment in segments:
                # Copied from larnd-sim's quenching.py
                # NOTE This seems to be a simplified compared to FD IonAndScint module
                dE = segment.GetEnergyDeposit()
                xd = (segment.GetStop().X() - segment.GetStart().X()) * EDEP2CM
                yd = (segment.GetStop().Y() - segment.GetStart().Y()) * EDEP2CM
                zd = (segment.GetStop().Z() - segment.GetStart().Z()) * EDEP2CM
                dx = math.sqrt(xd**2 + yd**2 + zd**2)
                dEdx = dE / dx if dx > 0 else 0
                recomb = BIRKS_AB / (1 + BIRKS_KB * dEdx / (E_FIELD * LAR_DENSITY))
                n_electrons = recomb * dE / W_ION

                depo = ROOT.vector("double")(12)
                depo[0] = -1 # track id from greatest contributor (not using)
                depo[1] = -1 # pdg code from greatest contributor (not using)
                depo[2] = segment.GetStart().X() * EDEP2CM
                depo[3] = segment.GetStop().X() * EDEP2CM
                depo[4] = segment.GetStart().Y() * EDEP2CM
                depo[5] = segment.GetStop().Y() * EDEP2CM
                depo[6] = segment.GetStart().Z() * EDEP2CM
                depo[7] = segment.GetStop().Z() * EDEP2CM
                depo[8] = segment.GetStart().T() * EDEP2US
                depo[9] = segment.GetStop().T() * EDEP2US
                depo[10] = n_electrons
                depo[11] = dE
                depos.push_back(depo)

        out_t.Fill()

    out_f.Write()
    out_f.Close()


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("input_file", type=str)
    parser.add_argument("output_file", type=str)

    parser.add_argument(
        "--muoncontained_ids_file",
        type=str, default="",
        help="File containing eventIDs of events with lepton contained in LAr"
    )

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    main(parse_arguments())

