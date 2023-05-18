"""
Merge CAF files with nd param reco + truth information with recodumps of FD reconstruction on the
same events.
Expects a directory of ND CAFs names FHC.{number}.caf.root and a directory of fd reco dumps with
names FHC.{number}.LoadedDepos.detsim_reco_recodump.root. Within each file events have an eventID
data product, the ND CAFs can have events without a pair.
"""

import argparse, os, subprocess, re, warnings
from array import array

import ROOT


def main(args):
    fd_files = subprocess.check_output(
        "pnfs2xrootd {}".format(os.path.join(args.fd_dir, "*")), shell=True
    )
    fd_files = fd_files.decode().strip("\n").split(" ")

    f_out = ROOT.TFile.Open(
        args.output_path if args.output_path else "nd_fd_reco.root", "RECREATE"
    )
    t_out = ROOT.TTree("nd_fd_reco", "nd_fd_reco")

    run_num = array("i", [0])
    t_out.Branch("run", run_num, "run/I")
    event_num = array("i", [0])
    t_out.Branch("eventID", event_num, "eventID/I")

    # Truth and param reco from ND CAF
    is_cc = array("i", [0])
    t_out.Branch("isCC", is_cc, "isCC/I")
    nu_pdg = array("i", [0])
    t_out.Branch("nuPDG", nu_pdg, "nuPDG/I")
    nu_mom_x = array("f", [0])
    t_out.Branch("NuMomX", nu_mom_x, "nuMomX/F")
    nu_mom_y = array("f", [0])
    t_out.Branch("NuMomY", nu_mom_y, "NuMomY/F")
    nu_mom_z = array("f", [0])
    t_out.Branch("NuMomZ", nu_mom_z, "NuMomZ/F")
    nu_E = array("f", [0])
    t_out.Branch("Ev", nu_E, "Ev/F")
    mode = array("i", [0])
    t_out.Branch("mode", mode, "mode/I")
    lep_pdg = array("i", [0])
    t_out.Branch("LepPDG", lep_pdg, "LepPDG/I")
    lep_mom_x = array("f", [0])
    t_out.Branch("LepMomX", lep_mom_x, "LepMomX/F")
    lep_mom_y = array("f", [0])
    t_out.Branch("LepMomY", lep_mom_y, "LepMomY/F")
    lep_mom_z = array("f", [0])
    t_out.Branch("LepMomZ", lep_mom_z, "LepMomZ/F")
    lep_E = array("f", [0])
    t_out.Branch("LepE", lep_E, "lepE/F")
    lep_nu_angle = array("f", [0])
    t_out.Branch("LepNuAngle", lep_nu_angle, "LepNuAngle/F")
    num_prtn = array("i", [0])
    t_out.Branch("nP", num_prtn, "nP/I")
    num_ntrn = array("i", [0])
    t_out.Branch("nN", num_ntrn, "nN/I")
    num_pip = array("i", [0])
    t_out.Branch("nPip", num_pip, "nPip/I")
    num_pim = array("i", [0])
    t_out.Branch("nPim", num_pim, "nPim/I")
    num_pi0 = array("i", [0])
    t_out.Branch("nPi0", num_pi0, "nPi0/I")
    num_kaonp = array("i", [0])
    t_out.Branch("nKp", num_kaonp, "nKp/I")
    num_kaonm = array("i", [0])
    t_out.Branch("nKm", num_kaonm, "nKm/I")
    num_kaon0 = array("i", [0])
    t_out.Branch("nK0", num_kaon0, "nK0/I")
    num_e = array("i", [0])
    t_out.Branch("ne", num_e, "ne/I")
    num_other = array("i", [0])
    t_out.Branch("nOther", num_other, "nOther/I")
    num_nuc = array("i", [0])
    t_out.Branch("nNucleus", num_nuc, "nNucleus/I")
    num_unkwn = array("i", [0])
    t_out.Branch("nUNKOWN", num_unkwn, "nUNKOWN/I")
    prtn_E = array("f", [0])
    t_out.Branch("eP", prtn_E, "eP/F")
    ntrn_E = array("f", [0])
    t_out.Branch("eN", ntrn_E, "eN/F")
    pip_E = array("f", [0])
    t_out.Branch("ePip", pip_E, "ePip/F")
    pim_E = array("f", [0])
    t_out.Branch("ePim", pim_E, "ePim/F")
    pi0_E = array("f", [0])
    t_out.Branch("ePi0", pi0_E, "ePi0/F")
    other_E = array("f", [0])
    t_out.Branch("eOther", other_E, "eOther/F")
    prtn_E_reco = array("f", [0])
    t_out.Branch("eRecoP", prtn_E_reco, "eRecoP/F")
    ntrn_E_reco = array("f", [0])
    t_out.Branch("eRecoN", ntrn_E_reco, "eRecoN/F")
    pip_E_reco = array("f", [0])
    t_out.Branch("eRecoPip", pip_E_reco, "eRecoPip/F")
    pim_E_reco = array("f", [0])
    t_out.Branch("eRecoPim", pim_E_reco, "eRecoPim/F")
    pi0_E_reco = array("f", [0])
    t_out.Branch("eRecoPi0", pi0_E_reco, "eRecoPi0/F")
    other_E_reco = array("f", [0])
    t_out.Branch("eRecoOther", other_E_reco, "eRecoOther/F")
    vtx_x = array("f", [0])
    t_out.Branch("vtxX", vtx_x, "vtxX/F")
    vtx_y = array("f", [0])
    t_out.Branch("vtxY", vtx_y, "vtxY/F")
    vtx_z = array("f", [0])
    t_out.Branch("vtxZ", vtx_z, "vtxZ/F")
    nu_E_reco = array("f", [0])
    t_out.Branch("Ev_reco", nu_E_reco, "Ev_reco/F")
    lep_E_reco = array("f", [0])
    t_out.Branch("Elep_reco", lep_E_reco, "Elep_reoo/F")
    theta_reco = array("f", [0])
    t_out.Branch("theta_reco", theta_reco, "theta_reco/F")
    reco_numu = array("i", [0])
    t_out.Branch("reco_numu", reco_numu, "reco_numu/I")
    reco_nue = array("i", [0])
    t_out.Branch("reco_nue", reco_nue, "reco_nue/I")
    reco_nc = array("i", [0])
    t_out.Branch("reco_nc", reco_nc, "reco_nc/I")
    reco_q = array("i", [0])
    t_out.Branch("reco_q", reco_q, "reco_q/I")
    muon_contained = array("i", [0])
    t_out.Branch("muon_contained", muon_contained, "muon_contained/I")
    muon_tracker = array("i", [0])
    t_out.Branch("muon_tracker", muon_tracker, "muon_tracker/I")
    muon_ecal = array("i", [0])
    t_out.Branch("muon_ecal", muon_ecal, "muon_ecal/I")
    muon_exit = array("i", [0])
    t_out.Branch("muon_exit", muon_exit, "muon_exit/I")
    lep_pdg_reco = array("i", [0])
    t_out.Branch("reco_lepton_pdg", lep_pdg_reco, "reco_lepton_pdg/I")
    muon_endpnt_x = array("f", [0])
    t_out.Branch("muon_endpntX", muon_endpnt_x, "muon_endpntX/F")
    muon_endpnt_y = array("f", [0])
    t_out.Branch("muon_endpntY", muon_endpnt_y, "muon_endpntY/F")
    muon_endpnt_z = array("f", [0])
    t_out.Branch("muon_endpntZ", muon_endpnt_z, "muon_endpntZ/F")
    hadE_veto = array("f", [0])
    t_out.Branch("Ehad_veto", hadE_veto, "Ehad_veto/F")

    # FD Reco
    numu_score = array("f", [0])
    t_out.Branch("numu_score", numu_score, "fd_numu_score/F")
    nue_score = array("f", [0])
    t_out.Branch("nue_score", nue_score, "fd_nue_score/F")
    nc_score = array("f", [0])
    t_out.Branch("nc_score", nc_score, "fd_nc_score/F")
    nutau_score = array("f", [0])
    t_out.Branch("nutau_score", nutau_score, "fd_nutau_score/F")
    antinu_score = array("f", [0])
    t_out.Branch("antinu_score", antinu_score, "fd_antinu_score/F")
    proton0_score = array("f", [0])
    t_out.Branch("proton0_score", proton0_score, "fd_proton0_score/F")
    proton1_score = array("f", [0])
    t_out.Branch("proton1_score", proton1_score, "fd_proton1_score/F")
    proton2_score = array("f", [0])
    t_out.Branch("proton2_score", proton2_score, "fd_proton2_score/F")
    protonN_score = array("f", [0])
    t_out.Branch("protonN_score", protonN_score, "fd_protonN_score/F")
    pion0_score = array("f", [0])
    t_out.Branch("pion0_score", pion0_score, "fd_pion0_score/F")
    pion1_score = array("f", [0])
    t_out.Branch("pion1_score", pion1_score, "fd_pion1_score/F")
    pion2_score = array("f", [0])
    t_out.Branch("pion2_score", pion2_score, "fd_pion2_score/F")
    pionN_score = array("f", [0])
    t_out.Branch("pionN_score", pionN_score, "fd_pionN_score/F")
    pionzero0_score = array("f", [0])
    t_out.Branch("pionzero0_score", pionzero0_score, "fd_pionzero0_score/F")
    pionzero1_score = array("f", [0])
    t_out.Branch("pionzero1_score", pionzero1_score, "fd_pionzero1_score/F")
    pionzero2_score = array("f", [0])
    t_out.Branch("pionzero2_score", pionzero2_score, "fd_pionzero2_score/F")
    pionzeroN_score = array("f", [0])
    t_out.Branch("pionzeroN_score", pionzeroN_score, "fd_pionzeroN_score/F")
    neutron0_score = array("f", [0])
    t_out.Branch("neutron0_score", neutron0_score, "fd_neutron0_score/F")
    neutron1_score = array("f", [0])
    t_out.Branch("neutron1_score", neutron1_score, "fd_neutron1_score/F")
    neutron2_score = array("f", [0])
    t_out.Branch("neutron2_score", neutron2_score, "fd_neutron2_score/F")
    neutronN_score = array("f", [0])
    t_out.Branch("neutronN_score", neutronN_score, "fd_neutronN_score/F")
    numu_nu_E = array("f", [0])
    t_out.Branch("numu_nu_E", numu_nu_E, "fd_numu_nu_E/F")
    numu_had_E = array("f", [0])
    t_out.Branch("numu_had_E", numu_had_E, "fd_numu_had_E/F")
    numu_lep_E = array("f", [0])
    t_out.Branch("numu_lep_E", numu_lep_E, "fd_numu_lep_E/F")
    numu_recomethod = array("i", [0])
    t_out.Branch("numu_recomethod", numu_recomethod, "fd_numu_recomethod/I")
    numu_longesttrackcontained = array("i", [0])
    t_out.Branch(
        "fd_numu_longesttrackcontained", numu_longesttrackcontained, "numu_longesttrackcontained/I"
    )
    numu_trackmommethod = array("i", [0])
    t_out.Branch("numu_trackmommethod", numu_trackmommethod, "fd_numu_trackmommethod/I")
    nue_nu_E = array("f", [0])
    t_out.Branch("nue_nu_E", nue_nu_E, "fd_nue_nu_E/F")
    nue_had_E = array("f", [0])
    t_out.Branch("nue_had_E", nue_had_E, "fd_nue_had_E/F")
    nue_lep_E = array("f", [0])
    t_out.Branch("nue_lep_E", nue_lep_E, "fd_nue_lep_E/F")
    nue_recomethod = array("i", [0])
    t_out.Branch("nue_recomethod", nue_recomethod, "fd_nue_recomethod/I")
    nc_nu_E = array("f", [0])
    t_out.Branch("nc_nu_E", nc_nu_E, "fd_nc_nu_E/F")
    nc_had_E = array("f", [0])
    t_out.Branch("nc_had_E", nc_had_E, "fd_nc_had_E/F")
    nc_lep_E = array("f", [0])
    t_out.Branch("nc_lep_E", nc_lep_E, "fd_nc_lep_E/F")
    nc_recomethod = array("i", [0])
    t_out.Branch("nc_recomethod", nc_recomethod, "fd_nc_recomethod/I")

    event_cnt = 0

    for i_file, fd_filepath in enumerate(fd_files):
        if i_file + 1 % 100 == 0:
            print("{} files processed -- {} events merged".format(i_file, event_cnt))

        # XXX test
        if i_file > 20:
            break

        f_fd = ROOT.TFile.Open(fd_filepath)
        t_fd = f_fd.Get("recodump/FDReco")

        num = re.findall("FHC\.[0-9]*\.", os.path.basename(fd_filepath))[0][4:-1]

        # Annoying genie warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            f_nd = ROOT.TFile.Open(
                subprocess.check_output(
                    "pnfs2xrootd {}".format(
                        os.path.join(args.nd_dir, "FHC.{}.CAF.root".format(num))
                    ),
                    shell=True
                ).decode().strip("\n")
            )
        # Some caf files did not get the caf tree, grid job failure I missed I suppose
        if not f_nd.GetListOfKeys().Contains("caf"):
            print("FHC.{}.CAF.root missing caf tree".format(num))
            # f_nd.Close()
            # f_fd.Close()
            continue
        t_nd = f_nd.Get("caf")

        for e_fd in t_fd:
            event_num[0] = int(e_fd.eventid)
            print(event_num[0])
            for e_nd in t_nd:
                if int(e_nd.event) != event_num[0]:
                    continue

                run_num[0] = int(e_nd.run)

                is_cc[0] = int(e_nd.isCC)
                nu_pdg[0] = int(e_nd.nuPDG)
                nu_mom_x[0] = float(e_nd.NuMomX)
                nu_mom_y[0] = float(e_nd.NuMomY)
                nu_mom_z[0] = float(e_nd.NuMomZ)
                nu_E[0] = float(e_nd.Ev)
                mode[0] = int(e_nd.mode)
                lep_pdg[0] = int(e_nd.LepPDG)
                lep_mom_x[0] = float(e_nd.LepMomX)
                lep_mom_y[0] = float(e_nd.LepMomY)
                lep_mom_z[0] = float(e_nd.LepMomZ)
                lep_E[0] = float(e_nd.LepE)
                lep_nu_angle[0] = float(e_nd.LepNuAngle)
                num_prtn[0] = int(e_nd.nP)
                num_ntrn[0] = int(e_nd.nN)
                num_pip[0] = int(e_nd.nipip)
                num_pim[0] = int(e_nd.nipim)
                num_pi0[0] = int(e_nd.nipi0)
                num_kaonp[0] = int(e_nd.nikp)
                num_kaonm[0] = int(e_nd.nikm)
                num_kaon0[0] = int(e_nd.nik0)
                num_e[0] = int(e_nd.niem)
                num_other[0] = int(e_nd.niother)
                num_nuc[0] = int(e_nd.nNucleus)
                num_unkwn[0] = int(e_nd.nUNKNOWN)
                prtn_E[0] = float(e_nd.eP)
                ntrn_E[0] = float(e_nd.eN)
                pip_E[0] = float(e_nd.ePip)
                pim_E[0] = float(e_nd.ePim)
                pi0_E[0] = float(e_nd.ePi0)
                other_E[0] = float(e_nd.eOther)
                prtn_E_reco[0] = float(e_nd.eRecoP)
                ntrn_E_reco[0] = float(e_nd.eRecoN)
                pip_E_reco[0] = float(e_nd.eRecoPip)
                pim_E_reco[0] = float(e_nd.eRecoPim)
                pi0_E_reco[0] = float(e_nd.eRecoPi0)
                other_E_reco[0] = float(e_nd.eRecoOther)
                vtx_x[0] = float(e_nd.vtx_x)
                vtx_y[0] = float(e_nd.vtx_y)
                vtx_z[0] = float(e_nd.vtx_z)
                nu_E_reco[0] = float(e_nd.Ev_reco)
                lep_E_reco[0] = float(e_nd.Elep_reco)
                theta_reco[0] = float(e_nd.theta_reco)
                reco_numu[0] = int(e_nd.reco_numu)
                reco_nue[0] = int(e_nd.reco_nue)
                reco_nc[0] = int(e_nd.reco_nc)
                reco_q[0] = int(e_nd.reco_q)
                muon_contained[0] = int(e_nd.muon_contained)
                muon_tracker[0] = int(e_nd.muon_tracker)
                muon_ecal[0] = int(e_nd.muon_ecal)
                muon_exit[0] = int(e_nd.muon_exit)
                lep_pdg_reco[0] = int(e_nd.reco_lepton_pdg)
                muon_endpnt_x[0] = float(e_nd.muon_endpoint[0])
                muon_endpnt_y[0] = float(e_nd.muon_endpoint[1])
                muon_endpnt_z[0] = float(e_nd.muon_endpoint[2])
                hadE_veto[0] = float(e_nd.Ehad_veto)

                numu_score[0] = float(e_fd.numuScore)
                nue_score[0] = float(e_fd.nueScore)
                nc_score[0] = float(e_fd.NCScore)
                nutau_score[0] = float(e_fd.nutauScore)
                antinu_score[0] = float(e_fd.antiNuScore)
                proton0_score[0] = float(e_fd.proton0Score)
                proton1_score[0] = float(e_fd.proton1Score)
                proton2_score[0] = float(e_fd.proton2Score)
                protonN_score[0] = float(e_fd.protonNScore)
                pion0_score[0] = float(e_fd.pion0Score)
                pion1_score[0] = float(e_fd.pion1Score)
                pion2_score[0] = float(e_fd.pion2Score)
                pionN_score[0] = float(e_fd.pionNScore)
                pionzero0_score[0] = float(e_fd.pionZero0Score)
                pionzero1_score[0] = float(e_fd.pionZero1Score)
                pionzero2_score[0] = float(e_fd.pionZero2Score)
                pionzeroN_score[0] = float(e_fd.pionZeroNScore)
                neutron0_score[0] = float(e_fd.neutron0Score)
                neutron1_score[0] = float(e_fd.neutron1Score)
                neutron2_score[0] = float(e_fd.neutron2Score)
                neutronN_score[0] = float(e_fd.neutronNScore)
                numu_nu_E[0] = float(e_fd.numuNuE)
                numu_had_E[0] = float(e_fd.numuHadE)
                numu_lep_E[0] = float(e_fd.numuLepE)
                numu_recomethod[0] = int(e_fd.numuRecoMethod)
                numu_longesttrackcontained[0] = int(e_fd.numuLongestTrackContained)
                numu_trackmommethod[0] = int(e_fd.numuTrackMomMethod)
                nue_nu_E[0] = float(e_fd.nueNuE)
                nue_had_E[0] = float(e_fd.nueHadE)
                nue_lep_E[0] = float(e_fd.nueLepE)
                nue_recomethod[0] = int(e_fd.nueRecoMethod)
                nc_nu_E[0] = float(e_fd.NCNuE)
                nc_had_E[0] = float(e_fd.NCHadE)
                nc_lep_E[0] = float(e_fd.NCLepE)
                nc_recomethod[0] = int(e_fd.NCRecoMethod)

                event_cnt += 1

                break

            t_out.Fill()

        # f_fd.Close()
        # f_fd.Close()

    print("Finished: {} files processed -- {} events merged".format(i_file + 1, event_cnt))

    f_out.Write()
    # f_out.Close()


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("nd_dir", type=str)
    parser.add_argument("fd_dir", type=str)

    parser.add_argument("--output_path", default="", type=str)

    return parser.parse_args()


if __name__ == '__main__':
    main(parse_arguments())

