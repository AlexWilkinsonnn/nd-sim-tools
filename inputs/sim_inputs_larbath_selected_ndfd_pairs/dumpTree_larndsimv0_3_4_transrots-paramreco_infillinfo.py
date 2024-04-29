"""
Converts ROOT file created by edep-sim and passed through translations and rotations
into HDF5 format. Optionally adds ND parameterised reco (this would be produced with the same genie
info processed through edep-sim with the ND geom rather than larbath)
"""
from math import sqrt

import numpy as np
import h5py, fire

from ROOT import TFile

# Output array datatypes
segments_dtype = np.dtype([("eventID", "u4"), ("trackID", "u4"), ("uniqID", "u4"), ("pdgID", "i4"),
                           ("x_start", "f4"), ("x_end", "f4"), ("x", "f4"),
                           ("y_start", "f4"), ("y_end", "f4"), ("y", "f4"),
                           ("z_start", "f4"), ("z_end", "f4"), ("z", "f4"),
                           ("t_start", "f4"), ("t_end", "f4"), ("t", "f4"),
                           ("t0_start", "f4"), ("t0_end", "f4"), ("t0", "f4"),
                           ("n_electrons", "u4"), ("n_photons", "u4"),
                           ("tran_diff", "f4"), ("long_diff", "f4"),
                           ("dx", "f4"), ("dEdx", "f4"), ("dE", "f4"),
                           ("pixel_plane", "i4")])

vertices_dtype = np.dtype([("eventID","u4"),("x_vert","f4"),("y_vert","f4"),("z_vert","f4")])

depos_dtype = np.dtype([("eventID", "u4"), ("uniqID", "u4"),
                        ("x_start", "f4"), ("x_end", "f4"), ("x", "f4"),
                        ("y_start", "f4"), ("y_end", "f4"), ("y", "f4"),
                        ("z_start", "f4"), ("z_end", "f4"), ("z", "f4"),
                        ("t0_start", "f8"), ("t0_end", "f8"), ("t0", "f8"),
                        ("dx", "f4"), ("dEdx", "f4"), ("dE", "f4"),
                        ("outsideNDLAr", "i4")])

paramreco_dtype = np.dtype([("eventID", "u4"), ("cafTree_event", "u4"),
                            ("isFHC", "u4"), ("isCC", "u4"),
                            ("nuPDG", "i4"),
                            ("NuMomX", "f4"), ("NuMomY", "f4"), ("NuMomZ", "f4"),
                            ("Ev", "f4"),
                            ("mode", "i4"),
                            ("LepPDG", "i4"),
                            ("LepMomX", "f4"), ("LepMomY", "f4"), ("LepMomZ", "f4"),
                            ("LepE", "f4"), ("LepNuAngle", "f4"),
                            ("Q2", "f4"), ("W", "f4"), ("X", "f4"), ("Y", "f4"),
                            ("nP", "u4"), ("nN", "u4"),
                            ("nipip", "u4"), ("nipim", "u4"), ("nipi0", "u4"),
                            ("nikp", "u4"), ("nikm", "u4"), ("nik0", "u4"),
                            ("niem", "u4"), ("niother", "u4"),
                            ("nNucleus", "u4"), ("nUNKNOWN", "u4"),
                            ("eP", "f4"), ("eN", "f4"),
                            ("ePip", "f4"), ("ePim", "f4"), ("ePi0", "f4"), ("eOther", "f4"),
                            ("eRecoP", "f4"), ("eRecoN", "f4"),
                            ("eRecoPip", "f4"), ("eRecoPim", "f4"), ("eRecoPi0", "f4"),
                            ("eRecoOther", "f4"),
                            ("det_x", "i4"),
                            ("vtx_x", "f4"), ("vtx_y", "f4"), ("vtx_z", "f4"),
                            ("Ev_reco", "f4"), ("Elep_reco", "f4"),
                            ("theta_reco", "f4"),
                            ("reco_numu", "i4"), ("reco_nue", "i4"),
                            ("reco_nc", "i4"), ("reco_q", "i4"),
                            ("muon_contained", "u4"), ("muon_tracker", "u4"),
                            ("muon_ecal", "u4"), ("muon_exit", "u4"),
                            ("reco_lepton_pdg", "u4"),
                            ("muon_endpoint_x", "f4"), ("muon_endpoint_y", "f4"),
                            ("muon_endpoint_z", "f4"),
                            ("Ehad_veto", "f4"),
                            ("vtxInGap", "i4"), ("hadEFracInGap", "f4"), ("lepEFracInGap", "f4")])

primaries_dtype = np.dtype([("eventID", "u4"),
                            ("pdg", "u4"),
                            ("p0_MeV", "f4"), ("p1_MeV", "f4"), ("p2_MeV", "f4"), ("p3_MeV", "f4")])

lepton_dtype = np.dtype([("eventID", "u4"),
                         ("pdg", "u4"),
                         ("p0_MeV", "f4"), ("p1_MeV", "f4"), ("p2_MeV", "f4"), ("p3_MeV", "f4"),
                         ("nn_lep_contained_prob", "f4"), ("nn_lep_tracker_prob", "f4"),
                         ("lep_ke_outside_ndlar_MeV", "f4")])

# Prep HDF5 file for writing
def initHDF5File(output_file, paramreco):
    with h5py.File(output_file, 'w') as f:
        f.create_dataset('segments', (0,), dtype=segments_dtype, maxshape=(None,))
        f.create_dataset('vertices', (0,), dtype=vertices_dtype, maxshape=(None,))
        f.create_dataset('fd_deps', (0,), dtype=depos_dtype, maxshape=(None,))
        f.create_dataset('fd_vertices', (0,), dtype=vertices_dtype, maxshape=(None,))
        if paramreco:
            f.create_dataset('nd_paramreco', (0,), dtype=paramreco_dtype, maxshape=(None,))
        f.create_dataset('primaries', (0,), dtype=primaries_dtype, maxshape=(None,))
        f.create_dataset('lepton', (0,), dtype=lepton_dtype, maxshape=(None,))

# Resize HDF5 file and save output arrays
def updateHDF5File(
    output_file, segments, vertices, fd_deps, fd_vertices, nd_paramreco, primaries, lepton
):
    if any([
        len(segments), len(vertices), len(fd_deps), len(fd_vertices), len(nd_paramreco),
        len(primaries), len(lepton)
    ]):
        with h5py.File(output_file, 'a') as f:
            if len(segments):
                nseg = len(f['segments'])
                f['segments'].resize((nseg+len(segments),))
                f['segments'][nseg:] = segments

            if len(vertices):
                nvert = len(f['vertices'])
                f['vertices'].resize((nvert+len(vertices),))
                f['vertices'][nvert:] = vertices

            if len(fd_deps):
                ndeps = len(f['fd_deps'])
                f['fd_deps'].resize((ndeps+len(fd_deps),))
                f['fd_deps'][ndeps:] = fd_deps

            if len(fd_vertices):
                nvert = len(f['fd_vertices'])
                f['fd_vertices'].resize((nvert+len(fd_vertices),))
                f['fd_vertices'][nvert:] = fd_vertices

            if len(nd_paramreco):
                nprec = len(f['nd_paramreco'])
                f['nd_paramreco'].resize((nvert+len(nd_paramreco),))
                f['nd_paramreco'][nprec:] = nd_paramreco

            if len(primaries):
                nprim = len(f['primaries'])
                f['primaries'].resize((nvert+len(primaries),))
                f['primaries'][nprim:] = primaries

            if len(lepton):
                nlep = len(f['lepton'])
                f['lepton'].resize((nvert+len(lepton),))
                f['lepton'][nlep:] = lepton

# Read a file and dump it.
def dump(input_file, output_file, param_reco_file=None, min_nEdeps=20):
    # The input file is generated in a previous test (100TestTree.sh).
    inputFile = TFile(input_file)

    # Get the input tree out of the file.
    inputTree = inputFile.Get("myEvents")

    # Get the caf tree for the param reco if it exists
    if param_reco_file is not None:
        paramrecoFile = TFile(param_reco_file)
        paramrecoTree = paramrecoFile.Get("caf")
        if paramrecoTree.GetEntries() != inputTree.GetEntries():
            raise ValueError(
                "param reco tree and edep-sim tree have different number of events: {} and {}".format(
                    paramrecoTree.GetEntries(), inputTree.GetEntries()
                )
            )
        paramrecoTree_itr = iter(paramrecoTree)

    # Prep output file
    initHDF5File(output_file, param_reco_file is not None)

    segments_list = list()
    fd_depos_list = list()
    vertices_list = list()
    fd_vertices_list = list()
    param_reco_list = list()
    primaries_list = list()
    lepton_list = list()

    for i_event, event in enumerate(inputTree):
        if ((i_event) % 50 == 0):
            print(i_event)

        if not event.nd_fd_throws_passed:
            print("Throws for event {} failed, skipping".format(i_event))
            if param_reco_file is not None:
                next(paramrecoTree_itr)
            continue

        # Events that result in generate no packets in larnd-sim stage cause problems so try to
        # minimise these. Shouldn't be selecting these events anyway
        if event.nEdeps < min_nEdeps:
            print(
                "Event has too few energy deposition {} < {}, skipping".format(
                    event.nEdeps, min_nEdeps
                )
            )
            if param_reco_file is not None:
                next(paramrecoTree_itr)
            continue

        # Get the ND vertex
        vertex = np.empty(1, dtype=vertices_dtype)
        vertex["eventID"] = i_event
        vertex["x_vert"] = event.nd_vtx_cm_nonecc[0]
        vertex["y_vert"] = event.nd_vtx_cm_nonecc[1]
        vertex["z_vert"] = event.nd_vtx_cm_nonecc[2]
        vertices_list.append(vertex)

        # Dump the ND segments that have gone through geoEff translations and rotations
        segment = np.empty(event.nEdeps, dtype=segments_dtype)
        for (
            i_seg, (
                dep_E, start_t, stop_t, start_x, stop_x, start_y, stop_y, start_z, stop_z, trackID
            )
        ) in enumerate(
            zip(
                event.deps_E_MeV,
                event.deps_start_t_us, event.deps_stop_t_us,
                event.nd_deps_start_x_cm_nonecc, event.nd_deps_stop_x_cm_nonecc,
                event.nd_deps_start_y_cm_nonecc, event.nd_deps_stop_y_cm_nonecc,
                event.nd_deps_start_z_cm_nonecc, event.nd_deps_stop_z_cm_nonecc,
                event.deps_trkID
            )
        ):
            segment[i_seg]["eventID"] = i_event
            segment[i_seg]["x_start"] = start_x
            segment[i_seg]["y_start"] = start_y
            segment[i_seg]["z_start"] = start_z
            segment[i_seg]["t0_start"] = start_t
            segment[i_seg]["t_start"] = 0
            segment[i_seg]["x_end"] = stop_x
            segment[i_seg]["y_end"] = stop_y
            segment[i_seg]["z_end"] = stop_z
            segment[i_seg]["t0_end"] = stop_t
            segment[i_seg]["t_end"] = 0
            segment[i_seg]["dE"] = dep_E
            xd = segment[i_seg]["x_end"] - segment[i_seg]["x_start"]
            yd = segment[i_seg]["y_end"] - segment[i_seg]["y_start"]
            zd = segment[i_seg]["z_end"] - segment[i_seg]["z_start"]
            dx = sqrt(xd**2 + yd**2 + zd**2)
            segment[i_seg]["dx"] = dx
            segment[i_seg]["x"] = (
                (segment[i_seg]["x_start"] + segment[i_seg]["x_end"]) / 2.
            )
            segment[i_seg]["y"] = (
                (segment[i_seg]["y_start"] + segment[i_seg]["y_end"]) / 2.
            )
            segment[i_seg]["z"] = (
                (segment[i_seg]["z_start"] + segment[i_seg]["z_end"]) / 2.
            )
            segment[i_seg]["t0"] = (
                (segment[i_seg]["t0_start"] + segment[i_seg]["t0_end"]) / 2.
            )
            segment[i_seg]["t"] = 0
            segment[i_seg]["dEdx"] = segment[i_seg]["dE"] / dx if dx > 0 else 0
            segment[i_seg]["n_electrons"] = 0
            segment[i_seg]["long_diff"] = 0
            segment[i_seg]["tran_diff"] = 0
            segment[i_seg]["pixel_plane"] = 0
            segment[i_seg]["n_photons"] = 0
            segment[i_seg]["trackID"] = trackID
            segment[i_seg]["uniqID"] = i_seg
        segments_list.append(segment)

        # Dump the FD deps + vertex.
        # The deps do not need to be in the segment format required to larnd-sim
        vertex_fd = np.empty(1, dtype=vertices_dtype)
        vertex_fd["eventID"] = i_event
        vertex_fd["x_vert"] = event.fd_vtx_cm_pair_nd_nonecc[0]
        vertex_fd["y_vert"] = event.fd_vtx_cm_pair_nd_nonecc[1]
        vertex_fd["z_vert"] = event.fd_vtx_cm_pair_nd_nonecc[2]
        fd_vertices_list.append(vertex_fd)

        dep = np.empty(event.nEdeps, dtype=depos_dtype)
        for (
            i_dep, (dep_E, start_t, stop_t, start_x, stop_x, start_y, stop_y, start_z, stop_z)
        ) in enumerate(
            zip(
                event.deps_E_MeV,
                event.deps_start_t_us, event.deps_stop_t_us,
                event.fd_deps_start_x_cm_pair_nd_nonecc, event.fd_deps_stop_x_cm_pair_nd_nonecc,
                event.fd_deps_start_y_cm_pair_nd_nonecc, event.fd_deps_stop_y_cm_pair_nd_nonecc,
                event.fd_deps_start_z_cm_pair_nd_nonecc, event.fd_deps_stop_z_cm_pair_nd_nonecc,
            )
        ):
            dep[i_dep]["eventID"] = i_event
            dep[i_dep]["x_start"] = start_x
            dep[i_dep]["y_start"] = start_y
            dep[i_dep]["z_start"] = start_z
            dep[i_dep]["t0_start"] = start_t
            dep[i_dep]["x_end"] = stop_x
            dep[i_dep]["y_end"] = stop_y
            dep[i_dep]["z_end"] = stop_z
            dep[i_dep]["t0_end"] = stop_t
            dep[i_dep]["dE"] = dep_E
            xd = dep[i_dep]["x_end"] - dep[i_dep]["x_start"]
            yd = dep[i_dep]["y_end"] - dep[i_dep]["y_start"]
            zd = dep[i_dep]["z_end"] - dep[i_dep]["z_start"]
            dx = sqrt(xd**2 + yd**2 + zd**2)
            dep[i_dep]["dx"] = dx
            dep[i_dep]["x"] = (
                (dep[i_dep]["x_start"] + dep[i_dep]["x_end"]) / 2.
            )
            dep[i_dep]["y"] = (
                (dep[i_dep]["y_start"] + dep[i_dep]["y_end"]) / 2.
            )
            dep[i_dep]["z"] = (
                (dep[i_dep]["z_start"] + dep[i_dep]["z_end"]) / 2.
            )
            dep[i_dep]["t0"] = (
                (dep[i_dep]["t0_start"] + dep[i_dep]["t0_end"]) / 2.
            )
            dep[i_dep]["dEdx"] = dep[i_dep]["dE"] / dx if dx > 0 else 0
            dep[i_dep]["uniqID"] = i_dep
            dep[i_dep]["outsideNDLAr"] = -1
        fd_depos_list.append(dep)

        # Dump genie primaries
        primaries = np.empty(event.Genie_nParts, dtype=primaries_dtype)
        for i_prim, (pdg, p0, p1, p2, p3) in enumerate(
            zip(
                event.GenieParts_pdg,
                event.GenieParts_p0_MeV, event.GenieParts_p1_MeV,
                event.GenieParts_p2_MeV, event.GenieParts_p3_MeV
            )
        ):
            primaries[i_prim]["eventID"] = i_event
            primaries[i_prim]["pdg"] = pdg
            primaries[i_prim]["p0_MeV"] = p0
            primaries[i_prim]["p1_MeV"] = p1
            primaries[i_prim]["p2_MeV"] = p2
            primaries[i_prim]["p3_MeV"] = p3
        primaries_list.append(primaries)

        # Dump lepton kinematic at point of exiting ND-LAr
        lepton = np.empty(1, dtype=lepton_dtype)
        lepton["eventID"] = i_event
        lepton["pdg"] = event.Genie_final_lep_pdg
        lepton["p0_MeV"] = event.Genie_final_lep_p0_MeV
        lepton["p1_MeV"] = event.Genie_final_lep_p1_MeV
        lepton["p2_MeV"] = event.Genie_final_lep_p2_MeV
        lepton["p3_MeV"] = event.Genie_final_lep_p3_MeV
        lepton["nn_lep_contained_prob"] = event.nd_lep_contained_prob_nonecc
        lepton["nn_lep_tracker_prob"] = event.nd_lep_tracker_prob_nonecc
        lepton["lep_ke_outside_ndlar_MeV"] = event.nd_lep_ke_MeV_exit_ndlar_nonecc
        lepton_list.append(lepton)

        # Dump ND param reco results if provided
        if param_reco_file is not None:
            prec_event = next(paramrecoTree_itr)
            prec = np.empty(1, dtype=paramreco_dtype)
            prec["eventID"] = i_event
            prec["cafTree_event"] = prec_event.event
            prec["isFHC"] = prec_event.isFHC
            prec["isCC"] = prec_event.isCC
            prec["nuPDG"] = prec_event.nuPDG
            prec["NuMomX"] = prec_event.NuMomX
            prec["NuMomY"] = prec_event.NuMomY
            prec["NuMomZ"] = prec_event.NuMomZ
            prec["Ev"] = prec_event.Ev
            prec["mode"] = prec_event.mode
            prec["LepPDG"] = prec_event.LepPDG
            prec["LepMomX"] = prec_event.LepMomX
            prec["LepMomY"] = prec_event.LepMomY
            prec["LepMomZ"] = prec_event.LepMomZ
            prec["LepE"] = prec_event.LepE
            prec["LepNuAngle"] = prec_event.LepNuAngle
            prec["Q2"] = prec_event.Q2
            prec["W"] = prec_event.W
            prec["X"] = prec_event.X
            prec["Y"] = prec_event.Y
            prec["nP"] = prec_event.nP
            prec["nN"] = prec_event.nN
            prec["nipip"] = prec_event.nipip
            prec["nipim"] = prec_event.nipim
            prec["nipi0"] = prec_event.nipi0
            prec["nikp"] = prec_event.nikp
            prec["nikm"] = prec_event.nikm
            prec["nik0"] = prec_event.nik0
            prec["niem"] = prec_event.niem
            prec["niother"] = prec_event.niother
            prec["nNucleus"] = prec_event.nNucleus
            prec["nUNKNOWN"] = prec_event.nUNKNOWN
            prec["eP"] = prec_event.eP
            prec["eN"] = prec_event.eN
            prec["ePip"] = prec_event.ePip
            prec["ePim"] = prec_event.ePim
            prec["ePi0"] = prec_event.ePi0
            prec["eOther"] = prec_event.eOther
            prec["eRecoP"] = prec_event.eRecoP
            prec["eRecoN"] = prec_event.eRecoN
            prec["eRecoPip"] = prec_event.eRecoPip
            prec["eRecoPim"] = prec_event.eRecoPim
            prec["eRecoPi0"] = prec_event.eRecoPi0
            prec["eRecoOther"] = prec_event.eRecoOther
            prec["det_x"] = prec_event.det_x
            prec["vtx_x"] = prec_event.vtx_x
            prec["vtx_y"] = prec_event.vtx_y
            prec["vtx_z"] = prec_event.vtx_z
            prec["Ev_reco"] = prec_event.Ev_reco
            prec["Elep_reco"] = prec_event.Elep_reco
            prec["theta_reco"] = prec_event.theta_reco
            prec["reco_numu"] = prec_event.reco_numu
            prec["reco_nue"] = prec_event.reco_nue
            prec["reco_nc"] = prec_event.reco_nc
            prec["reco_q"] = prec_event.reco_q
            prec["muon_contained"] = prec_event.muon_contained
            prec["muon_tracker"] = prec_event.muon_tracker
            prec["muon_ecal"] = prec_event.muon_ecal
            prec["muon_exit"] = prec_event.muon_exit
            prec["reco_lepton_pdg"] = prec_event.reco_lepton_pdg
            prec["muon_endpoint_x"] = prec_event.muon_endpoint[0]
            prec["muon_endpoint_y"] = prec_event.muon_endpoint[1]
            prec["muon_endpoint_z"] = prec_event.muon_endpoint[2]
            prec["Ehad_veto"] = prec_event.Ehad_veto
            prec["vtxInGap"] = prec_event.vtxInGap
            prec["hadEFracInGap"] = prec_event.hadEFracInGap
            prec["lepEFracInGap"] = prec_event.lepEFracInGap
            param_reco_list.append(prec)

    # save any lingering data not written to file
    updateHDF5File(
        output_file,
        np.concatenate(segments_list, axis=0) if segments_list else np.empty((0,)),
        np.concatenate(vertices_list, axis=0) if vertices_list else np.empty((0,)),
        np.concatenate(fd_depos_list, axis=0) if fd_depos_list else np.empty((0,)),
        np.concatenate(fd_vertices_list, axis=0) if fd_vertices_list else np.empty((0,)),
        np.concatenate(param_reco_list, axis=0) if param_reco_list else np.empty((0,)),
        np.concatenate(primaries_list, axis=0) if primaries_list else np.empty((0,)),
        np.concatenate(lepton_list, axis=0) if lepton_list else np.empty((0,)))

if __name__ == "__main__":
    fire.Fire(dump)

