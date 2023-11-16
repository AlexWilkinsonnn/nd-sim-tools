#! /usr/bin/env python
"""
Converts ROOT file created by edep-sim and passed through translations and rotations
into HDF5 format
"""

from math import sqrt

import numpy as np
import fire
import h5py

from array import array

from ROOT import TFile

# Output array datatypes
segments_dtype = np.dtype([("eventID", "u4"), ("z_end", "f4"),
                           ("trackID", "u4"), ("tran_diff", "f4"),
                           ("z_start", "f4"), ("x_end", "f4"),
                           ("y_end", "f4"), ("n_electrons", "u4"),
                           ("pdgId", "i4"), ("x_start", "f4"),
                           ("y_start", "f4"), ("t_start", "f4"),
                           ("t0_start", "f8"), ("t0_end", "f8"), ("t0", "f8"),
                           ("dx", "f4"), ("long_diff", "f4"),
                           ("pixel_plane", "i4"), ("t_end", "f4"),
                           ("dEdx", "f4"), ("dE", "f4"), ("t", "f4"),
                           ("y", "f4"), ("x", "f4"), ("z", "f4"),
                           ("n_photons","f4")])

trajectories_dtype = np.dtype([("eventID", "u4"), ("trackID", "u4"),
                               ("parentID", "i4"),
                               ("pxyz_start", "f4", (3,)),
                               ("xyz_start", "f4", (3,)), ("t_start", "f4"),
                               ("pxyz_end", "f4", (3,)),
                               ("xyz_end", "f4", (3,)), ("t_end", "f4"),
                               ("pdgId", "i4"), ("start_process", "u4"),
                               ("start_subprocess", "u4"),
                               ("end_process", "u4"),
                               ("end_subprocess", "u4")])

vertices_dtype = np.dtype([("eventID","u4"),("x_vert","f4"),("y_vert","f4"),("z_vert","f4")])

depos_dtype = np.dtype([("eventID", "u4"), ("z_end", "f4"),
                        ("z_start", "f4"), ("x_end", "f4"),
                        ("y_end", "f4"), ("x_start", "f4"),
                        ("y_start", "f4"), ("t_start", "f4"),
                        ("t0_start", "f8"), ("t0_end", "f8"), ("t0", "f8"),
                        ("dx", "f4"), ("dEdx", "f4"),
                        ("dE", "f4"), ("y", "f4"),
                        ("x", "f4"), ("z", "f4")])

# Print the fields in a TG4PrimaryParticle object
def printPrimaryParticle(depth, primaryParticle):
    print(depth,"Class: ", primaryParticle.ClassName())
    print(depth,"Track Id:", primaryParticle.GetTrackId())
    print(depth,"Name:", primaryParticle.GetName())
    print(depth,"PDG Code:",primaryParticle.GetPDGCode())
    print(depth,"Momentum:",primaryParticle.GetMomentum().X(),
          primaryParticle.GetMomentum().Y(),
          primaryParticle.GetMomentum().Z(),
          primaryParticle.GetMomentum().E(),
          primaryParticle.GetMomentum().P(),
          primaryParticle.GetMomentum().M())

# Print the fields in an TG4PrimaryVertex object
def printPrimaryVertex(depth, primaryVertex):
    print(depth,"Class: ", primaryVertex.ClassName())
    print(depth,"Position:", primaryVertex.GetPosition().X(),
          primaryVertex.GetPosition().Y(),
          primaryVertex.GetPosition().Z(),
          primaryVertex.GetPosition().T())
    print(depth,"Generator:",primaryVertex.GetGeneratorName())
    print(depth,"Reaction:",primaryVertex.GetReaction())
    print(depth,"Filename:",primaryVertex.GetFilename())
    print(depth,"InteractionNumber:",primaryVertex.GetInteractionNumber())
    depth = depth + ".."
    for infoVertex in primaryVertex.Informational:
        printPrimaryVertex(depth,infoVertex)
    for primaryParticle in primaryVertex.Particles:
        printPrimaryParticle(depth,primaryParticle)

# Print the fields in a TG4TrajectoryPoint object
def printTrajectoryPoint(depth, trajectoryPoint):
    print(depth,"Class: ", trajectoryPoint.ClassName())
    print(depth,"Position:", trajectoryPoint.GetPosition().X(),
          trajectoryPoint.GetPosition().Y(),
          trajectoryPoint.GetPosition().Z(),
          trajectoryPoint.GetPosition().T())
    print(depth,"Momentum:", trajectoryPoint.GetMomentum().X(),
          trajectoryPoint.GetMomentum().Y(),
          trajectoryPoint.GetMomentum().Z(),
          trajectoryPoint.GetMomentum().Mag())
    print(depth,"Process",trajectoryPoint.GetProcess())
    print(depth,"Subprocess",trajectoryPoint.GetSubprocess())

# Print the fields in a TG4Trajectory object
def printTrajectory(depth, trajectory):
    print(depth,"Class: ", trajectory.ClassName())
    depth = depth + ".."
    print(depth,"Track Id/Parent Id:",
          trajectory.GetTrackId(),
          trajectory.GetParentId())
    print(depth,"Name:",trajectory.GetName())
    print(depth,"PDG Code",trajectory.GetPDGCode())
    print(depth,"Initial Momentum:",trajectory.GetInitialMomentum().X(),
          trajectory.GetInitialMomentum().Y(),
          trajectory.GetInitialMomentum().Z(),
          trajectory.GetInitialMomentum().E(),
          trajectory.GetInitialMomentum().P(),
          trajectory.GetInitialMomentum().M())
    for trajectoryPoint in trajectory.Points:
        printTrajectoryPoint(depth,trajectoryPoint)

# Print the fields in a TG4HitSegment object
def printHitSegment(depth, hitSegment):
    print(depth,"Class: ", hitSegment.ClassName())
    print(depth,"Primary Id:", hitSegment.GetPrimaryId())
    print(depth,"Energy Deposit:",hitSegment.GetEnergyDeposit())
    print(depth,"Secondary Deposit:", hitSegment.GetSecondaryDeposit())
    print(depth,"Track Length:",hitSegment.GetTrackLength())
    print(depth,"Start:", hitSegment.GetStart().X(),
          hitSegment.GetStart().Y(),
          hitSegment.GetStart().Z(),
          hitSegment.GetStart().T())
    print(depth,"Stop:", hitSegment.GetStop().X(),
          hitSegment.GetStop().Y(),
          hitSegment.GetStop().Z(),
          hitSegment.GetStop().T())
    print(depth,"Contributor:", [contributor for contributor in hitSegment.Contrib])

# Print the fields in a single element of the SegmentDetectors map.
# The container name is the key, and the hitSegments is the value (a
# vector of TG4HitSegment objects).
def printSegmentContainer(depth, containerName, hitSegments):
    print(depth,"Detector: ", containerName, hitSegments.size())
    depth = depth + ".."
    for hitSegment in hitSegments: printHitSegment(depth, hitSegment)

# Prep HDF5 file for writing
def initHDF5File(output_file):
    with h5py.File(output_file, 'w') as f:
        # f.create_dataset('trajectories', (0,), dtype=trajectories_dtype, maxshape=(None,))
        f.create_dataset('segments', (0,), dtype=segments_dtype, maxshape=(None,))
        f.create_dataset('vertices', (0,), dtype=vertices_dtype, maxshape=(None,))
        f.create_dataset('fd_deps', (0,), dtype=depos_dtype, maxshape=(None,))
        f.create_dataset('fd_vertices', (0,), dtype=vertices_dtype, maxshape=(None,))

# Resize HDF5 file and save output arrays
def updateHDF5File(output_file, trajectories, segments, vertices, fd_deps, fd_vertices):
    if any([len(trajectories), len(segments), len(vertices), len(fd_deps), len(fd_vertices)]):
        with h5py.File(output_file, 'a') as f:
            if len(trajectories):
                ntraj = len(f['trajectories'])
                f['trajectories'].resize((ntraj+len(trajectories),))
                f['trajectories'][ntraj:] = trajectories

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

# Read a file and dump it.
def dump(input_file, output_file, nd_ecc=False):

    # The input file is generated in a previous test (100TestTree.sh).
    inputFile = TFile(input_file)

    # Get the input tree out of the file.
    inputTree = inputFile.Get("myEvents")

    # Prep output file
    initHDF5File(output_file)

    segments_list = list()
    fd_depos_list = list()
    trajectories_list = list()
    vertices_list = list()
    fd_vertices_list = list()

    for i_event, event in enumerate(inputTree):
        if ((i_event) % 10 == 0):
            print(i_event)

        # Dump the primary vertices
        vertex = np.empty(1, dtype=vertices_dtype)
        vertex["eventID"] = i_event
        if nd_ecc:
            vertex["x_vert"] = event.throwVtx_nd_ecc_cm[0]
            vertex["y_vert"] = event.throwVtx_nd_ecc_cm[1]
            vertex["z_vert"] = event.throwVtx_nd_ecc_cm[2]
        else:
            vertex["x_vert"] = event.throwVtx_nd_cm[0]
            vertex["y_vert"] = event.throwVtx_nd_cm[1]
            vertex["z_vert"] = event.throwVtx_nd_cm[2]
        vertices_list.append(vertex)

        # Dump the segments that have gone through geoEff translations and rotations
        segment = np.empty(event.nEdeps, dtype=segments_dtype)
        if nd_ecc:
            data_zip = zip(
                event.deps_E_MeV,
                event.deps_start_t_us, event.deps_stop_t_us,
                event.ndthrow_ecc_deps_start_x_cm, event.ndthrow_ecc_deps_stop_x_cm,
                event.ndthrow_ecc_deps_start_y_cm, event.ndthrow_ecc_deps_stop_y_cm,
                event.ndthrow_ecc_deps_start_z_cm, event.ndthrow_ecc_deps_stop_z_cm
            )
        else:
            data_zip = zip(
                event.deps_E_MeV,
                event.deps_start_t_us, event.deps_stop_t_us,
                event.ndthrow_deps_start_x_cm, event.ndthrow_deps_stop_x_cm,
                event.ndthrow_deps_start_y_cm, event.ndthrow_deps_stop_y_cm,
                event.ndthrow_deps_start_z_cm, event.ndthrow_deps_stop_z_cm
            )
        for (
            i_seg, (dep_E, start_t, stop_t, start_x, stop_x, start_y, stop_y, start_z, stop_z)
        ) in enumerate(data_zip):
            segment[i_seg]["eventID"] = i_event
            # segment[i_seg]["trackID"] = trajectories[hitSegment.Contrib[0]]["trackID"]
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
            # segment[i_seg]["pdgId"] = trajectories[hitSegment.Contrib[0]]["pdgId"]
            segment[i_seg]["n_electrons"] = 0
            segment[i_seg]["long_diff"] = 0
            segment[i_seg]["tran_diff"] = 0
            segment[i_seg]["pixel_plane"] = 0
            segment[i_seg]["n_photons"] = 0

        segments_list.append(segment)

        # Dump the FD deps. These do not need to be in the segment format required to larnd-sim
        vertex_fd = np.empty(1, dtype=vertices_dtype)
        vertex_fd["eventID"] = i_event
        vertex_fd["x_vert"] = event.throwVtx_fd_cm[0]
        vertex_fd["y_vert"] = event.throwVtx_fd_cm[1]
        vertex_fd["z_vert"] = event.throwVtx_fd_cm[2]
        fd_vertices_list.append(vertex_fd)

        dep = np.empty(event.nEdeps, dtype=depos_dtype)
        for (
            i_dep, (dep_E, start_t, stop_t, start_x, stop_x, start_y, stop_y, start_z, stop_z)
        ) in enumerate(
            zip(
                event.deps_E_MeV,
                event.deps_start_t_us, event.deps_stop_t_us,
                event.fdthrow_deps_start_x_cm, event.fdthrow_deps_stop_x_cm,
                event.fdthrow_deps_start_y_cm, event.fdthrow_deps_stop_y_cm,
                event.fdthrow_deps_start_z_cm, event.fdthrow_deps_stop_z_cm
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

        fd_depos_list.append(dep)

    # save any lingering data not written to file
    updateHDF5File(
        output_file,
        np.concatenate(trajectories_list, axis=0) if trajectories_list else np.empty((0,)),
        np.concatenate(segments_list, axis=0) if segments_list else np.empty((0,)),
        np.concatenate(vertices_list, axis=0) if vertices_list else np.empty((0,)),
        np.concatenate(fd_depos_list, axis=0) if fd_depos_list else np.empty((0,)),
        np.concatenate(fd_vertices_list, axis=0) if fd_vertices_list else np.empty((0,)))

if __name__ == "__main__":
    fire.Fire(dump)
