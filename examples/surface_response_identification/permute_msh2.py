#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Renumber a MSH 2.2 mesh (a seeded permutation of the node tags and of the element order within
every type, the corner nodes of every element rotated cyclically; optionally mirrored in x or
rotated about z in the plan view) so that the identification's tie-break freedom can be checked:
the manifest content of the permuted mesh must equal the canonical one (decision 82
infrastructure: canonical numbering by sorted quantized coordinates), the mirrored mesh must give
the same signatures with the opposite chirality, and the rotated mesh the same signatures and
chirality (the closed-loop start rule of the arc scan: the seed vertex of a loop is
coordinate-dependent; review fix-4 m-C).

    python3 -m surface_response_identification.permute_msh2 IN.msh2 OUT.msh2 [--seed N] [--mirror-x] [--rotate-degrees A]
"""
import argparse
import json
import sys

import numpy as np

from .msh2 import read_msh2
from .refine_msh2 import write_binary_msh2


def permute_file(mesh_path, output_path, seed=1, mirror_x=False, rotate_degrees=0.0):
    mesh = read_msh2(mesh_path)
    rng = np.random.default_rng(seed)
    coordinates = np.array(mesh.coordinates, dtype=float)
    if mirror_x:
        coordinates[:, 0] = -coordinates[:, 0]
    if rotate_degrees:
        # Rotation about z (the plan-view normal) by an arbitrary angle: the metal planes stay
        # planes, every in-plane direction and the canonical numbering's coordinate order change.
        angle = np.radians(rotate_degrees)
        x, y = coordinates[:, 0].copy(), coordinates[:, 1].copy()
        coordinates[:, 0] = np.cos(angle) * x - np.sin(angle) * y
        coordinates[:, 1] = np.sin(angle) * x + np.cos(angle) * y
    n = len(coordinates)
    order = rng.permutation(n)  # new position of every old node
    new_index = np.empty(n, dtype=np.int64)
    new_index[order] = np.arange(n)
    elements = {}
    for t in mesh.elements:
        physical = np.asarray(mesh.physical_tags(t))
        nodes = new_index[mesh.corner_indices(t)]
        # Rotate the corner nodes of every element by a seeded amount (the orientation of a
        # mirrored mesh is restored by swapping two nodes of every 2D / 3D element).
        k = nodes.shape[1]
        shift = rng.integers(0, k, size=len(physical))
        rotated = np.empty_like(nodes)
        for s in range(k):
            rows = shift == s
            rotated[rows] = np.roll(nodes[rows], -s, axis=1)
        if mirror_x and k >= 3:
            rotated[:, [0, 1]] = rotated[:, [1, 0]]
        perm = rng.permutation(len(physical))
        elements[t] = (physical[perm], rotated[perm])
    write_binary_msh2(output_path, coordinates[order], elements, mesh.physical_names)
    return {"Nodes": int(n), "Elements": {int(t): int(len(v[0])) for t, v in elements.items()}, "Seed": seed, "MirrorX": mirror_x,
            "RotateDegrees": rotate_degrees}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mesh")
    parser.add_argument("output")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--mirror-x", action="store_true")
    parser.add_argument("--rotate-degrees", type=float, default=0.0, help="rotate the plan view about z by this angle")
    args = parser.parse_args(argv)
    print(json.dumps(permute_file(args.mesh, args.output, seed=args.seed, mirror_x=args.mirror_x, rotate_degrees=args.rotate_degrees)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
