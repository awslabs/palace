# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Derive the slab mesh from the grating mesh; run from the repository root.

Requires meshio 5.3.5 and NumPy. The original topology and periodic node mappings
are retained; only physical volume tags in the middle layer are changed.
"""

from pathlib import Path

import meshio
import numpy as np

mesh = meshio.read("examples/dielectric_grating/mesh/dielectric_grating.msh")
assert mesh.gmsh_periodic is not None
slab_volume = 0.0
for block, tags in zip(mesh.cells, mesh.cell_data["gmsh:physical"]):
    if block.type != "hexahedron":
        continue
    vertices = mesh.points[block.data]
    center = vertices.mean(axis=1)
    selected = np.abs(center[:, 2]) < 0.25
    lower, upper = vertices.min(axis=1), vertices.max(axis=1)
    assert np.all(lower[selected, 2] >= -0.25 - 1e-12)
    assert np.all(upper[selected, 2] <= 0.25 + 1e-12)
    tags[:] = np.where(selected, 2, 1)
    slab_volume += np.prod(upper[selected] - lower[selected], axis=1).sum()
assert abs(slab_volume - 2.0) < 1e-10  # cm³
meshio.write(
    Path(__file__).with_name("slab.msh"), mesh, file_format="gmsh22", binary=False
)
