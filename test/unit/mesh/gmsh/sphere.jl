# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# This script uses Gmsh to generate a spherical mesh.

# The mesh has two attributes:
# - 1 (`domain`): the 3D domain
# - 2 (`boundary`): the 2D boundary
#
# Set up (only the first time):
#   julia --project -e 'using Pkg; Pkg.instantiate()'
#
# Run with
#   julia --project -e 'include("sphere.jl")'

using Gmsh: gmsh

radius = 1.5
mesh_size = 0.25  # Controls mesh resolution (smaller = finer mesh)
filename = "sphere.msh"

gmsh.initialize()
x = y = z = 0.0
boundary = gmsh.model.occ.addSphere(x, y, z, radius)

gmsh.model.occ.synchronize()

gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size)

gmsh.model.mesh.generate(3)     # 3D mesh
gmsh.model.mesh.setOrder(3)     # Cubically curved elements

extract_tag(object) = last(object) # Find tag from tuple of (ndim, tag)

gmsh.model.addPhysicalGroup(3, extract_tag.(gmsh.model.occ.getEntities(3)), -1, "domain")
gmsh.model.addPhysicalGroup(2, extract_tag.(gmsh.model.occ.getEntities(2)), -1, "boundary")

gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write(joinpath(@__DIR__, filename))
gmsh.finalize()
