# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# This script uses Gmsh to generate a spherical mesh split in two top and
# bottom hemispheres.

# The mesh has three attributes:
# - 1 (`domain`): the 3D domain
# - 2 (`north`): the top 2D boundary
# - 3 (`south`): the bottom 2D boundary
#
# Set up (only the first time):
#   julia --project -e 'using Pkg; Pkg.instantiate()'
#
# Run with
#   julia --project -e 'include("two_hemispheres.jl")'

using Gmsh: gmsh

radius = 1.5
mesh_size = 0.125  # Controls mesh resolution (smaller = finer mesh)
filename = "two_hemispheres.msh"

gmsh.initialize()

# Create sphere
x = y = z = 0.0
sphere = gmsh.model.occ.addSphere(x, y, z, radius)

# Create cutting plane at z=0 to split hemispheres
plane = gmsh.model.occ.addPlaneSurface([
    gmsh.model.occ.addRectangle(-radius, -radius, 0, 2*radius, 2*radius)
])

# Fragment the sphere with the plane to create separate surfaces
gmsh.model.occ.fragment([(3, sphere)], [(2, plane)])
gmsh.model.occ.synchronize()

# Set mesh size
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size)

# Generate mesh
gmsh.model.mesh.generate(3)     # 3D mesh
gmsh.model.mesh.setOrder(3)     # Cubically curved elements

# Get all entities
volumes = gmsh.model.getEntities(3)
surfaces = gmsh.model.getEntities(2)

# Create physical groups
# Domain (all volumes)
gmsh.model.addPhysicalGroup(3, [tag for (dim, tag) in volumes], 1, "domain")

# Separate north and south hemispheres by z-coordinate
north_surfaces = []
south_surfaces = []

for (dim, tag) in surfaces
    # Get center of mass of surface
    com = gmsh.model.occ.getCenterOfMass(dim, tag)
    z_coord = com[3]

    if z_coord > 1e-10  # Northern hemisphere (z > 0)
        push!(north_surfaces, tag)
    elseif z_coord < -1e-10  # Southern hemisphere (z < 0)
        push!(south_surfaces, tag)
    end
    # Skip surfaces at z ≈ 0 (cutting plane)
end

# Add physical groups for boundaries
gmsh.model.addPhysicalGroup(2, north_surfaces, 2, "north")
gmsh.model.addPhysicalGroup(2, south_surfaces, 3, "south")

gmsh.option.setNumber("Mesh.Algorithm", 6)
gmsh.option.setNumber("Mesh.Algorithm3D", 1)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.option.setNumber("Mesh.Binary", 1)
gmsh.write(joinpath(@__DIR__, filename))
gmsh.finalize()
