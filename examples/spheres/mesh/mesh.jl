# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using Gmsh: gmsh

"""
    generate_spheres_mesh(;
        filename::AbstractString,
        radius_a::Real   = 1.0,
        radius_b::Real   = 2.0,
        center_d::Real   = 5.0,
        verbose::Integer = 5,
        gui::Bool        = false
    )

Generate a mesh for the spheres example using Gmsh

# Arguments

  - filename - the filename to use for the generated mesh
  - radius_a - radius of the first sphere
  - radius_b - radius of the second sphere
  - center_d - distance between sphere centers
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
  - gui - whether to launch the Gmsh GUI on mesh generation
"""
function generate_spheres_mesh(;
    filename::AbstractString,
    radius_a::Real   = 1.0,
    radius_b::Real   = 2.0,
    center_d::Real   = 5.0,
    verbose::Integer = 5,
    gui::Bool        = false
)
    kernel = gmsh.model.occ

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    # Add model
    if "spheres" in gmsh.model.list()
        gmsh.model.setCurrent("spheres")
        gmsh.model.remove()
    end
    gmsh.model.add("spheres")

    # Geometry parameters (in cm)
    radius_farfield = 15.0 * center_d

    # Mesh parameters
    n_sphere = 16
    l_farfield = 20.0

    # Geometry
    sphere_a = kernel.addSphere(-0.5 * center_d, 0.0, 0.0, radius_a)
    sphere_b = kernel.addSphere(0.5 * center_d, 0.0, 0.0, radius_b)
    sphere_farfield = kernel.add_sphere(0.0, 0.0, 0.0, radius_farfield)

    cut_dimtags, cut_dimtags_map =
        kernel.cut([(3, sphere_farfield)], [(3, sphere_a), (3, sphere_b)])
    @assert length(cut_dimtags) == 1 && first(cut_dimtags)[1] == 3
    domain = first(cut_dimtags)[2]

    eps = 1.0e-3
    sphere_a = kernel.getEntitiesInBoundingBox(
        -0.5 * center_d - radius_a - eps,
        -radius_a - eps,
        -radius_a - eps,
        -0.5 * center_d + radius_a + eps,
        radius_a + eps,
        radius_a + eps,
        2
    )
    sphere_b = kernel.getEntitiesInBoundingBox(
        0.5 * center_d - radius_b - eps,
        -radius_b - eps,
        -radius_b - eps,
        0.5 * center_d + radius_b + eps,
        radius_b + eps,
        radius_b + eps,
        2
    )
    @assert length(sphere_a) == 1 && length(sphere_b) == 1
    sphere_a = first(sphere_a)[2]
    sphere_b = first(sphere_b)[2]

    sphere_farfield =
        filter(x -> x != (2, sphere_a) && x != (2, sphere_b), kernel.getEntities(2))
    @assert length(sphere_farfield) == 1
    sphere_farfield = first(sphere_farfield)[2]

    kernel.synchronize()

    # Add physical groups
    domain_group = gmsh.model.addPhysicalGroup(3, [domain], -1, "domain")

    farfield_group = gmsh.model.addPhysicalGroup(2, [sphere_farfield], -1, "farfield")

    sphere_a_group = gmsh.model.addPhysicalGroup(2, [sphere_a], -1, "sphere_a")
    sphere_b_group = gmsh.model.addPhysicalGroup(2, [sphere_b], -1, "sphere_b")

    # Generate mesh
    gmsh.option.setNumber("Mesh.MinimumCurveNodes", 2)
    gmsh.option.setNumber("Mesh.MinimumCircleNodes", 0)

    gmsh.option.setNumber("Mesh.MeshSizeMin", 2.0 * pi * radius_a / n_sphere / 2.0)
    gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", n_sphere)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    gmsh.model.mesh.field.add("Extend", 1)
    gmsh.model.mesh.field.setNumbers(1, "SurfacesList", [sphere_a, sphere_b])
    gmsh.model.mesh.field.setNumber(1, "Power", 1.0)
    gmsh.model.mesh.field.setNumber(1, "DistMax", radius_farfield)
    gmsh.model.mesh.field.setNumber(1, "SizeMax", l_farfield)

    gmsh.model.mesh.field.add("Min", 101)
    gmsh.model.mesh.field.setNumbers(101, "FieldsList", [1])
    gmsh.model.mesh.field.setAsBackgroundMesh(101)

    gmsh.option.setNumber("Mesh.Algorithm3D", 10)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(3)

    # Save mesh
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(joinpath(@__DIR__, filename))

    # Print some information
    println("\nFinished generating mesh. Physical group tags:")
    println("Domain: ", domain_group)
    println("Farfield boundary: ", farfield_group)
    println("Sphere A: ", sphere_a_group)
    println("Sphere B: ", sphere_b_group)
    println()

    # Optionally launch GUI
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
