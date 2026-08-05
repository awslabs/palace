# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using DeviceLayout

const gmsh = DeviceLayout.SolidModels.gmsh
const output = joinpath(@__DIR__, "singular_finite_metal_fixed.msh2")

gmsh.initialize()
try
    gmsh.model.add("singular_finite_metal_fixed")
    outer = gmsh.model.occ.add_box(-2.0, -2.0, -2.0, 4.0, 4.0, 4.0)
    metal = gmsh.model.occ.add_box(-0.5, -0.5, -0.5, 1.0, 1.0, 1.0)
    volumes, _ = gmsh.model.occ.cut([(3, outer)], [(3, metal)], -1, true, true)
    gmsh.model.occ.synchronize()
    volume_tags = Int32[tag for (dimension, tag) in volumes if dimension == 3]
    only(volume_tags)

    boundary = gmsh.model.get_boundary([(3, only(volume_tags))], false, false, false)
    exterior = Int32[]
    conductor = Int32[]
    for (dimension, tag) in boundary
        bounds = gmsh.model.get_bounding_box(dimension, tag)
        if any(isapprox(abs(value), 2.0; atol=1.0e-6) for value in bounds)
            push!(exterior, tag)
        else
            push!(conductor, tag)
        end
    end
    length(exterior) == 6 || error("Expected six exterior surfaces")
    length(conductor) == 6 || error("Expected six conductor surfaces")

    gmsh.model.add_physical_group(3, volume_tags, 1, "vacuum")
    gmsh.model.add_physical_group(2, exterior, 3, "exterior_boundary")
    gmsh.model.add_physical_group(2, conductor, 5, "metal")

    conductor_curves = unique(
        gmsh.model.get_boundary([(2, tag) for tag in conductor], false, false, false)
    )
    for (_, tag) in conductor_curves
        gmsh.model.mesh.set_transfinite_curve(tag, 2)
    end
    for tag in conductor
        gmsh.model.mesh.set_transfinite_surface(tag)
    end
    gmsh.option.set_number("General.NumThreads", 1)
    gmsh.option.set_number("Mesh.MaxNumThreads1D", 1)
    gmsh.option.set_number("Mesh.MaxNumThreads2D", 1)
    gmsh.option.set_number("Mesh.MaxNumThreads3D", 1)
    gmsh.option.set_number("Mesh.Algorithm3D", 10)
    gmsh.option.set_number("Mesh.MeshSizeMin", 0.2)
    gmsh.option.set_number("Mesh.MeshSizeMax", 0.5)
    gmsh.model.mesh.generate(3)
    gmsh.write(output)
    println(output)
finally
    gmsh.finalize()
end
