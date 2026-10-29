# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate the mesh for the parallel-plate London tests (from the parallel_plate_london
# directory):
#   julia -e 'include("mesh/mesh.jl"); generate_parallel_plate_mesh()'

"""
    generate_parallel_plate_mesh(; filename, t, d, W, ℓ, nx, nt, nd, grading)

Generate a hexahedral mesh of a parallel-plate line: two slabs of thickness `t` separated
by a gap `d`, of width `W` and length `ℓ` along x. The mesh is a tensor product, graded
geometrically in z toward the slab faces bounding the gap, and is written directly in the
Gmsh 2.2 ASCII format.

# Attributes

  - Volume 1: the gap
  - Volume 2: the lower slab, 0 ≤ z ≤ t
  - Volume 3: the upper slab, t + d ≤ z ≤ 2t + d
  - Boundary 2: the slab end faces at x = 0
  - Boundary 3: the gap end face at x = 0
  - Boundary 4: the far end x = ℓ
  - Boundary 5: the other exterior faces
  - Boundary 6: the interior face z = t + d between the gap and the upper slab

Units: μm (Palace L0 = 1e-6).
"""
function generate_parallel_plate_mesh(;
    filename::AbstractString="mesh/parallel_plate.msh",
    t::Real=0.2,
    d::Real=0.5,
    W::Real=1.0,
    ℓ::Real=4.0,
    nx::Integer=2,
    nt::Integer=8,
    nd::Integer=4,
    grading::Real=1.3
)
    # Nodes in [z0, z1] with sizes growing geometrically away from the end at fine.
    function graded(z0, z1, n, fine_at_end)
        h = [grading^(i - 1) for i = 1:n]
        fine_at_end && reverse!(h)
        return z0 .+ (z1 - z0) .* cumsum(vcat(0.0, h)) ./ sum(h)
    end
    xs = collect(range(0.0, ℓ; length=nx + 1))
    ys = [0.0, W]
    zs = vcat(
        graded(0.0, t, nt, true),
        collect(range(t, t + d; length=nd + 1))[2:end],
        graded(t + d, 2t + d, nt, false)[2:end]
    )
    nxp, nyp, nzp = length(xs), length(ys), length(zs)
    node(i, j, k) = i + (j - 1) * nxp + (k - 1) * nxp * nyp
    zmid(k) = (zs[k] + zs[k + 1]) / 2
    region(k) = zmid(k) < t ? 2 : (zmid(k) > t + d ? 3 : 1)

    # Hexahedra (Gmsh type 5) and boundary quadrilaterals (type 3).
    elements = Tuple{Int, Int, Vector{Int}}[]
    for k = 1:(nzp - 1), j = 1:(nyp - 1), i = 1:(nxp - 1)
        push!(
            elements,
            (
                5,
                region(k),
                [
                    node(i, j, k),
                    node(i + 1, j, k),
                    node(i + 1, j + 1, k),
                    node(i, j + 1, k),
                    node(i, j, k + 1),
                    node(i + 1, j, k + 1),
                    node(i + 1, j + 1, k + 1),
                    node(i, j + 1, k + 1)
                ]
            )
        )
    end
    for k = 1:(nzp - 1), j = 1:(nyp - 1)
        x0 = [node(1, j, k), node(1, j + 1, k), node(1, j + 1, k + 1), node(1, j, k + 1)]
        x1 = [
            node(nxp, j, k),
            node(nxp, j + 1, k),
            node(nxp, j + 1, k + 1),
            node(nxp, j, k + 1)
        ]
        push!(elements, (3, region(k) == 1 ? 3 : 2, x0))
        push!(elements, (3, 4, x1))
    end
    for k = 1:(nzp - 1), i = 1:(nxp - 1), j in (1, nyp)
        push!(
            elements,
            (
                3,
                5,
                [
                    node(i, j, k),
                    node(i + 1, j, k),
                    node(i + 1, j, k + 1),
                    node(i, j, k + 1)
                ]
            )
        )
    end
    kf = findfirst(z -> isapprox(z, t + d), zs)
    for j = 1:(nyp - 1), i = 1:(nxp - 1), (k, attr) in ((1, 5), (nzp, 5), (kf, 6))
        push!(
            elements,
            (
                3,
                attr,
                [
                    node(i, j, k),
                    node(i + 1, j, k),
                    node(i + 1, j + 1, k),
                    node(i, j + 1, k)
                ]
            )
        )
    end

    open(joinpath(@__DIR__, "..", filename), "w") do io
        println(io, "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat")
        println(io, "\$Nodes\n", nxp * nyp * nzp)
        for k = 1:nzp, j = 1:nyp, i = 1:nxp
            println(io, node(i, j, k), " ", xs[i], " ", ys[j], " ", zs[k])
        end
        println(io, "\$EndNodes\n\$Elements\n", length(elements))
        for (e, (type, attr, nodes)) in enumerate(elements)
            println(io, e, " ", type, " 2 ", attr, " ", attr, " ", join(nodes, " "))
        end
        return println(io, "\$EndElements")
    end
    return nothing
end
