# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate the mesh for the parallel-plate London tests (from the parallel_plate_london
# directory):
#   julia -e 'include("mesh/mesh.jl"); generate_parallel_plate_mesh()'

"""
    generate_parallel_plate_mesh(; filename, t, d, W, ℓ, nx, nt, nd, grading)

Generate a hexahedral mesh of a two-port parallel-plate line: three slabs of thickness `t`
separated by two gaps `d`, of width `W` and length `ℓ` along x. The mesh is a tensor
product, graded geometrically in z toward the slab faces bounding the gaps, and is written
directly in the binary Gmsh 2.2 format.

# Attributes

  - Volume 1: the gaps
  - Volume 2: the lower slab, 0 ≤ z ≤ t
  - Volume 3: the middle slab, t + d ≤ z ≤ 2t + d
  - Volume 4: the upper slab, 2t + 2d ≤ z ≤ 3t + 2d
  - Boundary 2: the slab end faces at x = 0
  - Boundary 3: the lower gap end face at x = 0
  - Boundary 4: the far end x = ℓ
  - Boundary 5: the other exterior faces
  - Boundary 6: the interior face z = 2t + 2d between the upper gap and the upper slab
  - Boundary 7: the upper gap end face at x = 0

Units: μm (Palace L0 = 1e-6).
"""
function generate_parallel_plate_mesh(;
    filename::AbstractString="mesh/parallel_plate.msh",
    t::Real=0.2,
    d::Real=0.5,
    W::Real=1.0,
    ℓ::Real=4.0,
    nx::Integer=8,
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
    gap(z0) = collect(range(z0, z0 + d; length=nd + 1))[2:end]
    zm = 1.5t + d
    xs = collect(range(0.0, ℓ; length=nx + 1))
    ys = [0.0, W]
    zs = vcat(
        graded(0.0, t, nt, true),
        gap(t),
        graded(t + d, zm, nt ÷ 2, false)[2:end],
        graded(zm, 2t + d, nt ÷ 2, true)[2:end],
        gap(2t + d),
        graded(2t + 2d, 3t + 2d, nt, false)[2:end]
    )
    nxp, nyp, nzp = length(xs), length(ys), length(zs)
    node(i, j, k) = i + (j - 1) * nxp + (k - 1) * nxp * nyp
    zmid(k) = (zs[k] + zs[k + 1]) / 2
    function region(k)
        z = zmid(k)
        return z < t ? 2 : z < t + d ? 1 : z < 2t + d ? 3 : z < 2t + 2d ? 1 : 4
    end
    end_face(k) = region(k) != 1 ? 2 : (zmid(k) < 2t + d ? 3 : 7)

    # Hexahedra (Gmsh type 5) and boundary quadrilaterals (type 3).
    hexes = Tuple{Int, Vector{Int}}[]
    quads = Tuple{Int, Vector{Int}}[]
    for k = 1:(nzp - 1), j = 1:(nyp - 1), i = 1:(nxp - 1)
        push!(
            hexes,
            (
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
        push!(quads, (end_face(k), x0))
        push!(quads, (4, x1))
    end
    for k = 1:(nzp - 1), i = 1:(nxp - 1), j in (1, nyp)
        push!(
            quads,
            (
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
    kf = findfirst(z -> isapprox(z, 2t + 2d), zs)
    for j = 1:(nyp - 1), i = 1:(nxp - 1), (k, attr) in ((1, 5), (nzp, 5), (kf, 6))
        push!(
            quads,
            (
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

    # Binary Gmsh 2.2: a header per element type, then the tag, two attribute tags, and
    # nodes of each element.
    open(joinpath(@__DIR__, "..", filename), "w") do io
        print(io, "\$MeshFormat\n2.2 1 8\n")
        write(io, Int32(1))
        print(io, "\n\$EndMeshFormat\n\$Nodes\n", nxp * nyp * nzp, "\n")
        for k = 1:nzp, j = 1:nyp, i = 1:nxp
            write(io, Int32(node(i, j, k)), Float64(xs[i]), Float64(ys[j]), Float64(zs[k]))
        end
        print(io, "\n\$EndNodes\n\$Elements\n", length(hexes) + length(quads), "\n")
        tag = 0
        for (type, elements) in ((5, hexes), (3, quads))
            write(io, Int32(type), Int32(length(elements)), Int32(2))
            for (attr, nodes) in elements
                write(io, Int32(tag += 1), Int32(attr), Int32(attr), Int32.(nodes)...)
            end
        end
        return print(io, "\n\$EndElements\n")
    end
    return nothing
end
