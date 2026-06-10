# Single-transmon sensitivity chain: df/d(cap_length) on the real transmon geometry.
#
# Uses the existing Palace mesh sensitivity output (mesh_for_sensitivity.mesh +
# mesh_sensitivity.solb) and computes dx/dp via finite-difference geometry
# perturbation: render at cap_length and cap_length+δ, match Gmsh mesh vertices
# by nearest-neighbor, compute displacement.
#
# Prerequisites:
#   - Run `generate_transmon()` first (creates mesh/transmon.msh2 + transmon.json)
#   - Run Palace with MeshSensitivity enabled → postpro/mesh_for_sensitivity.mesh
#
# Usage: julia --project transmon/sensitivity_chain.jl

import DeviceLayout
import JSON
import Gmsh: gmsh

const ST_PATH = joinpath(pkgdir(DeviceLayout), "examples", "SingleTransmon")
include(joinpath(ST_PATH, "SingleTransmon.jl"))

using LinearAlgebra, Printf
using DeviceLayout.PreferredUnits

include(joinpath(@__DIR__, "..", "shared", "palace_io.jl"))

# ─── Generate a Gmsh mesh and return vertices ──────────────────────────────

"""
Generate the transmon mesh at given cap_length, return CORNER vertex positions.
Uses DeviceLayout's SingleTransmon pipeline (Gmsh/OCC kernel).
Extracts only 1st-order (corner) nodes to match Palace's mesh_for_sensitivity.mesh.
"""
function generate_mesh_vertices(; cap_length=620e-3, kwargs...)
    # Generate the solid model (geometry only), then mesh with deterministic settings
    sm = SingleTransmon.single_transmon(; save_mesh=false, cap_length=cap_length, kwargs...)
    # Force single-threaded meshing for deterministic vertex ordering
    gmsh.option.setNumber("General.NumThreads", 1)
    gmsh.model.mesh.generate(3)

    # Extract corner nodes only: for 2nd-order tets (type 11), first 4 of 10 nodes are corners
    elem_types, elem_tags, elem_nodes = gmsh.model.mesh.getElements(3)  # dim=3 (volumes)

    corner_tags = Set{Int}()
    for (etype, enodes) in zip(elem_types, elem_nodes)
        # Determine how many corner nodes per element
        _, _, _, nnodes, _ = gmsh.model.mesh.getElementProperties(etype)
        ncorner = etype == 4 ? 4 : etype == 11 ? 4 : etype == 5 ? 8 : etype == 17 ? 8 : nnodes
        ne = div(length(enodes), nnodes)
        for e in 1:ne
            offset = (e - 1) * nnodes
            for k in 1:ncorner
                push!(corner_tags, Int(enodes[offset + k]))
            end
        end
    end

    # Get coordinates for corner nodes
    all_tags, all_coords, _ = gmsh.model.mesh.getNodes()
    tag_to_coord = Dict{Int, Vector{Float64}}()
    for (i, t) in enumerate(all_tags)
        tag_to_coord[Int(t)] = all_coords[3*(i-1)+1:3*i]
    end

    sorted_tags = sort(collect(corner_tags))
    nv = length(sorted_tags)
    vertices = zeros(Float64, nv, 3)
    for (i, t) in enumerate(sorted_tags)
        vertices[i, :] .= tag_to_coord[t]
    end

    return (vertices=Matrix{Float64}(vertices), node_tags=sorted_tags,
            solid_model=sm)
end

# ─── Compute dx/dp via FD on Gmsh vertices ─────────────────────────────────

"""
Compute dx/dp at Palace mesh vertices via FD on the Gmsh geometry.

If Gmsh produces the same vertex count for both cap_length values (same mesh
topology), dx/dp is computed by direct index matching. This is the expected case
for small δ. Falls back to spatial hashing if counts differ.
"""
function compute_dxdp_fd(palace_mesh, is_boundary; cap_length=620e-3, delta=0.1e-3, delta_num=nothing)
    println("  Generating baseline Gmsh mesh...")
    baseline = generate_mesh_vertices(; cap_length=cap_length)

    println("  Generating perturbed Gmsh mesh (δ = $(delta))...")
    perturbed = generate_mesh_vertices(; cap_length=cap_length + delta)

    nv_base = size(baseline.vertices, 1)
    nv_pert = size(perturbed.vertices, 1)
    nv_palace = size(palace_mesh.vertices, 1)

    @printf("  Gmsh baseline: %d verts, perturbed: %d verts, Palace: %d verts\n",
            nv_base, nv_pert, nv_palace)

    # Gmsh→Gmsh: compute FD displacement via spatial hash matching
    # Use the stripped (unitless) delta for division — vertices are plain numbers in μm
    δ = isnothing(delta_num) ? Float64(delta) : delta_num

    gmsh_dxdp = zeros(Float64, nv_base, 3)
    if nv_base == nv_pert
        max_disp = maximum(norm(baseline.vertices[i,:] - perturbed.vertices[i,:])
                           for i in 1:nv_base)
        @printf("  Same vertex count → direct index match (max displacement: %.4e μm)\n", max_disp)
        gmsh_dxdp .= (perturbed.vertices .- baseline.vertices) ./ δ
    else
        @printf("  Vertex count mismatch (%d vs %d) — spatial hash matching...\n", nv_base, nv_pert)
        _spatial_hash_match_fd!(gmsh_dxdp, baseline.vertices, perturbed.vertices, δ)
    end

    # Palace→Gmsh mapping via spatial hash (counts will differ due to high-order nodes)
    dxdp = zeros(Float64, nv_palace, 3)
    if nv_palace == nv_base
        spot_err = maximum(norm(palace_mesh.vertices[i,:] - baseline.vertices[i,:])
                           for i in 1:min(100, nv_palace))
        @printf("  Palace→Gmsh: same count, spot check max err: %.4e μm\n", spot_err)
        if spot_err < 1.0  # vertices within 1 μm → direct mapping
            dxdp .= gmsh_dxdp
        else
            println("  Positions don't match — spatial hash mapping...")
            _spatial_hash_map!(dxdp, palace_mesh.vertices, baseline.vertices, gmsh_dxdp)
        end
    else
        @printf("  Palace (%d) ≠ Gmsh corners (%d) — spatial hash mapping...\n", nv_palace, nv_base)
        _spatial_hash_map!(dxdp, palace_mesh.vertices, baseline.vertices, gmsh_dxdp)
    end

    n_nonzero = count(maximum(abs, dxdp; dims=2)[:] .> 1e-12)
    @printf("  Vertices with |dx/dp| > 0: %d / %d\n", n_nonzero, nv_palace)

    return dxdp
end

"""Spatial hash FD: for each baseline vertex, find nearest perturbed vertex and compute FD."""
function _spatial_hash_match_fd!(dxdp, base_verts, pert_verts, delta)
    nb = size(base_verts, 1)
    np = size(pert_verts, 1)
    bbox_min = minimum(pert_verts; dims=1)
    bbox_max = maximum(pert_verts; dims=1)
    cell_size = max(norm(bbox_max - bbox_min) / sqrt(np), 1e-10)

    grid = Dict{Tuple{Int,Int,Int}, Vector{Int}}()
    for i in 1:np
        cx = floor(Int, pert_verts[i, 1] / cell_size)
        cy = floor(Int, pert_verts[i, 2] / cell_size)
        cz = floor(Int, pert_verts[i, 3] / cell_size)
        push!(get!(grid, (cx, cy, cz), Int[]), i)
    end

    for bi in 1:nb
        qx = floor(Int, base_verts[bi, 1] / cell_size)
        qy = floor(Int, base_verts[bi, 2] / cell_size)
        qz = floor(Int, base_verts[bi, 3] / cell_size)
        best_idx = 0
        best_dist = Inf
        for dx in -1:1, dy in -1:1, dz in -1:1
            key = (qx+dx, qy+dy, qz+dz)
            haskey(grid, key) || continue
            for pi in grid[key]
                d = norm(@view(base_verts[bi, :]) .- @view(pert_verts[pi, :]))
                if d < best_dist
                    best_dist = d
                    best_idx = pi
                end
            end
        end
        if best_idx > 0
            dxdp[bi, :] .= (pert_verts[best_idx, :] .- base_verts[bi, :]) ./ delta
        end
    end
end

"""Spatial hash map: for each query vertex, find nearest source vertex and copy its dxdp."""
function _spatial_hash_map!(dxdp, query_verts, source_verts, source_dxdp)
    nq = size(query_verts, 1)
    ns = size(source_verts, 1)

    # Build grid: cell_size = typical mesh spacing
    bbox_min = minimum(source_verts; dims=1)
    bbox_max = maximum(source_verts; dims=1)
    cell_size = max(norm(bbox_max - bbox_min) / sqrt(ns), 1e-10)

    # Hash source vertices into grid cells
    grid = Dict{Tuple{Int,Int,Int}, Vector{Int}}()
    for i in 1:ns
        cx = floor(Int, source_verts[i, 1] / cell_size)
        cy = floor(Int, source_verts[i, 2] / cell_size)
        cz = floor(Int, source_verts[i, 3] / cell_size)
        key = (cx, cy, cz)
        push!(get!(grid, key, Int[]), i)
    end

    # For each query vertex, search nearby cells
    for qi in 1:nq
        qx = floor(Int, query_verts[qi, 1] / cell_size)
        qy = floor(Int, query_verts[qi, 2] / cell_size)
        qz = floor(Int, query_verts[qi, 3] / cell_size)

        best_idx = 0
        best_dist = Inf
        for dx in -1:1, dy in -1:1, dz in -1:1
            key = (qx+dx, qy+dy, qz+dz)
            haskey(grid, key) || continue
            for si in grid[key]
                d = norm(@view(query_verts[qi, :]) .- @view(source_verts[si, :]))
                if d < best_dist
                    best_dist = d
                    best_idx = si
                end
            end
        end
        if best_idx > 0
            dxdp[qi, :] .= source_dxdp[best_idx, :]
        end
    end
end

# ─── Main ───────────────────────────────────────────────────────────────────

function main()
    cd(@__DIR__)

    palace = joinpath(@__DIR__, "..", "..", "build", "palace-build", "palace-arm64.bin")
    @assert isfile(palace) "Palace binary not found"

    println("═══ Single Transmon Sensitivity Chain ═══\n")

    # ─── Step 1: Read Palace mesh sensitivity ──────────────────────────────

    println("Step 1: Reading Palace mesh sensitivity output...")
    mesh_file = "postpro/mesh_for_sensitivity.mesh"
    solb_file = "postpro/mesh_sensitivity.solb"
    @assert isfile(mesh_file) "Run Palace first with MeshSensitivity enabled"
    @assert isfile(solb_file) "mesh_sensitivity.solb not found"

    palace_mesh = read_mesh(mesh_file)
    sol = read_solb(solb_file)
    dfdx = sol.data
    is_boundary = classify_boundary_vertices(palace_mesh)
    nv = size(palace_mesh.vertices, 1)

    f0 = read_eigenfrequency("postpro/eig.csv"; mode=1)
    @printf("  f₀ = %.6f GHz\n", f0)
    @printf("  Palace mesh: %d verts, %d boundary\n", nv, count(is_boundary))

    # ─── Step 2: Compute dx/dp via FD geometry perturbation ────────────────

    println("\nStep 2: Computing dx/d(cap_length) via FD...")
    # SingleTransmon.single_transmon uses Unitful: cap_length=620μm
    # Gmsh mesh coordinates are in μm (dimensionless numbers).
    # delta must be a Unitful quantity (for generate_mesh_vertices) but FD division
    # uses the stripped value.
    cap_length_val = 620μm
    delta_val = 10μm  # Larger δ improves signal-to-noise when Gmsh remeshes

    dxdp_all = compute_dxdp_fd(palace_mesh, is_boundary;
                                cap_length=cap_length_val, delta=delta_val,
                                delta_num=Float64(DeviceLayout.ustrip(delta_val)))

    # ─── Step 3: Laplacian propagation + chain contraction ─────────────────
    # If Palace→Gmsh mapping was direct, dxdp_all has values everywhere.
    # If spatial hash was used, interior vertices may be zero → propagate.
    # Either way, Laplacian propagation fills any gaps smoothly.

    println("\nStep 3: Laplacian propagation + chain contraction...")
    # Extract boundary-only values, propagate to interior
    dxdp_bdy = zeros(Float64, nv, 3)
    for vi in findall(is_boundary)
        dxdp_bdy[vi, :] .= dxdp_all[vi, :]
    end
    dxdp_full = propagate_boundary_sensitivity(palace_mesh, dxdp_bdy, is_boundary)

    dfdp = dot(vec(dfdx), vec(dxdp_full))
    @printf("  df/d(cap_length) = %.6e GHz/μm\n", dfdp)

    # ─── Step 4: Mesh-perturbation Taylor test ─────────────────────────────

    println("\nStep 4: Mesh-perturbation Taylor test...")
    h_values = [5e-2, 2e-2, 1e-2, 5e-3]
    taylor_errs = Float64[]

    for (i, h) in enumerate(h_values)
        # Write perturbed mesh
        perturbed_mesh_file = "mesh/transmon_taylor_$(i).mesh"
        mkpath("mesh")
        open(perturbed_mesh_file, "w") do io
            println(io, "MeshVersionFormatted 2\n\nDimension 3\n")
            println(io, "Vertices\n$nv")
            for vi in 1:nv
                x = palace_mesh.vertices[vi, 1] + h * dxdp_full[vi, 1]
                y = palace_mesh.vertices[vi, 2] + h * dxdp_full[vi, 2]
                z = palace_mesh.vertices[vi, 3] + h * dxdp_full[vi, 3]
                @printf(io, "%.15e %.15e %.15e 0\n", x, y, z)
            end
            _write_elems(io, "Tetrahedra", palace_mesh.tetrahedra, palace_mesh.tet_refs, 4)
            _write_elems(io, "Hexahedra", palace_mesh.hexahedra, palace_mesh.hex_refs, 8)
            _write_elems(io, "Triangles", palace_mesh.triangles, palace_mesh.tri_refs, 3)
            _write_elems(io, "Quadrilaterals", palace_mesh.quadrilaterals, palace_mesh.quad_refs, 4)
            println(io, "\nEnd")
        end

        # Write config for perturbed mesh
        config = JSON.parsefile("transmon.json")
        config["Model"]["Mesh"] = perturbed_mesh_file
        config["Problem"]["Output"] = "postpro/taylor_$(i)"
        config["Solver"]["Eigenmode"]["Sensitivity"] = false
        config["Solver"]["Eigenmode"]["MeshSensitivity"] = false
        taylor_config = "transmon_taylor_$(i).json"
        open(taylor_config, "w") do f
            JSON.print(f, config, 2)
            println(f)
        end

        mkpath("postpro/taylor_$(i)")
        run(pipeline(`$palace $taylor_config`; stdout=devnull, stderr=devnull))

        f_h = read_eigenfrequency("postpro/taylor_$(i)/eig.csv"; mode=1)
        linear_pred = f0 + h * dfdp
        err = abs(f_h - linear_pred)
        push!(taylor_errs, err)
        @printf("  ε=%.0e: f=%.8f, pred=%.8f, err=%.2e, err/ε²=%.2e\n",
                h, f_h, linear_pred, err, err / h^2)

        if i > 1
            rate = log(taylor_errs[i-1] / taylor_errs[i]) / log(h_values[i-1] / h_values[i])
            @printf("  → rate = %.2f\n", rate)
        end

        # Cleanup
        rm("postpro/taylor_$(i)"; recursive=true, force=true)
        rm(taylor_config; force=true)
        rm(perturbed_mesh_file; force=true)
    end

    # ─── Summary ───────────────────────────────────────────────────────────

    rates = Float64[]
    for i in 2:length(taylor_errs)
        push!(rates, log(taylor_errs[i-1] / taylor_errs[i]) / log(h_values[i-1] / h_values[i]))
    end

    println("\n═══ Summary ═══")
    @printf("  f₀ = %.6f GHz at cap_length = 620 μm\n", f0)
    @printf("  df/d(cap_length) = %.6e GHz/μm\n", dfdp)
    if !isempty(rates)
        @printf("  Taylor rates = %s\n", join([@sprintf("%.2f", r) for r in rates], ", "))
        max_rate = maximum(rates)
        @printf("  Status: %s (max rate %.2f)\n", max_rate > 1.5 ? "PASS ✓" : "FAIL ✗", max_rate)
    end
    println("\nDone.")
end

function _write_elems(io, keyword, elems, refs, ncols)
    ne = size(elems, 1)
    ne > 0 || return
    println(io, "\n$keyword\n$ne")
    for e in 1:ne
        for c in 1:ncols
            print(io, elems[e, c], " ")
        end
        println(io, refs[e])
    end
end

main()
