#!/usr/bin/env julia
#
# End-to-end sensitivity chain orchestration script.
#
# Computes the full design sensitivity dJ/dp by contracting Palace's mesh sensitivity
# dJ/dx (from P6) with DeviceLayout/EGADS's geometry sensitivity dx/dp (from DL6):
#
#   dJ/dp = sum_i (dJ/dx_i) * (dx_i/dp)
#
# where the sum runs over all mesh boundary nodes in each spatial direction.
#
# Usage:
#   julia sensitivity_chain.jl <palace_output_dir> <geometry_sens_dir> [--mapping <type>]
#
# Arguments:
#   palace_output_dir    Directory containing Palace outputs, including:
#                        - mesh_sensitivity.solb  (dJ/dx vector field, all mesh vertices)
#                        - metric0000.meshb       (mesh in .meshb format, for node mapping)
#
#   geometry_sens_dir    Directory containing EGADS geometry sensitivities:
#                        - dxdp_<param_name>.solb files (one per design parameter)
#                        - Optional: egads_tess.meshb   (tessellation mesh for node mapping)
#
# Options:
#   --mapping identity   Assume identical vertex ordering (default, for EGADS-generated mesh)
#   --mapping nearest    Use nearest-neighbor matching between vertex sets
#
# ─── Verification test case ──────────────────────────────────────────────────
#
# Box resonator: a rectangular cavity with dimensions W x H x D.
#
# The TE_mnp eigenfrequency of a perfectly conducting rectangular cavity is:
#
#   f_mnp = c / (2 * sqrt(eps * mu)) * sqrt((m/W)^2 + (n/H)^2 + (p/D)^2)
#
# where m, n, p are mode indices and c is the speed of light.
#
# Analytical sensitivity with respect to width W:
#
#   df/dW = -f * (m/W)^2 / ((m/W)^2 + (n/H)^2 + (p/D)^2) / W
#
# For the fundamental TE_101 mode (m=1, n=0, p=1) of a cube W = H = D = a:
#
#   f_101 = c / (2*a) * sqrt(2)
#   df/da = -c * sqrt(2) / (2 * a^2)  =  -f / a
#
# Verification procedure:
#   1. Run Palace eigenmode solve on a box with mesh sensitivity enabled
#   2. Run this script to contract with dx/dp from EGADS (parameter = width W)
#   3. Compare computed dJ/dW against the analytical formula above
#   4. Also verify via finite differences: (f(W+h) - f(W)) / h
#
# Expected agreement: adjoint chain should match FD to within O(h) (FD truncation
# error) + O(h_mesh) (mesh sensitivity FD error). For a well-resolved mesh with
# h_FD ~ 1e-7 * char_length, agreement within 1-5% is typical.
#
# ─── End-to-end pipeline ─────────────────────────────────────────────────────
#
# The full derivative propagation pipeline that this script connects:
#
#   1. Design parameters p (width, height, gap, ...)
#           |
#           v  [EGADS _dot API via DeviceLayout DL6]
#   2. Mesh node positions x(p) — via EGADS tessellation velocity dx/dp
#           |
#           v  [Palace P6 — mesh sensitivity via FD or analytic]
#   3. Objective J(x) — eigenfrequency, S-parameter, quality factor, ...
#           |
#           v  [Chain rule contraction — this script]
#   4. Design gradient dJ/dp = (dJ/dx) . (dx/dp)
#
# The gradient dJ/dp enables gradient-based optimization of EM device geometry
# using standard optimization algorithms (L-BFGS, trust-region, etc.).
#
# ─────────────────────────────────────────────────────────────────────────────

using LinearAlgebra

# ─── Standalone .solb reader ─────────────────────────────────────────────────
#
# A minimal .solb reader that does not require DeviceLayout. This allows the script
# to run standalone with only the Julia standard library. For production use with
# DeviceLayout loaded, the extension's read_sol is preferred.

const SOL_SCALAR = 1
const SOL_VECTOR = 2
const SOL_SYMMAT = 3

struct SolField
    field_type::Int
    data::Matrix{Float64}  # npts x ncomp
end

struct SolFile
    version::Int
    dimension::Int
    fields::Vector{SolField}
end

function sol_ncomp(field_type::Int, dimension::Int)
    field_type == SOL_SCALAR && return 1
    field_type == SOL_VECTOR && return dimension
    field_type == SOL_SYMMAT && return dimension * (dimension + 1) ÷ 2
    return error("Unknown field type: $field_type")
end

function read_solb(filename::AbstractString)
    lines = readlines(filename)
    version = 2
    dimension = 3
    fields = SolField[]

    i = 1
    while i <= length(lines)
        line = strip(lines[i])
        isempty(line) && (i += 1; continue)
        tokens = split(line)
        keyword = tokens[1]

        if keyword == "MeshVersionFormatted"
            version =
                length(tokens) > 1 ? parse(Int, tokens[2]) :
                parse(Int, strip(lines[i += 1]))
            i += 1
        elseif keyword == "Dimension"
            dimension =
                length(tokens) > 1 ? parse(Int, tokens[2]) :
                parse(Int, strip(lines[i += 1]))
            i += 1
        elseif keyword == "SolAtVertices"
            i += 1
            header = split(strip(lines[i]))
            npts = parse(Int, header[1])
            nfields = parse(Int, header[2])
            ftypes = [parse(Int, header[2 + f]) for f = 1:nfields]
            ncomps = [sol_ncomp(ft, dimension) for ft in ftypes]
            field_data = [zeros(Float64, npts, nc) for nc in ncomps]

            i += 1
            for p = 1:npts
                row_tokens = split(strip(lines[i]))
                col = 1
                for (fi, nc) in enumerate(ncomps)
                    for c = 1:nc
                        field_data[fi][p, c] = parse(Float64, row_tokens[col])
                        col += 1
                    end
                end
                i += 1
            end

            fields = [SolField(ft, fd) for (ft, fd) in zip(ftypes, field_data)]
        elseif keyword == "End"
            break
        else
            i += 1
        end
    end

    return SolFile(version, dimension, fields)
end

# ─── Standalone .meshb vertex reader ─────────────────────────────────────────

function read_meshb_vertices(filename::AbstractString)
    lines = readlines(filename)
    dimension = 3
    vertices = zeros(Float64, 0, 3)

    i = 1
    while i <= length(lines)
        line = strip(lines[i])
        isempty(line) && (i += 1; continue)
        tokens = split(line)
        keyword = tokens[1]

        if keyword == "Dimension"
            dimension =
                length(tokens) > 1 ? parse(Int, tokens[2]) :
                parse(Int, strip(lines[i += 1]))
            i += 1
        elseif keyword == "Vertices"
            i += 1
            n = parse(Int, strip(lines[i]))
            i += 1
            vertices = zeros(Float64, n, dimension)
            for j = 1:n
                row = split(strip(lines[i]))
                for d = 1:dimension
                    vertices[j, d] = parse(Float64, row[d])
                end
                i += 1
            end
        elseif keyword == "End"
            break
        else
            i += 1
        end
    end

    return vertices
end

# ─── Node mapping ────────────────────────────────────────────────────────────

"""
    nearest_neighbor_mapping(src::Matrix{Float64}, dst::Matrix{Float64}) -> Vector{Int}

For each row in `src`, find the index of the nearest row in `dst`.
Returns a mapping vector of length size(src, 1).
"""
function nearest_neighbor_mapping(src::Matrix{Float64}, dst::Matrix{Float64})
    n_src = size(src, 1)
    n_dst = size(dst, 1)
    dim = size(src, 2)
    size(dst, 2) == dim || error("Dimension mismatch")

    mapping = Vector{Int}(undef, n_src)
    max_dist = 0.0

    for i = 1:n_src
        best_j = 1
        best_d2 = Inf
        for j = 1:n_dst
            d2 = 0.0
            for d = 1:dim
                diff = src[i, d] - dst[j, d]
                d2 += diff * diff
            end
            if d2 < best_d2
                best_d2 = d2
                best_j = j
            end
        end
        mapping[i] = best_j
        max_dist = max(max_dist, sqrt(best_d2))
    end

    if max_dist > 1e-6
        @warn "Nearest-neighbor mapping: max distance = $max_dist (threshold 1e-6). " *
              "Meshes may be from different sources."
    end

    return mapping
end

# ─── Chain rule contraction ──────────────────────────────────────────────────

"""
    contract(dJdx::Matrix{Float64}, dxdp::Matrix{Float64}) -> Float64

Compute dJ/dp = sum(dJ/dx .* dx/dp) via dot product.
Both matrices must be Nx3 (same vertex count and dimension).
"""
function contract(dJdx::Matrix{Float64}, dxdp::Matrix{Float64})
    size(dJdx) == size(dxdp) ||
        error("Shape mismatch: dJ/dx is $(size(dJdx)), dx/dp is $(size(dxdp))")
    return dot(vec(dJdx), vec(dxdp))
end

"""
    run_sensitivity_chain(palace_dir, geom_sens_dir; mapping_type="identity")

Run the full sensitivity chain and print results.

# Arguments

  - `palace_dir`: Palace output directory (contains mesh_sensitivity.solb)
  - `geom_sens_dir`: Directory with dxdp_*.solb files from DeviceLayout
  - `mapping_type`: "identity" or "nearest"
"""
function run_sensitivity_chain(
    palace_dir::AbstractString,
    geom_sens_dir::AbstractString;
    mapping_type::AbstractString="identity"
)
    # Read Palace's dJ/dx
    dJdx_file = joinpath(palace_dir, "mesh_sensitivity.solb")
    isfile(dJdx_file) || error("Palace dJ/dx file not found: $dJdx_file")
    dJdx_sol = read_solb(dJdx_file)
    isempty(dJdx_sol.fields) && error("dJ/dx .solb has no fields")
    dJdx_sol.fields[1].field_type == SOL_VECTOR ||
        error("dJ/dx field is not a vector (type=$(dJdx_sol.fields[1].field_type))")
    dJdx = dJdx_sol.fields[1].data  # Nx3

    println("Palace dJ/dx: $(size(dJdx, 1)) vertices, $(dJdx_sol.dimension)D")
    println("  ||dJ/dx|| = $(norm(vec(dJdx)))")
    println("  nonzero entries: $(count(!iszero, dJdx)) / $(length(dJdx))")

    # Build node mapping if needed
    node_map = nothing
    if mapping_type == "nearest"
        # Need vertex coordinates from both meshes
        palace_meshb = joinpath(palace_dir, "metric0000.meshb")
        egads_meshb = joinpath(geom_sens_dir, "egads_tess.meshb")
        if !isfile(palace_meshb)
            @warn "Palace .meshb not found at $palace_meshb, trying mesh.meshb"
            palace_meshb = joinpath(palace_dir, "mesh.meshb")
        end
        isfile(palace_meshb) || error("No Palace .meshb file found for node mapping")
        isfile(egads_meshb) ||
            error("No EGADS .meshb file found at $egads_meshb for node mapping")

        palace_verts = read_meshb_vertices(palace_meshb)
        egads_verts = read_meshb_vertices(egads_meshb)
        node_map = nearest_neighbor_mapping(palace_verts, egads_verts)
        println(
            "Node mapping: nearest-neighbor, $(length(node_map)) Palace -> $(size(egads_verts, 1)) EGADS vertices"
        )
    else
        println("Node mapping: identity (same mesh source)")
    end

    # Find all dx/dp files
    param_files = filter(readdir(geom_sens_dir; join=true)) do f
        return endswith(f, ".solb") && startswith(basename(f), "dxdp_")
    end

    if isempty(param_files)
        println("\nNo dxdp_*.solb files found in $geom_sens_dir")
        return Dict{String, Float64}()
    end

    println("\nFound $(length(param_files)) design parameter(s):")

    results = Dict{String, Float64}()

    for pfile in param_files
        param_name = replace(basename(pfile), "dxdp_" => "", ".solb" => "")
        dxdp_sol = read_solb(pfile)
        isempty(dxdp_sol.fields) && continue
        dxdp_sol.fields[1].field_type == SOL_VECTOR || continue
        dxdp = dxdp_sol.fields[1].data  # Mx3

        # Apply node mapping: rearrange dJ/dx to match dx/dp ordering
        if !isnothing(node_map)
            npts_palace = size(dJdx, 1)
            npts_geom = size(dxdp, 1)
            dim = size(dJdx, 2)
            dJdx_mapped = zeros(Float64, npts_geom, dim)
            for i = 1:npts_palace
                j = node_map[i]
                if 1 <= j <= npts_geom
                    for d = 1:dim
                        dJdx_mapped[j, d] += dJdx[i, d]
                    end
                end
            end
            dJdp = contract(dJdx_mapped, dxdp)
        else
            dJdp = contract(dJdx, dxdp)
        end

        results[param_name] = dJdp
        println("  $param_name: dJ/dp = $dJdp  (||dx/dp|| = $(norm(vec(dxdp))))")
    end

    # Write results to file
    output_file = joinpath(palace_dir, "sensitivity_chain.txt")
    open(output_file, "w") do io
        println(io, "# Full sensitivity chain: dJ/dp = (dJ/dx) . (dx/dp)")
        println(io, "# Palace dJ/dx vertices: $(size(dJdx, 1))")
        println(io, "# Node mapping: $mapping_type")
        println(io, "#")
        println(io, "# Parameter  dJ/dp")
        for (param, val) in sort(collect(results))
            println(io, "$param  $val")
        end
    end
    println("\nResults written to $output_file")

    return results
end

# ─── Main ────────────────────────────────────────────────────────────────────

function main()
    if length(ARGS) < 2
        println("""
        Usage: julia sensitivity_chain.jl <palace_output_dir> <geometry_sens_dir> [options]

        Options:
          --mapping identity   Assume identical vertex ordering (default)
          --mapping nearest    Use nearest-neighbor matching between vertex sets

        Palace output dir should contain:
          - mesh_sensitivity.solb  (dJ/dx from Palace P6)

        Geometry sensitivity dir should contain:
          - dxdp_<param>.solb     (dx/dp from DeviceLayout/EGADS DL6)
          - (optional) egads_tess.meshb  (for nearest-neighbor node mapping)
        """)
        return
    end

    palace_dir = ARGS[1]
    geom_sens_dir = ARGS[2]
    mapping_type = "identity"

    # Parse optional arguments
    i = 3
    while i <= length(ARGS)
        if ARGS[i] == "--mapping" && i + 1 <= length(ARGS)
            mapping_type = ARGS[i + 1]
            i += 2
        else
            @warn "Unknown argument: $(ARGS[i])"
            i += 1
        end
    end

    return run_sensitivity_chain(palace_dir, geom_sens_dir; mapping_type=mapping_type)
end

# Run if executed as a script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
