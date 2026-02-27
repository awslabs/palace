# Shared utilities for reading Palace mesh sensitivity output (.mesh + .solb)
# and propagating boundary sensitivities to interior nodes via Laplacian smoothing.
#
# Used by: chain_test.jl, chain_test_egads.jl, egads_pipeline_test.jl,
#          boolean_sensitivity_test.jl

using LinearAlgebra, SparseArrays

# ─── Minimal .mesh reader ────────────────────────────────────────────────────

struct SimpleMesh
    vertices::Matrix{Float64}   # Nx3
    hexahedra::Matrix{Int}      # Mx8 (1-based)
    hex_refs::Vector{Int}
    quadrilaterals::Matrix{Int}  # Mx4 (1-based)
    quad_refs::Vector{Int}
    tetrahedra::Matrix{Int}     # Mx4 (1-based)
    tet_refs::Vector{Int}
    triangles::Matrix{Int}      # Mx3 (1-based)
    tri_refs::Vector{Int}
end

function read_mesh(filename::AbstractString)
    lines = readlines(filename)
    vertices = zeros(Float64, 0, 3)
    hexahedra = zeros(Int, 0, 8);
    hex_refs = Int[]
    quadrilaterals = zeros(Int, 0, 4);
    quad_refs = Int[]
    tetrahedra = zeros(Int, 0, 4);
    tet_refs = Int[]
    triangles = zeros(Int, 0, 3);
    tri_refs = Int[]

    i = 1
    while i <= length(lines)
        line = strip(lines[i])
        isempty(line) && (i += 1; continue)
        keyword = first(split(line))

        if keyword == "Vertices"
            i += 1
            n = parse(Int, strip(lines[i]));
            i += 1
            vertices = zeros(Float64, n, 3)
            for j = 1:n
                tokens = split(strip(lines[i]))
                for d = 1:3
                    vertices[j, d] = parse(Float64, tokens[d])
                end
                i += 1
            end
        elseif keyword == "Hexahedra"
            hexahedra, hex_refs, i = _read_elem_section(lines, i, 8)
        elseif keyword == "Quadrilaterals"
            quadrilaterals, quad_refs, i = _read_elem_section(lines, i, 4)
        elseif keyword == "Tetrahedra"
            tetrahedra, tet_refs, i = _read_elem_section(lines, i, 4)
        elseif keyword == "Triangles"
            triangles, tri_refs, i = _read_elem_section(lines, i, 3)
        elseif keyword == "End"
            break
        else
            # Skip unknown keyword: try count + that many lines
            i += 1
            if i <= length(lines)
                count = tryparse(Int, strip(lines[i]))
                if !isnothing(count) && count >= 0
                    i += 1 + count
                end
            end
        end
    end
    return SimpleMesh(
        vertices,
        hexahedra,
        hex_refs,
        quadrilaterals,
        quad_refs,
        tetrahedra,
        tet_refs,
        triangles,
        tri_refs
    )
end

function _read_elem_section(lines, i, ncols)
    i += 1  # past keyword
    n = parse(Int, strip(lines[i]));
    i += 1
    elems = zeros(Int, n, ncols)
    refs = zeros(Int, n)
    for j = 1:n
        tokens = split(strip(lines[i]))
        for c = 1:ncols
            elems[j, c] = parse(Int, tokens[c])
        end
        if length(tokens) > ncols
            refs[j] = parse(Int, tokens[ncols + 1])
        end
        i += 1
    end
    return elems, refs, i
end

# ─── Minimal .solb reader ────────────────────────────────────────────────────

struct SolField
    npts::Int
    dim::Int
    data::Matrix{Float64}  # npts × dim
end

function read_solb(filename::AbstractString)
    lines = readlines(filename)
    i = 1
    dimension = 3
    while i <= length(lines)
        line = strip(lines[i])
        isempty(line) && (i += 1; continue)
        keyword = first(split(line))

        if keyword == "Dimension"
            rest = strip(line[(length("Dimension") + 1):end])
            dimension = isempty(rest) ? parse(Int, strip(lines[i += 1])) : parse(Int, rest)
            i += 1
        elseif keyword == "SolAtVertices"
            i += 1
            header = split(strip(lines[i]))
            npts = parse(Int, header[1])
            nfields = parse(Int, header[2])
            field_type = parse(Int, header[3])
            # field_type 2 = vector → dimension components
            ncomp =
                field_type == 1 ? 1 :
                field_type == 2 ? dimension : dimension * (dimension + 1) ÷ 2
            i += 1
            data = zeros(Float64, npts, ncomp)
            for p = 1:npts
                tokens = split(strip(lines[i]))
                for c = 1:ncomp
                    data[p, c] = parse(Float64, tokens[c])
                end
                i += 1
            end
            return SolField(npts, ncomp, data)
        elseif keyword == "End"
            break
        else
            i += 1
        end
    end
    return error("No SolAtVertices section found in $filename")
end

# ─── Read eigenfrequency from Palace eig.csv ─────────────────────────────────

function read_eigenfrequency(eig_csv::AbstractString; mode::Int=1)
    lines = readlines(eig_csv)
    for line in lines
        stripped = strip(line)
        startswith(stripped, "m,") && continue
        isempty(stripped) && continue
        tokens = split(stripped, ",")
        m = parse(Float64, strip(tokens[1]))
        if round(Int, m) == mode
            freq_ghz = parse(Float64, strip(tokens[2]))
            return freq_ghz  # GHz (same units as mesh sensitivity)
        end
    end
    return error("Mode $mode not found in $eig_csv")
end

# ─── Laplacian boundary → interior sensitivity propagation ────────────────────

"""
    build_adjacency(mesh::SimpleMesh) -> Vector{Set{Int}}

Build vertex adjacency graph from element connectivity.
"""
function build_adjacency(mesh::SimpleMesh)
    nv = size(mesh.vertices, 1)
    adj = [Set{Int}() for _ = 1:nv]

    function add_element_edges!(conn::Matrix{Int})
        for e = 1:size(conn, 1)
            verts = conn[e, :]
            for i in eachindex(verts), j in eachindex(verts)
                i == j && continue
                push!(adj[verts[i]], verts[j])
            end
        end
    end

    size(mesh.hexahedra, 1) > 0 && add_element_edges!(mesh.hexahedra)
    size(mesh.tetrahedra, 1) > 0 && add_element_edges!(mesh.tetrahedra)

    return adj
end

"""
    classify_boundary_vertices(mesh::SimpleMesh) -> BitVector

Returns a BitVector where true = boundary vertex (appears in surface elements).
"""
function classify_boundary_vertices(mesh::SimpleMesh)
    nv = size(mesh.vertices, 1)
    is_boundary = falses(nv)

    for elems in (mesh.quadrilaterals, mesh.triangles)
        for e = 1:size(elems, 1)
            for c = 1:size(elems, 2)
                is_boundary[elems[e, c]] = true
            end
        end
    end
    return is_boundary
end

"""
    propagate_boundary_sensitivity(mesh::SimpleMesh, boundary_dxdp::Matrix{Float64},
                                   is_boundary::BitVector) -> Matrix{Float64}

Propagate boundary dx/dp values to interior vertices via Laplacian smoothing.

The Laplacian system preserves linear functions on convex domains, so for the
cylinder radius change (dx/dR = x/R, dy/dR = y/R), interior values should match
the analytical solution.
"""
function propagate_boundary_sensitivity(
    mesh::SimpleMesh,
    boundary_dxdp::Matrix{Float64},
    is_boundary::BitVector
)
    nv = size(mesh.vertices, 1)
    dim = size(boundary_dxdp, 2)
    adj = build_adjacency(mesh)

    # Build sparse Laplacian with boundary rows as identity
    I_idx = Int[];
    J_idx = Int[];
    V_val = Float64[]

    for i = 1:nv
        if is_boundary[i]
            # Boundary: identity row
            push!(I_idx, i);
            push!(J_idx, i);
            push!(V_val, 1.0)
        else
            # Interior: Laplacian row — L[i,i] = degree, L[i,j] = -1
            neighbors = adj[i]
            deg = length(neighbors)
            push!(I_idx, i);
            push!(J_idx, i);
            push!(V_val, Float64(deg))
            for j in neighbors
                push!(I_idx, i);
                push!(J_idx, j);
                push!(V_val, -1.0)
            end
        end
    end

    L = sparse(I_idx, J_idx, V_val, nv, nv)

    # Solve for each spatial dimension
    result = zeros(Float64, nv, dim)
    for d = 1:dim
        rhs = zeros(Float64, nv)
        for i = 1:nv
            if is_boundary[i]
                rhs[i] = boundary_dxdp[i, d]
            end
        end
        result[:, d] = L \ rhs
    end
    return result
end
