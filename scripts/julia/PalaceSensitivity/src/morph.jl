# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Morphing of MFEM meshes: the nodes of a mesh moved by a displacement, such as a
# combination of design velocities, and the rest of the file kept.

"""
    struct MorphQuality

The change in element quality from a [`morph_mesh`](@ref):

  - `num_elements`: the number of elements checked: the volume elements, or the surface
    elements of a mesh without volume elements. Each is checked on the simplices at its
    corners, which its vertices define.
  - `num_inverted`: the number of checked elements that the morph inverted at some corner.
  - `min_jacobian_ratio`, `max_jacobian_ratio`: the extreme ratios of the volume of a corner
    simplex (the area, for surface elements) after the morph to its volume before.
  - `max_strain`: the largest change in an edge of a corner simplex relative to its length,
    an estimate of the largest displacement gradient.
  - `max_displacement`: the largest node displacement, in mesh units.

A small `min_jacobian_ratio` (below about 0.3) or a large `max_strain` (above about 0.3)
indicates that the morph has degraded the mesh, which should be replaced by a new one.
"""
struct MorphQuality
    num_elements::Int
    num_inverted::Int
    min_jacobian_ratio::Float64
    max_jacobian_ratio::Float64
    max_strain::Float64
    max_displacement::Float64
end

"""
    morph_mesh(src, dst, displacement; length_unit=1e-6, velocity_length_unit=1e-6,
               check=true)
    morph_mesh(src, dst, velocities, δp::NamedTuple; kwargs...)

Write to `dst` the MFEM mesh in the file `src` with each node `x` moved to
`x + displacement(x)`, and return the [`MorphQuality`](@ref) of the result. Only the node
coordinates change.

`displacement` maps a `3 × N` matrix of node coordinates to a `3 × N` matrix of
displacements, both with lengths in `velocity_length_unit` meters, and `length_unit` is the
mesh unit in meters (`Model.L0`). The second form moves the nodes by
`Σₖ δp[k] velocities[k](x)`, for velocities as in [`jacobian`](@ref) and parameter changes
`δp` in the units of the parameters, with lengths in `velocity_length_unit` meters.

The mesh must be an uncompressed MFEM mesh v1.x in three dimensions, with vertices or with
nodes. With `check=true`, an error is thrown if the morph inverts any element at any corner,
and nothing is written. The quality of a mesh with nodes is checked at its vertices, which
requires the nodes to be an `H1` grid function, whose first degrees of freedom are the
vertices; otherwise the quality is `nothing`, and `check=true` is an error.
"""
function morph_mesh(
    src::AbstractString,
    dst::AbstractString,
    displacement;
    length_unit=1e-6,
    velocity_length_unit=1e-6,
    check=true
)
    mesh = _read_mfem(src)
    scale = length_unit / velocity_length_unit # mesh unit in velocity units
    D = displacement(mesh.nodes .* scale) ./ scale
    size(D) == size(mesh.nodes) || throw(
        DimensionMismatch(
            "the displacement of $(size(mesh.nodes, 2)) nodes has size $(size(D))"
        )
    )
    newnodes = mesh.nodes .+ D
    quality = _morph_quality(mesh, mesh.nodes, newnodes)
    if check
        isnothing(quality) && error(
            "Cannot check the quality of the morphed mesh $src, whose nodes are not an H1 grid function; pass check=false to morph without checking"
        )
        quality.num_inverted > 0 && error(
            "The morph inverts $(quality.num_inverted) of $(quality.num_elements) elements of $src (minimum Jacobian ratio $(quality.min_jacobian_ratio)); the displacement is too large for this mesh"
        )
    end
    write(dst, _morphed_bytes(mesh, newnodes))
    return quality
end

function morph_mesh(
    src::AbstractString,
    dst::AbstractString,
    velocities,
    δp::NamedTuple;
    kwargs...
)
    for (k, v) in pairs(δp)
        haskey(velocities, k) || throw(ArgumentError("no velocity for parameter :$k"))
        v isa Real || throw(
            ArgumentError(
                "the change of :$k must be a number, in the unit of the parameter (lengths in velocity_length_unit meters), got $v"
            )
        )
    end
    displacement(X) = sum(δp[k] .* velocities[k](X) for k in keys(δp); init=zero(X))
    return morph_mesh(src, dst, displacement; kwargs...)
end

"""
    morph_quality(reference, mesh)

The [`MorphQuality`](@ref) of the MFEM mesh file `mesh` relative to the mesh file
`reference`, which must have the same nodes and elements: for example, a mesh after a
sequence of morphs relative to the mesh before the first, to decide when to replace it.
"""
function morph_quality(reference::AbstractString, mesh::AbstractString)
    a, b = _read_mfem(reference), _read_mfem(mesh)
    (size(a.nodes) == size(b.nodes) && a.cells == b.cells) ||
        error("$mesh does not have the nodes and elements of $reference")
    quality = _morph_quality(a, a.nodes, b.nodes)
    isnothing(quality) && error(
        "Cannot check the quality of $reference, whose nodes are not an H1 grid function"
    )
    return quality
end

######## Reading and writing

# An MFEM mesh file: its bytes, its node coordinates (3 × N), where they are in the file,
# and the corners of the elements to check.
struct _MFEMMesh
    bytes::Vector{UInt8}
    nodes::Matrix{Float64}
    # The byte range of the coordinates of each node, or of each value for nodes by node.
    locations::Vector{UnitRange{Int}}
    by_nodes::Bool
    cells::Vector{Tuple{Symbol, Vector{Int}}}
    checkable::Bool # false for nodes that are not an H1 grid function
end

# The number of vertices of each MFEM geometry type, and the shapes checked for quality.
const _MFEM_VERTICES = Dict(0 => 1, 1 => 2, 2 => 3, 3 => 4, 4 => 4, 5 => 8, 6 => 6, 7 => 5)
const _MFEM_SHAPES =
    Dict(2 => :tri, 3 => :quad, 4 => :tet, 5 => :hex, 6 => :prism, 7 => :pyramid)

# A cursor over the lines of a file's bytes.
mutable struct _Cursor
    bytes::Vector{UInt8}
    pos::Int
end

function _readline!(c::_Cursor)
    i = findnext(==(UInt8('\n')), c.bytes, c.pos)
    stop = isnothing(i) ? length(c.bytes) : i
    r = (c.pos):(isnothing(i) ? stop : stop - 1)
    c.pos = stop + 1
    # Tolerate Windows line endings.
    if !isempty(r) && c.bytes[last(r)] == UInt8('\r')
        r = first(r):(last(r) - 1)
    end
    return r
end

_line(c::_Cursor, r) = String(c.bytes[r])

function _read_mfem(path::AbstractString)
    bytes = read(path)
    length(bytes) >= 2 &&
        bytes[1] == 0x1f &&
        bytes[2] == 0x8b &&
        error(
            "$path is gzip-compressed; decompress it (for example with gunzip) to morph it"
        )
    startswith(String(bytes[1:min(end, 64)]), "MFEM mesh v1.") ||
        error("Unsupported mesh format in $path (expected an MFEM mesh v1.x)")
    c = _Cursor(bytes, 1)
    cells = Tuple{Symbol, Vector{Int}}[]
    nodes = Matrix{Float64}(undef, 3, 0)
    locations = UnitRange{Int}[]
    by_nodes, h1, has_nodes = false, true, false
    while c.pos <= length(bytes)
        line = strip(_line(c, _readline!(c)))
        if line == "elements"
            for _ = 1:parse(Int, strip(_line(c, _readline!(c))))
                fields = parse.(Int, split(_line(c, _readline!(c))))
                geom = fields[2]
                haskey(_MFEM_VERTICES, geom) ||
                    error("Unsupported MFEM geometry type $geom in $path")
                haskey(_MFEM_SHAPES, geom) &&
                    push!(cells, (_MFEM_SHAPES[geom], fields[3:end] .+ 1))
            end
        elseif line == "vertices"
            nv = parse(Int, strip(_line(c, _readline!(c))))
            # Without nodes, the vertex coordinates follow their dimension.
            save = c.pos
            l = strip(_line(c, _readline!(c)))
            if l == "nodes" || isempty(l)
                c.pos = save
            else
                sdim = parse(Int, l)
                sdim == 3 || error("Unsupported MFEM mesh space dimension $sdim in $path")
                nodes = Matrix{Float64}(undef, 3, nv)
                locations = Vector{UnitRange{Int}}(undef, nv)
                for j = 1:nv
                    r = _readline!(c)
                    fields = split(_line(c, r))
                    for d = 1:3
                        nodes[d, j] = parse(Float64, fields[d])
                    end
                    locations[j] = r
                end
            end
        elseif line == "nodes"
            has_nodes = true
            vdim, ordering = 0, -1
            while true
                l = strip(_line(c, _readline!(c)))
                isempty(l) && vdim > 0 && ordering >= 0 && break
                if startswith(l, "FiniteElementCollection:")
                    h1 = startswith(
                        strip(l[(length("FiniteElementCollection:") + 1):end]),
                        "H1"
                    )
                elseif startswith(l, "VDim:")
                    vdim = parse(Int, strip(l[6:end]))
                elseif startswith(l, "Ordering:")
                    ordering = parse(Int, strip(l[10:end]))
                end
            end
            vdim == 3 || error("Unsupported MFEM nodes dimension $vdim in $path")
            by_nodes = ordering == 0
            values, ranges = Float64[], UnitRange{Int}[]
            while c.pos <= length(bytes)
                save = c.pos
                r = _readline!(c)
                fields = split(_line(c, r))
                isempty(fields) && continue
                v = tryparse(Float64, fields[1])
                if isnothing(v)
                    c.pos = save
                    break
                end
                by_nodes ? push!(values, v) : append!(values, parse.(Float64, fields[1:3]))
                push!(ranges, r)
            end
            nodes = by_nodes ? permutedims(reshape(values, :, 3)) : reshape(values, 3, :)
            locations = ranges
        end
    end
    isempty(locations) && error("No vertex or node coordinates in $path")
    # The volume elements are checked, or the surface elements if there are none.
    any(c -> _isvolume(c[1]), cells) && filter!(c -> _isvolume(c[1]), cells)
    return _MFEMMesh(bytes, nodes, locations, by_nodes, cells, !has_nodes || h1)
end

# The file's bytes with the coordinate text replaced by that of `newnodes`.
function _morphed_bytes(mesh::_MFEMMesh, newnodes)
    io = IOBuffer()
    pos = 1
    n = size(newnodes, 2)
    for (k, r) in enumerate(mesh.locations)
        write(io, view(mesh.bytes, pos:(first(r) - 1)))
        if mesh.by_nodes
            d, j = divrem(k - 1, n)
            print(io, newnodes[d + 1, j + 1])
        else
            print(io, newnodes[1, k], ' ', newnodes[2, k], ' ', newnodes[3, k])
        end
        pos = last(r) + 1
    end
    write(io, view(mesh.bytes, pos:length(mesh.bytes)))
    return take!(io)
end

######## Quality

# The simplices at the corners of each element shape: each corner with its neighbors, by
# local vertex index in MFEM's order. An element is inverted when the volume (area) of any
# of them changes sign.
const _CORNER_SIMPLICES = Dict{Symbol, Vector{Vector{Int}}}(
    :tet => [[1, 2, 3, 4]],
    :hex => [
        [1, 2, 4, 5],
        [2, 3, 1, 6],
        [3, 4, 2, 7],
        [4, 1, 3, 8],
        [5, 8, 6, 1],
        [6, 5, 7, 2],
        [7, 6, 8, 3],
        [8, 7, 5, 4]
    ],
    :prism => [
        [1, 2, 3, 4],
        [2, 3, 1, 5],
        [3, 1, 2, 6],
        [4, 6, 5, 1],
        [5, 4, 6, 2],
        [6, 5, 4, 3]
    ],
    :pyramid => [[1, 2, 4, 5], [2, 3, 1, 5], [3, 4, 2, 5], [4, 1, 3, 5]],
    :tri => [[1, 2, 3]],
    :quad => [[1, 2, 4], [2, 3, 1], [3, 4, 2], [4, 1, 3]]
)
_isvolume(shape) = shape in (:tet, :hex, :prism, :pyramid)

function _morph_quality(mesh::_MFEMMesh, X, Y)
    mesh.checkable || return nothing
    dmax = maximum(j -> _norm3(_sub(_point(Y, j), _point(X, j))), axes(X, 2); init=0.0)
    rmin, rmax, smax, ninv = Inf, -Inf, 0.0, 0
    for (shape, corners) in mesh.cells
        inverted = false
        for simplex in _CORNER_SIMPLICES[shape]
            c = view(corners, simplex)
            r = length(c) == 4 ? _volume_ratio(X, Y, c) : _area_ratio(X, Y, c)
            isnan(r) && continue
            rmin, rmax = min(rmin, r), max(rmax, r)
            inverted |= r <= 0
            for a = 1:length(c), b = (a + 1):length(c)
                e = _sub(_point(X, c[b]), _point(X, c[a]))
                f = _sub(_point(Y, c[b]), _point(Y, c[a]))
                L = _norm3(e)
                L > 0 && (smax = max(smax, _norm3(_sub(f, e)) / L))
            end
        end
        ninv += inverted
    end
    isempty(mesh.cells) && ((rmin, rmax) = (1.0, 1.0))
    return MorphQuality(length(mesh.cells), ninv, rmin, rmax, smax, dmax)
end

function _volume_ratio(X, Y, c)
    volume(X) = _dot3(
        _sub(_point(X, c[2]), _point(X, c[1])),
        _cross3(
            _sub(_point(X, c[3]), _point(X, c[1])),
            _sub(_point(X, c[4]), _point(X, c[1]))
        )
    )
    v0 = volume(X)
    iszero(v0) && return NaN
    return volume(Y) / v0
end

function _area_ratio(X, Y, c)
    normal(X) = _cross3(
        _sub(_point(X, c[2]), _point(X, c[1])),
        _sub(_point(X, c[3]), _point(X, c[1]))
    )
    n0 = normal(X)
    a0 = _dot3(n0, n0)
    iszero(a0) && return NaN
    return _dot3(normal(Y), n0) / a0
end

_point(X, i) = (X[1, i], X[2, i], X[3, i])
_sub(a, b) = (a[1] - b[1], a[2] - b[2], a[3] - b[3])
_dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
_cross3(a, b) =
    (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
_norm3(a) = sqrt(_dot3(a, a))
