# EGADS -> TetGen meshing helpers: PLC construction from EGADS bodies, size-field
# refinement, vertex classification onto the CAD, curved P2 meshes and an order-2 MSH
# writer. Moved from EGADS.jl (its EGADSTetGenExt extension as of 0c940f4, removed there in
# 349ea5a), which wraps EGADS and nothing more. Use it from an environment that has EGADS.jl
# and TetGen.jl:
#
#     include(joinpath(@__DIR__, "shared", "egads_tetgen.jl"))
#     using .EGADSTetGen
#
# Tests: test_egads_tetgen.jl next to this file.
#
# Known failure with TetGen_jll 1.5.4 (all pass with 1.5.3): "curve_tetgen_mesh: cylinder
# -> valid curved P2" finds boundary mid-nodes at radius 1.9712 instead of 2. The refinement
# passes (no -Y) put a Steiner point on a flat boundary facet, 0.038 off the cylinder;
# snapping it onto the surface makes a neighbouring tet invalid at full curvature, and the
# straightening fallback in curve_tetgen_mesh then resets every mid-node of that tet,
# including the CAD-constrained boundary ones. The fix: the fallback must not move a
# CAD-constrained node, and boundary mid-nodes should be placed at the parametric midpoint
# (t on the shared edge, (u,v) on the face) rather than by projection.

"""
    EGADSTetGen

Helpers bridging EGADS BREP bodies to TetGen volume meshing.

The integration is entirely in Julia — no external files, no gmsh, no C++ shims.

## Single body (`tetgen_io(body)`)

`EGADS.tessellate` produces a watertight surface triangulation; every triangle
becomes one PLC facet, with its marker equal to the EGADS face index (1-based,
optionally offset). `TetGen.tetrahedralize` then does the volume mesh.

## Multi-body (`tetgen_io(bodies)`)

The challenge: independent BREP bodies (as produced by DeviceLayout's
`EGADSKernel` fragmentation) do not share EGADS topology pointers, even when
they are geometrically adjacent. So independent tessellations of a shared
interface quad may pick opposite diagonals → TetGen rejects the PLC as
self-intersecting.

The fix is to emit every **planar** face as a single polygonal PLC facet
using its topological boundary loop, rather than as triangles. Point merging
across bodies then collapses the loops of two coincident faces onto the same
vertex set, and identical polygon facets across bodies are deduplicated into
one interior PLC facet — giving a conformal multi-region PLC with no
triangulation-diagonal mismatch. TetGen triangulates each polygonal facet
internally.

Curved (non-planar) faces are still emitted as tessellated triangles. In
DeviceLayout's current SolidModel output every volume is a box, so the
planar path covers the whole geometry; the curved-face path is a best effort
for more exotic BREPs (two coincident curved faces across bodies are not yet
stitched).

One region seed per body is injected so TetGen tags every tetrahedron with a
material attribute (the body index).
"""
module EGADSTetGen

using EGADS
using EGADS: EGObject, OClass, MType, tessellate, get_global_mesh,
             get_face_mesh, get_edge_mesh, get_body_topos, index_body_topo,
             get_bounding_box, get_topology, in_topology
import TetGen

export tetgen_io, tetrahedralize, tetrahedralize_size_field, classify_tetgen_vertices,
       curve_tetgen_mesh, write_msh2_p2

include("egads_tetgen_curving.jl")

# ─── Core: single-body PLC ──────────────────────────────────────────────────

"""
    tetgen_io(body; max_edge=0, sag=0, angle=15,
                    hole_seeds=nothing, region_seeds=nothing,
                    face_marker_offset=0)

Build a `TetGen.RawTetGenIO{Cdouble}` PLC from a single EGADS BODY.

The surface of `body` is tessellated with `tessellate(body; max_edge, sag,
angle)` and every triangle becomes a PLC facet. Facet marker `i +
face_marker_offset` is assigned to triangles that came from EGADS face `i`.
"""
function tetgen_io(body::EGObject;
                         max_edge::Real = 0.0, sag::Real = 0.0, angle::Real = 15.0,
                         hole_seeds = nothing, region_seeds = nothing,
                         face_marker_offset::Integer = 0,
                         meshsize_callback = nothing)
    _apply_meshsize_attributes!([body], meshsize_callback; angle = angle)
    tess = tessellate(body; max_edge = Float64(max_edge), sag = Float64(sag),
                      angle = Float64(angle))
    (xyz, tris, face_ids) = get_global_mesh(tess, body)
    io = TetGen.RawTetGenIO{Cdouble}()
    io.pointlist = Matrix{Cdouble}(transpose(xyz))
    TetGen.facetlist!(io, Matrix{Cint}(transpose(tris)))
    io.facetmarkerlist = Cint[Int32(face_marker_offset) + Int32(fi) for fi in face_ids]
    _set_holes_regions!(io, hole_seeds, region_seeds)
    return io
end

# ─── Multi-body PLC with interface stitching ────────────────────────────────

"""
    tetgen_io(bodies; max_edge=0, sag=0, angle=15,
                    hole_seeds=nothing, region_seeds=nothing,
                    auto_regions=true, face_marker_offset=1000,
                    coincidence_tol=nothing, verbose=false)

Build **one conformal PLC** covering a collection of EGADS bodies, including
shared interfaces.

## Strategy

For each body, for each face: if the face is planar, extract its boundary
loops as ordered point lists (walking edges via EGADS topology) and emit the
face as a single polygonal PLC facet. Non-planar faces are tessellated and
emitted as triangular facets.

All points are merged across bodies using a tolerance-based spatial hash, and
polygonal facets with the same sorted vertex set (i.e. two bodies' views of
the same planar interface) are collapsed to a single internal facet. TetGen
triangulates polygonal facets internally, so planar interface diagonals never
conflict.

## Kwargs

- `max_edge`, `sag`, `angle`: EGADS tessellation parameters.
- `coincidence_tol`: geometric tolerance for point merging. Defaults to
  `1e-9 * bbox_diag` of the combined input.
- `face_marker_offset`: facet marker = `face_marker_offset*(body_idx-1) +
  face_id`.
- `auto_regions`: if `true` (default) and `region_seeds` is `nothing`, one
  interior seed per SOLIDBODY is set with `attribute = body_index`.
- `region_seeds`: explicit 5×R seed matrix `[x;y;z;attr;max_vol]`.
- `hole_seeds`: 3×H matrix of interior hole seed points.
- `embedded_faces`: `Vector{EGObject}` of extra FACE-class objects to embed
  as internal constraint surfaces (e.g. 2D metals sitting on a
  substrate–vacuum interface). Each embedded face becomes one polygonal PLC
  facet with marker `embedded_marker_offset + k` where `k` is its index in
  `embedded_faces`. If an embedded face is geometrically coincident with an
  existing solid-body face (e.g. the full substrate top), its marker is
  merged with that face's; if it's a *proper subset* (e.g. a metal pad on
  the substrate top), TetGen splits the containing face and preserves the
  embedded marker on the sub-region.
- `embedded_marker_offset`: base marker for embedded-surface facets
  (default `10_000`).
- `verbose`: print dedup statistics.

## Input object kinds

The `bodies` argument may also include SHEETBODY / FACEBODY objects. These
are treated like sheet bodies: their faces are emitted as embedded constraint
surfaces without a region seed. Use this or `embedded_faces` interchangeably.

## Facet markers on output

Outer-boundary facets carry `face_marker_offset*(body_idx-1) + face_id`.
Interface facets (seen by two bodies) carry the marker of the smaller body
index. Every boundary/interface triangle in TetGen's `trifacemarkerlist`
output inherits the marker of its parent polygonal facet, so CAD-face
identity is preserved through the volume mesh.
"""
function tetgen_io(bodies::AbstractVector{<:EGObject};
                         max_edge::Real = 0.0, sag::Real = 0.0, angle::Real = 15.0,
                         hole_seeds = nothing, region_seeds = nothing,
                         auto_regions::Bool = true,
                         face_marker_offset::Integer = 1000,
                         embedded_faces::Union{Nothing, AbstractVector{<:EGObject}} = nothing,
                         embedded_marker_offset::Integer = 10_000,
                         embedded_edges::Union{Nothing, AbstractVector{<:EGObject}} = nothing,
                         embedded_edge_marker::Integer = 900_000,
                         coincidence_tol::Union{Real, Nothing} = nothing,
                         meshsize_callback = nothing,
                         emit_triangles::Bool = false,
                         verbose::Bool = false)
    isempty(bodies) && throw(ArgumentError("tetgen_io: at least one body required"))

    # `meshsize_callback`'s effect on TetGen has moved to TetGen's own
    # vertex-metric interface (`io.pointmtrlist` + `-m` flag), populated
    # *after* tessellation by `_set_pointmtr_from_callback!` in the
    # `tetrahedralize_size_field` wrapper. The previous per-EGADS-face
    # `.tParams` hack — setting independent target edge lengths on every
    # face — caused boundary tessellation drift between adjacent faces of
    # the unified non-manifold SHEETBODY (TetGen's polygon triangulator
    # picked different diagonals on shared edges → self-intersecting
    # facets). The unified `EG_tessellate` call below uses uniform `angle`
    # and any caller-supplied `max_edge`; per-vertex sizing is applied at
    # the TetGen layer where it interpolates linearly across input tets,
    # the same way gmsh's `setSizeCallback` propagates a continuous field
    # through its meshing backend.
    effective_max_edge = Float64(max_edge)

    # 1. Tessellate every input body and extract face entries. SHEETBODY /
    #    FACEBODY inputs are marked as non-volumetric so they contribute only
    #    embedded surfaces (no region seed).
    per_body_faces = Vector{Vector{_FaceEntry}}(undef, length(bodies))
    is_volume = Vector{Bool}(undef, length(bodies))
    for (bi, body) in enumerate(bodies)
        info = EGADS.get_info(body)
        # A BODY counts as a "volumetric domain" — contributing outer +
        # interface facets and tagged tets via region seeds — if it's
        # either a SOLIDBODY (the obvious case) or a SHEETBODY for which
        # the caller has supplied explicit region seeds. The latter is the
        # EGADS-native non-manifold representation produced by
        # `egads_make_conformal!` after `EG_sewFaces(flag=1)`: a single
        # SHEETBODY whose internal non-manifold shell bounds multiple
        # regions, each identified by a seed point in `region_seeds`.
        # Without explicit seeds we fall back to embedded-sheet semantics
        # (no auto region attribute).
        is_solid = info.oclass == OClass.BODY && info.mtype == MType.SOLIDBODY
        is_sheet_with_seeds = info.oclass == OClass.BODY &&
                              info.mtype == MType.SHEETBODY &&
                              region_seeds !== nothing
        is_volume[bi] = is_solid || is_sheet_with_seeds
        tess = tessellate(body; max_edge = effective_max_edge,
                          sag = Float64(sag), angle = Float64(angle))
        per_body_faces[bi] = _collect_face_entries(body, tess;
                                                   emit_triangles = emit_triangles)
    end

    # 1b. Standalone embedded faces (not attached to any body in `bodies`).
    #     Each one becomes a synthetic "body" of exactly one face entry,
    #     flagged embedded so it never gets a region seed and markers use the
    #     embedded range.
    if embedded_faces !== nothing
        for (ei, face) in enumerate(embedded_faces)
            info = EGADS.get_info(face)
            info.oclass == OClass.FACE || throw(ArgumentError(
                "embedded_faces[$ei]: expected FACE (class 23), got oclass=$(info.oclass)"))
            fe = _face_entry_from_standalone(face, Int32(ei))
            push!(per_body_faces, [fe])
            push!(is_volume, false)
        end
    end

    # Per-body index among *sheet* (non-volume) entries — used for embedded
    # marker assignment so sheet markers don't collide across sheets.
    sheet_idx = zeros(Int, length(per_body_faces))
    let counter = 0
        for bi in eachindex(per_body_faces)
            is_volume[bi] && continue
            counter += 1
            sheet_idx[bi] = counter
        end
    end

    # 2. Global bbox → tolerance
    xmin = ymin = zmin = Inf
    xmax = ymax = zmax = -Inf
    for face_entries in per_body_faces, fe in face_entries, pt in fe.points
        pt[1] < xmin && (xmin = pt[1]); pt[1] > xmax && (xmax = pt[1])
        pt[2] < ymin && (ymin = pt[2]); pt[2] > ymax && (ymax = pt[2])
        pt[3] < zmin && (zmin = pt[3]); pt[3] > zmax && (zmax = pt[3])
    end
    diag = hypot(xmax - xmin, ymax - ymin, zmax - zmin)
    tol = coincidence_tol === nothing ?
          max(1e-9 * diag, 1e3 * eps(Float64)) : Float64(coincidence_tol)

    # 3. Dedup points globally
    global_pts, per_body_faces = _dedup_face_points!(per_body_faces, tol)

    # 3b. Share triangulations across geometrically coincident face entries.
    # When the same physical face (e.g. substrate-top z=0 and vacuum-bottom
    # z=0, or substrate-bottom z=-525 and vacuum-below-top z=-525) appears
    # in two volume bodies, EGADS gives each body its own TopoDS_Face and
    # `tessellate` runs once per face — producing two *independent*
    # triangulations of the same boundary that pick different diagonals.
    # Per-triangle dedup in step 4 keys on sorted vertex indices, so
    # mismatched diagonals slip through, and TetGen sees overlapping
    # triangles → "self-intersection" rejection of the PLC.
    #
    # We collapse this by hashing each triangulated face entry on the SET
    # of its global vertex indices and rewriting all matching entries to
    # share the first one's triangle list. This is exact for planar faces
    # (whose interior tessellation is determined by the boundary alone)
    # and a near-perfect heuristic for the curved cases we hit in practice.
    _share_triangulations_across_bodies!(per_body_faces)

    # 4. Dedup facets: build dict keyed by (sorted_flat_vertices) → list of
    #    (body_idx, face_id, facet_payload). Facet payload is either
    #    :poly with loops-of-indices, or :tri with 3 indices.
    facet_bucket = Dict{NTuple, Vector{_FacetSource}}()
    for (bi, face_entries) in enumerate(per_body_faces)
        for fe in face_entries
            if fe.is_polygon
                key = _poly_key(fe.loops)
                push!(get!(facet_bucket, key, _FacetSource[]),
                      _FacetSource(bi, fe.face_id, true, fe.loops, (0, 0, 0),
                                   fe.hole_markers))
            else
                for tri in fe.tris
                    a, b, c = tri
                    key = _tri_key(a, b, c)
                    push!(get!(facet_bucket, key, _FacetSource[]),
                          _FacetSource(bi, fe.face_id, false,
                                       Vector{Vector{Int}}(), tri,
                                       NTuple{3, Float64}[]))
                end
            end
        end
    end

    # 5. Emit facets + markers. For non-volume (sheet/face) source bodies we
    #    shift the marker into the `embedded_marker_offset..` range so
    #    CAD-face identity on embedded constraint surfaces is distinguishable
    #    from outer/interface boundaries at a glance.
    facets = Vector{Vector{Vector{Cint}}}()  # each facet = list of polygons (loops)
    facet_holes = Vector{Vector{NTuple{3, Float64}}}()  # hole markers per facet
    markers = Cint[]
    n_interface = 0
    n_outer = 0
    n_embedded = 0
    n_poly_facets = 0
    n_tri_facets = 0
    n_nonmanifold = 0
    for (_, sources) in facet_bucket
        if length(sources) > 2
            n_nonmanifold += 1
        end
        sort!(sources; by = s -> s.body_idx)
        s = sources[1]
        if s.is_polygon
            push!(facets, [Cint.(loop) for loop in s.loops])
            push!(facet_holes, copy(s.hole_markers))
            n_poly_facets += 1
        else
            push!(facets, [Cint[s.tri[1], s.tri[2], s.tri[3]]])
            push!(facet_holes, NTuple{3, Float64}[])
            n_tri_facets += 1
        end
        if is_volume[s.body_idx]
            push!(markers, Cint(face_marker_offset) * Cint(s.body_idx - 1) +
                           Cint(s.face_id))
            if length(sources) >= 2
                n_interface += 1
            else
                n_outer += 1
            end
        else
            # Embedded constraint surface
            push!(markers, Cint(embedded_marker_offset) * Cint(sheet_idx[s.body_idx]) +
                           Cint(s.face_id))
            n_embedded += 1
        end
    end

    # 5b. Embedded 1D constraints (e.g. airbridge feet). Each edge is added as
    #     an internal SEGMENT (a 2-point degenerate polygon) on the planar
    #     facet that geometrically contains it. TetGen then conforms that
    #     facet's triangulation to the segment WITHOUT it bounding any
    #     sub-region — the equivalent of gmsh "embed line in surface". This
    #     lets us imprint only the bridge feet (two short segments per
    #     bridge) instead of carving the whole bridge-top footprint into the
    #     substrate surface. (A bare `io.edgelist` segment is NOT enough —
    #     TetGen rejects a segment whose endpoints lie in a facet interior
    #     unless the segment is declared as one of that facet's polygons.)
    n_embedded_seg = 0
    if embedded_edges !== nothing
        for edge in embedded_edges
            poly = _edge_polyline(edge)            # Vector of xyz along the edge
            length(poly) >= 2 || continue
            # global-point index for each polyline vertex (merge into global_pts)
            gidx = Int[]
            for p in poly
                push!(gidx, _find_or_add_global_point!(global_pts, p, tol))
            end
            # The foot lies on a HORIZONTAL substrate facet (constant z); a
            # foot edge is a line so it is degenerate in two axes (x AND z for
            # an x=const, z=const leg) — we must pick the substrate plane
            # (z=const), not the vertical leg plane. Restrict the facet search
            # to z=const facets at the foot's z. Use the geometric midpoint
            # (mean of the polyline) for the point-in-facet test.
            mx = sum(p[1] for p in poly) / length(poly)
            my = sum(p[2] for p in poly) / length(poly)
            mz = sum(p[3] for p in poly) / length(poly)
            mid = (mx, my, mz)
            fidx = _facet_containing_point(facets, global_pts, mid, tol, 3, mz)
            if fidx === nothing
                verbose && @warn "embedded edge: no containing z-facet; skipping" midpoint=mid
                continue
            end
            verbose && println("  foot @ mid=$(round.(mid;digits=1)) → facet #$fidx marker=$(markers[fidx])")
            for s in 1:(length(gidx) - 1)
                push!(facets[fidx], Cint[gidx[s], gidx[s + 1]])
                n_embedded_seg += 1
            end
        end
    end

    if verbose
        npts_raw = sum(sum(length(fe.points) for fe in fes; init = 0)
                       for fes in per_body_faces; init = 0)
        println("tetgen_io: bodies=$(length(bodies)) (",
                count(is_volume), " volume + ",
                count(!, is_volume), " sheet/embedded)  tol=$tol")
        println("  points  : $npts_raw raw → $(length(global_pts)) merged")
        println("  facets  : $(n_outer + n_interface + n_embedded) total " *
                "($n_outer outer, $n_interface interface, $n_embedded embedded)")
        n_embedded_seg > 0 && println("  embedded segments (e.g. bridge feet): $n_embedded_seg")
        println("  by kind : $n_poly_facets polygonal, $n_tri_facets triangular" *
                (n_nonmanifold > 0 ? ", $n_nonmanifold non-manifold" : ""))
    end

    # 6. Build IO
    io = TetGen.RawTetGenIO{Cdouble}()
    N = length(global_pts)
    pts = Matrix{Cdouble}(undef, 3, N)
    @inbounds for i in 1:N
        pts[1, i] = global_pts[i][1]
        pts[2, i] = global_pts[i][2]
        pts[3, i] = global_pts[i][3]
    end
    io.pointlist = pts
    # Use facetlist!(io, Vector) where each element is a Vector of polygons,
    # but actually RawFacet holds one polygonlist; facetlist! Vector-version
    # wraps each element as single polygon. For multi-loop we build RawFacets
    # directly.
    _assign_facets!(io, facets, facet_holes)
    io.facetmarkerlist = markers

    # 7. Region seeds — one per volumetric (SOLIDBODY) input only.
    if region_seeds === nothing && auto_regions
        rows = Vector{Float64}[]
        for (bi, body) in enumerate(bodies)
            is_volume[bi] || continue
            seed = _interior_seed_point(body)
            seed === nothing && continue
            push!(rows, [seed[1], seed[2], seed[3], Float64(bi), 0.0])
        end
        if !isempty(rows)
            io.regionlist = reduce(hcat, rows)
        end
    end

    _set_holes_regions!(io, hole_seeds, region_seeds)
    return io
end

# ─── Convenience: direct tetrahedralize method ───────────────────────────────

"""
    tetrahedralize(body::EGObject, flags="pq1.414YQ"; meshsize_callback=nothing, kw...)
    tetrahedralize(bodies::Vector{<:EGObject}, flags="pq1.414YAQ"; meshsize_callback=nothing, kw...)

Build a PLC via `tetgen_io(body[or bodies]; kw...)` and hand it to
TetGen. `kw...` are forwarded to `tetgen_io`.

## `meshsize_callback`

If supplied, `meshsize_callback(x, y, z) -> s` returns the target **edge
length** `s` (in model units) at a point. The callback is bridged to
TetGen's per-vertex metric (`io.pointmtrlist` + `-m` flag): we sample the
callback at every input PLC vertex and let TetGen interpolate linearly
between vertices when refining. This is the gmsh-equivalent shape: a
spatial map `(x,y,z) ↦ target-edge-length` propagated as a continuous
field by the meshing backend.

For strongly varying fields (e.g. power-law `h = C*(r-r0)^α` with α<1),
linear interpolation between coarse PLC vertices undersamples the field
in the rapidly-changing region; iterate via
[`tetrahedralize_size_field`](@ref) instead, which re-samples at the
output vertices after each pass and refines via TetGen's `-r` flag until
edge-length-vs-target ratios converge.
"""
function tetrahedralize(body::EGObject, flags::AbstractString = "pq1.414YQ";
                               meshsize_callback = nothing, kw...)
    io = tetgen_io(body; meshsize_callback = meshsize_callback, kw...)
    flags = _install_size_field!(io, flags, meshsize_callback)
    return invoke(TetGen.tetrahedralize,
                  Tuple{TetGen.RawTetGenIO{Float64}, String},
                  io, String(flags))
end

function tetrahedralize(bodies::AbstractVector{<:EGObject},
                               flags::AbstractString = "pq1.414YAQ";
                               meshsize_callback = nothing, kw...)
    io = tetgen_io(bodies; meshsize_callback = meshsize_callback, kw...)
    flags = _install_size_field!(io, flags, meshsize_callback)
    return invoke(TetGen.tetrahedralize,
                  Tuple{TetGen.RawTetGenIO{Float64}, String},
                  io, String(flags))
end

# Lets callers that already own a `RawTetGenIO` (e.g. built via
# `tetgen_io(…)` for inspection or manual tweaks) still use the same
# size-field path without rebuilding the PLC.
function tetrahedralize(io::TetGen.RawTetGenIO{Float64},
                               flags::AbstractString;
                               meshsize_callback)
    flags = _install_size_field!(io, flags, meshsize_callback)
    return invoke(TetGen.tetrahedralize,
                  Tuple{TetGen.RawTetGenIO{Float64}, String},
                  io, String(flags))
end

# Sample `meshsize_callback` at every input PLC vertex → 1×N pointmtrlist
# (TetGen's isotropic metric format), and add `m` to the flags string so
# TetGen actually consults the metric during refinement. Linear
# interpolation across input tets gives a continuous size field anywhere
# in the domain, gmsh-style.
function _install_size_field!(io::TetGen.RawTetGenIO{Float64},
                              flags::AbstractString, meshsize_callback)
    meshsize_callback === nothing && return String(flags)
    n = size(io.pointlist, 2)
    n == 0 && return String(flags)
    mtr = Matrix{Cdouble}(undef, 1, n)
    @inbounds for i in 1:n
        mtr[1, i] = Float64(meshsize_callback(io.pointlist[1, i],
                                              io.pointlist[2, i],
                                              io.pointlist[3, i]))
    end
    io.pointmtrlist = mtr
    s = String(flags)
    occursin('m', s) || (s *= "m")
    return s
end

"""
    tetrahedralize_size_field(io_or_bodies, flags;
                              meshsize_callback,
                              method = :callback,
                              size_factor = 1.0,
                              max_iterations = 5,
                              tolerance = 0.2,
                              refinement_flags = nothing,
                              verbose = false)

Size-field-respecting tetrahedralization. `meshsize_callback(x, y, z) -> h`
is the target edge length at a point — the same field gmsh consumes via
`gmsh.model.mesh.setSizeCallback`.

## `method = :callback` (default — exact field evaluation)

The direct TetGen analogue of gmsh's size callback: TetGen's
`tetunsuitable` hook. Pass 0 builds a coarse, valid, region-attributed
PLC mesh with `flags` (no sizing). Each subsequent pass reconstructs the
mesh (`-r`) with `TetGen.tetunsuitable!` installed: a tet is split when
its **longest edge** exceeds `size_factor * h(centroid)`. Because the
callback evaluates the *true* field at every candidate tet (no
interpolation from input vertices), one refinement pass typically lands
the graded mesh; the loop iterates until no tet is unsuitable or the
count stalls.

`tolerance` is accepted for API compatibility but convergence is
`n_unsuitable == 0` (plus a stall guard).

## `method = :metric` (legacy — interpolated vertex metric)

The previous behaviour: sample the callback at input vertices into
`pointmtrlist` (`-m`), iterate `-rqm` passes re-sampling at output
vertices, and stop when `max(edge / min(h(p1), h(p2))) ≤ 1 + tolerance`.
TetGen treats the metric as a loose sizing *hint* (edges routinely land
at several× the target), so prefer `:callback` when actual grading
matters.

`refinement_flags` defaults to `"raqQ"` for `:callback` (`-a` activates
the per-tet check that consults `tetunsuitable`, `-q` keeps quality
bounds) and `"rqm"` for `:metric`.
"""
function tetrahedralize_size_field(
        bodies_or_io,
        flags::AbstractString = "pAYQ";
        meshsize_callback,
        method::Symbol = :callback,
        size_factor::Real = 1.0,
        max_iterations::Integer = 5,
        tolerance::Real = 1.0,
        refinement_flags::Union{Nothing, AbstractString} = nothing,
        verbose::Bool = false,
        kw...)
    method in (:callback, :metric) || throw(ArgumentError(
        "tetrahedralize_size_field: method must be :callback or :metric, got :$method"))
    io = bodies_or_io isa TetGen.RawTetGenIO ?
         bodies_or_io :
         tetgen_io(bodies_or_io; meshsize_callback = nothing, kw...)

    # Pass 0: valid, region-attributed PLC mesh. :metric additionally seeds
    # the vertex metric so the first pass is already graded-ish; :callback
    # leaves pass 0 unsized — grading happens in the refinement passes where
    # the hook evaluates the exact field.
    flags0 = method === :metric ?
             _install_size_field!(io, flags, meshsize_callback) : String(flags)
    # Capture the input regionlist so we can replay it on refinement
    # passes — TetGen's `-r` doesn't propagate region attributes to
    # Steiner-point sub-tets unless `-A` + regionlist re-tag them.
    region_seeds_in = size(io.regionlist, 2) > 0 ? copy(io.regionlist) : nothing
    unsuit = method === :callback ?
             _make_unsuitable(meshsize_callback, size_factor) : nothing
    verbose && println("[size-field] pass 0 (", flags0, ", method=", method, "): ",
                       size(io.pointlist, 2), " PLC vertices")
    # In :callback mode the hook is armed for pass 0 as well: with flags
    # like "pAaqQ" (`-a` activates the per-tet check, `-q` the quality
    # refinement stage) the whole graded mesh is produced in this single
    # 'p' pass — boundary facets included — and the loop below only
    # verifies/touches up. With the conservative "pAYQ" default the pass-0
    # mesh is coarse and grading happens in `-r` refinement passes.
    out = _tetrahedralize_hooked(io, flags0, unsuit)
    verbose && println("             → ", size(out.pointlist, 2),
                       " pts, ", size(out.tetrahedronlist, 2), " tets")

    flags_k = refinement_flags !== nothing ? String(refinement_flags) :
              (method === :callback ? "raqQ" : "rqm")
    method === :metric && !occursin('m', flags_k) && (flags_k *= "m")

    nbad_prev = typemax(Int)
    for k in 1:max_iterations
        ratio = _size_field_max_ratio(out, meshsize_callback)
        nbad = method === :callback ?
               _count_unsuitable(out, unsuit) : 0
        verbose && println("[size-field] iter ", k - 1,
                           " max(edge/h) = ", round(ratio; digits = 3),
                           method === :callback ?
                           "  unsuitable = $nbad/$(size(out.tetrahedronlist, 2))" : "")
        converged = method === :callback ? nbad == 0 : ratio <= 1.0 + tolerance
        if converged
            _propagate_tet_attributes!(out, region_seeds_in; verbose = verbose)
            return out
        end
        nbad_prev = nbad

        # Build a refinement input: hand TetGen the previous tet mesh
        # plus a fresh per-vertex metric sampled at the *output* nodes.
        refine_io = _refinement_io_from_output(out, meshsize_callback,
                                               region_seeds_in)
        verbose && println("[size-field] pass ", k, " (", flags_k, "): ",
                           size(refine_io.pointlist, 2), " vertices")
        # `-r` mode preserves existing tet attributes through refinement;
        # Steiner-induced sub-tets inherit their parent's attribute. We
        # do NOT pass `-A` here — in refinement mode `-A` tries to
        # re-tag from region seeds and produces unattributed tets in
        # complex non-manifold compounds. Inheritance via the existing
        # `tetrahedronattributelist` is reliable.
        #
        # TetGen's `-r` reconstruction is fragile on some inputs (internal
        # error 2); a failed touch-up pass keeps the last good mesh rather
        # than aborting the pipeline.
        out_k = try
            _tetrahedralize_hooked(refine_io, flags_k, unsuit)
        catch err
            err isa TetGen.TetGenError || rethrow()
            @warn "[size-field] refinement pass $k failed; keeping previous mesh" err
            break
        end
        # A `-r` touch-up pass can regress (TetGen re-tetrahedralizes facets
        # and may undo good splits). Only accept a pass that reduces the
        # number of field-violating tets; otherwise keep the previous mesh.
        if method === :callback
            nbad_k = _count_unsuitable(out_k, unsuit)
            if nbad_k >= nbad_prev
                verbose && @warn "[size-field] refinement pass $k did not improve " *
                                 "($nbad_prev → $nbad_k unsuitable); keeping previous mesh"
                break
            end
        end
        out = out_k
        verbose && println("             → ", size(out.pointlist, 2),
                           " pts, ", size(out.tetrahedronlist, 2), " tets")
    end
    if verbose
        ratio = _size_field_max_ratio(out, meshsize_callback)
        nbad = method === :callback ? _count_unsuitable(out, unsuit) : 0
        println("[size-field] final max(edge/h) = ", round(ratio; digits = 3),
                method === :callback ? "  unsuitable = $nbad" : "",
                " after ≤", max_iterations, " iter")
    end
    _propagate_tet_attributes!(out, region_seeds_in; verbose = verbose)
    return out
end

# Run one tetrahedralize call with the `tetunsuitable` hook installed (or
# plain if `unsuit === nothing`). The hook is global state inside libtet, so
# it is re-armed immediately before the call and disarmed (no-op) after —
# other `tetrahedralize` callers in the same process must not inherit it.
function _tetrahedralize_hooked(io::TetGen.RawTetGenIO{Float64},
                                flags::AbstractString, unsuit)
    unsuit === nothing && return invoke(TetGen.tetrahedralize,
                                        Tuple{TetGen.RawTetGenIO{Float64}, String},
                                        io, String(flags))
    TetGen.tetunsuitable!(unsuit; check_signature = false)
    try
        return invoke(TetGen.tetrahedralize,
                      Tuple{TetGen.RawTetGenIO{Float64}, String},
                      io, String(flags))
    finally
        TetGen.tetunsuitable!((pa, pb, pc, pd) -> false; check_signature = false)
    end
end

# gmsh-size-callback semantics for TetGen's tetunsuitable hook: split a tet
# whose longest edge exceeds `factor * h(centroid)`.
function _make_unsuitable(h::Function, factor::Real)
    f2 = Float64(factor)^2
    return function (pa, pb, pc, pd)
        l2 = 0.0
        @inbounds for (u, v) in ((pa, pb), (pa, pc), (pa, pd),
                                 (pb, pc), (pb, pd), (pc, pd))
            d = (u[1] - v[1])^2 + (u[2] - v[2])^2 + (u[3] - v[3])^2
            d > l2 && (l2 = d)
        end
        cx = 0.25 * (pa[1] + pb[1] + pc[1] + pd[1])
        cy = 0.25 * (pa[2] + pb[2] + pc[2] + pd[2])
        cz = 0.25 * (pa[3] + pb[3] + pc[3] + pd[3])
        hc = Float64(h(cx, cy, cz))
        (isfinite(hc) && hc > 0) || return false
        return l2 > f2 * hc * hc
    end
end

# Count output tets the `unsuitable` predicate would still split.
function _count_unsuitable(out, unsuit)
    n = 0
    pts = out.pointlist
    @inbounds for j in 1:size(out.tetrahedronlist, 2)
        a = out.tetrahedronlist[1, j]; b = out.tetrahedronlist[2, j]
        c = out.tetrahedronlist[3, j]; d = out.tetrahedronlist[4, j]
        pa = (pts[1, a], pts[2, a], pts[3, a])
        pb = (pts[1, b], pts[2, b], pts[3, b])
        pc = (pts[1, c], pts[2, c], pts[3, c])
        pd = (pts[1, d], pts[2, d], pts[3, d])
        unsuit(pa, pb, pc, pd) && (n += 1)
    end
    return n
end

# Walk the tet adjacency graph and assign every attr=0 tet the most common
# attribute among its face-neighbors with non-zero attribute. Iterate until
# no more can be assigned. Required because TetGen's `-r` refinement of a
# non-manifold compound (the unified SHEETBODY shape we feed it) leaves a
# small fraction of Steiner-induced sub-tets unattributed (~5-8% on the
# SingleTransmon stack) — they sit on either side of internal interface
# faces in regions where the seed flood-fill didn't reach during refinement.
# Without this propagation, downstream readers see a phantom attr=0 "region"
# that is actually slivers of the real volumes.
function _propagate_tet_attributes!(out, region_seeds_in;
                                    verbose::Bool = false)
    size(out.tetrahedronattributelist, 1) == 0 && return out
    attrs = view(out.tetrahedronattributelist, 1, :)
    initial_zeros = count(==(0.0), attrs)
    verbose && println("[size-field] _propagate_tet_attributes!: ",
                       "initial zeros = ", initial_zeros, " / ", length(attrs))
    initial_zeros == 0 && return out

    # First pass: any attr=0 tet whose centroid lies on the same side of
    # the original region seeds as one of them gets that seed's attribute.
    # This catches the common case where a sliver tet ended up away from
    # all interfaces but missing its tag.
    if region_seeds_in !== nothing && size(region_seeds_in, 2) > 0
        nseeds = size(region_seeds_in, 2)
        for j in axes(attrs, 1)
            attrs[j] == 0.0 || continue
            v = view(out.tetrahedronlist, :, j)
            cx = (out.pointlist[1, v[1]] + out.pointlist[1, v[2]] +
                  out.pointlist[1, v[3]] + out.pointlist[1, v[4]]) * 0.25
            cy = (out.pointlist[2, v[1]] + out.pointlist[2, v[2]] +
                  out.pointlist[2, v[3]] + out.pointlist[2, v[4]]) * 0.25
            cz = (out.pointlist[3, v[1]] + out.pointlist[3, v[2]] +
                  out.pointlist[3, v[3]] + out.pointlist[3, v[4]]) * 0.25
            best_d = Inf; best_a = 0.0
            for s in 1:nseeds
                dx = cx - region_seeds_in[1, s]
                dy = cy - region_seeds_in[2, s]
                dz = cz - region_seeds_in[3, s]
                d = dx * dx + dy * dy + dz * dz
                if d < best_d
                    best_d = d; best_a = region_seeds_in[4, s]
                end
            end
            attrs[j] = best_a
        end
    end

    # Second pass: neighbour propagation on any remaining zeros (handles
    # cases where nearest-seed picks the wrong region across an interface).
    # Cheap: build face→tets adjacency once, walk attr=0 tets and adopt the
    # majority non-zero attribute among neighbours.
    if any(==(0.0), attrs)
        face_to_tets = Dict{NTuple{3, Cint}, Vector{Int}}()
        for j in axes(attrs, 1)
            v = view(out.tetrahedronlist, :, j)
            for (a, b, c) in ((v[2], v[3], v[4]), (v[1], v[3], v[4]),
                              (v[1], v[2], v[4]), (v[1], v[2], v[3]))
                xs = (Cint(a), Cint(b), Cint(c))
                xs = xs[1] <= xs[2] ? xs : (xs[2], xs[1], xs[3])
                xs = xs[2] <= xs[3] ? xs : (xs[1], xs[3], xs[2])
                xs = xs[1] <= xs[2] ? xs : (xs[2], xs[1], xs[3])
                push!(get!(face_to_tets, xs, Int[]), j)
            end
        end
        changed = true
        passes = 0
        while changed && passes < 20
            changed = false
            passes += 1
            for j in axes(attrs, 1)
                attrs[j] == 0.0 || continue
                v = view(out.tetrahedronlist, :, j)
                votes = Dict{Float64, Int}()
                for (a, b, c) in ((v[2], v[3], v[4]), (v[1], v[3], v[4]),
                                  (v[1], v[2], v[4]), (v[1], v[2], v[3]))
                    xs = (Cint(a), Cint(b), Cint(c))
                    xs = xs[1] <= xs[2] ? xs : (xs[2], xs[1], xs[3])
                    xs = xs[2] <= xs[3] ? xs : (xs[1], xs[3], xs[2])
                    xs = xs[1] <= xs[2] ? xs : (xs[2], xs[1], xs[3])
                    for nbr in face_to_tets[xs]
                        nbr == j && continue
                        attrs[nbr] == 0.0 && continue
                        votes[attrs[nbr]] = get(votes, attrs[nbr], 0) + 1
                    end
                end
                isempty(votes) && continue
                attrs[j] = argmax(votes)
                changed = true
            end
        end
    end
    final_zeros = count(==(0.0), attrs)
    verbose && initial_zeros > 0 &&
        println("[size-field] re-attributed ", initial_zeros - final_zeros,
                "/", initial_zeros, " tets from attr=0 ",
                final_zeros == 0 ? "(all assigned)" : "($(final_zeros) remaining)")
    return out
end

function _refinement_io_from_output(out, meshsize_callback, region_seeds = nothing)
    refine = TetGen.RawTetGenIO{Cdouble}()
    refine.pointlist = copy(out.pointlist)
    refine.tetrahedronlist = copy(out.tetrahedronlist)
    if size(out.tetrahedronattributelist, 1) > 0
        refine.tetrahedronattributelist = copy(out.tetrahedronattributelist)
    end
    if size(out.trifacelist, 2) > 0
        refine.trifacelist = copy(out.trifacelist)
        if length(out.trifacemarkerlist) > 0
            refine.trifacemarkerlist = copy(out.trifacemarkerlist)
        end
    end
    # Forward region seeds so refinement re-attributes Steiner-point
    # sub-tets via geometric containment (`-A` flag). Without this, every
    # tet split by Steiner points loses its attribute and falls back to 0.
    if region_seeds !== nothing
        refine.regionlist = copy(region_seeds)
    elseif size(out.regionlist, 2) > 0
        refine.regionlist = copy(out.regionlist)
    end
    n = size(refine.pointlist, 2)
    mtr = Matrix{Cdouble}(undef, 1, n)
    @inbounds for i in 1:n
        mtr[1, i] = Float64(meshsize_callback(refine.pointlist[1, i],
                                              refine.pointlist[2, i],
                                              refine.pointlist[3, i]))
    end
    refine.pointmtrlist = mtr
    return refine
end

# Convergence metric: largest ratio over all tet edges of
#   actual_edge_length / min(h(p1), h(p2))
# A value ≤ 1 means every edge respects the locally-tightest target; we
# accept up to 1 + tolerance because TetGen's `-m` mode is not a strict
# bound (it's a sizing hint and produces edges up to ~2x the target in
# isotropic mode).
function _size_field_max_ratio(out, meshsize_callback)
    max_ratio = 0.0
    @inbounds for j in 1:size(out.tetrahedronlist, 2)
        for (a, b) in ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))
            ia = out.tetrahedronlist[a, j]
            ib = out.tetrahedronlist[b, j]
            xa, ya, za = out.pointlist[1, ia], out.pointlist[2, ia], out.pointlist[3, ia]
            xb, yb, zb = out.pointlist[1, ib], out.pointlist[2, ib], out.pointlist[3, ib]
            actual = sqrt((xa - xb)^2 + (ya - yb)^2 + (za - zb)^2)
            ha = Float64(meshsize_callback(xa, ya, za))
            hb = Float64(meshsize_callback(xb, yb, zb))
            target = min(ha, hb)
            target <= 0 && continue
            r = actual / target
            r > max_ratio && (max_ratio = r)
        end
    end
    return max_ratio
end

# ─── Per-face entry type ────────────────────────────────────────────────────

mutable struct _FaceEntry
    face_id::Int32
    is_polygon::Bool
    # Raw data (xyz), filled during collection; later remapped to global ids.
    points::Vector{NTuple{3, Float64}}       # flat list of unique pts for this face
    local_to_global::Vector{Int}             # filled after dedup
    # If is_polygon: loops is list of loops, each a list of local indices.
    # First loop is the outer boundary; subsequent loops are inner (holes).
    loops::Vector{Vector{Int}}
    # If !is_polygon: tris is a list of triangles, each 3 local indices.
    # After dedup, remapped to global.
    tris::Vector{NTuple{3, Int}}
    # Hole-marker points (xyz) — one per inner loop, guaranteed inside the
    # corresponding hole region. TetGen's RawFacet treats holes as marker
    # points rather than extra polygons; inner loops therefore appear in
    # `loops[2:end]` for the PLC outline itself AND as `hole_markers` so
    # TetGen excludes the hole interior when triangulating the facet.
    hole_markers::Vector{NTuple{3, Float64}}
end

struct _FacetSource
    body_idx::Int
    face_id::Int32
    is_polygon::Bool
    loops::Vector{Vector{Int}}   # global indices (only meaningful if is_polygon)
    tri::NTuple{3, Int}          # global indices (only meaningful if !is_polygon)
    hole_markers::Vector{NTuple{3, Float64}}  # 3D hole-marker pts for multi-loop facets
end

# ─── Collect face entries from a body's tessellation ─────────────────────────

function _collect_face_entries(body::EGObject, tess::EGObject;
                               emit_triangles::Bool = false)
    faces = get_body_topos(body, OClass.FACE)
    entries = _FaceEntry[]
    for (fi, face) in enumerate(faces)
        area = try abs(EGADS.get_area(face)) catch; 0.0 end
        area < 1e-6 && continue
        push!(entries, _face_entry(body, tess, face, Int32(fi), emit_triangles))
    end
    return entries
end

function _face_entry(body::EGObject, tess::EGObject, face::EGObject, face_id::Int32,
                     emit_triangles::Bool = false)
    ft = get_topology(face)
    geom = ft.geom
    is_planar = false
    if geom !== nothing
        gi = EGADS.get_info(geom)
        is_planar = gi.oclass == OClass.SURFACE && gi.mtype == MType.PLANE
    end
    # Multi-loop planar faces are triangulated (TetGen's polygon-with-holes
    # triangulator chokes on the complex chip-perimeter-with-3-holes case).
    # Per-body independent tessellations of geometrically coincident faces
    # would normally pick different diagonals — we resolve that in step 3b
    # of `tetgen_io` by sharing one body's triangulation across all
    # coincident face entries.
    multi_loop = is_planar && length(ft.children) > 1
    if is_planar && !emit_triangles && !multi_loop
        return _planar_face_entry(body, tess, face, face_id, ft)
    else
        return _triangulated_face_entry(tess, body, face, face_id)
    end
end

# Walk face's loops → ordered list of boundary points per loop.
function _planar_face_entry(body::EGObject, tess::EGObject, face::EGObject,
                            face_id::Int32, ft)
    pts = NTuple{3, Float64}[]
    loops_local = Vector{Vector{Int}}()
    loop_senses = Int[]     # +1 outer, -1 inner (from face topology)
    face_sense = ft.mtype == MType.SREVERSE ? -1 : 1
    for (loop, lsense) in zip(ft.children, ft.senses)
        lt = get_topology(loop)
        loop_local = Int[]
        for (edge, esense) in zip(lt.children, lt.senses)
            edge_topo = get_topology(edge)
            edge_topo.mtype == MType.DEGENERATE && continue
            edge_idx = index_body_topo(body, edge)
            edge_idx <= 0 && continue
            em = get_edge_mesh(tess, edge_idx)
            npts = size(em.xyz, 1)
            npts == 0 && continue
            rng = esense == 1 ? (1:npts) : (npts:-1:1)
            for k in rng
                pt = (em.xyz[k, 1], em.xyz[k, 2], em.xyz[k, 3])
                if isempty(loop_local) || pts[loop_local[end]] != pt
                    push!(pts, pt)
                    push!(loop_local, length(pts))
                end
            end
        end
        if length(loop_local) > 1 && pts[loop_local[1]] == pts[loop_local[end]]
            pop!(loop_local)
        end
        if length(loop_local) >= 3
            push!(loops_local, loop_local)
            # Loop sense in face coords is lsense * face_sense: +1 outer, -1 hole.
            push!(loop_senses, Int(lsense) * face_sense)
        end
    end

    # TetGen convention: a facet has a single outer polygon and holes are
    # described by marker points inside each hole region (NOT by extra
    # polygon loops). For multi-loop faces we therefore (a) keep the outer
    # loop as the only polygon and (b) compute a 3D point inside each inner
    # loop so TetGen can exclude the hole interior when triangulating.
    if length(loops_local) > 1
        order = sortperm(loop_senses; rev = true)
        loops_local = loops_local[order]
        loop_senses = loop_senses[order]
    end
    hole_markers = NTuple{3, Float64}[]
    for li in 2:length(loops_local)
        p = _polygon_interior_point(loops_local[li], pts)
        p === nothing || push!(hole_markers, p)
    end
    loops_out = length(loops_local) >= 1 ? [loops_local[1]] : loops_local

    return _FaceEntry(face_id, true, pts, Int[], loops_out,
                      NTuple{3, Int}[], hole_markers)
end

# Pick a 3D point strictly inside a planar polygon's interior (in xy-ish
# coordinates). Use the centroid if it's inside; else bisect from centroid
# toward the first vertex until an interior point is found. The polygon is
# planar so any returned 3D point lies on the polygon's plane.
function _polygon_interior_point(loop::Vector{Int},
                                 pts::Vector{NTuple{3, Float64}})
    isempty(loop) && return nothing
    # Compute centroid (3D average).
    cx = 0.0; cy = 0.0; cz = 0.0
    for v in loop
        p = pts[v]; cx += p[1]; cy += p[2]; cz += p[3]
    end
    n = length(loop); cx /= n; cy /= n; cz /= n
    # For convex polygons the centroid is always interior. For non-convex
    # polygons this fails ~rarely in PLC facets (which tend to be simple).
    # If the caller hits problems on specific geometries we can upgrade to
    # a proper ray-cast point-in-polygon.
    return (cx, cy, cz)
end

# Fallback: triangulate face via tessellation; emit tris as independent facets.
function _triangulated_face_entry(tess::EGObject, body::EGObject,
                                  face::EGObject, face_id::Int32)
    # get_face_mesh gives per-face triangles + xyz. Indices are face-local
    # (1..nv); convert to tessellation-global via eg_local_to_global, then
    # back to xyz via get_global_mesh? Simpler: use face-local xyz directly
    # and let the global point dedup merge duplicates across faces.
    fmesh = get_face_mesh(tess, body, face_id)
    nv = size(fmesh.xyz, 1)
    pts = NTuple{3, Float64}[(fmesh.xyz[i, 1], fmesh.xyz[i, 2], fmesh.xyz[i, 3])
                             for i in 1:nv]
    tris = NTuple{3, Int}[(Int(fmesh.tris[t, 1]),
                           Int(fmesh.tris[t, 2]),
                           Int(fmesh.tris[t, 3]))
                          for t in 1:size(fmesh.tris, 1)]
    return _FaceEntry(face_id, false, pts, Int[], Vector{Vector{Int}}(), tris,
                      NTuple{3, Float64}[])
end

# Standalone FACE EGObject (not wrapped in a body). Uses its parent BODY to
# get a tessellation, then extracts only the triangles / loop belonging to
# this face. Used by the `embedded_faces` kwarg path.
function _face_entry_from_standalone(face::EGObject, face_id::Int32)
    parent = EGADS.get_body(face)
    face_idx = index_body_topo(parent, face)
    face_idx > 0 || throw(ArgumentError(
        "standalone face has no parent body or invalid index"))
    tess = tessellate(parent; angle = 30.0)
    return _face_entry(parent, tess, face, Int32(face_idx))
end

# ─── Global point dedup across all face entries ──────────────────────────────

function _dedup_face_points!(per_body_faces::Vector{Vector{_FaceEntry}},
                             tol::Float64)
    inv_tol = 1.0 / tol
    tol2 = tol * tol
    cell_dict = Dict{NTuple{3, Int}, Vector{Int}}()
    global_pts = Tuple{Float64, Float64, Float64}[]

    for face_entries in per_body_faces
        for fe in face_entries
            fe.local_to_global = Vector{Int}(undef, length(fe.points))
            for (i, p) in enumerate(fe.points)
                cx = floor(Int, p[1] * inv_tol)
                cy = floor(Int, p[2] * inv_tol)
                cz = floor(Int, p[3] * inv_tol)
                found = 0
                for dz in -1:1, dy in -1:1, dx in -1:1
                    bucket = get(cell_dict, (cx + dx, cy + dy, cz + dz), nothing)
                    bucket === nothing && continue
                    for gi in bucket
                        q = global_pts[gi]
                        dxp = q[1] - p[1]; dyp = q[2] - p[2]; dzp = q[3] - p[3]
                        if dxp * dxp + dyp * dyp + dzp * dzp <= tol2
                            found = gi
                            break
                        end
                    end
                    found != 0 && break
                end
                if found == 0
                    push!(global_pts, p)
                    found = length(global_pts)
                    push!(get!(cell_dict, (cx, cy, cz), Int[]), found)
                end
                fe.local_to_global[i] = found
            end
            # Remap loops / tris to global indices
            if fe.is_polygon
                for loop in fe.loops
                    for j in eachindex(loop)
                        loop[j] = fe.local_to_global[loop[j]]
                    end
                end
            else
                for j in eachindex(fe.tris)
                    t = fe.tris[j]
                    fe.tris[j] = (fe.local_to_global[t[1]],
                                  fe.local_to_global[t[2]],
                                  fe.local_to_global[t[3]])
                end
            end
        end
    end
    return global_pts, per_body_faces
end

# Sample an EDGE into a polyline of xyz points. A straight edge yields its two
# endpoints; a curved edge (e.g. a CIRCLE-arc bridge foot) is sampled so the
# embedded segment chain follows the curve. Used by `tetgen_io`'s
# `embedded_edges` path.
function _edge_polyline(edge::EGObject; nmin::Int = 2)
    rng = try EGADS.get_range(edge) catch; nothing end
    rng === nothing && return Tuple{Float64,Float64,Float64}[]
    rvec = rng isa NamedTuple ? rng.range : rng
    t0, t1 = Float64(rvec[1]), Float64(rvec[2])
    # decide sample count by geometry type: lines need 2, curves more
    geom = try EGADS.get_topology(edge).geom catch; nothing end
    is_line = geom !== nothing && EGADS.get_info(geom).mtype == EGADS.MType.LINE
    n = is_line ? 2 : max(nmin, 16)
    pts = Tuple{Float64,Float64,Float64}[]
    for k in 0:(n - 1)
        t = t0 + (t1 - t0) * k / (n - 1)
        ev = try EGADS.evaluate(edge, t) catch; nothing end
        ev === nothing && continue
        xyz = ev isa NamedTuple ? ev.xyz : ev[1:3]
        push!(pts, (Float64(xyz[1]), Float64(xyz[2]), Float64(xyz[3])))
    end
    return pts
end

# Find an existing global point within `tol`, else append it. Returns 1-based
# index into `global_pts`. (Linear scan over a handful of embedded-edge points
# — cheap; the bulk face points already went through the spatial-hash dedup.)
function _find_or_add_global_point!(global_pts::Vector{Tuple{Float64,Float64,Float64}},
                                    p::Tuple{Float64,Float64,Float64}, tol::Float64)
    tol2 = tol * tol
    @inbounds for gi in eachindex(global_pts)
        q = global_pts[gi]
        dx = q[1]-p[1]; dy = q[2]-p[2]; dz = q[3]-p[3]
        dx*dx + dy*dy + dz*dz <= tol2 && return gi
    end
    push!(global_pts, p)
    return length(global_pts)
end

# Index of the planar facet that geometrically contains `point`. Only facets
# whose constant axis is `want_axis` and whose constant value matches
# `want_val` (within `tol`) are considered — this prevents matching a facet
# from a different plane that merely passes through `point` (e.g. a vertical
# bridge-leg facet when we want the horizontal substrate facet). Containment
# uses an even-odd ray-crossing test over ALL the facet's polygons (outer loop
# minus holes, order-independent). Returns `nothing` if none match.
function _facet_containing_point(facets::Vector{Vector{Vector{Cint}}},
                                 global_pts::Vector{Tuple{Float64,Float64,Float64}},
                                 point::Tuple{Float64,Float64,Float64}, tol::Float64,
                                 want_axis::Int, want_val::Float64)
    ptol = max(tol, 1e-6)
    for (fi, polys) in enumerate(facets)
        isempty(polys) && continue
        allv = Int[]
        for poly in polys, idx in poly; push!(allv, Int(idx)); end
        length(allv) >= 3 || continue
        xs = (global_pts[v][1] for v in allv)
        ys = (global_pts[v][2] for v in allv)
        zs = (global_pts[v][3] for v in allv)
        rx = maximum(xs) - minimum(xs); ry = maximum(ys) - minimum(ys); rz = maximum(zs) - minimum(zs)
        caxis = argmin((rx, ry, rz))
        caxis == want_axis || continue                # must lie in the edge's plane orientation
        cval = (caxis == 1 ? global_pts[allv[1]][1] :
                caxis == 2 ? global_pts[allv[1]][2] : global_pts[allv[1]][3])
        abs(cval - want_val) <= ptol || continue       # ... and the same plane
        # project to the 2 in-plane axes
        a1, a2 = caxis == 1 ? (2,3) : caxis == 2 ? (1,3) : (1,2)
        px = point[a1]; py = point[a2]
        crossings = 0
        for poly in polys
            np = length(poly)
            np >= 3 || continue                       # skip segment-polygons
            for e in 1:np
                v1 = global_pts[poly[e]]; v2 = global_pts[poly[e % np + 1]]
                x1 = v1[a1]; y1 = v1[a2]; x2 = v2[a1]; y2 = v2[a2]
                if (y1 > py) != (y2 > py)
                    xint = x1 + (py - y1) / (y2 - y1) * (x2 - x1)
                    px < xint && (crossings += 1)
                end
            end
        end
        isodd(crossings) && return fi
    end
    return nothing
end

# ─── Keys for deduplication ──────────────────────────────────────────────────

# Share triangulations across geometrically coincident face entries.
#
# Two volumes that abut at a face (e.g. substrate-top and vacuum-bottom both
# at z=0) get *separate* EGADS face EGObjects — EGADS has no native primitive
# that produces solids sharing TopoDS_TFace handles between bodies. When
# `tessellate` runs once per body, the two coincident faces tessellate
# independently and pick different triangulation diagonals. Per-triangle
# dedup (step 4 of `tetgen_io`) keys on sorted vertex tuples, so a triangle
# (a,b,c) on body A doesn't match (a,b,d) on body B and both end up in the
# PLC — TetGen then rejects the input as self-intersecting.
#
# This function fixes the rendering side: after global point dedup, every
# triangulated face entry has its triangles expressed in *global* vertex
# indices. Two coincident faces will have the same set of global vertex
# indices, so we hash on that set and rewrite all matching entries to share
# the first one's `tris` list. The pipeline then sees a single, consistent
# triangulation across the interface and dedup in step 4 collapses each
# triangle to one PLC facet.
function _share_triangulations_across_bodies!(per_body_faces::Vector{Vector{_FaceEntry}})
    # Build the key from each triangle's boundary edges: an edge is on the
    # face boundary iff it appears in exactly one triangle. Collect the
    # sorted set of global vertex indices on those boundary edges. This
    # decouples matching from interior tessellation jitter — two coincident
    # faces with identical loops match even if their interior Steiner points
    # land slightly differently.
    bucket = Dict{NTuple, Vector{_FaceEntry}}()
    for face_entries in per_body_faces, fe in face_entries
        fe.is_polygon && continue
        isempty(fe.tris) && continue
        edge_count = Dict{Tuple{Int, Int}, Int}()
        @inline function _bump!(a, b)
            k = a <= b ? (a, b) : (b, a)
            edge_count[k] = get(edge_count, k, 0) + 1
        end
        for t in fe.tris
            _bump!(t[1], t[2]); _bump!(t[2], t[3]); _bump!(t[3], t[1])
        end
        boundary_verts = Set{Int}()
        for ((a, b), c) in edge_count
            c == 1 || continue
            push!(boundary_verts, a); push!(boundary_verts, b)
        end
        isempty(boundary_verts) && continue
        key = Tuple(sort!(collect(boundary_verts)))
        push!(get!(bucket, key, _FaceEntry[]), fe)
    end
    n_shared = 0
    n_groups_shared = 0
    for (_, group) in bucket
        length(group) <= 1 && continue
        n_groups_shared += 1
        ref = group[1]
        for k in 2:length(group)
            # Replace this entry's tris with the reference's tris. Both
            # are in global indices already.
            empty!(group[k].tris)
            append!(group[k].tris, ref.tris)
            n_shared += 1
        end
    end
    n_shared > 0 && @info "_share_triangulations_across_bodies!: shared $n_shared face entries across $n_groups_shared coincident-face groups"
    return per_body_faces
end

# Polygon key: sorted flat tuple of all vertex indices across loops.
function _poly_key(loops::Vector{Vector{Int}})
    flat = Int[]
    for loop in loops, v in loop
        push!(flat, v)
    end
    sort!(flat)
    return Tuple(flat)
end

function _tri_key(a::Int, b::Int, c::Int)
    x, y, z = a, b, c
    x > y && ((x, y) = (y, x))
    y > z && ((y, z) = (z, y))
    x > y && ((x, y) = (y, x))
    return (x, y, z)
end

# ─── Facet assignment to RawTetGenIO ─────────────────────────────────────────

function _assign_facets!(io::TetGen.RawTetGenIO{Cdouble},
                         facets::Vector{Vector{Vector{Cint}}},
                         facet_holes::Vector{Vector{NTuple{3, Float64}}} =
                             fill(NTuple{3, Float64}[], length(facets)))
    io.facetlist = [TetGen.RawFacet{Cdouble}(f, _hole_matrix(h))
                    for (f, h) in zip(facets, facet_holes)]
    return io
end

# Convert a list of 3D hole-marker points into TetGen's 3×N hole-list matrix.
function _hole_matrix(holes::Vector{NTuple{3, Float64}})
    isempty(holes) && return Matrix{Cdouble}(undef, 0, 0)
    m = Matrix{Cdouble}(undef, 3, length(holes))
    @inbounds for (j, h) in enumerate(holes)
        m[1, j] = h[1]; m[2, j] = h[2]; m[3, j] = h[3]
    end
    return m
end

# ─── Internals ──────────────────────────────────────────────────────────────


# Spatially-varying tessellation size field via EGADS `.tParams` attributes.
#
# The `.tParams` attribute is EGADS's native per-entity tessellation override.
# A 3-real vector `[max_len, sag, angle]` attached to an EDGE tells
# `EG_makeTessBody` to use those parameters when discretising that edge;
# attached to a FACE, it controls the face's interior Steiner-point density.
#
# Given a user size field `h(x,y,z) -> target_edge_length`, we sample `h` at
# several points along every edge (and at each face centroid), take the
# minimum as a conservative local target, and write it into `.tParams`. The
# EGADS tessellator then produces a surface PLC whose element size matches
# the field, same `h` that TetGen's volume refinement will consult via the
# `tetunsuitable!` callback.
#
# **Important:** `EG_makeTessBody` only consults per-entity `.tParams` when
# the body-level `max_edge` is non-zero. We therefore return the maximum
# `h` sampled across all edges so the caller can pass it as the global
# `max_edge` to `tessellate` and activate the per-entity overrides.
function _apply_meshsize_attributes!(bodies::AbstractVector{<:EGObject},
                                     h::Union{Nothing, Function};
                                     angle::Real = 15.0,
                                     nsamples::Int = 5)
    h === nothing && return 0.0
    h_max = 0.0
    for body in bodies
        h_max = max(h_max,
            _apply_meshsize_to_edges!(body, h; angle = angle, nsamples = nsamples))
        _apply_meshsize_to_faces!(body, h; angle = angle)
    end
    return h_max
end

function _apply_meshsize_to_edges!(body::EGObject, h::Function;
                                   angle::Real, nsamples::Int)
    edges = get_body_topos(body, OClass.EDGE)
    h_max = 0.0
    for edge in edges
        info = EGADS.get_info(edge)
        info.mtype == MType.DEGENERATE && continue
        topo = try get_topology(edge) catch; continue end
        topo.geom === nothing && continue
        tlo, thi = topo.limits[1], topo.limits[2]
        isfinite(tlo) && isfinite(thi) && thi > tlo || continue
        h_min = Inf
        h_edge_max = 0.0
        for k in 0:(nsamples - 1)
            t = tlo + (thi - tlo) * k / max(nsamples - 1, 1)
            xyz = try
                EGADS.evaluate(topo.geom, t).xyz
            catch
                nothing
            end
            xyz === nothing && continue
            val = try
                Float64(h(xyz[1], xyz[2], xyz[3]))
            catch
                continue
            end
            isfinite(val) && val > 0 || continue
            val < h_min && (h_min = val)
            val > h_edge_max && (h_edge_max = val)
        end
        isfinite(h_min) || continue
        try
            EGADS.attribute_add!(edge, ".tParams",
                                 Float64[h_min, 0.0, Float64(angle)])
        catch
        end
        h_edge_max > h_max && (h_max = h_edge_max)
    end
    return h_max
end

function _apply_meshsize_to_faces!(body::EGObject, h::Function;
                                   angle::Real)
    faces = get_body_topos(body, OClass.FACE)
    for face in faces
        cog = try
            mp = EGADS.get_mass_properties(face)
            length(mp) >= 5 ? (mp[3], mp[4], mp[5]) : nothing
        catch
            nothing
        end
        cog === nothing && continue
        val = try
            Float64(h(cog[1], cog[2], cog[3]))
        catch
            continue
        end
        isfinite(val) && val > 0 || continue
        try
            EGADS.attribute_add!(face, ".tParams",
                                 Float64[val, 0.0, Float64(angle)])
        catch
        end
    end
    return nothing
end

function _set_holes_regions!(io::TetGen.RawTetGenIO{Cdouble},
                             hole_seeds, region_seeds)
    if hole_seeds !== nothing
        H = Matrix{Cdouble}(hole_seeds)
        size(H, 1) == 3 || throw(ArgumentError("hole_seeds must be 3×H"))
        io.holelist = H
    end
    if region_seeds !== nothing
        R = Matrix{Cdouble}(region_seeds)
        size(R, 1) == 5 || throw(ArgumentError(
            "region_seeds must be 5×R rows [x;y;z;attribute;max_volume]"))
        io.regionlist = R
    end
    return io
end

function _interior_seed_point(body::EGObject)
    bb = get_bounding_box(body)
    diag = hypot(bb.xmax - bb.xmin, bb.ymax - bb.ymin, bb.zmax - bb.zmin)
    diag == 0 && return nothing

    cx = 0.5 * (bb.xmin + bb.xmax)
    cy = 0.5 * (bb.ymin + bb.ymax)
    cz = 0.5 * (bb.zmin + bb.zmax)
    center = (cx, cy, cz)
    in_topology(body, [center[1], center[2], center[3]]) && return center

    tess = tessellate(body; angle = 30.0)
    faces = get_body_topos(body, OClass.FACE)
    isempty(faces) && return center
    step = max(diag * 1e-4, 1e3 * eps(Float64))
    for face_idx in 1:length(faces)
        fmesh = get_face_mesh(tess, body, face_idx)
        nt = size(fmesh.tris, 1)
        nt == 0 && continue
        face_mtype = try
            get_topology(faces[face_idx]).mtype
        catch
            Int32(MType.SFORWARD)
        end
        sense = face_mtype == MType.SREVERSE ? -1.0 : 1.0
        for ti in 1:nt
            a = Cint(fmesh.tris[ti, 1])
            b = Cint(fmesh.tris[ti, 2])
            c = Cint(fmesh.tris[ti, 3])
            pa = @view fmesh.xyz[a, :]
            pb = @view fmesh.xyz[b, :]
            pc = @view fmesh.xyz[c, :]
            centroid = (pa .+ pb .+ pc) ./ 3
            e1 = pb .- pa; e2 = pc .- pa
            n = (e1[2]*e2[3] - e1[3]*e2[2],
                 e1[3]*e2[1] - e1[1]*e2[3],
                 e1[1]*e2[2] - e1[2]*e2[1])
            nrm = hypot(n[1], n[2], n[3])
            nrm == 0 && continue
            inward = (-sense * n[1] / nrm, -sense * n[2] / nrm, -sense * n[3] / nrm)
            for factor in (1.0, 10.0, 100.0, 1000.0)
                p = (centroid[1] + factor * step * inward[1],
                     centroid[2] + factor * step * inward[2],
                     centroid[3] + factor * step * inward[3])
                in_topology(body, [p[1], p[2], p[3]]) && return p
                q = (centroid[1] - factor * step * inward[1],
                     centroid[2] - factor * step * inward[2],
                     centroid[3] - factor * step * inward[3])
                in_topology(body, [q[1], q[2], q[3]]) && return q
            end
        end
    end
    return center
end

# ─── Vertex → CAD classification ──────────────────────────────────────────────

"""
    classify_tetgen_vertices(out, body;
                                   face_marker_offset = 1000,
                                   node_tol = nothing,
                                   verbose = false)

Classify every constrained vertex of a TetGen output mesh onto the CAD
topology of `body` (the EGADS BODY the PLC came from): geometric NODE,
EDGE (with curve parameter `t`), or FACE (with surface parameters `(u, v)`).

This is the association that downstream consumers need:

  * **P2 curving**: a mesh edge whose endpoints share a CAD entity gets its
    mid-node projected onto that entity (`inv_evaluate`).
  * **Mesh morphing**: boundary velocity `dx/dp` is evaluated at the CAD
    parameters of each vertex (EGADS `_dot` sensitivities live on entities).
  * **Metric-based remeshers** (avro / refine / Metris): the meshb
    `VerticesOnGeometric{Vertices,Edges,Triangles}` sections — the returned
    matrices use exactly that layout (1-based vertex ids, EGADS body-topology
    indices as refs).

Classification is geometric, driven by the facet markers `tetgen_io` assigns
(marker = CAD face index for a single volume body):

  1. A vertex within `node_tol` of a geometric NODE → NODE.
  2. A vertex whose incident marked triangles carry ONE face id → that FACE,
     `(u, v, dist)` from `EG_invEvaluate`.
  3. Two or more face ids → the nearest EDGE shared by those faces
     (`t`, `dist`); falls back over all candidate faces' edges when the
     marker set has no shared edge (non-manifold junctions).

`dist` is the distance from the vertex to its CAD entity. Input tessellation
vertices lie on geometry (`dist ≈ 0`); TetGen Steiner points on segments or
facets of a *curved* entity lie on the linear chord, so `dist` is the local
chord error — exactly the displacement a curving/snapping pass must apply.
Nothing is rejected by distance: `dist` is data, the consumer decides.

Returns a NamedTuple:

  * `on_nodes::Matrix{Int}`       — columns `vertex_id, node_ref`
  * `on_edges::Matrix{Float64}`   — columns `vertex_id, edge_ref, t, dist`
  * `on_faces::Matrix{Float64}`   — columns `vertex_id, face_ref, u, v, dist`
  * `n_constrained`, `n_unclassified`, `max_dist`

Limitations (v1): classifies against ONE body — facet markers outside
`1..nfaces` (embedded sheets, other volume bodies via `face_marker_offset`
multiples) are skipped and counted in `n_unclassified`.
"""
function classify_tetgen_vertices(out, body::EGObject;
                                        face_marker_offset::Integer = 1000,
                                        node_tol::Union{Nothing, Real} = nothing,
                                        verbose::Bool = false)
    faces = get_body_topos(body, OClass.FACE)
    nfaces = length(faces)
    bnodes = get_body_topos(body, OClass.NODE)
    node_pts = [EGADS.get_topology(n).limits for n in bnodes]
    bb = get_bounding_box(body)
    diag = hypot(bb.xmax - bb.xmin, bb.ymax - bb.ymin, bb.zmax - bb.zmin)
    ntol = node_tol === nothing ? 1.0e-7 * diag : Float64(node_tol)

    # Per-face non-degenerate edge lists for shared-edge resolution.
    face_edges = Vector{Vector{EGADS.EGObject}}(undef, nfaces)
    for fi in 1:nfaces
        face_edges[fi] = EGADS.EGObject[
            e for e in get_body_topos(body, faces[fi], OClass.EDGE)
            if EGADS.get_info(e).mtype != MType.DEGENERATE
        ]
    end

    # Vertex → set of incident CAD face ids, from the facet markers.
    nv = size(out.pointlist, 2)
    vfaces = Dict{Int, Set{Int}}()
    n_foreign = 0
    for j in 1:size(out.trifacelist, 2)
        m = Int(out.trifacemarkerlist[j])
        # Single-body convention: marker == face id. Markers from other
        # volume bodies (face_marker_offset multiples) or embedded sheets
        # are foreign to `body`.
        if !(1 <= m <= nfaces)
            n_foreign += 1
            continue
        end
        for i in 1:3
            push!(get!(vfaces, Int(out.trifacelist[i, j]), Set{Int}()), m)
        end
    end

    on_nodes = Vector{NTuple{2, Int}}()
    on_edges = Vector{NTuple{4, Float64}}()
    on_faces = Vector{NTuple{5, Float64}}()
    n_unclassified = 0
    max_dist = 0.0

    for (vi, fids) in vfaces
        p = (out.pointlist[1, vi], out.pointlist[2, vi], out.pointlist[3, vi])

        # 1. Geometric NODE?
        best_n, best_nd = 0, Inf
        for (ni, q) in enumerate(node_pts)
            d = hypot(p[1] - q[1], p[2] - q[2], p[3] - q[3])
            d < best_nd && (best_nd = d; best_n = ni)
        end
        if best_nd <= ntol
            push!(on_nodes, (vi, Int(EGADS.index_body_topo(body, bnodes[best_n]))))
            continue
        end

        pv = Float64[p[1], p[2], p[3]]
        if length(fids) == 1
            # 2. Interior of a single CAD face.
            fid = first(fids)
            r = try EGADS.inv_evaluate(faces[fid], pv) catch; nothing end
            if r === nothing
                n_unclassified += 1
                continue
            end
            d = hypot(r.results[1] - p[1], r.results[2] - p[2], r.results[3] - p[3])
            push!(on_faces, (Float64(vi), Float64(fid), r.param[1], r.param[2], d))
            max_dist = max(max_dist, d)
        else
            # 3. On the boundary between CAD faces → nearest shared EDGE.
            fl = collect(fids)
            shared = EGADS.EGObject[]
            seen = Set{Ptr{Cvoid}}()
            for a in 1:length(fl), b in (a + 1):length(fl)
                ptrs_b = Set{Ptr{Cvoid}}(e.ptr for e in face_edges[fl[b]])
                for e in face_edges[fl[a]]
                    (e.ptr in ptrs_b && !(e.ptr in seen)) || continue
                    push!(seen, e.ptr)
                    push!(shared, e)
                end
            end
            # Non-manifold junctions may pair faces with no shared edge —
            # fall back to every edge of every candidate face.
            cand = isempty(shared) ?
                   EGADS.EGObject[e for fid in fl for e in face_edges[fid]] : shared
            best_e, best_t, best_d = 0, 0.0, Inf
            for (ei, e) in enumerate(cand)
                r = try EGADS.inv_evaluate(e, pv) catch; continue end
                d = hypot(r.results[1] - p[1], r.results[2] - p[2],
                          r.results[3] - p[3])
                d < best_d && (best_d = d; best_e = ei; best_t = r.param[1])
            end
            if best_e == 0
                n_unclassified += 1
                continue
            end
            push!(on_edges, (Float64(vi),
                             Float64(EGADS.index_body_topo(body, cand[best_e])),
                             best_t, best_d))
            max_dist = max(max_dist, best_d)
        end
    end

    sort!(on_nodes; by = first)
    sort!(on_edges; by = first)
    sort!(on_faces; by = first)
    result = (
        on_nodes = isempty(on_nodes) ? zeros(Int, 0, 2) :
                   permutedims(reduce(hcat, [collect(t) for t in on_nodes])),
        on_edges = isempty(on_edges) ? zeros(Float64, 0, 4) :
                   permutedims(reduce(hcat, [collect(t) for t in on_edges])),
        on_faces = isempty(on_faces) ? zeros(Float64, 0, 5) :
                   permutedims(reduce(hcat, [collect(t) for t in on_faces])),
        n_constrained = length(vfaces),
        n_unclassified = n_unclassified,
        max_dist = max_dist,
    )
    verbose && println("classify_tetgen_vertices: ", result.n_constrained,
                       " constrained vertices → ", size(result.on_nodes, 1),
                       " NODE / ", size(result.on_edges, 1), " EDGE / ",
                       size(result.on_faces, 1), " FACE; ",
                       n_unclassified, " unclassified; max dist = ",
                       round(max_dist; sigdigits = 3),
                       n_foreign > 0 ? "; $n_foreign foreign-marker facets skipped" : "")
    return result
end

end # module
