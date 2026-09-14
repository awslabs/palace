# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Standalone defaults are diagnostic: no scalar trace sizing. Positive sizes need
# an explicit scope AND CAD/size segment mode; legacy is never inferred.
function trace_size_policy(environment=ENV)
    maximum_size = parse(Float64, get(environment,"TET_TRACE_SIZE","0"))
    relative_size = parse(Float64, get(environment,"TET_TRACE_RELATIVE_SIZE","0"))
    scope = get(environment,"TET_TRACE_SIZE_SCOPE","off")
    mode = get(environment,"TET_TRACE_CONSTRAINT_MODE","all")
    scope in ("off","matching","matching-and-volume","legacy-global") ||
        error("Invalid trace-size scope")
    mode in ("all","sides","levels","none") || error("Invalid trace constraint mode")
    isfinite(maximum_size) && maximum_size >= 0 || error("Invalid trace size")
    isfinite(relative_size) && relative_size >= 0 || error("Invalid relative trace size")
    maximum_size == 0 || scope != "off" ||
        error("Positive TET_TRACE_SIZE requires explicit TET_TRACE_SIZE_SCOPE")
    maximum_size > 0 || scope == "off" || error("A trace-size scope requires a positive TET_TRACE_SIZE")
    maximum_size == 0 || haskey(environment,"TET_TRACE_CONSTRAINT_MODE") ||
        error("Positive TET_TRACE_SIZE requires explicit TET_TRACE_CONSTRAINT_MODE")
    mode != "none" || scope == "off" || error("Trace mode none requires scalar trace sizing off")
    relative_size == 0 || maximum_size > 0 || error("Relative trace sizing requires a positive maximum trace size")
    surface_growth = parse(Float64, get(environment,"TET_TRACE_SURFACE_GROWTH",
                                       scope == "legacy-global" ? "0.5" : "4"))
    isfinite(surface_growth) && surface_growth > 0 || error("Invalid trace-surface growth")
    scope != "legacy-global" || surface_growth == 0.5 ||
        error("Legacy-global replay requires trace surface growth 0.5")
    return (;maximum_size,relative_size,scope=Symbol(replace(scope,"-"=>"_")),mode,surface_growth)
end

# Trace topology and optional scalar sizing are separate. Per-edge altitude is a
# geometric heuristic, NOT a function-resolution or accuracy bound. Long edges of
# skinny triangles can still request tiny isotropic sizes. Legacy replay assigns
# the triangle minimum altitude to EVERY edge and ignores CAD segment selection.
function trace_segment_sizes(triangles, maximum_size, relative_size; mode="all", legacy=false)
    isfinite(maximum_size) && maximum_size > 0 || error("Invalid maximum trace size")
    isfinite(relative_size) && relative_size >= 0 || error("Invalid relative trace size")
    mode in ("all", "sides", "levels") || error("Unknown matching trace constraint mode")
    areas = trace_triangle_areas(triangles)
    vertices = [point for triangle in values(triangles) for point in triangle]
    lower = ntuple(d -> minimum(point[d] for point in vertices), 3)
    upper = ntuple(d -> maximum(point[d] for point in vertices), 3)
    tolerance = 1.0e-10 * max(maximum(upper[d] - lower[d] for d in 1:3), 1.0)
    if mode == "levels" && !legacy
        levels = sort!(unique(point[3] for point in vertices))
        edges = Tuple{NTuple{3,Float64},NTuple{3,Float64}}[]
        for z in levels
            lower[3] + tolerance < z < upper[3] - tolerance || continue
            ring = [(lower[1],lower[2],z), (upper[1],lower[2],z),
                    (upper[1],upper[2],z), (lower[1],upper[2],z)]
            append!(edges, [(ring[i],ring[mod1(i+1,4)]) for i in 1:4])
        end
        isempty(edges) && error("Trace levels contain no interior matching-box rings")
        return NTuple{6,Float64}[(first...,last...) for (first,last) in edges],
               fill(maximum_size,length(edges))
    end
    scales = Dict{Tuple{NTuple{3,Float64},NTuple{3,Float64}},Float64}()
    dimensions = mode == "sides" ? (1:2) : (1:3)
    for (index,triangle) in triangles
        edge_lengths = [sqrt(sum((triangle[mod1(i+1,3)][d]-triangle[i][d])^2
                                for d in 1:3)) for i in 1:3]
        all(h -> isfinite(h) && h > 0, edge_lengths) || error("Invalid trace edge")
        minimum_altitude = areas[index] / maximum(edge_lengths)
        for i in 1:3
            first, last = triangle[i], triangle[mod1(i+1,3)]
            on_box = any((abs(first[d]-lower[d])<=tolerance && abs(last[d]-lower[d])<=tolerance) ||
                         (abs(first[d]-upper[d])<=tolerance && abs(last[d]-upper[d])<=tolerance)
                         for d in dimensions)
            legacy || on_box || continue
            key = isless(first, last) ? (first, last) : (last, first)
            scale = legacy ? minimum_altitude : areas[index] / edge_lengths[i]
            scales[key] = min(get(scales, key, Inf), scale)
        end
    end
    isempty(scales) && error("Trace geometry has no edges for mode $mode")
    edges = sort!(collect(keys(scales)))
    segments = NTuple{6,Float64}[(first..., last...) for (first, last) in edges]
    sizes = [relative_size == 0 ? maximum_size : min(maximum_size, relative_size * scales[edge])
             for edge in edges]
    all(size -> isfinite(size) && size > 0, sizes) || error("Invalid trace feature scale")
    return segments, sizes
end
