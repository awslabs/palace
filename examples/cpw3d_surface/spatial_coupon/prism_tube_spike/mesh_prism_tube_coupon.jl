# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Prism edge tube SPIKE (supervisor decision 37): the four-edge spatial coupon with
# a localized prism tube around every straight metal edge (top and bottom edge of
# each metal segment), Gmsh only (no MMG), the production size fields and corner
# grading everywhere else. NOT production: the canonical pipeline
# (mesh_spatial_coupon.jl, manifests, gates) is untouched; this script reuses its
# CAD helpers and reproduces its fabricated sharp-sidewall CAD path for the
# single-slot single-conductor coupon with an explicit etch footprint.
#
# usage: julia --project=test/examples mesh_prism_tube_coupon.jl <case directory> output.msh
#        [--tube-inner-size 0.00025] [--tube-ratio 2] [--tube-rings 7] [--tube-sector-degrees 30]
#        [--pyramid-height 0.008] [--census census.json]

include(joinpath(@__DIR__, "..", "mesh_spatial_coupon.jl"))
include(joinpath(@__DIR__, "prism_tube.jl"))
include(joinpath(@__DIR__, "hybrid_census.jl"))
using TOML

# Straight metal edge segments of the plan-view boundary loops: sides not lying on
# the outer box, with the outward horizontal normal (away from the metal) and the
# tube interval shrunk by the corner ball radius at semantic corners (0 at box
# continuation vertices). Returns (start, stop, normal, s_start, s_end) in the plane.
function metal_edge_segments(loops, corners, radius, lower, upper, tolerance)
    segments = NamedTuple[]
    on_box(p) = any(abs(p[d] - lower[d]) <= tolerance || abs(p[d] - upper[d]) <= tolerance
                    for d in 1:2)
    same_box_face(p, q) = any((abs(p[d] - lower[d]) <= tolerance && abs(q[d] - lower[d]) <= tolerance) ||
                              (abs(p[d] - upper[d]) <= tolerance && abs(q[d] - upper[d]) <= tolerance)
                              for d in 1:2)
    is_corner(p, plane) = any(norm(collect(corner) .- [p[1], p[2], plane]) <= tolerance
                              for corner in corners)
    for loop in loops
        loop.hole && error("Prism tube spike supports exterior conductor loops only")
        n = length(loop.points)
        for i in 1:n
            p = loop.points[i]
            q = loop.points[i % n + 1]
            same_box_face(p, q) && continue
            direction = [q[1] - p[1], q[2] - p[2]]
            span = norm(direction)
            direction ./= span
            normal = [direction[2], -direction[1]]
            midpoint = [0.5 * (p[1] + q[1]), 0.5 * (p[2] + q[2])]
            probe = midpoint .+ 1.0e-3 .* normal
            if point_in_polygon((probe[1], probe[2]), loop.points, tolerance)
                normal .*= -1.0
            end
            probe = midpoint .- 1.0e-3 .* normal
            point_in_polygon((probe[1], probe[2]), loop.points, tolerance) ||
                error("Unable to orient the metal edge $p -> $q")
            shrink_start = is_corner(p, loop.plane) ? radius : (on_box(p) ? 0.0 :
                error("Metal edge end $p is neither a semantic corner nor on the box"))
            shrink_end = is_corner(q, loop.plane) ? radius : (on_box(q) ? 0.0 :
                error("Metal edge end $q is neither a semantic corner nor on the box"))
            push!(segments, (start=[p[1], p[2]], stop=[q[1], q[2]], direction=direction,
                             normal=normal, span=span, s_start=shrink_start,
                             s_end=span - shrink_end, plane=loop.plane, conductor=loop.conductor))
        end
    end
    isempty(segments) && error("No straight metal edges found")
    return segments
end

# The etch footprint must carry the metal edge (trench wall under the sidewall) so
# that the bottom tube's substrate quadrant / vacuum half split is exact.
function assert_etch_carries_edge(etch_loops, segment, tolerance)
    for loop in etch_loops
        n = length(loop.points)
        for i in 1:n
            p = loop.points[i]
            q = loop.points[i % n + 1]
            for (a, b) in ((p, q), (q, p))
                norm([a[1], a[2]] .- segment.start) <= tolerance &&
                    norm([b[1], b[2]] .- segment.stop) <= tolerance && return true
            end
        end
    end
    error("Etch footprint has no side coincident with the metal edge $(segment.start) -> $(segment.stop)")
end

function generate_prism_tube_coupon(; case_directory::String, filename::String,
                                    tube_inner_size::Float64=0.00025, tube_ratio::Float64=2.0,
                                    tube_rings::Int=7, tube_sector_degrees::Float64=30.0,
                                    pyramid_height::Float64=0.008,
                                    lc_tangent::Float64=0.05, corner_size::Float64=0.004,
                                    edge_growth_ratio::Float64=2.0,
                                    trace_basis_size_ratio::Float64=1.0,
                                    max_elements::Int=4_000_000,
                                    census_path::Union{Nothing, String}=nothing)
    inputs(name) = joinpath(case_directory, name)
    process = TOML.parsefile(inputs("process.toml"))
    process["Units"] == "um" || error("Prism tube spike expects um")
    recipe = parse_json(read(joinpath(@__DIR__, "..", "testdata", "generality-mesh-recipe.json"), String))
    radius = Float64(process["Radius"])
    metal_thickness = Float64(process["MetalThickness"])
    overetch = Float64(process["Overetch"])
    Float64(process["SidewallAngle"]) == 90.0 && Float64(process["TopRounding"]) == 0.0 &&
        Float64(process["TrenchRounding"]) == 0.0 ||
        error("Prism tube spike requires sharp vertical fabricated geometry")
    lc_fine = Float64(recipe["NormalSizeOverThickness"]) * metal_thickness
    lc_far = Float64(recipe["FarSizeOverRadius"]) * radius
    corner_isotropy_radius = Float64(recipe["TangentialSizeOverNormalSize"]) * lc_fine
    process_core_width = max(2metal_thickness, 4overetch, 8lc_fine)
    process_fine_width = 0.0
    corner_grading_slope = (lc_far - lc_fine) / (process_core_width - process_fine_width)
    corner_grading = CornerGrading(corner_size, edge_growth_ratio, lc_fine, corner_isotropy_radius)
    corner_grading_reach(corner_grading) <= corner_isotropy_radius ||
        error("The corner grading must reach the normal size inside the corner ball")

    edges = read_edges(inputs("mesh-signature.csv"))
    facets = read_mask(inputs("plan-view-mask.csv"))
    boundary_loops = read_boundary(inputs("plan-view-boundary.csv"))
    etch_loops = read_boundary(inputs("retained-etch.csv"))
    transform = copy(IDENTITY_RIGID_TRANSFORM)
    semantic_corners = read_semantic_corners(inputs("semantic-contract.json"), transform)
    trace_basis = read_trace_basis(inputs("basis-contract.json"), inputs("trace-vertices.csv"),
                                   inputs("trace-triangles.csv"), inputs("process-library.json"))
    lower, upper = coupon_bounds(edges, radius, metal_thickness, overetch)
    tolerance = 1.0e-7 * radius
    outer_tolerance = 1.0e-4 * radius
    all(abs(trace_basis.lower[d] - lower[d]) <= tolerance &&
        abs(trace_basis.upper[d] - upper[d]) <= tolerance for d in 1:3) ||
        error("Trace basis box differs from the coupon box")
    validate_plan_view_geometry(edges, radius, tolerance, facets)
    layers = layer_groups(edges, tolerance)
    length(layers) == 1 && layers[1].sign > 0 || error("Prism tube spike supports one upward layer")
    layer = layers[1]
    length(unique(edge.slot for edge in edges)) == 1 || error("single slot only")

    segments = metal_edge_segments(boundary_loops, semantic_corners, corner_isotropy_radius,
                                   lower, upper, tolerance)
    for segment in segments
        assert_etch_carries_edge(etch_loops, segment, tolerance)
    end
    sectors = round(Int, 270.0 / tube_sector_degrees)
    abs(sectors * tube_sector_degrees - 270.0) <= 1.0e-9 || error("sector angle must divide 270")
    per_quadrant = round(Int, 90.0 / tube_sector_degrees)
    abs(per_quadrant * tube_sector_degrees - 90.0) <= 1.0e-9 || error("sector angle must divide 90")
    # Top edge: vacuum from the sidewall ray (-90) through the outward normal to the
    # top face ray (180). Bottom edge: substrate from the metal bottom face ray (180)
    # to the trench wall ray (270), vacuum from the trench wall to the sidewall (450).
    top_section = TubeSection(tube_inner_size, tube_ratio, tube_rings,
                              [-90.0 + tube_sector_degrees * j for j in 0:sectors], fill(2, sectors))
    bottom_section = TubeSection(tube_inner_size, tube_ratio, tube_rings,
                                 [180.0 + tube_sector_degrees * j for j in 0:sectors],
                                 vcat(fill(1, per_quadrant), fill(2, sectors - per_quadrant)))
    tube_radius(top_section) < overetch ||
        error("tube radius $(tube_radius(top_section)) must stay above the trench floor ($overetch)")
    tube_radius(top_section) + pyramid_height < overetch ||
        error("pyramid apexes must stay above the trench floor")

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)
    gmsh.model.add("prism_tube_coupon")
    occ = gmsh.model.occ
    outer = (3, occ.addBox(lower[1], lower[2], lower[3], upper[1] - lower[1],
                           upper[2] - lower[2], upper[3] - lower[3]))
    slab = [(3, occ.addBox(lower[1], lower[2], lower[3], upper[1] - lower[1],
                           upper[2] - lower[2], layer.plane - lower[3]))]
    footprint_polygons = Dict{String, Any}[]
    footprint = [loop for loop in etch_loops if abs(loop.plane - layer.plane) <= tolerance]
    isempty(footprint) && error("Explicit etch footprint is missing the process layer")
    trench = loft_mask(occ, footprint, layer.plane, layer.plane - layer.sign * overetch, 0.0,
                       tolerance; simplify=true, footprint=footprint_polygons)
    substrates, _ = occ.cut(slab, trench)
    metal = Tuple{Int32, Int32}[]
    for conductor in sort!(unique(edge.conductor for edge in layer.edges))
        conductor_loops = [loop for loop in boundary_loops if loop.conductor == conductor &&
                           abs(loop.plane - layer.plane) <= tolerance]
        isempty(conductor_loops) && error("Plan-view boundary is missing conductor $conductor")
        append!(metal, loft_mask(occ, conductor_loops, layer.plane,
                                 layer.plane + layer.sign * metal_thickness, 0.0, tolerance))
    end
    field, _ = occ.cut([outer], metal, -1, true, true)
    vacuum, _ = occ.cut(field, substrates, -1, true, false)
    objects = vcat(substrates, vacuum)
    # Tubes: one per straight metal edge and ridge (top z = plane + thickness, bottom
    # z = plane), extruded along e = n x z from the segment end that makes e point
    # along the segment.
    tubes = Tuple{EdgeTube, TubeSection}[]
    records = TubeRecord[]
    for segment in segments
        n = [segment.normal[1], segment.normal[2], 0.0]
        b = [0.0, 0.0, 1.0]
        e = cross(n, b)
        along = dot(e[1:2], segment.direction)
        abs(abs(along) - 1.0) <= 1.0e-12 || error("tube frame is not aligned with the edge")
        for (z, section) in ((layer.plane + metal_thickness, top_section), (layer.plane, bottom_section))
            if along > 0.0
                origin = [segment.start[1], segment.start[2], z]
                tube = EdgeTube(origin, n, b, segment.s_start, segment.s_end, lc_tangent)
            else
                origin = [segment.stop[1], segment.stop[2], z]
                tube = EdgeTube(origin, n, b, segment.span - segment.s_end,
                                segment.span - segment.s_start, lc_tangent)
            end
            push!(tubes, (tube, section))
            for (tool, group) in add_tube_volumes!(occ, tube, section)
                push!(records, TubeRecord(tube, section, group, tool))
            end
        end
    end
    tools = [record.tool for record in records]
    domains, fragment_map = occ.fragment(objects, tools)
    occ.synchronize()
    substrate_tags = Int32[]
    vacuum_tags = Int32[]
    for (i, object) in enumerate(objects)
        descendants = [tag for (dim, tag) in fragment_map[i] if dim == 3]
        append!(object in substrates ? substrate_tags : vacuum_tags, descendants)
    end
    unique!(substrate_tags)
    unique!(vacuum_tags)
    tube_volumes = Int32[tube_volume_after_fragment(record, index, fragment_map, length(objects))
                         for (index, record) in enumerate(records)]
    all(volume in substrate_tags || volume in vacuum_tags for volume in tube_volumes) ||
        error("A tube volume is neither substrate nor vacuum")
    for (index, record) in enumerate(records)
        expected = record.group[3] == 1 ? substrate_tags : vacuum_tags
        tube_volumes[index] in expected ||
            error("Tube volume $(tube_volumes[index]) material differs from its section material")
    end
    corner_point_tags = semantic_corner_points(semantic_corners, tolerance)
    substrate_set = Set(substrate_tags)
    vacuum_set = Set(vacuum_tags)

    # Surface labels as the production mesher, with one addition: a face between
    # two volumes of the same material (tube / outside) is an interior face.
    matching = Int32[]
    boundary_groups = Dict{Int, Vector{Int32}}()
    interface_surfaces = Int32[]
    for (dim, tag) in gmsh.model.getEntities(2)
        up, _ = gmsh.model.getAdjacencies(dim, tag)
        adjacent_substrate = [volume for volume in up if volume in substrate_set]
        adjacent_vacuum = [volume for volume in up if volume in vacuum_set]
        isempty(adjacent_substrate) && isempty(adjacent_vacuum) && continue
        bounds = gmsh.model.getBoundingBox(dim, tag)
        if on_outer_box(bounds, lower, upper, outer_tolerance)
            push!(matching, tag)
            continue
        end
        point = point_on_surface(tag)
        edge = nearest_edge(edges, point, radius)
        attribute = 0
        if !isempty(adjacent_substrate) && !isempty(adjacent_vacuum)
            _, _, zmin, _, _, zmax = bounds
            attribute = abs(zmin - edge.point[3]) < tolerance && abs(zmax - edge.point[3]) < tolerance ?
                        3000 + edge.slot : 3100 + edge.slot
            push!(interface_surfaces, tag)
        elseif length(up) == 1 && !isempty(adjacent_substrate)
            owner = nearest_metal_edge(edges, facets, point, radius, tolerance)
            attribute = metal_surface_attribute(5000, owner.slot, owner.conductor)
        elseif length(up) == 1 && !isempty(adjacent_vacuum)
            owner = nearest_metal_edge(edges, facets, point, radius, tolerance)
            attribute = metal_surface_attribute(6000, owner.slot, owner.conductor)
        end
        attribute > 0 && push!(get!(boundary_groups, attribute, Int32[]), tag)
    end
    isempty(matching) && error("No matching surface was generated")
    gmsh.model.addPhysicalGroup(2, unique(matching), 1, "matching_surface")
    for (attribute, surfaces) in sort(collect(boundary_groups))
        unique!(surfaces)
        gmsh.model.addPhysicalGroup(2, surfaces, attribute, "surface_$attribute")
    end

    # Feature curves exactly as the production mesher (coplanar same-family seams
    # discarded, junction curves added).
    curve_surfaces = Dict{Int32, Vector{Tuple{Int, Int32}}}()
    candidate_curves = Set{Int32}()
    for (attribute, surfaces) in boundary_groups, surface in surfaces
        for (curve_dim, curve) in gmsh.model.getBoundary([(2, surface)], false, false, false)
            curve_dim == 1 || continue
            on_outer_box(gmsh.model.getBoundingBox(curve_dim, curve), lower, upper,
                         outer_tolerance) && continue
            push!(candidate_curves, curve)
            push!(get!(curve_surfaces, curve, Tuple{Int, Int32}[]), (attribute, surface))
        end
    end
    feature_curves = Int32[]
    discarded_seams = Int32[]
    for curve in candidate_curves
        records_of_curve = unique(curve_surfaces[curve])
        attributes = unique(first(record) for record in records_of_curve)
        surfaces = unique(last(record) for record in records_of_curve)
        families = unique(surface_family(attribute) for attribute in attributes)
        artificial_seam = length(surfaces) >= 2 && length(families) == 1 &&
                          coplanar_surfaces(surfaces, tolerance)
        push!(artificial_seam ? discarded_seams : feature_curves, curve)
    end
    junction_curves = Int32[]
    for surface in interface_surfaces
        for (curve_dim, curve) in gmsh.model.getBoundary([(2, surface)], false, false, false)
            curve_dim == 1 || continue
            on_outer_box(gmsh.model.getBoundingBox(curve_dim, curve), lower, upper,
                         outer_tolerance) || continue
            push!(junction_curves, curve)
        end
    end
    sort!(unique!(junction_curves))
    isempty(junction_curves) && error("No junction curves were generated")
    append!(feature_curves, junction_curves)
    sort!(unique!(feature_curves))

    # Explicit tube meshes, phase 1 (points, curves), before the other explicit
    # curve meshes so that the shared ridge / ball-boundary points are registered.
    next_node = Ref(0)
    point_nodes = Dict{Int32, Int}()
    states = TubeMesh[]
    tube_curves = Set{Int32}()
    for (tube, section) in tubes
        volumes = Tuple{Int32, Tuple{Int, Int, Int}, Dict}[]
        for (index, record) in enumerate(records)
            record.tube === tube || continue
            matched = match_tube_entities(tube_volumes[index], tube, section, record.group,
                                          1.0e-6 * radius)
            push!(volumes, (tube_volumes[index], record.group, matched))
            for ((kind, _), tag) in matched
                kind in (:edge_line, :outer_line, :cap_polygon, :cap_ray) && push!(tube_curves, tag)
            end
        end
        state = TubeMesh(tube, section, volumes; pyramid_height=pyramid_height)
        install_tube_curves!(state, next_node, point_nodes)
        push!(states, state)
    end
    # Longitudinal feature curves outside the tubes: the production placement
    # (corner-graded explicit nodes inside the ball law's reach, transfinite
    # lc_tangent grid otherwise).
    longitudinal_curves = Int32[]
    explicit_curves = Int32[]
    transfinite_count = 0
    for curve in feature_curves
        curve in tube_curves && continue
        lower_parameter, upper_parameter = gmsh.model.getParametrizationBounds(1, curve)
        parameter = 0.5 * (lower_parameter[1] + upper_parameter[1])
        tangent = gmsh.model.getDerivative(1, curve, [parameter])[1:3]
        tangent ./= norm(tangent)
        any(abs(dot(tangent, edge.tangent)) >= 1.0 - 1.0e-6 for edge in edges) || continue
        push!(longitudinal_curves, curve)
        placed = corner_isotropic_curve_nodes(curve, semantic_corners, corner_grading, lc_tangent,
                                              corner_grading_slope)
        if placed === nothing
            curve_length = gmsh.model.occ.getMass(1, curve)
            gmsh.model.mesh.setTransfiniteCurve(curve, max(2, ceil(Int, curve_length / lc_tangent) + 1))
            transfinite_count += 1
        else
            add_explicit_curve_mesh!(curve, placed[1], placed[2], next_node, point_nodes)
            push!(explicit_curves, curve)
        end
    end
    gmsh.option.setNumber("Mesh.MeshOnlyEmpty", 1)
    gmsh.option.setNumber("Mesh.Renumber", 0)

    # Production size fields: anisotropic curve attractor, corner law, trace basis.
    trace_basis_record = prepare_trace_basis_sizing!(trace_basis, trace_basis_size_ratio, lc_far,
                                                     corner_grading_slope)
    mesh_size_minimum = min(lc_fine, trace_basis_record["MinimumRequestedSize"], corner_size)
    install_trace_basis_callback!()
    gmsh.model.mesh.field.add("AttractorAnisoCurve", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", Float64.(feature_curves))
    gmsh.model.mesh.field.setNumber(1, "DistMin", process_fine_width)
    gmsh.model.mesh.field.setNumber(1, "DistMax", process_core_width)
    gmsh.model.mesh.field.setNumber(1, "SizeMinNormal", lc_fine)
    gmsh.model.mesh.field.setNumber(1, "SizeMaxNormal", lc_far)
    gmsh.model.mesh.field.setNumber(1, "SizeMinTangent", lc_tangent)
    gmsh.model.mesh.field.setNumber(1, "SizeMaxTangent", lc_far)
    gmsh.model.mesh.field.setNumber(1, "Sampling", 100)
    gmsh.model.mesh.field.add("Distance", 2)
    gmsh.model.mesh.field.setNumbers(2, "PointsList", Float64.(corner_point_tags))
    gmsh.model.mesh.field.add("MathEval", 3)
    gmsh.model.mesh.field.setString(3, "F", "min(F1," * corner_size_expression(
        "F2", corner_grading, lc_far, process_core_width - process_fine_width) * ")")
    gmsh.model.mesh.field.setAsBackgroundMesh(3)
    for (name, value) in [
        ("Mesh.MeshSizeMin", mesh_size_minimum), ("Mesh.MeshSizeMax", lc_far),
        ("Mesh.Algorithm3D", 1), ("Mesh.MeshSizeExtendFromBoundary", 0),
        ("Mesh.MeshSizeFromPoints", 0), ("Mesh.MeshSizeFromCurvature", 0),
        ("Mesh.MinimumCirclePoints", 24), ("Mesh.MinimumCurvePoints", 3),
        ("Mesh.MshFileVersion", 2.2), ("Mesh.Binary", 1)]
        gmsh.option.setNumber(name, value)
    end
    println("Prism tube spike: tubes=$(length(tubes)) (volumes $(length(records))), " *
            "rings=$(ring_sizes(top_section)), radius=$(tube_radius(top_section)), " *
            "sectors=$sectors x $(tube_sector_degrees) deg, pyramid_height=$pyramid_height, " *
            "layers=$(sum(tube.layers for (tube, _) in tubes)), " *
            "feature curves=$(length(feature_curves)) (junction $(length(junction_curves)), " *
            "seams discarded $(length(discarded_seams))), longitudinal outside tubes " *
            "$(length(longitudinal_curves)) (explicit $(length(explicit_curves)), " *
            "transfinite $transfinite_count), mesh size minimum=$mesh_size_minimum")
    timings = Dict{String, Float64}()
    timings["Generate2D"] = @elapsed gmsh.model.mesh.generate(2)
    timings["TubeFaces"] = @elapsed tube_faces = [install_tube_faces!(state, next_node) for state in states]
    remove_tube_volumes!(states)
    timings["Generate3D"] = @elapsed gmsh.model.mesh.generate(3)
    gmsh.model.mesh.removeSizeCallback()
    duplicates = remove_duplicate_volume_elements!()
    tube_discrete, tube_census = finalize_tube_volumes!(states)
    outside_substrate = setdiff(substrate_tags, tube_volumes)
    outside_vacuum = setdiff(vacuum_tags, tube_volumes)
    gmsh.model.addPhysicalGroup(3, vcat(outside_substrate, get(tube_discrete, 1, Int32[])), 1, "substrate")
    gmsh.model.addPhysicalGroup(3, vcat(outside_vacuum, get(tube_discrete, 2, Int32[])), 2, "vacuum")
    census = hybrid_volume_census()
    census["Total"] <= max_elements ||
        error("Prism tube coupon exceeds the element budget: $(census["Total"]) > $max_elements")
    census["Areas"] = physical_surface_areas()
    census["Tubes"] = tube_census
    census["TubeFaces"] = tube_faces
    census["Timings"] = timings
    census["DuplicateVolumeElementsRemoved"] = duplicates
    census["Anisotropy"] = tube_anisotropy(tubes[1][1], top_section, pyramid_height)
    census["TubeSpanLength"] = sum(tube.s_end - tube.s_start for (tube, _) in tubes)
    census["Sections"] = Dict{String, Any}(
        "Top" => Dict("Angles" => top_section.angles, "Materials" => top_section.materials),
        "Bottom" => Dict("Angles" => bottom_section.angles, "Materials" => bottom_section.materials))
    census["Production"] = Dict{String, Any}(
        "NormalSize" => lc_fine, "TangentialSize" => lc_tangent, "FarSize" => lc_far,
        "CornerIsotropyRadius" => corner_isotropy_radius, "CornerSize" => corner_size,
        "GrowthRatio" => edge_growth_ratio, "TraceBasisSizeRatio" => trace_basis_size_ratio,
        "ProcessCoreWidth" => process_core_width, "EdgeLayer" => "replaced by the prism tubes")
    census["MaxNodeTag"] = Int(gmsh.model.mesh.getMaxNodeTag())
    gmsh.write(filename)
    gmsh.finalize()
    if census_path !== nothing
        open(census_path, "w") do io
            write_census_json(io, census)
        end
    end
    return census
end

function parse_spike_options(args)
    options = Dict{Symbol, Any}()
    index = 3
    names = Dict("--tube-inner-size" => (:tube_inner_size, Float64), "--tube-ratio" => (:tube_ratio, Float64),
                 "--tube-rings" => (:tube_rings, Int), "--tube-sector-degrees" => (:tube_sector_degrees, Float64),
                 "--pyramid-height" => (:pyramid_height, Float64), "--max-elements" => (:max_elements, Int))
    while index <= length(args)
        if args[index] == "--census"
            options[:census_path] = abspath(args[index + 1])
        elseif haskey(names, args[index])
            symbol, type = names[args[index]]
            options[symbol] = parse(type, args[index + 1])
        else
            error("unknown option $(args[index])")
        end
        index += 2
    end
    return options
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 || error("usage: mesh_prism_tube_coupon.jl <case directory> output.msh [options]")
    isfile(abspath(ARGS[2])) && error("Refuse to overwrite $(ARGS[2])")
    census = generate_prism_tube_coupon(; case_directory=abspath(ARGS[1]), filename=abspath(ARGS[2]),
                                        parse_spike_options(ARGS)...)
    for (name, value) in sort(collect(census); by=first)
        name in ("Areas", "Tubes", "TubeFaces", "Sections") && continue
        println(name, ": ", value)
    end
    println("Areas: ", census["Areas"])
end
