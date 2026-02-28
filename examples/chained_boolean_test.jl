# Chained boolean sensitivity test: substrate - cavity1 - cavity2
#
# Tests sensitivity through TWO sequential boolean subtractions.
# A real transmon has multiple subtractions (island, resonator, claw).
# This validates that provenance matching works through chained operations.
#
# Geometry (cm):
#   Outer box: 4 × 4 × 2 cm
#   Cavity 1: 1.5 × 1.5 × 1.5 cm centered, Z-extruded from z=0.25
#   Cavity 2: 0.5 × w2 × 1.5 cm offset, Z-extruded from z=0.25
#   Boolean: (outer - cavity1) - cavity2
#   Design parameter: w2 (width of cavity2 in Y)
#
# Usage: julia --project=path/to/DeviceLayout.jl-ad-sensitivity chained_boolean_test.jl

using LinearAlgebra, Printf
using EGADS
using TetGen

include(joinpath(@__DIR__, "shared", "palace_io.jl"))

# ─── Inlined helpers ────────────────────────────────────────────────────────

function face_plane(face)
    ri = EGADS.get_range(face)
    um = (ri.range[1] + ri.range[2]) / 2
    vm = (ri.range[3] + ri.range[4]) / 2
    res = EGADS.evaluate(face, [um, vm])
    n = cross(res.du, res.dv)
    nrm = norm(n)
    nrm < 1e-30 && return (normal=[0.0, 0.0, 0.0], offset=0.0, point=res.xyz)
    n = n / nrm
    d = dot(n, res.xyz)
    return (normal=n, offset=d, point=res.xyz)
end

function match_provenance(result_body, parent_bodies; tol=1e-6)
    result_faces = EGADS.get_body_topos(result_body, EGADS.OClass.FACE)
    parent_data = Tuple{Int, Int, EGADS.EGObject, NamedTuple}[]
    for (bi, body) in enumerate(parent_bodies)
        pfaces = EGADS.get_body_topos(body, EGADS.OClass.FACE)
        for (fi, pf) in enumerate(pfaces)
            push!(parent_data, (bi, fi, pf, face_plane(pf)))
        end
    end
    provenance = Dict{Int, NamedTuple{(:body_idx, :face_idx), Tuple{Int, Int}}}()
    for (ri, rf) in enumerate(result_faces)
        rp = face_plane(rf)
        for (bi, fi, _, pp) in parent_data
            ndot = abs(dot(rp.normal, pp.normal))
            if ndot > 1 - tol &&
               (abs(rp.offset - pp.offset) < tol || abs(rp.offset + pp.offset) < tol)
                provenance[ri] = (body_idx=bi, face_idx=fi)
                break
            end
        end
    end
    return provenance
end

function find_extrusion_base(faces, dir)
    min_proj = Inf
    base_pos = [0.0, 0.0, 0.0]
    for f in faces
        rng = EGADS.get_range(f)
        um = (rng.range[1] + rng.range[2]) / 2
        vm = (rng.range[3] + rng.range[4]) / 2
        center = EGADS.evaluate(f, [um, vm]).xyz
        proj = dot(center, dir)
        if proj < min_proj
            min_proj = proj
            base_pos = center
        end
    end
    return base_pos
end

struct SourceFaceVelocity
    vertices::Matrix{Float64}
    velocity::Matrix{Float64}
end

function interpolate_source_velocity(sfv::SourceFaceVelocity, x::Float64, y::Float64)
    nv = size(sfv.vertices, 1)
    nv == 0 && return [0.0, 0.0, 0.0]
    min_dist_sq = Inf
    best_s = 0.0
    best_i = 1
    for i in 1:nv
        j = mod1(i + 1, nv)
        ax, ay = sfv.vertices[i, 1], sfv.vertices[i, 2]
        bx, by = sfv.vertices[j, 1], sfv.vertices[j, 2]
        dx_ab, dy_ab = bx - ax, by - ay
        len_sq = dx_ab^2 + dy_ab^2
        if len_sq < 1e-30
            s = 0.0
            dist_sq = (x - ax)^2 + (y - ay)^2
        else
            s = clamp(((x - ax) * dx_ab + (y - ay) * dy_ab) / len_sq, 0.0, 1.0)
            px, py = ax + s * dx_ab, ay + s * dy_ab
            dist_sq = (x - px)^2 + (y - py)^2
        end
        if dist_sq < min_dist_sq
            min_dist_sq = dist_sq
            best_s = s
            best_i = i
        end
    end
    j = mod1(best_i + 1, nv)
    vx = (1 - best_s) * sfv.velocity[best_i, 1] + best_s * sfv.velocity[j, 1]
    vy = (1 - best_s) * sfv.velocity[best_i, 2] + best_s * sfv.velocity[j, 2]
    return [vx, vy, 0.0]
end

# ─── Geometry ────────────────────────────────────────────────────────────────

function make_box_z(ctx, x0, y0, z0, dx, dy, dist)
    n1 = EGADS.make_node(ctx, x0,      y0,      z0)
    n2 = EGADS.make_node(ctx, x0 + dx, y0,      z0)
    n3 = EGADS.make_node(ctx, x0 + dx, y0 + dy, z0)
    n4 = EGADS.make_node(ctx, x0,      y0 + dy, z0)
    e1 = EGADS.make_edge_from_nodes(ctx, n1, n2)
    e2 = EGADS.make_edge_from_nodes(ctx, n2, n3)
    e3 = EGADS.make_edge_from_nodes(ctx, n3, n4)
    e4 = EGADS.make_edge_from_nodes(ctx, n4, n1)
    loop = EGADS.make_wire_loop(ctx, [e1, e2, e3, e4])
    face = EGADS.make_plane_face(ctx, loop)
    src = EGADS.make_face_body(ctx, face)
    body = EGADS.extrude(src, dist, [0.0, 0.0, 1.0])
    return body
end

function build_chained_geometry(ctx; w2=1.0)
    # Outer box
    outer = make_box_z(ctx, 0.0, 0.0, 0.0, 4.0, 4.0, 2.0)

    # Cavity 1: centered 1.5×1.5×1.5
    c1 = make_box_z(ctx, 1.25, 1.25, 0.25, 1.5, 1.5, 1.5)

    # First boolean: outer - c1
    result1_vec = EGADS.solid_boolean(outer, c1, EGADS.BoolOp.SUBTRACTION)
    @assert length(result1_vec) >= 1
    result1 = result1_vec[1]

    # Cavity 2: offset, width w2 in Y
    c2_x0 = 3.0  # right side
    c2_y0 = 2.0 - w2 / 2  # centered in Y
    c2 = make_box_z(ctx, c2_x0, c2_y0, 0.25, 0.5, w2, 1.5)

    # Second boolean: result1 - c2
    result2_vec = EGADS.solid_boolean(result1, c2, EGADS.BoolOp.SUBTRACTION)
    @assert length(result2_vec) >= 1

    return (outer=outer, c1=c1, c2=c2, result1=result1, result=result2_vec[1])
end

# ─── Mesh + Palace ────────────────────────────────────────────────────────────

function mesh_and_write(body, mesh_file; max_edge=0.35, tetgen_opts="pq1.2a0.02")
    tess = EGADS.tessellate(body; max_edge=max_edge, sag=0.01, angle=15.0)
    npts = EGADS.tess_status(tess).npts
    faces = EGADS.get_body_topos(body, EGADS.OClass.FACE)

    pts = zeros(Float64, 3, npts)
    for fi in eachindex(faces)
        fm = EGADS.get_face_mesh(tess, body, fi)
        for j in 1:size(fm.xyz, 1)
            ig = EGADS.eg_local_to_global(tess.ptr, fi, j)
            pts[:, ig] .= fm.xyz[j, :]
        end
    end

    tris = Int32[]
    tri_markers = Int32[]
    for fi in eachindex(faces)
        fm = EGADS.get_face_mesh(tess, body, fi)
        for t in 1:size(fm.tris, 1)
            for k in 1:3
                push!(tris, Int32(EGADS.eg_local_to_global(tess.ptr, fi, fm.tris[t, k])))
            end
            push!(tri_markers, Int32(fi))
        end
    end

    input = RawTetGenIO{Float64}()
    input.pointlist = pts
    facetlist!(input, reshape(tris, 3, length(tri_markers)))
    input.facetmarkerlist = tri_markers
    output = tetrahedralize(input, tetgen_opts)

    nvol = numberofpoints(output)
    ntet = numberoftetrahedra(output)
    ntri = numberoftrifaces(output)

    mkpath(dirname(mesh_file))
    open(mesh_file, "w") do io
        println(io, "MeshVersionFormatted 2\n\nDimension 3\n")
        pts_out = output.pointlist
        println(io, "Vertices\n$nvol")
        for i in 1:nvol
            @printf(io, "%.15e %.15e %.15e 0\n", pts_out[1, i], pts_out[2, i], pts_out[3, i])
        end
        tets = output.tetrahedronlist
        println(io, "\nTetrahedra\n$ntet")
        for i in 1:ntet
            @printf(io, "%d %d %d %d 1\n", tets[1, i], tets[2, i], tets[3, i], tets[4, i])
        end
        tfo = output.trifacelist
        tmo = output.trifacemarkerlist
        println(io, "\nTriangles\n$ntri")
        for i in 1:ntri
            @printf(io, "%d %d %d %d\n", tfo[1, i], tfo[2, i], tfo[3, i], tmo[i])
        end
        println(io, "\nEnd")
    end
    return (nvol=nvol, ntet=ntet, nfaces=length(faces))
end

function write_config(filename, mesh_path, nfaces; outdir="postpro/chained")
    config = """{
  "Problem": { "Type": "Eigenmode", "Verbose": 2, "Output": "$outdir" },
  "Model": { "Mesh": "$mesh_path", "L0": 1.0e-2 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0 }] },
  "Boundaries": { "PEC": { "Attributes": $(collect(1:nfaces)) } },
  "Solver": {
    "Order": 1, "Device": "CPU",
    "Eigenmode": { "N": 2, "Tol": 1.0e-8, "Target": 1.5, "Sensitivity": true, "MeshSensitivity": true },
    "Linear": { "Type": "Default", "Tol": 1.0e-8, "MaxIts": 200 }
  }
}"""
    write(filename, config)
end

# ─── Chain contraction for d/d(w2) ────────────────────────────────────────────

function compute_chain(geo, palace_mesh, dfdx, is_boundary; w2=1.0)
    # For the chained boolean, provenance matching against the SECOND boolean's parents:
    # result2 = result1 - c2. So parent bodies are [result1, c2].
    parent_bodies = [geo.result1, geo.c2]
    provenance = match_provenance(geo.result, parent_bodies)

    # Cavity 2 extrusion parameters
    dir = [0.0, 0.0, 1.0]
    dist = 1.5
    c2_faces = EGADS.get_body_topos(geo.c2, EGADS.OClass.FACE)
    source_origin = find_extrusion_base(c2_faces, dir)

    # Source face velocity for d/d(w2): centered square expanding in Y
    c2_x0, c2_y0 = 3.0, 2.0 - w2 / 2
    sfv = SourceFaceVelocity(
        Float64[c2_x0 c2_y0; c2_x0+0.5 c2_y0; c2_x0+0.5 c2_y0+w2; c2_x0 c2_y0+w2],
        Float64[0.0 -0.5; 0.0 -0.5; 0.0 0.5; 0.0 0.5],
    )

    result_faces = EGADS.get_body_topos(geo.result, EGADS.OClass.FACE)
    nv = size(palace_mesh.vertices, 1)
    dxdp_bdy = zeros(Float64, nv, 3)
    n_mapped = 0

    for vi in findall(is_boundary)
        pt = palace_mesh.vertices[vi, :]
        for fi in eachindex(result_faces)
            inv_res = EGADS.inv_evaluate(result_faces[fi], pt)
            if norm(inv_res.results[1:3] - pt) < 1e-6
                if haskey(provenance, fi) && provenance[fi].body_idx == 2
                    # Cavity 2 face: source face velocity
                    α = clamp(dot(pt - source_origin, dir) / dist, 0.0, 1.0)
                    pos_source = pt .- (α * dist) .* dir
                    vel = interpolate_source_velocity(sfv, pos_source[1], pos_source[2])
                    dxdp_bdy[vi, :] .= vel
                end
                # result1 faces (body_idx==1): zero velocity (only c2 changes)
                n_mapped += 1
                break
            end
        end
    end

    dxdp_full = propagate_boundary_sensitivity(palace_mesh, dxdp_bdy, is_boundary)
    dfdp = dot(vec(dfdx), vec(dxdp_full))
    return (dfdp=dfdp, n_mapped=n_mapped, dxdp_full=dxdp_full)
end

# ─── Main ────────────────────────────────────────────────────────────────────

function main()
    cd(@__DIR__)
    w2 = 1.0
    palace = joinpath(@__DIR__, "..", "build", "palace-build", "palace-arm64.bin")
    @assert isfile(palace)

    println("═══ Chained Boolean Sensitivity Test ═══")
    println("  Geometry: 4×4×2 box - 1.5³ cavity1 - 0.5×$(w2)×1.5 cavity2")
    println("  Parameter: w2 (cavity2 Y-width, source face shape)")
    println()

    ctx = EGADS.EGADSContext()
    EGADS.set_out_level!(ctx, 0)
    geo = build_chained_geometry(ctx; w2=w2)

    nfaces_result = length(EGADS.get_body_topos(geo.result, EGADS.OClass.FACE))
    println("  Result faces: $nfaces_result")

    mesh_file = "mesh/chained.mesh"
    minfo = mesh_and_write(geo.result, mesh_file)
    @printf("  Mesh: %d verts, %d tets, %d faces\n", minfo.nvol, minfo.ntet, minfo.nfaces)

    config_file = "chained.json"
    outdir = "postpro/chained"
    write_config(config_file, mesh_file, minfo.nfaces; outdir=outdir)
    mkpath(outdir)

    println("  Running Palace...")
    flush(stdout)
    run(pipeline(`$palace $config_file`; stdout=devnull, stderr=devnull))

    freq = read_eigenfrequency(joinpath(outdir, "eig.csv"); mode=1)
    @printf("  f₀ = %.6f GHz\n", freq)

    palace_mesh = read_mesh(joinpath(outdir, "mesh_for_sensitivity.mesh"))
    sol = read_solb(joinpath(outdir, "mesh_sensitivity.solb"))
    dfdx = sol.data
    is_boundary = classify_boundary_vertices(palace_mesh)
    nv = size(palace_mesh.vertices, 1)

    println("\n  Chain contraction...")
    chain = compute_chain(geo, palace_mesh, dfdx, is_boundary; w2=w2)
    @printf("  Mapped: %d / %d boundary\n", chain.n_mapped, count(is_boundary))
    @printf("  df/d(w2) = %.6e GHz/cm\n", chain.dfdp)

    EGADS.close(ctx)

    # Taylor test
    println("\n  Mesh-perturbation Taylor test:")
    epsilons = [0.1, 0.05, 0.025, 0.0125]
    remainders = Float64[]

    for (i, ε) in enumerate(epsilons)
        pert_mesh = "mesh/chained_pert_$(i).mesh"
        open(pert_mesh, "w") do io
            println(io, "MeshVersionFormatted 2\n\nDimension 3\n")
            println(io, "Vertices\n$nv")
            for vi in 1:nv
                @printf(io, "%.15e %.15e %.15e 0\n",
                        palace_mesh.vertices[vi, 1] + ε * chain.dxdp_full[vi, 1],
                        palace_mesh.vertices[vi, 2] + ε * chain.dxdp_full[vi, 2],
                        palace_mesh.vertices[vi, 3] + ε * chain.dxdp_full[vi, 3])
            end
            ntet = size(palace_mesh.tetrahedra, 1)
            println(io, "\nTetrahedra\n$ntet")
            for e in 1:ntet
                @printf(io, "%d %d %d %d %d\n",
                        palace_mesh.tetrahedra[e, 1], palace_mesh.tetrahedra[e, 2],
                        palace_mesh.tetrahedra[e, 3], palace_mesh.tetrahedra[e, 4],
                        palace_mesh.tet_refs[e])
            end
            ntri = size(palace_mesh.triangles, 1)
            println(io, "\nTriangles\n$ntri")
            for e in 1:ntri
                @printf(io, "%d %d %d %d\n",
                        palace_mesh.triangles[e, 1], palace_mesh.triangles[e, 2],
                        palace_mesh.triangles[e, 3], palace_mesh.tri_refs[e])
            end
            println(io, "\nEnd")
        end

        pert_config = "chained_pert_$(i).json"
        pert_outdir = "postpro/chained_pert_$(i)"
        write_config(pert_config, pert_mesh, minfo.nfaces; outdir=pert_outdir)
        mkpath(pert_outdir)
        run(pipeline(`$palace $pert_config`; stdout=devnull, stderr=devnull))
        f_pert = read_eigenfrequency(joinpath(pert_outdir, "eig.csv"); mode=1)

        rem = abs(f_pert - freq - ε * chain.dfdp)
        push!(remainders, rem)
        rate = length(remainders) >= 2 ?
            log(remainders[end-1] / remainders[end]) / log(epsilons[length(remainders)-1] / ε) : NaN
        @printf("  ε=%.4f: f=%.8f, rem=%.2e, rate=%.2f\n", ε, f_pert, rem, rate)

        rm(pert_outdir; recursive=true, force=true)
        rm(pert_config; force=true)
        rm(pert_mesh; force=true)
    end

    rates = [log(remainders[i-1] / remainders[i]) / log(epsilons[i-1] / epsilons[i])
             for i in 2:length(remainders)]
    max_rate = isempty(rates) ? 0.0 : maximum(rates)
    pass = max_rate > 1.5
    @printf("\n  Rates: %s\n", join([@sprintf("%.2f", r) for r in rates], ", "))
    @printf("  Status: %s (max rate %.2f)\n", pass ? "PASS ✓" : "FAIL ✗", max_rate)

    rm(outdir; recursive=true, force=true)
    rm(config_file; force=true)
    rm(mesh_file; force=true)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
