# Source face velocity test: sensitivity w.r.t. a parameter that changes the 2D cross-section.
#
# Extends transmon_sensitivity_test.jl to cover the case where the design parameter
# changes the SOURCE FACE SHAPE (not just the extrusion distance). This is the bridge
# between ForwardDiff AD through DeviceLayout and the EGADS sensitivity pipeline.
#
# Geometry (cm, L0=1e-2):
#   Outer PEC box: 3 × 3 × 1.5 cm (Z-extruded)
#   Cavity: w × w × 1.3 cm centered hole (Z-extruded from z=0.1)
#   Boolean: box - cavity → re-entrant structure
#   Design parameter: w (cavity width, changes source face shape)
#
# The source face velocity for d/d(w) is:
#   - Left edge vertices: vx = -0.5, vy = 0
#   - Right edge vertices: vx = +0.5, vy = 0
#   - Bottom edge vertices: vx = 0, vy = -0.5
#   - Top edge vertices: vx = 0, vy = +0.5
#   All vertices move outward by ±0.5 in their respective directions.
#
# Usage: julia --project=path/to/DeviceLayout.jl-ad-sensitivity transmon_source_face_test.jl

using LinearAlgebra, Printf
using EGADS
using TetGen

include(joinpath(@__DIR__, "shared", "palace_io.jl"))

# ─── Inlined from DeviceLayoutEGADSExt (self-contained, no extension needed) ─

# ─── face_plane (from boolean_sensitivity_test.jl) ────────────────────────────

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

# ─── match_provenance (from boolean_sensitivity_test.jl) ─────────────────────

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

# ─── find_extrusion_base (from sensitivity_extract.jl) ──────────────────────

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

# ─── SourceFaceVelocity (from sensitivity_2d.jl) ────────────────────────────

struct SourceFaceVelocity
    vertices::Matrix{Float64}  # Nx2
    velocity::Matrix{Float64}  # Nx2
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
            dist_sq = (x - ax)^2 + (y - ay)^2
            s = 0.0
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

# ─── Geometry: Z-extruded box with Z-extruded cavity ─────────────────────────

"""Extrude a rectangle in XY at z=z0 by `dist` in +Z."""
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
    return body, src
end

const GEOM = (
    Lx = 3.0, Ly = 3.0, Lz = 1.5,  # outer box (cm)
    cavity_h = 1.3,                   # cavity height in Z (cm)
    cavity_z0 = 0.1,                  # cavity bottom offset from floor (cm)
)

function build_geometry(ctx; w=1.0)
    xc, yc = GEOM.Lx / 2, GEOM.Ly / 2

    # Outer PEC box
    box, _ = make_box_z(ctx, 0.0, 0.0, 0.0, GEOM.Lx, GEOM.Ly, GEOM.Lz)

    # Cavity: centered square w×w in XY, extruded in Z
    cavity_x0 = xc - w / 2
    cavity_y0 = yc - w / 2
    cavity, cavity_src = make_box_z(ctx, cavity_x0, cavity_y0, GEOM.cavity_z0,
                                     w, w, GEOM.cavity_h)

    # Boolean: result = box - cavity
    result_vec = EGADS.solid_boolean(box, cavity, EGADS.BoolOp.SUBTRACTION)
    @assert length(result_vec) >= 1 "Boolean produced no result"

    return (box=box, cavity=cavity, cavity_src=cavity_src, result=result_vec[1])
end

# ─── Source face velocity for d/d(w) ─────────────────────────────────────────

"""
Build SourceFaceVelocity for d/d(w) on the cavity source face (XY rectangle).

The source face is a w×w square centered at (xc, yc) in the XY plane.
For d/d(w), each corner moves outward:
  (xc-w/2, yc-w/2) → vel = (-0.5, -0.5)
  (xc+w/2, yc-w/2) → vel = (+0.5, -0.5)
  (xc+w/2, yc+w/2) → vel = (+0.5, +0.5)
  (xc-w/2, yc+w/2) → vel = (-0.5, +0.5)
"""
function make_source_face_velocity(w)
    xc, yc = GEOM.Lx / 2, GEOM.Ly / 2
    x0, y0 = xc - w / 2, yc - w / 2
    x1, y1 = xc + w / 2, yc + w / 2

    vertices = Float64[
        x0 y0
        x1 y0
        x1 y1
        x0 y1
    ]
    velocity = Float64[
        -0.5 -0.5
        +0.5 -0.5
        +0.5 +0.5
        -0.5 +0.5
    ]

    return SourceFaceVelocity(vertices, velocity)
end

# ─── Mesh + Palace ────────────────────────────────────────────────────────────

function mesh_geometry(body, mesh_file; max_edge=0.25, tetgen_opts="pq1.2a0.005")
    tess = EGADS.tessellate(body; max_edge=max_edge, sag=0.01, angle=15.0)
    status = EGADS.tess_status(tess)
    npts = status.npts
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

    mkpath(dirname(mesh_file))
    nvol = numberofpoints(output)
    ntet = numberoftetrahedra(output)
    ntri = numberoftrifaces(output)

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
        tfaces_out = output.trifacelist
        tmarkers_out = output.trifacemarkerlist
        println(io, "\nTriangles\n$ntri")
        for i in 1:ntri
            @printf(io, "%d %d %d %d\n", tfaces_out[1, i], tfaces_out[2, i], tfaces_out[3, i], tmarkers_out[i])
        end
        println(io, "\nEnd")
    end

    return (nvol=nvol, ntet=ntet, ntri=ntri, nfaces=length(faces))
end

function write_config(filename, mesh_path, nfaces; output_dir="postpro/sfv_test")
    config = """{
  "Problem": { "Type": "Eigenmode", "Verbose": 2, "Output": "$output_dir" },
  "Model": { "Mesh": "$mesh_path", "L0": 1.0e-2 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0 }] },
  "Boundaries": { "PEC": { "Attributes": $(collect(1:nfaces)) } },
  "Solver": {
    "Order": 1, "Device": "CPU",
    "Eigenmode": { "N": 2, "Tol": 1.0e-8, "Target": 2.0, "Sensitivity": true, "MeshSensitivity": true },
    "Linear": { "Type": "Default", "Tol": 1.0e-8, "MaxIts": 200 }
  }
}"""
    write(filename, config)
end

# ─── Chain contraction with source face velocity ─────────────────────────────

"""
Compute df/d(w) using provenance + analytical extrude formula with source face velocity.

For cavity faces (body_idx==2), the velocity at each point is:
  v(x,y,z) = v_source(x₀, y₀)
where (x₀, y₀) is the projection of (x,y,z) onto the source face (XY plane at z=cavity_z0).

Since dist_dot=0 and dir_dot=0 (cavity height and direction don't change with w),
the velocity is purely from the source face shape change.
"""
function compute_chain_sfv(geo, palace_mesh, dfdx, is_boundary; w=1.0)
    parent_bodies = [geo.box, geo.cavity]
    provenance = match_provenance(geo.result, parent_bodies)

    # Cavity extrusion parameters
    cavity_dir = [0.0, 0.0, 1.0]
    cavity_dist = GEOM.cavity_h

    # Find extrusion base
    cavity_faces = EGADS.get_body_topos(geo.cavity, EGADS.OClass.FACE)
    source_origin = find_extrusion_base(cavity_faces, cavity_dir)

    # Source face velocity for d/d(w)
    sfv = make_source_face_velocity(w)

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
                    # Cavity face: analytical extrude with source face velocity
                    α = clamp(dot(pt - source_origin, cavity_dir) / cavity_dist, 0.0, 1.0)
                    # dist_dot = 0, dir_dot = 0, so only source face velocity contributes
                    pos_source = pt .- (α * cavity_dist) .* cavity_dir
                    vel = interpolate_source_velocity(sfv, pos_source[1], pos_source[2])
                    dxdp_bdy[vi, :] .= vel
                end
                # Box faces (body_idx==1): zero velocity (box doesn't change)
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
    w = 1.0  # cavity width (cm)

    println("═══ Source Face Velocity Sensitivity Test ═══")
    println("  Geometry: 3×3×1.5 cm PEC box, $(w)×$(w)×1.3 cm cavity (Z-extruded)")
    println("  Parameter: w (cavity width, changes 2D source face shape)")
    println()

    palace = joinpath(@__DIR__, "..", "build", "palace-build", "palace-arm64.bin")
    @assert isfile(palace) "Palace binary not found: $palace"

    # ─── Step 1: Build, mesh, solve at baseline ──────────────────────────────

    println("Step 1: Baseline geometry + Palace...")
    ctx = EGADS.EGADSContext()
    EGADS.set_out_level!(ctx, 0)
    geo = build_geometry(ctx; w=w)

    mesh_file = "mesh/sfv_baseline.mesh"
    minfo = mesh_geometry(geo.result, mesh_file)

    config_file = "sfv_baseline.json"
    outdir = "postpro/sfv_baseline"
    write_config(config_file, mesh_file, minfo.nfaces; output_dir=outdir)
    mkpath(outdir)

    flush(stdout)
    t0 = time()
    run(pipeline(`$palace $config_file`; stdout=devnull, stderr=devnull))
    @printf("  Palace completed in %.1fs\n", time() - t0)

    freq = read_eigenfrequency(joinpath(outdir, "eig.csv"); mode=1)
    @printf("  Mesh: %d verts, %d tets, %d faces\n", minfo.nvol, minfo.ntet, minfo.nfaces)
    @printf("  Mode 1: f = %.6f GHz\n", freq)

    # Read mesh sensitivity
    palace_mesh = read_mesh(joinpath(outdir, "mesh_for_sensitivity.mesh"))
    sol = read_solb(joinpath(outdir, "mesh_sensitivity.solb"))
    dfdx = sol.data
    is_boundary = classify_boundary_vertices(palace_mesh)
    nv = size(palace_mesh.vertices, 1)
    @printf("  Palace mesh: %d verts, %d boundary\n", nv, count(is_boundary))

    # ─── Step 2: Chain contraction with source face velocity ─────────────────

    println("\nStep 2: Chain contraction (source face velocity for d/d(w))...")
    chain = compute_chain_sfv(geo, palace_mesh, dfdx, is_boundary; w=w)
    @printf("  Mapped: %d / %d boundary vertices\n", chain.n_mapped, count(is_boundary))
    @printf("  df/d(w) = %.6e GHz/cm\n", chain.dfdp)

    EGADS.close(ctx)

    # ─── Step 3: Mesh-perturbation Taylor test ───────────────────────────────

    println("\nStep 3: Mesh-perturbation Taylor test...")
    epsilons = [0.1, 0.05, 0.025, 0.0125]
    remainders = Float64[]

    println("  ε          f(p+ε)       |Δf|            |Δf-ε·df/dp|     rate")

    for (i, ε) in enumerate(epsilons)
        # Perturbed mesh
        pert_mesh = "mesh/sfv_pert_$(i).mesh"
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
            if ntet > 0
                println(io, "\nTetrahedra\n$ntet")
                for e in 1:ntet
                    @printf(io, "%d %d %d %d %d\n",
                            palace_mesh.tetrahedra[e, 1], palace_mesh.tetrahedra[e, 2],
                            palace_mesh.tetrahedra[e, 3], palace_mesh.tetrahedra[e, 4],
                            palace_mesh.tet_refs[e])
                end
            end
            ntri = size(palace_mesh.triangles, 1)
            if ntri > 0
                println(io, "\nTriangles\n$ntri")
                for e in 1:ntri
                    @printf(io, "%d %d %d %d\n",
                            palace_mesh.triangles[e, 1], palace_mesh.triangles[e, 2],
                            palace_mesh.triangles[e, 3], palace_mesh.tri_refs[e])
                end
            end
            println(io, "\nEnd")
        end

        pert_config = "sfv_pert_$(i).json"
        pert_outdir = "postpro/sfv_pert_$(i)"
        write_config(pert_config, pert_mesh, minfo.nfaces; output_dir=pert_outdir)
        mkpath(pert_outdir)

        run(pipeline(`$palace $pert_config`; stdout=devnull, stderr=devnull))
        f_pert = read_eigenfrequency(joinpath(pert_outdir, "eig.csv"); mode=1)

        delta_f = f_pert - freq
        rem = abs(delta_f - ε * chain.dfdp)
        push!(remainders, rem)

        rate = length(remainders) >= 2 ?
            log(remainders[end-1] / remainders[end]) /
            log(epsilons[length(remainders)-1] / ε) : NaN

        @printf("  %.4f   %10.6f   %12.6e   %12.6e   %.2f\n",
                ε, f_pert, abs(delta_f), rem, rate)

        # Clean up
        rm(pert_outdir; recursive=true, force=true)
        rm(pert_config; force=true)
        rm(pert_mesh; force=true)
    end

    # ─── Summary ──────────────────────────────────────────────────────────────

    rates = Float64[]
    for i in 2:length(remainders)
        push!(rates, log(remainders[i-1] / remainders[i]) / log(epsilons[i-1] / epsilons[i]))
    end

    println("\n═══ Summary ═══")
    @printf("  f₀ = %.6f GHz at w = %.2f cm\n", freq, w)
    @printf("  df/d(w) = %.6e GHz/cm\n", chain.dfdp)
    if !isempty(rates)
        @printf("  Taylor rates = %s\n", join([@sprintf("%.2f", r) for r in rates], ", "))
        max_rate = maximum(rates)
        pass = max_rate > 1.5
        @printf("  Status: %s (max rate %.2f, need > 1.5 for O(ε²))\n",
                pass ? "PASS ✓" : "FAIL ✗", max_rate)
    end

    # Clean up baseline
    rm(outdir; recursive=true, force=true)
    rm(config_file; force=true)
    rm(mesh_file; force=true)

    println("\nDone.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
