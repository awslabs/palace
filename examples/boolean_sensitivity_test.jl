# Boolean sensitivity test: box with a pocket subtracted.
#
# Step 1: Entity provenance — identify which post-boolean faces came from which pre-boolean body
# Step 2: Face velocities from parent bodies
# Step 3: Edge velocities via implicit function theorem at intersections
# Step 4: Full pipeline: EGADS boolean → TetGen → Palace → chain contraction → Taylor test
#
# Usage: julia --project=path/to/DeviceLayout.jl-egads-dev boolean_sensitivity_test.jl

using LinearAlgebra, Printf, SparseArrays, DelimitedFiles
using EGADS
using TetGen

# Shared Palace mesh/solb readers and boundary utilities
include(joinpath(@__DIR__, "shared", "palace_io.jl"))

# ─── Face plane extraction ────────────────────────────────────────────────
# For a planar face, returns (normal, offset) defining the plane n·x = d.
# Uses du × dv from EGADS evaluate at the face midpoint.

function face_plane(face)
    ri = EGADS.get_range(face)
    um = (ri.range[1] + ri.range[2]) / 2
    vm = (ri.range[3] + ri.range[4]) / 2
    res = EGADS.evaluate(face, [um, vm])
    n = cross(res.du, res.dv)
    n = n / norm(n)
    d = dot(n, res.xyz)
    return (normal=n, offset=d, point=res.xyz)
end

# ─── Entity provenance matching ───────────────────────────────────────────
# Match each post-boolean face to a pre-boolean face from one of the source bodies.
# Returns Dict{Int, (body_idx, face_idx)} where body_idx ∈ {1=src, 2=tool}.

function match_provenance(result_body, src_body, tool_body; tol=1e-6)
    result_faces = EGADS.get_body_topos(result_body, EGADS.OClass.FACE)
    src_faces = EGADS.get_body_topos(src_body, EGADS.OClass.FACE)
    tool_faces = EGADS.get_body_topos(tool_body, EGADS.OClass.FACE)

    # Compute plane equations for all pre-boolean faces
    src_planes = [face_plane(f) for f in src_faces]
    tool_planes = [face_plane(f) for f in tool_faces]

    provenance = Dict{Int, NamedTuple{(:body_idx, :face_idx, :body_name), Tuple{Int, Int, String}}}()

    for (ri, rf) in enumerate(result_faces)
        rp = face_plane(rf)

        # Try matching against src faces
        matched = false
        for (si, sp) in enumerate(src_planes)
            # Planes match if normals are parallel and offsets agree
            ndot = abs(dot(rp.normal, sp.normal))
            if ndot > 1 - tol && abs(rp.offset - sp.offset) < tol
                provenance[ri] = (body_idx=1, face_idx=si, body_name="src")
                matched = true
                break
            end
            # Also check anti-parallel normals (face orientation may flip)
            if ndot > 1 - tol && abs(rp.offset + sp.offset) < tol
                provenance[ri] = (body_idx=1, face_idx=si, body_name="src")
                matched = true
                break
            end
        end
        matched && continue

        # Try matching against tool faces
        for (ti, tp) in enumerate(tool_planes)
            ndot = abs(dot(rp.normal, tp.normal))
            if ndot > 1 - tol && abs(rp.offset - tp.offset) < tol
                provenance[ri] = (body_idx=2, face_idx=ti, body_name="tool")
                matched = true
                break
            end
            if ndot > 1 - tol && abs(rp.offset + tp.offset) < tol
                provenance[ri] = (body_idx=2, face_idx=ti, body_name="tool")
                matched = true
                break
            end
        end

        if !matched
            @warn "Face $ri: no provenance match" normal=rp.normal offset=rp.offset
        end
    end

    return provenance
end

# ─── Geometry builder ─────────────────────────────────────────────────────
# Build outer box + pocket at given parameter, optionally set _dot.

function build_boolean_geometry(ctx; pocket_dx=2.0, set_dot=false, dx_dot=0.0)
    outer_data = Float64[0.0, 0.0, 0.0, 4.0, 4.0, 4.0]
    pocket_data = Float64[1.0, 1.0, 3.0, pocket_dx, 2.0, 1.0]

    outer = EGADS.make_box(ctx, outer_data...)
    pocket = EGADS.make_box(ctx, pocket_data...)

    if set_dot
        EGADS.make_solid_body_dot!(outer, EGADS.SolidType.BOX, outer_data, zeros(6))
        pocket_dot = Float64[0.0, 0.0, 0.0, dx_dot, 0.0, 0.0]
        EGADS.make_solid_body_dot!(pocket, EGADS.SolidType.BOX, pocket_data, pocket_dot)
    end

    result_bodies = EGADS.solid_boolean(outer, pocket, EGADS.BoolOp.SUBTRACTION)
    @assert length(result_bodies) == 1 "Expected 1 result body, got $(length(result_bodies))"

    return (outer=outer, pocket=pocket, result=result_bodies[1],
            outer_data=outer_data, pocket_data=pocket_data)
end

# ─── Test 1: Provenance matching ─────────────────────────────────────────

function test_boolean_provenance()
    println("═══ Test 1: Provenance matching ═══")
    ctx = EGADS.EGADSContext()

    geo = build_boolean_geometry(ctx)
    result_faces = EGADS.get_body_topos(geo.result, EGADS.OClass.FACE)
    println("Result: $(length(result_faces)) faces")

    prov = match_provenance(geo.result, geo.outer, geo.pocket)
    println("Matched $(length(prov)) / $(length(result_faces)) faces:")
    for ri in sort(collect(keys(prov)))
        p = prov[ri]
        rp = face_plane(result_faces[ri])
        @printf("  Result face %2d → %4s face %d  (n=(%.2f,%.2f,%.2f), d=%.2f)\n",
                ri, p.body_name, p.face_idx, rp.normal..., rp.offset)
    end

    unmatched = setdiff(1:length(result_faces), keys(prov))
    pass = isempty(unmatched)
    println(pass ? "  PASS ✓" : "  FAIL ✗ — unmatched: $unmatched")

    EGADS.close(ctx)
    return pass
end

# ─── Test 2: Face velocity via provenance ─────────────────────────────────
# Set _dot on pocket body (d/d(pocket_dx)), compute face velocities on the
# post-boolean body by projecting to parent faces via provenance.
# Verify against finite differences at face midpoints.

function test_boolean_face_velocity()
    println("\n═══ Test 2: Face velocities via provenance + evaluate_dot ═══")

    pocket_dx = 2.0
    h = 1e-6  # FD step

    # ── Analytical: build at p, set _dot, evaluate_dot on parents ──
    ctx = EGADS.EGADSContext()
    geo = build_boolean_geometry(ctx; pocket_dx=pocket_dx, set_dot=true, dx_dot=1.0)

    result_faces = EGADS.get_body_topos(geo.result, EGADS.OClass.FACE)
    outer_faces = EGADS.get_body_topos(geo.outer, EGADS.OClass.FACE)
    pocket_faces = EGADS.get_body_topos(geo.pocket, EGADS.OClass.FACE)
    parent_faces = [outer_faces, pocket_faces]

    prov = match_provenance(geo.result, geo.outer, geo.pocket)

    # Compute analytical velocity at each result face midpoint
    analytical = Dict{Int, Vector{Float64}}()
    midpoints = Dict{Int, Vector{Float64}}()

    for (ri, rf) in enumerate(result_faces)
        rng = EGADS.get_range(rf)
        uv_mid = [(rng.range[1] + rng.range[2]) / 2, (rng.range[3] + rng.range[4]) / 2]
        xyz = EGADS.evaluate(rf, uv_mid).xyz
        midpoints[ri] = xyz

        p = prov[ri]
        pf = parent_faces[p.body_idx][p.face_idx]

        # Project to parent face
        inv = EGADS.inv_evaluate(pf, xyz)
        dist = norm(inv.results[1:3] - xyz)
        @assert dist < 1e-8 "inv_evaluate projection failed: dist=$dist on face $ri"

        res = EGADS.evaluate_dot(pf, inv.param, zeros(2))
        analytical[ri] = res.result_dot.xyz
    end
    EGADS.close(ctx)

    # ── FD: build at p+h, evaluate at same UV on corresponding parent face ──
    ctx2 = EGADS.EGADSContext()
    geo_h = build_boolean_geometry(ctx2; pocket_dx=pocket_dx + h)
    outer_faces_h = EGADS.get_body_topos(geo_h.outer, EGADS.OClass.FACE)
    pocket_faces_h = EGADS.get_body_topos(geo_h.pocket, EGADS.OClass.FACE)
    parent_faces_h = [outer_faces_h, pocket_faces_h]

    # For FD, we evaluate at the same UV on the parent face of the perturbed geometry.
    # The parent face at p+h has the same surface type, so same UV → different xyz.
    fd_vel = Dict{Int, Vector{Float64}}()

    for (ri, rf) in enumerate(result_faces)
        p = prov[ri]
        pf = parent_faces[p.body_idx][p.face_idx]  # parent at p (already freed!)
        pf_h = parent_faces_h[p.body_idx][p.face_idx]  # parent at p+h

        # Re-project the base point to the perturbed parent to get consistent UV
        # (For planes the UV is nearly identical, but be safe)
        xyz_base = midpoints[ri]
        inv_h = EGADS.inv_evaluate(pf_h, xyz_base)
        xyz_h = EGADS.evaluate(pf_h, inv_h.param).xyz

        fd_vel[ri] = (xyz_h - xyz_base) / h
    end
    EGADS.close(ctx2)

    # ── Compare ──
    println("Face midpoint velocities (d/d(pocket_dx)):")
    max_err = 0.0
    for ri in sort(collect(keys(analytical)))
        a = analytical[ri]
        fd = fd_vel[ri]
        err = norm(a - fd)
        max_err = max(max_err, err)
        body = prov[ri].body_name
        @printf("  Face %2d (%4s): analytical=(%.4f, %.4f, %.4f)  FD=(%.4f, %.4f, %.4f)  err=%.2e\n",
                ri, body, a..., fd..., err)
    end

    println("\n  Max error: $(max_err)")
    pass = max_err < 1e-4
    println("  Status: $(pass ? "PASS ✓" : "FAIL ✗")")
    return pass
end

# ─── Test 3: Edge velocity via implicit function theorem ──────────────────
# For edges at the intersection of faces from different parent bodies,
# solve the 3×3 system: [n₁; n₂; t] · v = [n₁·v₁; n₂·v₂; 0]

function edge_velocity_ift(face1, face2, edge_point, vel_face1, vel_face2)
    # Get surface normals at the edge point
    inv1 = EGADS.inv_evaluate(face1, edge_point)
    inv2 = EGADS.inv_evaluate(face2, edge_point)
    res1 = EGADS.evaluate(face1, inv1.param)
    res2 = EGADS.evaluate(face2, inv2.param)

    n1 = normalize(cross(res1.du, res1.dv))
    n2 = normalize(cross(res2.du, res2.dv))

    # Edge tangent: intersection of two planes
    t = normalize(cross(n1, n2))

    # 3×3 system: [n₁; n₂; t] · v = [n₁·v₁; n₂·v₂; 0]
    A = [n1'; n2'; t']
    b = [dot(n1, vel_face1); dot(n2, vel_face2); 0.0]
    return A \ b
end

function test_boolean_edge_velocity()
    println("\n═══ Test 3: Edge velocity via IFT ═══")

    pocket_dx = 2.0
    h = 1e-6

    # Build geometry with _dot
    ctx = EGADS.EGADSContext()
    geo = build_boolean_geometry(ctx; pocket_dx=pocket_dx, set_dot=true, dx_dot=1.0)

    result_faces = EGADS.get_body_topos(geo.result, EGADS.OClass.FACE)
    result_edges = EGADS.get_body_topos(geo.result, EGADS.OClass.EDGE)
    outer_faces = EGADS.get_body_topos(geo.outer, EGADS.OClass.FACE)
    pocket_faces = EGADS.get_body_topos(geo.pocket, EGADS.OClass.FACE)
    parent_faces_list = [outer_faces, pocket_faces]

    prov = match_provenance(geo.result, geo.outer, geo.pocket)

    println("Result edges: $(length(result_edges))")

    # Find intersection edges: edges bordered by faces from different parent bodies.
    # Use get_body_topos(body, edge, FACE) to get the 2 adjacent faces per edge.
    intersection_edges = Int[]
    edge_parent_faces = Dict{Int, Tuple{Int, Int}}()  # edge_idx → (face_idx1, face_idx2)

    for (ei, edge) in enumerate(result_edges)
        adj_faces = EGADS.get_body_topos(geo.result, edge, EGADS.OClass.FACE)
        if length(adj_faces) == 2
            # Map adjacent faces to result face indices
            fi1 = findfirst(f -> EGADS.is_same(f, adj_faces[1]), result_faces)
            fi2 = findfirst(f -> EGADS.is_same(f, adj_faces[2]), result_faces)
            if !isnothing(fi1) && !isnothing(fi2) &&
               haskey(prov, fi1) && haskey(prov, fi2) &&
               prov[fi1].body_idx != prov[fi2].body_idx
                push!(intersection_edges, ei)
                edge_parent_faces[ei] = (fi1, fi2)
            end
        end
    end

    println("Intersection edges: $(length(intersection_edges))")

    # Compute edge midpoint velocities via IFT and compare to FD
    max_err = 0.0
    for ei in intersection_edges
        edge = result_edges[ei]
        rng = EGADS.get_range(edge)
        t_mid = (rng.range[1] + rng.range[2]) / 2
        res = EGADS.evaluate(edge, [t_mid])
        edge_pt = res.xyz

        fi1, fi2 = edge_parent_faces[ei]
        p1, p2 = prov[fi1], prov[fi2]
        pf1 = parent_faces_list[p1.body_idx][p1.face_idx]
        pf2 = parent_faces_list[p2.body_idx][p2.face_idx]

        # Face velocities at the edge point (from parent evaluate_dot)
        inv1 = EGADS.inv_evaluate(pf1, edge_pt)
        inv2 = EGADS.inv_evaluate(pf2, edge_pt)
        vel1 = EGADS.evaluate_dot(pf1, inv1.param, zeros(2)).result_dot.xyz
        vel2 = EGADS.evaluate_dot(pf2, inv2.param, zeros(2)).result_dot.xyz

        # IFT edge velocity
        v_ift = edge_velocity_ift(pf1, pf2, edge_pt, vel1, vel2)

        @printf("  Edge %2d at (%.2f,%.2f,%.2f): v_IFT=(%.4f,%.4f,%.4f)  faces=%d(%s)+%d(%s)\n",
                ei, edge_pt..., v_ift..., fi1, p1.body_name, fi2, p2.body_name)
    end

    # FD validation: project base physical point onto perturbed edge.
    # Using parametric midpoint is wrong because edge extent changes.
    EGADS.close(ctx)
    ctx2 = EGADS.EGADSContext()
    geo_h = build_boolean_geometry(ctx2; pocket_dx=pocket_dx + h)
    result_edges_h = EGADS.get_body_topos(geo_h.result, EGADS.OClass.EDGE)
    outer_faces_h = EGADS.get_body_topos(geo_h.outer, EGADS.OClass.FACE)
    pocket_faces_h = EGADS.get_body_topos(geo_h.pocket, EGADS.OClass.FACE)
    parent_faces_h = [outer_faces_h, pocket_faces_h]

    # Re-build at p for the base positions
    ctx3 = EGADS.EGADSContext()
    geo_base = build_boolean_geometry(ctx3; pocket_dx=pocket_dx, set_dot=true, dx_dot=1.0)
    result_edges_base = EGADS.get_body_topos(geo_base.result, EGADS.OClass.EDGE)
    outer_faces_base = EGADS.get_body_topos(geo_base.outer, EGADS.OClass.FACE)
    pocket_faces_base = EGADS.get_body_topos(geo_base.pocket, EGADS.OClass.FACE)
    parent_faces_base = [outer_faces_base, pocket_faces_base]
    prov_base = match_provenance(geo_base.result, geo_base.outer, geo_base.pocket)

    println("\nFD verification (project base point to perturbed parent faces):")
    for ei in intersection_edges
        edge = result_edges_base[ei]
        rng = EGADS.get_range(edge)
        t_mid = (rng.range[1] + rng.range[2]) / 2
        xyz_base = EGADS.evaluate(edge, [t_mid]).xyz

        # IFT at base
        fi1, fi2 = edge_parent_faces[ei]
        p1, p2 = prov_base[fi1], prov_base[fi2]
        pf1 = parent_faces_base[p1.body_idx][p1.face_idx]
        pf2 = parent_faces_base[p2.body_idx][p2.face_idx]
        inv1 = EGADS.inv_evaluate(pf1, xyz_base)
        inv2 = EGADS.inv_evaluate(pf2, xyz_base)
        vel1 = EGADS.evaluate_dot(pf1, inv1.param, zeros(2)).result_dot.xyz
        vel2 = EGADS.evaluate_dot(pf2, inv2.param, zeros(2)).result_dot.xyz
        v_ift = edge_velocity_ift(pf1, pf2, xyz_base, vel1, vel2)

        # FD: evaluate perturbed parent faces at base UV → xyz(p+h)
        pf1_h = parent_faces_h[p1.body_idx][p1.face_idx]
        pf2_h = parent_faces_h[p2.body_idx][p2.face_idx]
        xyz1_h = EGADS.evaluate(pf1_h, inv1.param).xyz
        xyz2_h = EGADS.evaluate(pf2_h, inv2.param).xyz
        # Both surfaces should give same point at intersection; use average
        # Actually, use IFT approach: find intersection of perturbed surfaces near base
        # For planes, simpler: project base point onto each perturbed face, then
        # the IFT point is where both constraints are satisfied.
        # Simple FD: average the two face FD velocities projected via IFT
        fd1 = (xyz1_h - xyz_base) / h  # velocity according to face 1
        fd2 = (xyz2_h - xyz_base) / h  # velocity according to face 2
        # Use IFT on the FD velocities
        fd_ift = edge_velocity_ift(pf1, pf2, xyz_base, fd1, fd2)

        err = norm(v_ift - fd_ift)
        max_err = max(max_err, err)
        @printf("  Edge %2d: IFT=(%.4f,%.4f,%.4f)  FD_IFT=(%.4f,%.4f,%.4f)  err=%.2e\n",
                ei, v_ift..., fd_ift..., err)
    end

    EGADS.close(ctx2)
    EGADS.close(ctx3)

    println("\n  Max error: $(max_err)")
    pass = max_err < 1e-3
    println("  Status: $(pass ? "PASS ✓" : "FAIL ✗")")
    return pass
end

# ─── Boolean velocity field computation ───────────────────────────────────
# Compute dx/dp at every tessellation surface node of a boolean result body.
# Handles face-interior nodes, same-body edge nodes, and intersection edge nodes (IFT).

function compute_boolean_velocities(geo, tess, prov, parent_faces)
    result_body = geo.result
    result_faces = EGADS.get_body_topos(result_body, EGADS.OClass.FACE)
    result_edges = EGADS.get_body_topos(result_body, EGADS.OClass.EDGE)
    status = EGADS.tess_status(tess)
    npts = status.npts

    # ── Step 1: Build global coordinate array and node→faces map ──
    global_xyz = zeros(Float64, 3, npts)
    node_faces = Dict{Int, Vector{Int}}()  # global_id → [face_ids that contain this node]

    for (fi, face) in enumerate(result_faces)
        fm = EGADS.get_face_mesh(tess, result_body, fi)
        for j in 1:size(fm.xyz, 1)
            ig = EGADS.eg_local_to_global(tess.ptr, fi, j)
            global_xyz[:, ig] .= fm.xyz[j, :]
            if !haskey(node_faces, ig)
                node_faces[ig] = Int[]
            end
            if !(fi in node_faces[ig])
                push!(node_faces[ig], fi)
            end
        end
    end

    # ── Step 2: Find intersection edges ──
    intersection_edge_set = Set{Int}()  # result edge indices
    edge_adj_faces = Dict{Int, Tuple{Int, Int}}()

    for (ei, edge) in enumerate(result_edges)
        adj = EGADS.get_body_topos(result_body, edge, EGADS.OClass.FACE)
        if length(adj) == 2
            fi1 = findfirst(f -> EGADS.is_same(f, adj[1]), result_faces)
            fi2 = findfirst(f -> EGADS.is_same(f, adj[2]), result_faces)
            if !isnothing(fi1) && !isnothing(fi2) &&
               haskey(prov, fi1) && haskey(prov, fi2) &&
               prov[fi1].body_idx != prov[fi2].body_idx
                push!(intersection_edge_set, ei)
                edge_adj_faces[ei] = (fi1, fi2)
            end
        end
    end

    # ── Step 3: Compute velocity for each global node ──
    # For face-interior nodes (1 face): use evaluate_dot on parent face.
    # For shared nodes (2+ faces): gather surface normal constraints from ALL
    # adjacent faces and solve the (possibly overdetermined) linear system.
    # This handles same-body edges (e.g., pocket wall + pocket floor) and
    # cross-body edges (e.g., pocket wall + outer top) uniformly.
    xyz_dot = zeros(Float64, 3 * npts)
    n_face_vel = 0
    n_shared_vel = 0
    n_zero = 0

    for ig in 1:npts
        !haskey(node_faces, ig) && continue
        pt = global_xyz[:, ig]
        faces_list = node_faces[ig]

        if length(faces_list) == 1
            # Face-interior node: unambiguous, use evaluate_dot
            fi = faces_list[1]
            !haskey(prov, fi) && continue
            p = prov[fi]
            pf = parent_faces[p.body_idx][p.face_idx]
            inv = EGADS.inv_evaluate(pf, pt)
            vel = EGADS.evaluate_dot(pf, inv.param, zeros(2)).result_dot.xyz
        else
            # Shared node (edge or vertex): gather constraints from ALL adjacent faces.
            # Each face gives one constraint: nᵢ · v = nᵢ · vᵢ
            normals = Vector{Float64}[]
            rhs_vals = Float64[]
            seen_normals = Set{Tuple{Float64,Float64,Float64}}()

            for fi in faces_list
                !haskey(prov, fi) && continue
                p = prov[fi]
                pf = parent_faces[p.body_idx][p.face_idx]
                inv = EGADS.inv_evaluate(pf, pt)
                res = EGADS.evaluate(pf, inv.param)
                n = normalize(cross(res.du, res.dv))
                vel_face = EGADS.evaluate_dot(pf, inv.param, zeros(2)).result_dot.xyz

                # Deduplicate by normal direction (parallel or anti-parallel)
                n_key = (round(n[1]; digits=8), round(n[2]; digits=8), round(n[3]; digits=8))
                n_key_neg = (round(-n[1]; digits=8), round(-n[2]; digits=8), round(-n[3]; digits=8))
                if n_key in seen_normals || n_key_neg in seen_normals
                    continue
                end
                push!(seen_normals, n_key)
                push!(normals, n)
                push!(rhs_vals, dot(n, vel_face))
            end

            if length(normals) >= 3
                # Vertex: 3+ independent constraints → solve exactly or least-squares
                A = reduce(hcat, normals)'
                vel = A \ rhs_vals
            elseif length(normals) == 2
                # Edge: 2 constraints + tangent constraint (no motion along edge)
                n1, n2 = normals[1], normals[2]
                t = normalize(cross(n1, n2))
                A = [n1'; n2'; t']
                vel = A \ [rhs_vals[1]; rhs_vals[2]; 0.0]
            else
                # Single constraint (degenerate): fall back to face evaluate_dot
                fi = first(fi for fi in faces_list if haskey(prov, fi))
                p = prov[fi]
                pf = parent_faces[p.body_idx][p.face_idx]
                inv = EGADS.inv_evaluate(pf, pt)
                vel = EGADS.evaluate_dot(pf, inv.param, zeros(2)).result_dot.xyz
            end
            n_shared_vel += 1
        end

        if norm(vel) > 1e-10
            n_face_vel += 1
        else
            n_zero += 1
        end

        idx = 3 * (ig - 1)
        xyz_dot[idx+1] = vel[1]
        xyz_dot[idx+2] = vel[2]
        xyz_dot[idx+3] = vel[3]
    end

    return (xyz_dot=xyz_dot, npts=npts, n_nonzero=n_face_vel, n_shared=n_shared_vel, n_zero=n_zero)
end

# ─── Test 4: Full tessellation velocity field ─────────────────────────────
# Compare dV/dp from analytical velocity field against FD on volume.

function test_tessellation_velocity()
    println("\n═══ Test 4: Full tessellation velocity field ═══")

    pocket_dx = 2.0
    ctx = EGADS.EGADSContext()
    geo = build_boolean_geometry(ctx; pocket_dx=pocket_dx, set_dot=true, dx_dot=1.0)

    tess = EGADS.tessellate(geo.result; max_edge=0.4, sag=0.01, angle=15.0)
    npts = EGADS.tess_status(tess).npts
    println("  Tessellation: $npts surface vertices")

    result_faces = EGADS.get_body_topos(geo.result, EGADS.OClass.FACE)
    outer_faces = EGADS.get_body_topos(geo.outer, EGADS.OClass.FACE)
    pocket_faces = EGADS.get_body_topos(geo.pocket, EGADS.OClass.FACE)
    parent_faces = [outer_faces, pocket_faces]
    prov = match_provenance(geo.result, geo.outer, geo.pocket)

    vel = compute_boolean_velocities(geo, tess, prov, parent_faces)
    println("  Nonzero velocity: $(vel.n_nonzero), shared nodes: $(vel.n_shared), zero: $(vel.n_zero)")

    # Analytical dV/dp via tess_mass_props_dot
    mass_dot = EGADS.tess_mass_props_dot(tess, vel.xyz_dot)
    vol = mass_dot.props[1]
    dvol_dp = mass_dot.props_dot[1]
    println("  Volume = $(@sprintf("%.6f", vol))")
    println("  dV/dp (analytical) = $(@sprintf("%.6e", dvol_dp))")

    # Analytical expectation: V = 4³ - pocket_dx × 2 × 1 = 64 - 2*dx → dV/dp = -2.0
    println("  dV/dp (exact)      = $(@sprintf("%.6e", -2.0))")

    # FD: rebuild at p+h, tessellate, compute volume
    EGADS.close(ctx)
    h = 1e-5
    ctx2 = EGADS.EGADSContext()
    geo_h = build_boolean_geometry(ctx2; pocket_dx=pocket_dx + h)
    tess_h = EGADS.tessellate(geo_h.result; max_edge=0.4, sag=0.01, angle=15.0)
    vol_h = EGADS.tess_mass_props(tess_h)[1]
    EGADS.close(ctx2)

    dvol_dp_fd = (vol_h - vol) / h
    err_fd = abs(dvol_dp - dvol_dp_fd) / abs(dvol_dp_fd + 1e-30)
    err_exact = abs(dvol_dp - (-2.0)) / 2.0
    @printf("  dV/dp (FD)         = %.6e\n", dvol_dp_fd)
    @printf("  Relative error (vs FD)    = %.2e\n", err_fd)
    @printf("  Relative error (vs exact) = %.2e\n", err_exact)

    pass = err_exact < 0.05
    println("  Status: $(pass ? "PASS ✓" : "FAIL ✗")")
    return pass
end

# ─── Volume meshing (EGADS tessellation → TetGen) ────────────────────────

function make_boolean_volume_mesh(body; max_edge=0.3, sag=0.01, angle=15.0, max_vol=0.05)
    tess = EGADS.tessellate(body; max_edge=max_edge, sag=sag, angle=angle)
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
                ig = EGADS.eg_local_to_global(tess.ptr, fi, fm.tris[t, k])
                push!(tris, Int32(ig))
            end
            push!(tri_markers, Int32(fi))
        end
    end

    input = RawTetGenIO{Float64}()
    input.pointlist = pts
    facetlist!(input, reshape(tris, 3, length(tri_markers)))
    input.facetmarkerlist = tri_markers

    output = tetrahedralize(input, "pq1.2a$(max_vol)")
    return output, tess, npts
end

function write_boolean_mesh(filename, tetgen_out)
    pts = tetgen_out.pointlist
    nvol = size(pts, 2)
    tets = tetgen_out.tetrahedronlist
    ntet = size(tets, 2)
    trifaces = tetgen_out.trifacelist
    trimarkers = tetgen_out.trifacemarkerlist
    ntri = size(trifaces, 2)

    open(filename, "w") do io
        println(io, "MeshVersionFormatted 2\n\nDimension 3\n")
        println(io, "Vertices\n$nvol")
        for i in 1:nvol
            @printf(io, "%.15e %.15e %.15e 0\n", pts[1, i], pts[2, i], pts[3, i])
        end
        println(io, "\nTetrahedra\n$ntet")
        for i in 1:ntet
            @printf(io, "%d %d %d %d 1\n", tets[1, i], tets[2, i], tets[3, i], tets[4, i])
        end
        println(io, "\nTriangles\n$ntri")
        for i in 1:ntri
            @printf(io, "%d %d %d %d\n", trifaces[1, i], trifaces[2, i], trifaces[3, i], trimarkers[i])
        end
        println(io, "\nEnd")
    end
end

function write_boolean_palace_config(filename, mesh_path, nfaces)
    config = """{
  "Problem": { "Type": "Eigenmode", "Verbose": 2, "Output": "postpro/boolean" },
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

# ─── Boolean dx/dp at Palace mesh boundary nodes ─────────────────────────
# For each Palace boundary node, project to EGADS parent faces and compute
# velocity via the multi-face constraint solver.

function compute_palace_boolean_dxdp(palace_mesh, parent_faces)
    nv = size(palace_mesh.vertices, 1)
    is_boundary = classify_boundary_vertices(palace_mesh)
    dxdp = zeros(Float64, nv, 3)

    all_pf = vcat([(bi, fi, parent_faces[bi][fi])
                   for bi in 1:length(parent_faces) for fi in eachindex(parent_faces[bi])]...)

    n_projected = 0
    for vi in findall(is_boundary)
        pt = palace_mesh.vertices[vi, :]

        # Find all parent faces within projection tolerance
        close = Tuple{Int, Int, Float64, Vector{Float64}, Vector{Float64}}[]
        for (bi, fi, pf) in all_pf
            inv = EGADS.inv_evaluate(pf, pt)
            dist = norm(inv.results[1:3] - pt)
            if dist < 1e-5
                res = EGADS.evaluate(pf, inv.param)
                n = normalize(cross(res.du, res.dv))
                vel = EGADS.evaluate_dot(pf, inv.param, zeros(2)).result_dot.xyz
                push!(close, (bi, fi, dist, n, vel))
            end
        end
        isempty(close) && continue

        # Deduplicate by normal direction and solve
        normals = Vector{Float64}[]
        rhs_vals = Float64[]
        seen = Set{Tuple{Float64,Float64,Float64}}()

        for (_, _, _, n, vel) in close
            n_key = (round(n[1]; digits=8), round(n[2]; digits=8), round(n[3]; digits=8))
            n_neg = (round(-n[1]; digits=8), round(-n[2]; digits=8), round(-n[3]; digits=8))
            (n_key in seen || n_neg in seen) && continue
            push!(seen, n_key)
            push!(normals, n)
            push!(rhs_vals, dot(n, vel))
        end

        if length(normals) >= 3
            A = reduce(hcat, normals)'
            vel = A \ rhs_vals
        elseif length(normals) == 2
            t = normalize(cross(normals[1], normals[2]))
            A = [normals[1]'; normals[2]'; t']
            vel = A \ [rhs_vals[1]; rhs_vals[2]; 0.0]
        else
            vel = close[1][5]
        end

        dxdp[vi, :] .= vel
        n_projected += 1
    end

    return dxdp, is_boundary, n_projected
end

# ─── Test 5: Full pipeline Taylor test ────────────────────────────────────
# EGADS boolean → TetGen → Palace eigensolve → chain contraction → Taylor test

function run_boolean_pipeline(pocket_dx; max_edge=0.3, max_vol=0.05)
    ctx = EGADS.EGADSContext()
    geo = build_boolean_geometry(ctx; pocket_dx=pocket_dx)
    nfaces = length(EGADS.get_body_topos(geo.result, EGADS.OClass.FACE))

    tetgen_out, _, _ = make_boolean_volume_mesh(geo.result; max_edge=max_edge, max_vol=max_vol)

    dir = joinpath(@__DIR__, "boolean_pipeline")
    mkpath(joinpath(dir, "mesh"))
    mkpath(joinpath(dir, "postpro", "boolean"))

    mesh_file = joinpath(dir, "mesh", "boolean.mesh")
    write_boolean_mesh(mesh_file, tetgen_out)

    config_file = joinpath(dir, "boolean.json")
    write_boolean_palace_config(config_file, "mesh/boolean.mesh", nfaces)

    palace = joinpath(@__DIR__, "..", "build", "palace-build", "palace-arm64.bin")
    @assert isfile(palace) "Palace binary not found: $palace"

    run(pipeline(Cmd(`$palace $config_file`; dir=dir); stdout=devnull, stderr=devnull))

    eig_csv = joinpath(dir, "postpro", "boolean", "eig.csv")
    freq = read_eigenfrequency(eig_csv; mode=1)
    EGADS.close(ctx)
    return freq
end

function test_boolean_taylor()
    println("\n═══ Test 5: Full pipeline Taylor test ═══")
    println("  Parameter: pocket_dx (pocket x-width)")

    pocket_dx = 2.0
    max_edge = 0.3
    max_vol = 0.05

    # ── Step 1: Build geometry, mesh, solve at base parameter ──
    println("\n  Building base geometry (pocket_dx = $pocket_dx)...")
    ctx = EGADS.EGADSContext()
    geo = build_boolean_geometry(ctx; pocket_dx=pocket_dx, set_dot=true, dx_dot=1.0)

    result_faces = EGADS.get_body_topos(geo.result, EGADS.OClass.FACE)
    nfaces = length(result_faces)
    outer_faces = EGADS.get_body_topos(geo.outer, EGADS.OClass.FACE)
    pocket_faces = EGADS.get_body_topos(geo.pocket, EGADS.OClass.FACE)
    parent_faces = [outer_faces, pocket_faces]

    tetgen_out, _, _ = make_boolean_volume_mesh(geo.result; max_edge=max_edge, max_vol=max_vol)

    dir = joinpath(@__DIR__, "boolean_pipeline")
    mkpath(joinpath(dir, "mesh"))
    mkpath(joinpath(dir, "postpro", "boolean"))

    mesh_file = joinpath(dir, "mesh", "boolean.mesh")
    write_boolean_mesh(mesh_file, tetgen_out)

    config_file = joinpath(dir, "boolean.json")
    write_boolean_palace_config(config_file, "mesh/boolean.mesh", nfaces)

    palace = joinpath(@__DIR__, "..", "build", "palace-build", "palace-arm64.bin")
    @assert isfile(palace) "Palace binary not found: $palace"

    println("  Running Palace eigensolve...")
    flush(stdout)
    t0 = time()
    run(pipeline(Cmd(`$palace $config_file`; dir=dir); stdout=devnull, stderr=devnull))
    @printf("  Palace completed in %.1fs\n", time() - t0)

    eig_csv = joinpath(dir, "postpro", "boolean", "eig.csv")
    freq = read_eigenfrequency(eig_csv; mode=1)
    @printf("  f(p) = %.6f GHz\n", freq)

    # ── Step 2: Read mesh sensitivity and compute chain ──
    solb_file = joinpath(dir, "postpro", "boolean", "mesh_sensitivity.solb")
    mesh_sens_file = joinpath(dir, "postpro", "boolean", "mesh_for_sensitivity.mesh")
    @assert isfile(solb_file) "No mesh sensitivity output"

    palace_mesh = read_mesh(mesh_sens_file)
    sol = read_solb(solb_file)
    dfdx = sol.data  # N×3

    println("  Computing dx/dp on boundary via provenance + IFT...")
    dxdp, is_boundary, n_proj = compute_palace_boolean_dxdp(palace_mesh, parent_faces)
    println("  Projected $n_proj / $(count(is_boundary)) boundary nodes")

    # Propagate boundary dx/dp to interior via Laplacian smoothing.
    # The Taylor test perturbs ALL mesh nodes, so the gradient must use the same
    # full-field perturbation: df/dp = dot(dfdx, dxdp_full).
    dxdp_full = propagate_boundary_sensitivity(palace_mesh, dxdp, is_boundary)
    dfdp = dot(vec(dfdx), vec(dxdp_full))
    @printf("  df/dp = %.6e GHz/cm\n", dfdp)

    nv_palace = size(palace_mesh.vertices, 1)
    n_nonzero = count(i -> norm(dxdp[i, :]) > 1e-10, 1:nv_palace)
    @printf("  Nonzero velocity nodes: %d / %d boundary\n", n_nonzero, count(is_boundary))
    EGADS.close(ctx)

    # ── Step 3: Mesh-perturbation Taylor test ──
    # Perturb node coordinates of the SAME mesh to avoid mesh noise.
    # For each ε: x_new = x_old + ε * dxdp_full
    # EGADS and mesh are both in cm — no unit conversion needed.
    println("\n  Mesh-perturbation Taylor test:")
    println("  ε          f(p+ε)       |Δf|            |Δf-ε·df/dp|     rate")

    dir_base = joinpath(@__DIR__, "boolean_pipeline")
    mesh_pert = joinpath(dir_base, "mesh", "boolean_pert.mesh")
    config_pert = joinpath(dir_base, "boolean_pert.json")
    write_boolean_palace_config(config_pert, "mesh/boolean_pert.mesh", nfaces)

    epsilons = [0.1, 0.05, 0.025, 0.0125]
    remainders = Float64[]

    for ε in epsilons
        # Write perturbed mesh
        pts_pert = copy(palace_mesh.vertices)
        for i in 1:nv_palace
            pts_pert[i, :] .+= ε .* dxdp_full[i, :]
        end

        open(mesh_pert, "w") do io
            println(io, "MeshVersionFormatted 2\n\nDimension 3\n")
            println(io, "Vertices\n$nv_palace")
            for i in 1:nv_palace
                @printf(io, "%.15e %.15e %.15e 0\n", pts_pert[i, 1], pts_pert[i, 2], pts_pert[i, 3])
            end
            println(io, "\nTetrahedra\n$(size(palace_mesh.tetrahedra, 1))")
            for i in 1:size(palace_mesh.tetrahedra, 1)
                t = palace_mesh.tetrahedra[i, :]
                @printf(io, "%d %d %d %d %d\n", t[1], t[2], t[3], t[4], palace_mesh.tet_refs[i])
            end
            println(io, "\nTriangles\n$(size(palace_mesh.triangles, 1))")
            for i in 1:size(palace_mesh.triangles, 1)
                t = palace_mesh.triangles[i, :]
                @printf(io, "%d %d %d %d\n", t[1], t[2], t[3], palace_mesh.tri_refs[i])
            end
            println(io, "\nEnd")
        end

        # Solve perturbed mesh
        mkpath(joinpath(dir_base, "postpro", "boolean"))
        run(pipeline(Cmd(`$palace $config_pert`; dir=dir_base); stdout=devnull, stderr=devnull))
        f_pert = read_eigenfrequency(joinpath(dir_base, "postpro", "boolean", "eig.csv"); mode=1)

        delta_f = f_pert - freq
        rem = abs(delta_f - ε * dfdp)
        push!(remainders, rem)

        rate = length(remainders) >= 2 ?
            log(remainders[end-1] / remainders[end]) /
            log(epsilons[length(remainders)-1] / ε) : NaN

        @printf("  %.3f   %10.6f   %12.6e   %12.6e   %.2f\n",
                ε, f_pert, abs(delta_f), rem, rate)
    end

    final_rate = log(remainders[end-1] / remainders[end]) /
                 log(epsilons[length(remainders)-1] / epsilons[end])

    pass = final_rate > 1.5
    @printf("\n  Final convergence rate: %.2f (expect ~2.0 for correct gradient)\n", final_rate)
    println("  Status: $(pass ? "PASS ✓" : "FAIL ✗")")
    return pass
end

# ─── Main ─────────────────────────────────────────────────────────────────

test_boolean_provenance()
test_boolean_face_velocity()
test_boolean_edge_velocity()
test_tessellation_velocity()
test_boolean_taylor()
