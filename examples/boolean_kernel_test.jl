# Boolean sensitivity via EGADSKernel tape integration test.
#
# Builds a box-with-pocket via EGADS extrude + boolean subtraction, manually
# constructs tape entries, and verifies that the tape-based sensitivity pipeline
# produces correct dV/dp through the boolean.
#
# Tests two parameters:
# 1. slab_height: affects only the outer slab extrusion → dV/dp = 16.0
# 2. pocket_height: affects only the pocket extrusion → dV/dp = -4.0
#
# Verification uses tess_mass_props_dot on the EGADS surface tessellation.
#
# Usage: julia --project=path/to/DeviceLayout.jl-egads-dev boolean_kernel_test.jl

using LinearAlgebra, Printf
using EGADS
using DeviceLayout

# Get the extension module for accessing internal functions
const DLE = Base.get_extension(DeviceLayout, :DeviceLayoutEGADSExt)

# ─── Build extruded rectangle ────────────────────────────────────────────────
# Create a rectangular face at (x0,y0,z0) with size (dx,dy), then extrude by
# `dist` along `dir`. Returns the extruded solid body and the source face body.

function make_extruded_rect(ctx, x0, y0, z0, dx, dy, dist, dir)
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
    src_body = EGADS.make_face_body(ctx, face)

    body = EGADS.extrude(src_body, dist, dir)
    return body, src_body
end

# ─── Build boolean geometry (slab - pocket) using extrusions ─────────────────

function build_boolean_geometry(ctx; slab_height=4.0, pocket_height=1.0)
    # Outer slab: rectangle at z=0, size 4×4, extruded by slab_height in z
    slab_body, slab_src = make_extruded_rect(ctx, 0.0, 0.0, 0.0, 4.0, 4.0,
                                              slab_height, [0.0, 0.0, 1.0])

    # Pocket: rectangle at z=2 (fixed, well-interior to slab), size 2×2, extruded by pocket_height in z
    pocket_z = 2.0
    pocket_body, pocket_src = make_extruded_rect(ctx, 1.0, 1.0, pocket_z, 2.0, 2.0,
                                                  pocket_height, [0.0, 0.0, 1.0])

    # Boolean subtraction: result = slab - pocket
    result_bodies = EGADS.solid_boolean(slab_body, pocket_body, EGADS.BoolOp.SUBTRACTION)
    @assert length(result_bodies) >= 1 "Boolean produced no result bodies"
    result = result_bodies[1]

    return (slab=slab_body, pocket=pocket_body, result=result,
            slab_src=slab_src, pocket_src=pocket_src)
end

# ─── Manually construct tape (mimics what kernel_extrude + kernel_cut record) ─

function build_tape(geo; slab_height=4.0, pocket_height=1.0)
    tape = DLE.GeometryTape()
    DLE.start_recording!(tape)

    # Entry 1: slab extrusion
    DLE.record!(tape, DLE.GeometryTapeEntry(
        :extrude,
        Dict{Symbol, Any}(
            :dist => slab_height,
            :dir => [0.0, 0.0, 1.0],
            :offset => [0.0, 0.0, slab_height],
        ),
        geo.slab,
        geo.slab_src,
        Set{Symbol}([:slab_height]),
    ))

    # Entry 2: pocket extrusion
    DLE.record!(tape, DLE.GeometryTapeEntry(
        :extrude,
        Dict{Symbol, Any}(
            :dist => pocket_height,
            :dir => [0.0, 0.0, 1.0],
            :offset => [0.0, 0.0, pocket_height],
        ),
        geo.pocket,
        geo.pocket_src,
        Set{Symbol}([:pocket_height]),
    ))

    # Entry 3: boolean subtraction
    DLE.record!(tape, DLE.GeometryTapeEntry(
        :boolean,
        Dict{Symbol, Any}(
            :op => EGADS.BoolOp.SUBTRACTION,
            :obj_dt => Tuple{Int32, Int32}[],
            :tool_dt => Tuple{Int32, Int32}[],
            :obj_bodies => [geo.slab],
            :tool_bodies => [geo.pocket],
        ),
        geo.result,
        nothing,
        Set{Symbol}(),  # deps inferred during replay
    ))

    DLE.stop_recording!(tape)
    return tape
end

# ─── Run the sensitivity pipeline for a given parameter ──────────────────────

function compute_boolean_sensitivity(tape, result_body, param; seed_val=1.0,
                                      max_edge=0.0, sag=0.0, angle=15.0)
    # Step 1: Replay tape with sensitivity
    DLE.replay_with_sensitivity!(tape, param, seed_val)

    # Step 2: Build extrude fallback dict
    fallback_bodies = Dict{Ptr{Cvoid}, Any}()
    for entry in tape.entries
        if entry.func === :extrude && get(entry.args, :_extrude_analytical, false)
            fallback_bodies[entry.result_body.ptr] = (
                dist = entry.args[:dist]::Float64,
                dir = entry.args[:dir]::Vector{Float64},
                dist_dot = get(entry.args, :_dist_dot, 0.0)::Float64,
                dir_dot = get(entry.args, :_dir_dot, [0.0, 0.0, 0.0])::Vector{Float64},
            )
        end
    end

    # Step 3: Build boolean fallback dict
    boolean_bodies = Dict{Ptr{Cvoid}, Any}()
    for entry in tape.entries
        if entry.func === :boolean && get(entry.args, :_boolean_analytical, false)
            all_parent_bodies = vcat(
                get(entry.args, :obj_bodies, EGADS.EGObject[]),
                get(entry.args, :tool_bodies, EGADS.EGObject[]),
            )

            function _parent_vel_fn(body_idx, parent_face, point)
                parent_body = all_parent_bodies[body_idx]
                fb = get(fallback_bodies, parent_body.ptr, nothing)
                if isnothing(fb)
                    return [0.0, 0.0, 0.0]
                end
                parent_faces = EGADS.get_body_topos(parent_body, EGADS.OClass.FACE)
                source_origin = DLE._find_extrusion_base(parent_faces, fb.dir)
                α = clamp(dot(point - source_origin, fb.dir) / fb.dist, 0.0, 1.0)
                return (α * fb.dist_dot) .* fb.dir .+ (α * fb.dist) .* fb.dir_dot
            end

            boolean_bodies[entry.result_body.ptr] = (
                parent_bodies = all_parent_bodies,
                parent_velocity_fn = _parent_vel_fn,
            )
        end
    end

    # Step 4: Tessellate and extract sensitivity
    tess = EGADS.tessellate(result_body; max_edge=max_edge, sag=sag, angle=angle)
    status = EGADS.tess_status(tess)
    @assert status.state == 1 "Tessellation failed"

    fb_boolean = get(boolean_bodies, result_body.ptr, nothing)
    fb_extrude = get(fallback_bodies, result_body.ptr, nothing)

    if !isnothing(fb_boolean)
        xyz_dot_mat = DLE._boolean_sensitivity(tess, result_body, fb_boolean)
    elseif !isnothing(fb_extrude)
        xyz_dot_mat = DLE._analytical_extrude_sensitivity(tess, result_body, fb_extrude)
    else
        xyz_dot_mat = DLE._evaluate_dot_sensitivity(tess, result_body)
    end

    # Flatten to 3*N vector for tess_mass_props_dot
    npts = status.npts
    xyz_dot = zeros(Float64, 3 * npts)
    for g in 1:npts
        xyz_dot[3*(g-1)+1] = xyz_dot_mat[g, 1]
        xyz_dot[3*(g-1)+2] = xyz_dot_mat[g, 2]
        xyz_dot[3*(g-1)+3] = xyz_dot_mat[g, 3]
    end

    return (tess=tess, xyz_dot=xyz_dot, xyz_dot_mat=xyz_dot_mat, npts=npts)
end

# ─── Reset tape _dot state for a new parameter ──────────────────────────────

function reset_tape_dot_state!(tape)
    for entry in tape.entries
        delete!(entry.args, :_extrude_analytical)
        delete!(entry.args, :_dist_dot)
        delete!(entry.args, :_dir_dot)
        delete!(entry.args, :_boolean_analytical)
        delete!(entry.args, :_parent_entry_indices)
    end
end

# ─── Test: dV/d(slab_height) ────────────────────────────────────────────────

function test_slab_height_sensitivity()
    println("═══ Test 1: dV/d(slab_height) ═══")

    slab_height = 4.0
    pocket_height = 1.0

    ctx = EGADS.EGADSContext()
    geo = build_boolean_geometry(ctx; slab_height=slab_height, pocket_height=pocket_height)
    tape = build_tape(geo; slab_height=slab_height, pocket_height=pocket_height)

    result = compute_boolean_sensitivity(tape, geo.result, :slab_height;
                                          max_edge=0.4, sag=0.01, angle=15.0)

    # Compute dV/dp via tess_mass_props_dot
    mass_dot = EGADS.tess_mass_props_dot(result.tess, result.xyz_dot)
    vol = mass_dot.props[1]
    dvol_dp = mass_dot.props_dot[1]

    # Analytical: V = 4 × 4 × h - 2 × 2 × pocket_height = 16h - 4
    # dV/dh = 16.0
    analytical = 16.0

    @printf("  Volume          = %.4f (expected %.4f)\n", vol, 16 * slab_height - 4)
    @printf("  dV/d(slab_h)    = %.4e (analytical: %.4e)\n", dvol_dp, analytical)

    rel_err = abs(dvol_dp - analytical) / abs(analytical)
    @printf("  Relative error  = %.2e\n", rel_err)

    # FD verification
    h = 1e-5
    ctx2 = EGADS.EGADSContext()
    geo_h = build_boolean_geometry(ctx2; slab_height=slab_height + h, pocket_height=pocket_height)
    tess_h = EGADS.tessellate(geo_h.result; max_edge=0.4, sag=0.01, angle=15.0)
    vol_h = EGADS.tess_mass_props(tess_h)[1]
    dvol_fd = (vol_h - vol) / h
    @printf("  dV/dp (FD)      = %.4e\n", dvol_fd)
    @printf("  FD rel error    = %.2e\n", abs(dvol_dp - dvol_fd) / abs(dvol_fd + 1e-30))
    EGADS.close(ctx2)

    EGADS.close(ctx)

    pass = rel_err < 0.05
    println("  Status: $(pass ? "PASS ✓" : "FAIL ✗")")
    return pass
end

# ─── Test: dV/d(pocket_height) ──────────────────────────────────────────────

function test_pocket_height_sensitivity()
    println("\n═══ Test 2: dV/d(pocket_height) ═══")

    slab_height = 4.0
    pocket_height = 1.0

    ctx = EGADS.EGADSContext()
    geo = build_boolean_geometry(ctx; slab_height=slab_height, pocket_height=pocket_height)
    tape = build_tape(geo; slab_height=slab_height, pocket_height=pocket_height)

    result = compute_boolean_sensitivity(tape, geo.result, :pocket_height;
                                          max_edge=0.4, sag=0.01, angle=15.0)

    mass_dot = EGADS.tess_mass_props_dot(result.tess, result.xyz_dot)
    vol = mass_dot.props[1]
    dvol_dp = mass_dot.props_dot[1]

    # Analytical: V = 16*slab_height - 2*2*pocket_height = 64 - 4*pocket_height
    # dV/d(pocket_height) = -4.0
    analytical = -4.0

    @printf("  Volume              = %.4f\n", vol)
    @printf("  dV/d(pocket_h)      = %.4e (analytical: %.4e)\n", dvol_dp, analytical)

    rel_err = abs(dvol_dp - analytical) / abs(analytical)
    @printf("  Relative error      = %.2e\n", rel_err)

    # FD verification
    h = 1e-5
    ctx2 = EGADS.EGADSContext()
    geo_h = build_boolean_geometry(ctx2; slab_height=slab_height, pocket_height=pocket_height + h)
    tess_h = EGADS.tessellate(geo_h.result; max_edge=0.4, sag=0.01, angle=15.0)
    vol_h = EGADS.tess_mass_props(tess_h)[1]
    dvol_fd = (vol_h - vol) / h
    @printf("  dV/dp (FD)          = %.4e\n", dvol_fd)
    @printf("  FD rel error        = %.2e\n", abs(dvol_dp - dvol_fd) / abs(dvol_fd + 1e-30))
    EGADS.close(ctx2)

    EGADS.close(ctx)

    pass = rel_err < 0.05
    println("  Status: $(pass ? "PASS ✓" : "FAIL ✗")")
    return pass
end

# ─── Test: Provenance coverage ──────────────────────────────────────────────

function test_provenance_coverage()
    println("\n═══ Test 3: Provenance coverage ═══")

    ctx = EGADS.EGADSContext()
    geo = build_boolean_geometry(ctx)

    result_faces = EGADS.get_body_topos(geo.result, EGADS.OClass.FACE)
    println("  Result body: $(length(result_faces)) faces")

    parent_bodies = [geo.slab, geo.pocket]
    prov = DLE._match_provenance(geo.result, parent_bodies)

    println("  Matched: $(length(prov)) / $(length(result_faces)) faces")
    for ri in sort(collect(keys(prov)))
        p = prov[ri]
        parent_name = p.body_idx == 1 ? "slab" : "pocket"
        plane = DLE._face_plane(result_faces[ri])
        @printf("    Face %2d → %6s body, face %d  (n=(%.2f,%.2f,%.2f), d=%.2f)\n",
                ri, parent_name, p.face_idx, plane.normal..., plane.offset)
    end

    unmatched = setdiff(1:length(result_faces), keys(prov))
    pass = isempty(unmatched)
    if !pass
        println("  Unmatched faces: $unmatched")
    end

    EGADS.close(ctx)
    println("  Status: $(pass ? "PASS ✓" : "FAIL ✗")")
    return pass
end

# ─── Main ────────────────────────────────────────────────────────────────────

println("Boolean Kernel Tape Integration Test")
println("====================================\n")

p1 = test_provenance_coverage()
p2 = test_slab_height_sensitivity()
p3 = test_pocket_height_sensitivity()

println("\n═══ Summary ═══")
println("  Provenance coverage:       $(p1 ? "PASS ✓" : "FAIL ✗")")
println("  dV/d(slab_height) = 16.0:  $(p2 ? "PASS ✓" : "FAIL ✗")")
println("  dV/d(pocket_height) = -4:  $(p3 ? "PASS ✓" : "FAIL ✗")")

all_pass = p1 && p2 && p3
all_pass || @warn "Some tests failed"
println("\nDone.")
