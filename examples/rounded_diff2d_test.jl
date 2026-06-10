# Gate test: ForwardDiff through Rounded(radius)(polygon) and difference2d
#
# Tests whether ForwardDiff.Dual numbers propagate correctly through DeviceLayout's
# Rounded polygon style. This is a prerequisite for adding island_rounding as param 11
# in the transmon sensitivity pipeline.
#
# Test 1: Rounded rectangle — differentiate vertex coords w.r.t. (width, height, radius)
# Test 2: Rounded difference2d — Rounded(rnd)(difference2d(outer, inner))
#          (exact pattern from ExampleRectangleIsland._geometry!)
#
# Gate criteria:
#   - Stable vertex count under Dual perturbation
#   - All Jacobian entries finite
#   - Max error vs central FD < 1e-4
#
# If PASS → island_rounding can be added as param 11
# If FAIL → investigate vertex count instability or non-differentiable operations

using DeviceLayout, DeviceLayout.PreferredUnits
using ForwardDiff
using Unitful
using Printf, LinearAlgebra

# ─── Test 1: Rounded rectangle ───────────────────────────────────────────────

function rounded_rect_coords(params)
    # params = [width_μm, height_μm, rounding_radius_μm]
    w = params[1] * 1μm
    h = params[2] * 1μm
    rnd = params[3] * 1μm

    r = Rectangle(Point(zero(w), zero(w)), Point(w, h))
    styled = Rounded(rnd)(r)
    result = to_polygons(styled)
    # to_polygons returns Polygon (not vector) for simple shapes
    polys = result isa AbstractVector ? result : [result]

    coords = eltype(params)[]
    for poly in polys
        for pt in points(poly)
            push!(coords, Unitful.ustrip(μm, getx(pt)))
            push!(coords, Unitful.ustrip(μm, gety(pt)))
        end
    end
    return coords
end

# ─── Test 2: Rounded difference2d (transmon island pattern) ──────────────────

function rounded_island_coords(params)
    # params = [cap_width_μm, cap_length_μm, cap_gap_μm, junction_gap_μm, rounding_μm]
    # Mirrors ExampleRectangleIsland._geometry! from DeviceLayout
    cap_width = params[1] * 1μm
    cap_length = params[2] * 1μm
    cap_gap = params[3] * 1μm
    junction_gap = params[4] * 1μm
    rnd = params[5] * 1μm

    # Inner rectangle (island metal)
    r1 = Rectangle(
        Point(-cap_width / 2, junction_gap),
        Point(cap_width / 2, junction_gap + cap_length)
    )
    # Outer rectangle (gap boundary)
    r2 = Rectangle(
        Point(-cap_width / 2 - cap_gap, zero(cap_gap)),
        Point(cap_width / 2 + cap_gap, junction_gap + cap_length + cap_gap)
    )

    diff = difference2d(r2, r1)
    styled = Rounded(rnd)(diff)
    polys = to_polygons(styled)

    # Flatten all polygon vertices into a single coordinate vector
    coords = eltype(params)[]
    for poly in polys
        for pt in points(poly)
            push!(coords, Unitful.ustrip(μm, getx(pt)))
            push!(coords, Unitful.ustrip(μm, gety(pt)))
        end
    end
    return coords
end

# ─── Jacobian verification ────────────────────────────────────────────────────

function check_vertex_count_stability(f, p0, name; δ=1e-7)
    n0 = length(f(p0))
    stable = true
    for j in 1:length(p0)
        pp = copy(p0); pp[j] += δ
        pm = copy(p0); pm[j] -= δ
        np = length(f(pp))
        nm = length(f(pm))
        if np != n0 || nm != n0
            @printf("  ⚠ Vertex count unstable for param %d: %d / %d / %d (−δ/0/+δ)\n",
                    j, nm ÷ 2, n0 ÷ 2, np ÷ 2)
            stable = false
        end
    end
    if stable
        @printf("  Vertex count stable: %d vertices across all perturbations\n", n0 ÷ 2)
    end
    return stable
end

function test_jacobian_vs_fd(f, p0, names; δ=1e-7, tol=1e-4)
    # ForwardDiff Jacobian
    println("  Computing ForwardDiff Jacobian...")
    local J
    try
        J = ForwardDiff.jacobian(f, p0)
    catch e
        println("  ✗ ForwardDiff.jacobian FAILED: $e")
        return false
    end
    ncoords = size(J, 1)
    nparams = size(J, 2)
    @printf("  Jacobian size: %d × %d\n", ncoords, nparams)

    # Check for NaN/Inf
    nnan = count(isnan, J)
    ninf = count(isinf, J)
    if nnan > 0 || ninf > 0
        @printf("  ✗ Non-finite entries: %d NaN, %d Inf\n", nnan, ninf)
        return false
    end
    println("  All entries finite ✓")

    # Central FD comparison
    println("  Comparing to central FD (δ=$δ)...")
    all_pass = true
    for j in 1:nparams
        pp = copy(p0); pp[j] += δ
        pm = copy(p0); pm[j] -= δ
        fp = f(pp)
        fm = f(pm)
        if length(fp) != ncoords || length(fm) != ncoords
            @printf("    %-15s: SKIP (vertex count changed: %d vs %d vs %d)\n",
                    names[j], length(fm) ÷ 2, ncoords ÷ 2, length(fp) ÷ 2)
            all_pass = false
            continue
        end
        fd_col = (fp - fm) / (2δ)
        err = maximum(abs.(J[:, j] - fd_col))
        pass = err < tol
        status = pass ? "✓" : "✗"
        @printf("    %-15s: max|J-FD| = %.2e  %s\n", names[j], err, status)
        all_pass &= pass
    end
    return all_pass
end

# ─── Main ────────────────────────────────────────────────────────────────────

function main()
    println("═══ Gate Test: ForwardDiff through Rounded polygons ═══\n")

    # --- Test 1: Rounded rectangle ---
    println("── Test 1: Rounded Rectangle ──")
    p1 = [100.0, 80.0, 5.0]  # width, height, rounding_radius (μm)
    names1 = ["width", "height", "rnd_radius"]

    println("  Nominal: width=100, height=80, rnd=5 μm")
    stable1 = check_vertex_count_stability(rounded_rect_coords, p1, "rounded_rect")
    pass1 = stable1 && test_jacobian_vs_fd(rounded_rect_coords, p1, names1)
    println("  Result: $(pass1 ? "PASS ✓" : "FAIL ✗")\n")

    # --- Test 2: Rounded difference2d (island pattern) ---
    println("── Test 2: Rounded difference2d (island pattern) ──")
    # Matches ExampleRectangleIsland defaults
    p2 = [24.0, 620.0, 30.0, 12.0, 3.0]  # cap_width, cap_length, cap_gap, junction_gap, rounding
    names2 = ["cap_width", "cap_length", "cap_gap", "junction_gap", "rounding"]

    println("  Nominal: cap_w=24, cap_l=620, gap=30, jg=12, rnd=3 μm")
    stable2 = check_vertex_count_stability(rounded_island_coords, p2, "rounded_island")
    pass2 = stable2 && test_jacobian_vs_fd(rounded_island_coords, p2, names2)
    println("  Result: $(pass2 ? "PASS ✓" : "FAIL ✗")\n")

    # --- Summary ---
    println("═══ Gate Test Summary ═══")
    println("  Rounded rectangle:      $(pass1 ? "PASS ✓" : "FAIL ✗")")
    println("  Rounded difference2d:   $(pass2 ? "PASS ✓" : "FAIL ✗")")

    overall = pass1 && pass2
    println()
    if overall
        println("  ✓ All tests pass. island_rounding is safe to add as param 11.")
        println("    Next: add rounding param to island_gap_coords_μm in the full pipeline.")
    elseif pass1 && !pass2
        println("  △ Rounded works alone but not with difference2d.")
        println("    Known issue: difference2d produces self-touching polygons.")
        println("    Workaround: use EGADS boolean subtraction (separate bodies) instead.")
    else
        println("  ✗ Rounded polygon rounding is not ForwardDiff-compatible.")
        println("    Investigate: vertex count instability or non-differentiable code paths.")
        println("    Proceeding with 10 params (no rounding) for now.")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
