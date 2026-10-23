# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generates docs/src/assets/examples/cavity-4.svg, the schematic of the Lorentz-filled
# cylindrical cavity in docs/src/examples/cylinder.md, during docs builds. Panels: (a) the
# axial section of the cavity with the TM010 field, drawn from the exact profiles
# E_z ∝ J₀(p₀₁ r/a) and H_φ ∝ J₁(p₀₁ r/a); (b) the graphical solution of the mode condition
# ω² ε_r(ω) = ε∞ ω_c² for a lossless medium resonant with the cavity mode, with the coupling
# exaggerated so that the two roots and the gap between them are visible. In this lossless,
# resonant case the splitting is exactly ω₀ √(Δε/ε∞).
#
# Sub- and superscripts use dy offsets, since Firefox ignores baseline-shift, and each reset
# is attached to the text that follows, since Chrome drops dy on an empty tspan, as in
# generate_coaxial_schematic.jl.

using Printf

const CAVITY_SCHEMATIC_PATH =
    joinpath(@__DIR__, "src", "assets", "examples", "cavity-4.svg")

# Bessel functions of the first kind by their power series, which converges quickly on the
# interval [0, p₀₁] needed here.
function besselj_series(n::Integer, x::Real)
    term = (x / 2)^n / factorial(n)
    total = term
    for k = 1:40
        term *= -(x / 2)^2 / (k * (k + n))
        total += term
    end
    return total
end

function generate_cavity_schematic(; path::String=CAVITY_SCHEMATIC_PATH)
    W, H = 1200, 560
    BLUE, RED, INK, GREY = "#1f5fa8", "#c0392b", "#1a1a1a", "#666666"
    SERIF = "\"STIX Two Text\", \"Latin Modern Roman\", \"Times New Roman\", Times, serif"
    PEC, MEDIUM, GAP = "url(#pec)", "#eef3e6", "#f4e1de"
    out = String[]
    add(s) = push!(out, s)
    f1(x) = @sprintf("%.1f", x)

    # Primitives
    attrs(kw) = join(("$(replace(string(k), "_" => "-"))='$v'" for (k, v) in kw), " ")
    it(s) = "<tspan font-style='italic'>$s</tspan>"
    function sub(s, then; italic=false)
        style = italic ? " font-style='italic'" : ""
        return "<tspan dy='5' font-size='0.7em'$style>$s</tspan>" *
               "<tspan dy='-5'>$then</tspan>"
    end
    sup(s, then) = "<tspan dy='-7' font-size='0.7em'>$s</tspan><tspan dy='7'>$then</tspan>"
    function text(x, y, body, size=17; anchor="middle", color=INK, kw...)
        return add(
            "<text x='$(f1(x))' y='$(f1(y))' font-size='$size' text-anchor='$anchor' " *
            "fill='$color' $(attrs(kw))>$body</text>"
        )
    end
    function line(x1, y1, x2, y2, color=INK, width=1.5; kw...)
        return add(
            "<line x1='$(f1(x1))' y1='$(f1(y1))' x2='$(f1(x2))' y2='$(f1(y2))' " *
            "stroke='$color' stroke-width='$width' $(attrs(kw))/>"
        )
    end
    arrow(x1, y1, x2, y2, color=INK, width=1.5, head="ink") =
        line(x1, y1, x2, y2, color, width; marker_end="url(#head-$head)")
    dim(x1, y1, x2, y2) = line(
        x1,
        y1,
        x2,
        y2,
        INK,
        1.0;
        marker_start="url(#dim-start)",
        marker_end="url(#dim-end)"
    )
    function rect(x, y, w, h, fill; stroke=INK, width=1.4)
        return add(
            "<rect x='$(f1(x))' y='$(f1(y))' width='$(f1(w))' height='$(f1(h))' " *
            "fill='$fill' stroke='$stroke' stroke-width='$width'/>"
        )
    end
    circle(cx, cy, r, color; fill="white") = add(
        "<circle cx='$(f1(cx))' cy='$(f1(cy))' r='$(f1(r))' fill='$fill' " *
        "stroke='$color' stroke-width='1.4'/>"
    )
    function cross(cx, cy, r, color)  # Vector into the page
        circle(cx, cy, r, color)
        k = 0.6 * r
        line(cx - k, cy - k, cx + k, cy + k, color, 1.4)
        return line(cx - k, cy + k, cx + k, cy - k, color, 1.4)
    end
    function dot(cx, cy, r, color)  # Vector out of the page
        circle(cx, cy, r, color)
        return add(
            "<circle cx='$(f1(cx))' cy='$(f1(cy))' r='$(f1(max(1.5, 0.3 * r)))' " *
            "fill='$color'/>"
        )
    end
    function polyline(pts, color, width=2.4; kw...)
        p = join(("$(f1(x)),$(f1(y))" for (x, y) in pts), " ")
        return add(
            "<polyline points='$p' fill='none' stroke='$color' stroke-width='$width' " *
            "stroke-linejoin='round' stroke-linecap='round' $(attrs(kw))/>"
        )
    end
    function panel_title(x, y, tag, caption)
        text(x, y, tag, 19; anchor="start", font_weight="bold")
        return text(x + 36, y, caption, 17; anchor="start", color=GREY)
    end

    # Document
    add(
        "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 $W $H' width='$W' " *
        "height='$H' font-family='$SERIF' role='img' aria-labelledby='title desc'>"
    )
    add("<title id='title'>Cylindrical cavity filled with a Lorentz medium</title>")
    add(
        "<desc id='desc'>Axial section of a PEC cylindrical cavity filled with a Lorentz " *
        "medium, with the TM010 electric and magnetic field profiles, and the graphical " *
        "solution of the mode condition omega squared times the relative permittivity " *
        "equals epsilon infinity times the cavity frequency squared, which has one root " *
        "on each side of a gap of negative permittivity between the resonance and the " *
        "longitudinal frequency.</desc>"
    )
    add("<defs>")
    for (name, color) in (("ink", INK), ("red", RED), ("blue", BLUE))
        add(
            "<marker id='head-$name' viewBox='0 0 10 10' refX='8.5' refY='5' " *
            "markerWidth='6.5' markerHeight='6.5' orient='auto'><path " *
            "d='M0,0.5 L10,5 L0,9.5 L2.8,5 z' fill='$color'/></marker>"
        )
    end
    for (name, refx, d) in
        (("end", 10, "M0,1.8 L10,5 L0,8.2 z"), ("start", 0, "M10,1.8 L0,5 L10,8.2 z"))
        add(
            "<marker id='dim-$name' viewBox='0 0 10 10' refX='$refx' refY='5' " *
            "markerWidth='8' markerHeight='8' orient='auto'><path d='$d' " *
            "fill='$INK'/></marker>"
        )
    end
    add(
        "<pattern id='pec' width='7' height='7' patternUnits='userSpaceOnUse' " *
        "patternTransform='rotate(45)'><rect width='7' height='7' fill='white'/>" *
        "<line x1='0' y1='0' x2='0' y2='7' stroke='$INK' stroke-width='1.5'/></pattern>"
    )
    add("</defs>")
    add("<rect width='$W' height='$H' fill='white'/>")

    # (a) Axial section with the TM010 field
    panel_title(40, 36, "(a)", "Axial section, TM" * sub("010", "") * " field")
    cx, ytop, R, wall = 290.0, 100.0, 180.0, 12.0  # Axis, top end, radius (px)
    ybot = ytop + 2R  # d = 2a
    p01 = 2.404825557695773
    rect(cx - R, ytop, 2R, ybot - ytop, MEDIUM; stroke="none")
    rect(cx - R - wall, ytop - wall, 2R + 2wall, wall, PEC)
    rect(cx - R - wall, ybot, 2R + 2wall, wall, PEC)
    rect(cx - R - wall, ytop, wall, ybot - ytop, PEC)
    rect(cx + R, ytop, wall, ybot - ytop, PEC)
    line(cx, ytop - 22, cx, ybot + 22, GREY, 0.9; stroke_dasharray="16 4 3 4")  # Axis

    # E_z ∝ J₀(p₀₁ r/a): uniform in z, largest on the axis, zero on the side wall.
    ymid = 0.5 * (ytop + ybot)
    for r in range(-0.88, 0.88; length=9)
        len = 0.5 * (ybot - ytop - 64) * besselj_series(0, p01 * abs(r))
        x = cx + r * R
        arrow(x, ymid + len, x, ymid - len, RED, 2.2, "red")
    end
    # H_φ ∝ J₁(p₀₁ r/a), circulating about the axis: into the page on one side, out of it
    # on the other, with symbol size ∝ |H_φ|.
    for r in (0.25, 0.5, 0.77), y in (ytop + 30, ybot - 30)
        size = 2.5 + 7.5 * besselj_series(1, p01 * r) / besselj_series(1, 1.8412)
        cross(cx + r * R, y, size, BLUE)
        dot(cx - r * R, y, size, BLUE)
    end
    text(cx + 0.47R, ymid + 6, it("E") * sub("z", ""; italic=true), 19; color=RED)
    text(
        cx - R + 8,
        ytop + 66,
        it("H") * sub("φ", ""; italic=true),
        19;
        anchor="start",
        color=BLUE
    )

    # Dimensions and labels.
    yd = ybot + wall + 30
    dim(cx + 1, yd, cx + R - 1, yd)
    text(cx + 0.5R, yd - 8, it("a"), 18)
    xd = cx + R + wall + 26
    dim(xd, ytop + 1, xd, ybot - 1)
    text(xd + 12, ymid + 6, it("d") * " = 2" * it("a"), 18; anchor="start")

    # Legend
    lx = 40.0
    for (k, (fill, label)) in enumerate((
        (MEDIUM, "Lorentz medium, " * it("ε") * sub("r", "(") * it("ω") * ")"),
        (PEC, "PEC")
    ))
        ly = yd - 4 + 24 * (k - 1)
        rect(lx, ly - 11, 26, 13, fill; stroke=INK, width=1.0)
        text(lx + 36, ly, label, 15; anchor="start")
    end

    # (b) Graphical solution of the mode condition
    panel_title(560, 36, "(b)", "Mode condition, resonant medium")
    κ = 0.25  # Exaggerated Δε/ε∞ for a visible splitting; ω₀ = ω_c
    xL = sqrt(1 + κ)
    xm = (sqrt(κ + 4) - sqrt(κ)) / 2  # Roots of (x² - 1)² = κ x², x = ω/ω_c
    xp = (sqrt(κ + 4) + sqrt(κ)) / 2
    gx0, gx1, gy0, gy1 = 640.0, 1150.0, 92.0, 470.0  # Plot box (px)
    xmin, xmax, ymin, ymax = 0.55, 1.6, -0.35, 2.6
    X(x) = gx0 + (x - xmin) / (xmax - xmin) * (gx1 - gx0)
    Y(y) = gy1 - (y - ymin) / (ymax - ymin) * (gy1 - gy0)

    # Gap ω₀ < ω < ω_L where ε_r < 0 and no transverse mode exists.
    rect(X(1.0), gy0, X(xL) - X(1.0), gy1 - gy0, GAP; stroke="none")
    text(0.5 * (X(1.0) + X(xL)), gy0 + 22, it("ε") * sub("r", " &lt; 0"), 15; color=RED)
    # Axes
    line(gx0, Y(0.0), gx1, Y(0.0), INK, 1.2)
    arrow(gx0, gy1, gx0, gy0 - 10, INK, 1.2)
    arrow(gx1 - 30, Y(0.0), gx1 + 10, Y(0.0), INK, 1.2)
    text(gx1 + 4, Y(0.0) + 26, it("ω") * " / " * it("ω") * sub("c", ""; italic=true), 17)
    text(
        gx0 + 8,
        gy0 - 2,
        it("ω") *
        sup("2", "") *
        it("ε") *
        sub("r", "(") *
        it("ω") *
        ") / " *
        it("ε") *
        sub("∞", "") *
        it("ω") *
        sub("c", ""; italic=true) *
        sup("2", ""),
        17;
        anchor="start"
    )

    # Without the medium resonance: ω²/ω_c², crossing 1 at the cavity mode.
    bare = [(X(x), Y(x^2)) for x in range(xmin, 1.6; length=60) if x^2 <= ymax]
    polyline(bare, GREY, 1.6; stroke_dasharray="6 4")
    # With it: x² (1 + κ/(1 - x²)), with a pole at the resonance and a zero at ω_L.
    g(x) = x^2 * (1 + κ / (1 - x^2))
    left = [(X(x), Y(g(x))) for x in range(xmin, 1 - 1e-3; length=200) if g(x) <= ymax]
    right = [(X(x), Y(g(x))) for x in range(1 + 1e-3, xmax; length=300) if g(x) >= ymin]
    polyline(left, BLUE, 2.6)
    polyline(right, BLUE, 2.6)
    line(X(1.0), gy0, X(1.0), gy1, GREY, 1.0; stroke_dasharray="3 3")  # Pole at ω₀

    # Mode condition = 1 and its two roots.
    line(gx0, Y(1.0), gx1 - 20, Y(1.0), INK, 1.4)
    text(gx1 - 16, Y(1.0) + 5, "1", 16; anchor="start")
    for (x, label) in ((xm, "−"), (xp, "+"))
        line(X(x), Y(1.0), X(x), Y(0.0), INK, 0.9; stroke_dasharray="2 3")
        circle(X(x), Y(1.0), 5.5, INK; fill=INK)
        text(X(x), Y(0.0) + 24, it("ω") * sub(label, ""), 17)
    end
    text(X(1.0), Y(0.0) + 24, it("ω") * sub("0", ""), 17; color=GREY)
    text(
        X(xL) + 12,
        Y(0.0) + 24,
        it("ω") * sub("L", ""; italic=true),
        17;
        anchor="start",
        color=GREY
    )
    # Splitting, below the axis labels
    ys = Y(0.0) + 52
    for x in (xm, xp)
        line(X(x), Y(0.0) + 32, X(x), ys + 6, INK, 0.9; stroke_dasharray="2 3")
    end
    dim(X(xm) + 1, ys, X(xp) - 1, ys)
    text(
        0.5 * (X(xm) + X(xp)),
        ys + 24,
        "Δ" *
        it("ω") *
        " = " *
        it("ω") *
        sub("0", " √(") *
        "Δ" *
        it("ε") *
        " / " *
        it("ε") *
        sub("∞", ")"),
        16
    )
    # Curve labels
    text(X(0.66), Y(0.66^2) + 26, "no resonance", 15; color=GREY)
    text(X(1.42), Y(g(1.42)) - 14, "with resonance", 15; anchor="end", color=BLUE)

    add("</svg>")
    mkpath(dirname(path))
    write(path, join(out, "\n") * "\n")
    return path
end
