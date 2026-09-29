# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generates docs/src/assets/examples/coaxial-4.svg, the schematic of the superconducting
# coaxial line in docs/src/examples/coaxial.md, during docs builds. Panels: (a) the axial
# section, with the driven end cap, PEC short, current loop and B_φ; (b) the cross-section
# seen from the driven end; (c) the field and current across a one-sided London film, drawn
# from the exact profiles B ∝ sinh((d - x)/λ) and J ∝ cosh((d - x)/λ).
#
# Signs follow the right-hand rule for the current loop: radially out across the driven
# cap, +z along the outer conductor, -z along the inner conductor. In the axial section (z
# to the right, r up) B_φ points into the page above the axis and out of it below; seen
# from the driven end the inner current comes toward the viewer and B circulates
# counterclockwise.
#
# Sub- and superscripts use dy offsets, since Firefox ignores baseline-shift, and each reset
# is attached to the text that follows, since Chrome drops dy on an empty tspan. The result
# renders the same in Chrome and Firefox.

using Printf

const COAX_SCHEMATIC_PATH = joinpath(@__DIR__, "src", "assets", "examples", "coaxial-4.svg")

function generate_coaxial_schematic(; path::String=COAX_SCHEMATIC_PATH)
    W, H = 1200, 720
    BLUE, RED, INK, GREY = "#1f5fa8", "#c0392b", "#1a1a1a", "#666666"
    SERIF = "\"STIX Two Text\", \"Latin Modern Roman\", \"Times New Roman\", Times, serif"
    FILM, PEC, CORE = "url(#film)", "url(#pec)", "#f1f1f1"
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
    mu0(then) = "μ" * sub("0", then)
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
    circle(cx, cy, r, color) = add(
        "<circle cx='$(f1(cx))' cy='$(f1(cy))' r='$(f1(r))' fill='white' " *
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
    polar(cx, cy, r, deg) = (cx + r * cosd(deg), cy - r * sind(deg))

    # Document
    add(
        "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 $W $H' width='$W' " *
        "height='$H' font-family='$SERIF' role='img' aria-labelledby='title desc'>"
    )
    add("<title id='title'>Shorted superconducting coaxial line</title>")
    add(
        "<desc id='desc'>Axial section and cross-section of a shorted coaxial line " *
        "driven by a radial surface current on one end cap, with the field and current " *
        "profiles across a one-sided London film of thickness d and penetration depth " *
        "lambda.</desc>"
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
        "<pattern id='film' width='5' height='5' patternUnits='userSpaceOnUse' " *
        "patternTransform='rotate(45)'><rect width='5' height='5' fill='#dde6ef'/>" *
        "<line x1='0' y1='0' x2='0' y2='5' stroke='#7b90a6' stroke-width='1.1'/></pattern>"
    )
    add(
        "<pattern id='pec' width='7' height='7' patternUnits='userSpaceOnUse' " *
        "patternTransform='rotate(45)'><rect width='7' height='7' fill='white'/>" *
        "<line x1='0' y1='0' x2='0' y2='7' stroke='$INK' stroke-width='1.5'/></pattern>"
    )
    add("</defs>")
    add("<rect width='$W' height='$H' fill='white'/>")

    # (a) Axial section
    panel_title(40, 36, "(a)", "Axial section (not to scale)")
    x0, x1, y0 = 200.0, 880.0, 196.0  # Driven cap, short, axis
    ra, rb, tf = 34.0, 112.0, 7.0     # Radii and exaggerated film thickness (px)

    # Field-free metal: the inner core (r < a - d) and the region outside the outer film.
    rect(x0, y0 - ra + tf, x1 - x0, 2 * (ra - tf), CORE; stroke="none")
    for s in (-1, 1)
        rect(x0, y0 - s * (rb + tf) - (s > 0 ? 14 : 0), x1 - x0, 14, CORE; stroke="none")
    end
    # Films: outer r ∈ [b, b + d], inner r ∈ [a - d, a].
    for s in (-1, 1)
        rect(x0, s > 0 ? y0 - rb - tf : y0 + rb, x1 - x0, tf, FILM)
        rect(x0, s > 0 ? y0 - ra : y0 + ra - tf, x1 - x0, tf, FILM)
    end
    line(x0 - 36, y0, x1 + 50, y0, GREY, 0.9; stroke_dasharray="16 4 3 4")  # Axis

    # Driven end cap: "+R" surface current, radially out from the inner to the outer.
    for s in (-1, 1)
        ya, yb = y0 - s * ra, y0 - s * rb
        line(x0, ya, x0, yb, RED, 2.8)
        for f in (0.3, 0.74)
            yy = ya + (yb - ya) * f
            arrow(x0, yy + s * 10, x0, yy - s * 10, RED, 2.4, "red")
        end
    end
    text(x0 - 16, y0 - rb + 20, "Driven end cap", 16; anchor="end")
    text(x0 - 16, y0 - rb + 40, "SurfaceCurrent, +R", 14.5; anchor="end", color=GREY)

    # Far end: PEC short across the annulus.
    rect(x1, y0 - rb - tf, 13, 2 * (rb + tf), PEC)
    text(x1 + 24, y0 - 0.5 * (ra + rb) + 6, "PEC short", 16; anchor="start")

    # Current: +z just outside the outer film, -z inside the inner core, so each arrow is
    # on the metal side of its film.
    for s in (-1, 1), xc in (340, 540, 740)
        yo = y0 - s * (rb + tf + 7)
        yi = y0 - s * (ra - tf - 9)
        arrow(xc - 28, yo, xc + 28, yo, RED, 2.4, "red")
        arrow(xc + 28, yi, xc - 28, yi, RED, 2.4, "red")
    end
    text(640, y0 - rb - tf - 1, it("I"), 18; color=RED)
    text(640, y0 - (ra - tf - 9) + 6, it("I"), 18; color=RED)

    # Azimuthal B in the gap, with symbol size ∝ 1/r.
    for r in (48.0, 70.0, 94.0), xc in (270, 440, 840)
        cross(xc, y0 - r, 9.0 * 48.0 / r, BLUE)
        dot(xc, y0 + r, 9.0 * 48.0 / r, BLUE)
    end
    text(
        640,
        y0 - 0.5 * (ra + rb) + 7,
        it("B") * sub("φ", " = "; italic=true) * mu0(it(" I")) * " / 2π" * it("r"),
        18;
        color=BLUE
    )
    text(640, y0 + 0.5 * (ra + rb) + 6, "vacuum", 16; color=GREY)

    # Legend under the radius dimensions.
    legend = (
        (FILM, INK, "Superconductor film (" * it("d") * ", " * it("λ") * ")"),
        (CORE, "#bdbdbd", "Field-free metal"),
        (PEC, INK, "PEC")
    )
    for (k, (fill, edge, label)) in enumerate(legend)
        ly = y0 + rb + 40 + 24 * (k - 1)
        rect(x1 + 60, ly - 11, 26, 13, fill; stroke=edge, width=1.0)
        text(x1 + 96, ly, label, 15; anchor="start")
    end

    # Dimensions.
    xa, xb = 1010.0, 1060.0
    for (yy, xe) in ((y0, xb), (y0 + ra, xa), (y0 + rb, xb))
        line(x1 + 18, yy, xe + 12, yy, GREY, 0.8; stroke_dasharray="3 3")
    end
    dim(xa, y0 + 2, xa, y0 + ra - 2)
    text(xa - 12, y0 + ra / 2 + 6, it("a"), 18; anchor="end")
    dim(xb, y0 + 2, xb, y0 + rb - 2)
    text(xb + 12, y0 + rb / 2 + 6, it("b"), 18; anchor="start")
    yl = y0 + rb + tf + 44
    line(x0, y0 + rb + tf + 18, x0, yl + 8, GREY, 0.8)
    line(x1, y0 + rb + tf + 4, x1, yl + 8, GREY, 0.8)
    dim(x0 + 1, yl, x1 - 1, yl)
    text(0.5 * (x0 + x1), yl - 8, it("ℓ"), 19)
    text(x0, yl + 26, it("z") * " = 0", 15; color=GREY)
    text(x1, yl + 26, it("z") * " = " * it("ℓ"), 15; color=GREY)

    # (b) Cross-section
    panel_title(40, 432, "(b)", "Cross-section, seen from the driven end")
    cx, cy = 245.0, 580.0
    Ra, Rb, Tf = 34.0, 112.0, 8.0
    function ring(ro, ri, fill)
        return add(
            "<path d='M$(cx - ro),$cy a$ro,$ro 0 1,0 $(2 * ro),0 " *
            "a$ro,$ro 0 1,0 $(-2 * ro),0 M$(cx - ri),$cy a$ri,$ri 0 1,0 $(2 * ri),0 " *
            "a$ri,$ri 0 1,0 $(-2 * ri),0 z' fill='$fill' fill-rule='evenodd' " *
            "stroke='$INK' stroke-width='1.4'/>"
        )
    end
    add("<circle cx='$cx' cy='$cy' r='$(Ra - Tf)' fill='$CORE'/>")
    ring(Rb + Tf, Rb, FILM)
    ring(Ra, Ra - Tf, FILM)
    for r in (54.0, 74.0, 94.0)  # Counterclockwise field lines
        add(
            "<circle cx='$cx' cy='$cy' r='$r' fill='none' stroke='$BLUE' " *
            "stroke-width='1.3' stroke-dasharray='6 4'/>"
        )
        for deg in (60.0, 180.0, 300.0)
            px, py = polar(cx, cy, r, deg)
            tx, ty = -sind(deg), -cosd(deg)  # Counterclockwise tangent, screen coordinates
            arrow(px - 5 * tx, py - 5 * ty, px + 6 * tx, py + 6 * ty, BLUE, 1.7, "blue")
        end
    end
    for deg in (90.0, 210.0, 330.0)  # Inner current toward the viewer
        dot(polar(cx, cy, Ra - Tf / 2, deg)..., 3.8, RED)
    end
    for deg in (30.0, 150.0, 270.0)  # Outer current away from the viewer
        cross(polar(cx, cy, Rb + Tf / 2, deg)..., 4.2, RED)
    end
    add("<circle cx='$cx' cy='$cy' r='2' fill='$INK'/>")
    # Radii along directions clear of the field-line arrows.
    arrow(cx, cy, polar(cx, cy, Ra - 1, 135.0)..., INK, 1.1)
    lx, ly = polar(cx, cy, 44, 135.0)
    text(lx, ly + 6, it("a"), 17)
    arrow(cx, cy, polar(cx, cy, Rb - 1, 225.0)..., INK, 1.1)
    lx, ly = polar(cx, cy, Rb + Tf + 14, 225.0)
    text(lx, ly + 6, it("b"), 17)
    # Callouts on the right.
    tx0 = 420.0
    line(polar(cx, cy, Rb + Tf, 18.0)..., tx0 - 6, 494, INK, 0.9)
    text(
        tx0,
        499,
        it("K") * " = " * it("I") * " / 2π" * it("b"),
        16;
        anchor="start",
        color=RED
    )
    line(polar(cx, cy, 94, -8.0)..., tx0 - 6, 580, INK, 0.9)
    text(tx0, 585, it("B") * sub("φ", ""; italic=true), 17; anchor="start", color=BLUE)
    line(polar(cx, cy, Ra, -40.0)..., tx0 - 6, 664, INK, 0.9)
    text(
        tx0,
        669,
        it("K") * " = " * it("I") * " / 2π" * it("a"),
        16;
        anchor="start",
        color=RED
    )

    # (c) One-sided London film
    panel_title(590, 432, "(c)", "One-sided London film")
    gx0, fx0, fx1, gx1 = 610.0, 790.0, 910.0, 1000.0  # Vacuum | film | field-free
    ybase, ytop = 660.0, 520.0
    λ_d = 0.45  # λ/d
    rect(fx0, ytop - 26, fx1 - fx0, ybase - ytop + 26, FILM)
    line(gx0, ybase, gx1, ybase, INK, 1.4)
    for (xc, label) in (
        (0.5 * (gx0 + fx0), "vacuum"),
        (0.5 * (fx0 + fx1), "film"),
        (0.5 * (fx1 + gx1) + 4, "field-free")
    )
        text(xc, ytop - 36, label, 15.5; color=GREY)
    end

    # |B| is uniform in the gap and falls as sinh((d - x)/λ)/sinh(d/λ) through the film; the
    # current density, normalized to its value at the surface, is cosh((d - x)/λ)/cosh(d/λ).
    n, d_px = 80, fx1 - fx0
    ts = range(0, 1; length=n + 1)
    b_profile(t) = sinh((1 - t) / λ_d) / sinh(1 / λ_d)  # t = x/d
    j_profile(t) = cosh((1 - t) / λ_d) / cosh(1 / λ_d)
    B = [
        (gx0 + 8, ytop);
        (fx0, ytop);
        [(fx0 + t * d_px, ybase - b_profile(t) * (ybase - ytop)) for t in ts];
        (gx1 - 6, ybase)
    ]
    J = [(fx0 + t * d_px, ybase - 0.72 * j_profile(t) * (ybase - ytop)) for t in ts]
    polyline(B, BLUE, 2.6)
    polyline(J, RED, 2.3; stroke_dasharray="7 4")
    text(0.5 * (gx0 + fx0), ytop - 10, "|" * it("B") * "|", 17; color=BLUE)
    y_j = ybase - 0.72 / cosh(1 / λ_d) * (ybase - ytop) - 12
    text(fx1 - 8, y_j, it("J"), 17; anchor="end", color=RED)

    x_λ = fx0 + λ_d * d_px
    line(fx0, ybase, fx0, ybase + 44, INK, 0.9)
    line(x_λ, ybase, x_λ, ybase + 22, INK, 0.9)
    line(fx1, ybase, fx1, ybase + 44, INK, 0.9)
    dim(fx0 + 1, ybase + 16, x_λ - 1, ybase + 16)
    text(0.5 * (fx0 + x_λ), ybase + 34, it("λ"), 17)
    dim(fx0 + 1, ybase + 40, fx1 - 1, ybase + 40)
    text(fx1 + 12, ybase + 46, it("d"), 17; anchor="start")

    fx = 1016.0
    ksq =
        it("L") *
        sub("ksq", " = ") *
        mu0(it("λ")) *
        " coth(" *
        it("d") *
        "/" *
        it("λ") *
        ")"
    text(fx, 530, ksq, 17.5; anchor="start")
    text(
        fx,
        574,
        "thin, " * it("d") * " ≪ " * it("λ") * ":",
        15;
        anchor="start",
        color=GREY
    )
    text(fx + 14, 598, "→ " * mu0(it("λ")) * sup("2", "/") * it("d"), 17; anchor="start")
    text(
        fx,
        632,
        "thick, " * it("d") * " ≫ " * it("λ") * ":",
        15;
        anchor="start",
        color=GREY
    )
    text(fx + 14, 656, "→ " * mu0(it("λ")), 17; anchor="start")

    add("</svg>")
    mkpath(dirname(path))
    write(path, join(out, "\n") * "\n")
    return path
end
