# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generates docs/src/assets/examples/microstrip-1.svg, the schematic of the dispersive
# microstrip in docs/src/examples/microstrip.md, during docs builds. Panels: (a) the section
# along the line, with the two wave ports and the scattering parameters; (b) the
# cross-section at a port plane, with the field between the plates and around their edges;
# (c) the Djordjevic-Sarkar permittivity of the substrate, drawn from the exact model
# ε_r(f) = ε∞ + K log((f_hi + i f) / (f_lo + i f)) with the parameters of
# examples/microstrip/microstrip_driven.json.
#
# The air box of the mesh (10 mm × 44 mm × 40.2 mm) is far larger than the line, and the
# substrate (0.2 mm) far thinner than it is wide (4 mm), so (a) and (b) are not to scale.
# The field lines of (b) are a sketch: they leave the strip, at the higher potential, and
# wrap around the plate edges to the ground plane, which is as wide as the strip.
#
# Sub- and superscripts use dy offsets, since Firefox ignores baseline-shift, and each reset
# is attached to the text that follows, since Chrome drops dy on an empty tspan, as in
# generate_coaxial_schematic.jl.

using Printf

const MICROSTRIP_SCHEMATIC_PATH =
    joinpath(@__DIR__, "src", "assets", "examples", "microstrip-1.svg")

function generate_microstrip_schematic(; path::String=MICROSTRIP_SCHEMATIC_PATH)
    W, H = 1200, 740
    BLUE, RED, INK, GREY = "#1f5fa8", "#c0392b", "#1a1a1a", "#666666"
    SERIF = "\"STIX Two Text\", \"Latin Modern Roman\", \"Times New Roman\", Times, serif"
    PEC, SUBSTRATE, ABSORBING, BAND = "url(#pec)", "#eef3e6", "url(#absorbing)", "#f4e1de"
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
    # Air box of (a) and (b): its faces carry the absorbing boundary, drawn as a band
    # outside the domain.
    function air_box(x0, y0, x1, y1; sides=(:top, :bottom, :left, :right), band=8)
        rect(x0, y0, x1 - x0, y1 - y0, "white"; stroke="none")
        # The top and bottom bands cover the corners when the side bands are present.
        xa = :left in sides ? x0 - band : x0
        xb = :right in sides ? x1 + band : x1
        :top in sides && rect(xa, y0 - band, xb - xa, band, ABSORBING; width=1.0)
        :bottom in sides && rect(xa, y1, xb - xa, band, ABSORBING; width=1.0)
        :left in sides && rect(x0 - band, y0, band, y1 - y0, ABSORBING; width=1.0)
        :right in sides && rect(x1, y0, band, y1 - y0, ABSORBING; width=1.0)
        return
    end

    # Document
    add(
        "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 $W $H' width='$W' " *
        "height='$H' font-family='$SERIF' role='img' aria-labelledby='title desc'>"
    )
    add("<title id='title'>Microstrip line on a dispersive substrate</title>")
    add(
        "<desc id='desc'>Section along a PEC strip and ground plane of equal width on a " *
        "thin dispersive substrate, in an air box with absorbing faces and a wave port at " *
        "each end; the cross-section at a port plane with the electric field between the " *
        "plates and around their edges; and the real and imaginary parts of the " *
        "Djordjevic-Sarkar permittivity against frequency on a logarithmic axis.</desc>"
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
    add(
        "<pattern id='absorbing' width='5' height='5' patternUnits='userSpaceOnUse'>" *
        "<rect width='5' height='5' fill='#ececec'/>" *
        "<circle cx='2.5' cy='2.5' r='0.9' fill='#7a7a7a'/></pattern>"
    )
    add("</defs>")
    add("<rect width='$W' height='$H' fill='white'/>")

    # (a) Section along the line, through the middle of the strip
    panel_title(40, 36, "(a)", "Section along the line (not to scale)")
    x0, x1 = 200.0, 780.0         # Port planes x = -ℓ/2 and x = ℓ/2
    yt, yb, yc = 84.0, 316.0, 200.0  # Air box faces and the middle of the substrate
    hg, tp = 34.0, 7.0            # Exaggerated substrate and plate thickness (px)
    air_box(x0, yt, x1, yb; sides=(:top, :bottom))
    rect(x0, yc - hg / 2, x1 - x0, hg, SUBSTRATE; stroke="none")
    rect(x0, yc - hg / 2 - tp, x1 - x0, tp, PEC; width=1.2)  # Strip
    rect(x0, yc + hg / 2, x1 - x0, tp, PEC; width=1.2)       # Ground plane
    text(0.5 * (x0 + x1), yc + 5, "FR4, " * it("ε") * sub("r", "(") * it("ω") * ")", 15)
    text(0.5 * (x0 + x1), yc - hg / 2 - tp - 10, "strip", 15)
    text(0.5 * (x0 + x1), yc + hg / 2 + tp + 21, "ground plane", 15)
    text(x0 + 70, yt + 34, "air", 15; color=GREY)
    text(x1 - 12, yt - 16, "absorbing", 14; anchor="end", color=GREY)

    # Wave ports on the two end planes, which are the whole cross-section of the box. The
    # VoltagePath of each runs from the strip to the ground plane.
    for (x, k, side) in ((x0, 1, -1), (x1, 2, 1))
        line(x, yt - 8, x, yb + 8, BLUE, 2.6)
        text(x, yb + 32, "port $k", 16; color=BLUE)
        text(
            x,
            yb + 52,
            it("x") * (k == 1 ? " = −" : " = ") * it("ℓ") * "/2",
            15;
            color=GREY
        )
        xv = x - 16 * side
        arrow(xv, yc - hg / 2 + 1, xv, yc + hg / 2 - 1, INK, 1.4)
        text(xv - 9 * side, yc + 5, it("V"), 15; anchor=(side < 0 ? "start" : "end"))
    end
    # Incident, reflected and transmitted waves
    arrow(x0 - 120, yc - 13, x0 - 10, yc - 13, RED, 2.2, "red")
    text(x0 - 65, yc - 23, it("a") * sub("1", ""), 17; color=RED)
    arrow(x0 - 10, yc + 13, x0 - 120, yc + 13, RED, 2.2, "red")
    text(x0 - 65, yc + 38, it("S") * sub("11", "") * it("a") * sub("1", ""), 17; color=RED)
    arrow(x1 + 10, yc, x1 + 120, yc, RED, 2.2, "red")
    text(x1 + 65, yc - 12, it("S") * sub("21", "") * it("a") * sub("1", ""), 17; color=RED)
    # Length between the port planes
    yd = yb + 78
    dim(x0 + 1, yd, x1 - 1, yd)
    text(0.5 * (x0 + x1), yd - 8, it("ℓ"), 17)

    # Legend
    lx = 950.0
    for (k, (fill, label)) in enumerate((
        (PEC, "PEC"),
        (SUBSTRATE, "Djordjevic-Sarkar substrate"),
        (ABSORBING, "absorbing boundary"),
        ("", "wave port")
    ))
        ly = 120.0 + 26 * (k - 1)
        if isempty(fill)
            line(lx, ly - 4.5, lx + 26, ly - 4.5, BLUE, 2.6)
        else
            rect(lx, ly - 11, 26, 13, fill; stroke=INK, width=1.0)
        end
        text(lx + 36, ly, label, 15; anchor="start")
    end

    # (b) Cross-section at a port plane
    panel_title(40, 436, "(b)", "Cross-section at a port plane (not to scale)")
    bx0, bx1, by0, by1 = 90.0, 470.0, 470.0, 712.0  # Air box
    cx, cy = 0.5 * (bx0 + bx1), 590.0                # Middle of the substrate
    w = 190.0                                        # Plate width (px)
    xl, xr = cx - w / 2, cx + w / 2
    ys, yg = cy - hg / 2, cy + hg / 2                # Substrate faces
    air_box(bx0, by0, bx1, by1)
    rect(xl, ys, w, hg, SUBSTRATE; stroke="none")
    # Field between the plates, uniform away from the edges, from the strip to the ground
    for x in (xl + 22, xl + 52, xr - 52, xr - 22)
        arrow(x, ys + 2, x, yg - 2, RED, 1.8, "red")
    end
    # Fringing field: arcs of concentric ellipses about each plate edge, from the back of
    # the strip around the edge to the back of the ground plane, with the arrow where they
    # cross the plane of the substrate
    hb = hg / 2 + tp  # Half the distance between the plate backs
    for (xe, side) in ((xl, -1), (xr, 1)),
        (rx, ry) in ((16.0, 31.0), (30.0, 45.0), (46.0, 61.0))

        θ0 = π - asin(hb / ry)
        pts =
            [(xe + side * rx * cos(θ), cy - ry * sin(θ)) for θ in range(θ0, -θ0; length=80)]
        polyline(pts, RED, 1.3)
        arrow(xe + side * rx, cy - 3, xe + side * rx, cy + 3, RED, 1.3, "red")
    end
    rect(xl, ys - tp, w, tp, PEC; width=1.2)  # Strip
    rect(xl, yg, w, tp, PEC; width=1.2)       # Ground plane
    text(xl + 30, cy + 5, it("E"), 16; color=RED, anchor="start")
    dim(cx, ys + 1, cx, yg - 1)
    text(cx + 8, cy + 6, it("h"), 17; anchor="start")
    yw = cy - 82
    dim(xl + 1, yw, xr - 1, yw)
    text(cx, yw - 8, it("w"), 17)
    text(cx, ys - tp - 8, "strip", 14; color=GREY)
    text(cx, yg + tp + 19, "ground plane", 14; color=GREY)
    text(bx0 + 14, by0 + 26, "air", 15; anchor="start", color=GREY)
    text(bx1 - 12, by1 - 14, "absorbing on all four faces", 14; anchor="end", color=GREY)

    # (c) Djordjevic-Sarkar permittivity on a logarithmic frequency axis
    panel_title(560, 436, "(c)", "Djordjevic-Sarkar substrate")
    ε_inf, K, f_lo, f_hi = 4.0, 0.056, 1.0e-3, 1.0e3  # GHz, as in the example inputs
    ε(f) = ε_inf + K * log(complex(f_hi, f) / complex(f_lo, f))
    gx0, gx1 = 660.0, 1120.0
    lmin, lmax = -6.0, 5.0  # log10(f / GHz): 1 kHz to 100 THz
    X(f) = gx0 + (log10(f) - lmin) / (lmax - lmin) * (gx1 - gx0)
    # Upper graph Re ε_r, lower graph -Im ε_r
    u0, u1, re_min, re_max = 488.0, 588.0, 3.92, 4.86
    v0, v1, im_max = 612.0, 690.0, 0.105
    Yre(e) = u1 - (e - re_min) / (re_max - re_min) * (u1 - u0)
    Yim(e) = v1 - e / im_max * (v1 - v0)

    # Sweep of the example, 2-14 GHz, and the corner frequencies
    rect(X(2.0), u0 - 8, X(14.0) - X(2.0), v1 - u0 + 8, BAND; stroke="none")
    text(0.5 * (X(2.0) + X(14.0)), u0 - 14, "2–14 GHz", 14; color=RED)
    for f in (f_lo, f_hi)
        line(X(f), u0 - 8, X(f), v1, GREY, 1.0; stroke_dasharray="3 3")
    end
    # Axes
    for (y0, y1) in ((u0 - 14, u1), (v0 - 14, v1))
        arrow(gx0, y1, gx0, y0, INK, 1.2)
        line(gx0, y1, gx1, y1, INK, 1.2)
    end
    arrow(gx1 - 20, v1, gx1 + 14, v1, INK, 1.2)
    text(gx1 + 18, v1 + 5, it("f"), 17; anchor="start")
    for l = (lmin + 1):(lmax - 1)
        line(X(10.0^l), v1, X(10.0^l), v1 + 5, INK, 1.0)
    end
    text(X(f_lo), v1 + 24, it("f") * sub("lo", "\u00a0= 1 MHz"), 15)
    text(X(1.0), v1 + 24, "1 GHz", 15)
    text(X(f_hi), v1 + 24, it("f") * sub("hi", "\u00a0= 1 THz"), 15)
    text(
        gx0 - 8,
        Yre(4.4) + 6,
        "Re " * it("ε") * sub("r", ""; italic=false),
        16;
        anchor="end"
    )
    text(gx0 - 8, Yim(0.05) + 6, "−Im " * it("ε") * sub("r", ""), 16; anchor="end")

    # The model: Re ε_r falls by K ln 10 per decade between the corners, and -Im ε_r is
    # close to π K / 2 there.
    fs = exp10.(range(lmin, lmax; length=400))
    polyline([(X(f), Yre(real(ε(f)))) for f in fs], BLUE, 2.4)
    polyline([(X(f), Yim(-imag(ε(f)))) for f in fs], BLUE, 2.4)
    text(
        X(2.0e-6),
        Yre(ε_inf + K * log(f_hi / f_lo)) - 9,
        it("ε") *
        sub("∞", "\u00a0+ ") *
        it("K") *
        " ln(" *
        it("f") *
        sub("hi", "\u00a0/ ") *
        it("f") *
        sub("lo", ")"),
        14;
        anchor="start"
    )
    text(X(2.0e3), Yre(ε_inf) - 9, it("ε") * sub("∞", ""), 15; anchor="start")
    text(X(0.5), Yre(4.32), "slope −" * it("K") * " ln 10", 14; anchor="end", color=BLUE)
    text(X(0.5), Yre(4.32) + 17, "per decade", 14; anchor="end", color=BLUE)
    text(X(3.0e-2), Yim(π * K / 2) - 9, "≈ π" * it("K") * "/2", 15; color=BLUE)

    add("</svg>")
    mkpath(dirname(path))
    write(path, join(out, "\n") * "\n")
    return path
end
