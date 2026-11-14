# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generates docs/src/assets/examples/sic_reststrahlen_schematic.svg, the schematic of the
# silicon carbide half-space in docs/src/examples/dielectric_grating.md, during docs builds:
# a section of the periodic unit cell along the port normal, with the two Floquet ports, the
# waves, and the lengths of examples/dielectric_grating/mesh/sic_halfspace.msh, to scale
# along the normal. Port 2 lies inside the silicon carbide, so the silicon carbide is drawn
# beyond it, fading out, as the half-space which the port terminates.
#
# Sub- and superscripts use dy offsets, since Firefox ignores baseline-shift, and each reset
# is attached to the text that follows, since Chrome drops dy on an empty tspan, as in
# generate_coaxial_schematic.jl.

using Printf

const SIC_SCHEMATIC_PATH =
    joinpath(@__DIR__, "src", "assets", "examples", "sic_reststrahlen_schematic.svg")

function generate_sic_schematic(; path::String=SIC_SCHEMATIC_PATH)
    W, H = 900, 400
    RED, INK, GREY = "#c0392b", "#1a1a1a", "#666666"
    SERIF = "\"STIX Two Text\", \"Latin Modern Roman\", \"Times New Roman\", Times, serif"
    SIC = "url(#sic)"
    out = String[]
    add(s) = push!(out, s)
    f1(x) = @sprintf("%.1f", x)

    # Primitives
    attrs(kw) = join(("$(replace(string(k), "_" => "-"))='$v'" for (k, v) in kw), " ")
    it(s) = "<tspan font-style='italic'>$s</tspan>"
    # A space after a sub- or superscript is a dx offset: renderers drop a leading space
    # in a tspan.
    function sub(s, then; italic=false)
        style = italic ? " font-style='italic'" : ""
        dx = startswith(then, " ") ? " dx='4'" : ""
        return "<tspan dy='5' font-size='0.7em'$style>$s</tspan>" *
               "<tspan dy='-5'$dx>$(lstrip(then))</tspan>"
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
    arrow(x1, y1, x2, y2, color=INK, width=1.5, head="ink"; kw...) =
        line(x1, y1, x2, y2, color, width; marker_end="url(#head-$head)", kw...)
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
    function rect(x, y, w, h, fill; stroke=INK, width=1.4, kw...)
        return add(
            "<rect x='$(f1(x))' y='$(f1(y))' width='$(f1(w))' height='$(f1(h))' " *
            "fill='$fill' stroke='$stroke' stroke-width='$width' $(attrs(kw))/>"
        )
    end
    ε_r(then) = it("ε") * sub("r", then)

    # Document
    add(
        "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 $W $H' width='$W' " *
        "height='$H' font-family='$SERIF' role='img' aria-labelledby='title desc'>"
    )
    add("<title id='title'>Silicon carbide half-space between Floquet ports</title>")
    add(
        "<desc id='desc'>Section of the periodic unit cell along the normal of the ports: " *
        "Floquet port 1 in vacuum half a micrometer above the surface of the silicon " *
        "carbide, and Floquet port 2 a quarter of a micrometer below it, inside the " *
        "silicon carbide, which continues beyond port 2 as a half-space. A plane wave is " *
        "incident on the surface at normal incidence, and is reflected to port 1 and " *
        "transmitted to port 2.</desc>"
    )
    add("<defs>")
    for (name, color) in (("ink", INK), ("red", RED))
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
        "<pattern id='sic' width='7' height='7' patternUnits='userSpaceOnUse' " *
        "patternTransform='rotate(45)'><rect width='7' height='7' fill='#f3ece2'/>" *
        "<line x1='0' y1='0' x2='0' y2='7' stroke='#8c7b67' stroke-width='1.0'/></pattern>"
    )
    add(
        "<linearGradient id='fade' x1='0' y1='0' x2='1' y2='0'>" *
        "<stop offset='0' stop-color='white' stop-opacity='0'/>" *
        "<stop offset='1' stop-color='white' stop-opacity='1'/></linearGradient>"
    )
    add("</defs>")
    add("<rect width='$W' height='$H' fill='white'/>")

    # Geometry along the port normal z, to scale: 600 px per μm.
    scale = 600.0
    xp1 = 140.0                    # Port 1, z = -0.5 μm
    xs = xp1 + 0.5 * scale         # Surface, z = 0
    xp2 = xs + 0.25 * scale        # Port 2, z = 0.25 μm
    xend = xp2 + 170.0             # The half-space, drawn fading out
    ytop, ybot = 80.0, 300.0       # One period across the cell
    ymid = 0.5 * (ytop + ybot)

    # Silicon carbide: the meshed layer, then the half-space beyond port 2.
    rect(xs, ytop, xend - xs, ybot - ytop, SIC; stroke="none")
    rect(xp2 + 40, ytop - 1, xend - xp2 - 39, ybot - ytop + 2, "url(#fade)"; stroke="none")
    line(xs, ytop, xs, ybot, INK, 1.6)  # Surface
    # Periodic boundaries of the cell, from port 1 to port 2
    for y in (ytop, ybot)
        line(xp1, y, xp2, y, GREY, 1.0; stroke_dasharray="12 4 3 4")
    end
    # Floquet ports
    for (x, label) in ((xp1, "port 1"), (xp2, "port 2"))
        line(x, ytop - 18, x, ybot + 18, INK, 2.6)
        text(x, ytop - 28, "Floquet " * label, 17)
    end

    # Waves at normal incidence
    arrow(xp1 + 30, ymid - 44, xs - 30, ymid - 44, RED, 2.2, "red")
    text(0.5 * (xp1 + xs), ymid - 56, "incident", 16; color=RED)
    arrow(xs - 30, ymid + 44, xp1 + 30, ymid + 44, RED, 2.2, "red")
    text(
        0.5 * (xp1 + xs),
        ymid + 72,
        it("S") *
        sub("11", " = " * it("r e")) *
        "<tspan dy='-8' font-size='0.7em'>−2" *
        it("ik") *
        "</tspan><tspan dy='3' font-size='0.5em'>0</tspan>" *
        "<tspan dy='-3' font-size='0.7em'>" *
        it("d") *
        "</tspan><tspan dy='3' font-size='0.5em'>1</tspan>",
        17;
        color=RED
    )
    arrow(xs + 18, ymid, xp2 - 18, ymid, RED, 2.2, "red")
    text(
        0.5 * (xs + xp2),
        ymid - 12,
        it("S") * "<tspan dy='5' font-size='0.7em'>21</tspan>",
        17;
        color=RED
    )
    # The wave which port 2 absorbs, as the half-space beyond it would
    arrow(xp2 + 18, ymid, xp2 + 110, ymid, RED, 2.2, "red"; opacity="0.35")

    # Media
    text(0.5 * (xp1 + xs), ytop + 26, "vacuum, " * ε_r(" = 1"), 16; color=GREY)
    text(
        xp2 + 14,
        ytop + 26,
        "6H-SiC, " * ε_r("(" * it("iω") * ")"),
        16;
        anchor="start",
        color=GREY
    )
    text(xp2 + 14, ybot - 14, "half-space", 16; anchor="start", color=GREY)

    # Dimensions
    yd = ybot + 48
    for x in (xp1, xs, xp2)
        line(x, ybot + 22, x, yd + 8, INK, 0.8)
    end
    dim(xp1 + 1, yd, xs - 1, yd)
    text(0.5 * (xp1 + xs), yd + 24, it("d") * sub("1", " = 0.5 μm"), 17)
    dim(xs + 1, yd, xp2 - 1, yd)
    text(0.5 * (xs + xp2), yd + 24, it("d") * sub("2", " = 0.25 μm"), 17)
    xd = xp1 - 34
    dim(xd, ytop + 1, xd, ybot - 1)
    text(xd - 10, ymid + 6, it("a") * " = 1 μm", 17; anchor="end")
    text(xs, ytop - 10, it("z") * " = 0", 15; color=GREY)

    # Axes
    ox, oy = 800.0, 362.0
    arrow(ox, oy, ox + 60, oy, INK, 1.2)
    arrow(ox, oy, ox, oy - 60, INK, 1.2)
    text(ox + 70, oy + 6, it("z"), 17; anchor="start")
    text(ox, oy - 70, it("x"), 17)

    add("</svg>")
    mkpath(dirname(path))
    write(path, join(out, "\n") * "\n")
    return path
end
