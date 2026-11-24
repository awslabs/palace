# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generates docs/src/assets/examples/pml-extrusion.svg, the diagrams of the extruded PML in
# docs/src/guide/boundaries.md, during docs builds. Each panel is a cross-section of a
# physical region (blue) with the layers of the PML extruded from its boundary (orange),
# the center of the rays (red) and a few arrows along the extrusion, with the entry of
# config["Boundaries"]["PML"] under it.
#
# The layers are drawn by moving the points of the boundary to the depths k t / N, as the mesh
# extrusion does with the nodes, so the shapes follow from the rule for the directions: along
# the normal of a flat part; along the ray d = (x - p) / |x - p| from the closest point p of
# the center on a curved part; and on a box, component by component, by t (x_i - p_i) / g_i
# with g_i the distance from p_i to the face of the box on that side, so that the outside of
# the layer is the box grown by t. The colors and the gray text read on the light and the
# dark theme of the documentation.

using Printf

const PML_DIAGRAMS_PATH =
    joinpath(@__DIR__, "src", "assets", "examples", "pml-extrusion.svg")

function generate_pml_diagrams(; path::String=PML_DIAGRAMS_PATH)
    PW, PH = 300, 330             # Panel size
    NX, NY = 4, 2                 # Panels per row and column
    W, H = NX * PW, NY * PH + 50  # With a legend row
    S = 55.0                      # Pixels per unit length in the panels
    GRAY, PHYS, PML, CORE = "#8c8c8c", "#4a90d9", "#e0892e", "#d6455d"
    SUB = "#9a7fd1"
    SERIF = "\"STIX Two Text\", \"Latin Modern Roman\", \"Times New Roman\", Times, serif"
    MONO = "Menlo, Consolas, \"DejaVu Sans Mono\", monospace"
    out = String[]
    add(s) = push!(out, s)
    f1(x) = @sprintf("%.1f", x)

    # Primitives, in pixels
    attrs(kw) = join(("$(replace(string(k), "_" => "-"))='$v'" for (k, v) in kw), " ")
    function text(x, y, body, size=15; anchor="middle", color=GRAY, kw...)
        return add(
            "<text x='$(f1(x))' y='$(f1(y))' font-size='$size' text-anchor='$anchor' " *
            "fill='$color' $(attrs(kw))>$body</text>"
        )
    end
    function line(x1, y1, x2, y2, color=GRAY, width=1.5; kw...)
        return add(
            "<line x1='$(f1(x1))' y1='$(f1(y1))' x2='$(f1(x2))' y2='$(f1(y2))' " *
            "stroke='$color' stroke-width='$width' $(attrs(kw))/>"
        )
    end
    pts(p) = join(("$(f1(x)),$(f1(y))" for (x, y) in p), " ")
    polyline(p, color, width=1.0; kw...) = add(
        "<polyline points='$(pts(p))' fill='none' stroke='$color' stroke-width='$width' " *
        "stroke-linejoin='round' $(attrs(kw))/>"
    )
    polygon(p, fill, opacity; kw...) = add(
        "<polygon points='$(pts(p))' fill='$fill' fill-opacity='$opacity' $(attrs(kw))/>"
    )
    # A region between an inner and an outer closed curve.
    function ring(inner, outer, fill, opacity)
        d = "M " * pts(outer) * " Z M " * pts(reverse(inner)) * " Z"
        return add(
            "<path d='$(replace(d, "," => " "))' fill='$fill' fill-opacity='$opacity' " *
            "fill-rule='evenodd'/>"
        )
    end
    dot(x, y, r=4.0, color=CORE) =
        add("<circle cx='$(f1(x))' cy='$(f1(y))' r='$(f1(r))' fill='$color'/>")
    # An arrow with its head drawn as a triangle, which every SVG renderer shows.
    function arrow(x1, y1, x2, y2, color=CORE, width=1.4; head=6.0)
        l = hypot(x2 - x1, y2 - y1)
        ux, uy = (x2 - x1) / l, (y2 - y1) / l
        bx, by = x2 - head * ux, y2 - head * uy
        line(x1, y1, bx, by, color, width)
        tip = [
            (x2, y2),
            (bx - 0.45 * head * uy, by + 0.45 * head * ux),
            (bx + 0.45 * head * uy, by - 0.45 * head * ux)
        ]
        return polygon(tip, color, 1.0)
    end

    # Panels, in units of S with v up, centered on the drawing area of panel (i, j).
    function panel(i, j, title, config, draw)
        x0, y0 = (i - 1) * PW, (j - 1) * PH
        cx, cy = x0 + PW / 2, y0 + 140
        P(u, v) = (cx + S * u, cy - S * v)
        text(x0 + 12, y0 + 24, title, 16; anchor="start", font_family=SERIF)
        draw(P)
        for (k, l) in enumerate(config)
            text(
                x0 + 14,
                y0 + 262 + 15 * (k - 1),
                replace(l, " " => "&#160;"),
                11.5;
                anchor="start",
                font_family=MONO
            )
        end
    end

    # The layers of a closed boundary, given as points with their directions.
    function layers(P, boundary, dirs, t, N; closed=true)
        curves = [
            [P((b .+ s * t .* d)...) for (b, d) in zip(boundary, dirs)] for s in (0:N) ./ N
        ]
        if closed
            ring(curves[1], curves[end], PML, 0.28)
        else
            polygon(vcat(curves[1], reverse(curves[end])), PML, 0.28)
        end
        for c in curves[2:(end - 1)]
            polyline(closed ? vcat(c, [c[1]]) : c, PML, 1.0)
        end
        c = curves[end]
        polyline(closed ? vcat(c, [c[1]]) : c, GRAY, 1.6)  # The back of the layer
        return curves
    end
    unit(v) = v ./ hypot(v...)
    closest_box(x, lo, hi) = clamp.(x, lo, hi)
    # The move per unit thickness of a point b of the box [lo, hi] for the center point p.
    box_dir(b, p, lo, hi) = Tuple(
        b[i] > p[i] ? (b[i] - p[i]) / (hi[i] - p[i]) :
        (b[i] < p[i] ? (b[i] - p[i]) / (p[i] - lo[i]) : 0.0) for i = 1:2
    )
    function square(a, n)  # Points along the boundary of the square |u|, |v| ≤ a
        s = range(-a, a; length=n)
        return vcat(
            [(a, v) for v in s[1:(end - 1)]],
            [(-u, a) for u in s[1:(end - 1)]],
            [(-a, -v) for v in s[1:(end - 1)]],
            [(u, -a) for u in s[1:(end - 1)]]
        )
    end
    # Rays from the points b of the boundary: dashed from the center to b, an arrow across
    # the layer.
    function rays(P, starts, boundary, dirs, t)
        for (p, b, d) in zip(starts, boundary, dirs)
            if p !== nothing && hypot((b .- p)...) > 1.0e-9
                line(
                    P(p...)...,
                    P(b...)...,
                    CORE,
                    1.0;
                    stroke_dasharray="3 3",
                    stroke_opacity="0.7"
                )
            end
            arrow(P(b...)..., P((b .+ 0.9 * t .* d)...)...)
        end
    end

    # Document
    add(
        "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 $W $H' width='$W' " *
        "height='$H' font-family='$SERIF' role='img' aria-labelledby='title desc'>"
    )
    add("<title id='title'>Extrusion of the PML</title>")
    add(
        "<desc id='desc'>Cross-sections of physical regions with the layers of a PML " *
        "extruded from their boundaries, along the normals of flat parts or along rays " *
        "from a point, a segment or a box, with the configuration of each.</desc>"
    )
    # (a) Sphere: rays from its center.
    panel(
        1,
        1,
        "(a) Sphere, radial from its center",
        ["{ \"Attributes\": [4], \"Thickness\": 1.0,", "  \"Layers\": 4 }"],
        P -> begin
            R, t, N = 1.25, 0.6, 4
            θ = range(0, 2π; length=181)[1:(end - 1)]
            b = [(R * cos(x), R * sin(x)) for x in θ]
            polygon([P(p...) for p in b], PHYS, 0.16; stroke=PHYS, stroke_width="1.2")
            layers(P, b, [unit(p) for p in b], t, N)
            ϕ = (0:7) .* (π / 4) .+ π / 8
            bb = [(R * cos(x), R * sin(x)) for x in ϕ]
            rays(P, fill((0.0, 0.0), 8), bb, [unit(p) for p in bb], t)
            dot(P(0, 0)...)
        end
    )

    # (b) Waveguide end: flat, along its normal; the walls continue through the layer.
    function guide(P; left::Bool)
        L, h, t, N = 2.0, 0.7, 0.7, 6
        xl, xr = left ? -L + t : -L, L - t
        polygon([P(xl, -h), P(xr, -h), P(xr, h), P(xl, h)], PHYS, 0.16)
        for (x0, s) in (left ? ((xl, -1), (xr, 1)) : ((xr, 1),))
            polygon([P(x0, -h), P(x0 + s * t, -h), P(x0 + s * t, h), P(x0, h)], PML, 0.28)
            for k = 1:N
                line(P(x0 + s * k * t / N, -h)..., P(x0 + s * k * t / N, h)..., PML, 1.0)
            end
            line(P(x0, -h)..., P(x0, h)..., PML, 1.6)
            line(P(x0 + s * t, -h)..., P(x0 + s * t, h)..., GRAY, 3.0)
            for v in (-0.45, 0.0, 0.45)
                arrow(P(x0, v)..., P(x0 + s * 0.9 * t, v)...)
            end
        end
        xa, xb = left ? (xl - t, xr + t) : (xl, xr + t)
        line(P(xa, h)..., P(xb, h)..., GRAY, 3.0)
        return line(P(xa, -h)..., P(xb, -h)..., GRAY, 3.0)
    end
    panel(
        2,
        1,
        "(b) Waveguide end, along its normal",
        ["{ \"Attributes\": [3], \"Thickness\": 0.1,", "  \"Layers\": 6 }"],
        P -> guide(P; left=false)
    )

    # (c) Pipe with both ends: two flat parts, each along its own normal.
    panel(
        3,
        1,
        "(c) Pipe, both ends",
        ["{ \"Attributes\": [5, 6], \"Thickness\": 0.1,", "  \"Layers\": 6 }"],
        P -> guide(P; left=true)
    )

    # (d) Cube: radial from its center, one closed layer whose outside is a larger cube.
    panel(
        4,
        1,
        "(d) Cube, radial from its center",
        ["{ \"Attributes\": [7], \"Thickness\": 0.5,", "  \"Layers\": 4 }"],
        P -> begin
            a, t, N = 1.0, 0.5, 4
            b = square(a, 61)
            polygon([P(p...) for p in b], PHYS, 0.16; stroke=PHYS, stroke_width="1.2")
            dir(p) = box_dir(p, (0.0, 0.0), (-a, -a), (a, a))
            layers(P, b, dir.(b), t, N)
            bb = [
                (a, 0.0),
                (a, a),
                (0.0, a),
                (-a, a),
                (-a, 0.0),
                (-a, -a),
                (0.0, -a),
                (a, -a),
                (a, 0.5a),
                (0.5a, a),
                (-0.5a, a),
                (-a, 0.5a),
                (-a, -0.5a),
                (-0.5a, -a),
                (0.5a, -a),
                (a, -0.5a)
            ]
            rays(P, fill((0.0, 0.0), length(bb)), bb, dir.(bb), t)
            dot(P(0, 0)...)
        end
    )

    # (e) Cube: radial from a box inside it, the same outside, with straight columns over the
    # box and fans at the edges.
    panel(
        1.5,
        2,
        "(e) Cube, radial from a box",
        [
            "{ \"Attributes\": [7], \"Thickness\": 0.5,",
            "  \"Layers\": 4,",
            "  \"Center\": [[-0.7, -0.7, -0.7],",
            "             [0.7, 0.7, 0.7]] }"
        ],
        P -> begin
            a, c, t, N = 1.0, 0.7, 0.5, 4
            b = square(a, 121)
            dir(p) = box_dir(p, closest_box(p, (-c, -c), (c, c)), (-a, -a), (a, a))
            polygon([P(p...) for p in b], PHYS, 0.16; stroke=PHYS, stroke_width="1.2")
            layers(P, b, dir.(b), t, N)
            polygon(
                [P(-c, -c), P(c, -c), P(c, c), P(-c, c)],
                CORE,
                0.12;
                stroke=CORE,
                stroke_width="1.4"
            )
            bb = [
                (a, 0.0),
                (0.0, a),
                (-a, 0.0),
                (0.0, -a),
                (a, 0.4),
                (-0.4, a),
                (-a, -0.4),
                (0.4, -a)
            ]
            for (sx, sy) in ((1, 1), (-1, 1), (-1, -1), (1, -1))
                append!(bb, [(sx * a, sy * 0.85), (sx * a, sy * a), (sx * 0.85, sy * a)])
            end
            rays(P, [closest_box(p, (-c, -c), (c, c)) for p in bb], bb, dir.(bb), t)
        end
    )

    # (f) Floquet cell with a PML on top; the periodic sides continue through the layer.
    panel(
        2.5,
        2,
        "(f) Floquet cell, PML on top",
        ["{ \"Attributes\": [6], \"Thickness\": 0.1,", "  \"Layers\": 6 }"],
        P -> begin
            w, h, t, N = 0.75, 1.1, 0.8, 6
            v0 = -1.0
            polygon([P(-w, v0), P(w, v0), P(w, v0 + h), P(-w, v0 + h)], PHYS, 0.16)
            polygon(
                [P(-w, v0 + h), P(w, v0 + h), P(w, v0 + h + t), P(-w, v0 + h + t)],
                PML,
                0.28
            )
            for k = 1:N
                line(P(-w, v0 + h + k * t / N)..., P(w, v0 + h + k * t / N)..., PML, 1.0)
            end
            line(P(-w, v0 + h)..., P(w, v0 + h)..., PML, 1.6)
            line(P(-w, v0 + h + t)..., P(w, v0 + h + t)..., GRAY, 3.0)
            for u in (-w, w)
                line(P(u, v0)..., P(u, v0 + h + t)..., GRAY, 1.6; stroke_dasharray="6 4")
            end
            line(P(-w, v0)..., P(w, v0)..., PHYS, 3.0)
            text(P(0, v0)[1], P(0, v0)[2] + 18, "Floquet port", 13)
            for u in (-w, w)
                x, y = P(u + (u < 0 ? -0.12 : 0.26), v0 + h / 2 + t / 2)
                text(x, y, "periodic", 13; transform="rotate(-90 $(f1(x)) $(f1(y)))")
            end
            for u in (-0.4, 0.0, 0.4)
                arrow(P(u, v0 + h)..., P(u, v0 + h + 0.9 * t)...)
            end
        end
    )

    # (g) Patch antenna: one closed layer on the sides and the top, radial from a segment
    # through the substrate on the ground plane, so that the substrate stays flat through
    # the side layers and the layer ends flush with the ground plane.
    panel(
        3.5,
        2,
        "(g) Patch antenna, radial from a segment",
        [
            "{ \"Attributes\": [7], \"Thickness\": 4.0,",
            "  \"Layers\": 4,",
            "  \"Center\": [[0, 0, 0], [0, 0, 0.794]] }"
        ],
        P -> begin
            a, zt, hs, t, N = 1.2, 1.3, 0.3, 0.5, 4
            v0 = -0.8
            seg(p) = (0.0, clamp(p[2], 0.0, hs))
            dir(p) = box_dir(p, seg(p), (-a, 0.0), (a, zt))
            z = range(0.0, zt; length=60)
            u = range(a, -a; length=81)
            b = vcat(
                [(a, v) for v in z],
                [(x, zt) for x in u[2:end]],
                [(-a, v) for v in reverse(z)[2:end]]
            )
            Q(p) = P(p[1], p[2] + v0)
            polygon([Q(p) for p in vcat(b, [(-a, 0.0), (a, 0.0)])], PHYS, 0.16)
            polygon([Q((-a, 0.0)), Q((a, 0.0)), Q((a, hs)), Q((-a, hs))], SUB, 0.35)
            curves = [[Q(p .+ s * t .* dir(p)) for p in b] for s in (0:N) ./ N]
            polygon(vcat(curves[1], reverse(curves[end])), PML, 0.28)
            for (x0, s) in ((a, 1), (-a, -1))  # The substrate continues through the layer
                polygon(
                    [Q((x0, 0.0)), Q((x0 + s * t, 0.0)), Q((x0 + s * t, hs)), Q((x0, hs))],
                    SUB,
                    0.35
                )
            end
            for c in curves[2:(end - 1)]
                polyline(c, PML, 1.0)
            end
            polyline(curves[end], GRAY, 1.6)
            line(Q((-a - t, 0.0))..., Q((a + t, 0.0))..., GRAY, 3.0)
            line(Q((-0.55, hs))..., Q((0.55, hs))..., GRAY, 3.0)
            text(Q((0.0, 0.0))[1], Q((0.0, 0.0))[2] + 18, "ground plane", 13)
            text(Q((0.6, hs))[1], Q((0.6, hs))[2] - 6, "patch", 13; anchor="start")
            bb = [
                (a, 0.15),
                (a, 0.7),
                (a, zt),
                (0.6, zt),
                (0.0, zt),
                (-0.6, zt),
                (-a, zt),
                (-a, 0.7),
                (-a, 0.15)
            ]
            rays(
                P,
                [seg(p) .+ (0.0, v0) for p in bb],
                [p .+ (0.0, v0) for p in bb],
                dir.(bb),
                t
            )
            line(Q((0.0, 0.0))..., Q((0.0, hs))..., CORE, 3.0)
        end
    )

    # Legend
    y = NY * PH + 25
    items = [
        (
            "Physical region",
            x -> polygon(
                [(x, y - 7), (x + 22, y - 7), (x + 22, y + 7), (x, y + 7)],
                PHYS,
                0.3
            )
        ),
        (
            "PML layers",
            x -> polygon(
                [(x, y - 7), (x + 22, y - 7), (x + 22, y + 7), (x, y + 7)],
                PML,
                0.45
            )
        ),
        (
            "Substrate",
            x -> polygon(
                [(x, y - 7), (x + 22, y - 7), (x + 22, y + 7), (x, y + 7)],
                SUB,
                0.5
            )
        ),
        ("Center", x -> dot(x + 11, y, 4.5)),
        ("Extrusion", x -> arrow(x, y, x + 24, y)),
        ("PEC", x -> line(x, y, x + 24, y, GRAY, 3.0))
    ]
    x = 110
    for (label, draw) in items
        draw(x)
        text(x + 32, y + 5, label, 15; anchor="start")
        x += 52 + 9 * length(label)
    end
    add("</svg>")

    mkpath(dirname(path))
    open(path, "w") do f
        return write(f, join(out, "\n") * "\n")
    end
    return path
end
