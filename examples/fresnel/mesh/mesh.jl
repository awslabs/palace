# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate mesh with:
# julia -e 'include("mesh/mesh.jl"); generate_fresnel_mesh(filename="mesh/fresnel.msh")'

using Gmsh: gmsh

"""
    generate_fresnel_mesh(;
        filename::AbstractString = "fresnel.msh",
        a::Real = 1.0,
        L::Real = 5.0,
        n_xy::Integer = 2,
        n_z::Integer = 5,
        order::Integer = 1,
        verbose::Integer = 5,
        gui::Bool = false
    )

Generate a structured hex mesh for the Fresnel reflection test case.

The geometry is a box [0,a] x [0,a] x [-L,L] with:

  - Two material regions: vacuum (z < 0, attr 1) and dielectric (z > 0, attr 2)
  - Periodic BCs on x-faces (attr 1,2) and y-faces (attr 3,4)
  - Floquet port at z = -L (attr 5, excitation) and z = L (attr 6, transmission)

Uses transfinite structured hex mesh to avoid OCC fragmentation issues.
Units: cm (Palace L0 = 1e-2 for cm -> m conversion).
"""
function generate_fresnel_mesh(;
    filename::AbstractString="fresnel.msh",
    a::Real=1.0,
    L::Real=5.0,
    n_xy::Integer=2,
    n_z::Integer=5,
    order::Integer=1,
    verbose::Integer=5,
    gui::Bool=false
)
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    if "fresnel" in gmsh.model.list()
        gmsh.model.setCurrent("fresnel")
        gmsh.model.remove()
    end
    gmsh.model.add("fresnel")

    geo = gmsh.model.geo

    # Define points for the two stacked boxes.
    # Bottom box: [0,a] x [0,a] x [-L, 0] (vacuum)
    # Top box:    [0,a] x [0,a] x [0, L]  (dielectric)

    # Layer z = -L (bottom port)
    p1 = geo.addPoint(0.0, 0.0, -L)
    p2 = geo.addPoint(a, 0.0, -L)
    p3 = geo.addPoint(a, a, -L)
    p4 = geo.addPoint(0.0, a, -L)

    # Layer z = 0 (interface)
    p5 = geo.addPoint(0.0, 0.0, 0.0)
    p6 = geo.addPoint(a, 0.0, 0.0)
    p7 = geo.addPoint(a, a, 0.0)
    p8 = geo.addPoint(0.0, a, 0.0)

    # Layer z = L (top port)
    p9  = geo.addPoint(0.0, 0.0, L)
    p10 = geo.addPoint(a, 0.0, L)
    p11 = geo.addPoint(a, a, L)
    p12 = geo.addPoint(0.0, a, L)

    # Lines — bottom face (z = -L)
    l1 = geo.addLine(p1, p2)
    l2 = geo.addLine(p2, p3)
    l3 = geo.addLine(p3, p4)
    l4 = geo.addLine(p4, p1)

    # Lines — interface (z = 0)
    l5 = geo.addLine(p5, p6)
    l6 = geo.addLine(p6, p7)
    l7 = geo.addLine(p7, p8)
    l8 = geo.addLine(p8, p5)

    # Lines — top face (z = L)
    l9  = geo.addLine(p9, p10)
    l10 = geo.addLine(p10, p11)
    l11 = geo.addLine(p11, p12)
    l12 = geo.addLine(p12, p9)

    # Vertical lines — bottom to interface
    l13 = geo.addLine(p1, p5)
    l14 = geo.addLine(p2, p6)
    l15 = geo.addLine(p3, p7)
    l16 = geo.addLine(p4, p8)

    # Vertical lines — interface to top
    l17 = geo.addLine(p5, p9)
    l18 = geo.addLine(p6, p10)
    l19 = geo.addLine(p7, p11)
    l20 = geo.addLine(p8, p12)

    # Curve loops and surfaces for bottom box (vacuum).
    # Bottom face (z = -L)
    cl_bot = geo.addCurveLoop([l1, l2, l3, l4])
    s_bot = geo.addPlaneSurface([cl_bot])

    # Interface face (z = 0)
    cl_int = geo.addCurveLoop([l5, l6, l7, l8])
    s_int = geo.addPlaneSurface([cl_int])

    # Top face (z = L)
    cl_top = geo.addCurveLoop([l9, l10, l11, l12])
    s_top = geo.addPlaneSurface([cl_top])

    # Side faces of bottom box.
    cl_s1 = geo.addCurveLoop([l1, l14, -l5, -l13])  # x-z face at y=0, bottom
    s_s1 = geo.addPlaneSurface([cl_s1])
    cl_s2 = geo.addCurveLoop([l2, l15, -l6, -l14])  # y-z face at x=a, bottom
    s_s2 = geo.addPlaneSurface([cl_s2])
    cl_s3 = geo.addCurveLoop([l3, l16, -l7, -l15])  # x-z face at y=a, bottom
    s_s3 = geo.addPlaneSurface([cl_s3])
    cl_s4 = geo.addCurveLoop([l4, l13, -l8, -l16])  # y-z face at x=0, bottom
    s_s4 = geo.addPlaneSurface([cl_s4])

    # Side faces of top box.
    cl_s5 = geo.addCurveLoop([l5, l18, -l9, -l17])  # x-z face at y=0, top
    s_s5 = geo.addPlaneSurface([cl_s5])
    cl_s6 = geo.addCurveLoop([l6, l19, -l10, -l18])  # y-z face at x=a, top
    s_s6 = geo.addPlaneSurface([cl_s6])
    cl_s7 = geo.addCurveLoop([l7, l20, -l11, -l19])  # x-z face at y=a, top
    s_s7 = geo.addPlaneSurface([cl_s7])
    cl_s8 = geo.addCurveLoop([l8, l17, -l12, -l20])  # y-z face at x=0, top
    s_s8 = geo.addPlaneSurface([cl_s8])

    # Surface loops and volumes.
    sl_vac = geo.addSurfaceLoop([s_bot, s_int, s_s1, s_s2, s_s3, s_s4])
    v_vac = geo.addVolume([sl_vac])

    sl_die = geo.addSurfaceLoop([s_int, s_top, s_s5, s_s6, s_s7, s_s8])
    v_die = geo.addVolume([sl_die])

    # Physical groups.
    geo.synchronize()

    gmsh.model.addPhysicalGroup(3, [v_vac], 1, "vacuum")
    gmsh.model.addPhysicalGroup(3, [v_die], 2, "dielectric")

    # Boundary attributes:
    # 1 = x_min (x=0), 2 = x_max (x=a)  [periodic pair in x]
    # 3 = y_min (y=0), 4 = y_max (y=a)  [periodic pair in y]
    # 5 = z_min (z=-L, port 1 excitation)
    # 6 = z_max (z=L, port 2 transmission)
    gmsh.model.addPhysicalGroup(2, [s_s4, s_s8], 1, "x_min")   # x=0
    gmsh.model.addPhysicalGroup(2, [s_s2, s_s6], 2, "x_max")   # x=a
    gmsh.model.addPhysicalGroup(2, [s_s1, s_s5], 3, "y_min")   # y=0
    gmsh.model.addPhysicalGroup(2, [s_s3, s_s7], 4, "y_max")   # y=a
    gmsh.model.addPhysicalGroup(2, [s_bot], 5, "port1")         # z=-L
    gmsh.model.addPhysicalGroup(2, [s_top], 6, "port2")         # z=L

    # Transfinite constraints for structured hex mesh.
    for l in [l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12]
        gmsh.model.mesh.setTransfiniteCurve(l, n_xy + 1)
    end
    for l in [l13, l14, l15, l16, l17, l18, l19, l20]
        gmsh.model.mesh.setTransfiniteCurve(l, n_z + 1)
    end

    for s in [s_bot, s_int, s_top, s_s1, s_s2, s_s3, s_s4, s_s5, s_s6, s_s7, s_s8]
        gmsh.model.mesh.setTransfiniteSurface(s)
        gmsh.model.mesh.setRecombine(2, s)
    end

    gmsh.model.mesh.setTransfiniteVolume(v_vac)
    gmsh.model.mesh.setTransfiniteVolume(v_die)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)

    # Write in Gmsh MSH 2.2 binary format (required by MFEM).
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(filename)

    if gui
        gmsh.fltk.run()
    end

    gmsh.finalize()
    return nothing
end
