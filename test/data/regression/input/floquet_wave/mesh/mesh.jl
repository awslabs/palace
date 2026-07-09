# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate mesh for the Floquet + wave port test (from the floquet_wave directory):
#   julia -e 'include("mesh/mesh.jl"); generate_floquet_wave_mesh()'

using Gmsh: gmsh

"""
    generate_floquet_wave_mesh(; filename, ...)

Generate a structured hex mesh for an open-ended waveguide array unit cell: a doubly
periodic cell with a Floquet port on top, fed by a rectangular waveguide from below with
a wave port at its base. Used for testing mixed wave/Floquet port S-parameter power
normalization.

# Geometry

  - Periodic cell: `a` x `b` in xy, z in [0, L_top]. Floquet port at z = L_top.
  - Waveguide: `w` x `h` centered in xy, z in [-L_wg, 0]. Wave port at z = -L_wg.
  - Ground-plane flange (PEC) at z = 0 outside the waveguide aperture.
  - Waveguide side walls are PEC.

With w = 2, h = 1 cm the TE10 cutoff is 7.5 GHz; the first Floquet grating lobe onsets
at c/a = 12 GHz for a = 2.5 cm. An 8-11 GHz sweep is single-mode on both sides, so for
the lossless structure |S_wave_refl|^2 + |S_floquet<-wave|^2 = 1 exactly (and the
reciprocal relation for the Floquet-driven column). TE10 is y-polarized, matching the
Floquet TE polarization at normal incidence.

# Attributes

  - Volume 1: cell above the ground plane, Volume 2: waveguide
  - Boundary 1,2: x periodic pair
  - Boundary 3,4: y periodic pair
  - Boundary 5: Floquet port (z = L_top)
  - Boundary 6: wave port (z = -L_wg)
  - Boundary 7: PEC (flange at z = 0 and waveguide walls)

Units: cm (Palace L0 = 1e-2 for cm -> m conversion).
"""
function generate_floquet_wave_mesh(;
    filename::AbstractString="mesh/floquet_wave.msh",
    a::Real=2.5,
    b::Real=2.5,
    w::Real=2.0,
    h::Real=1.0,
    L_top::Real=2.5,
    L_wg::Real=1.5,
    n_ap_x::Integer=6,
    n_ap_y::Integer=4,
    n_rim::Integer=1,
    n_z_top::Integer=6,
    n_z_wg::Integer=4,
    order::Integer=1,
    verbose::Integer=5,
    gui::Bool=false
)
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    if "floquet_wave" in gmsh.model.list()
        gmsh.model.setCurrent("floquet_wave")
        gmsh.model.remove()
    end
    gmsh.model.add("floquet_wave")

    geo = gmsh.model.geo

    # Grid: 4 x-positions, 4 y-positions, 3 z-levels. The central xy block
    # (ix,iy = 2..3) matches the waveguide cross-section.
    xs = [0.0, (a - w) / 2.0, (a + w) / 2.0, a]
    ys = [0.0, (b - h) / 2.0, (b + h) / 2.0, b]
    zs = [-L_wg, 0.0, L_top]

    nxs = [n_rim, n_ap_x, n_rim]
    nys = [n_rim, n_ap_y, n_rim]

    # Points: full 4x4 grid at z-levels 2 (z=0) and 3 (z=L_top); central 2x2 at level 1.
    p = Dict{Tuple{Int,Int,Int},Int}()
    for iz in 2:3, iy in 1:4, ix in 1:4
        p[(ix, iy, iz)] = geo.addPoint(xs[ix], ys[iy], zs[iz])
    end
    for iy in 2:3, ix in 2:3
        p[(ix, iy, 1)] = geo.addPoint(xs[ix], ys[iy], zs[1])
    end

    # Lines along x
    lx = Dict{Tuple{Int,Int,Int},Int}()
    for iz in 2:3, iy in 1:4, ix in 1:3
        lx[(ix, iy, iz)] = geo.addLine(p[(ix, iy, iz)], p[(ix + 1, iy, iz)])
    end
    for iy in 2:3
        lx[(2, iy, 1)] = geo.addLine(p[(2, iy, 1)], p[(3, iy, 1)])
    end

    # Lines along y
    ly = Dict{Tuple{Int,Int,Int},Int}()
    for iz in 2:3, iy in 1:3, ix in 1:4
        ly[(ix, iy, iz)] = geo.addLine(p[(ix, iy, iz)], p[(ix, iy + 1, iz)])
    end
    for ix in 2:3
        ly[(ix, 2, 1)] = geo.addLine(p[(ix, 2, 1)], p[(ix, 3, 1)])
    end

    # Lines along z: upper (level 2->3) full grid; lower (level 1->2) central only.
    lz = Dict{Tuple{Int,Int,Int},Int}()
    for iy in 1:4, ix in 1:4
        lz[(ix, iy, 2)] = geo.addLine(p[(ix, iy, 2)], p[(ix, iy, 3)])
    end
    for iy in 2:3, ix in 2:3
        lz[(ix, iy, 1)] = geo.addLine(p[(ix, iy, 1)], p[(ix, iy, 2)])
    end

    function make_surface(l1, l2, l3, l4)
        cl = geo.addCurveLoop([l1, l2, l3, l4])
        return geo.addPlaneSurface([cl])
    end

    # Z-faces (normal +z): 3x3 at levels 2 and 3; central at level 1.
    sz = Dict{Tuple{Int,Int,Int},Int}()
    for iz in 2:3, iy in 1:3, ix in 1:3
        sz[(ix, iy, iz)] = make_surface(lx[(ix, iy, iz)], ly[(ix + 1, iy, iz)],
                                        -lx[(ix, iy + 1, iz)], -ly[(ix, iy, iz)])
    end
    sz[(2, 2, 1)] = make_surface(lx[(2, 2, 1)], ly[(3, 2, 1)],
                                 -lx[(2, 3, 1)], -ly[(2, 2, 1)])

    # X-faces (normal +x): upper region full; lower region central.
    sx = Dict{Tuple{Int,Int,Int},Int}()
    for iy in 1:3, ix in 1:4
        sx[(ix, iy, 2)] = make_surface(ly[(ix, iy, 2)], lz[(ix, iy + 1, 2)],
                                       -ly[(ix, iy, 3)], -lz[(ix, iy, 2)])
    end
    for ix in 2:3
        sx[(ix, 2, 1)] = make_surface(ly[(ix, 2, 1)], lz[(ix, 3, 1)],
                                      -ly[(ix, 2, 2)], -lz[(ix, 2, 1)])
    end

    # Y-faces (normal -y): upper region full; lower region central.
    sy = Dict{Tuple{Int,Int,Int},Int}()
    for iy in 1:4, ix in 1:3
        sy[(ix, iy, 2)] = make_surface(lx[(ix, iy, 2)], lz[(ix + 1, iy, 2)],
                                       -lx[(ix, iy, 3)], -lz[(ix, iy, 2)])
    end
    for iy in 2:3
        sy[(2, iy, 1)] = make_surface(lx[(2, iy, 1)], lz[(3, iy, 1)],
                                      -lx[(2, iy, 2)], -lz[(2, iy, 1)])
    end

    # Volumes: 9 upper + 1 waveguide.
    upper_vols = Int[]
    for iy in 1:3, ix in 1:3
        sl = geo.addSurfaceLoop([
            -sz[(ix, iy, 2)], sz[(ix, iy, 3)],
            -sx[(ix, iy, 2)], sx[(ix + 1, iy, 2)],
            sy[(ix, iy, 2)], -sy[(ix, iy + 1, 2)]
        ])
        push!(upper_vols, geo.addVolume([sl]))
    end
    sl_wg = geo.addSurfaceLoop([
        -sz[(2, 2, 1)], sz[(2, 2, 2)],
        -sx[(2, 2, 1)], sx[(3, 2, 1)],
        sy[(2, 2, 1)], -sy[(2, 3, 1)]
    ])
    wg_vol = geo.addVolume([sl_wg])

    geo.synchronize()

    # Physical groups
    gmsh.model.addPhysicalGroup(3, upper_vols, 1, "cell")
    gmsh.model.addPhysicalGroup(3, [wg_vol], 2, "waveguide")

    xmin_surfs = [sx[(1, iy, 2)] for iy in 1:3]
    xmax_surfs = [sx[(4, iy, 2)] for iy in 1:3]
    ymin_surfs = [sy[(ix, 1, 2)] for ix in 1:3]
    ymax_surfs = [sy[(ix, 4, 2)] for ix in 1:3]
    gmsh.model.addPhysicalGroup(2, xmin_surfs, 1, "x_min")
    gmsh.model.addPhysicalGroup(2, xmax_surfs, 2, "x_max")
    gmsh.model.addPhysicalGroup(2, ymin_surfs, 3, "y_min")
    gmsh.model.addPhysicalGroup(2, ymax_surfs, 4, "y_max")
    gmsh.model.addPhysicalGroup(2, [sz[(ix, iy, 3)] for iy in 1:3 for ix in 1:3], 5,
                                "floquet_port")
    gmsh.model.addPhysicalGroup(2, [sz[(2, 2, 1)]], 6, "wave_port")
    # PEC: ground-plane flange (all z=0 faces except the aperture) + waveguide walls.
    pec_surfs = [sz[(ix, iy, 2)] for iy in 1:3 for ix in 1:3 if !(ix == 2 && iy == 2)]
    append!(pec_surfs, [sx[(2, 2, 1)], sx[(3, 2, 1)], sy[(2, 2, 1)], sy[(2, 3, 1)]])
    gmsh.model.addPhysicalGroup(2, pec_surfs, 7, "pec")

    # Transfinite constraints
    for iz in 2:3, iy in 1:4, ix in 1:3
        gmsh.model.mesh.setTransfiniteCurve(lx[(ix, iy, iz)], nxs[ix] + 1)
    end
    for iy in 2:3
        gmsh.model.mesh.setTransfiniteCurve(lx[(2, iy, 1)], nxs[2] + 1)
    end
    for iz in 2:3, iy in 1:3, ix in 1:4
        gmsh.model.mesh.setTransfiniteCurve(ly[(ix, iy, iz)], nys[iy] + 1)
    end
    for ix in 2:3
        gmsh.model.mesh.setTransfiniteCurve(ly[(ix, 2, 1)], nys[2] + 1)
    end
    for iy in 1:4, ix in 1:4
        gmsh.model.mesh.setTransfiniteCurve(lz[(ix, iy, 2)], n_z_top + 1)
    end
    for iy in 2:3, ix in 2:3
        gmsh.model.mesh.setTransfiniteCurve(lz[(ix, iy, 1)], n_z_wg + 1)
    end

    for (_, s) in sz
        gmsh.model.mesh.setTransfiniteSurface(s)
        gmsh.model.mesh.setRecombine(2, s)
    end
    for (_, s) in sx
        gmsh.model.mesh.setTransfiniteSurface(s)
        gmsh.model.mesh.setRecombine(2, s)
    end
    for (_, s) in sy
        gmsh.model.mesh.setTransfiniteSurface(s)
        gmsh.model.mesh.setRecombine(2, s)
    end
    for v in upper_vols
        gmsh.model.mesh.setTransfiniteVolume(v)
    end
    gmsh.model.mesh.setTransfiniteVolume(wg_vol)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)

    # Periodic mesh matching for the upper cell side walls.
    gmsh.model.mesh.setPeriodic(2, xmax_surfs, xmin_surfs,
                                [1, 0, 0, a, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
    gmsh.model.mesh.setPeriodic(2, ymax_surfs, ymin_surfs,
                                [1, 0, 0, 0, 0, 1, 0, b, 0, 0, 1, 0, 0, 0, 0, 1])

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(filename)

    if gui
        gmsh.fltk.run()
    end

    gmsh.finalize()
    return nothing
end
