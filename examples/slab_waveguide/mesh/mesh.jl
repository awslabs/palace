# Dielectric-slab-loaded rectangular waveguide, structured hex mesh.
# Cross-section [0,a]x[0,b] in x-y; slab eps_r fills x in [0,d], air fills [d,a];
# propagation along z, length L. Transverse inhomogeneity -> hybrid LSM modes with
# E_z != 0 and strongly frequency-dependent (rotating) transverse shape.
#
# Physical groups: dielectric(vol), air(vol), port_lo(z=0), port_hi(z=L), walls(sides).
using Gmsh: gmsh

function gen_slab_guide(; filename, a=2.0, b=1.0, L=2.0, d=1.0, eps_r=10.0,
                          nx_diel=6, nx_air=6, ny=6, nz=8, order=2, verbose=2)
    kernel = gmsh.model.occ
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)
    if "slab" in gmsh.model.list()
        gmsh.model.setCurrent("slab"); gmsh.model.remove()
    end
    gmsh.model.add("slab")

    diel = kernel.addBox(0.0, 0.0, 0.0, d, b, L)
    air  = kernel.addBox(d, 0.0, 0.0, a - d, b, L)
    # glue the shared interface so meshes are conformal
    out, _ = kernel.fragment([(3, diel)], [(3, air)])
    kernel.synchronize()

    # identify the two volumes by centroid x
    diel_vol = -1; air_vol = -1
    for (dim, tag) in out
        dim == 3 || continue
        x, _, _ = gmsh.model.occ.getCenterOfMass(3, tag)
        if x < d
            diel_vol = tag
        else
            air_vol = tag
        end
    end
    @assert diel_vol > 0 && air_vol > 0

    gmsh.model.addPhysicalGroup(3, [diel_vol], -1, "dielectric")
    gmsh.model.addPhysicalGroup(3, [air_vol], -1, "air")

    eps = 1e-6 * max(a, b, L)
    faces_at(zlo, zhi) = last.(gmsh.model.getEntitiesInBoundingBox(
        -eps, -eps, zlo - eps, a + eps, b + eps, zhi + eps, 2))
    port_lo = faces_at(0.0, 0.0)
    port_hi = faces_at(L, L)
    # walls = the four OUTER side faces only (x=0, x=a, y=0, y=b). The fragment leaves an
    # internal dielectric/air interface at x=d spanning the guide; it must NOT be swept into
    # the PEC walls group (that would put a conducting sheet through the guide). Select walls
    # by outer-plane bounding boxes so the interface (centroid x=d, interior) is excluded.
    on_plane(x0, y0, z0, x1, y1, z1) = last.(gmsh.model.getEntitiesInBoundingBox(
        x0 - eps, y0 - eps, z0 - eps, x1 + eps, y1 + eps, z1 + eps, 2))
    walls = unique(vcat(
        on_plane(0.0, 0.0, 0.0, 0.0, b, L),   # x = 0
        on_plane(a, 0.0, 0.0, a, b, L),        # x = a
        on_plane(0.0, 0.0, 0.0, a, 0.0, L),    # y = 0
        on_plane(0.0, b, 0.0, a, b, L)))       # y = b

    gmsh.model.addPhysicalGroup(2, port_lo, -1, "port_lo")
    gmsh.model.addPhysicalGroup(2, port_hi, -1, "port_hi")
    gmsh.model.addPhysicalGroup(2, walls, -1, "walls")

    # structured hex: transfinite curves along x (per region), y, z + recombine
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    for (dim, tag) in gmsh.model.getEntities(1)
        x0, y0, z0 = gmsh.model.occ.getCenterOfMass(1, tag)
        bb = gmsh.model.getBoundingBox(1, tag)
        dx = bb[4] - bb[1]; dy = bb[5] - bb[2]; dz = bb[6] - bb[3]
        if dz > 0.5L                      # z-direction edge
            gmsh.model.mesh.setTransfiniteCurve(tag, nz + 1)
        elseif dy > 0.5b                  # y-direction edge
            gmsh.model.mesh.setTransfiniteCurve(tag, ny + 1)
        else                              # x-direction edge (within one region)
            n = (x0 < d) ? nx_diel : nx_air
            gmsh.model.mesh.setTransfiniteCurve(tag, n + 1)
        end
    end
    for (dim, tag) in gmsh.model.getEntities(2)
        gmsh.model.mesh.setTransfiniteSurface(tag)
        gmsh.model.mesh.setRecombine(2, tag)
    end
    for (dim, tag) in gmsh.model.getEntities(3)
        gmsh.model.mesh.setTransfiniteVolume(tag)
    end

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(filename)
    println("a=$a b=$b L=$L d=$d eps_r=$eps_r  diel=$diel_vol air=$air_vol")
    println("groups: dielectric air port_lo port_hi walls")
    gmsh.finalize()
end

gen_slab_guide(filename=get(ARGS, 1, joinpath(@__DIR__, "slab_waveguide.msh")),
               L=parse(Float64, get(ARGS, 2, "2.0")),
               eps_r=parse(Float64, get(ARGS, 3, "10.0")))
