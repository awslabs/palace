# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Iris-coupled, dielectric-slab-loaded rectangular waveguide filter (two-port). Cross-section
# [0,a]x[0,b]: slab eps_r fills x in [0,d], air fills [d,a], over the full length L in z. The
# full-length slab keeps both wave-port modes hybrid LSM (frequency-rotating shape), so the
# modal correction W is rank>=2 -- unlike an air-filled guide whose TE10 port mode gives W=0.
# Two inductive PEC posts at z1, z2 form a half-wave cavity whose resonance is a strongly-
# coupled |S21| transmission peak, so real-axis greedy PROM sampling captures the pole.
#
# Physical groups: dielectric(vol), air(vol), port_lo(z=0), port_hi(z=L),
#                  walls(all remaining boundary faces = outer sides + post surfaces = PEC).
using Gmsh: gmsh

function gen_iris_filter(;
    filename,
    a=2.0,
    b=1.0,
    L=5.0,
    d=1.0,
    eps_r=10.0,
    z1=1.5,          # first post center (z)
    z2=3.5,          # second post center (z) -> cavity length Lc = z2 - z1
    ap=0.3,          # iris aperture: guide stays open only for x in [0, ap] (a slit in the
    # high-field slab region); the post blocks x in [ap, a] over full height
    post_t=0.2,      # post thickness along z
    h=0.22,          # target mesh size
    order=2,
    verbose=2
)
    kernel = gmsh.model.occ
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)
    if "iris" in gmsh.model.list()
        gmsh.model.setCurrent("iris")
        gmsh.model.remove()
    end
    gmsh.model.add("iris")

    diel = kernel.addBox(0.0, 0.0, 0.0, d, b, L)
    air = kernel.addBox(d, 0.0, 0.0, a - d, b, L)
    # glue slab<->air first (conformal interface at x=d), then carve the irises out of the
    # whole cross-section so the aperture is a narrow slit x in [0, ap] in the high-field slab
    frag0, _ = kernel.fragment([(3, diel)], [(3, air)])
    guide_vols = [(dim, tag) for (dim, tag) in frag0 if dim == 3]
    # two inductive irises: PEC blocks spanning x in [ap, a], full height y, thin in z
    p1 = kernel.addBox(ap, 0.0, z1 - post_t / 2, a - ap, b, post_t)
    p2 = kernel.addBox(ap, 0.0, z2 - post_t / 2, a - ap, b, post_t)
    frag, _ = kernel.cut(guide_vols, [(3, p1), (3, p2)])
    kernel.synchronize()

    # identify volumes by centroid x: slab (x<d) vs air (x>d). The irises may split each
    # column into several z-pieces, so collect lists on both sides.
    diel_vols = Int[]
    air_vols = Int[]
    for (dim, tag) in frag
        dim == 3 || continue
        x, _, _ = gmsh.model.occ.getCenterOfMass(3, tag)
        if x < d
            push!(diel_vols, tag)
        else
            push!(air_vols, tag)
        end
    end
    @assert !isempty(diel_vols) && !isempty(air_vols)

    gmsh.model.addPhysicalGroup(3, diel_vols, -1, "dielectric")
    gmsh.model.addPhysicalGroup(3, air_vols, -1, "air")

    eps = 1e-6 * max(a, b, L)
    faces_in(x0, y0, z0, x1, y1, z1_) = last.(
        gmsh.model.getEntitiesInBoundingBox(
            x0 - eps,
            y0 - eps,
            z0 - eps,
            x1 + eps,
            y1 + eps,
            z1_ + eps,
            2
        )
    )
    port_lo = faces_in(0.0, 0.0, 0.0, a, b, 0.0)
    port_hi = faces_in(0.0, 0.0, L, a, b, L)

    # walls (PEC) = every boundary face that is not a port. This sweeps in the four outer
    # side walls AND the exposed post surfaces automatically. The internal slab/air interface
    # at x=d is shared between two volumes (not a boundary face) so it is excluded.
    all_faces = last.(gmsh.model.getEntities(2))
    port_set = Set(vcat(port_lo, port_hi))
    walls = [f for f in all_faces if !(f in port_set) && is_boundary_face(f)]

    gmsh.model.addPhysicalGroup(2, port_lo, -1, "port_lo")
    gmsh.model.addPhysicalGroup(2, port_hi, -1, "port_hi")
    gmsh.model.addPhysicalGroup(2, walls, -1, "walls")

    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(filename)
    println("a=$a b=$b L=$L d=$d eps_r=$eps_r  Lc=$(z2 - z1) ap=$ap post_t=$post_t")
    println(
        "diel=$diel_vols air=$air_vols  #port_lo=$(length(port_lo)) " *
        "#port_hi=$(length(port_hi)) #walls=$(length(walls))"
    )
    println("groups: dielectric air port_lo port_hi walls")
    return gmsh.finalize()
end

# A dim-2 entity is a boundary face iff it bounds exactly one volume.
function is_boundary_face(ftag)
    up, _ = gmsh.model.getAdjacencies(2, ftag)
    return length(up) == 1
end

gen_iris_filter(
    filename=get(ARGS, 1, joinpath(@__DIR__, "iris_filter.msh")),
    L=parse(Float64, get(ARGS, 2, "5.0")),
    eps_r=parse(Float64, get(ARGS, 3, "10.0"))
)
