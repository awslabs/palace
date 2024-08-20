using Gmsh: gmsh

#=
  make_box_mesh(;hex_mesh = true, verbose = 10, periodic =[false,false,false], filename="")

Generate a mesh of a box for testing periodicity
  Arguments:
  - hex_mesh - Whether the generated mesh should be hexahedral
  - verbose - Verbosity setting for gmsh
  - length - length of the box for each direction, [lx, ly, lz]
  - res - element size at each direction [rez_x, rez_y, rez_z]
  - periodic - which faces to be periodic, [x periodic, y periodic, z periodic]
  - filename - filename to save the generated mesh, if empty the mesh is not saved
  - gui - Opens gui if true
=#
function make_box_mesh(;
    hex_mesh=true,
    verbose=10,
    length=[1., 1., 1.],
    res=[2/3, 2/3, 2/3],
    periodic=[false, false, false],
    filename="",
    gui=false
)
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    # Add model
    if "box" in gmsh.model.list()
        gmsh.model.setCurrent("box")
        gmsh.model.remove()
    end
    gmsh.model.add("box")



    kernel = gmsh.model.occ
    @show box = kernel.add_box(0.0, 0.0, 0.0, length[1], length[2], length[3])
    @show _, box_boundary = kernel.getSurfaceLoops(box)

    kernel.synchronize()

    gmsh.option.setNumber("Mesh.MeshSizeMin", minimum(res)) # TODO: This doesn't set the resolution for each direction
    gmsh.option.setNumber("Mesh.MeshSizeMax", maximum(res))

    # Makes the surfaces of the hex transfinite
    if hex_mesh
        gmsh.model.mesh.set_transfinite_automatic()
    end

    gmsh.model.add_physical_group(3, [box], -1, "Volume")
    gmsh.model.add_physical_group(2, [1], 1, "Xmin")
    gmsh.model.add_physical_group(2, [2], 2, "Xmax")
    gmsh.model.add_physical_group(2, [3], 3, "Ymin")
    gmsh.model.add_physical_group(2, [4], 4, "Ymax")
    gmsh.model.add_physical_group(2, [5], 5, "Zmin")
    gmsh.model.add_physical_group(2, [6], 6, "Zmax")

    gmsh.model.mesh.generate(3)

    # Add in periodicity -- the row vector is a 4x4 affine matrix in row format -- see
    # https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_13_1/tutorials/julia/t18.jl

    # Hard coded surface pairings for the box.
    if periodic[1]
        gmsh.model.mesh.set_periodic(
            2,
            [2],
            [1],
            [1, 0, 0, length[1], 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
        )
    end
    if periodic[2]
        gmsh.model.mesh.setPeriodic(
            2,
            [4],
            [3],
            [1, 0, 0, 0, 0, 1, 0, length[2], 0, 0, 1, 0, 0, 0, 0, 1]
        )
    end
    if periodic[3]
        gmsh.model.mesh.setPeriodic(
            2,
            [6],
            [5],
            [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, length[3], 0, 0, 0, 1]
        )
    end

    if gui
        gmsh.fltk.run()
    end

    if !isempty(filename)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.write(joinpath(@__DIR__, filename))
    end
    return gmsh.finalize()
end
