using Gmsh: gmsh

#=
  make_cube_mesh(;hex_mesh = true, verbose = 10, periodic =[false,false,false], filename="")

Generate a mesh of a cube for testing periodicity
  Arguments:
  - hex_mesh - Whether the generated mesh should be hexahedral
  - verbose - Verbosity setting for gmsh
  - periodic - which faces to be periodic, [x periodic, y periodic, z periodic]
  - filename - filename to save the generated mesh, if empty the mesh is not saved
=#
function make_cube_mesh(;
    hex_mesh=true,
    verbose=10,
    periodic=[false, false, false],
    filename=""
)
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    # Add model
    if "cube" in gmsh.model.list()
        gmsh.model.setCurrent("cube")
        gmsh.model.remove()
    end
    gmsh.model.add("cube")

    length = 1.0

    kernel = gmsh.model.occ
    @show cube = kernel.add_box(0.0, 0.0, 0.0, length, length, length)
    @show _, cube_boundary = kernel.getSurfaceLoops(cube)

    kernel.synchronize()

    gmsh.option.setNumber("Mesh.MeshSizeMin", length / 1.5)
    gmsh.option.setNumber("Mesh.MeshSizeMax", length / 1.5)

    # Makes the surfaces of the hex transfinite
    if hex_mesh
        gmsh.model.mesh.set_transfinite_automatic()
    end

    # Add in periodicity -- the row vector is a 4x4 affine matrix in row format -- see
    # https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_13_1/tutorials/julia/t18.jl

    # Hard coded surface pairings for the box.
    if periodic[1]
        gmsh.model.mesh.set_periodic(
            2,
            [2],
            [1],
            [1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
        )
    end
    if periodic[2]
        gmsh.model.mesh.setPeriodic(
            2,
            [4],
            [3],
            [1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1]
        )
    end
    if periodic[3]
        gmsh.model.mesh.setPeriodic(
            2,
            [6],
            [5],
            [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1]
        )
    end

    gmsh.model.add_physical_group(3, [cube], -1, "Volume")
    gmsh.model.add_physical_group(2, [1], 1, "Xmin")
    gmsh.model.add_physical_group(2, [2], 2, "Xmax")
    gmsh.model.add_physical_group(2, [3], 3, "Ymin")
    gmsh.model.add_physical_group(2, [4], 4, "Ymax")
    gmsh.model.add_physical_group(2, [5], 5, "Zmin")
    gmsh.model.add_physical_group(2, [6], 6, "Zmax")

    gmsh.model.mesh.generate(3)

    gmsh.fltk.run()

    if isempty(filename)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 0)
        gmsh.write(joinpath(@__DIR__, "cube.msh"))
    end
    return gmsh.finalize()
end
