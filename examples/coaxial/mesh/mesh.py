# Automatic translation of `mesh.jl` made with Google Gemini.

import gmsh
from pathlib import Path
from math import pi

def generate_coaxial_mesh(
    filename: Path,
    refinement: int = 2,
    order: int = 2,
    inner_diameter_mm: float = 1.6383,
    outer_diameter_mm: float = 5.461,
    length_mm: float = 40.0,
    verbose: int = 5,
    gui: bool = False
):
    """
    Generate a mesh for the coaxial cable example using Gmsh

    :param filename: the filename to use for the generated mesh
    :param refinement: measure of how many elements to include, 0 is least
    :param order: the polynomial order of the approximation, minimum 1
    :param inner_diameter_mm: the inner diameter of the cable, in millimeters
    :param outer_diameter_mm: the outer diameter of the cable, in millimeters
    :param length_mm: the length of the cable, in millimeters
    :param verbose: flag to dictate the level of print to console, passed to Gmsh
    :param gui: whether to launch the Gmsh GUI on mesh generation
    """
    assert outer_diameter_mm > inner_diameter_mm > 0, "Outer diameter must be greater than inner diameter, and both must be positive."
    assert length_mm > 0, "Length must be positive."
    assert refinement >= 0, "Refinement must be non-negative."
    assert order > 0, "Order must be positive."

    kernel = gmsh.model.occ

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    if "coaxial" in gmsh.model.list():
        gmsh.model.setCurrent("coaxial")
        gmsh.model.remove() # Removes current model and all its entities
    gmsh.model.add("coaxial")

    # Geometry parameters (in mm)
    ri = inner_diameter_mm / 2
    ro = outer_diameter_mm / 2

    # Mesh parameters
    n_circum = 4 * (2**refinement) # min 4 elements in round
    n_length = 4 * (2**refinement) # min 4 elements on length

    # Geometry
    p0 = kernel.addPoint(ri, 0.0, 0.0)
    p1 = kernel.addPoint(ro, 0.0, 0.0)

    l0 = kernel.addLine(p0, p1)

    base_face_0_list = kernel.revolve(
        [(1, l0)],  # Entities to revolve: list of (dim, tag)
        0.0, 0.0, 0.0,  # Center of revolution
        0.0, 0.0, 1.0,  # Axis of revolution (z-axis)
        pi,             # Angle (pi radians)
        [n_circum // 2], # Mesh element count
        [1.0],          # Mesh element ratios
        True            # Remove duplicates
    )
    base_face_1_list = kernel.revolve(
        [(1, l0)],
        0.0, 0.0, 0.0,
        0.0, 0.0, 1.0,
        -pi,
        [n_circum // 2],
        [1.0],
        True
    )

    base_face_0 = [e for e in base_face_0_list if e[0] == 2]
    base_face_1 = [e for e in base_face_1_list if e[0] == 2]
    assert len(base_face_0) == 1 and len(base_face_1) == 1

    cylinder_0_entities = kernel.extrude(
        base_face_0,
        0.0, 0.0, length_mm,
        [n_length],
        [1.0],
        True
    )
    cylinder_1_entities = kernel.extrude(
        base_face_1,
        0.0, 0.0, length_mm,
        [n_length],
        [1.0],
        True
    )

    # Extract the tags (the second element of the (dim, tag) tuple)
    base_face_0_tag = base_face_0[0][1]
    base_face_1_tag = base_face_1[0][1]

    # Find the 2D and 3D entities in the results
    far_face_0_tag = [e[1] for e in cylinder_0_entities if e[0] == 2 and e[1] != base_face_0_tag][0]
    cylinder_0_tag = [e[1] for e in cylinder_0_entities if e[0] == 3][0]

    far_face_1_tag = [e[1] for e in cylinder_1_entities if e[0] == 2 and e[1] != base_face_1_tag][0]
    cylinder_1_tag = [e[1] for e in cylinder_1_entities if e[0] == 3][0]

    # Remove duplicates but preserves tags for non-removed objects
    kernel.fragment(kernel.getEntities(), [])
    kernel.synchronize()

    # 1. Get all boundaries of the two cylinders
    boundaries = []

    _, local_boundaries_0 = gmsh.model.getAdjacencies(3, cylinder_0_tag)
    _, local_boundaries_1 = gmsh.model.getAdjacencies(3, cylinder_1_tag)

    # Combine and find unique boundaries
    all_boundaries = list(local_boundaries_0) + list(local_boundaries_1)

    # 2. Identify the faces that are the shared *internal* boundary (due to the fragment)
    all_faces = gmsh.model.getEntities(2)
    all_face_tags = [tag for dim, tag in all_faces]

    # Port faces to exclude from the 'boundaries' group
    port_faces_tags = [
        base_face_0_tag,
        base_face_1_tag,
        far_face_0_tag,
        far_face_1_tag
    ]

    boundaries_tags = [
        tag for tag in all_face_tags if tag not in port_faces_tags
    ]

    # Add physical groups
    cylinder_group = gmsh.model.addPhysicalGroup(
        3, [cylinder_0_tag, cylinder_1_tag], name="cylinder"
    )
    boundary_group = gmsh.model.addPhysicalGroup(
        2, boundaries_tags, name="boundaries"
    )

    port1_group = gmsh.model.addPhysicalGroup(
        2, [base_face_0_tag, base_face_1_tag], name="port1"
    )
    port2_group = gmsh.model.addPhysicalGroup(
        2, [far_face_0_tag, far_face_1_tag], name="port2"
    )

    # Generate mesh
    gmsh.option.setNumber("Mesh.MinimumCurveNodes", 2)
    gmsh.option.setNumber("Mesh.MinimumCircleNodes", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    # Set meshing algorithms
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)

    # Generate the 3D mesh
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)

    filepath = Path(filename)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(str(filepath))

    # Print some information
    if verbose > 0:
        print("\nFinished generating mesh. Physical group tags:")
        print(f"Cylinder: {cylinder_group}")
        print(f"Boundaries: {boundary_group}")
        print(f"Port 1: {port1_group}")
        print(f"Port 2: {port2_group}")
        print()

    # Optionally launch GUI
    if gui:
        gmsh.fltk.run()

    return gmsh.finalize()

if __name__ == '__main__':
    generate_coaxial_mesh(filename="coaxial_py.msh", verbose=5, gui=True)
