#=
    Render real single transmon → mesh → (Palace later)

    Uses OCC kernel (same geometry as EGADS, but meshes correctly via Gmsh
    with proper physical groups so the mesh conforms to island, resonator,
    feedline boundaries).

    Usage:
        cd examples/single_transmon
        julia --project=/path/to/DeviceLayout.jl-ad-sensitivity render_and_solve.jl
=#

using DeviceLayout, DeviceLayout.SchematicDrivenLayout, DeviceLayout.PreferredUnits
import DeviceLayout.SolidModels
import DeviceLayout.SolidModels: SolidModel
import DeviceLayout.SchematicDrivenLayout: ExamplePDK
import .ExamplePDK: LayerVocabulary, add_bridges!
using .ExamplePDK.Transmons, .ExamplePDK.ReadoutResonators
import .ExamplePDK.SimpleJunctions: ExampleSimpleJunction
import DeviceLayout: save

# ─── Build the same schematic as SingleTransmon.jl ───────────────────────────

reset_uniquename!()

cpw_width = 10μm; cpw_gap = 6μm
PATH_STYLE = Paths.SimpleCPW(cpw_width, cpw_gap)
BRIDGE_STYLE = ExamplePDK.bridge_geometry(PATH_STYLE)

cap_length = 620μm; cap_gap = 30μm; cap_width = 24μm
coupling_gap = 5μm; bend_radius = 50μm; n_meander_turns = 5
hanger_length = 500μm; arm_length = 428μm; total_length = 5000μm
w_grasp = cap_width + 2 * cap_gap
total_height = arm_length + coupling_gap + Paths.extent(PATH_STYLE) +
               hanger_length + (3 + n_meander_turns * 2) * bend_radius

qubit = ExampleRectangleTransmon(;
    jj_template=ExampleSimpleJunction(), name="qubit",
    cap_length, cap_gap, cap_width)

rres = ExampleClawedMeanderReadout(;
    name="rres", coupling_length=400μm, coupling_gap, total_length,
    w_shield=2μm, w_claw=34μm, l_claw=121μm, claw_gap=6μm,
    w_grasp, n_meander_turns, total_height, hanger_length,
    bend_radius, bridge=BRIDGE_STYLE)

readout_length = 2700μm
p_readout = Path(Point(0μm, 0μm); α0=π/2, name="p_ro",
                 metadata=LayerVocabulary.METAL_NEGATIVE)
straight!(p_readout, readout_length / 2, PATH_STYLE)
straight!(p_readout, readout_length / 2, PATH_STYLE)

csport = CoordinateSystem(uniquename("port"), nm)
render!(csport, only_simulated(centered(Rectangle(cpw_width, cpw_width))),
        LayerVocabulary.PORT)
attach!(p_readout, sref(csport), cpw_width, i=1)
attach!(p_readout, sref(csport), readout_length / 2 - cpw_width, i=2)

g = SchematicGraph("single-transmon")
qubit_node = add_node!(g, qubit)
rres_node = fuse!(g, qubit_node, rres)
p_readout_node = add_node!(g, p_readout)
attach!(g, p_readout_node, rres_node => :feedline, 0mm, location=1)

floorplan = plan(g)
# Skip bridges for now — they cause PLC mesh errors and don't affect eigenmode
# add_bridges!(floorplan, BRIDGE_STYLE, spacing=300μm)

substrate_x = 4mm; substrate_y = 3.7mm
center_xyz = DeviceLayout.center(floorplan)
chip = centered(Rectangle(substrate_x, substrate_y), on_pt=center_xyz)
sim_area = centered(Rectangle(substrate_x, substrate_y), on_pt=center_xyz)
render!(floorplan.coordinate_system, sim_area, LayerVocabulary.SIMULATED_AREA)
render!(floorplan.coordinate_system, sim_area, LayerVocabulary.WRITEABLE_AREA)
render!(floorplan.coordinate_system, chip, LayerVocabulary.CHIP_AREA)
check!(floorplan)

println("Schematic built.")

# ─── Render + Mesh via OCC/Gmsh ──────────────────────────────────────────────

# Bridge-free target: the AD sensitivity branch bridge postrender ops are broken
# (bridge/bridge_base physical groups not generated). Use a target that skips them.
tech = SolidModelTarget(
    ExamplePDK.EXAMPLE_SINGLECHIP_TECHNOLOGY;
    simulation=true,
    bounding_layers=[:simulated_area],
    substrate_layers=[:chip_area],
    indexed_layers=[:port, :lumped_element, :integration, :wave_port],
    wave_port_layers=[:wave_port],
    ignored_layers=[:bridge_base, :bridge],
    postrender_ops=Any[
        ("metal_negative", SolidModels.union_geom!,
         ("metal_negative", "metal_negative", 2, 2), :remove_object => true),
        ("metal", SolidModels.difference_geom!,
         ("writeable_area", "metal_negative", 2, 2), :remove_object => true),
        ("metal", SolidModels.union_geom!,
         ("metal", "metal_positive", 2, 2), :remove_tool => true),
        ("substrate", SolidModels.union_geom!,
         ("chip_area_extrusion", "chip_area_extrusion", 3, 3),
         :remove_object => true, :remove_tool => true),
        ("vacuum", SolidModels.difference_geom!,
         ("simulated_area_extrusion", "substrate", 3, 3)),
        ("metal", SolidModels.difference_geom!,
         ("metal", "port"), :remove_object => true),
        ("metal", SolidModels.difference_geom!,
         ("metal", ["port_1", "port_2", "lumped_element"]), :remove_object => true),
    ],
    retained_physical_groups=[
        ("vacuum", 3), ("substrate", 3), ("metal", 2), ("exterior_boundary", 2),
        ("port_1", 2), ("port_2", 2), ("lumped_element", 2),
    ]
)
sm = SolidModel("single_transmon", overwrite=true)
SolidModels.set_gmsh_option("General.Verbosity", 1)
SolidModels.mesh_order(1)
println("Rendering (OCC, bridge-free)...")
render!(sm, floorplan, tech)
println("Render complete.")

println("Generating mesh...")
SolidModels.gmsh.model.mesh.generate(3)

outdir = @__DIR__
meshfile = joinpath(outdir, "single_transmon.msh2")
save(meshfile, sm)
println("Saved: $meshfile")

attrs = SolidModels.attributes(sm)
println("Attributes: $attrs")
