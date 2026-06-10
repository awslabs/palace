# Inspect the EGADSKernel tape for the real SingleTransmon geometry.
# Shows what operations are recorded during postrender.

push!(LOAD_PATH, joinpath(ENV["HOME"], "AWS", "EGADS.jl"))
using EGADS
lib = joinpath(ENV["HOME"], "AWS", "EngineeringSketchPad", "EngSketchPad", "lib", "libegads.dylib")
isfile(lib) && EGADS.set_library_path!(lib)

using DeviceLayout, DeviceLayout.SchematicDrivenLayout, DeviceLayout.PreferredUnits
import DeviceLayout.SolidModels
import DeviceLayout.SolidModels: SolidModel, kernel
import DeviceLayout.SchematicDrivenLayout: SolidModelTarget
import DeviceLayout.SchematicDrivenLayout: ExamplePDK
import .ExamplePDK: LayerVocabulary
using .ExamplePDK.Transmons, .ExamplePDK.ReadoutResonators
import .ExamplePDK.SimpleJunctions: ExampleSimpleJunction

Base.retry_load_extensions()
ext = Base.get_extension(DeviceLayout, :DeviceLayoutEGADSExt)
isnothing(ext) && error("DeviceLayoutEGADSExt not loaded")
const EK = ext.EGADSKernel

# ─── Build floorplan (same as SingleTransmon_EGADS.jl) ──────────────────────

function build_floorplan(; cap_length=620μm, cap_gap=30μm, cap_width=24μm)
    reset_uniquename!()
    cpw_width = 10μm; cpw_gap = 6μm
    PATH_STYLE = Paths.SimpleCPW(cpw_width, cpw_gap)

    qubit = ExampleRectangleTransmon(;
        jj_template=ExampleSimpleJunction(), name="qubit",
        cap_length=cap_length, cap_gap=cap_gap, cap_width=cap_width)

    coupling_gap = 5μm
    w_grasp = cap_width + 2 * cap_gap
    arm_length = 428μm
    n_meander_turns = 5
    bend_radius = 50μm
    hanger_length = 500μm
    total_height = arm_length + coupling_gap + Paths.extent(PATH_STYLE) +
                   hanger_length + (3 + n_meander_turns * 2) * bend_radius

    rres = ExampleClawedMeanderReadout(;
        name="rres", coupling_length=400μm, coupling_gap=coupling_gap,
        total_length=5000μm, w_shield=2μm, w_claw=34μm, l_claw=121μm,
        claw_gap=6μm, w_grasp=w_grasp, n_meander_turns=n_meander_turns,
        total_height=total_height, hanger_length=hanger_length,
        bend_radius=bend_radius,
        bridge=ExamplePDK.bridge_geometry(PATH_STYLE))

    readout_length = 2700μm
    p_ro = Path(Point(0μm, 0μm); α0=π/2, name="p_ro",
                metadata=LayerVocabulary.METAL_NEGATIVE)
    straight!(p_ro, readout_length / 2, PATH_STYLE)
    straight!(p_ro, readout_length / 2, PATH_STYLE)

    csport = CoordinateSystem(uniquename("port"), nm)
    render!(csport, only_simulated(centered(Rectangle(cpw_width, cpw_width))),
            LayerVocabulary.PORT)
    attach!(p_ro, sref(csport), cpw_width, i=1)
    attach!(p_ro, sref(csport), readout_length / 2 - cpw_width, i=2)

    g = SchematicGraph("single-transmon-tape")
    qn = add_node!(g, qubit)
    rn = fuse!(g, qn, rres)
    pn = add_node!(g, p_ro)
    attach!(g, pn, rn => :feedline, 0mm, location=1)

    fp = plan(g)
    c = DeviceLayout.center(fp)
    sim_area = centered(Rectangle(4mm, 3.7mm), on_pt=c)
    render!(fp.coordinate_system, sim_area, LayerVocabulary.SIMULATED_AREA)
    render!(fp.coordinate_system, sim_area, LayerVocabulary.WRITEABLE_AREA)
    render!(fp.coordinate_system, centered(Rectangle(4mm, 3.7mm), on_pt=c),
            LayerVocabulary.CHIP_AREA)
    check!(fp)
    return fp
end

# ─── Render via EGADSKernel with tape ───────────────────────────────────────

function render_with_tape(fp; cap_length=620μm)
    # Bridge-free target (OCC can't dispatch egads_staple_bridge!)
    tech = SolidModelTarget(ExamplePDK.EXAMPLE_SINGLECHIP_TECHNOLOGY;
        simulation=true,
        bounding_layers=[:simulated_area],
        substrate_layers=[:chip_area],
        indexed_layers=[:port, :lumped_element, :integration, :wave_port],
        wave_port_layers=Symbol[],
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
        ],
        retained_physical_groups=[
            ("vacuum", 3), ("substrate", 3), ("metal", 2), ("exterior_boundary", 2)])

    sm = SolidModel("transmon_tape", EK(enable_tape=true))
    render!(sm, fp, tech)
    return sm
end

# ─── Main ───────────────────────────────────────────────────────────────────

println("Building floorplan...")
fp = build_floorplan()

println("Rendering via EGADSKernel with tape...")
sm = render_with_tape(fp)

println("\n=== Geometry Tape ===")
println(ext.tape_summary(sm))

# Count bodies
k = kernel(sm)
println("\n=== Bodies ===")
bodies = ext._bodies(k)
println("  Total bodies: $(length(bodies))")
for (i, b) in enumerate(bodies)
    info = EGADS.get_info(b)
    faces = EGADS.get_body_topos(b, EGADS.OClass.FACE)
    println("  Body $i: mtype=$(info.mtype), $(length(faces)) faces")
end
