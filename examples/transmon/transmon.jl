
import DeviceLayout
import JSON

const ST_PATH = joinpath(pkgdir(DeviceLayout), "examples", "SingleTransmon")

# Include SingleTransmon module from the DeviceLayout examples.
include(joinpath(ST_PATH, "SingleTransmon.jl"))

# Now `Main.SingleTransmon` is available.
"""

"""
function generate_transmon(;
    mesh_filename::AbstractString = "transmon.msh2",
    config_filename::AbstractString = "transmon.json",
    kwargs...
)
    solid_model = SingleTransmon.single_transmon(; save_mesh = true, kwargs...)
    mesh_path = joinpath(ST_PATH, "single_transmon.msh2")
    new_mesh_path = joinpath(@__DIR__, "mesh", mesh_filename)
    Base.mv(mesh_path, new_mesh_path, force = true)

    config_path = joinpath(@__DIR__, config_filename)

    # solver_order=1 to keep the simulation fast enough for CI.
    config = SingleTransmon.configfile(solid_model, solver_order=1, amr=1)

    # Fix paths (by default, they are within pkgdir(DeviceLayout)).
    config["Model"]["Mesh"] = new_mesh_path
    config["Problem"]["Output"] = "postpro"
    # Reduce verbosity.
    config["Problem"]["Verbose"] = 1
    open(joinpath(config_path), "w") do f
        indent = 2
        return JSON.print(f, config, indent)
    end
end
