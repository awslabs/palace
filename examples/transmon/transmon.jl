
import DeviceLayout
import JSON
import Gmsh: gmsh

const ST_PATH = joinpath(pkgdir(DeviceLayout), "examples", "SingleTransmon")

# Include SingleTransmon module from the DeviceLayout examples.
include(joinpath(ST_PATH, "SingleTransmon.jl"))

# Now `Main.SingleTransmon` is available.
"""

Physical details and properties can be customized by passing the keywords
accepted by `SingleTransmon.single_transmon()`. For instance, the increase the
number of meander turns to 7, call
```
generate_transmon(; n_meander_turns=7)
```

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
    config = SingleTransmon.configfile(solid_model, solver_order=1)

    # Fix paths (by default, they are within pkgdir(DeviceLayout)).
    config["Model"]["Mesh"] = joinpath("mesh", mesh_filename)
    config["Problem"]["Output"] = "postpro"
    # Tighten tolerances to improve reproduciability.
    config["Solver"]["Eigenmode"]["Tol"] = 1e-8
    config["Solver"]["Linear"]["Tol"] = 1e-10
    open(joinpath(config_path), "w") do f
        indent = 2
        return JSON.print(f, config, indent)
    end
end
