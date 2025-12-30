# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script uses DeviceLayout.jl to create a mesh and Palace configuration for a
superconducting transmon qubit coupled to a readout resonator.

The generated system contains:
1. A transmon qubit (superconducting capacitor pads connected by a Josephson junction)
2. A quarter-wave coplanar resonator (meandering transmission line for readout)
3. A feedline (straight coplanar waveguide for signal input/output)
4. Lumped ports (boundary conditions for external connections and Josephson junction)

## Prerequisites

This script requires DeviceLayout.jl and its dependencies. If you don't already have it
installed, you can install it with (from this directory)

```bash
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## How to run

From this directory, run:
```bash
julia --project -e 'include("transmon.jl"); generate_transmon()'
```
This generates the mesh and configuration files used in the example.

To customize the device geometry, pass parameters to the function:
```bash
julia --project -e 'include("transmon.jl"); generate_transmon(cap_length=800μm, n_meander_turns=7)'
```

To generate a plot of the mesh for visualization (used for the documentation):
```bash
julia --project -e 'include("transmon.jl"); plot_transmon_mesh()'
```

The script will generate mesh and configuration files and print the output file locations.
These files are ready for use with Palace simulations.
=#

import DeviceLayout
import JSON
import Gmsh: gmsh

const ST_PATH = joinpath(pkgdir(DeviceLayout), "examples", "SingleTransmon")

# Include SingleTransmon module from the DeviceLayout examples.
include(joinpath(ST_PATH, "SingleTransmon.jl"))
# Now the `SingleTransmon` module is available.

"""
    generate_transmon(;
        mesh_filename::AbstractString = "transmon.msh2",
        config_filename::AbstractString = "transmon.json",
        solver_order = 2,
        kwargs...
    )

Generate a mesh file and Palace configuration for a transmon-resonator system.

This function creates a complete simulation setup by:

 1. Generating the 3D geometry and mesh using DeviceLayout.jl
 2. Creating the corresponding Palace configuration file

# Arguments

  - `mesh_filename`: Name of the output mesh file (default: "transmon.msh2")
  - `config_filename`: Name of the output configuration file (default: "transmon.json")
  - `solver_order`: Finite element order (1 or 2). Higher order gives better accuracy
    but increases computational cost significantly (default: 2)

# Keyword Arguments

All keyword arguments are passed to `SingleTransmon.single_transmon()` to
customize the device geometry. Parameters include:

  - `cap_length=620μm`: Length of transmon capacitor pads
  - `cap_width=24μm`: Width of transmon capacitor pads
  - `cap_gap=30μm`: Gap between transmon capacitor pads
  - `total_length=5000μm`: Total electrical length of readout resonator
  - `n_meander_turns=5`: Number of meander turns in resonator
  - `claw_gap=6μm`: Gap between resonator claw and transmon
  - `w_claw=34μm`: Width of resonator claw fingers
  - `l_claw=121μm`: Length of resonator claw fingers

For a full list of parameters, consult the keyword arguments for
`SingleTransmon.single_transmon()`.

# Output Files

  - `mesh/transmon.msh2`: Gmsh mesh file for Palace
  - `transmon.json`: Palace configuration file
"""
function generate_transmon(;
    mesh_filename::AbstractString="transmon.msh2",
    config_filename::AbstractString="transmon.json",
    solver_order=2,
    kwargs...
)
    # Generate the solid model and mesh using DeviceLayout.jl.
    solid_model = SingleTransmon.single_transmon(; save_mesh=true, kwargs...)

    # Move mesh file to local mesh directory.
    mesh_path = joinpath(ST_PATH, "single_transmon.msh2")
    new_mesh_path = joinpath(@__DIR__, "mesh", mesh_filename)
    Base.mv(mesh_path, new_mesh_path, force=true)

    # Generate Palace configuration.
    config_path = joinpath(@__DIR__, config_filename)
    config = SingleTransmon.configfile(solid_model; solver_order)

    # Update paths for local directory structure.
    config["Model"]["Mesh"] = joinpath("mesh", mesh_filename)
    config["Problem"]["Output"] = "postpro"

    # Set tight tolerances for reproducible results (important for CI/testing).
    config["Solver"]["Eigenmode"]["Tol"] = 1e-8
    config["Solver"]["Linear"]["Tol"] = 1e-12

    # Get closer to expected values to speed up convergence.
    config["Solver"]["Eigenmode"]["Target"] = 4

    # Enable GridFunction output for GLVis visualization (Used in the
    # documentation).
    config["Problem"]["OutputFormats"] = Dict("GridFunction" => true)

    # Write configuration file with pretty formatting.
    open(config_path, "w") do f
        indent = 2
        JSON.print(f, config, indent)
        println(f)
        return nothing
    end

    println("Generated files:")
    println("  Mesh: $(new_mesh_path)")
    println("  Config: $(config_path)")

    return nothing
end

"""
    plot_transmon_mesh(; 
                         mesh_file = "mesh/transmon.msh2",
                         filename = "transmon-1.svg"
                         )

Generate a plot of the transmon mesh to use in the documentation.

This function opens the generated mesh file in Gmsh, configures the
visualization settings to highlight the device structure, adds annotations for
key components, and exports the result as an SVG file. The output is used in the
Palace documentation.

You may convert the output to PNG using `imagemagick`.

# Requirements

Requires a mesh file at `mesh/transmon.msh2`. Run `generate_transmon()` first
if the mesh file doesn't exist.
"""
function plot_transmon_mesh(; mesh_file="mesh/transmon.msh2", filename="transmon-1.svg")
    # Initialize Gmsh in batch mode (no interactive GUI).
    gmsh.initialize(["-batch"])
    gmsh.fltk.initialize()

    # Load the mesh file.
    if !isfile(mesh_file)
        error("Mesh file not found: $mesh_file. Run generate_transmon() first.")
    end

    # Check that filename has .svg extension.
    if !endswith(filename, ".svg")
        error("Output filename must have .svg extension, got: $filename")
    end

    gmsh.open(mesh_file)

    # Configure visualization to show only surface mesh.
    gmsh.option.setNumber("Mesh.SurfaceEdges", 1)  # Show surface edges.
    gmsh.option.setNumber("Mesh.SurfaceFaces", 0)  # Hide surface faces.
    gmsh.option.setNumber("Mesh.VolumeEdges", 0)   # Hide volume edges.
    gmsh.option.setNumber("Mesh.VolumeFaces", 0)   # Hide volume faces.

    # Hide coordinate axes for cleaner appearance.
    gmsh.option.setNumber("General.Axes", 0)
    gmsh.option.setNumber("General.SmallAxes", 0)
    gmsh.option.setNumber("General.AxesAutoPosition", 0)
    gmsh.option.setNumber("General.GraphicsWidth", 1600)
    gmsh.option.setNumber("General.GraphicsHeight", 800)

    # Add text annotations for key components.
    gmsh.view.add("Annotations")

    # Component labels with positions and font sizes.
    annotations = [
        ("Transmon", 725.0, 525.0, 30),
        ("Read-out line", 680.0, 210.0, 30),
        ("Coplanar resonator", 770.0, 400.0, 30),
        ("Port", 500.0, 225.0, 26),
        ("Port", 910.0, 225.0, 26),
        ("Port", 700.0, 580.0, 26)
    ]

    for (text, x, y, size) in annotations
        gmsh.plugin.setString("Annotate", "Text", text)
        gmsh.plugin.setNumber("Annotate", "FontSize", size)
        gmsh.plugin.setNumber("Annotate", "X", x)
        gmsh.plugin.setNumber("Annotate", "Y", y)
        gmsh.plugin.setNumber("Annotate", "View", 0)
        gmsh.plugin.run("Annotate")
    end

    # Hide exterior boundary for cleaner visualization.
    physical_groups = gmsh.model.getPhysicalGroups(2)
    for (dim, tag) in physical_groups
        name = gmsh.model.getPhysicalName(dim, tag)
        if name == "exterior_boundary"
            entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
            gmsh.model.setVisibility([(dim, e) for e in entities], 0)
        end
    end

    # NOTE: If there's an easier way to customize a specific color and to crop
    # the whitespace from the figure directly in Gmsh, @gbozzola could not find
    # it. But if you do, please fix this up.

    # Export initial SVG.
    gmsh.write(filename)

    # Post-process SVG.
    _postprocess_svg(filename)

    # Clean up Gmsh.
    gmsh.finalize()

    return println("Generated mesh plot: $filename")
end

"""
    _postprocess_svg(filename::String)

Post-process the exported SVG file to crop the image to the content and increase
contrast of the lines.
"""
function _postprocess_svg(filename::String)
    svg_content = read(filename, String)

    # Change stroke color for better contrast (blue instead of default).
    svg_content = replace(
        svg_content,
        r"stroke=\"#[0-9a-f]{2}ff[0-9a-f]{2}\"" => "stroke=\"#27B4F5\""
    )

    # Extract coordinates to determine bounding box for cropping.
    coords = Float64[]
    for m in eachmatch(r"points=\"([^\"]+)\"", svg_content)
        for pair in split(m.captures[1], r"\s+")
            if occursin(',', pair)
                x, y = split(pair, ',')
                push!(coords, parse(Float64, x), parse(Float64, y))
            end
        end
    end

    # Crop to geometry bounds with padding.
    if !isempty(coords)
        xs = coords[1:2:end]
        ys = coords[2:2:end]

        min_x, max_x = minimum(xs), maximum(xs)
        min_y, max_y = minimum(ys), maximum(ys)

        # Add 5% padding around geometry.
        padding = 0.05 * max(max_x - min_x, max_y - min_y)
        viewbox = "$(min_x - padding) $(min_y - padding) $(max_x - min_x + 2*padding) $(max_y - min_y + 2*padding)"
        svg_content = replace(svg_content, r"viewBox=\"[^\"]*\"" => "viewBox=\"$viewbox\"")
    end

    # Write processed SVG.
    return write(filename, svg_content)
end
