```@raw html
<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
```

```@setup json
using JSON
using OrderedCollections # Needed to preserve printing order as in the original file

# JSON.jl doesn't support // comments
strip_json_comments(json_str) = replace(json_str, r"//.*?($|\n)" => s"\1")

# JSON.jl expands lists on multiple lines, but we prefer preserving the attributes on a single line
fix_attributes_format(json_str) = replace(replace(json_str, r"\"Attributes\":\s*\[\s*(\d+)\s*\]" => s"\"Attributes\": [\1]"))

# Print just one section of the entire JSON file and prepend the correct header
print_section(json, section) = println("\"$section\":", fix_attributes_format(JSON.json(json[section], 4)))

# Move mesh and json file to build directory, so that we don't pollute the example directory
spheres_basepath = joinpath(@__DIR__, "..", "..", "examples", "spheres")
json_name = "spheres.json"
mesh_name = "spheres.msh"
json_basepath = joinpath(spheres_basepath, json_name)

Base.cp(json_basepath, json_name, force = true)
Base.mkpath("mesh")
Base.cp(joinpath(spheres_basepath, "mesh", mesh_name), joinpath("mesh", mesh_name), force = true)

json_text = strip_json_comments(read(json_basepath, String))
spheres_json = JSON.parse(json_text, dicttype=OrderedDict)

function logrun(cmd)
    try
        run(cmd)
    catch e
        println(e)
    end
    return nothing
end

```

# Run your first simulation with Palace

Welcome to your first tutorial with *Palace*!

In this tutorial, we will:

 1. [Install *Palace* using Spack](#Installing-Palace)
 2. [Set up a simulation using a provided mesh](#The-config-file)
 3. [Run the simulation and visualize results with ParaView](#Running-the-simulation-and-inspecting-the-output)

By the end of this page, you'll understand the basic workflow of electromagnetic
simulations with *Palace*. You will be able to follow the
[examples](examples/examples.md) and start experimenting with setting up and
running your own simulations.

As a motivating example, we'll solve for the electrostatic potential and
capacitance matrix of two conducting spheres inside a larger grounded sphere.
Our system consists of:

  - Two inner spherical conductors
  - A larger outer sphere at zero potential
  - Vacuum in between

This is [one of the examples](examples/spheres.md) included with *Palace*. For
more details on the physics and comparisons against analytical results, see the
[full example page](examples/spheres.md).

## A bird's-eye view of *Palace*

*Palace* is a finite-element code for **PA**rallel **LA**rge-scale
**C**omputational **E**lectromagnetics. It is a command-line executable that
runs on single or multiple nodes, supporting both [CPUs and
GPUs](guide/parallelism.md).

A *Palace* simulation requires two main inputs:

  - **Mesh file**: The mesh file describes the target geometry. *Palace* does not
    construct meshes, so you must supply one. A large number of formats are
    supported, allowing you to use your preferred CAD or meshing software. See the
    [supported mesh formats](guide/model.md#Supported-mesh-formats) for more
    details.

  - **Config file**: A JSON file that defines what problem to solve and how.
    It specifies material properties, boundary conditions, problem type,
    solver parameters, and output settings. See the [configuration
    documentation](config/config.md) for complete details.

*Palace* can solve Maxwell's equations under a few different set of assumptions,
leading to different electromagnetic [problem types](guide/problem.md):

  - **Electrostatic**: Computes static electric fields and capacitance matrices (used in this tutorial)
  - **Magnetostatic**: Computes static magnetic fields and inductance matrices
  - **Eigenmode**: Computes eigenmodes and eigenvalues in the frequency domain
  - **Driven**: Computes frequency-domain response to harmonic excitations
  - **Transient**: Computes time-domain response to boundary excitations

*Palace* produces two types of primary output:

 1. CSV files with post-processed quantities (capacitance matrices, field values
    at probe points, etc.)
 2. PVD files for visualizing fields with [ParaView](https://www.paraview.org)
    or compatible software

The full list of problem types and their outputs is available in the
[problem configuration guide](config/problem.md).

In this tutorial, we'll use a mesh generated with [Gmsh](https://gmsh.info/),
create a configuration file for an `Electrostatic` problem, and visualize the
resulting electric field with ParaView.

## Installing Palace

As a user, the simplest way to install *Palace* is with
[Spack](https://spack.readthedocs.io), a package manager designed for
high-performance computing applications.

Follow the [instructions on the official
website](https://spack.readthedocs.io/en/latest/getting_started.html) to install
Spack. This involves cloning a repository and sourcing a `setup-env` shell
script. Come back here once you are done with that.

Let's check that Spack is correctly installed on your system. This can be
accomplished by running:

```bash
spack --version
```

You should see your Spack version (yours may differ, that's okay):

```@example json
logrun(`spack --version`); # hide
nothing # hide
```

!!! note "spack command not found"
    
    If you get a `command not found` error, revisit the [Spack instructions](https://spack.readthedocs.io/en/latest/getting_started.html)
    and ensure you've completed all steps, including sourcing the setup script (the command starting with `.`).
    Consider adding this to your shell initialization file (typically, `.bashrc` or `.zshrc`).

With Spack installed, we can now move to *Palace*:

```bash
spack install palace
```

This will download and compile the latest release of *Palace* and its
dependencies. This step may take tens of minutes depending on your system.

Once installed, load *Palace* and verify it works:

```bash
spack load palace
palace --help
```

You should see:

```@example json
logrun(`palace --help`); # hide
nothing # hide
```

!!! tip "Loading Palace"
    
    You need to load *Palace* with `spack load palace` in each new shell session.
    For convenience, add this command to your shell initialization file if you are a frequent *Palace* user.

## (Optional) Install ParaView

*Palace* optionally saves electromagnetic field data in the PVD format, which is
immediately accessible by ParaView or ParaView-compatible software. You can
download ParaView from the [official
website](https://www.
paraview.org/download/) or using your package manager
(`dnf`, `apt`, `brew`, ...). ParaView is not required for running simulations,
but we will use it in this tutorial to visualize our simulated fields.

## The mesh

The mesh describes the geometry over which we want to solve the problem.
Constructing a mesh is not *Palace*'s responsibility and we'll use a pre-made
mesh from the *Palace* examples for this tutorial.

Create a `mesh` directory and download the [mesh
file](https://raw.githubusercontent.com/awslabs/palace/refs/heads/main/examples/spheres/mesh/spheres.msh):

```bash
mkdir -p mesh
curl -o mesh/spheres.msh https://raw.githubusercontent.com/awslabs/palace/refs/heads/main/examples/spheres/mesh/spheres.msh
```

This mesh was created using the [Julia](https://julialang.org/) interface for
Gmsh. If you're interested in how it was created, see the
[`mesh.jl`](https://raw.githubusercontent.com/awslabs/palace/refs/heads/main/examples/spheres/mesh/mesh.jl)
file.

### Understanding mesh attributes

To set up a simulation, you need to identify regions in the mesh (collections of
volumes or surfaces) to assign material properties and boundary conditions. The
mesh formats used within *Palace* all support this via a notion of `Attributes`.
An attribute is a 1-based index within a mesh that indicates a subset of
elements. The particulars of how these attributes are identified are specific to
each mesh format.

Our mesh follows the `msh2` format and contains four distinct regions:

 1. `domain` (3D volume, attribute: 1) - the vacuum region between spheres
 2. `farfield` (2D surface, attribute: 2) - the outer boundary surface
 3. `sphere_a` (2D surface, attribute: 3) - the first conductor surface
 4. `sphere_b` (2D surface, attribute: 4) - the second conductor surface

We'll reference these attributes in our configuration file to specify boundary
conditions and material properties.

## The config file

The configuration file defines the electromagnetic problem: what to solve for,
material properties and boundary conditions, the details of the algorithm to be
employed, and so on.

*Palace* config files contain five sections:

 1. **`Problem`**: Defines the physics type and output settings
 2. **`Model`**: Specifies the mesh file and geometric parameters
 3. **`Domains`**: Defines material properties for 3D regions and related postprocessing operations
 4. **`Boundaries`**: Sets boundary conditions for 2D surfaces and related postprocessing operations
 5. **`Solver`**: Controls numerical parameters and solution methods

Here's a complete configuration for our electrostatic problem:

```@raw html
<div style="max-height: 300px; overflow-y: auto;">
```

```@example json
println(fix_attributes_format(JSON.json(spheres_json, 4))) # hide
```

```@raw html
</div>
```

You can save this as `spheres.json` or download it directly:

```bash
curl -O https://raw.githubusercontent.com/awslabs/palace/refs/heads/main/examples/spheres/spheres.json
```

Let's examine each section in detail. For complete documentation on all available
options, see the [configuration reference](config/config.md).

If you wish to skip the explanation and jump directly to running your
simulations, go to [Running the simulation and inspecting the
output](#Running-the-simulation-and-inspecting-the-output).

#### Section 1: `Problem` definition

The `Problem` section identifies the problem type and the output directory. In
this case, we choose `Electrostatic`. This means that *Palace* solves Poisson's
equation for electric potential sequentially activating all the `Terminals` on
the mesh while setting the non-activated terminals to ground. All simulation
types in *Palace* have some form of iteration (over frequencies, times, mode
numbers, or terminals). The output is saved to the `"Output"` folder specified
in the `"Problem"` section in the JSON file, `postpro` in this example. If the
output already exists, it will be overwritten. See [`config["Problem"]`](config/problem.md)
for details on all available problem types and their outputs.

```@example json
print_section(spheres_json, "Problem") # hide
```

#### Section 2: `Model` specification

The `Model` section specifies the desired geometry. In addition to defining the
mesh, it specifies how to convert mesh units to physical units using the `L0`
parameter. For example, `L0` of `1e-2` means that one mesh unit corresponds to
one centimeter. The `Model` section can also include settings for adaptive mesh
refinement. See [`config["Model"]`](config/model.md) for more information.

```@example json
print_section(spheres_json, "Model") # hide
```

#### Section 3: `Domain` properties

The `Domains` section defines properties for the 3D regions in the geometry.

Each 3D region (identified by its `Attribute`) must have a `Material` definition
specifying its physical properties. In our mesh, we have just one 3D region (the
vacuum between the spheres and outer boundary) identified by attribute 1. While
vacuum properties are applied by default, you can specify various material
characteristics as detailed in
[`config["Domains"]["Materials"]`](config/domains.md#domains%5B%22Materials%22%5D).

The `Domains` section also includes a `Postprocessing` subsection for
calculating specific quantities. In this example, we add:

  - `Energy`, which activates integrated energy calculations in the 3D domain
  - `Probe`, which requests field values at specific coordinates defined by a `Center`
    (in mesh units)

When configuring `Postprocessing`, you must specify an `Index` that determines
the suffix for column headers in the output CSV files. For example, with `Index: 1`, the probe output will show headers like `E_x[1]`.

!!! note "What is the difference between `Attributes` and `Index`?"
    
    `Attributes` identify mesh regions and come from the mesh file. In our example,
    attributes 1-4 identify the vacuum region, outer boundary, and two spheres.
    
    `Index` is used only for postprocessing and defines a notation used in the output CSV files. It has no relation to mesh attributes and can
    be any positive integer.
    
    Note how `Attributes` is an array and `Index` an integer: multiple attributes might
    be needed to specify a given region in the mesh that corresponds to a single output.

```@example json
print_section(spheres_json, "Domains") # hide
```

#### Section 4: `Boundary` conditions

The `Boundaries` section maps 2D surfaces in your mesh to their physical
boundary conditions. *Palace* offers numerous boundary condition types, all
documented in [`config["Boundaries"]`](config/boundaries.md).

Unlike 3D regions, which all require `Material` specifications, 2D surfaces have
default behavior: any external surface without an explicit boundary condition is
treated as a Perfect Magnetic Conductor (PMC), where the tangential component of
the magnetic field is zero, and no conditions are imposed on internal surfaces
(since terms from either sides cancel out on such boundaries).

For our electrostatic problem, we define:

  - The outer boundary as `Ground` (zero potential)
  - Two `Terminal` surfaces (one for each sphere)

`Terminals` are particularly important for `Electrostatic` simulations. *Palace*
activates each terminal sequentially (applying a unit of potential and grounding
all the other ones) and solves Maxwell's equations. Each of these steps adds a
new row to the output CSV files.

Like the `Domains` section, `Boundaries` also includes a `Postprocessing`
subsection for calculating quantities such as surface fluxes across 2D regions.
Here, we compute the fluxes of electric fields across the spherical conductors.
See [`config["Boundaries"]`](config/boundaries.md) for all available
postprocessing options.

```@example json
print_section(spheres_json, "Boundaries") # hide
```

#### Section 5: Solver settings

Finally, the `Solver` section prescribes properties of the problem and the
numerical algorithm, what device to use for the solution, and how much to save
as PVD files. For this problem, we run on CPU, specify third-order finite
elements, and save the fields for both terminal activations. The details of the
linear solver parameters in `"Linear"` are not essential for this tutorial.

Other problem types typically have more extensive `Solver` configurations,
including excitation parameters and frequency sweep settings. For complete
details on all solver options, see [`config["Solver"]`](config/solver.md).

```@example json
print_section(spheres_json, "Solver") # hide
```

### Running the simulation and inspecting the output

If you've followed along, you should now have two files:

```
├── mesh
│   └── spheres.msh
└── spheres.json
```

!!! note "Do not have the files?"
    
    If you need to download the files, run:
    
    ```bash
    mkdir -p mesh
    curl -o mesh/spheres.msh https://raw.githubusercontent.com/awslabs/palace/refs/heads/main/examples/spheres/mesh/spheres.msh
    curl -O https://raw.githubusercontent.com/awslabs/palace/refs/heads/main/examples/spheres/spheres.json
    ```

Before running your simulation, it's a good idea to validate the configuration:

```bash
palace --dry-run spheres.json
```

This checks for syntax errors and basic configuration issues. The validator should return:

```@example json
logrun(`palace --dry-run spheres.json`); # hide
nothing # hide
```

Finally, we are ready to run the simulation:

```bash
palace -np 1 spheres.json
```

`-np 1` instruct *Palace* to run with a single MPI process.

You'll see output including mesh details, solver progress, and timing
information. The amount of information can be controlled with the `Verbose`
configuration option in the JSON file.

```@raw html
<div style="max-height: 300px; overflow-y: auto;">
```

```@example json
logrun(`palace -np 1 spheres.json`); # hide
nothing # hide
```

```@raw html
</div>
```

Notice that *Palace* ran two iterations, one for each `Terminal`. Different
problem types will have different iteration patterns (e.g., `Driven` iterating
over frequencies or `Transient` over time steps) and many of the output CSV
files are organized along these iterations.

### Understanding the output

Once the simulation is completed, you'll find a `postpro` directory containing:

```@example json
logrun(`ls -R postpro`); # hide
nothing # hide
```

In addition to the `palace.json`, which contains metadata about the simulation
(including timing information and counts), the output consists of CSV and PVD
files. You can safely ignore all the `Cycle` directories as their content is
accessed through the corresponding PVD file. For more details on output files
and formats, see the [output documentation](guide/postprocessing.md).

#### CSV files

CSV files contain post-processed quantities and depend on the specific problem
type chosen. Let's look at two examples:

`probe-E.csv` shows electric field values at the probe point we defined in the
`Postprocessing` section in `"Domains"`:

```@example json
logrun(`cat postpro/probe-E.csv`); # hide
nothing # hide
```

The first column `i` indicates the iteration and corresponds to the `Index`
associated to each `Terminal`, whereas the `[1]` in column headers corresponds
to the `Index` we specified in the `Probe` section.

One of the key outputs of the `Electrostatics` problem type is the capacitance
matrix, saved in `terminal-C.csv`:

```@example json
logrun(`cat postpro/terminal-C.csv`); # hide
nothing # hide
```

Here, both rows and columns correspond to `Terminal` indices. As expected, the
matrix is symmetric.

#### Visualizing with ParaView

In this final step, we'll create a visualization of our simulation results using
[ParaView](https://www.paraview.org). We'll work with both the volume field data
(`electrostatic.pvd`) and the boundary surface data
(`electrostatic_boundaries.pvd`) to reproduce the figures in [the example
page](examples/spheres.md).

 1. Launch ParaView and navigate to your `postpro/paraview` directory

 2. Open the volume data:
    
      + Click File → Open → Navigate to
        `postpro/paraview/electrostatic/electrostatic.pvd`, nothing should be
        rendered so far
      + Click the Apply button in the Properties panel (left side), a sphere should
        appear
      + In the Coloring section, select `V`, the sphere should now be colored
        according to the potential values
 3. Create a slice to see inside:
    
      + From the menu bar, select Filters → Common → Slice
      + In the Properties panel (left side), set the Origin to (0, 0, 0)
      + Set the Normal to (0, 1, 0) for a vertical slice along the Y-axis
      + Click Apply
      + Use the mouse to rotate and zoom until you can see the outlines of both inner spheres
 4. Add the boundary surfaces:
    
      + Click File → Open → Navigate to `postpro/paraview/electrostatic_boundaries/electrostatic_boundaries.pvd`
      + Click Apply in the Properties panel
      + In the Coloring section, select `V`
      + The two inner spheres should now appear with their surface potentials displayed

Notice the time slider at the top of the ParaView window:

  - Frame 0: First terminal activated (first sphere at unit potential)
  - Frame 1: Second terminal activated (second sphere at unit potential)
  - Frame 99: Error estimates

You can save the visualization as an image with File → Save Screenshot, or save
the entire ParaView state with File → Save State (allows reopening your complete
setup later). The result should look more or less like the images below:

```@raw html
<br/><p align="center">
  <img src="../assets/examples/spheres-3.png" width="45%" />
  <img src="../assets/examples/spheres-4.png" width="45%" />
</p>
```

ParaView offers many more advanced features for data analysis and visualization.
For more details, refer to the [official ParaView
documentation](https://docs.paraview.org/en/latest/).

## Where to go next

Congratulations! You've completed your first *Palace* simulation. To continue
learning:

 1. Try looking at the other output files in this simulation
 2. Try modifying this example with different materials or boundary conditions
 3. Explore the [examples](examples/examples.md) to see different problem types
    and more complex geometries
 4. Read the [configuration reference](config/config.md) to understand all
    available options

If you encounter any issues or have questions, please report them to our
[GitHub issue tracker](https://github.com/awslabs/palace/issues).
