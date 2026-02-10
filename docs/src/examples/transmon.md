```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

```@setup include_example
function include_example_file(example_path, filename, max_lines = 4)
    content = read(joinpath(@__DIR__, "..", "..", "..", "test", "examples", "ref", example_path, filename), String)
    lines = split(content, '\n')
    println(join(lines[1:min(max_lines, length(lines))], '\n'))
end
```

# Single transmon with read-out resonator from [`DeviceLayout.jl`](https://aws-cqc.github.io/DeviceLayout.jl/)

!!! note

    The files for this example can be found in the
    [`examples/transmon/`](https://github.com/awslabs/palace/blob/main/examples/transmon)
    directory of the *Palace* source code. The configuration and mesh files
    were generated with [`DeviceLayout.jl`](https://aws-cqc.github.io/DeviceLayout.jl).

In this example, we simulate a superconducting
[transmon](https://en.wikipedia.org/wiki/Transmon) qubit coupled to a readout
resonator, a fundamental building block of superconducting quantum processors.

The configuration and mesh for this example are generated with
[`DeviceLayout.jl`](https://aws-cqc.github.io/DeviceLayout.jl), a Julia package
for computer-aided design of quantum integrated circuits. DeviceLayout.jl
supports the generation of 2D layouts and 3D models using both low-level
geometry interfaces and high-level schematic-driven workflows, enabling
automated design and optimization of quantum devices. For a detailed tutorial on
using DeviceLayout.jl to generate this geometry, see [this
example](https://aws-cqc.github.io/DeviceLayout.jl/stable/examples/singletransmon/).

### Overview and Configuration

The system consists of three main components: a transmon qubit, a quarter-wave
coplanar waveguide resonator, and a feedline for signal input and output. The
transmon is capacitively coupled to the resonator through a "claw" structure
that extends around the transmon cutout, while the resonator couples to the
feedline to enable readout of the qubit state.

The image below shows the meshed geometry of the system. The metal structures
are fabricated on a dielectric substrate and enclosed within a simulation box
filled with vacuum to model the actual device environment.

```@raw html
<!-- Plot generated with plot_transmon_mesh() in transmon/transmon.jl -->
<br/><p align="center">
  <img src="../../assets/examples/transmon-1.png" width="100%" />
</p><br/>
```

The metal conductors are modeled as perfect electric conductors (PEC), with
three lumped ports defining the electromagnetic boundary conditions. Two
resistive ports at the feedline ends represent 50-Ohm terminations to the
external world, while a third port models the Josephson junction as a linearized
LC circuit with both inductance and capacitance. The substrate is sapphire,
chosen for its low dielectric loss and well-characterized anisotropic material
properties (see,
[`config["Domains"]["Materials"]`](../config/domains.md#domains%5B%22Materials%22%5D)
for more information on how to configure anisotropic materials).
[`config["Boundaries"]["Absorbing"]`](../config/boundaries.md#boundaries%5B%22Absorbing%22%5D)
boundary conditions on the external boundaries of the simulation box capture electromagnetic
radiation and prevent artificial reflections.

The JSON configuration file and mesh are generated automatically by the
[`SingleTransmon`
module](https://aws-cqc.github.io/DeviceLayout.jl/stable/examples/singletransmon/)
in DeviceLayout.jl. This module demonstrates how to programmatically assign
domain and boundary attributes from a schematic design, enabling semi-automated
configuration generation for electromagnetic simulations.

The system design is fully parametric, allowing customization of all geometric
and material properties. Physical details can be modified by passing keyword
arguments to the `generate_transmon` function in the `transmon.jl` file. To see
the available design parameters, run the following command in the transmon
example directory:

```bash
julia --project -E 'include("transmon.jl"); @doc SingleTransmon.single_transmon'
```

This displays the available design parameters, including geometric dimensions
and processing options:

```julia
single_transmon(
    w_shield=2μm,
    claw_gap=6μm,
    w_claw=34μm,
    l_claw=121μm,
    cap_width=24μm,
    cap_length=620μm,
    cap_gap=30μm,
    n_meander_turns=5,
    hanger_length=500μm,
    bend_radius=50μm,
    save_mesh::Bool=false,
    save_gds::Bool=false
)
```

Beyond the physical dimensions, the `save_mesh` and `save_gds` options enable
output of a mesh file for simulation and a [GDS
file](https://en.wikipedia.org/wiki/GDSII) for fabrication, respectively. This
integration between design, simulation, and fabrication workflows demonstrates
the power of automated quantum device design tools.

We perform an eigenmode analysis to determine the system's resonant frequencies
and quality factors. This simulation identifies the two fundamental modes
corresponding to the transmon and readout resonator, providing essential
information for device characterization including mode frequencies, decay rates,
and electromagnetic field distributions. Such analysis is crucial in the quantum
device design process, as it informs design decisions and enables optimization
before fabrication. The DeviceLayout.jl
[example](https://aws-cqc.github.io/DeviceLayout.jl/stable/examples/singletransmon/)
extends this approach by demonstrating an automated optimization loop that tunes
device parameters to achieve target frequencies.

### Results

The simulation completes in approximately 10 minutes on a modern computer using
6 CPU cores.

The primary results are the eigenfrequencies and quality factors, stored in the
`eig.csv` file:

```@example include_example
include_example_file("transmon/transmon_coarse", "eig.csv") #hide
```

The columns represent the mode number, real and imaginary parts of the
eigenfrequency, and computed quality factor, respectively. The results
successfully identify two distinct eigenmodes with predominantly real
frequencies and small imaginary components representing dissipation.

To verify the physical nature of these modes, we can examine the electromagnetic
field distributions. Below, we embed an interactive visualization of the
magnetic energy density of the second eigenmode. This visualization leverages
[GLVis](https://glvis.org). For a quick introduction to this tool, refer to the
box below. We encourage you to play with the visualization and inspect the mesh,
looking at how `DeviceLayout` refines certain regions of interest (e.g., the
claw).

!!! tip "Quick Start to GLVis"

    [GLVis](https://glvis.org) provides fast visualization of Palace output when
    [`config["Problem"]["Output"]["OutputFormats"]["GridFunction"]`](../config/problem.md#problem%5B%OutputFormats%22%5D) is enabled.
    While less quantitative than [ParaView](https://www.paraview.org/download/), GLVis excels at rapid solution inspection and is particularly well-suited
    for interactive exploration. The interface is primarily keyboard-driven with extensive
    [keybindings](https://github.com/glvis/glvis/blob/v4.4/README.md).

    Essential commands include:

      - `i` cut through the solution
      - `x/X`, `y/Y` rotate the cutting plane
      - `z/Z` translate the cutting plane
      - `R` cycle through 2D viewpoints
      - `c` display the colorbar
      - `m` show the mesh
      - `right click + mouse movement` zoom in/out
      - `center click + mouse movement` translate viewpoint

```@raw html
<div id="glvis-container" style="width: 100%; height: 500px;">
  <div id="glvis-div" style="width: 100%; height: 100%;"></div>
</div>

<!-- This was generated with order=1 and without mesh cracking. -->
<!-- We disable mesh cracking because it looks glitchy otherwise. -->

<!-- Note, the snippet below only works with one MPI process because we are
    manually composing a stream file. -->
<script type="text/javascript">
  var div = document.getElementById("glvis-div");
  require(["../../assets/js/glvis/index.js"], function (glvis) {
    var glv = new glvis.State(div);

    Promise.all([
      fetch('../../assets/examples/transmon-mesh.000000').then(r => r.text()),
      fetch('../../assets/examples/transmon-U_m_000002.gf.000000').then(r => r.text())
    ]).then(function(results) {
      var stream = "solution\n" + results[0] + results[1] + "keys OOOOOOOOOiXXXXXXXXXXXXXXXXXX1Lac\n";
      var originalTitle = document.title;
      glv.display(stream).then(function() {
        document.title = originalTitle;
        // Enable keyboard controls on click
        div.addEventListener('click', function() {
          this.setAttribute('tabindex', '0');
          this.focus();
        }, { once: true });
      });
    }).catch(function(e) {
      console.error('Failed to load GLVis data:', e);
    });
  });
</script>
```

The ParaView visualization below shows the magnetic energy density distribution,
providing clear confirmation that the computed modes correspond to the expected
transmon and resonator physics.

```@raw html
<!-- Plot generated with transmon/plot_slice_docs.py -->
<br/><p align="center">
  <img src="../../assets/examples/transmon-2.png" width="100%" />
</p><br/>
```

Additional simulation outputs provide detailed characterization of the system.
The `port-Q.csv` file contains port coupling rates and quality factors for the
resistive terminations at the feedline ends, while `port-EPR.csv` reports energy
participation ratios for inductive lumped elements (the Josephson junction). The
`domain-E.csv` file provides bulk dielectric loss and participation ratios for
the sapphire substrate, enabling comprehensive analysis of loss mechanisms
throughout the device.

#### Balancing Speed and Accuracy: Solver Order Considerations and

While this simulation provides accurate results, computational time can be
significantly reduced by adjusting the finite element order. Changing from
second order to first order decreases runtime by a factor of 30, though this
comes at the cost of reduced accuracy. This trade-off can be controlled by
modifying the solver order in the JSON configuration file or by passing the
`solver_order` parameter to the `generate_transmon` function.

The `error-indicators.csv` file provides quantitative measures of numerical
accuracy (see, the reference section on [error
estimation](../reference.md#Error-estimation-and-adaptive-mesh-refinement-(AMR))).
For a second-order solution:

```@example include_example
include_example_file("transmon/transmon_coarse", "error-indicators.csv") #hide
```

These error indicators provide metrics for estimating numerical accuracy across
the computational domain, with the `Norm` column representing a measure related
to the H(curl) seminorm of the solution error. When reducing to first order, the
error norm increases by a factor of two, and the computed eigenfrequencies
deviate noticeably from the converged values. While first-order simulations are
valuable for rapid prototyping and debugging, second order (or higher) is
generally recommended for accurate results.

#### Adaptive Mesh Refinement

Adaptive Mesh Refinement (AMR) automatically refines the computational mesh in
regions where the solution error is large, providing an efficient way to improve
accuracy without uniformly refining the entire domain. Palace uses a
flux-recovery [error
estimator](../reference.md#Error-estimation-and-adaptive-mesh-refinement-(AMR))
to identify elements that contribute most to the global error.

AMR works by:

 1. Solving the problem on the current mesh
 2. Computing element-wise error indicators
 3. Marking elements with large errors for refinement
 4. Refining the mesh
 5. Repeating

AMR is particularly effective for problems with localized features where uniform
refinement would be wasteful.

In *Palace*, AMR continues until one or more configurable conditions are met:

  - `Tol`: Error norm drops below tolerance (the `Norm` column in `error-indicators.csv`)
  - `MaxIts`: Maximum number of refinement iterations (default: `0`, disabled)
  - `MaxSize`: Maximum degrees of freedom (default: `0`, unlimited)

See
[`config["Model"]["Refinement"]`](../config/model.md#config%5B%22Model%22%5D)
for the various AMR parameters.

For a more realistic simulation of a transmon, we enable AMR with:

```@example
import JSON #hide
json_amr = JSON.parsefile(  #hide
    joinpath(@__DIR__, "..", "..", "..", "examples", "transmon", "transmon_amr.json")  #hide
) #hide
println(  #hide
    "config[\"Model\"][\"Refinement\"][\"MaxIts\"] = $(json_amr["Model"]["Refinement"]["MaxIts"])"  #hide
)  #hide
println("config[\"Solver\"][\"Order\"] = $(json_amr["Solver"]["Order"])") #hide
```

After running the simulation, the `postpro` folder contains an `iterationN`
subfolder for each AMR cycle. The convergence history shows the error decreasing
with each refinement:

```@example
import JSON #hide
# Path to AMR simulation results #hide
amr_dir = joinpath(@__DIR__, "..", "..", "..", "examples", "transmon", "postpro", "amr") #hide
# Find all iteration folders and sort them #hide
iter_folders = filter(x -> startswith(x, "iteration"), readdir(amr_dir)) #hide
iter_numbers = [parse(Int, replace(f, "iteration" => "")) for f in iter_folders] #hide
sorted_indices = sortperm(iter_numbers) #hide
iter_folders = iter_folders[sorted_indices] #hide
# Collect data for each iteration #hide
all_dofs = Int[] #hide
all_errors = Float64[] #hide
for iter_name in [iter_folders; ""] #hide
    iter_dir = isempty(iter_name) ? amr_dir : joinpath(amr_dir, iter_name) #hide
    palace_json = JSON.parsefile(joinpath(iter_dir, "palace.json")) #hide
    push!(all_dofs, palace_json["Problem"]["DegreesOfFreedom"]) #hide
    error_file = joinpath(iter_dir, "error-indicators.csv") #hide
    lines = readlines(error_file) #hide
    push!(all_errors, parse(Float64, strip(split(lines[2], ',')[1]))) #hide
end #hide
# Print convergence table #hide
for idx = 1:length(all_dofs) #hide
    println(  #hide
        "iteration = $(idx-1), DOFs = $(all_dofs[idx]), error = $(round(all_errors[idx], sigdigits=4))"  #hide
    ) #hide
end #hide
```
