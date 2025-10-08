```@raw html
<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
```

# `config["Problem"]`

```json
"Problem":
{
    "Type": <string>,
    "Verbose": <int>,
    "Output": <string>,
    "OutputFormats":
    {
      ...
    }
}
```

with

`"Type" [None]` :  Controls the simulation type. The available options are:

  - `"Eigenmode"` :  Perform a undamped or damped eigenfrequency analysis.
  - `"Driven"` :  Perform a frequency response simulation.
  - `"Transient"` :  Perform a time domain excitation response simulation.
  - `"Electrostatic"` :  Perform an electrostatic analysis to compute the capacitance matrix
    for a set of voltage terminals.
  - `"Magnetostatic"` :  Perform a magnetostatic analysis to compute the inductance matrix
    for a set of current sources.

`"Verbose" [1]` :  Controls the level of log file printing.

`"Output" [None]` :  Directory path for saving postprocessing outputs.

`"OutputFormats"` :  Top-level object for configuring the field output formats.

## `problem["OutputFormats"]`

```json
"OutputFormats":
{
    "Paraview": <bool>,
    "GridFunction": <bool>
}
```

with

`"Paraview" [true]` :  Set to true to output fields in Paraview format.

`"GridFunction" [false]` :  Set to true to output fields in MFEM grid function format to visualize with GLVis.
