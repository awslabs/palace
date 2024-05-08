```@raw html
<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
```

# `config["Boundaries"]`

```json
"Boundaries":
{
    "PEC":
    {
        ...
    },
    "PMC":
    {
        ...
    },
    "Impedance":
    [
        ...
    ],
    "Absorbing":
    {
        ...
    },
    "Conductivity":
    [
        ...
    ],
    "LumpedPort":
    [
        ...
    ],
    "WavePort":
    [
        ...
    ],
    "WavePortPEC":
    {
        ...
    },
    "SurfaceCurrent":
    [
        ...
    ],
    "Ground":
    {
        ...
    },
    "ZeroCharge":
    {
        ...
    },
    "Terminal":
    [
        ...
    ],
    "Postprocessing":
    {
        "Side": <string>,
        "SurfaceFlux":
        [
            ...
        ],
        "Dielectric":
        [
            ...
        ]
    }
}
```

with

`"PEC"` :  Top-level object for configuring perfect electric conductor (PEC) boundary
conditions (zero tangential electric field).

`"PMC"` :  Top-level object for configuring perfect magnetic conductor (PMC) boundary
conditions (zero tangential magnetic field). Also imposes symmetry of the electric field
across the boundary surface.

`"Impedance"` :  Array of objects for configuring surface impedance boundary conditions. A
surface impedance boundary relates the tangential electric and magnetic fields on the
boundary using a user specified surface impedance.

`"Absorbing"` : Top-level object for configuring absorbing boundary conditions. These are
artificial scattering boundary conditions at farfield boundaries.

`"Conductivity"` :  Array of objects for configuring finite conductivity surface impedance
boundary conditions. Finite conductivity boundaries are only available for the frequency
domain driven simulation type.

`"LumpedPort"` :  Array of objects for configuring lumped port boundary conditions. Lumped
ports can be specified on boundaries which are internal to the computational domain.

`"WavePort"` :  Array of objects for configuring numeric wave port boundary conditions. Wave
ports can only be specified on boundaries which are on the true boundary of the
computational domain. Additionally, wave port boundaries are only available for the
frequency domain driven simulation type.

`"WavePortPEC"` :  Top-level object for configuring PEC boundary conditions for boundary
mode analysis performed on the wave port boundaries. Thus, this object is only relevant
when wave port boundaries are specified under
[`config["Boundaries"]["WavePort"]`](#boundaries%5B%22WavePort%22%5D).

`"SurfaceCurrent"` :  Array of objects for configuring surface current boundary conditions.
This boundary prescribes a unit source surface current excitation on the given boundary in
order to excite a frequency or time domain driven simulation or magnetostatic simulation.
For the magnetostatic simulation type, entries of the inductance matrix are extracted
corresponding to each surface current boundary.

`"Ground"` :  Top-level object for specifying ground, or zero voltage, boundary conditions
for for electrostatic simulations.

`"ZeroCharge"` :  Top-level object for specifying zero charge boundary conditions for for
electrostatic simulations. Also imposes symmetry of the electric field across the boundary
surface.

`"Terminal"` :  Array of objects for configuring terminal boundary conditions for
electrostatic simulations. Entries of the capacitance matrix are extracted corresponding to
each terminal boundary.

`"Postprocessing"` :  Top-level object for configuring boundary postprocessing.

`"Side" ["SmallerRefractiveIndex"]` :  Defines the postprocessing side for internal
boundary surfaces where the fields are in general double-valued. This is only relevant for
output for [boundary visualization with ParaView](../guide/postprocessing.md#Visualization).
The available options are:

  - `"SmallerRefractiveIndex"` :  Take the value from the side where the material index of
    refraction is smaller (speed of light is larger). Typically this selects the vacuum
    side. For anisotropic materials, the index of refraction associated with the principal
    direction with the smallest value is used.
  - `"LargerRefractiveIndex"` :  Take the value from the side where the material index of
    refraction is larger (speed of light is smaller). Typically this selects the non-vacuum
    side.

`"SurfaceFlux"` :  Array of objects for postprocessing surface flux.

`"Dielectric"` :  Array of objects for postprocessing surface interface dielectric loss.

## `boundaries["PEC"]`

```json
"PEC":
{
    "Attributes": [<int array>]
}
```

with

`"Attributes" [None]` :  Integer array of mesh boundary attributes at which to apply the PEC
boundary condition.

## `boundaries["PMC"]`

```json
"PMC":
{
    "Attributes": [<int array>]
}
```

with

`"Attributes" [None]` :  Integer array of mesh boundary attributes at which to apply the
PMC boundary condition.

## `boundaries["Impedance"]`

```json
"Impedance":
[
    {
        "Attributes": [<int array>],
        "Rs": <float>,
        "Ls": <float>,
        "Cs": <float>
    },
    ...
]
```

with

`"Attributes" [None]` :  Integer array of mesh boundary attributes for this surface
impedance boundary.

`"Rs" [0.0]` :  Surface resistance used for computing this surface impedance boundary's
impedance per square, ``\Omega``/sq.

`"Ls" [0.0]` :  Surface inductance used for computing this surface impedance boundary's
impedance per square, H/sq.

`"Cs" [0.0]` :  Surface capacitance used computing this surface impedance boundary's
impedance per square, F/sq.

## `boundaries["Absorbing"]`

```json
"Absorbing":
{
    "Attributes": [<int array>],
    "Order": <int>
}
```

with

`"Attributes" [None]` :  Integer array of mesh boundary attributes at which to apply
farfield absorbing boundary conditions.

`"Order" [1]` :  Specify a first- or second-order approximation for the farfield absorbing
boundary condition. Second-order absorbing boundary conditions are only available for the
frequency domain driven simulation type.

## `boundaries["Conductivity"]`

```json
"Conductivity":
[
    {
        "Attributes": [<int array>],
        "Conductivity": <float>,
        "Permeability": <float>,
        "Thickness": <float>
    },
    ...
]
```

with

`"Attributes" [None]` :  Integer array of mesh boundary attributes for this finite
conductivity boundary.

`"Conductivity" [None]` :  Electrical conductivity for this finite conductivity boundary,
S/m.

`"Permeability" [1.0]` :  Relative permeability for this finite conductivity boundary.

`"Thickness" [None]` :  Optional conductor thickness for this finite conductivity boundary
specified in mesh length units. Activates a finite conductivity boundary condition which
accounts for nonzero metal thickness.

## `boundaries["LumpedPort"]`

```json
"LumpedPort":
[
    {
        "Index": <int>,
        "Attributes": [<int array>],
        "Direction": <string> or [<float array>],
        "CoordinateSystem": <string>,
        "Excitation": <bool>,
        "Active": <bool>,
        "R": <float>,
        "L": <float>,
        "C": <float>,
        "Rs": <float>,
        "Ls": <float>,
        "Cs": <float>,
        "Elements":
        [
            {
                "Attributes": <string>,
                "Direction": <string> or [<float array>],
                "CoordinateSystem": <string>
            },
            ...
        ]
    },
    ...
]
```

with

`"Index" [None]` :  Index of this lumped port, used in postprocessing output files.

`"Attributes" [None]` :  Integer array of mesh boundary attributes for this lumped port
boundary. If this port is to be a multielement lumped port with more than a single lumped
element, use the `"Elements"` array described below.

`"Direction" [None]` :  Direction to define the polarization direction of the port field
mode on this lumped port boundary. Axis aligned lumped ports can be specified using
keywords: `"+X"`, `"-X"`, `"+Y"`, `"-Y"`, `"+Z"`, `"-Z"`, while coaxial lumped ports can be
specified using `"+R"`, `"-R"`. The direction can alternatively be specified as a
normalized array of three values, for example `[0.0, 1.0, 0.0]`. If a vector direction is
specified, the `"CoordinateSystem"` value specifies the coordinate system it is expressed
in. If this port is to be a multielement lumped port with more than a single lumped
element, use the `"Elements"` array described below.

`"CoordinateSystem" ["Cartesian"]` : Coordinate system used to express the `"Direction"`
vector, the options are `"Cartesian"` and `"Cylindrical"`. If a keyword argument is used
for `"Direction"` this value is ignored, and the appropriate coordinate system is used
instead.

`"Excitation" [false]` :  Turns on or off port excitation for this lumped port boundary for
driven or transient simulation types.

`"Active" [true]` :  Turns on or off damping boundary condition for this lumped port
boundary for driven or transient simulation types.

`"R" [0.0]` :  Circuit resistance used for computing this lumped port boundary's impedance,
``\Omega``. This option should only be used along with the corresponding `"L"` and `"C"`
parameters, and not with any of the surface parameters `"Rs"`, `"Ls"`, or `"Cs"`.

`"L" [0.0]` :  Circuit inductance used for computing this lumped port boundary's impedance,
H. This option should only be used along with the corresponding `"R"` and `"C"` parameters,
and not with any of the surface parameters `"Rs"`, `"Ls"`, or `"Cs"`.

`"C" [0.0]` :  Circuit capacitance used for computing this lumped port boundary's impedance,
F. This option should only be used along with the corresponding `"R"` and `"L"` parameters,
and not with any of the surface parameters `"Rs"`, `"Ls"`, or `"Cs"`.

`"Rs" [0.0]` :  Surface resistance used for computing this lumped port boundary's impedance,
``\Omega``/sq. This option should only be used along with the corresponding `"Ls"` and
`"Cs"` parameters, and not with any of the circuit parameters `"R"`, `"L"`, or `"C"`.

`"Ls" [0.0]` :  Surface inductance used for computing this lumped port boundary's impedance,
H/sq. This option should only be used along with the corresponding `"Rs"` and `"Cs"`
parameters, and not with any of the circuit parameters `"R"`, `"L"`, or `"C"`.

`"Cs" [0.0]` :  Surface capacitance used for computing this lumped port boundary's
impedance, F/sq. This option should only be used along with the corresponding `"Rs"` and
`"Ls"` parameters, and not with any of the circuit parameters `"R"`, `"L"`, or `"C"`.

`"Elements"[]["Attributes"] [None]` :  This option is for multielement lumped ports and
should not be combined with the `"Attributes"` field described above. Each element of a
multielement lumped port can be described by its own unique integer array of mesh boundary
attributes, which are specified here. The elements of a multielement port add in parallel.

`"Elements"[]["Direction"] [None]` :  This option is for multielement lumped ports and
should not be combined with the `"Direction"` field described above. Each element of a
multielement lumped port can be described by its own unique direction, which is specified
here. The elements of a multielement port add in parallel.

`"Elements"[]["CoordinateSystem"] ["Cartesian"]` :  This option is for multielement lumped
ports and should not be combined with the `"CoordinateSystem"` field described above. Each
element of a multielement lumped port can be described by its own unique direction, and
corresponding coordinate system.

## `boundaries["WavePort"]`

```json
"WavePort":
[
    {
        "Index": <int>,
        "Attributes": [<int array>],
        "Excitation": <bool>,
        "Active": <bool>,
        "Mode": <int>,
        "Offset": <float>
    },
    ...
]
```

with

`"Index" [None]` :  Index of this wave port boundary, used in postprocessing output files.

`"Attributes" [None]` :  Integer array of mesh boundary attributes for this wave port
boundary.

`"Excitation" [false]` :  Turns on or off port excitation for this wave port boundary for
driven simulation types.

`"Active" [true]` :  Turns on or off damping boundary condition for this wave port boundary
for driven simulation types.

`"Mode" [1]` :  Mode index (1-based) for the characteristic port mode of this wave port
boundary. Ranked in order of decreasing wave number.

`"Offset" [0.0]` :  Offset distance used for scattering parameter de-embedding for this wave
port boundary, specified in mesh length units.

## `boundaries["WavePortPEC"]`

```json
"WavePortPEC":
{
    "Attributes": [<int array>]
}
```

with

`"Attributes" [None]` :  Integer array of mesh boundary attributes to consider along with
those specified under
[`config["Boundaries"]["PEC"]["Attributes"]`](#boundaries%5B%22PEC%22%5D) as PEC when
performing wave port boundary mode analysis.

## `boundaries["SurfaceCurrent"]`

```json
"SurfaceCurrent":
[
    {
        "Index": <int>,
        "Attributes": [<int array>],
        "Direction": <string> or [<float array>],
        "CoordinateSystem": <string>,
        "Elements":
        [
            {
                "Attributes": [<int array>],
                "Direction": <string> or [<float array>],
                "CoordinateSystem": <string>,
            },
            ...
        ]
    },
    ...
]
```

with

`"Index" [None]` :  Index of this surface current boundary, used in postprocessing output
files.

`"Attributes" [None]` :  Integer array of mesh boundary attributes for this surface current
boundary. If this source is to be a multielement source which distributes the source
across more than a single lumped element, use the `"Elements"` array described below.

`"Direction" [None]` :  Defines the source current direction for this surface current
boundary. The available options are the same as under
[`config["Boundaries"]["LumpedPort"]["Direction"]`](#boundaries%5B%22LumpedPort%22%5D). If
this source is to be a multielement source which distributes the source across more than a
single lumped element, use the `"Elements"` array described below.

`"CoordinateSystem" ["Cartesian"]` :  Defines the coordinate system for the source current
direction for this surface current boundary. The available options are the same as under
[`config["Boundaries"]["LumpedPort"]["CoordinateSystem"]`](#boundaries%5B%22LumpedPort%22%5D).
If this source is to be a multielement source which distributes the source across more than
a single lumped element, use the `"Elements"` array described below.

`"Elements"[]["Attributes"] [None]` :  This option is for multielement surface current
boundaries should not be combined with the `"Attributes"` field described above. Each
element of a multielement current source can be described by its own unique integer array of
mesh boundary attributes, which are specified here. The elements of a multielement source
add in parallel to give the same total current as a single-element source.

`"Elements"[]["Direction"] [None]` :  This option is for multielement surface current
boundaries and should not be combined with the `"Direction"` field described above. Each
element of a multielement current source can be described by its own unique direction,
which is specified here. The elements of a multielement source add in parallel to give the
same total current as a single-element source.

`"Elements"[]["CoordinateSystem"] ["Cartesian"]` :  This option is for multielement surface
current boundaries and should not be combined with the `"CoordinateSystem"` field described
above. Each element of a multielement current source can be described by its own unique
direction, and corresponding coordinate system.

## `boundaries["Ground"]`

```json
"Ground":
{
    "Attributes": [<int array>]
}
```

with

`"Attributes" [None]` :  Integer array of mesh boundary attributes at which to apply the
ground boundary condition.

## `boundaries["ZeroCharge"]`

```json
"ZeroCharge":
{
    "Attributes": [<int array>]
}
```

with

`"Attributes" [None]` :  Integer array of mesh boundary attributes at which to apply the
zero-charge boundary condition.

## `boundaries["Terminal"]`

```json
"Terminal":
[
    {
        "Index": <int>
        "Attributes": [<int array>],
    },
    ...
]
```

with

`"Index" [None]` :  Index of this terminal boundary, used in postprocessing output files and
to index the computed capacitance matrix.

`"Attributes" [None]` :  Integer array of mesh boundary attributes for this terminal
boundary.

## `boundaries["Postprocessing"]["SurfaceFlux"]`

```json
"Postprocessing":
{
    "SurfaceFlux":
    [
        {
            "Index": <int>
            "Attributes": [<int array>],
            "Type": <string>,
            "TwoSided": <bool>,
            "Center": [<float array>]
        },
        ...
    ]
}
```

with

`"Index" [None]` :  Index of this surface flux postprocessing boundary, used in
postprocessing output files.

`"Attributes" [None]` :  Integer array of mesh boundary attributes for this surface flux
postprocessing boundary.

`"Type" [None]` :  Specifies the type of surface flux to calculate for this postprocessing
boundary. The available options are:

  - `"Electric"` :  Integrate the electric flux density over the boundary surface.
  - `"Magnetic"` :  Integrate the magnetic flux density over the boundary surface.
  - `"Power"` :  Integrate the energy flux density, given by the Poynting vector, over the
    boundary surface.

`"TwoSided" [false]` :  Specifies how to account for internal boundary surfaces with a
possible discontinuous field on either side. When set to `false`, the flux on either side of
an internal boundary surface is averaged. When `true`, it is summed with an opposite normal
direction.

`"Center" [None]` :  Floating point array of length equal to the model spatial dimension
specfiying the coordinates of a central point used to compute the outward flux. The true
surface normal is used in the calculation, and this point is only used to ensure the correct
orientation of the normal. Specified in mesh length units, and only relevant when
`"TwoSided"` is `false`. If not specified, the point will be computed as the centroid of the
axis-aligned bounding box for all elements making up the postprocessing boundary.

## `boundaries["Postprocessing"]["Dielectric"]`

```json
"Postprocessing":
{
    "Dielectric":
    [
        {
            "Index": <int>,
            "Attributes": [<int array>],
            "Type": <string>,
            "Thickness": <float>,
            "Permittivity": <float>,
            "LossTan": <float>,
            "Side": <string>
        },
        ...
    ]
}
```

with

`"Index" [None]` :  Index of this dielectric interface, used in postprocessing output files.

`"Attributes" [None]` :  Integer array of mesh boundary attributes for this dielectric
interface. If the interface consists of multiple elements with different `"Side"` values,
use the `"Elements"` array described below.

`"Type" [None]` :  Specifies the type of dielectric interface for this postprocessing
boundary. See also [this page](../reference.md#Bulk-and-interface-dielectric-loss).
Available options are:

  - `"Default"` :  Use the full electric field evaulated at the boundary to compute the
    energy participation ratio (EPR) of this dielectric interface and estimate loss.
  - `"MA"` :  Use the boundary conditions assuming a metal-air interface to compute the EPR
    of this dielectric interface.
  - `"MS"` :  Use the boundary conditions assuming a metal-substrate interface to compute
    the EPR of this dielectric interface.
  - `"SA"` :  Use the boundary conditions assuming a substrate-air interface to compute the
    EPR of this dielectric interface.

`"Thickness" [None]` :  Thickness of this dielectric interface, specified in mesh length
units.

`"Permittivity" [None]` :  Relative permittivity for this dielectric interface. This should
be the interface layer permittivity for the specific `"Type"` of interface specified.

`"LossTan" [0.0]` :  Loss tangent for this lossy dielectric interface.

`"Side" ["SmallerRefractiveIndex"]` :  Defines the postprocessing side when this dielectric
interface is an internal boundary surface (and thus the electric field on the boundary is in
general double-valued). The available options are:

  - `"SmallerRefractiveIndex"` :  Take the value from the side where the material index of
    refraction is smaller (speed of light is larger). Typically this selects the vacuum
    side. For anisotropic materials, the index of refraction associated with the principal
    direction with the smallest value is used.
  - `"LargerRefractiveIndex"` :  Take the value from the side where the material index of
    refraction is larger (speed of light is smaller). Typically this selects the non-vacuum
    side.
