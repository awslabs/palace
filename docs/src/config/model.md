```@raw html
<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
```

# `config["Model"]`

```json
"Model":
{
    "Mesh": <string>
    "L0": <float>,
    "Lc": <float>,
    "Refinement":
    {
        ...
    }
}
```

with

`"Mesh" [None]` :  Input mesh file path, an absolute path is recommended.

`"L0" [1.0e-6]` :  Mesh vertex coordinate length unit, m.

`"Lc" [0.0]` :  Characteristic length scale used for nondimensionalization, specified in
mesh length units. A value less than or equal to zero uses an internally calculated length
scale based on the bounding box of the computational domain.

`"Refinement"` : Top-level object for configuring mesh refinement.

## `model["Refinement"]`

```json
"Refinement":
{
    "UniformLevels": <int>,
    "Boxes":
    [
        {
            "Levels": <int>,
            "XLimits": [<float array>],
            "YLimits": [<float array>],
            "ZLimits": [<float array>]
        },
        ...
    ],
    "Spheres":
    [
        {
            "Levels": <int>,
            "Center": [<float array>],
            "Radius": float
        },
        ...
    ],
    "Adaptation":
    {
        "Tol": <float>,
        "MaxIts": <int>,
        "UpdateFraction": <float>,
        "DOFLimit": <int>,
        "SaveStep": <int>,
        "MaximumImbalance": <float>,
        "Nonconformal": <bool>,
        "UseCoarsening": <bool>,
        "WriteSerialMesh": <bool>
    }
}
```

with

`"UniformLevels" [0]` :  Levels of uniform parallel mesh refinement to be performed on the
input mesh.

`"Boxes"` :  Array of box region refinement objects. All elements with a node inside the box
region will be marked for refinement.

`"Spheres"` :  Array of sphere region refinement objects. All elements with a node inside
the sphere region will be marked for refinement.

`"Levels" [0]` : Levels of parallel mesh refinement inside the specified refinement region.

`"XLimits" [None]` : Floating point array of length 2 specifying the limits in the
``x``-direction of the axis-aligned bounding box for this box refinement region. Specified
in mesh length units.

`"YLimits" [None]` : Floating point array of length 2 specifying the limits in the
``y``-direction of the axis-aligned bounding box for this box refinement region. Specified
in mesh length units.

`"ZLimits" [None]` : Floating point array of length 2 specifying the limits in the
``z``-direction of the axis-aligned bounding box for this box refinement region. Specified
in mesh length units.

`"Center" [None]` : Floating point array of length equal to the model spatial dimension
specfiying the center coordinates of the sphere for this sphere refinement region.
Specified in mesh length units.

`"Radius" [None]` : The radius of the sphere for this sphere refinement region, specified in
mesh length units.

`"Tol" [1e-3]` : Relative error convergence tolerance for mesh adaptation.

`"MaxIts" [0]` : Maximum number of iterations to perform of adaptive mesh refinement.

`"UpdateFraction" [0.4]` : Dörfler marking fraction used to specify which elements to
refine. This marking strategy will mark the smallest number of elements that make up
`"UpdateFraction"` of the total error in the mesh.

`"DOFLimit" [0]` : The maximum allowable number of DOF within an AMR loop, if an adapted
mesh exceeds this value no further adaptation will occur.

`"SaveStep" [1]` : Specify which adaptive iterations to save off within the post processing
output folder. `"SaveStep"` of 1 specifies to save all iterations, while for example
`"SaveStep"` of 3 would save every third iteration.

`"MaximumImbalance" [1.0]` : The ratio between maximum number of elements on a processor and
minimum number of elements on a processor before a rebalance is performed. A value of `2.0`
would result in rebalancing occurring only if one processor had more than double the number
of elements on another.

`"Nonconformal" [false]` : Whether the adaptation should use nonconformal refinement.
`"Nonconformal"` is necessary to enable `"UseCoarsening"`.

`"UseCoarsening" [false]` : Whether or not to perform coarsening if the total number of DOF
exceeds the `"DOFLimit"`. Coarsening may be useful to improve the efficiency of the mesh if
a large value of `"UpdateFraction"` is used. Requires `"Nonconformal"` mesh refinement to be
enabled.

`"WriteSerialMesh" [true]` : Whether to write a serialized mesh to file after adaptation.

### Advanced model options

  - `"MaxNCLevels" [0]`
  - `"Partition" [""]`
  - `"ReorientTetMesh" [false]`
  - `"RemoveCurvature" [false]`
  - `"WritePostBalanceMesh" [false]`
  - `"WritePreBalanceMesh" [false]`
