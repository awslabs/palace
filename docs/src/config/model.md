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
    ]
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

### Advanced model options

  - `"Partition" [""]`
  - `"ReorientTetMesh" [false]`
  - `"RemoveCurvature" [false]`
