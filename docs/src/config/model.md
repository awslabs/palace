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

`"Mesh" [None]` :  Input mesh file path, an absolute path is recommended. If the provided
mesh is nonconformal, it is assumed that it comes from a previous *Palace* solve using AMR,
and all mesh preprocessing checks and modifications (for example
[`model["Refinement"]["CrackInternalBoundaryElements"]`](#model%5B%22Refinement%22%5D)), are
skipped .

`"L0" [1.0e-6]` :  Unit, relative to m, for mesh vertex coordinates. For example, a value
of `1.0e-6` implies the mesh coordinates are in μm.

`"Lc" [0.0]` :  Characteristic length scale used for nondimensionalization, specified in
mesh length units. This keyword should typically not be specified by the user. A value less
than or equal to zero uses an internally calculated length scale based on the bounding box
of the computational domain. A value of 1.0 will disable nondimensionalization of lengths
in the model and all computations will take place in the same units as the mesh.

`"Refinement"` : Top-level object for configuring mesh refinement.

## `model["Refinement"]`

```json
"Refinement":
{
    "Tol": <float>,
    "MaxIts": <int>,
    "MaxSize": <int>,
    "Nonconformal": <bool>,
    "UpdateFraction": <float>,
    "UniformLevels": <int>,
    "Boxes":
    [
        {
            "Levels": <int>,
            "BoundingBoxMin": [<float array>],
            "BoundingBoxMax": [<float array>]
        },
        ...
    ],
    "Spheres":
    [
        {
            "Levels": <int>,
            "Center": [<float array>],
            "Radius": <float>
        },
        ...
    ]
}
```

with

`"Tol" [1.0e-2]` : Relative error convergence tolerance for adaptive mesh refinement (AMR).

`"MaxIts" [0]` : Maximum number of iterations of AMR to perform.

`"MaxSize" [0]` : The maximum allowable number of degrees of freedom for AMR. If an adapted
mesh exceeds this value no further adaptation will occur. A value less than 1 means that no
maximum size constraint will be imposed.

`"Nonconformal" [true]` : Chose whether the adaptation should use nonconformal refinement.
Nonconformal refinement is required for non-simplex meshes.

`"UpdateFraction" [0.7]` : Dörfler marking fraction used to specify which elements to
refine. This marking strategy will mark the smallest number of elements that make up
`"UpdateFraction"` of the total error in the mesh. A larger value will refine more elements
per iteration, at the cost of the final mesh being less efficient.

`"UniformLevels" [0]` :  Levels of uniform parallel mesh refinement to be performed on the
input mesh. If not performing AMR, these may be used as levels within a geometric multigrid
scheme. If performing AMR the most refined mesh is used as the initial mesh and the coarser
meshes cannot be used in a geometric multigrid scheme.

`"Boxes"` :  Array of box region refinement objects. All elements with a node inside the box
region will be marked for refinement.

`"Spheres"` :  Array of sphere region refinement objects. All elements with a node inside
the sphere region will be marked for refinement.

`"Levels" [0]` : Levels of parallel mesh refinement inside the specified refinement region.

`"BoundingBoxMin" [None]` : Floating point array of length equal to the model spatial
dimension specifying the minimum coordinates of the axis-aligned bounding box for this
refinement region. Specified in mesh length units.

`"BoundingBoxMax" [None]` : Floating point array of length equal to the model spatial
dimension specifying the maximum coordinates of the axis-aligned bounding box for this
refinement region. Specified in mesh length units.

`"Center" [None]` : Floating point array of length equal to the model spatial dimension
specifying the center coordinates of the sphere for this sphere refinement region.
Specified in mesh length units.

`"Radius" [None]` : The radius of the sphere for this sphere refinement region, specified in
mesh length units.

### Advanced model options

  - `"RemoveCurvature" [false]`
  - `"MakeSimplex" [false]`
  - `"MakeHexahedral" [false]`
  - `"ReorderElements" [false]`
  - `"CleanUnusedElements" [true]`
  - `"CrackInternalBoundaryElements" [true]`
  - `"RefineCrackElements" [true]`
  - `"CrackDisplacementFactor" [1.0e-12]`
  - `"AddInterfaceBoundaryElements" [true]`
  - `"ExportPrerefinedMesh" [false]`
  - `"ReorientTetMesh" [false]`
  - `"Partitioning" [""]`
  - `"MaxNCLevels" [1]`
  - `"MaximumImbalance" [1.1]`
  - `"SaveAdaptIterations" [true]`
  - `"SaveAdaptMesh" [false]`
  - `"SerialUniformLevels" [0]`
