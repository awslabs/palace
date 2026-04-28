```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# [Configuration File Reference](@id config-reference)

```@raw html
<dl class="config-toc">
  <dt><a href="#config-problem">Problem</a></dt>
  <dd><a href="#config-problem-outputformats">Output Formats</a></dd>
  <dt><a href="#config-model">Model</a></dt>
  <dd><a href="#config-model-refinement">Mesh Refinement</a></dd>
  <dt><a href="#config-domains">Domains</a></dt>
  <dd><a href="#config-domains-materials">Materials</a> · <a href="#config-domains-currentdipole">CurrentDipole</a> · <a href="#config-domains-postprocessing">Domain Postprocessing</a></dd>
  <dt><a href="#config-boundaries">Boundaries</a></dt>
  <dd><a href="#config-boundaries-pec">PEC Boundary</a> · <a href="#config-boundaries-pmc">PMC Boundary</a> · <a href="#config-boundaries-impedance">Impedance</a> · <a href="#config-boundaries-absorbing">Absorbing Boundary</a> · <a href="#config-boundaries-conductivity">Conductivity</a> · <a href="#config-boundaries-lumpedport">Lumped Port</a> · <a href="#config-boundaries-waveport">Wave Port</a> · <a href="#config-boundaries-waveportpec">Wave Port PEC</a> · <a href="#config-boundaries-surfacecurrent">Surface Current</a> · <a href="#config-boundaries-ground">Ground Boundary</a> · <a href="#config-boundaries-zerocharge">Zero Charge Boundary</a> · <a href="#config-boundaries-terminal">Terminal</a> · <a href="#config-boundaries-periodic">Periodic Boundary</a> · <a href="#config-boundaries-postprocessing">Boundary Postprocessing</a></dd>
  <dt><a href="#config-solver">Solver</a></dt>
  <dd><a href="#config-solver-eigenmode">Eigenmode Solver</a> · <a href="#config-solver-driven">Driven Solver</a> · <a href="#config-solver-transient">Transient Solver</a> · <a href="#config-solver-electrostatic">Electrostatic Solver</a> · <a href="#config-solver-magnetostatic">Magnetostatic Solver</a> · <a href="#config-solver-linear">Linear Solver</a></dd>
</dl>
```

## [Problem](@id config-problem)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-problem">Problem</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-required">required</span></span>
```

Top-level configuration for the simulation type and output.

```@raw html
<dl class="palace-config">
  <dt id="config-problem-type"><a href="#config-problem-type"><code>"Type"</code></a> <span class="config-type">string</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Controls the simulation type.</p>
    <dl class="config-enum">
      <dt><code>"Eigenmode"</code></dt>
      <dd>Perform an undamped or damped eigenfrequency analysis.</dd>
      <dt><code>"Driven"</code></dt>
      <dd>Perform a frequency-domain driven simulation.</dd>
      <dt><code>"Transient"</code></dt>
      <dd>Perform a time-domain excitation response simulation.</dd>
      <dt><code>"Electrostatic"</code></dt>
      <dd>Perform an electrostatic analysis to compute the capacitance matrix for a set of voltage terminals.</dd>
      <dt><code>"Magnetostatic"</code></dt>
      <dd>Perform a magnetostatic analysis to compute the inductance matrix for a set of current sources.</dd>
    </dl>
  </dd>
  <dt id="config-problem-verbose"><a href="#config-problem-verbose"><code>"Verbose"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>1</code></span> <span class="config-constraint"><code>≥ 0</code></span></dt>
  <dd>
    <p>Controls the level of log file printing.</p>
  </dd>
  <dt id="config-problem-output"><a href="#config-problem-output"><code>"Output"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;&quot;</code></span></dt>
  <dd>
    <p>Directory path for saving postprocessing outputs.</p>
  </dd>
  <dt><code>"OutputFormats"</code> <span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-problem-outputformats">See full reference ↓</a></p>
  </dd>
</dl>
```

### [Output Formats](@id config-problem-outputformats)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-problem">Problem</a>/<a href="#config-problem-outputformats">OutputFormats</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></span>
```

Configures the field output formats.

```@raw html
<dl class="palace-config">
  <dt id="config-problem-outputformats-paraview"><a href="#config-problem-outputformats-paraview"><code>"Paraview"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>true</code></span></dt>
  <dd>
    <p>Set to <code>true</code> to output fields in <a href="https://www.paraview.org/">ParaView</a> format.</p>
  </dd>
  <dt id="config-problem-outputformats-gridfunction"><a href="#config-problem-outputformats-gridfunction"><code>"GridFunction"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span></dt>
  <dd>
    <p>Set to <code>true</code> to output fields in MFEM grid function format for visualization with <a href="https://glvis.org/">GLVis</a>.</p>
  </dd>
</dl>
```

## [Model](@id config-model)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-model">Model</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-required">required</span></span>
```

Mesh and model configuration.

```@raw html
<dl class="palace-config">
  <dt id="config-model-mesh"><a href="#config-model-mesh"><code>"Mesh"</code></a> <span class="config-type">string</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Input mesh file path. An absolute path is recommended. If the provided mesh is nonconformal, it is assumed to come from a previous <em>Palace</em> AMR solve, and all mesh preprocessing checks and modifications (for example <a href="#config-model-crackinternalboundaryelements">CrackInternalBoundaryElements</a>) are skipped.</p>
  </dd>
  <dt id="config-model-l0"><a href="#config-model-l0"><code>"L0"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>1.0e-6</code></span> <span class="config-constraint"><code>&gt; 0.0</code></span></dt>
  <dd>
    <p>Unit, relative to meters, for mesh vertex coordinates. For example, a value of <code>1.0e-6</code> means the mesh coordinates are in μm.</p>
  </dd>
  <dt id="config-model-lc"><a href="#config-model-lc"><code>"Lc"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.0</code></span></dt>
  <dd>
    <p>Characteristic length scale used for nondimensionalization, specified in mesh length units. This keyword should typically not be specified by the user. A value less than or equal to zero uses an internally calculated length scale based on the bounding box of the computational domain. A value of <code>1.0</code> will disable nondimensionalization, so that all computations will take place in the same units as the mesh.</p>
  </dd>
  <dt id="config-model-removecurvature"><a href="#config-model-removecurvature"><code>"RemoveCurvature"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Project high-order nodes to the mesh surface, removing all curvature before the simulation.</p>
  </dd>
  <dt id="config-model-makesimplex"><a href="#config-model-makesimplex"><code>"MakeSimplex"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Convert all mesh elements to simplices (tetrahedra/triangles).</p>
  </dd>
  <dt id="config-model-makehexahedral"><a href="#config-model-makehexahedral"><code>"MakeHexahedral"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Convert all mesh elements to hexahedra.</p>
  </dd>
  <dt id="config-model-reorderelements"><a href="#config-model-reorderelements"><code>"ReorderElements"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Reorder mesh elements to improve cache efficiency.</p>
  </dd>
  <dt id="config-model-cleanunusedelements"><a href="#config-model-cleanunusedelements"><code>"CleanUnusedElements"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>true</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Remove elements not connected to any domain material.</p>
  </dd>
  <dt id="config-model-crackinternalboundaryelements"><a href="#config-model-crackinternalboundaryelements"><code>"CrackInternalBoundaryElements"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>true</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Duplicate nodes along internal boundary elements to create a crack.</p>
  </dd>
  <dt id="config-model-refinecrackelements"><a href="#config-model-refinecrackelements"><code>"RefineCrackElements"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>true</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Refine elements adjacent to cracked internal boundaries.</p>
  </dd>
  <dt id="config-model-crackdisplacementfactor"><a href="#config-model-crackdisplacementfactor"><code>"CrackDisplacementFactor"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>1.0e-12</code></span> <span class="config-constraint"><code>≥ 0.0</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Displacement factor applied to cracked nodes as a fraction of the local element size.</p>
  </dd>
  <dt id="config-model-addinterfaceboundaryelements"><a href="#config-model-addinterfaceboundaryelements"><code>"AddInterfaceBoundaryElements"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>true</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Add boundary elements at interfaces between domains that lack them.</p>
  </dd>
  <dt id="config-model-exportprerefinedmesh"><a href="#config-model-exportprerefinedmesh"><code>"ExportPrerefinedMesh"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Export the mesh after preprocessing but before AMR refinement.</p>
  </dd>
  <dt id="config-model-reorienttetmesh"><a href="#config-model-reorienttetmesh"><code>"ReorientTetMesh"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Reorient tetrahedral elements to ensure positive Jacobians.</p>
  </dd>
  <dt id="config-model-partitioning"><a href="#config-model-partitioning"><code>"Partitioning"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;&quot;</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Path to a mesh partitioning file. If empty, partitioning is computed automatically.</p>
  </dd>
  <dt><code>"Refinement"</code> <span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-model-refinement">See full reference ↓</a></p>
  </dd>
</dl>
```

### [Mesh Refinement](@id config-model-refinement)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-model">Model</a>/<a href="#config-model-refinement">Refinement</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></span>
```

Configuration for adaptive and uniform mesh refinement.

```@raw html
<dl class="palace-config">
  <dt id="config-model-refinement-tol"><a href="#config-model-refinement-tol"><code>"Tol"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.01</code></span> <span class="config-constraint"><code>&gt; 0.0</code></span></dt>
  <dd>
    <p>Stop adaptive mesh refinement (AMR) when the norm of the estimated error falls below this value. The error is reported in <code>error-indicators.csv</code>.</p>
  </dd>
  <dt id="config-model-refinement-maxits"><a href="#config-model-refinement-maxits"><code>"MaxIts"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>0</code></span> <span class="config-constraint"><code>≥ 0</code></span></dt>
  <dd>
    <p>Maximum number of AMR iterations to perform.</p>
  </dd>
  <dt id="config-model-refinement-maxsize"><a href="#config-model-refinement-maxsize"><code>"MaxSize"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0</code></span> <span class="config-constraint"><code>≥ 0</code></span></dt>
  <dd>
    <p>The maximum allowable number of degrees of freedom for AMR. If an adapted mesh exceeds this value no further adaptation will occur. A value less than 1 means that no maximum size constraint will be imposed.</p>
  </dd>
  <dt id="config-model-refinement-updatefraction"><a href="#config-model-refinement-updatefraction"><code>"UpdateFraction"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.7</code></span> <span class="config-constraint"><code>&gt; 0.0, &lt; 1.0</code></span></dt>
  <dd>
    <p>Dörfler marking fraction used to specify which elements to refine. This marking strategy will mark the smallest number of elements that make up &quot;UpdateFraction&quot; of the total error in the mesh. A larger value will refine more elements per iteration, at the cost of the final mesh being less efficient.</p>
  </dd>
  <dt id="config-model-refinement-nonconformal"><a href="#config-model-refinement-nonconformal"><code>"Nonconformal"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>true</code></span></dt>
  <dd>
    <p>Use nonconformal refinement in adaptation. Required for non-simplex meshes.</p>
  </dd>
  <dt id="config-model-refinement-maxnclevels"><a href="#config-model-refinement-maxnclevels"><code>"MaxNCLevels"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>1</code></span> <span class="config-constraint"><code>≥ 0</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Maximum number of nonconformal refinement levels. <code>0</code> means no limit.</p>
  </dd>
  <dt id="config-model-refinement-maximumimbalance"><a href="#config-model-refinement-maximumimbalance"><code>"MaximumImbalance"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>1.1</code></span> <span class="config-constraint"><code>≥ 1.0</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Maximum ratio of elements between the most- and least-loaded MPI ranks before repartitioning.</p>
  </dd>
  <dt id="config-model-refinement-saveadaptiterations"><a href="#config-model-refinement-saveadaptiterations"><code>"SaveAdaptIterations"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>true</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Save postprocessing results from each AMR iteration in a subdirectory <code>iterationX</code>.</p>
  </dd>
  <dt id="config-model-refinement-saveadaptmesh"><a href="#config-model-refinement-saveadaptmesh"><code>"SaveAdaptMesh"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Save the final adapted mesh to disk.</p>
  </dd>
  <dt id="config-model-refinement-uniformlevels"><a href="#config-model-refinement-uniformlevels"><code>"UniformLevels"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>0</code></span> <span class="config-constraint"><code>≥ 0</code></span></dt>
  <dd>
    <p>Levels of uniform parallel mesh refinement to be performed on the input mesh. If not performing AMR, these may be used as levels within a geometric multigrid scheme. If performing AMR the most refined mesh is used as the initial mesh and the coarser meshes cannot be used in a geometric multigrid scheme.</p>
  </dd>
  <dt id="config-model-refinement-serialuniformlevels"><a href="#config-model-refinement-serialuniformlevels"><code>"SerialUniformLevels"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>0</code></span> <span class="config-constraint"><code>≥ 0</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Levels of uniform serial mesh refinement applied before parallel distribution.</p>
  </dd>
  <dt><code>"Boxes"</code> <span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-model-refinement-boxes">See full reference ↓</a></p>
  </dd>
  <dt><code>"Spheres"</code> <span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-model-refinement-spheres">See full reference ↓</a></p>
  </dd>
</dl>
```

#### [Boxes](@id config-model-refinement-boxes)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-model">Model</a>/<a href="#config-model-refinement">Refinement</a>/<a href="#config-model-refinement-boxes">Boxes</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></span>
```

Array of axis-aligned box refinement regions. All elements with a node inside the box are marked for refinement.

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-model">Model</a>/<a href="#config-model-refinement">Refinement</a>/<a href="#config-model-refinement-boxes">Boxes</a>/0</code></p>
```

```@raw html
<dl class="palace-config">
  <dt id="config-model-refinement-boxes-levels"><a href="#config-model-refinement-boxes-levels"><code>"Levels"</code></a> <span class="config-type">integer</span> <span class="config-required">required</span> <span class="config-constraint"><code>≥ 0</code></span></dt>
  <dd>
    <p>Levels of parallel mesh refinement inside this box region.</p>
  </dd>
  <dt id="config-model-refinement-boxes-boundingboxmin"><a href="#config-model-refinement-boxes-boundingboxmin"><code>"BoundingBoxMin"</code></a> <span class="config-type">[number × 3]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Minimum coordinates <code>[x, y, z]</code> of the axis-aligned bounding box for this refinement region, in mesh length units.</p>
  </dd>
  <dt id="config-model-refinement-boxes-boundingboxmax"><a href="#config-model-refinement-boxes-boundingboxmax"><code>"BoundingBoxMax"</code></a> <span class="config-type">[number × 3]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Maximum coordinates <code>[x, y, z]</code> of the axis-aligned bounding box for this refinement region, in mesh length units.</p>
  </dd>
</dl>
```

#### [Spheres](@id config-model-refinement-spheres)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-model">Model</a>/<a href="#config-model-refinement">Refinement</a>/<a href="#config-model-refinement-spheres">Spheres</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></span>
```

Array of sphere refinement regions. All elements with a node inside the sphere are marked for refinement.

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-model">Model</a>/<a href="#config-model-refinement">Refinement</a>/<a href="#config-model-refinement-spheres">Spheres</a>/0</code></p>
```

```@raw html
<dl class="palace-config">
  <dt id="config-model-refinement-spheres-levels"><a href="#config-model-refinement-spheres-levels"><code>"Levels"</code></a> <span class="config-type">integer</span> <span class="config-required">required</span> <span class="config-constraint"><code>≥ 0</code></span></dt>
  <dd>
    <p>Levels of parallel mesh refinement inside this sphere region.</p>
  </dd>
  <dt id="config-model-refinement-spheres-radius"><a href="#config-model-refinement-spheres-radius"><code>"Radius"</code></a> <span class="config-type">number</span> <span class="config-required">required</span> <span class="config-constraint"><code>&gt; 0.0</code></span></dt>
  <dd>
    <p>Radius of the sphere, in mesh length units.</p>
  </dd>
  <dt id="config-model-refinement-spheres-center"><a href="#config-model-refinement-spheres-center"><code>"Center"</code></a> <span class="config-type">[number × 3]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Center coordinates <code>[x, y, z]</code> of the sphere, in mesh length units.</p>
  </dd>
</dl>
```

## [Domains](@id config-domains)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-domains">Domains</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-required">required</span></span>
```

Material and domain configuration.

```@raw html
<dl class="palace-config">
  <dt><code>"Materials"</code> <span class="config-type">[object, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p><a href="#config-domains-materials">See full reference ↓</a></p>
  </dd>
  <dt><code>"CurrentDipole"</code> <span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-domains-currentdipole">See full reference ↓</a></p>
  </dd>
  <dt><code>"Postprocessing"</code> <span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-domains-postprocessing">See full reference ↓</a></p>
  </dd>
</dl>
```

### [Materials](@id config-domains-materials)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-domains">Domains</a>/<a href="#config-domains-materials">Materials</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">[object, ...]</span> <span class="config-required">required</span></span>
```

Array of material property objects.

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-domains">Domains</a>/<a href="#config-domains-materials">Materials</a>/0</code></p>
```

```@raw html
<dl class="palace-config">
  <dt id="config-domains-materials-attributes"><a href="#config-domains-materials-attributes"><code>"Attributes"</code></a> <span class="config-type">[integer, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Integer array of mesh domain this object applies to.</p>
  </dd>
  <dt id="config-domains-materials-permeability"><a href="#config-domains-materials-permeability"><code>"Permeability"</code></a> <span class="config-type">number or [number × 3]</span> <span class="config-default">default: <code>1.0</code></span></dt>
  <dd>
    <p>Relative permeability for this material. Scalar or vector of 3 coefficients corresponding to each of <code>&quot;MaterialAxes&quot;</code>.</p>
  </dd>
  <dt id="config-domains-materials-permittivity"><a href="#config-domains-materials-permittivity"><code>"Permittivity"</code></a> <span class="config-type">number or [number × 3]</span> <span class="config-default">default: <code>1.0</code></span></dt>
  <dd>
    <p>Relative permittivity for this material. Scalar or vector of 3 coefficients corresponding to each of <code>&quot;MaterialAxes&quot;</code>.</p>
  </dd>
  <dt id="config-domains-materials-losstan"><a href="#config-domains-materials-losstan"><code>"LossTan"</code></a> <span class="config-type">number or [number × 3]</span> <span class="config-default">default: <code>0.0</code></span></dt>
  <dd>
    <p>Loss tangent for this material. Scalar or vector of 3 coefficients corresponding to each of <code>&quot;MaterialAxes&quot;</code>.</p>
  </dd>
  <dt id="config-domains-materials-conductivity"><a href="#config-domains-materials-conductivity"><code>"Conductivity"</code></a> <span class="config-type">number or [number × 3]</span> <span class="config-default">default: <code>0.0</code></span></dt>
  <dd>
    <p>Electrical conductivity for this material, S/m. Activates the Ohmic loss model in this domain. Scalar or vector of 3 coefficients corresponding to each of <code>&quot;MaterialAxes&quot;</code>.</p>
  </dd>
  <dt id="config-domains-materials-londondepth"><a href="#config-domains-materials-londondepth"><code>"LondonDepth"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.0</code></span></dt>
  <dd>
    <p>London penetration depth for this material, specified in mesh length units. Activates the London equations-based model relating superconducting current and electromagnetic fields in this domain.</p>
  </dd>
  <dt id="config-domains-materials-materialaxes"><a href="#config-domains-materials-materialaxes"><code>"MaterialAxes"</code></a> <span class="config-type">[[number × 3] × 3]</span> <span class="config-default">default: <code>[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]</code></span></dt>
  <dd>
    <p>Axes directions for specification of anisotropic material properties. Required to be unit length and orthogonal.</p>
  </dd>
</dl>
```

### [CurrentDipole](@id config-domains-currentdipole)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-domains">Domains</a>/<a href="#config-domains-currentdipole">CurrentDipole</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></span>
```

Array of current dipole source excitations.

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-domains">Domains</a>/<a href="#config-domains-currentdipole">CurrentDipole</a>/0</code></p>
```

```@raw html
<dl class="palace-config">
  <dt id="config-domains-currentdipole-index"><a href="#config-domains-currentdipole-index"><code>"Index"</code></a> <span class="config-type">integer</span> <span class="config-required">required</span> <span class="config-constraint"><code>&gt; 0</code></span></dt>
  <dd>
    <p>Index of this current dipole source, used in postprocessing output files.</p>
  </dd>
  <dt id="config-domains-currentdipole-moment"><a href="#config-domains-currentdipole-moment"><code>"Moment"</code></a> <span class="config-type">number</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Current dipole moment magnitude, A·m.</p>
  </dd>
  <dt id="config-domains-currentdipole-center"><a href="#config-domains-currentdipole-center"><code>"Center"</code></a> <span class="config-type">[number × 3]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Coordinates of the dipole center position <code>[x, y, z]</code>, in mesh length units.</p>
  </dd>
  <dt id="config-domains-currentdipole-direction"><a href="#config-domains-currentdipole-direction"><code>"Direction"</code></a> <span class="config-type">string or [number × 3]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Direction of the Dirac current source specifying the dipole. Axis-aligned directions can be specified using keywords: <code>&quot;+X&quot;</code>, <code>&quot;-X&quot;</code>, <code>&quot;+Y&quot;</code>, <code>&quot;-Y&quot;</code>, <code>&quot;+Z&quot;</code>, <code>&quot;-Z&quot;</code>. The direction can alternatively be specified as a normalized array of three values, for example <code>[0.0, 1.0, 0.0]</code>.</p>
  </dd>
</dl>
```

### [Domain Postprocessing](@id config-domains-postprocessing)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-domains">Domains</a>/<a href="#config-domains-postprocessing">Postprocessing</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></span>
```

Configuration for domain postprocessing.

```@raw html
<dl class="palace-config">
  <dt><code>"Energy"</code> <span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-domains-postprocessing-energy">See full reference ↓</a></p>
  </dd>
  <dt><code>"Probe"</code> <span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-domains-postprocessing-probe">See full reference ↓</a></p>
  </dd>
</dl>
```

#### [Energy](@id config-domains-postprocessing-energy)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-domains">Domains</a>/<a href="#config-domains-postprocessing">Postprocessing</a>/<a href="#config-domains-postprocessing-energy">Energy</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></span>
```

Array of objects for postprocessing domain energies. Postprocesses the electric and magnetic field energy inside a given domain. Results are written to `domain-E.csv` in the output directory.

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-domains">Domains</a>/<a href="#config-domains-postprocessing">Postprocessing</a>/<a href="#config-domains-postprocessing-energy">Energy</a>/0</code></p>
```

```@raw html
<dl class="palace-config">
  <dt id="config-domains-postprocessing-energy-index"><a href="#config-domains-postprocessing-energy-index"><code>"Index"</code></a> <span class="config-type">integer</span> <span class="config-required">required</span> <span class="config-constraint"><code>&gt; 0</code></span></dt>
  <dd>
    <p>Index of this energy postprocessing domain, used in output files.</p>
  </dd>
  <dt id="config-domains-postprocessing-energy-attributes"><a href="#config-domains-postprocessing-energy-attributes"><code>"Attributes"</code></a> <span class="config-type">[integer, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Integer array of mesh domain this object applies to.</p>
  </dd>
</dl>
```

#### [Probe](@id config-domains-postprocessing-probe)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-domains">Domains</a>/<a href="#config-domains-postprocessing">Postprocessing</a>/<a href="#config-domains-postprocessing-probe">Probe</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></span>
```

Array of probe points for evaluating field values at specified locations in space. The electric field **E** and magnetic flux density **B** are probed and written to `probe-E.csv` and `probe-B.csv` in the output directory.

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-domains">Domains</a>/<a href="#config-domains-postprocessing">Postprocessing</a>/<a href="#config-domains-postprocessing-probe">Probe</a>/0</code></p>
```

```@raw html
<dl class="palace-config">
  <dt id="config-domains-postprocessing-probe-index"><a href="#config-domains-postprocessing-probe-index"><code>"Index"</code></a> <span class="config-type">integer</span> <span class="config-required">required</span> <span class="config-constraint"><code>&gt; 0</code></span></dt>
  <dd>
    <p>Index of this probe, used in postprocessing output files.</p>
  </dd>
  <dt id="config-domains-postprocessing-probe-center"><a href="#config-domains-postprocessing-probe-center"><code>"Center"</code></a> <span class="config-type">[number × 3]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Coordinates of this probe <code>[x, y, z]</code>, in mesh length units.</p>
  </dd>
</dl>
```

## [Boundaries](@id config-boundaries)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-required">required</span></span>
```

Boundary condition configuration.

```@raw html
<dl class="palace-config">
  <dt><code>"PEC"</code> <span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-boundaries-pec">See full reference ↓</a></p>
  </dd>
  <dt><code>"PMC"</code> <span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-boundaries-pmc">See full reference ↓</a></p>
  </dd>
  <dt><code>"Impedance"</code> <span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-boundaries-impedance">See full reference ↓</a></p>
  </dd>
  <dt><code>"Absorbing"</code> <span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-boundaries-absorbing">See full reference ↓</a></p>
  </dd>
  <dt><code>"Conductivity"</code> <span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-boundaries-conductivity">See full reference ↓</a></p>
  </dd>
  <dt><code>"LumpedPort"</code> <span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-boundaries-lumpedport">See full reference ↓</a></p>
  </dd>
  <dt><code>"WavePort"</code> <span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-boundaries-waveport">See full reference ↓</a></p>
  </dd>
  <dt><code>"WavePortPEC"</code> <span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-boundaries-waveportpec">See full reference ↓</a></p>
  </dd>
  <dt><code>"SurfaceCurrent"</code> <span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-boundaries-surfacecurrent">See full reference ↓</a></p>
  </dd>
  <dt><code>"Ground"</code> <span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-boundaries-ground">See full reference ↓</a></p>
  </dd>
  <dt><code>"ZeroCharge"</code> <span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-boundaries-zerocharge">See full reference ↓</a></p>
  </dd>
  <dt><code>"Terminal"</code> <span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-boundaries-terminal">See full reference ↓</a></p>
  </dd>
  <dt><code>"Periodic"</code> <span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-boundaries-periodic">See full reference ↓</a></p>
  </dd>
  <dt><code>"Postprocessing"</code> <span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-boundaries-postprocessing">See full reference ↓</a></p>
  </dd>
</dl>
```

### [PEC Boundary](@id config-boundaries-pec)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-pec">PEC</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></span>
```

Perfect electric conductor (PEC) boundary condition: enforces zero tangential electric field. This is a homogeneous Dirichlet condition for frequency/time domain and magnetostatic formulations.

```@raw html
<dl class="palace-config">
  <dt id="config-boundaries-pec-attributes"><a href="#config-boundaries-pec-attributes"><code>"Attributes"</code></a> <span class="config-type">[integer, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Integer array of mesh boundary attributes this object applies to.</p>
  </dd>
</dl>
```

### [PMC Boundary](@id config-boundaries-pmc)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-pmc">PMC</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></span>
```

Perfect magnetic conductor (PMC) boundary condition: enforces zero tangential magnetic field. This is the natural (homogeneous Neumann) boundary condition; it also imposes symmetry of the electric field across the surface.

```@raw html
<dl class="palace-config">
  <dt id="config-boundaries-pmc-attributes"><a href="#config-boundaries-pmc-attributes"><code>"Attributes"</code></a> <span class="config-type">[integer, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Integer array of mesh boundary attributes this object applies to.</p>
  </dd>
</dl>
```

### [Impedance](@id config-boundaries-impedance)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-impedance">Impedance</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></span>
```

Array of surface impedance boundary conditions. The surface impedance relates the tangential electric and magnetic fields using the parallel combination of the specified resistance, inductance, and capacitance per square.

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-impedance">Impedance</a>/0</code></p>
```

```@raw html
<dl class="palace-config">
  <dt id="config-boundaries-impedance-attributes"><a href="#config-boundaries-impedance-attributes"><code>"Attributes"</code></a> <span class="config-type">[integer, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Integer array of mesh boundary attributes this object applies to.</p>
  </dd>
  <dt id="config-boundaries-impedance-rs"><a href="#config-boundaries-impedance-rs"><code>"Rs"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.0</code></span></dt>
  <dd>
    <p>Surface resistance for this impedance boundary, Ω/sq.</p>
  </dd>
  <dt id="config-boundaries-impedance-ls"><a href="#config-boundaries-impedance-ls"><code>"Ls"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.0</code></span></dt>
  <dd>
    <p>Surface inductance for this impedance boundary, H/sq.</p>
  </dd>
  <dt id="config-boundaries-impedance-cs"><a href="#config-boundaries-impedance-cs"><code>"Cs"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.0</code></span></dt>
  <dd>
    <p>Surface capacitance for this impedance boundary, F/sq.</p>
  </dd>
</dl>
```

### [Absorbing Boundary](@id config-boundaries-absorbing)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-absorbing">Absorbing</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></span>
```

Farfield absorbing (scattering) boundary conditions. These are artificial boundary conditions applied at farfield boundaries to minimize reflections.

```@raw html
<dl class="palace-config">
  <dt id="config-boundaries-absorbing-attributes"><a href="#config-boundaries-absorbing-attributes"><code>"Attributes"</code></a> <span class="config-type">[integer, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Integer array of mesh boundary attributes this object applies to.</p>
  </dd>
  <dt id="config-boundaries-absorbing-order"><a href="#config-boundaries-absorbing-order"><code>"Order"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>1</code></span> <span class="config-constraint"><code>≥ 1, ≤ 2</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Specify a first- or second-order approximation for the absorbing boundary condition. Second-order is only available for frequency domain driven simulations.</p>
  </dd>
</dl>
```

### [Conductivity](@id config-boundaries-conductivity)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-conductivity">Conductivity</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></span>
```

Array of finite conductivity surface impedance boundaries. Models the effect of a boundary with non-infinite conductivity for conductors with thickness much larger than the skin depth. Only available for frequency domain driven and eigenmode simulations.

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-conductivity">Conductivity</a>/0</code></p>
```

```@raw html
<dl class="palace-config">
  <dt id="config-boundaries-conductivity-attributes"><a href="#config-boundaries-conductivity-attributes"><code>"Attributes"</code></a> <span class="config-type">[integer, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Integer array of mesh boundary attributes this object applies to.</p>
  </dd>
  <dt id="config-boundaries-conductivity-conductivity"><a href="#config-boundaries-conductivity-conductivity"><code>"Conductivity"</code></a> <span class="config-type">number</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Electrical conductivity for this boundary, S/m.</p>
  </dd>
  <dt id="config-boundaries-conductivity-permeability"><a href="#config-boundaries-conductivity-permeability"><code>"Permeability"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>1.0</code></span></dt>
  <dd>
    <p>Relative permeability for this boundary.</p>
  </dd>
  <dt id="config-boundaries-conductivity-thickness"><a href="#config-boundaries-conductivity-thickness"><code>"Thickness"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.0</code></span></dt>
  <dd>
    <p>Optional conductor thickness in mesh length units. Activates a finite-thickness boundary condition for metal.</p>
  </dd>
  <dt id="config-boundaries-conductivity-external"><a href="#config-boundaries-conductivity-external"><code>"External"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span></dt>
  <dd>
    <p>Whether this boundary is on the exterior of the computational domain. Relevant for the thickness correction.</p>
  </dd>
</dl>
```

### [Lumped Port](@id config-boundaries-lumpedport)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-lumpedport">LumpedPort</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></span>
```

Array of lumped port boundary conditions. Lumped ports can be specified on boundaries internal to the computational domain.

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-lumpedport">LumpedPort</a>/0</code></p>
```

```@raw html
<dl class="palace-config">
  <dt id="config-boundaries-lumpedport-index"><a href="#config-boundaries-lumpedport-index"><code>"Index"</code></a> <span class="config-type">integer</span> <span class="config-required">required</span> <span class="config-constraint"><code>&gt; 0</code></span></dt>
  <dd>
    <p>Index of this lumped port, used in postprocessing output files. Must be unique across all port and source types.</p>
  </dd>
  <dt id="config-boundaries-lumpedport-attributes"><a href="#config-boundaries-lumpedport-attributes"><code>"Attributes"</code></a> <span class="config-type">[integer, ...]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p>Integer array of mesh boundary attributes for this lumped port boundary. If this port is to be a multielement lumped port with more than a single lumped element, use the &quot;Elements&quot; array described below.</p>
  </dd>
  <dt id="config-boundaries-lumpedport-direction"><a href="#config-boundaries-lumpedport-direction"><code>"Direction"</code></a> <span class="config-type">string or [number × 3]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p>Excitation direction. Axis-aligned Cartesian directions can be specified using keywords: <code>&quot;+X&quot;</code>, <code>&quot;-X&quot;</code>, <code>&quot;+Y&quot;</code>, <code>&quot;-Y&quot;</code>, <code>&quot;+Z&quot;</code>, <code>&quot;-Z&quot;</code>. Coaxial directions use <code>&quot;+R&quot;</code>, <code>&quot;-R&quot;</code>. Alternatively, specify a normalized 3-element array, e.g. <code>[0.0, 1.0, 0.0]</code>. The coordinate system is determined by <code>&quot;CoordinateSystem&quot;</code>.</p>
  </dd>
  <dt id="config-boundaries-lumpedport-coordinatesystem"><a href="#config-boundaries-lumpedport-coordinatesystem"><code>"CoordinateSystem"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;Cartesian&quot;</code></span></dt>
  <dd>
    <p>Coordinate system used to express the <code>&quot;Direction&quot;</code> vector. If a keyword argument is used for <code>&quot;Direction&quot;</code> this value is ignored.</p>
    <dl class="config-enum">
      <dt><code>"Cartesian"</code></dt>
      <dd>Standard Cartesian coordinates.</dd>
      <dt><code>"Cylindrical"</code></dt>
      <dd>Cylindrical coordinates (enables <code>+R</code>/<code>-R</code> directions).</dd>
    </dl>
  </dd>
  <dt id="config-boundaries-lumpedport-r"><a href="#config-boundaries-lumpedport-r"><code>"R"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.0</code></span></dt>
  <dd>
    <p>Circuit resistance, Ω. Use with <code>&quot;L&quot;</code> and <code>&quot;C&quot;</code>; do not mix with surface parameters <code>&quot;Rs&quot;</code>, <code>&quot;Ls&quot;</code>, <code>&quot;Cs&quot;</code>.</p>
  </dd>
  <dt id="config-boundaries-lumpedport-l"><a href="#config-boundaries-lumpedport-l"><code>"L"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.0</code></span></dt>
  <dd>
    <p>Circuit inductance, H. Use with <code>&quot;R&quot;</code> and <code>&quot;C&quot;</code>; do not mix with surface parameters <code>&quot;Rs&quot;</code>, <code>&quot;Ls&quot;</code>, <code>&quot;Cs&quot;</code>.</p>
  </dd>
  <dt id="config-boundaries-lumpedport-c"><a href="#config-boundaries-lumpedport-c"><code>"C"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.0</code></span></dt>
  <dd>
    <p>Circuit capacitance, F. Use with <code>&quot;R&quot;</code> and <code>&quot;L&quot;</code>; do not mix with surface parameters <code>&quot;Rs&quot;</code>, <code>&quot;Ls&quot;</code>, <code>&quot;Cs&quot;</code>.</p>
  </dd>
  <dt id="config-boundaries-lumpedport-rs"><a href="#config-boundaries-lumpedport-rs"><code>"Rs"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.0</code></span></dt>
  <dd>
    <p>Surface resistance, Ω/sq. Use with <code>&quot;Ls&quot;</code> and <code>&quot;Cs&quot;</code>; do not mix with circuit parameters <code>&quot;R&quot;</code>, <code>&quot;L&quot;</code>, <code>&quot;C&quot;</code>.</p>
  </dd>
  <dt id="config-boundaries-lumpedport-ls"><a href="#config-boundaries-lumpedport-ls"><code>"Ls"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.0</code></span></dt>
  <dd>
    <p>Surface inductance, H/sq. Use with <code>&quot;Rs&quot;</code> and <code>&quot;Cs&quot;</code>; do not mix with circuit parameters <code>&quot;R&quot;</code>, <code>&quot;L&quot;</code>, <code>&quot;C&quot;</code>.</p>
  </dd>
  <dt id="config-boundaries-lumpedport-cs"><a href="#config-boundaries-lumpedport-cs"><code>"Cs"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.0</code></span></dt>
  <dd>
    <p>Surface capacitance, F/sq. Use with <code>&quot;Rs&quot;</code> and <code>&quot;Ls&quot;</code>; do not mix with circuit parameters <code>&quot;R&quot;</code>, <code>&quot;L&quot;</code>, <code>&quot;C&quot;</code>.</p>
  </dd>
  <dt id="config-boundaries-lumpedport-excitation"><a href="#config-boundaries-lumpedport-excitation"><code>"Excitation"</code></a> <span class="config-type">boolean or integer</span> <span class="config-default">default: <code>false</code></span></dt>
  <dd>
    <p>Turns on or off port excitation for driven or transient simulations. Can be specified as a boolean or as a non-negative integer (excitation group index). See the <a href="../../guide/boundaries/#Lumped-and-wave-port-excitation">boundary conditions guide</a> for details.</p>
  </dd>
  <dt id="config-boundaries-lumpedport-active"><a href="#config-boundaries-lumpedport-active"><code>"Active"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>true</code></span></dt>
  <dd>
    <p>Turns on or off the damping boundary condition for this port for driven or transient simulations.</p>
  </dd>
  <dt><code>"Elements"</code> <span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-boundaries-lumpedport-elements">See full reference ↓</a></p>
  </dd>
</dl>
```

#### [Elements](@id config-boundaries-lumpedport-elements)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-lumpedport">LumpedPort</a>/<a href="#config-boundaries-lumpedport-elements">Elements</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></span>
```

Sub-elements for a multielement lumped port. Use this instead of the top-level `"Attributes"`/`"Direction"`/`"CoordinateSystem"` when the port spans multiple disjoint boundary surfaces. Elements add in parallel.

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-lumpedport">LumpedPort</a>/<a href="#config-boundaries-lumpedport-elements">Elements</a>/0</code></p>
```

```@raw html
<dl class="palace-config">
  <dt id="config-boundaries-lumpedport-elements-attributes"><a href="#config-boundaries-lumpedport-elements-attributes"><code>"Attributes"</code></a> <span class="config-type">[integer, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Integer array of mesh boundary attributes this object applies to.</p>
  </dd>
  <dt id="config-boundaries-lumpedport-elements-direction"><a href="#config-boundaries-lumpedport-elements-direction"><code>"Direction"</code></a> <span class="config-type">string or [number × 3]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Excitation direction. Axis-aligned Cartesian directions can be specified using keywords: <code>&quot;+X&quot;</code>, <code>&quot;-X&quot;</code>, <code>&quot;+Y&quot;</code>, <code>&quot;-Y&quot;</code>, <code>&quot;+Z&quot;</code>, <code>&quot;-Z&quot;</code>. Coaxial directions use <code>&quot;+R&quot;</code>, <code>&quot;-R&quot;</code>. Alternatively, specify a normalized 3-element array, e.g. <code>[0.0, 1.0, 0.0]</code>. The coordinate system is determined by <code>&quot;CoordinateSystem&quot;</code>.</p>
  </dd>
  <dt id="config-boundaries-lumpedport-elements-coordinatesystem"><a href="#config-boundaries-lumpedport-elements-coordinatesystem"><code>"CoordinateSystem"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;Cartesian&quot;</code></span></dt>
  <dd>
    <p>Coordinate system for this element's <code>&quot;Direction&quot;</code> vector.</p>
    <dl class="config-enum">
      <dt><code>"Cartesian"</code></dt>
      <dd>Standard Cartesian coordinates.</dd>
      <dt><code>"Cylindrical"</code></dt>
      <dd>Cylindrical coordinates.</dd>
    </dl>
  </dd>
</dl>
```

### [Wave Port](@id config-boundaries-waveport)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-waveport">WavePort</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></span>
```

Array of numeric wave port boundary conditions. Wave ports can only be specified on the true boundary of the computational domain (they must be "one-sided"). A 2D boundary mode eigenproblem is solved on each wave port to compute the port mode shape. Only available for frequency domain driven and eigenmode simulations.

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-waveport">WavePort</a>/0</code></p>
```

```@raw html
<dl class="palace-config">
  <dt id="config-boundaries-waveport-index"><a href="#config-boundaries-waveport-index"><code>"Index"</code></a> <span class="config-type">integer</span> <span class="config-required">required</span> <span class="config-constraint"><code>&gt; 0</code></span></dt>
  <dd>
    <p>Index of this wave port, used in postprocessing output files. Must be unique across all port and source types.</p>
  </dd>
  <dt id="config-boundaries-waveport-attributes"><a href="#config-boundaries-waveport-attributes"><code>"Attributes"</code></a> <span class="config-type">[integer, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Integer array of mesh boundary attributes this object applies to.</p>
  </dd>
  <dt id="config-boundaries-waveport-mode"><a href="#config-boundaries-waveport-mode"><code>"Mode"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>1</code></span> <span class="config-constraint"><code>&gt; 0</code></span></dt>
  <dd>
    <p>Mode index (1-based) for the characteristic port mode of this wave port, ranked in order of decreasing wave number.</p>
  </dd>
  <dt id="config-boundaries-waveport-offset"><a href="#config-boundaries-waveport-offset"><code>"Offset"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.0</code></span> <span class="config-constraint"><code>≥ 0.0</code></span></dt>
  <dd>
    <p>Offset distance used for S-parameter de-embedding for this wave port, specified in mesh length units.</p>
  </dd>
  <dt id="config-boundaries-waveport-solvertype"><a href="#config-boundaries-waveport-solvertype"><code>"SolverType"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;Default&quot;</code></span></dt>
  <dd>
    <p>Eigenvalue solver for computing the boundary mode. Accepts the same options as <a href="#config-solver-eigenmode-type"><code>/Solver/Eigenmode/Type</code></a>.</p>
  </dd>
  <dt id="config-boundaries-waveport-excitation"><a href="#config-boundaries-waveport-excitation"><code>"Excitation"</code></a> <span class="config-type">boolean or integer</span> <span class="config-default">default: <code>false</code></span></dt>
  <dd>
    <p>Turns on or off port excitation for driven simulations. Can be specified as a boolean or as a non-negative integer (excitation group index). See the <a href="../../guide/boundaries/#Lumped-and-wave-port-excitation">boundary conditions guide</a> for details.</p>
  </dd>
  <dt id="config-boundaries-waveport-active"><a href="#config-boundaries-waveport-active"><code>"Active"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>true</code></span></dt>
  <dd>
    <p>Turns on or off the damping boundary condition for this port for driven simulations.</p>
  </dd>
  <dt id="config-boundaries-waveport-maxits"><a href="#config-boundaries-waveport-maxits"><code>"MaxIts"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>30</code></span> <span class="config-constraint"><code>&gt; 0</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Maximum number of iterations for the GMRES solver used in the wave port boundary mode analysis.</p>
  </dd>
  <dt id="config-boundaries-waveport-ksptol"><a href="#config-boundaries-waveport-ksptol"><code>"KSPTol"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>1.0e-8</code></span> <span class="config-constraint"><code>&gt; 0.0</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Tolerance for the linear solver used in the wave port boundary mode analysis.</p>
  </dd>
  <dt id="config-boundaries-waveport-eigentol"><a href="#config-boundaries-waveport-eigentol"><code>"EigenTol"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>1.0e-6</code></span> <span class="config-constraint"><code>&gt; 0.0</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Tolerance for the eigenvalue solver used in the wave port boundary mode analysis.</p>
  </dd>
  <dt id="config-boundaries-waveport-verbose"><a href="#config-boundaries-waveport-verbose"><code>"Verbose"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>0</code></span> <span class="config-constraint"><code>≥ 0</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Verbosity level for the wave port linear and eigensolvers.</p>
  </dd>
</dl>
```

### [Wave Port PEC](@id config-boundaries-waveportpec)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-waveportpec">WavePortPEC</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></span>
```

Additional PEC boundary conditions for the 2D eigensolve used in wave port mode analysis, along with those already specified under [PEC](#config-boundaries-pec) and [Conductivity](#config-boundaries-conductivity). Only relevant when [WavePort](#config-boundaries-waveport) boundaries are present.

```@raw html
<dl class="palace-config">
  <dt id="config-boundaries-waveportpec-attributes"><a href="#config-boundaries-waveportpec-attributes"><code>"Attributes"</code></a> <span class="config-type">[integer, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Integer array of mesh boundary attributes this object applies to.</p>
  </dd>
</dl>
```

### [Surface Current](@id config-boundaries-surfacecurrent)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-surfacecurrent">SurfaceCurrent</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></span>
```

Array of surface current source boundaries. Prescribes a unit source surface current excitation on the given boundary to excite a driven, transient, or magnetostatic simulation. For magnetostatic simulations, inductance matrix entries are extracted for each surface current boundary.

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-surfacecurrent">SurfaceCurrent</a>/0</code></p>
```

```@raw html
<dl class="palace-config">
  <dt id="config-boundaries-surfacecurrent-index"><a href="#config-boundaries-surfacecurrent-index"><code>"Index"</code></a> <span class="config-type">integer</span> <span class="config-required">required</span> <span class="config-constraint"><code>&gt; 0</code></span></dt>
  <dd>
    <p>Index of this surface current source, used in postprocessing output files. Must be unique across all port and source types.</p>
  </dd>
  <dt id="config-boundaries-surfacecurrent-attributes"><a href="#config-boundaries-surfacecurrent-attributes"><code>"Attributes"</code></a> <span class="config-type">[integer, ...]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p>Integer array of mesh boundary attributes for this surface current boundary. If this is to be a object with more than a single element, use the &quot;Elements&quot; array described below.</p>
  </dd>
  <dt id="config-boundaries-surfacecurrent-direction"><a href="#config-boundaries-surfacecurrent-direction"><code>"Direction"</code></a> <span class="config-type">string or [number × 3]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p>Excitation direction. Axis-aligned Cartesian directions can be specified using keywords: <code>&quot;+X&quot;</code>, <code>&quot;-X&quot;</code>, <code>&quot;+Y&quot;</code>, <code>&quot;-Y&quot;</code>, <code>&quot;+Z&quot;</code>, <code>&quot;-Z&quot;</code>. Coaxial directions use <code>&quot;+R&quot;</code>, <code>&quot;-R&quot;</code>. Alternatively, specify a normalized 3-element array, e.g. <code>[0.0, 1.0, 0.0]</code>. The coordinate system is determined by <code>&quot;CoordinateSystem&quot;</code>.</p>
  </dd>
  <dt id="config-boundaries-surfacecurrent-coordinatesystem"><a href="#config-boundaries-surfacecurrent-coordinatesystem"><code>"CoordinateSystem"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;Cartesian&quot;</code></span></dt>
  <dd>
    <p>Coordinate system for the <code>&quot;Direction&quot;</code> vector. Same options as <a href="#config-boundaries-lumpedport-coordinatesystem"><code>/LumpedPort/CoordinateSystem</code></a>.</p>
    <dl class="config-enum">
      <dt><code>"Cartesian"</code></dt>
      <dd>Standard Cartesian coordinates.</dd>
      <dt><code>"Cylindrical"</code></dt>
      <dd>Cylindrical coordinates.</dd>
    </dl>
  </dd>
  <dt><code>"Elements"</code> <span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-boundaries-surfacecurrent-elements">See full reference ↓</a></p>
  </dd>
</dl>
```

#### [Elements](@id config-boundaries-surfacecurrent-elements)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-surfacecurrent">SurfaceCurrent</a>/<a href="#config-boundaries-surfacecurrent-elements">Elements</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></span>
```

Sub-elements for a multielement surface current source. Use this instead of the top-level `"Attributes"`/`"Direction"`/`"CoordinateSystem"` when the source spans multiple disjoint boundary surfaces. Elements add in parallel to give the same total current as a single-element source.

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-surfacecurrent">SurfaceCurrent</a>/<a href="#config-boundaries-surfacecurrent-elements">Elements</a>/0</code></p>
```

```@raw html
<dl class="palace-config">
  <dt id="config-boundaries-surfacecurrent-elements-attributes"><a href="#config-boundaries-surfacecurrent-elements-attributes"><code>"Attributes"</code></a> <span class="config-type">[integer, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Integer array of mesh boundary attributes this object applies to.</p>
  </dd>
  <dt id="config-boundaries-surfacecurrent-elements-direction"><a href="#config-boundaries-surfacecurrent-elements-direction"><code>"Direction"</code></a> <span class="config-type">string or [number × 3]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Excitation direction. Axis-aligned Cartesian directions can be specified using keywords: <code>&quot;+X&quot;</code>, <code>&quot;-X&quot;</code>, <code>&quot;+Y&quot;</code>, <code>&quot;-Y&quot;</code>, <code>&quot;+Z&quot;</code>, <code>&quot;-Z&quot;</code>. Coaxial directions use <code>&quot;+R&quot;</code>, <code>&quot;-R&quot;</code>. Alternatively, specify a normalized 3-element array, e.g. <code>[0.0, 1.0, 0.0]</code>. The coordinate system is determined by <code>&quot;CoordinateSystem&quot;</code>.</p>
  </dd>
  <dt id="config-boundaries-surfacecurrent-elements-coordinatesystem"><a href="#config-boundaries-surfacecurrent-elements-coordinatesystem"><code>"CoordinateSystem"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;Cartesian&quot;</code></span></dt>
  <dd>
    <p>Coordinate system for this element's <code>&quot;Direction&quot;</code> vector.</p>
    <dl class="config-enum">
      <dt><code>"Cartesian"</code></dt>
      <dd>Standard Cartesian coordinates.</dd>
      <dt><code>"Cylindrical"</code></dt>
      <dd>Cylindrical coordinates.</dd>
    </dl>
  </dd>
</dl>
```

### [Ground Boundary](@id config-boundaries-ground)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-ground">Ground</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></span>
```

Zero-voltage (ground) boundary condition for electrostatic simulations. Mutually exclusive with [PEC](#config-boundaries-pec).

```@raw html
<dl class="palace-config">
  <dt id="config-boundaries-ground-attributes"><a href="#config-boundaries-ground-attributes"><code>"Attributes"</code></a> <span class="config-type">[integer, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Integer array of mesh boundary attributes this object applies to.</p>
  </dd>
</dl>
```

### [Zero Charge Boundary](@id config-boundaries-zerocharge)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-zerocharge">ZeroCharge</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></span>
```

Zero surface charge (homogeneous Neumann) boundary condition for electrostatic simulations. Also imposes symmetry of the electric field across the surface. Mutually exclusive with [PMC](#config-boundaries-pmc).

```@raw html
<dl class="palace-config">
  <dt id="config-boundaries-zerocharge-attributes"><a href="#config-boundaries-zerocharge-attributes"><code>"Attributes"</code></a> <span class="config-type">[integer, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Integer array of mesh boundary attributes this object applies to.</p>
  </dd>
</dl>
```

### [Terminal](@id config-boundaries-terminal)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-terminal">Terminal</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></span>
```

Array of terminal boundaries for electrostatic simulations. Capacitance matrix entries are extracted for each terminal.

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-terminal">Terminal</a>/0</code></p>
```

```@raw html
<dl class="palace-config">
  <dt id="config-boundaries-terminal-index"><a href="#config-boundaries-terminal-index"><code>"Index"</code></a> <span class="config-type">integer</span> <span class="config-required">required</span> <span class="config-constraint"><code>&gt; 0</code></span></dt>
  <dd>
    <p>Index of this terminal, used in postprocessing output files and to index the computed capacitance matrix.</p>
  </dd>
  <dt id="config-boundaries-terminal-attributes"><a href="#config-boundaries-terminal-attributes"><code>"Attributes"</code></a> <span class="config-type">[integer, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Integer array of mesh boundary attributes this object applies to.</p>
  </dd>
</dl>
```

### [Periodic Boundary](@id config-boundaries-periodic)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-periodic">Periodic</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></span>
```

Periodic boundary conditions for surfaces whose meshes are identical after translation and/or rotation. Floquet periodic boundary conditions with a phase shift are also supported.

```@raw html
<dl class="palace-config">
  <dt><code>"BoundaryPairs"</code> <span class="config-type">[object, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p><a href="#config-boundaries-periodic-boundarypairs">See full reference ↓</a></p>
  </dd>
  <dt id="config-boundaries-periodic-floquetwavevector"><a href="#config-boundaries-periodic-floquetwavevector"><code>"FloquetWaveVector"</code></a> <span class="config-type">[number × 3]</span> <span class="config-default">default: <code>[0.0,0.0,0.0]</code></span></dt>
  <dd>
    <p>3-element Floquet wave vector <code>[kx, ky, kz]</code> defining the phase delay between periodic boundaries, in radians per mesh length unit.</p>
  </dd>
</dl>
```

#### [BoundaryPairs](@id config-boundaries-periodic-boundarypairs)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-periodic">Periodic</a>/<a href="#config-boundaries-periodic-boundarypairs">BoundaryPairs</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">[object, ...]</span> <span class="config-required">required</span></span>
```

Array of donor–receiver boundary pairs defining the periodic mapping.

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-periodic">Periodic</a>/<a href="#config-boundaries-periodic-boundarypairs">BoundaryPairs</a>/0</code></p>
```

```@raw html
<dl class="palace-config">
  <dt id="config-boundaries-periodic-boundarypairs-donorattributes"><a href="#config-boundaries-periodic-boundarypairs-donorattributes"><code>"DonorAttributes"</code></a> <span class="config-type">[integer, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Integer array of donor mesh boundary attributes for a periodic boundary pair.</p>
  </dd>
  <dt id="config-boundaries-periodic-boundarypairs-receiverattributes"><a href="#config-boundaries-periodic-boundarypairs-receiverattributes"><code>"ReceiverAttributes"</code></a> <span class="config-type">[integer, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Integer array of receiver mesh boundary attributes for a periodic boundary pair.</p>
  </dd>
  <dt id="config-boundaries-periodic-boundarypairs-translation"><a href="#config-boundaries-periodic-boundarypairs-translation"><code>"Translation"</code></a> <span class="config-type">[number × 3]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p>3-element translation vector <code>[dx, dy, dz]</code> from donor to receiver boundary, in mesh length units. If neither <code>&quot;Translation&quot;</code> nor <code>&quot;AffineTransformation&quot;</code> are specified, the transformation is detected automatically.</p>
  </dd>
  <dt id="config-boundaries-periodic-boundarypairs-affinetransformation"><a href="#config-boundaries-periodic-boundarypairs-affinetransformation"><code>"AffineTransformation"</code></a> <span class="config-type">[number × 16]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p>16-element row-major 4×4 affine transformation matrix from donor to receiver boundary, in mesh length units. If neither <code>&quot;Translation&quot;</code> nor <code>&quot;AffineTransformation&quot;</code> are specified, the transformation is detected automatically.</p>
  </dd>
</dl>
```

### [Boundary Postprocessing](@id config-boundaries-postprocessing)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-postprocessing">Postprocessing</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></span>
```

Configuration for boundary postprocessing.

```@raw html
<dl class="palace-config">
  <dt><code>"SurfaceFlux"</code> <span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-boundaries-postprocessing-surfaceflux">See full reference ↓</a></p>
  </dd>
  <dt><code>"Dielectric"</code> <span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-boundaries-postprocessing-dielectric">See full reference ↓</a></p>
  </dd>
  <dt><code>"FarField"</code> <span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-boundaries-postprocessing-farfield">See full reference ↓</a></p>
  </dd>
</dl>
```

#### [Surface Flux](@id config-boundaries-postprocessing-surfaceflux)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-postprocessing">Postprocessing</a>/<a href="#config-boundaries-postprocessing-surfaceflux">SurfaceFlux</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></span>
```

Array of surface flux postprocessing boundaries. Results are written to `surface-F.csv` in the output directory.

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-postprocessing">Postprocessing</a>/<a href="#config-boundaries-postprocessing-surfaceflux">SurfaceFlux</a>/0</code></p>
```

```@raw html
<dl class="palace-config">
  <dt id="config-boundaries-postprocessing-surfaceflux-index"><a href="#config-boundaries-postprocessing-surfaceflux-index"><code>"Index"</code></a> <span class="config-type">integer</span> <span class="config-required">required</span> <span class="config-constraint"><code>&gt; 0</code></span></dt>
  <dd>
    <p>Index of this surface flux postprocessing boundary, used in output files.</p>
  </dd>
  <dt id="config-boundaries-postprocessing-surfaceflux-attributes"><a href="#config-boundaries-postprocessing-surfaceflux-attributes"><code>"Attributes"</code></a> <span class="config-type">[integer, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Integer array of mesh boundary attributes this object applies to.</p>
  </dd>
  <dt id="config-boundaries-postprocessing-surfaceflux-type"><a href="#config-boundaries-postprocessing-surfaceflux-type"><code>"Type"</code></a> <span class="config-type">string</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Type of surface flux to integrate over the boundary.</p>
    <dl class="config-enum">
      <dt><code>"Electric"</code></dt>
      <dd>Integrate the electric flux density.</dd>
      <dt><code>"Magnetic"</code></dt>
      <dd>Integrate the magnetic flux density.</dd>
      <dt><code>"Power"</code></dt>
      <dd>Integrate the Poynting vector (energy flux).</dd>
    </dl>
  </dd>
  <dt id="config-boundaries-postprocessing-surfaceflux-twosided"><a href="#config-boundaries-postprocessing-surfaceflux-twosided"><code>"TwoSided"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span></dt>
  <dd>
    <p>For internal boundary surfaces: when <code>false</code>, the flux on both sides is averaged; when <code>true</code>, it is summed with opposite normal direction.</p>
  </dd>
  <dt id="config-boundaries-postprocessing-surfaceflux-center"><a href="#config-boundaries-postprocessing-surfaceflux-center"><code>"Center"</code></a> <span class="config-type">[number × 3]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p>Point used to determine the outward normal orientation, in mesh length units. Only used when <code>&quot;TwoSided&quot;</code> is <code>false</code>. If not specified, the point will be computed as the centroid of the axis-aligned bounding box for all elements making up the postprocessing boundary.</p>
  </dd>
</dl>
```

#### [Dielectric](@id config-boundaries-postprocessing-dielectric)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-postprocessing">Postprocessing</a>/<a href="#config-boundaries-postprocessing-dielectric">Dielectric</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></span>
```

Array of interface dielectric loss postprocessing boundaries. Computes energy participation ratios (EPR) and quality factors for dielectric interfaces. See also the [reference documentation](../reference.md#Bulk-and-interface-dielectric-loss). Results are written to `surface-Q.csv` in the output directory.

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-postprocessing">Postprocessing</a>/<a href="#config-boundaries-postprocessing-dielectric">Dielectric</a>/0</code></p>
```

```@raw html
<dl class="palace-config">
  <dt id="config-boundaries-postprocessing-dielectric-index"><a href="#config-boundaries-postprocessing-dielectric-index"><code>"Index"</code></a> <span class="config-type">integer</span> <span class="config-required">required</span> <span class="config-constraint"><code>&gt; 0</code></span></dt>
  <dd>
    <p>Index of this dielectric interface, used in output files.</p>
  </dd>
  <dt id="config-boundaries-postprocessing-dielectric-attributes"><a href="#config-boundaries-postprocessing-dielectric-attributes"><code>"Attributes"</code></a> <span class="config-type">[integer, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Integer array of mesh boundary attributes this object applies to.</p>
  </dd>
  <dt id="config-boundaries-postprocessing-dielectric-thickness"><a href="#config-boundaries-postprocessing-dielectric-thickness"><code>"Thickness"</code></a> <span class="config-type">number</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Thickness of this dielectric interface, in mesh length units.</p>
  </dd>
  <dt id="config-boundaries-postprocessing-dielectric-permittivity"><a href="#config-boundaries-postprocessing-dielectric-permittivity"><code>"Permittivity"</code></a> <span class="config-type">number</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Relative permittivity of this dielectric interface layer. This should be the interface layer permittivity for the specific &quot;Type&quot; of interface specified.</p>
  </dd>
  <dt id="config-boundaries-postprocessing-dielectric-type"><a href="#config-boundaries-postprocessing-dielectric-type"><code>"Type"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;Default&quot;</code></span></dt>
  <dd>
    <p>Interface type used to determine the boundary conditions for computing the EPR. See also <a href="../../reference/#Bulk-and-interface-dielectric-loss">theory reference</a>.</p>
    <dl class="config-enum">
      <dt><code>"Default"</code></dt>
      <dd>Use the full electric field evaluated at the boundary.</dd>
      <dt><code>"MA"</code></dt>
      <dd>Use metal-air interface boundary conditions.</dd>
      <dt><code>"MS"</code></dt>
      <dd>Use metal-substrate interface boundary conditions.</dd>
      <dt><code>"SA"</code></dt>
      <dd>Use substrate-air interface boundary conditions.</dd>
    </dl>
  </dd>
  <dt id="config-boundaries-postprocessing-dielectric-losstan"><a href="#config-boundaries-postprocessing-dielectric-losstan"><code>"LossTan"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.0</code></span></dt>
  <dd>
    <p>Loss tangent of this dielectric interface.</p>
  </dd>
</dl>
```

#### [Far Field](@id config-boundaries-postprocessing-farfield)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-boundaries">Boundaries</a>/<a href="#config-boundaries-postprocessing">Postprocessing</a>/<a href="#config-boundaries-postprocessing-farfield">FarField</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></span>
```

Far-field electric field extraction. The boundary attributes must enclose the system and be on an external boundary.

```@raw html
<dl class="palace-config">
  <dt id="config-boundaries-postprocessing-farfield-attributes"><a href="#config-boundaries-postprocessing-farfield-attributes"><code>"Attributes"</code></a> <span class="config-type">[integer, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Integer array of mesh boundary attributes this object applies to.</p>
  </dd>
  <dt id="config-boundaries-postprocessing-farfield-nsample"><a href="#config-boundaries-postprocessing-farfield-nsample"><code>"NSample"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>0</code></span> <span class="config-constraint"><code>≥ 0</code></span></dt>
  <dd>
    <p>Number of uniformly-spaced points used to discretize the far-field sphere.</p>
  </dd>
  <dt id="config-boundaries-postprocessing-farfield-thetaphis"><a href="#config-boundaries-postprocessing-farfield-thetaphis"><code>"ThetaPhis"</code></a> <span class="config-type">[[number × 2], ...]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p>Additional specific (θ, φ) angle pairs in degrees at which to evaluate the far field. θ ∈ [0°, 180°] is the polar angle, φ ∈ [0°, 360°] is the azimuthal angle.</p>
  </dd>
</dl>
```

## [Solver](@id config-solver)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-solver">Solver</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-required">required</span></span>
```

Solver configuration for all simulation types.

```@raw html
<dl class="palace-config">
  <dt id="config-solver-order"><a href="#config-solver-order"><code>"Order"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>1</code></span> <span class="config-constraint"><code>≥ 1</code></span></dt>
  <dd>
    <p>Finite element order (degree). Arbitrary high-order spaces are supported.</p>
  </dd>
  <dt id="config-solver-partialassemblyorder"><a href="#config-solver-partialassemblyorder"><code>"PartialAssemblyOrder"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>1</code></span> <span class="config-constraint"><code>≥ 1</code></span></dt>
  <dd>
    <p>Order at which to switch from full assembly of finite element operators to <a href="https://mfem.org/howto/assembly_levels/">partial assembly</a>. Setting to <code>1</code> fully activates partial assembly on all levels; a large value (greater than <a href="#config-solver-order">Order</a>) results in fully assembled sparse matrix operators.</p>
  </dd>
  <dt id="config-solver-quadratureorderjacobian"><a href="#config-solver-quadratureorderjacobian"><code>"QuadratureOrderJacobian"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Use the Jacobian-based quadrature order instead of the default.</p>
  </dd>
  <dt id="config-solver-quadratureorderextra"><a href="#config-solver-quadratureorderextra"><code>"QuadratureOrderExtra"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>0</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Extra quadrature order added on top of the default.</p>
  </dd>
  <dt id="config-solver-device"><a href="#config-solver-device"><code>"Device"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;CPU&quot;</code></span></dt>
  <dd>
    <p>The runtime device configuration passed to <a href="https://mfem.org/howto/assembly_levels/">MFEM</a> to activate different computation backends. When <em>Palace</em> is built with OpenMP support (<code>PALACE_WITH_OPENMP=ON</code>), <code>omp</code> is automatically added to the MFEM device list.</p>
    <dl class="config-enum">
      <dt><code>"CPU"</code></dt>
      <dd>Run on CPU.</dd>
      <dt><code>"GPU"</code></dt>
      <dd>Run on GPU via CUDA (<code>MFEM_USE_CUDA=ON</code>) or HIP (<code>MFEM_USE_HIP=ON</code>).</dd>
      <dt><code>"Debug"</code></dt>
      <dd>MFEM debug device, useful for diagnosing GPU-related issues.</dd>
    </dl>
  </dd>
  <dt id="config-solver-backend"><a href="#config-solver-backend"><code>"Backend"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;&quot;</code></span></dt>
  <dd>
    <p>Specifies the <a href="https://libceed.org/en/latest/gettingstarted/#backends">libCEED backend</a> to use. If empty, a suitable default is selected based on <a href="#config-solver-device">Device</a>.</p>
  </dd>
  <dt><code>"Eigenmode"</code> <span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-solver-eigenmode">See full reference ↓</a></p>
  </dd>
  <dt><code>"Driven"</code> <span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-solver-driven">See full reference ↓</a></p>
  </dd>
  <dt><code>"Transient"</code> <span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-solver-transient">See full reference ↓</a></p>
  </dd>
  <dt><code>"Electrostatic"</code> <span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-solver-electrostatic">See full reference ↓</a></p>
  </dd>
  <dt><code>"Magnetostatic"</code> <span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-solver-magnetostatic">See full reference ↓</a></p>
  </dd>
  <dt><code>"Linear"</code> <span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-solver-linear">See full reference ↓</a></p>
  </dd>
</dl>
```

### [Eigenmode Solver](@id config-solver-eigenmode)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-solver">Solver</a>/<a href="#config-solver-eigenmode">Eigenmode</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></span>
```

Configuration for the eigenvalue solver. Only relevant when [`/Problem/Type`](#config-problem-type) is `"Eigenmode"`.

```@raw html
<dl class="palace-config">
  <dt id="config-solver-eigenmode-target"><a href="#config-solver-eigenmode-target"><code>"Target"</code></a> <span class="config-type">number</span> <span class="config-required">required</span> <span class="config-constraint"><code>&gt; 0.0</code></span></dt>
  <dd>
    <p>(Nonzero) frequency target above which to search for eigenvalues, GHz.</p>
  </dd>
  <dt id="config-solver-eigenmode-tol"><a href="#config-solver-eigenmode-tol"><code>"Tol"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>1.0e-6</code></span> <span class="config-constraint"><code>≥ 0.0</code></span></dt>
  <dd>
    <p>Relative convergence tolerance for the eigenvalue solver.</p>
  </dd>
  <dt id="config-solver-eigenmode-maxits"><a href="#config-solver-eigenmode-maxits"><code>"MaxIts"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>0</code></span> <span class="config-constraint"><code>&gt; 0</code></span></dt>
  <dd>
    <p>Maximum number of iterations for the iterative eigenvalue solver. A value less than 1 uses the solver default.</p>
  </dd>
  <dt id="config-solver-eigenmode-maxsize"><a href="#config-solver-eigenmode-maxsize"><code>"MaxSize"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>0</code></span> <span class="config-constraint"><code>&gt; 0</code></span></dt>
  <dd>
    <p>Maximum subspace dimension for the eigenvalue solver. A value less than 1 uses the solver default.</p>
  </dd>
  <dt id="config-solver-eigenmode-n"><a href="#config-solver-eigenmode-n"><code>"N"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>1</code></span> <span class="config-constraint"><code>&gt; 0</code></span></dt>
  <dd>
    <p>Number of eigenvalues to compute.</p>
  </dd>
  <dt id="config-solver-eigenmode-save"><a href="#config-solver-eigenmode-save"><code>"Save"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>0</code></span></dt>
  <dd>
    <p>Number of computed field modes to save to disk for <a href="../../guide/postprocessing/#Visualization">visualization with ParaView</a>. Files are saved in the <code>paraview/</code> (and/or <code>gridfunction/</code>) directory under <a href="#config-problem-output"><code>/Problem/Output</code></a>.</p>
  </dd>
  <dt id="config-solver-eigenmode-type"><a href="#config-solver-eigenmode-type"><code>"Type"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;Default&quot;</code></span></dt>
  <dd>
    <p>Specifies the eigenvalue solver backend.</p>
    <dl class="config-enum">
      <dt><code>"SLEPc"</code></dt>
      <dd>Krylov-Schur solver from SLEPc.</dd>
      <dt><code>"ARPACK"</code></dt>
      <dd>ARPACK eigensolver.</dd>
      <dt><code>"Default"</code></dt>
      <dd>Use the default solver (currently SLEPc Krylov-Schur).</dd>
    </dl>
  </dd>
  <dt id="config-solver-eigenmode-nonlineartype"><a href="#config-solver-eigenmode-nonlineartype"><code>"NonlinearType"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;Hybrid&quot;</code></span></dt>
  <dd>
    <p>Specifies the nonlinear eigenvalue solver for problems with frequency-dependent boundary conditions.</p>
    <dl class="config-enum">
      <dt><code>"Hybrid"</code></dt>
      <dd>Hybrid algorithm: solve a polynomial (quadratic) approximation first, then refine with a quasi-Newton nonlinear eigensolver.</dd>
      <dt><code>"SLP"</code></dt>
      <dd>SLEPc Successive Linear Problem (SLP) nonlinear eigensolver.</dd>
    </dl>
  </dd>
  <dt id="config-solver-eigenmode-targetupper"><a href="#config-solver-eigenmode-targetupper"><code>"TargetUpper"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>-1</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Upper frequency bound for the eigenvalue search, GHz. Only used for nonlinear problems. A value less than or equal to zero uses <code>3 × Target</code> automatically. An inaccurate upper bound can negatively affect convergence of the nonlinear eigensolver.</p>
  </dd>
  <dt id="config-solver-eigenmode-peplinear"><a href="#config-solver-eigenmode-peplinear"><code>"PEPLinear"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>true</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Linearize the polynomial eigenvalue problem before solving.</p>
  </dd>
  <dt id="config-solver-eigenmode-scaling"><a href="#config-solver-eigenmode-scaling"><code>"Scaling"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>true</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Enable scaling of the eigenvalue problem.</p>
  </dd>
  <dt id="config-solver-eigenmode-startvector"><a href="#config-solver-eigenmode-startvector"><code>"StartVector"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>true</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Use a random start vector for the eigensolver.</p>
  </dd>
  <dt id="config-solver-eigenmode-startvectorconstant"><a href="#config-solver-eigenmode-startvectorconstant"><code>"StartVectorConstant"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Use a constant start vector instead of a random one.</p>
  </dd>
  <dt id="config-solver-eigenmode-massorthogonal"><a href="#config-solver-eigenmode-massorthogonal"><code>"MassOrthogonal"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Orthogonalize eigenvectors with respect to the mass matrix.</p>
  </dd>
  <dt id="config-solver-eigenmode-refinenonlinear"><a href="#config-solver-eigenmode-refinenonlinear"><code>"RefineNonlinear"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>true</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Refine nonlinear eigenvalue solutions after the initial solve.</p>
  </dd>
  <dt id="config-solver-eigenmode-lineartol"><a href="#config-solver-eigenmode-lineartol"><code>"LinearTol"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.001</code></span> <span class="config-constraint"><code>≥ 0.0</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Tolerance for the inner linear solve within the nonlinear eigensolver.</p>
  </dd>
  <dt id="config-solver-eigenmode-preconditionerlag"><a href="#config-solver-eigenmode-preconditionerlag"><code>"PreconditionerLag"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>10</code></span> <span class="config-constraint"><code>≥ 0</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Number of eigensolver iterations between preconditioner updates.</p>
  </dd>
  <dt id="config-solver-eigenmode-preconditionerlagtol"><a href="#config-solver-eigenmode-preconditionerlagtol"><code>"PreconditionerLagTol"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.0001</code></span> <span class="config-constraint"><code>≥ 0.0</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Residual tolerance below which preconditioner updates are skipped.</p>
  </dd>
  <dt id="config-solver-eigenmode-maxrestart"><a href="#config-solver-eigenmode-maxrestart"><code>"MaxRestart"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>2</code></span> <span class="config-constraint"><code>≥ 0</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Maximum number of restarts for the eigensolver.</p>
  </dd>
</dl>
```

### [Driven Solver](@id config-solver-driven)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-solver">Solver</a>/<a href="#config-solver-driven">Driven</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></span>
```

Configuration for the frequency domain driven solver. Only relevant when [`/Problem/Type`](#config-problem-type) is `"Driven"`.

```@raw html
<dl class="palace-config">
  <dt id="config-solver-driven-minfreq"><a href="#config-solver-driven-minfreq"><code>"MinFreq"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>—</code></span> <span class="config-deprecated">deprecated</span></dt>
  <dd>
    <p>Lower bound of the frequency sweep interval, GHz. Deprecated: use <a href="#config-solver-driven-samples-linear"><code>Linear Samples</code></a> interface instead.</p>
  </dd>
  <dt id="config-solver-driven-maxfreq"><a href="#config-solver-driven-maxfreq"><code>"MaxFreq"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>—</code></span> <span class="config-deprecated">deprecated</span></dt>
  <dd>
    <p>Upper bound of the frequency sweep interval, GHz. Deprecated: use <a href="#config-solver-driven-samples-linear"><code>Linear Samples</code></a> interface instead.</p>
  </dd>
  <dt id="config-solver-driven-freqstep"><a href="#config-solver-driven-freqstep"><code>"FreqStep"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>—</code></span> <span class="config-deprecated">deprecated</span></dt>
  <dd>
    <p>Frequency step size for the frequency sweep, GHz. Deprecated: use <a href="#config-solver-driven-samples-linear"><code>Linear Samples</code></a> interface instead.</p>
  </dd>
  <dt id="config-solver-driven-savestep"><a href="#config-solver-driven-savestep"><code>"SaveStep"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>0</code></span> <span class="config-deprecated">deprecated</span></dt>
  <dd>
    <p>Controls how often, in number of frequency steps, to save computed fields to disk for <a href="../../guide/postprocessing/#Visualization">visualization with ParaView</a>. Files are saved in the <code>paraview/</code> (and/or <code>gridfunction/</code>) directory under <a href="#config-problem-output"><code>/Problem/Output</code></a>. Deprecated: use <a href="#config-solver-driven-samples-linear"><code>Linear Samples</code></a> interface instead.</p>
  </dd>
  <dt><code>"Samples"</code> <span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p><a href="#config-solver-driven-samples">See full reference ↓</a></p>
  </dd>
  <dt id="config-solver-driven-save"><a href="#config-solver-driven-save"><code>"Save"</code></a> <span class="config-type">[number, ...]</span> <span class="config-default">default: <code>—</code></span></dt>
  <dd>
    <p>Additional frequencies at which to save computed fields to disk for <a href="../../guide/postprocessing/#Visualization">visualization with ParaView</a>, in addition to those specified by <code>&quot;SaveStep&quot;</code>. Files are saved in the <code>paraview/</code> (and/or <code>gridfunction/</code>) directory under <a href="#config-problem-output"><code>/Problem/Output</code></a>.</p>
  </dd>
  <dt id="config-solver-driven-restart"><a href="#config-solver-driven-restart"><code>"Restart"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>1</code></span> <span class="config-constraint"><code>&gt; 0</code></span></dt>
  <dd>
    <p>1-based sample index from which to restart a partial frequency sweep (i.e. <code>&quot;Restart&quot;: x</code> starts from the <em>x</em>-th sample of the combined sample set). Not valid for adaptive sweep.</p>
  </dd>
  <dt id="config-solver-driven-adaptivetol"><a href="#config-solver-driven-adaptivetol"><code>"AdaptiveTol"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.0</code></span> <span class="config-constraint"><code>≥ 0.0</code></span></dt>
  <dd>
    <p>Relative error convergence tolerance for adaptive fast frequency sweep. A value of <code>0</code> disables adaptive sweep and the full-order model is solved at each frequency step. A positive value ensures the reduced-order model is reliable relative to the full-order model in the frequency band of interest.</p>
  </dd>
  <dt id="config-solver-driven-adaptivemaxsamples"><a href="#config-solver-driven-adaptivemaxsamples"><code>"AdaptiveMaxSamples"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>20</code></span> <span class="config-constraint"><code>&gt; 0</code></span></dt>
  <dd>
    <p>Maximum number of frequency samples used to construct the reduced-order model for adaptive sweep, if the tolerance (<code>&quot;AdaptiveTol&quot;</code>) is not met first. In simulations with multiple excitations, this is the maximum per excitation.</p>
  </dd>
  <dt id="config-solver-driven-adaptiveconvergencememory"><a href="#config-solver-driven-adaptiveconvergencememory"><code>"AdaptiveConvergenceMemory"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>2</code></span> <span class="config-constraint"><code>&gt; 0</code></span></dt>
  <dd>
    <p>Number of consecutive samples satisfying the error tolerance required to declare convergence of the adaptive sampling algorithm.</p>
  </dd>
  <dt id="config-solver-driven-adaptivegsorthogonalization"><a href="#config-solver-driven-adaptivegsorthogonalization"><code>"AdaptiveGSOrthogonalization"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;CGS2&quot;</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Gram-Schmidt variant for orthogonalizing the adaptive reduced-order model basis. Same options as <a href="#config-solver-linear-gsorthogonalization"><code>/Linear/GSOrthogonalization</code></a>.</p>
    <dl class="config-enum">
      <dt><code>"MGS"</code></dt>
      <dd>Modified Gram-Schmidt.</dd>
      <dt><code>"CGS"</code></dt>
      <dd>Classical Gram-Schmidt.</dd>
      <dt><code>"CGS2"</code></dt>
      <dd>Two-step classical Gram-Schmidt with reorthogonalization.</dd>
    </dl>
  </dd>
  <dt id="config-solver-driven-adaptivecircuitsynthesis"><a href="#config-solver-driven-adaptivecircuitsynthesis"><code>"AdaptiveCircuitSynthesis"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span></dt>
  <dd>
    <p>Use the adaptive reduced-order model to print synthesized circuit-like matrices (L⁻¹, R⁻¹, C). Requires adaptive sweep to be enabled, all <code>LumpedPort</code> fields to be orthogonal, and only LRC-type frequency dependence (no <code>WavePort</code>, <code>Conductivity</code>, or second-order farfield BCs).</p>
  </dd>
  <dt id="config-solver-driven-adaptivecircuitsynthesisdomainorthogonalization"><a href="#config-solver-driven-adaptivecircuitsynthesisdomainorthogonalization"><code>"AdaptiveCircuitSynthesisDomainOrthogonalization"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;Energy&quot;</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Weight matrix type for domain orthogonalization when building synthesized circuit matrices.</p>
    <dl class="config-enum">
      <dt><code>"Energy"</code></dt>
      <dd>Use the energy-based domain mass matrix.</dd>
      <dt><code>"FEBasisIdentity"</code></dt>
      <dd>Use the identity matrix in the finite element basis.</dd>
      <dt><code>"SpaceOverlap"</code></dt>
      <dd>Use the physical-space field overlap.</dd>
    </dl>
  </dd>
</dl>
```

#### [Samples](@id config-solver-driven-samples)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-solver">Solver</a>/<a href="#config-solver-driven">Driven</a>/<a href="#config-solver-driven-samples">Samples</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">[object, ...]</span> <span class="config-default">default: <code>—</code></span></span>
```

Array of frequency sample specifications. Combined with `"MinFreq"`/`"MaxFreq"`/`"FreqStep"` to form a sorted, unique set of samples.

##### [Point Samples](@id config-solver-driven-samples-point)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-solver">Solver</a>/<a href="#config-solver-driven">Driven</a>/<a href="#config-solver-driven-samples">Samples</a>/<a href="#config-solver-driven-samples-point">0</a></code></p>
```

Explicit list of frequency sample points.

```@raw html
<dl class="palace-config">
  <dt id="config-solver-driven-samples-point-freq"><a href="#config-solver-driven-samples-point-freq"><code>"Freq"</code></a> <span class="config-type">[number, ...]</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Explicit frequencies to sample, GHz.</p>
  </dd>
  <dt id="config-solver-driven-samples-point-type"><a href="#config-solver-driven-samples-point-type"><code>"Type"</code></a> <span class="config-type">any</span> <span class="config-default">default: <code>&quot;Point&quot;</code></span></dt>
  <dd>
  </dd>
  <dt id="config-solver-driven-samples-point-savestep"><a href="#config-solver-driven-samples-point-savestep"><code>"SaveStep"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>0</code></span></dt>
  <dd>
    <p>Save fields every N steps within this sample. <code>0</code> disables saving.</p>
  </dd>
  <dt id="config-solver-driven-samples-point-addtoprom"><a href="#config-solver-driven-samples-point-addtoprom"><code>"AddToPROM"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span></dt>
  <dd>
    <p>Force inclusion of these points in the PROM for adaptive sweep (primarily a debugging tool).</p>
  </dd>
</dl>
```

##### [Linear Samples](@id config-solver-driven-samples-linear)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-solver">Solver</a>/<a href="#config-solver-driven">Driven</a>/<a href="#config-solver-driven-samples">Samples</a>/<a href="#config-solver-driven-samples-linear">1</a></code></p>
```

Linearly-spaced frequency samples.

```@raw html
<dl class="palace-config">
  <dt id="config-solver-driven-samples-linear-minfreq"><a href="#config-solver-driven-samples-linear-minfreq"><code>"MinFreq"</code></a> <span class="config-type">number</span> <span class="config-required">required</span> <span class="config-constraint"><code>≥ 0.0</code></span></dt>
  <dd>
    <p>Lower bound, GHz.</p>
  </dd>
  <dt id="config-solver-driven-samples-linear-maxfreq"><a href="#config-solver-driven-samples-linear-maxfreq"><code>"MaxFreq"</code></a> <span class="config-type">number</span> <span class="config-required">required</span> <span class="config-constraint"><code>≥ 0.0</code></span></dt>
  <dd>
    <p>Upper bound, GHz.</p>
  </dd>
  <dt id="config-solver-driven-samples-linear-type"><a href="#config-solver-driven-samples-linear-type"><code>"Type"</code></a> <span class="config-type">any</span> <span class="config-default">default: <code>&quot;Linear&quot;</code></span></dt>
  <dd>
  </dd>
  <dt id="config-solver-driven-samples-linear-freqstep"><a href="#config-solver-driven-samples-linear-freqstep"><code>"FreqStep"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>—</code></span> <span class="config-constraint"><code>≥ 0.0</code></span></dt>
  <dd>
    <p>Step size, GHz. Mutually exclusive with <code>&quot;NSample&quot;</code>.</p>
  </dd>
  <dt id="config-solver-driven-samples-linear-nsample"><a href="#config-solver-driven-samples-linear-nsample"><code>"NSample"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>—</code></span> <span class="config-constraint"><code>≥ 0</code></span></dt>
  <dd>
    <p>Number of samples. Mutually exclusive with <code>&quot;FreqStep&quot;</code>.</p>
  </dd>
  <dt id="config-solver-driven-samples-linear-savestep"><a href="#config-solver-driven-samples-linear-savestep"><code>"SaveStep"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>0</code></span></dt>
  <dd>
    <p>Save fields every N steps. <code>0</code> disables saving.</p>
  </dd>
  <dt id="config-solver-driven-samples-linear-addtoprom"><a href="#config-solver-driven-samples-linear-addtoprom"><code>"AddToPROM"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span></dt>
  <dd>
    <p>Force inclusion in the PROM for adaptive sweep.</p>
  </dd>
</dl>
```

##### [Log Samples](@id config-solver-driven-samples-log)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-solver">Solver</a>/<a href="#config-solver-driven">Driven</a>/<a href="#config-solver-driven-samples">Samples</a>/<a href="#config-solver-driven-samples-log">2</a></code></p>
```

Logarithmically-spaced frequency samples.

```@raw html
<dl class="palace-config">
  <dt id="config-solver-driven-samples-log-type"><a href="#config-solver-driven-samples-log-type"><code>"Type"</code></a> <span class="config-type">any</span> <span class="config-required">required</span></dt>
  <dd>
  </dd>
  <dt id="config-solver-driven-samples-log-minfreq"><a href="#config-solver-driven-samples-log-minfreq"><code>"MinFreq"</code></a> <span class="config-type">number</span> <span class="config-required">required</span> <span class="config-constraint"><code>&gt; 0.0</code></span></dt>
  <dd>
    <p>Lower bound, GHz.</p>
  </dd>
  <dt id="config-solver-driven-samples-log-maxfreq"><a href="#config-solver-driven-samples-log-maxfreq"><code>"MaxFreq"</code></a> <span class="config-type">number</span> <span class="config-required">required</span> <span class="config-constraint"><code>&gt; 0.0</code></span></dt>
  <dd>
    <p>Upper bound, GHz.</p>
  </dd>
  <dt id="config-solver-driven-samples-log-nsample"><a href="#config-solver-driven-samples-log-nsample"><code>"NSample"</code></a> <span class="config-type">integer</span> <span class="config-required">required</span> <span class="config-constraint"><code>≥ 1</code></span></dt>
  <dd>
    <p>Number of samples.</p>
  </dd>
  <dt id="config-solver-driven-samples-log-savestep"><a href="#config-solver-driven-samples-log-savestep"><code>"SaveStep"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>0</code></span></dt>
  <dd>
    <p>Save fields every N steps. <code>0</code> disables saving.</p>
  </dd>
  <dt id="config-solver-driven-samples-log-addtoprom"><a href="#config-solver-driven-samples-log-addtoprom"><code>"AddToPROM"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span></dt>
  <dd>
    <p>Force inclusion in the PROM for adaptive sweep.</p>
  </dd>
</dl>
```

### [Transient Solver](@id config-solver-transient)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-solver">Solver</a>/<a href="#config-solver-transient">Transient</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></span>
```

Configuration for the time domain driven solver. Only relevant when [`/Problem/Type`](#config-problem-type) is `"Transient"`. Simulations always start from rest at *t* = 0.

```@raw html
<dl class="palace-config">
  <dt id="config-solver-transient-excitation"><a href="#config-solver-transient-excitation"><code>"Excitation"</code></a> <span class="config-type">string</span> <span class="config-required">required</span></dt>
  <dd>
    <p>Controls the time dependence of the source excitation.</p>
    <dl class="config-enum">
      <dt><code>"Sinusoidal"</code></dt>
      <dd>Sinusoidal excitation at a user-specified frequency.</dd>
      <dt><code>"Gaussian"</code></dt>
      <dd>Gaussian pulse with a user-specified width (defines the bandwidth).</dd>
      <dt><code>"DifferentiatedGaussian"</code></dt>
      <dd>Differentiated Gaussian pulse with a user-specified width.</dd>
      <dt><code>"ModulatedGaussian"</code></dt>
      <dd>Modulated Gaussian pulse at a center frequency and width, with no DC component.</dd>
      <dt><code>"Ramp"</code></dt>
      <dd>Differentiable unit step function to model the ramp up to a DC signal.</dd>
      <dt><code>"SmoothStep"</code></dt>
      <dd>Smooth many-times differentiable unit step over a specified width.</dd>
    </dl>
  </dd>
  <dt id="config-solver-transient-maxtime"><a href="#config-solver-transient-maxtime"><code>"MaxTime"</code></a> <span class="config-type">number</span> <span class="config-required">required</span> <span class="config-constraint"><code>&gt; 0.0</code></span></dt>
  <dd>
    <p>End of simulation time interval, ns.</p>
  </dd>
  <dt id="config-solver-transient-timestep"><a href="#config-solver-transient-timestep"><code>"TimeStep"</code></a> <span class="config-type">number</span> <span class="config-required">required</span> <span class="config-constraint"><code>&gt; 0.0</code></span></dt>
  <dd>
    <p>Uniform time step size, ns.</p>
  </dd>
  <dt id="config-solver-transient-type"><a href="#config-solver-transient-type"><code>"Type"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;Default&quot;</code></span></dt>
  <dd>
    <p>Time integration scheme for the second-order ODE system.</p>
    <dl class="config-enum">
      <dt><code>"GeneralizedAlpha"</code></dt>
      <dd>Second-order implicit generalized-α method with $\rho_\infty = 1$. Unconditionally stable.</dd>
      <dt><code>"ARKODE"</code></dt>
      <dd>SUNDIALS ARKode implicit Runge-Kutta with adaptive time-stepping. Requires SUNDIALS support (see <a href="../../install/#Configuration-options">installation options</a>).</dd>
      <dt><code>"CVODE"</code></dt>
      <dd>SUNDIALS CVODE implicit multistep method with adaptive time-stepping. Requires SUNDIALS support (see <a href="../../install/#Configuration-options">installation options</a>).</dd>
      <dt><code>"RungeKutta"</code></dt>
      <dd>Two-stage singly diagonal implicit Runge-Kutta (SDIRK). Second-order, L-stable.</dd>
      <dt><code>"Default"</code></dt>
      <dd>Use the default <code>&quot;GeneralizedAlpha&quot;</code> scheme.</dd>
    </dl>
  </dd>
  <dt id="config-solver-transient-excitationfreq"><a href="#config-solver-transient-excitationfreq"><code>"ExcitationFreq"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.0</code></span></dt>
  <dd>
    <p>Center frequency for harmonic source excitations, GHz. Only relevant when <code>&quot;Excitation&quot;</code> is <code>&quot;Sinusoidal&quot;</code>, <code>&quot;Gaussian&quot;</code>, <code>&quot;DifferentiatedGaussian&quot;</code>, or <code>&quot;ModulatedGaussian&quot;</code>.</p>
  </dd>
  <dt id="config-solver-transient-excitationwidth"><a href="#config-solver-transient-excitationwidth"><code>"ExcitationWidth"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.0</code></span></dt>
  <dd>
    <p>Pulse width for Gaussian-type source excitations, ns. Only relevant when <code>&quot;Excitation&quot;</code> is <code>&quot;Gaussian&quot;</code>, <code>&quot;DifferentiatedGaussian&quot;</code>, <code>&quot;ModulatedGaussian&quot;</code>, or <code>&quot;SmoothStep&quot;</code>.</p>
  </dd>
  <dt id="config-solver-transient-savestep"><a href="#config-solver-transient-savestep"><code>"SaveStep"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>0</code></span></dt>
  <dd>
    <p>Controls how often, in number of time steps, to save computed fields to disk for <a href="../../guide/postprocessing/#Visualization">visualization with ParaView</a>. Files are saved in the <code>paraview/</code> (and/or <code>gridfunction/</code>) directory under <a href="#config-problem-output"><code>/Problem/Output</code></a>.</p>
  </dd>
  <dt id="config-solver-transient-order"><a href="#config-solver-transient-order"><code>"Order"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>2</code></span> <span class="config-constraint"><code>≥ 2, ≤ 5</code></span></dt>
  <dd>
    <p>Order of adaptive Runge-Kutta integrators or maximum multistep method order, must be within <code>[2, 5]</code>. Only relevant when <code>&quot;Type&quot;</code> is <code>&quot;ARKODE&quot;</code> or <code>&quot;CVODE&quot;</code>.</p>
  </dd>
  <dt id="config-solver-transient-reltol"><a href="#config-solver-transient-reltol"><code>"RelTol"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.0001</code></span> <span class="config-constraint"><code>&gt; 0.0</code></span></dt>
  <dd>
    <p>Relative tolerance for adaptive time-stepping. Only relevant when <code>&quot;Type&quot;</code> is <code>&quot;ARKODE&quot;</code> or <code>&quot;CVODE&quot;</code>.</p>
  </dd>
  <dt id="config-solver-transient-abstol"><a href="#config-solver-transient-abstol"><code>"AbsTol"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>1.0e-9</code></span> <span class="config-constraint"><code>&gt; 0.0</code></span></dt>
  <dd>
    <p>Absolute tolerance for adaptive time-stepping. Only relevant when <code>&quot;Type&quot;</code> is <code>&quot;ARKODE&quot;</code> or <code>&quot;CVODE&quot;</code>.</p>
  </dd>
</dl>
```

### [Electrostatic Solver](@id config-solver-electrostatic)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-solver">Solver</a>/<a href="#config-solver-electrostatic">Electrostatic</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></span>
```

Configuration for the electrostatic solver. Only relevant when [`/Problem/Type`](#config-problem-type) is `"Electrostatic"`.

```@raw html
<dl class="palace-config">
  <dt id="config-solver-electrostatic-save"><a href="#config-solver-electrostatic-save"><code>"Save"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>0</code></span> <span class="config-constraint"><code>≥ 0</code></span></dt>
  <dd>
    <p>Number of computed electric field solutions to save to disk for <a href="../../guide/postprocessing/#Visualization">visualization with ParaView</a>, ordered by the entries in the computed capacitance matrix. Files are saved in the <code>paraview/</code> (and/or <code>gridfunction/</code>) directory under <a href="#config-problem-output"><code>/Problem/Output</code></a>.</p>
  </dd>
</dl>
```

### [Magnetostatic Solver](@id config-solver-magnetostatic)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-solver">Solver</a>/<a href="#config-solver-magnetostatic">Magnetostatic</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></span>
```

Configuration for the magnetostatic solver. Only relevant when [`/Problem/Type`](#config-problem-type) is `"Magnetostatic"`.

```@raw html
<dl class="palace-config">
  <dt id="config-solver-magnetostatic-save"><a href="#config-solver-magnetostatic-save"><code>"Save"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>0</code></span> <span class="config-constraint"><code>≥ 0</code></span></dt>
  <dd>
    <p>Number of computed magnetic field solutions to save to disk for <a href="../../guide/postprocessing/#Visualization">visualization with ParaView</a>, ordered by the entries in the computed inductance matrix. Files are saved in the <code>paraview/</code> (and/or <code>gridfunction/</code>) directory under <a href="#config-problem-output"><code>/Problem/Output</code></a>.</p>
  </dd>
</dl>
```

### [Linear Solver](@id config-solver-linear)

```@raw html
<p class="config-keypath"><em>Path:</em> <code>/<a href="#config-solver">Solver</a>/<a href="#config-solver-linear">Linear</a></code></p>
```

```@raw html
<span class="config-section-badges"><span class="config-type">object</span> <span class="config-default">default: <code>—</code></span></span>
```

Configuration for the linear solver used by all simulation types.

```@raw html
<dl class="palace-config">
  <dt id="config-solver-linear-type"><a href="#config-solver-linear-type"><code>"Type"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;Default&quot;</code></span></dt>
  <dd>
    <p>Specifies the solver used for preconditioning the linear system.</p>
    <dl class="config-enum">
      <dt><code>"SuperLU"</code></dt>
      <dd><a href="https://github.com/xiaoyeli/superlu_dist">SuperLU_DIST</a> sparse direct solver in real double precision. For frequency domain problems uses a real approximation to the complex system matrix. Requires SuperLU_DIST support (see <a href="../../install/#Configuration-options">installation options</a>).</dd>
      <dt><code>"STRUMPACK"</code></dt>
      <dd><a href="https://portal.nersc.gov/project/sparse/strumpack">STRUMPACK</a> sparse direct solver in real double precision. Not compatible with magnetostatics (singular curl-curl operator); use <code>&quot;AMS&quot;</code> instead. Requires STRUMPACK support (see <a href="../../install/#Configuration-options">installation options</a>).</dd>
      <dt><code>"MUMPS"</code></dt>
      <dd><a href="http://mumps.enseeiht.fr/">MUMPS</a> sparse direct solver in real double precision. Requires MUMPS support (see <a href="../../install/#Configuration-options">installation options</a>).</dd>
      <dt><code>"AMS"</code></dt>
      <dd>Hypre's <a href="https://hypre.readthedocs.io/en/latest/solvers-ams.html">Auxiliary-space Maxwell Solver (AMS)</a>, an algebraic multigrid (AMG)-based preconditioner for curl-curl operators.</dd>
      <dt><code>"BoomerAMG"</code></dt>
      <dd>Hypre's <a href="https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html">BoomerAMG</a> algebraic multigrid solver.</dd>
      <dt><code>"Jacobi"</code></dt>
      <dd>Diagonal Jacobi preconditioner (not recommended in general).</dd>
      <dt><code>"Default"</code></dt>
      <dd>Use <code>&quot;AMS&quot;</code> for curl-curl and time-domain problems; a sparse direct solver (if available) for frequency domain; <code>&quot;BoomerAMG&quot;</code> for electrostatics.</dd>
    </dl>
  </dd>
  <dt id="config-solver-linear-ksptype"><a href="#config-solver-linear-ksptype"><code>"KSPType"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;Default&quot;</code></span></dt>
  <dd>
    <p>Specifies the iterative Krylov subspace solver.</p>
    <dl class="config-enum">
      <dt><code>"CG"</code></dt>
      <dd>Preconditioned conjugate gradient.</dd>
      <dt><code>"GMRES"</code></dt>
      <dd>GMRES.</dd>
      <dt><code>"FGMRES"</code></dt>
      <dd>Flexible GMRES.</dd>
      <dt><code>"Default"</code></dt>
      <dd>Use <code>&quot;GMRES&quot;</code> for frequency domain problems; <code>&quot;CG&quot;</code> for real symmetric positive-definite problems (transient, electrostatic, magnetostatic).</dd>
    </dl>
  </dd>
  <dt id="config-solver-linear-tol"><a href="#config-solver-linear-tol"><code>"Tol"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>1.0e-6</code></span> <span class="config-constraint"><code>≥ 0.0</code></span></dt>
  <dd>
    <p>Relative residual convergence tolerance for the iterative linear solver.</p>
  </dd>
  <dt id="config-solver-linear-maxits"><a href="#config-solver-linear-maxits"><code>"MaxIts"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>100</code></span> <span class="config-constraint"><code>&gt; 0</code></span></dt>
  <dd>
    <p>Maximum number of iterations for the iterative linear solver.</p>
  </dd>
  <dt id="config-solver-linear-maxsize"><a href="#config-solver-linear-maxsize"><code>"MaxSize"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>0</code></span> <span class="config-constraint"><code>&gt; 0</code></span></dt>
  <dd>
    <p>Maximum Krylov space size for GMRES/FGMRES. A value less than 1 defaults to <code>&quot;MaxIts&quot;</code>.</p>
  </dd>
  <dt id="config-solver-linear-initialguess"><a href="#config-solver-linear-initialguess"><code>"InitialGuess"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>true</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Use the previous solution as an initial guess for the iterative solver.</p>
  </dd>
  <dt id="config-solver-linear-mgmaxlevels"><a href="#config-solver-linear-mgmaxlevels"><code>"MGMaxLevels"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>100</code></span> <span class="config-constraint"><code>≥ 1</code></span></dt>
  <dd>
    <p>When greater than 1, enables <a href="https://en.wikipedia.org/wiki/Multigrid_method">geometric multigrid preconditioning</a> using p- and h-multigrid coarsening. The solver specified by <code>&quot;Type&quot;</code> is used on the coarsest level; fine levels use Chebyshev smoothing.</p>
  </dd>
  <dt id="config-solver-linear-mgcoarsentype"><a href="#config-solver-linear-mgcoarsentype"><code>"MGCoarsenType"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;Logarithmic&quot;</code></span></dt>
  <dd>
    <p>Coarsening strategy for constructing p-multigrid levels.</p>
    <dl class="config-enum">
      <dt><code>"Logarithmic"</code></dt>
      <dd>Logarithmic coarsening.</dd>
      <dt><code>"Linear"</code></dt>
      <dd>Linear coarsening.</dd>
    </dl>
  </dd>
  <dt id="config-solver-linear-mgusemesh"><a href="#config-solver-linear-mgusemesh"><code>"MGUseMesh"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>true</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Use the mesh hierarchy for h-multigrid coarsening.</p>
  </dd>
  <dt id="config-solver-linear-mgauxiliarysmoother"><a href="#config-solver-linear-mgauxiliarysmoother"><code>"MGAuxiliarySmoother"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>true</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Use an auxiliary space smoother on multigrid levels.</p>
  </dd>
  <dt id="config-solver-linear-mgcycleits"><a href="#config-solver-linear-mgcycleits"><code>"MGCycleIts"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>0</code></span> <span class="config-constraint"><code>&gt; 0</code></span></dt>
  <dd>
    <p>Number of V-cycle iterations per preconditioner application for multigrid preconditioners (when geometric multigrid is enabled or <code>&quot;Type&quot;</code> is <code>&quot;AMS&quot;</code> or <code>&quot;BoomerAMG&quot;</code>). A value less than 1 defaults to 2 for AMS frequency domain problems, 1 otherwise.</p>
  </dd>
  <dt id="config-solver-linear-mgsmoothits"><a href="#config-solver-linear-mgsmoothits"><code>"MGSmoothIts"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>1</code></span> <span class="config-constraint"><code>&gt; 0</code></span></dt>
  <dd>
    <p>Number of pre- and post-smooth iterations for multigrid preconditioners, when the geometric multigrid preconditioner is enabled.</p>
  </dd>
  <dt id="config-solver-linear-mgsmoothorder"><a href="#config-solver-linear-mgsmoothorder"><code>"MGSmoothOrder"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>0</code></span> <span class="config-constraint"><code>&gt; 0</code></span></dt>
  <dd>
    <p>Polynomial smoothing order for geometric multigrid. A value less than 1 defaults to <code>max(2 × solution order, 4)</code>.</p>
  </dd>
  <dt id="config-solver-linear-mgsmootheigscalemax"><a href="#config-solver-linear-mgsmootheigscalemax"><code>"MGSmoothEigScaleMax"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>1.0</code></span> <span class="config-constraint"><code>&gt; 0</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Maximum eigenvalue scale for the Chebyshev smoother.</p>
  </dd>
  <dt id="config-solver-linear-mgsmootheigscalemin"><a href="#config-solver-linear-mgsmootheigscalemin"><code>"MGSmoothEigScaleMin"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.0</code></span> <span class="config-constraint"><code>≥ 0</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Minimum eigenvalue scale for the Chebyshev smoother.</p>
  </dd>
  <dt id="config-solver-linear-mgsmoothchebyshev4th"><a href="#config-solver-linear-mgsmoothchebyshev4th"><code>"MGSmoothChebyshev4th"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>true</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Use 4th-kind Chebyshev polynomial smoother.</p>
  </dd>
  <dt id="config-solver-linear-pcmatreal"><a href="#config-solver-linear-pcmatreal"><code>"PCMatReal"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span></dt>
  <dd>
    <p>When <code>true</code>, builds the preconditioner for frequency domain problems using a real-valued approximation of the system matrix. Always used on the coarsest multigrid level regardless of this setting.</p>
  </dd>
  <dt id="config-solver-linear-pcmatshifted"><a href="#config-solver-linear-pcmatshifted"><code>"PCMatShifted"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span></dt>
  <dd>
    <p>When <code>true</code>, builds the preconditioner using a positive-definite approximation by flipping the sign of the mass matrix contribution. Can improve performance at high frequencies relative to the lowest nonzero eigenfrequencies of the model.</p>
  </dd>
  <dt id="config-solver-linear-complexcoarsesolve"><a href="#config-solver-linear-complexcoarsesolve"><code>"ComplexCoarseSolve"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span></dt>
  <dd>
    <p>When <code>true</code>, uses the true complex-valued system matrix for the coarse-level solver. When <code>false</code>, uses the real-valued approximation.</p>
  </dd>
  <dt id="config-solver-linear-dropsmallentries"><a href="#config-solver-linear-dropsmallentries"><code>"DropSmallEntries"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span></dt>
  <dd>
    <p>When <code>true</code>, drops entries smaller than double-precision machine epsilon from the system matrix used in the sparse direct solver.</p>
  </dd>
  <dt id="config-solver-linear-reorderingreuse"><a href="#config-solver-linear-reorderingreuse"><code>"ReorderingReuse"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>true</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Reuse the matrix reordering from the previous solve.</p>
  </dd>
  <dt id="config-solver-linear-pcside"><a href="#config-solver-linear-pcside"><code>"PCSide"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;Default&quot;</code></span></dt>
  <dd>
    <p>Preconditioning side. Not all options are available for all iterative solver choices.</p>
    <dl class="config-enum">
      <dt><code>"Left"</code></dt>
      <dd>Left preconditioning.</dd>
      <dt><code>"Right"</code></dt>
      <dd>Right preconditioning.</dd>
      <dt><code>"Default"</code></dt>
      <dd>Solver-default side.</dd>
    </dl>
  </dd>
  <dt id="config-solver-linear-columnordering"><a href="#config-solver-linear-columnordering"><code>"ColumnOrdering"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;Default&quot;</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Column ordering algorithm for sparse direct solvers: <code>&quot;METIS&quot;</code>, <code>&quot;ParMETIS&quot;</code>, <code>&quot;Scotch&quot;</code>, <code>&quot;PTScotch&quot;</code>, <code>&quot;PORD&quot;</code>, <code>&quot;AMD&quot;</code>, <code>&quot;RCM&quot;</code>, <code>&quot;Default&quot;</code>.</p>
  </dd>
  <dt id="config-solver-linear-strumpackcompressiontype"><a href="#config-solver-linear-strumpackcompressiontype"><code>"STRUMPACKCompressionType"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;None&quot;</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>STRUMPACK compression type: <code>&quot;None&quot;</code>, <code>&quot;BLR&quot;</code>, <code>&quot;HSS&quot;</code>, <code>&quot;HODLR&quot;</code>, <code>&quot;ZFP&quot;</code>, <code>&quot;BLR-HODLR&quot;</code>, <code>&quot;ZFP-BLR-HODLR&quot;</code>.</p>
  </dd>
  <dt id="config-solver-linear-strumpackcompressiontol"><a href="#config-solver-linear-strumpackcompressiontol"><code>"STRUMPACKCompressionTol"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>0.001</code></span> <span class="config-constraint"><code>≥ 0.0</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Tolerance for STRUMPACK lossy compression.</p>
  </dd>
  <dt id="config-solver-linear-strumpacklossyprecision"><a href="#config-solver-linear-strumpacklossyprecision"><code>"STRUMPACKLossyPrecision"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>16</code></span> <span class="config-constraint"><code>≥ 0</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Precision bits for STRUMPACK ZFP lossy compression.</p>
  </dd>
  <dt id="config-solver-linear-strumpackbutterflylevels"><a href="#config-solver-linear-strumpackbutterflylevels"><code>"STRUMPACKButterflyLevels"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>0</code></span> <span class="config-constraint"><code>≥ 0</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Number of butterfly levels for STRUMPACK HODLR compression.</p>
  </dd>
  <dt id="config-solver-linear-superlu3dcommunicator"><a href="#config-solver-linear-superlu3dcommunicator"><code>"SuperLU3DCommunicator"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Use a 3D process grid communicator for SuperLU_DIST.</p>
  </dd>
  <dt id="config-solver-linear-amsvectorinterpolation"><a href="#config-solver-linear-amsvectorinterpolation"><code>"AMSVectorInterpolation"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Use vector interpolation in AMS.</p>
  </dd>
  <dt id="config-solver-linear-amssingularoperator"><a href="#config-solver-linear-amssingularoperator"><code>"AMSSingularOperator"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Indicate to AMS that the operator is singular (e.g. magnetostatics with no essential boundary conditions).</p>
  </dd>
  <dt id="config-solver-linear-amgaggressivecoarsening"><a href="#config-solver-linear-amgaggressivecoarsening"><code>"AMGAggressiveCoarsening"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span> <span class="config-advanced">advanced</span></dt>
  <dd>
    <p>Enable aggressive coarsening in BoomerAMG.</p>
  </dd>
  <dt id="config-solver-linear-amsmaxits"><a href="#config-solver-linear-amsmaxits"><code>"AMSMaxIts"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>0</code></span> <span class="config-constraint"><code>≥ 0</code></span></dt>
  <dd>
    <p>Number of AMS iterations per preconditioner application when geometric multigrid is enabled (<code>&quot;MGMaxLevels&quot;</code> &gt; 1). A value less than 1 defaults to the solution order (<a href="#config-solver-order"><code>/Solver/Order</code></a>).</p>
  </dd>
  <dt id="config-solver-linear-divfreetol"><a href="#config-solver-linear-divfreetol"><code>"DivFreeTol"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>1.0e-12</code></span> <span class="config-constraint"><code>≥ 0.0</code></span></dt>
  <dd>
    <p>Relative tolerance for divergence-free cleaning used in the eigenmode simulation type. Ignored when a nonzero Floquet wave vector is specified or when a nonzero <a href="#config-domains-materials-londondepth">LondonDepth</a> is used.</p>
  </dd>
  <dt id="config-solver-linear-divfreemaxits"><a href="#config-solver-linear-divfreemaxits"><code>"DivFreeMaxIts"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>1000</code></span> <span class="config-constraint"><code>≥ 0</code></span></dt>
  <dd>
    <p>Maximum number of iterations for divergence-free cleaning. Ignored under the same conditions as <code>&quot;DivFreeTol&quot;</code>.</p>
  </dd>
  <dt id="config-solver-linear-estimatortol"><a href="#config-solver-linear-estimatortol"><code>"EstimatorTol"</code></a> <span class="config-type">number</span> <span class="config-default">default: <code>1.0e-6</code></span> <span class="config-constraint"><code>≥ 0.0</code></span></dt>
  <dd>
    <p>Relative tolerance for the flux projection solve used in the error estimate calculation.</p>
  </dd>
  <dt id="config-solver-linear-estimatormaxits"><a href="#config-solver-linear-estimatormaxits"><code>"EstimatorMaxIts"</code></a> <span class="config-type">integer</span> <span class="config-default">default: <code>10000</code></span> <span class="config-constraint"><code>≥ 0</code></span></dt>
  <dd>
    <p>Maximum number of iterations for the flux projection solve used in the error estimate calculation.</p>
  </dd>
  <dt id="config-solver-linear-estimatormg"><a href="#config-solver-linear-estimatormg"><code>"EstimatorMG"</code></a> <span class="config-type">boolean</span> <span class="config-default">default: <code>false</code></span></dt>
  <dd>
    <p>When <code>true</code>, uses a multigrid preconditioner with AMG coarse solve for the error estimator linear solver instead of Jacobi.</p>
  </dd>
  <dt id="config-solver-linear-gsorthogonalization"><a href="#config-solver-linear-gsorthogonalization"><code>"GSOrthogonalization"</code></a> <span class="config-type">string</span> <span class="config-default">default: <code>&quot;MGS&quot;</code></span></dt>
  <dd>
    <p>Gram-Schmidt variant used to orthogonalize vectors in Krylov subspace methods and other parts of the solver.</p>
    <dl class="config-enum">
      <dt><code>"MGS"</code></dt>
      <dd>Modified Gram-Schmidt.</dd>
      <dt><code>"CGS"</code></dt>
      <dd>Classical Gram-Schmidt.</dd>
      <dt><code>"CGS2"</code></dt>
      <dd>Two-step classical Gram-Schmidt with reorthogonalization.</dd>
    </dl>
  </dd>
</dl>
```
