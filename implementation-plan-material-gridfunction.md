# Implementation plan: spatially-varying material properties (GridFunction coefficients) for the 2D mode solver

Tracking issue: [awslabs/palace#925](https://github.com/awslabs/palace/issues/925)
Companion draft: [`issue-draft-material-gridfunction.md`](issue-draft-material-gridfunction.md)

## 1. Goal and non-goals

### Goal

Allow the `"BoundaryMode"` solver to take relative permittivity `ε_r(x)` and/or
permeability `μ_r(x)` as spatially-varying fields, supplied either

1. programmatically as an MFEM `GridFunction` on the solve mesh (used by tests), or
2. from a standard field file (Gmsh `$NodeData` / VTK-VTU point or cell data) whose
   named arrays are read into a `GridFunction` on the solve mesh.

The field is sampled once to the quadrature points already used by Palace's libCEED
operators and consumed by Q-functions as an extra inactive (`CEED_EVAL_NONE`) input,
mirroring the existing geometry-data mechanism. The default path (no field files) must
remain bit-identical.

### Non-goals (initial PR)

- General full-tensor / anisotropic fields. The first version is **scalar (isotropic)**
  fields, which avoids the 8–9× memory blow-up raised in the issue discussion. Full
  tensors are a follow-up.
- Fields for the 3D solvers. The machinery is designed to generalize, but only the 2D
  mode solver is wired up.
- Fields defined on a background mesh requiring cross-mesh projection.
- Fields combined with the 3D→2D submesh extraction path (`Solver.BoundaryMode.Attributes`)
  in the first cut (see §8). Reading the field on the parent 3D mesh and projecting it to
  the extracted cross-section is a Phase 2 item.
- Field-aware error estimation / domain-energy postprocessing (they continue to use the
  constant material operator; document this). Mode eigenvalues and fields are correct.

## 2. Key facts established from the code

These findings drive the design and should be re-verified while implementing.

1. **Where the mode operators are built.** `BoundaryModeOperator` (`palace/models/boundarymodeoperator.cpp`)
   assembles the fixed blocks `Atn`/`Btt` through `mode_assembly::AssembleAtn`/`AssembleBtt`
   and builds a `mode_assembly::ModeOperatorModel`, which assembles every frequency-independent
   component matrix. The p-multigrid preconditioner lives in
   `ModeEigenSolver::AssembleAttPreconditioner`/`AssembleAnnPreconditioner`
   (`palace/models/modeeigensolver.cpp`, ~lines 715–870).

2. **Integrators / Q-functions used (2D).** From `modeoperatorassembly.cpp` and
   `modeeigensolver.cpp`:

   | Use | Integrator | 2D Q-function (`f_apply_*`) |
   |---|---|---|
   | `μ⁻¹` in-plane ND mass (`Btt`) | `VectorFEMassIntegrator` | `hcurlmass_22` |
   | `μ⁻¹` mixed gradient (`Atn`) | `MixedVectorGradientIntegrator` (H_CURL) | `hcurlh1d_22` |
   | scalar `μ⁻¹` out-of-plane curl-curl | `CurlCurlIntegrator` | `l2_1` |
   | `-μ⁻¹` H1 diffusion (`Ann`) | `DiffusionIntegrator` | `hcurl_22` |
   | `ε` H1 mass (`Ann`) | `MassIntegrator` (`num_comp = 1`) | `h1_1` |
   | Preconditioner curl-curl + mass | `CurlCurlMassIntegrator` | `hdivmass_22` |
   | Preconditioner auxiliary diffusion + mass | `DiffusionMassIntegrator` | `hcurlmass_22` |

   The exact operator additionally has boundary mass integrators (`VectorFEMassIntegrator`,
   `MassIntegrator`) fed by surface-impedance/conductivity/rational-impedance coefficients;
   those stay constant.

3. **Geometry quadrature data is cached at mesh level.** `Mesh::GetCeedGeomFactorData(ceed)`
   lazily builds `ceed::CeedGeomFactorData` (attribute, `w|J|`, `adj(J)ᵀ/|J|`) per element
   geometry (`palace/fem/mesh.cpp:140–326`). `BilinearForm::PartialAssemble`
   (`palace/fem/bilinearform.cpp:27–107`) iterates this map and calls
   `integ->Assemble(ceed, trial_restr, test_restr, trial_basis, test_basis, data.geom_data,
   data.geom_data_restr, &sub_op)`.

4. **The quadrature rule is fixed for a run.** `FiniteElementSpace::BuildCeedBasis`
   (`palace/fem/fespace.cpp:134–160`) picks `IntRules.Get(geom, fem::DefaultIntegrationOrder::Get(T))`
   with `T` set to the **mesh nodal FE**, and `p_trial` = `solver.order`. So every p-level
   in a `FiniteElementSpaceHierarchy` uses the *same* quadrature rule, and
   `CeedOperatorMultigridLevelCreate` (`palace/fem/libceed/operator.cpp:1050–1110`)
   carries non-active fields (like `geom_data`, and our new field) to coarsened operators
   unchanged. **Per-quadrature field data is therefore level-independent and can be cached
   once per (ceed, geometry).**

5. **Two coefficient consumption paths exist.** The "apply" Q-functions read the
   attribute→material coefficient from `ctx` (`f_apply_*`), while the "build" Q-functions
   fold it into `q_data` once (`f_build_*`, driven by `assemble_q_data`). The mode solver
   currently uses the apply path.

6. **Field data flow.** `driver.cpp:66–78`: `mesh::Load` (serial) → `solver->Preprocess`
   (BoundaryMode submesh extraction) → `mesh::Partition` → `mesh::RefineMesh` (parallel
   uniform). `BoundaryModeSolver::Solve(mesh)` receives the **final** distributed solve
   mesh.

## 3. Recommended design

### 3.1 Core representation: `MaterialField`

Introduce a `palace::MaterialField` object (new files
`palace/fem/materialfield.hpp` / `.cpp`) that owns:

- a real scalar `mfem::ParGridFunction` on the solve mesh (P1, or higher order when built
  directly), and
- a lazily built cache of per-(Ceed, geometry) quadrature data:

```cpp
namespace palace {

// Per-geometry quadrature-point samples of a scalar material property, laid out like the
// geometry factor data so it can be attached to a libCEED operator as an inactive input.
struct CeedMaterialFieldData
{
  int ncomp;                        // 1 for scalar fields
  CeedVector data;                  // [ncomp, Q] per element, strided restriction
  CeedElemRestriction data_restr;
};

class MaterialField
{
public:
  // scaling multiplies the sampled values (e.g. 1.0 for ε_r, and we sample 1/μ_r for μ).
  MaterialField(const mfem::ParGridFunction &gf, double scaling = 1.0);

  // Build/return quadrature data for this (ceed, geometry) over `indices` (element order
  // must match CeedGeomFactorData::indices). Thread-safe-cached per ceed.
  const CeedMaterialFieldData &GetData(Ceed ceed, mfem::Geometry::Type geom,
                                       const std::vector<int> &indices) const;

  double GetMax() const;  // max sampled value, for the auto shift target

private:
  const mfem::ParGridFunction &gf;
  double scaling;
  mutable ceed::GeometryObjectMap<CeedMaterialFieldData> data;  // per ceed
};

}  // namespace palace
```

**Sampling.** For each element `e` in `indices` and each reference quadrature point `ip`
from `mfem::IntRules.Get(geom, fem::DefaultIntegrationOrder::Get(mesh, geom))` (the same
rule `BuildCeedBasis` uses), evaluate `gf.GetValue(T_e, ip) * scaling` and store it
element-contiguously, `[e*Q + q]`. Build a strided `CeedElemRestrictionCreateStrided(ceed,
num_elem, Q, ncomp, ...)`, exactly like `AssembleGeometryData` builds `geom_data_restr`.

Because the values are computed once per run, evaluate every `(ceed, geometry)` **eagerly**
before the OpenMP region in `BilinearForm::PartialAssemble` to avoid races on the cache
(generate them for all `ceed::internal::GetCeedObjects()` from the mesh geometry map).

For a permeability field the integrators need `μ⁻¹`; sample the field and store
`μ⁻¹(x)` (or store `μ_r` and let the Q-function invert — prefer precomputing `μ⁻¹` so the
Q-function stays a pure multiply, matching the existing coefficient semantics).

### 3.2 Integrator plumbing

Extend the integrator base class (`palace/fem/integrator.hpp`) with an optional field:

```cpp
class BilinearFormIntegrator
{
protected:
  const MaterialField *field = nullptr;
public:
  void SetMaterialField(const MaterialField *f) { field = f; }
  const MaterialField *GetMaterialField() const { return field; }
  ...
  virtual void Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                        CeedElemRestriction test_restr, CeedBasis trial_basis,
                        CeedBasis test_basis, CeedVector geom_data,
                        CeedElemRestriction geom_data_restr,
                        const CeedMaterialFieldData *field_data,  // new, may be null
                        CeedOperator *op) const = 0;
};
```

`BilinearForm::PartialAssemble` (domain and boundary branches) computes, for each
integrator, `field_data` from `integ->GetMaterialField()` (via
`GetData(ceed, geom, data.indices)`) and passes it. All existing integrator overrides add
the parameter; only the seven from §2.2 actually consume it.

`AssembleCeedOperator` (`palace/fem/libceed/integrator.cpp`) gains support for attaching
the field as an extra inactive input:

- If `field_data`, `CeedQFunctionAddInput(apply_qf, "coeff_data", ncomp, CEED_EVAL_NONE)`.
- `CeedOperatorSetField(*op, "coeff_data", field_data->data_restr, CEED_BASIS_NONE,
  field_data->data)`.

The **default path is untouched**: when no field is set, no input is added and the same
Q-function is selected.

### 3.3 Q-function variants

Add field-aware variants for the affected 2D Q-functions. Each variant is a small delta
on the existing one: read `coeff_data[i]` (pre-inverted for μ) and use it in place of the
attribute coefficient. Suggested file layout, mirroring the existing ones:

- `qfunctions/22/hcurlmass_field_22_qf.h` → ND mass and `DiffusionMass`
- `qfunctions/22/hdivmass_field_22_qf.h` → `CurlCurlMass`
- `qfunctions/22/hcurl_field_22_qf.h` → H1 diffusion
- `qfunctions/22/hcurlh1d_field_22_qf.h` → mixed gradient
- `qfunctions/22/l2_field_1_qf.h` → 2D curl-curl
- `qfunctions/22/h1_field_1_qf.h` → H1 mass

Selection happens in each integrator's `Assemble`:

```cpp
if (field_data) { info.apply_qf = f_apply_hcurlmass_field_22; ... }
else if (assemble_q_data) { ... existing build ... }
else { ... existing apply ... }
```

**Handling the preconditioner's combined coefficients.** `CurlCurlMassIntegrator` and
`DiffusionMassIntegrator` carry a coefficient that mixes `ε` and `μ⁻¹` (plus a frequency
scalar). Give those two field variants **two** inactive inputs, `"eps_field"` and
`"muinv_field"`, and have the Q-function apply each field to its own term (the frequency
scalars stay in `ctx`). A small `uint32_t` flag in the Q-function context records which
fields are present so a single variant can serve ε-only, μ-only, or both. If this proves
fiddly, fall back (documented) to keeping the **preconditioner** on constant coefficients
— it does not affect the converged eigenvalue, only the iteration count — and still ship
the field-aware exact operator.

The boundary integrators (surface impedance/conductivity/rational impedance) are left
field-unaware.

### 3.4 Assembly call sites

- `mode_assembly::AssembleAtn` / `AssembleBtt`: accept the `μ` field and call
  `SetMaterialField` on `MixedVectorGradientIntegrator` / `VectorFEMassIntegrator`.
- `ModeOperatorModel` constructor and its `muinv_cc`, `neg_muinv`, `eps_n`, `eps`-mass
  integrators: attach the `μ` / `ε` fields.
- `ModeEigenSolver::AssembleAttPreconditioner` / `AssembleAnnPreconditioner`: attach the
  fields to `CurlCurlMassIntegrator` / `DiffusionMassIntegrator` (if in scope).

### 3.5 Coefficient semantics

Initial semantics: **global replacement over the whole solve domain**. When a `μ` field
is present, every domain integrator that consumes `μ⁻¹` uses the sampled field and ignores
the per-material constant for that property; likewise for `ε`. All other properties
(London depth, conductivity, loss tangent) and all boundary coefficients keep using
`Domains.Materials`. `Domains.Materials` must still be present and define the remaining
properties.

Rationale: this exactly covers the target use case (a TCAD-derived profile over the whole
cross-section) and is the only well-defined semantics for a single global field without an
attribute mask. Per-region fields and field×constant weighting are deferred.

### 3.6 Field-file reader

New `palace/utils/gridfunctionio.hpp` / `.cpp`:

```cpp
namespace palace::fieldio {

// Parse a named scalar/vector array from a Gmsh (.msh v2.2 $NodeData) or VTK/VTU file
// into (points, values). Cell data is returned as element-centered values.
struct FieldData
{
  int dim;                                    // 2 or 3
  std::vector<std::array<double, 3>> points;  // node or cell-centroid coords
  std::vector<double> values;                 // scalar (vdim = 1) initially
  int vdim = 1;
};

FieldData ReadFieldFile(const std::string &filename, const std::string &name);

}  // namespace palace::fieldio
```

- **Gmsh v2.2**: parse `$PhysicalNames`/`$Nodes` and `$NodeData` blocks; `$NodeData`
  carries the string tag, real/integer tags, and per-node `(tag, value...)`. Match the
  requested `name` against the string tag. (MFEM's reader drops `$NodeData`; we only need
  the data section, and we can reuse MFEM to read the mesh itself.)
- **VTU**: parse `<PointData>` (and `<CellData>`) `<DataArray Name="..." format="ascii">`.
  Support ASCII first; base64/binary and appended data can be added if needed. Reuse
  MFEM's mesh reader for the topology.
- **Building the `GridFunction` on the solve mesh (MPI-safe).** Build a P1 H1 (nodal) or
  L2-order-0 (cell) `mfem::ParGridFunction` on `mesh.back()`. For each local vertex, match
  its coordinate to a file point with a spatial hash + tolerance (exact coordinate equality
  after nondimensionalization should hold when the file was written from the same mesh).
  Each rank does this independently — no communication is required and shared vertices get
  identical values. Abort with a clear message if any local vertex/cell is unmatched, or if
  the file's element count/topology is incompatible.
- **Refinement**: require the field file to describe the **final** solve mesh. If
  `Refinement.UniformLevels > 0` (or box/sphere refinement) is configured together with
  field files, abort with guidance pointing users to generate the field on the refined
  mesh. (Phase 2: read onto the input mesh and propagate through refinement with
  `mfem::InterpolationGridTransfer`.)
- Reject field files combined with the 3D→2D `Attributes` path in Phase 1 (clear error;
  Phase 2 handles the parent-mesh projection through the `mesh::Project3Dto2D` frame).

### 3.7 Driver wiring

In `BoundaryModeSolver::Solve` (before constructing `MaterialOperator`/`BoundaryModeOperator`):

1. If `iodata.solver.boundary_mode.material_fields` is set, read each file into a
   `GridFunction` (via the reader) or construct from the direct API.
2. Build the `ε` field (scaled values) and `μ⁻¹` field.
3. Pass them into `BoundaryModeOperator` (new constructor parameter, e.g.
   `const MaterialField *eps_field, const MaterialField *muinv_field`), which forwards
   them to `AssembleAtn`/`AssembleBtt`/`ModeOperatorModel`, and store them for
   `ModeEigenSolver`'s preconditioner.

A direct, file-free constructor (taking `GridFunction`s) is the seam used by unit tests.

### 3.8 Auto shift target

`BoundaryModeSolver::Solve` computes `kn_target = omega * sqrt(1.1 * GetMaxMuEpsilon())`
when `Target == 0`. Extend this to use the fields: compute the maximum of `ε_r(x)·μ_r(x)`
over the sampled field values where present, combined with the material constants where a
property has no field. Add `MaterialField::GetMax()` (max sampled value) and a small helper
that computes the effective `max μ_r · max ε_r`. Document the change.

## 4. Config schema and documentation

Add under `Solver.BoundaryMode` (localized to the mode solver for the initial PR):

```json
"Solver": {
  "BoundaryMode": {
    "Freq": 100.0,
    "MaterialFields": {
      "Permittivity": { "File": "eps.msh", "Name": "epsilon_r" },
      "Permeability": { "File": "mu.vtu",  "Name": "mu_r" }
    }
  }
}
```

- `palace/utils/configfile.hpp`: add a `MaterialFieldData` struct (file path + array name,
  with a default array name per property) and a
  `BoundaryModeSolverData::material_fields` member (e.g. `std::optional`-like struct with
  empty-means-unset, matching the codebase's convention).
- `palace/utils/configfile.cpp`: parse in `BoundaryModeSolverData::BoundaryModeSolverData`;
  validate path existence and (optionally) that Permittivity/Permeability are not both
  files without a `Domains.Materials` entry.
- `scripts/schema/config-schema.json`: add the `MaterialFields` object under the
  `BoundaryMode` definition (line ~2804) with `additionalProperties: false`, plus a
  `$defs/MaterialField` for `{File, Name}`. Bump SchemaVer if the project convention
  requires it (see the schema-versioning developer notes) and update the compatibility
  manifest if needed.
- `docs/src/guide/problem.md` (§"Mode analysis problems", ~line 238): document the new
  behavior, semantics, units (relative/dimensionless), and limitations.
- `CHANGELOG.md` "In progress → New Features": add an entry with the PR link placeholder.

## 5. Test plan

### 5.1 Regression / bit-identical default path

- Existing boundary-mode unit tests and the `cpw2d`/`cpw` regression solves must pass
  unchanged. No default-path code is exercised differently.
- Add a unit test asserting that a **constant** field equal to the material constant
  reproduces the default `kn` to machine precision (`WithinRel` ~1e-12), proving the
  field Q-function variant reduces correctly.

### 5.2 Mechanism unit tests

- `test-materialfield.cpp`: build a P1/P2 scalar field on a small `ParMesh`, verify
  `MaterialField::GetData` reproduces `gf.GetValue(T, ip)` at the quadrature points for
  quad and triangle meshes, serially and on 2 ranks.
- Field reader tests writing a Gmsh `$NodeData` and a VTU fixture in the shared temp dir
  and reading them back; verify values/topology mapping; verify clear failures on a name
  mismatch, topology mismatch, and missing vertices.
- Integrator equivalence: assemble e.g. the ND mass with a constant field via the new
  Q-function and compare the full matrix to the legacy attribute-based assembly
  (`CHECK` relative norm ~1e-13).

### 5.3 Convergence test (parabolic-index slab waveguide)

In `test-boundarymodeoperator.cpp`, following the existing `SolveRectangularModes` /
`SolveLondonSlabMode` structure (direct construction, not the config file reader):

- Geometry: `x ∈ [-W, W]`, `y ∈ [0, H]`, PEC at `y = 0, H` (mode uniform in `y`, `E`
  normal to the walls) and PEC/absorbing at `x = ±W` chosen wide enough that the Gaussian
  tail is below machine precision.
- `ε_r(x) = ε_c − (x/x₀)²`, represented exactly by an order-2 H1 `GridFunction`
  (`ProjectCoefficient`).
- Analytic full-vector TE solution: `k_n² = ω² ε_c − (2m+1)/x₀²` (Hermite–Gaussian
  profiles), for the lowest modes.
- Check `kn` against the analytic values and check the observed convergence rate under
  mesh/order refinement (`k_n − k_n^exact` at least super-linear; ideally close to the FE
  order). Because the material representation error is eliminated by the exact P2
  representation, the test isolates discretization error.
- Independent cross-check: reproduce the same profile with the per-element-attribute
  workaround on the same mesh and confirm the two agree to the workaround's
  piecewise-constant accuracy.
- MPI: run the same case on 2 ranks and check consistency.

### 5.4 Preconditioner test (if in scope)

Extend `"ModeEigenSolver p-multigrid preconditioning"` to a field case and check that the
multigrid result matches the direct result to the existing tolerance (and that the
iteration count is finite). If the preconditioner stays on constants, assert only
correctness, not iteration count.

## 6. Work breakdown (ordered)

> **Status:** Phases 1–2 complete. Phases 3–8 pending.

1. **`MaterialField` + integrator plumbing** (Q-functions, `AssembleCeedOperator`,
   `BilinearFormIntegrator`, `PartialAssemble`). Land the constant-field equivalence test
   here. **DONE** — `palace/fem/materialfield.{hpp,cpp}`, field-aware
   `f_apply_h1_field_1`, `AssembleCeedOperator` field input, per-integrator `field_data`
   parameter, eager cache warming in `PartialAssemble`; tests in
   `test/unit/test-materialfield.cpp` (elementwise-constant exact equivalence and smooth
   P2 sampling, serial and 2-rank) pass, and `[libceed]`/`[coefficient]` regressions pass.
2. **Exact mode operator wiring**: `AssembleAtn`/`AssembleBtt`/`ModeOperatorModel` and the
   direct-API seam. **DONE** — field-aware `f_apply_hcurl_field_22` (shared by ND mass, H1
   diffusion, and the mixed gradient) and `f_apply_l2_field_1` (2D curl-curl) wired into
   the four integrators; `AssembleAtn`/`AssembleBtt`/`ModeOperatorModel` replace the
   material constant with a sign-only unit coefficient and attach the field;
   `BoundaryModeOperator` accepts/exposes the ε and μ⁻¹ fields and `ModeEigenSolver`
   forwards them. A mode-integrator equivalence test and an end-to-end
   `BoundaryModeOperator` constant-field equivalence test pass.
3. **Convergence unit test** (P2 parabola, direct API) — validates the core feature.
   **DONE** — `BoundaryModeOperator parabolic-index slab convergence` solves
   `eps_r(x) = eps_c - ((x-W)/x0)^2` (exact P2 H1 field) against the analytic TE result
   `kn = sqrt(omega^2 eps_c - omega/x0)`. Errors are 2.2e-4, 1.5e-5, 9.3e-7 for
   nx = 16, 32, 64 (ratio ≈ 16, i.e. O(h^4) for order 2).
2. **Exact mode operator wiring**: `AssembleAtn`/`AssembleBtt`/`ModeOperatorModel` and the
   direct-API seam.
3. **Convergence unit test** (P2 parabola, direct API) — validates the core feature.
4. **Config + driver wiring + auto target**. **DONE** — `Solver.BoundaryMode.MaterialFields`
   (`Permittivity`/`Permeability` `{File, Name}`) parsed into `BoundaryModeSolverData`;
   `BoundaryModeSolver::Solve` reads the files, scales to the nondimensional mesh, inverts
   the permeability, builds `MaterialField`s, and passes them to `BoundaryModeOperator`;
   the auto shift target uses the field extrema. Schema updated (SchemaVer 2-3-0) with docs
   and a CHANGELOG entry.
5. **Field reader** (`gridfunctionio.*`) + reader tests. **DONE** — `ReadFieldFile` parses
   Gmsh v2.2 `$NodeData` (with `$Nodes`) and VTU ASCII point data; `BuildGridFunction` maps
   the points to a P1 grid function on the solve mesh by coordinate matching. Tests read
   both formats back and validate config parsing.
6. **Preconditioner** field-awareness (or documented fallback).
7. **Docs / schema / CHANGELOG.**
8. **Phase 2 follow-ups**: refinement propagation, 3D→2D parent projection, anisotropic
   and vector fields, generalizing to other solvers.

## 7. Risks and mitigations

- **Q-function proliferation.** Mitigate by scoping to 2D (`22`) variants and generating
  the field variants as minimal copies of the existing ones via the existing helpers
  (`coeff_*_qf.h`, `utils_22_qf.h`). Consider a helper macro for the field lookup and a
  build-gating CI check that the JIT sources list stays complete.
- **Thread-safety of the field cache.** Build all `(ceed, geometry)` data eagerly before
  the OpenMP region in `PartialAssemble`; this matches how geometry data is partitioned by
  thread and avoids concurrent map mutation.
- **Preconditioner mixing ε and μ.** Fallback: keep the preconditioner constant; it
  affects only convergence speed. Explicitly document.
- **Field/refinement and 3D→2D mismatches.** Fail loudly with actionable messages; keep
  the reader optional and orthogonal to the core mechanism so the core can land first.
- **GPU/device paths.** The new field is an inactive host/device `CeedVector` like
  `geom_data`; verify `.UseDevice` behavior on CUDA/HIP. Sampling into the vector should
  use `HostWrite` then copy to device, mirroring `InitCeedVector`.
- **Semantic surprises.** Document precisely that fields replace, not weight, the config
  material for that property, and that postprocessing/error estimation still use constants.

## 8. Open questions to confirm before coding

1. **Config location**: `Solver.BoundaryMode.MaterialFields` (localized) vs. a
   domain-level `Model`/`Domains` field (more general, larger surface). Recommendation:
   start in `Solver.BoundaryMode`.
2. **File format priority**: Gmsh `$NodeData`, VTU, or both in the first PR?
   Recommendation: both, ASCII, scalar, with a shared parser interface.
3. **Preconditioner scope** (§3.3): field-aware, or constant with a documented note?
4. **Semantics**: global replacement vs. per-attribute restriction. Recommendation:
   global replacement, as above.
5. **Refinement/3D→2D**: confirm deferring to Phase 2 is acceptable.
6. **Schema version bump**: confirm whether adding optional properties needs a SchemaVer
   bump and manifest update.

## 9. File-change checklist

New:

- `palace/fem/materialfield.hpp` / `.cpp`
- `palace/utils/gridfunctionio.hpp` / `.cpp`
- `palace/fem/qfunctions/22/*_field_*_qf.h` (six files)
- `test/unit/test-materialfield.cpp`
- `test/data/...` small Gmsh/VTU reader fixtures

Modified:

- `palace/fem/integrator.hpp` (field pointer, `Assemble` signature)
- `palace/fem/bilinearform.cpp` (fetch/pass field data)
- `palace/fem/libceed/integrator.hpp` / `.cpp` (`CeedQFunctionInfo`, `AssembleCeedOperator`)
- `palace/fem/integ/{vecfemass,mass,diffusion,curlcurl,curlcurlmass,diffusionmass,mixedvecgrad}.cpp`
- `palace/fem/qfunctions/*_qf.h` (include the new variants)
- `palace/models/modeoperatorassembly.hpp/.cpp`, `boundarymodeoperator.hpp/.cpp`,
  `modeeigensolver.hpp/.cpp`
- `palace/drivers/boundarymodesolver.cpp` (+ `.hpp` if needed)
- `palace/utils/configfile.hpp/.cpp`
- `scripts/schema/config-schema.json`
- `docs/src/guide/problem.md`, `CHANGELOG.md`
- `test/unit/CMakeLists.txt`, `test/unit/test-boundarymodeoperator.cpp`
