# Idea: finite NC face-map trace plan for boundary ParaView point fields

## Why this exists

The current nonconforming-AMR ParaView fallback is diagnostic, not a final design. We do **not** want a permanent `mesh.Nonconforming()` policy branch. The goal of this experiment is to make the opt-in libCEED boundary point-field path fast on NC AMR by operating at the correct abstraction: a boundary trace plan keyed by the finite NC face relation, not by individual boundary elements or floating-point point clouds.

The current bad path is exposed by running the AMR1 transmon ParaView-only benchmark with:

```bash
PALACE_CEED_NONCONFORMING_PARAVIEW=1 PALACE_SURFACE_PROFILE=1 ./.auto/measure.sh
```

Before the fallback, this path measured roughly:

```text
paraview_seconds ~1226-1246 s
operator_construction_seconds ~52 s
ParaView HWM ~60 GB
```

After the fallback, the normal benchmark is ~117 s, but that wins by avoiding the problematic libCEED NC boundary path. This experiment should be judged against the **opt-in no-fallback baseline**, not against the fallback best.

## Hypothesis

NC boundary output is slow because `SurfaceFunctional::AssembleLocal` loses the discrete NC face-map identity and treats mapped face point sets as arbitrary per-element mapped integration rules. In the non-AtPoints path, the grouping key intentionally includes:

```cpp
key.push_back(static_cast<long long>(mapped_group_id++));
```

which makes each mapped group unique. That is safe, but it prevents batching. For an NC slave face, the master-face submap is deterministic and comes from MFEM/NCMesh. With `MaxNCLevels=1`, there should be only a small finite set of slave placements (plus face orientation/local-face variants), not one operator per boundary element.

MFEM exposes this relation:

```cpp
mfem::Mesh::FaceInformation face = pmesh.GetFaceInformation(face_id);
face.point_matrix;                 // NC slave position within master face
face.ncface;
face.element[k].local_face_id;
face.element[k].orientation;
face.element[k].conformity;
face.tag;
```

and lower-level:

```cpp
mfem::Mesh::NCFaceInfo::PointMatrix
```

So the GPU path should group by this finite topology/map class, not by per-face floating-point points.

## Explicit functionals/fields involved

The regression is from private boundary visualization point fields using `SurfaceFunctional`, not the public reducing surface functionals. The implicated internal `KernelKind` values are:

```cpp
BDR_FIELD_E       // boundary E_real/E_imag or E
BDR_FIELD_B       // boundary B_real/B_imag or B
BDR_FLUX_Q        // Q_s_real/Q_s_imag or Q_s
BDR_CURRENT_J     // J_s_real/J_s_imag or J_s
BDR_ENERGY_E      // boundary U_e
BDR_ENERGY_M      // boundary U_m
BDR_POYNTING      // boundary S
```

Concrete qfunction families:

```cpp
f_eval_bdr_field_1_32 / f_eval_bdr_field_2_32
f_eval_bdr_flux_q_1_32 / f_eval_bdr_flux_q_2_32
f_eval_bdr_current_j_1_32 / f_eval_bdr_current_j_2_32
f_eval_bdr_energy_1_32 / f_eval_bdr_energy_2_32
f_eval_bdr_poynting_1_32 / f_eval_bdr_poynting_2_32
```

Autoresearch already showed `BDR_FIELD_E`/`BDR_FIELD_B` alone are a large part of the problem: re-enabling nonconforming boundary E/B pointfields while leaving the other boundary fields on legacy coefficients regressed from ~117 s to ~383 s.

## Concrete implementation plan

### 1. Instrument first

Run the no-fallback opt-in benchmark with `PALACE_SURFACE_PROFILE=1` and extend the existing profile output if needed. For each boundary buffer `KernelKind`, print:

- total groups, elements;
- AtPoints groups/elements;
- non-AtPoints groups/elements;
- two-sided group count;
- ghost/face-neighbor group count;
- number of unique NC map keys;
- max/median elements per group.

This confirms exactly which `BDR_*` kinds are exploding.

### 2. Introduce a finite NC trace-map key

In `SurfaceFunctional::AssembleLocal`, when building the `FaceConfigKey`, replace the per-element `mapped_group_id++` uniqueness for non-AtPoints groups with a deterministic key derived from MFEM face information.

Suggested key ingredients:

```text
bdr_geom
vol_geom_a / vol_geom_b
face geometry
side tag (a/b)
local face id on the volume element
face orientation
NC role: conforming / slave / ghost slave / master/subset
NC point-matrix class
nq
flip flag
side_scale * 2
normal_scale * 2
kernel kind or two-sided qfunction family if needed
```

For the NC point-matrix class:

- Prefer a stable canonical id if available from MFEM (`ncface` + `PointMatrix` identity is enough within one mesh instance only if it really maps to the unique matrix table).
- Otherwise hash the `DenseMatrix` dimensions and entries after quantization to the expected dyadic grid (for `MaxNCLevels=1`, entries are typically 0, 1/2, 1; use conservative exact/rounded hashing).
- Include face orientation/local-face id so oriented variants do not alias.

Do **not** use coordinate-fuzzy arbitrary point keys as the primary abstraction.

### 3. Prove key correctness while merging

When a key already exists and we are about to reuse the registered mapped integration rule, compare the newly computed mapped points against the stored representative points for that key. If they differ beyond a tiny tolerance, do not silently produce wrong output:

- either refine the key by adding the missing discrete face info;
- or abort in debug/profile mode with a clear message showing the mismatched key.

Avoid falling back to `mesh.Nonconforming()` or per-face uniqueness as the final path.

### 4. Prefer AtPoints batching where possible

For CUDA/libCEED AtPoints-capable cases, the operator should not need a separate mapped integration rule at all: mapped reference coordinates can be runtime point data. Check why the current code falls into the non-AtPoints path for the problematic NC boundary cases. If the blockers are ghost/two-sided handling, try splitting into one-sided local/ghost contributions that accumulate into the same output buffer so AtPoints remains usable.

### 5. Measurement target

Use the wrapper:

```bash
./.auto/measure_nc_trace.sh
```

which sets:

```bash
PALACE_CEED_NONCONFORMING_PARAVIEW=1
PALACE_SURFACE_PROFILE=1
```

A useful candidate is one that improves the opt-in no-fallback path substantially from ~1226 s and collapses group/operator counts, even if it does not immediately beat the current fallback result of ~117 s. The fallback result remains the performance reference, but this experiment is about removing the fallback by fixing the abstraction.

## Constraints

- Do not drop output fields.
- Do not optimize by changing AMR level, rank count, solver order, or disabling ParaView data.
- Do not ship a new data-dependent special case based on `mesh.Nonconforming()`.
- Preserve the existing scalar correctness checks and GPU/libCEED log checks.
- Use Spack only; package rebuilds through `.auto/measure*.sh` are okay.

## If this works

The final shape should be a reusable boundary trace-plan abstraction, not scattered conditionals:

```cpp
BoundaryTracePlan plan(mesh, bdr_marker, lod, semantics);
```

where each trace contribution knows:

```text
output slot
boundary element
source side / source element / ghost source
local face id and orientation
NC point-matrix/map key
side weight and normal sign
mapped reference points or AtPoints runtime coordinates
```

Then conforming and nonconforming meshes use the same path; NC only changes the finite trace-map key.
