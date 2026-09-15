<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Mesh-only geometry-independence suite

`geometry-independence-suite.json` is a fail-closed inventory and gate contract
for the mesh-only phase of the coupon-library acceleration plan. It does **not**
qualify response fields, trace accuracy, a production library, or all release
gates.

## Frozen inputs and required matrix

Every runnable case requires SHA-256-frozen `Signature`, `Boundary`, `Mask`,
`Process`, `SemanticContract`, and `MeshRecipe` roles. Expected materials,
labels/adjacency, corners, and protected supports come only from that contract.
Every case has identity and `rotate-z-0.63` variants with explicit transforms
and a fixed comparison pair. Concave/multislot, hole, rounded/filleted, and
opposed-layer controls are ordinary required cases. Feature-scaling and
CAD-subdivision-sensitivity comparisons are declared in the manifest.

The matrix contains 12 cases and 24 required case/variant entries. All 12 cases
now have hash-frozen local source contracts. A bounded read-only assessment on
`soca-green` copied only the approved source-contract files from campaign inputs
`07` and `09`; it copied no mesh, field, solution, response matrix, or source
bank. Campaign input `07` maps to
`spatialedgecluster_edgecount-4_9d2cb9bbb3fe` and has exactly four signature
rows. Input `09` maps to `spatialedgecluster_edgecount-10_6791f1c84123` and has
exactly ten. The process-library model names, local/remote SHA-256 values, trace
mesh references, slots, conductors, topology, and target-repository Apache-2.0
license were reviewed. Their semantic contracts derive label families and
slot/conductor coverage from the source process model and semantic corners from
plan-view vertices explicitly classified `Physical`.

## Preflight

```sh
python3 run_general_mesh_suite.py \
  --manifest geometry-independence-suite.json \
  --preflight-only \
  --root /tmp/coupon-generality-preflight
```

This now succeeds for all 12 source cases and verifies the four-/ten-edge row
counts as four and ten. `--input CASE=DIRECTORY` cannot bypass a mismatched
hash.

## Production evidence chain

A version-3 normalized record is assembled from five separate records:

1. `bounded-run`: `run_bounded_mesher.py` records the exact command, its
   absolute working directory (so relative argv paths resolve unambiguously),
   constrained environment, launcher digest, resource measurements, and output
   mesh digest.
2. `mesh-topology-quality`: the real Gmsh mesh is parsed and audited for
   topology, labels/material adjacency, corners, protected labels, achieved
   directional widths, diagonal bands, and Jacobian quality.
3. `mesh-complexity`: H1 DOFs are counted from mesh topology/order and semantic
   feature/CAD-subdivision counts are derived from the signature.
4. `mesh-invariants`: material volumes and labeled surface areas are integrated
   from the mesh.
5. `variant-transform`: candidate points, connectivity, and labels are checked
   against the exact 4x4 transform of the independently audited identity mesh.

Every record binds the same mesh digest, every manifest source digest, exact
transform and digest, case/variant, command/environment, and frozen producer
digest. The bounded record consumes five separately executed reports forming a
digest-linked DAG: seed generation, metric preparation, native adapter/MMG,
label restoration, and final Gmsh publication. A tool is accepted only when its
frozen path occurs in that stage's executed argv. The DAG binds the seed mesh,
MMG seed, tensor metric, pins, fixed triangles, restoration recipe, adapted
mesh, restored mesh, final candidate, and ownership partition. Merely declaring
an adaptor/MMG file cannot pass.

The topology audit reads Gmsh physical volume names, derives ownership from the
publisher's partition report, compares protected surface measures to the seed,
measures corner/subdivision/cut neighborhoods, and detects connected short-edge
surface bands rather than asserting constants. Feature and CAD-subdivision
counts are frozen in `FeatureTopology` and independently reproduced from finite
oriented signature segments plus boundary `Physical`/`Continuation` classes.
Records and meshes must be content-distinct by SHA-256 across matrix entries.

Run `general_mesh_audit_producer.py` once for each required kind, passing all
five `--stage-report STAGE=PATH` bindings to `bounded-run` and
`mesh-topology-quality`. Normalize only those five audit outputs:

```sh
python3 normalize_general_mesh_evidence.py \
  geometry-independence-suite.json CASE VARIANT \
  /large/run/CASE--VARIANT/candidate.msh \
  /large/audits/CASE--VARIANT.json \
  --audit-record bounded-run=/large/run/CASE--VARIANT/bounded.json \
  --audit-record mesh-topology-quality=/large/run/CASE--VARIANT/topology.json \
  --audit-record mesh-complexity=/large/run/CASE--VARIANT/complexity.json \
  --audit-record mesh-invariants=/large/run/CASE--VARIANT/invariants.json \
  --audit-record variant-transform=/large/run/CASE--VARIANT/transform.json
```

After all 24 records exist:

```sh
python3 run_general_mesh_suite.py \
  --manifest geometry-independence-suite.json \
  --audit-root /large/audits \
  --root /large/coupon-generality-summary
```

No release claim is permitted until all production records and scaling
comparisons pass.

## Trace tooling boundary

The historical 135-column/40-constrained/95-free tools remain explicitly
outside this gate. `trace_audit_contract.py` requires unique excitation names,
explicit kinds, smooth/global and localized-matching roles, exact coverage of
every declared conductor/interface class, and a fixed finite-coefficient
combination. It remains only a future arbitrary-source seam; no three/six/ten-
edge physics claim is made.

## Canonical build and mandatory placement publication

Every proper-rigid placement now uses the seven-stage contract in
`mesh_stage_contract.py`. The first six stages are one reusable source-local
canonical build: identity source validation, seed generation, metric preparation,
native MMG adaptation, label restoration, and canonical Gmsh publication. The
seventh stage, `proper-rigid-publication`, is mandatory even for identity.

`canonical_mesh_build.py` computes the cache key from all immutable source,
process, semantic, recipe, and gate hashes plus the exact canonical tool-role
set derived from the six stages (runtime, validator, mesher, metric preparer,
adaptation wrapper, adapter, runtime-resolved MMG library, restorer, publisher).
Its build record is derived from the six validated stage reports: it binds the
exact fifteen canonical artifact roles (every stage output, including the
adaptation receipt and the source-local restored mesh) and each stage report's
SHA-256, and rejects missing or extra roles. A placement may share the six
canonical reports only when the cache key, build hash, every canonical artifact
hash and every stage-report hash match. Final meshes, transform receipts,
ownership outputs, and audit records remain per-placement and content-distinct.
Canonical and placement seconds/RSS are reported separately.

Metric preparation classifies seed surface features by geometry, not by label:
`edge_volume_metric.surface_features` makes a shared triangle edge a feature
(MMG ridge, pin at graph turns/junctions, anisotropic `PhysicalSegment`,
protected band) only when its incident triangles are not coplanar within the
fixed dimensionless `COPLANAR_TOLERANCE` (sine of the normal angle and plane
offset relative to the local triangle diameter). A reference change between
coplanar triangles, such as a slot-ownership split of one plane, remains a
plain reference boundary that MMG preserves without ridge or corner
constraints; per-support relabeling and projection in the restorer are
unchanged. Label restoration compares repair displacements to the bound with a
roundoff-only relative allowance (`DISPLACEMENT_ROUNDOFF_TOLERANCE`, 1e-12) so
a vertex clamped onto the displacement ball is not misreported as an overshoot;
any true overshoot and any semantic-corner result that still misses the corner
target are rejected components (`RejectedQualityRepairComponents`,
`RejectedCornerRepairs`), never partial commits. The final scaled-Jacobian,
corner-aspect and displacement gates are unchanged.

Native adaptation must be launched through `run_native_mmg_adaptation.py`.
There is no caller-provided hmax: the wrapper reads
`FarFieldBudgetPolicy.EffectiveFarSize` from the bound restoration recipe,
requires it to equal the recipe's `FarSize`, passes that exact value to the
reviewed adapter, and emits an adaptation receipt containing the executed argv.
The wrapper also resolves the adapter's single MMG3D shared library from its
link table and rpath entries (`otool` on macOS, `readelf` elsewhere), requires
it to be the `--mmg-library` bound as the stage's `mmg-library` tool role,
removes dynamic-loader search overrides from the adapter environment, and
records the library path, link reference and SHA-256 in the receipt. The stage
contract requires the receipt's adapter and library hashes to equal the bound
tools, so replacing the loaded library changes `CanonicalBuildId`. Preflight
resolution:

```sh
python3 -c 'from run_native_mmg_adaptation import resolve_mmg_library; \
  print(resolve_mmg_library("$EDGE_METRIC_ADAPTER_EXE"))'
```

The normal, tangent, corner, protected-band targets and 4,000,000-element cap
remain unchanged.

Every stage command must name exactly its bound inputs and outputs: seed
generation binds the source signature/boundary/mask, canonical publication binds
the source process/signature/boundary (`--process/--signature/--boundary`), and
proper-rigid publication binds source semantic/signature/boundary/mask/process,
the canonical build record and every ownership/semantic/support output. Named
options must equal the bindings and positionally consumed paths must occur in
argv; `run_bounded_mesher.py` refuses to launch otherwise. No stage accepts a
source directory.

`publish_rigid_coupon_mesh.py` requires a fresh output and a finite orthonormal
homogeneous transform with determinant +1. It changes only Gmsh 2.2 node
coordinates, preserving node tags, cell-block order, connectivity, physical and
geometrical tags, physical names, and point/cell/field data. It reconstructs the
transformed semantic and support contracts from immutable source, checks exact
`R*x+t`, positive orientation, quality and invariant measures, and invokes
`audit_rigid_coupon_ownership.jl` with explicit `--process`, `--signature` and
`--boundary` files. That audit classifies quadrature points in the inverse
source-local frame and rejects any final owner-label mismatch or
nonpositive/nonclosing owner partition. The receipt records the ownership argv,
runtime/auditor hashes, output hashes and every source-input hash; the stage
contract requires the argv options to equal the placement bindings and the
hashes to equal the bound inputs/outputs. Identity therefore gets a fresh path
and receipt; reflection remains rejected.

A complete four-edge and ten-edge evidence run still requires an executable
reviewed build of `adapt_edge_metric.cpp` through `EDGE_METRIC_ADAPTER_EXE` and
its resolved MMG3D library bound as `mmg-library`. The fixture tests build a
tiny native adapter linked against a fixture `libmmg3d`
(`testdata/build_tiny_native_adapter.py`), so they require a C compiler.
Use fresh directories, 1800 seconds and 8 GiB for each bounded stage (16 GiB for
the aggregate audit), and the unchanged 4M element cap. Build the four-edge
canonical candidate once, publish identity and rotate-z through stage seven,
then build only ten-edge identity. Do not run ten-edge rotation, reflection,
physics, p4/p5, library, or HPC work in this mesh-only gate.
