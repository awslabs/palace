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

1. `bounded-run`: `run_bounded_mesher.py` records the exact command,
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

## Exact next complete mesh-only run

Run the identity four-edge chain below from the repository root after setting
`EDGE_METRIC_ADAPTER_EXE` to the reviewed native `adapt_edge_metric.cpp` build.
The current host has no such executable, so this command is exact but presently
blocked before execution. Use a fresh `ATTEMPT` directory.

```sh
set -eu
ATTEMPT=/tmp/coupon-generality-four-edge-complete-v1
ADAPTER="$EDGE_METRIC_ADAPTER_EXE"
test ! -e "$ATTEMPT" && test -x "$ADAPTER"
JULIA=/Users/simlap/.juliaup/bin/julia
PYTHON=/opt/homebrew/bin/python3
SOURCE=examples/cpw3d_surface/spatial_coupon/testdata/four-edge-9d2cb9bbb3fe
RUN=examples/cpw3d_surface/spatial_coupon/run_bounded_mesher.py
SEEDER=examples/cpw3d_surface/spatial_coupon/mesh_spatial_coupon.jl
PREPARE=examples/cpw3d_surface/spatial_coupon/prepare_edge_metric_scout.py
RESTORE=examples/cpw3d_surface/spatial_coupon/restore_planar_metric_mesh.py
PUBLISH=examples/cpw3d_surface/spatial_coupon/relabel_frozen_interface_mesh.jl
mkdir -p "$ATTEMPT"

$PYTHON "$RUN" --seconds 1800 --memory-gib 8 \
  --log "$ATTEMPT/seed.log" --stage seed-generation \
  --artifact "seed-mesh=$ATTEMPT/seed.msh" \
  --tool "runtime=$JULIA" --tool "mesher=$SEEDER" -- \
  "$JULIA" --startup-file=no --project=test/examples "$SEEDER" \
  "$SOURCE/mesh-signature.csv" fabricated "$ATTEMPT/seed.msh" \
  --mask "$SOURCE/plan-view-mask.csv" --boundary "$SOURCE/plan-view-boundary.csv" \
  --radius 2.0 --metal-thickness 0.1 --overetch 0.05 --sidewall-angle 90 \
  --top-radius 0 --bottom-radius 0 --lc-fine 0.025 --lc-tangent 0.1 \
  --lc-far 0.16 --mesh-order 1 --max-nodes 4000000 --max-elements 4000000

$PYTHON "$RUN" --seconds 1800 --memory-gib 8 \
  --log "$ATTEMPT/metric.log" --stage metric-preparation \
  --input "seed-mesh=$ATTEMPT/seed.msh" \
  --artifact "mmg-seed=$ATTEMPT/metric/seed.mesh" \
  --artifact "metric=$ATTEMPT/metric/metric.f64" \
  --artifact "pins=$ATTEMPT/metric/pins.txt" \
  --artifact "fixed-triangles=$ATTEMPT/metric/fixed-triangles.txt" \
  --artifact "restoration-recipe=$ATTEMPT/metric/recipe.json" \
  --tool "runtime=$PYTHON" --tool "metric-preparer=$PREPARE" -- \
  "$PYTHON" "$PREPARE" "$ATTEMPT/seed.msh" "$ATTEMPT/metric" \
  --normal 0.025 --tangent 0.1 --far 0.16 --protected-distance 0.05 \
  --protect-surface 0.05 --semantic-contract "$SOURCE/semantic-contract.json"

$PYTHON "$RUN" --seconds 1800 --memory-gib 8 \
  --log "$ATTEMPT/adapt.log" --stage native-adaptation-mmg \
  --input "mmg-seed=$ATTEMPT/metric/seed.mesh" \
  --input "metric=$ATTEMPT/metric/metric.f64" \
  --input "pins=$ATTEMPT/metric/pins.txt" \
  --input "fixed-triangles=$ATTEMPT/metric/fixed-triangles.txt" \
  --artifact "adapted-mesh=$ATTEMPT/adapted.meshb" \
  --tool "runtime=$ADAPTER" --tool "adapter-mmg=$ADAPTER" -- \
  "$ADAPTER" "$ATTEMPT/metric/seed.mesh" "$ATTEMPT/metric/metric.f64" \
  "$ATTEMPT/metric/pins.txt" "$ATTEMPT/adapted.meshb" 0.025 0.16 1.3 \
  freeze-selected "$ATTEMPT/metric/fixed-triangles.txt" 1e-8

$PYTHON "$RUN" --seconds 1800 --memory-gib 8 \
  --log "$ATTEMPT/restore.log" --stage label-restoration \
  --input "adapted-mesh=$ATTEMPT/adapted.meshb" \
  --input "restoration-recipe=$ATTEMPT/metric/recipe.json" \
  --artifact "restored-mesh=$ATTEMPT/restored.msh" \
  --tool "runtime=$PYTHON" --tool "label-restorer=$RESTORE" -- \
  "$PYTHON" "$RESTORE" "$ATTEMPT/adapted.meshb" \
  "$ATTEMPT/metric/recipe.json" "$ATTEMPT/restored.msh" --max-displacement 1e-8

$PYTHON "$RUN" --seconds 1800 --memory-gib 8 \
  --log "$ATTEMPT/publish.log" --stage final-gmsh-publication \
  --input "restored-mesh=$ATTEMPT/restored.msh" \
  --artifact "candidate-mesh=$ATTEMPT/candidate.msh" \
  --artifact "ownership-partition=$ATTEMPT/candidate.msh.interface-partition.csv" \
  --tool "runtime=$JULIA" --tool "publisher=$PUBLISH" -- \
  "$JULIA" --startup-file=no --project=test/examples "$PUBLISH" \
  "$SOURCE" fabricated "$ATTEMPT/restored.msh" "$ATTEMPT/candidate.msh"
```

This is mesh-only. Do not launch Palace response, p4/p5, library, HPC, or
physics jobs. Rotation still requires an explicit transformed source/producer;
curved and multilayer cases need corresponding CAD-preserving native-stage
evidence before any release claim.
