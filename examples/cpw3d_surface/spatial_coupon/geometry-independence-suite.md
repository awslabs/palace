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
digest. The bounded record additionally hashes the runtime, mesher, adaptor,
and MMG stages. Records and meshes must be content-distinct by SHA-256 across
matrix entries. The runner rereads every producer record, independently parses
the Gmsh mesh and rechecks topology/labels, and rejects normalized values that
differ from producer output.

Create the bounded mesh first:

```sh
python3 run_bounded_mesher.py \
  --seconds 1800 --memory-gib 8 \
  --log /large/run/CASE--VARIANT/mesh.log \
  --artifact /large/run/CASE--VARIANT/candidate.msh \
  --require-complete-toolchain \
  --tool runtime=/absolute/path/to/runtime \
  --tool mesher=/absolute/path/to/mesher \
  --tool adaptor=/absolute/path/to/adaptor \
  --tool mmg=/absolute/path/to/mmg3d \
  -- /absolute/path/to/runtime /absolute/path/to/mesher EXACT_FROZEN_ARGUMENTS
```

Then run `general_mesh_audit_producer.py` once for each required kind. Each
invocation takes `CASE VARIANT MESH INPUT_HASHES_JSON TRANSFORM_JSON OUTPUT` and
the kind-specific `--contract`, `--recipe`, `--signature`, `--identity-mesh`,
or `--bounded-report` arguments. Normalize only those five outputs:

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

## Next bounded mesh producer

The next bounded command is the identity four-edge Gmsh seed only. Run from the
repository root with fresh external output:

```sh
python3 examples/cpw3d_surface/spatial_coupon/run_bounded_mesher.py \
  --seconds 1800 --memory-gib 8 \
  --log /tmp/coupon-generality-four-edge-identity/mesh.log \
  --artifact /tmp/coupon-generality-four-edge-identity/candidate.msh \
  -- /Users/simlap/.juliaup/bin/julia --startup-file=no \
  --project=test/examples \
  examples/cpw3d_surface/spatial_coupon/mesh_graded_tet_experiment.jl \
  examples/cpw3d_surface/spatial_coupon/testdata/four-edge-9d2cb9bbb3fe \
  fabricated /tmp/coupon-generality-four-edge-identity/candidate.msh \
  1.0 0.02 0.25 --process \
  examples/cpw3d_surface/spatial_coupon/testdata/four-edge-9d2cb9bbb3fe/process.toml \
  --element-interface-slots
```

This does not yet create version-3 production evidence: no local MMG/adaptor
executable is available to freeze, and the rotated input still needs explicit
producer support. Do not launch Palace response, p4/p5, library, HPC, or physics
jobs. No mesh-only release claim is made.
