<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Mesh-only geometry-independence suite

`geometry-independence-suite.json` is the machine-readable inventory and gate
contract for the first phase of the coupon-library acceleration plan. It does
not qualify response fields or a production library.

The inventory freezes each available input file by SHA-256 and derives edge,
slot, and conductor counts from the signature CSV. The runner contains no
case-specific source count, edge count, box bound, or model hash. It records the
current local three-edge calibration without changing its source directory.
The known four- and ten-edge campaign cases are deliberately listed as required
but unavailable locally; they cannot silently fall back to a smaller fixture.
The repository six-edge fixture is likewise required and hash checked. Rotated
variants, transition/multislot coverage, and supplemental concave-slot and hole
fixtures are declared in data rather than selected by edge-count branches.

## Preflight

Write all output outside the repository and use a fresh path:

```sh
python3 run_general_mesh_suite.py \
  --manifest geometry-independence-suite.json \
  --preflight-only \
  --root /tmp/coupon-generality-preflight
```

The current local result is expected to fail closed before meshing because the
four- and ten-edge campaign input directories are not present. Available
1/2/3/6-edge inputs are still hash checked and inventoried in `summary.json`.
Supply a frozen external input only with `--input CASE=DIRECTORY`; the manifest
must already contain SHA-256 values for every consumed file, so an unrecorded
replacement remains rejected.

## Mesh evidence

The existing discovery and meshing path remains authoritative:
`mesh_graded_tet_experiment.jl` / `prepare_edge_metric_scout.py`, bounded by
`run_bounded_mesher.py`, with `audit_mesh_measures`, `audit_coupon_mesh`,
`audit_edge_metric_mesh.py`, and interface-partition output. Do not fork those
geometry or metric pipelines. Normalize their read-only results into one JSON
file per manifest case under an external audit directory, then run:

```sh
python3 run_general_mesh_suite.py \
  --manifest geometry-independence-suite.json \
  --audit-root /tmp/coupon-generality-audits \
  --root /tmp/coupon-generality-summary
```

Each evidence file must contain these sections:

- `ExactLabelsMaterials`: exact expected/actual volume and boundary attributes;
- `OwnershipClosure`: `UnmatchedPolicy: Error`, zero unmatched/overlaps, and
  exhaustive closure;
- `SemanticCorners`: expected/actual physical corner coordinates;
- `ProtectedSurfaces`: exact protected-support identities and zero changes;
- `AchievedAnisotropy`: sampled two-transverse and tangential widths, including
  target normal size and sample count;
- `RotationCovariance`: compared variant and maximum invariant error;
- `TraceDiagonal`: zero detected global diagonal refinement bands;
- `Resources`: exit status, elapsed seconds, peak process-tree RSS, and elements.

All tolerances and resource ceilings are declared in the manifest. Missing
inputs, hashes, evidence, fields, or samples fail closed. Mesh and audit outputs
remain external; existing AMR jobs, campaign results, libraries, and generated
user files are not read as mutable working state or modified.
