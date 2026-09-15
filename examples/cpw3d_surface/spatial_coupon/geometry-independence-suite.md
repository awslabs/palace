<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Mesh-only geometry-independence suite

`geometry-independence-suite.json` is a fail-closed inventory and gate contract
for the mesh-only phase of the coupon-library acceleration plan. It does **not**
qualify response fields, trace accuracy, a production library, or all release
gates.

## Frozen inputs and required matrix

Every runnable case requires SHA-256-frozen `Signature`, `Boundary`, `Mask`,
`Process`, `SemanticContract`, and `MeshRecipe` roles. A semantic contract must
contain nonempty materials, boundary labels and adjacency, semantic corners,
protected supports, metric/cut roles, and `UnmatchedPolicy: Error`. Expected
values come only from this contract, never from evidence under review.

Every case has identity and `rotate-z-0.63` variants with explicit 4x4
transforms and an exact comparison pair. Concave/multislot, hole,
rounded/filleted, and opposed-layer controls are ordinary required cases, not
ignored supplemental metadata. The one-edge CAD-subdivision control and the
one-to-six-edge feature comparison are declared scaling gates.

The checked-in matrix currently contains 12 cases and 24 required
case/variant entries. Ten cases (20 entries) have hash-frozen local source
contracts. Four-edge `008-spatialedgecluster-edgecount-4-9d2cb9bbb3fe` and the
unverified ten-edge candidate `exact-p5-6791f1c84123` remain required but
unavailable. They have null input hashes and cannot be overridden until their
actual signature edge counts, provenance, license, and semantics are verified.
No non-ten case is relabeled as ten-edge.

## Real preflight

Write output outside the repository and use a fresh path:

```sh
python3 run_general_mesh_suite.py \
  --manifest geometry-independence-suite.json \
  --preflight-only \
  --root /tmp/coupon-generality-preflight
```

On a host without the two `/data/home/simlap` libraries this must return nonzero.
It still hash-checks all available 1/2/3/6-edge and ordinary geometry controls;
only the four- and ten-edge records should fail. `--input CASE=DIRECTORY` cannot
bypass null or mismatched manifest hashes.

## Producing bound evidence

Use the existing bounded mesher and audits. Metric preparation and audit now
require a frozen semantic contract; numeric materials, matching labels,
multislot/multiconductor interfaces, and adjacency are data rather than
`3100/5001/6001` branches. The old one-slot convention remains available only
through the explicit fail-closed `--simple-sharp-contract` compatibility flag.

Normalize producer output only after the producer has written a case/variant
record, mesh, and independent audit record:

```sh
python3 normalize_general_mesh_evidence.py \
  geometry-independence-suite.json CASE VARIANT \
  /large/run/CASE--VARIANT-raw.json \
  /large/run/CASE--VARIANT.msh \
  /large/audits/CASE--VARIANT.json \
  --audit-record /large/run/CASE--VARIANT-raw.json
```

The normalizer binds `CaseId`, `Variant`, the canonical transform hash, every
input hash, process, semantic contract, mesh recipe, all frozen evidence tools,
audited mesh bytes, and the independent producer record. The gate recomputes artifact hashes and checks that
each producer record names the same case/variant. Audit records cannot be reused
between matrix entries.

A normalized record contains actual values only:

- exact volume attribute/material names, boundary labels, and material adjacency;
- exhaustive ownership with zero unmatched and overlap counts;
- semantic corners and protected supports;
- two-transverse/tangential achieved widths and zero global trace diagonals;
- nonnegative exit time, process-tree RSS, and element count;
- scaled-Jacobian and Jacobian-condition measurements;
- H1 DOFs, semantic feature count, and CAD-subdivision count;
- independently measured comparison invariants.

After every matrix record exists, run:

```sh
python3 run_general_mesh_suite.py \
  --manifest geometry-independence-suite.json \
  --audit-root /large/audits \
  --root /large/coupon-generality-summary
```

The runner computes rotation and scaling comparisons across independent records;
it does not accept a self-asserted covariance result.

## Trace tooling boundary

`prepare_trace_projection_audit.py` and
`summarize_trace_projection_audit.py` are historical model-specific diagnostics
for the 135-column, 40-constrained/95-free calibration. They are explicitly
outside this geometry-independence gate. `trace_audit_contract.py` is the tested
replacement seam for arbitrary source banks, constrained subsets, and semantic
excitation classes; no three/six/ten-edge physics claim is made here.

## Next bounded run

After verified four- and ten-edge source contracts are frozen, run only the
mesh matrix above with `run_bounded_mesher.py --seconds 1800 --memory-gib 8`, one
process per case/variant, and no more than six local workers. Keep meshes,
producer JSON, logs, and summaries under `/large` (or another external attempt
directory). Do not launch Palace response, HPC, p4/p5, or library solves in this
step. A release claim remains forbidden until all 24 evidence records and both
scaling comparisons pass.
