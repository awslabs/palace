<!--
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
-->

# Fixed-wedge finite-metal singular elements

Palace provides an opt-in approximation for three-dimensional finite-thickness PEC edges
following the wedge-superposition treatment of Elkin et al. (arXiv:2606.18140).

```json
"SingularElements": {
  "Attributes": [5, 6],
  "Order": 3,
  "FiniteMetalModel": "FixedWedgeEdgeSuperposition",
  "FixedExponent": 0.6666666666666666
}
```

The selected attributes identify candidate finite-metal faces. Palace walks the tetrahedral
material fan around each selected surface edge and:

- enriches dielectric openings within `1e-3` radians of 270 degrees;
- leaves 90-degree and smooth 180-degree openings to standard FEM;
- leaves other non-target openings to standard FEM and reports their count;
- reports edges where one selected face meets an unselected boundary, such as an excluded
  artificial port patch.

Every selected 270-degree segment receives the configured fixed real exponent. Collinear
segments remain one straight feature across material-sector changes. Three or more fixed
segments may meet at a `JUNCTION`; the incident edge features remain independent and share
the existing single node basis keyed by the physical mesh vertex.

This is a wedge superposition approximation, not an exact trihedral point solution. Runtime
output explicitly reports the extracted segment/junction counts and states that exact
trihedral point modes are absent.

The initial implementation is lossless. Bulk loss tangent, conductivity, impedance/far-field
boundaries, and resistive lumped-port terms are rejected. Interface dielectric loss tangents
remain available for SA/MS/MA participation and Q postprocessing. Low-loss complex-material
support requires a separate validated extension.

Tetrahedral singular matrices support two integration strategies. `Adaptive` uses entrywise
`AbsTol`, `RelTol`, and `MaxSubdivisions`. `FixedSubdivision` applies `QuadratureOrder`
uniformly after `Subdivisions` recursive one-to-eight refinements; order 8 and six
subdivisions reproduce the strategy described by Elkin et al. Fixed subdivision is 3D-only
and ignores adaptive tolerances.

For affine tetrahedra, fixed-subdivision electrostatic H1 assembly evaluates complete
standard-enrichment and enrichment-enrichment reference tensor blocks in one batched
quadrature traversal and caches them by the exact element-local singular basis pattern.
Every subsequent element with that pattern only contracts the cached tensor with its
Jacobian metric. Runtime diagnostics report memory and persistent hits, generated patterns,
tensor entries, and actual generated reference leaves.

Set `Solver.SingularElements.ReferenceCache` to a shared directory to persist these tensors
between Palace launches and AMR operator rebuilds. Cache keys include the MFEM and reference
convention versions, finite-element order, complete basis descriptors and exponent bit
patterns, and fixed quadrature order/depth. Files are validated before use and installed by
an atomic no-replace hard link.

For electrostatic H1 assembly, Palace collectively discovers persistent reference patterns
before element assembly, globally deduplicates their complete keys, and assigns each pattern
to one rank that reported a local use. A deterministic largest-work-first schedule balances
predicted tensor counts while avoiding duplicate cold-cache generation. The cache directory
must be coherently visible to every rank. During a cold prewarm, rank zero prints an immediate
summary and a heartbeat every 60 seconds with elapsed time and visible cache files. Set
`PALACE_REFERENCE_HEARTBEAT_SECONDS` to a positive value up to 3600 to change the interval.
The heartbeat thread performs no MPI calls or quadrature-loop instrumentation and is enabled
only when MPI provides at least `MPI_THREAD_FUNNELED`.

A cold depth-six generation remains expensive, but a warm depth-six cache performs no
recursive reference quadrature and has the same contraction cost as depth two.

`Subdivisions` controls matrix/reference-tensor construction only. Electrostatic total
energy, including the denominator of every surface participation ratio, is evaluated from
the algebraic quadratic form `V^T K V`. If per-domain energies are requested,
`PostprocessingSubdivisions` controls their independent direct-field integration and
defaults to two. Three-dimensional surface participation continues to use edge-aligned
Gauss-Jacobi quadrature and does not recursively follow either subdivision depth.

`TransmissionWedge` remains the default finite-metal model and preserves the existing
material-dependent exponent, loss-tangent validation, and strict junction behavior.
