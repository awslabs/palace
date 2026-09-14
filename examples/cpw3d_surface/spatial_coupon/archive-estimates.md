<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->
# Experimental read-only archive estimates

`estimate_archived_fields.py` launches Palace with
`PALACE_RESPONSE_ESTIMATE_ONLY=1`. This is a diagnostic reader, **not a new PDE
solve, AMR step, response-matrix reducer, or accuracy certificate**. Normal solver,
archive worker, and reducer paths remain unchanged. Use a separately built,
SHA256-bound executable; never replace the executable that produced the archives.

## Inputs and safety

Use the complete original prescribed-potential configuration, source files, mesh,
polynomial order, material/interface definitions, units, and MPI layout. Only
`Problem.Output` is changed by the launcher. Relative input paths are interpreted
relative to the invocation working directory, as in the Palace executable.
The launcher requires ordinary strict JSON, e.g. an existing resolved config.

The request has this structure (four-source synthetic example, not coupon data):

```json
{
  "Version": 1,
  "SourceIds": [1, 2, 3, 4],
  "SourceSHA256": ["<source1 SHA256>", "<source2 SHA256>", "<source3 SHA256>", "<source4 SHA256>"],
  "BasisContract": {"Path": "/absolute/path/basis-contract.json", "SHA256": "<contract SHA256>"},
  "ZeroTraceIndices": [3],
  "Validation": {
    "RelResidualTol": 1e-7,
    "AbsResidualTol": 2.220446049250313e-14,
    "BCAbsTolV": 1e-9,
    "MinEnergyJ": 1e-30,
    "MinCancellationRatio": 1e-8
  },
  "WriteElementIndicators": false,
  "Excitations": [
    {"Name": "source1", "Coefficients": [1, 0, 0, 0]},
    {"Name": "combination", "Coefficients": [2, -0.4, 0, 0]}
  ]
}
```

For the coupon assessment **all 135 source IDs/hashes and all 40 fixed constrained
indices must be present**, including unused columns. The complete dense coefficient
rows have exactly 135 entries, and every constrained coefficient must be exactly
zero. No magnitude cutoff is applied to nonzero coefficients. `BasisContract`
uses the preserved contract's `Sources`, `ZeroTraceIndices`, and
`OutputSourceSHA256` entries. The preserved `basis-NNNN.csv` filename convention
binds each original Index to its basename and hash, including unused sources; the
complete ordered mapping must match. Renaming source files is unsupported.
Contract and source hashes are verified and included in provenance.
Request coefficients themselves are hash-bound by the request file hash.

All validation parameters are explicit: the first three values above are the
approved initial **archive interpretation checks** for 1 V sources, not a solve or
discretization tolerance certificate. The energy/cancellation floors above are
illustrative numerical-reporting floors: choose and record them for the actual
probe. They never change the field, residual acceptance threshold, or raw indicator.
The launcher does not choose excitations, drop modes, automatically relax checks,
change recovery tolerance, or launch follow-up runs.

Example invocation after resource/probe approval:

```sh
python3 examples/cpw3d_surface/spatial_coupon/estimate_archived_fields.py CONFIG.json \
  --request REQUEST.json --archive /absolute/archive-view \
  --binary /absolute/palace-archive-estimate-SHA256.bin --binary-sha256 SHA256 \
  --output /absolute/new-run --ranks RANKS --seconds SECONDS --memory-gib GIB
```

`--output` must not exist, including as an empty directory, and must not overlap
any input. Both launcher and executable reject archive-only/reduce-only/recycle/
block-size flags. Both reject `Model.ExportPrerefinedMesh` (which otherwise writes
beside the input mesh), nonempty `Model.Partitioning`, and nonempty dielectric
`OwnershipDataFile`. These optional reconstruction/observable inputs are unsupported
in this experimental mode, not unrecorded exceptions. Revisit them explicitly before
any later NC/checkpoint workflow that needs them. The native mode also rejects AMR,
geometry-driven edge refinement, and surface response corrections. Ordinary a priori
mesh configuration is retained
for deterministic reconstruction. No archive is created or opened for writing.
Potential headers are checked for source, field kind, rank, rank count, width and
format; wrong widths are rejected **before allocation**, truncated/trailing payloads
and nonfinite values are rejected. Only used potential archives are read, never
archived D fields. Input hashes are compared before/after even when the run fails.

Use the launcher for assessed evidence: the lower-level environment mode performs
numerical validation but does not itself implement cryptographic SHA256 checks.
Native preflight runs before creating the solver output; all MPI ranks synchronize
before root creates it. There is no read-only guarantee against an unrelated
concurrent process modifying input files and restoring them between hash checks.
Use immutable archive storage and one writer for the output area.

## What is measured

For **each nonzero constituent** and then each final linear combination:

- Reproject the actual source definition using `LaplaceOperator::GetExcitationVector`.
- Compare archived essential true DOFs to these **imposed discrete Dirichlet values**
  using a global maximum in volts and in internal units. This is not the ideal-P1
  boundary representation error; retain the separate trace-projection audit.
  `BoundaryCategories` counts owned matching-only true DOFs, physical-ground true
  DOFs (including intersections), and trace/ground intersections. Separate maxima
  are reported for these categories, with null maxima if a category is absent.
  Physical ground takes precedence at intersections, as in the unchanged projection.
- Assemble and apply the original eliminated operator and test
  `||K V - RHS||_2 <= AbsResidualTol + RelResidualTol * ||RHS||_2`.
  Absolute numerator, RHS denominator, threshold, near-zero-RHS flag and pass/fail
  are recorded. Relative residual is null for an exactly zero denominator.
  Failed checks save their measured values then abort; later samples stay uncomputed.

The PDE KSP/preconditioner is **not constructed**. H1 operator/space hierarchy is
retained for residual reconstruction, and ND/RT operators are retained for flux
recovery. Only one excitation and one constituent are held at a time; full135
archived potentials are not loaded simultaneously.

`GradFluxErrorEstimator` recovers D from E using the existing global L2 projection
onto RT. Its unnormalized indicator is

`eta_raw^2 = integral |sqrt(epsilon) E - invsqrt(epsilon) D|^2 dOmega`.

The report gives `EtaRawSqrtJ`, domain energy `EnergyJ = U`, and
`NormalizationSqrtJ = sqrt(2 U)`, with normalized `Eta = eta_raw / sqrt(2 U)`.
This is exactly the existing estimator and energy normalization, not a replacement
quadrature/physics formula. The domain and configured MA/MS/SA whole/window energies
come from `PostOperator::GetElectrostaticEnergies`. For a coefficient vector c,
compare domain to `c^T Q_domain c`, each whole interface to `c^T Q_total c`, and
inside-window energy to `c^T Q_window c`. Off-diagonal symmetric CSV entries count
twice. Window energy is whole minus outside, as in existing postprocessing.
`InterfaceResponses` additionally reuses the unchanged single-basis response-matrix
integrator: `InsideJ`, `InsideNormalJ`, `InsideTangentialJ`, `TotalJ`, `TotalNormalJ`,
and `TotalTangentialJ` for each configured localized interface/window. This provides
all18 coupon surface quantities in addition to the domain energy. It requires the
same `LocalizeEdgeEnergy` and `EdgeDistances` configuration as the response matrices;
nonlocalized interfaces retain only the existing `Interfaces` totals, not invented
zero polarization results. `SurfaceUsesRecoveredFlux` records whether any configured
interface uses recovered flux. No surface observable or quadrature is switched.
For configurations using recovered D in surface postprocessing, these are freshly
recovered-D energies; do not expect equality to old low-tolerance D-archive energies
without separately assessing recovery convergence.

Recovery uses config `Solver.Linear.EstimatorTol`, `EstimatorMaxIts`, and
`EstimatorMG`. Actual per-excitation convergence, iteration counts, initial/final
**preconditioned CG** residual norms and their ratio are recorded along with all
controls, including the existing machine-epsilon absolute tolerance. Repeat with
explicitly tighter controls to assess stability. `0.5/5` is not an accuracy
certificate. An unconverged recovery is marked `recovery_unconverged`, not hidden or
silently retried. A completed report means all requested diagnostics were evaluated,
not that every recovery converged or any physical observable was certified.

Zero/small energies have null normalized eta but retain measured raw eta and energy.
`CoefficientCancellationRatio = ||sum c_i V_i||_2 / sum |c_i| ||V_i||_2` flags
cancellation separately from the energy floor. `NormalizationUsable` means only that
energy/cancellation floors and the configured recovery convergence test passed,
**not** reliable physics accuracy. Each excitation has `SampleCount: 1`; the top-level
count is the number actually completed. The legacy AMR console summary receives an
unnormalized aggregate; use the per-excitation JSON, not that summary, for assessment.

## Optional fixed-field surface quadrature sensitivity

A request may include `"SurfaceQuadratureExtras": [0, 4, 8, 12]` (sorted unique
integers in0–12, beginning at0). This re-evaluates only the existing interface-matrix
quadrature on the **same fixed E/D vectors**; it does not reassemble K, change BC,
recover flux again, or mutate global `DefaultIntegrationOrder` settings. Default0
preserves all ordinary calls. Ownership-rule overrides are explicitly unsupported.

Each excitation then has `SurfaceQuadratureSweeps` with all18 interface/window
components at each extra order, and `SurfaceQuadratureControls` records the trial
order and unchanged base/Jacobian-order settings. Rank-local
`archive-quadrature-case-NNNN-rank-NNNNNN.json` files record each interface/geometry,
requested and actual rule orders, point count, face count and minimum reference
weight. These files are needed because actual MFEM rules may serve several requested
orders. Nonpositive/nonfinite rule weights and invalid diagonal energies fail closed;
no higher-order triangle rule is silently replaced. All input/output hashes remain
in launcher provenance.

These are **quadrature sensitivity measurements**, not guaranteed convergence or
exact window cuts. The synthetic constant-normal-field test independently splits a
unit-square perimeter window to exact area0.64 and also checks each sampled GL rule;
it demonstrates that the hard window can remain inaccurate/nonmonotone even when
whole-surface polynomial energy is exact. Baseline0, domain/BC/residual/recovery and
subsequent excitations are regression checked against the no-sweep path. A frozen
baseline binary can additionally be bound with `ARCHIVE_DIAGNOSTIC_BASELINE_EXE`.

## Outputs and provenance

- `postpro/archive-estimates.json`: request, units, global FES counts, controls,
  checks, per-excitation indicators and energies, sample count and completion status.
- `postpro/archive-layout-rank-NNNNNN.json`: local H1/ND/RT widths, element count and
  native-endian FNV-1a64 fingerprints of local geometry/ordered FE dof maps. These are
  diagnostic fingerprints of the **reconstructed** spaces, not historical ordering
  hashes and not cryptographic proof of archive interpretation.
- Optional `archive-elements-case-NNNN-rank-NNNNNN.csv`: local element index,
  attribute, element-center coordinates in m, and raw squared indicator in J. Its sum across ranks equals raw eta
  squared. No replicated mesh is written. Element indices refer to the reconstructed
  rank-local mesh; keep mesh, partition configuration and layout records to localize.
- Launcher provenance hashes original config, request, basis contract, mesh, every
  source, executable, used archive payloads, run config, and all output records.
  `solver.log.json` records bounded process-tree RSS/time and exit status.
- Existing `palace.json` and resolved config provide normal executable/run metadata.
  Binary hash/source-build provenance takes precedence over a stale generated git
  banner from an existing dependency build.

Headers and hashes of reconstructed layouts cannot prove old DOF ordering. The
independent BC, assembled residual and old quadratic-energy checks are all necessary.
Even when they pass, the flux indicator is a **volume energy-norm heuristic**, not
an ideal-boundary or interface-observable error bound.

## Local regression

```sh
export ARCHIVE_DIAGNOSTIC_EXE=/absolute/new-binary
export ARCHIVE_DIAGNOSTIC_SCRATCH=/absolute/external-scratch
export TMPDIR="$ARCHIVE_DIAGNOSTIC_SCRATCH" PYTHONDONTWRITEBYTECODE=1
python3 examples/cpw3d_surface/spatial_coupon/run_bounded_mesher.py \
  --seconds 600 --memory-gib 6 --log "$ARCHIVE_DIAGNOSTIC_SCRATCH/tests-new.log" -- \
  python3 examples/cpw3d_surface/spatial_coupon/test_estimate_archived_fields.py -v
```

The tests preserve synthetic evidence in unique scratch directories. They use only
8 hexes/p4 (729 H1 DOFs), 1/2 ranks, analytic affine/zero/kink fields, explicit
combinations and cancellation, ordinary/archive agreement, domain and whole/window
MA/MS/SA quadratic checks (all19 quantities, nontrivial polarizations, raw and
recovered configurations, zero/tiny fields), localization sums, recovery tolerance
sweeps, corruption (including interior-only permutation with exactly unchanged BC
and ground-only corruption at a nonzero intended-trace intersection), separate
matching/ground/intersection counts on1/2 ranks, export sentinel preservation,
unsupported optional input rejection, unused-source mapping swap rejection,
flag/constraint/output-reuse rejection, input nonmutation and MPI layout agreement.

This implementation does not extend the separate p4 trace-projection helper's
quadrature assumptions to p5/p6. That requires a separately tested assessment-stage
change. No real135-source indicator or h/p convergence result is implied by these
synthetic tests. Review and a concrete resource/probe approval are required before
any production-sized archive diagnostic.
