<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->

# Changelog

> Note: *Palace* is under active initial development, pre-v1.0. Functionality and interfaces
> may change rapidly as development progresses.

The format of this changelog is based on
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to
[Semantic Versioning](https://semver.org/). Configuration files follow
[Schema Versioning](https://docs.snowplow.io/docs/pipeline-components-and-non-pipeline-components/iglu/common-architecture/schemaver/).

See the [developer notes on schema versioning](https://awslabs.github.io/palace/dev/developer/notes/#Schema-versioning) for how versions are bumped.

## In progress

#### Interface Changes

  - THE WORDS OF `surface-participation.csv` TELL THE TRUTH: each row reads ONE SUMMARY
    WORD in a new column `summary` (after `p corrected (fixed-trace)`), in this order:
    `NOT-TARGET`; `NOT-COMPUTED` in `"PreflightOnly"`; `DO-NOT-USE` when the value is raw
    where a correction is due (metal law or target not served), `UNDER-RESOLVED` or
    `NOT-EVALUATED`, not a number, or marked by `"MarkNotApplied"`; `OK` for an `SA` or `MS`
    value that is `RESOLVED` (or `RESOLVED-UNTESTED`) with coverage `OK`, no failed Maxwell
    gate and no tolerated library mismatch, and not a box (`region`) value while the bars'
    region scope is untested; `CHECK` otherwise (`NO-BARS`, every `MA` value:
    `MA not certified`; every box value under the built-in bars: no box value has read `OK`
    on the benchmark, and a box's census reads its own patches, optimistic on junction
    boxes). On the benchmark's 188 graded `SA` / `MS` cells: `OK` 98 within 5 %
    of thick and 1 beyond, `CHECK` 27 / 12, `DO-NOT-USE` 6 / 44; electrostatic and driven
    `Type` cells read `OK`, no eigenmode or box cell yet; a 5 % word only (44 of 99 within
    2 %), no chip.
    The columns `regime` and `confidence` are REMOVED (their gates failed on every served row
    and did not rank the error; the developer tables keep them unchanged). `flags` holds bad
    words only, now with the failed Maxwell gates of the value (`kR <v> above 0.1`,
    `loop residual <v> above 0.05`, `loop fraction <v> above 0.01`,
    `library distance <v> above 0.8`, `boundary law not verified`; the log names them in a
    Warning); a new column
    `information` (after `flags`) holds `Ls sheet served as PEC`, the two `S` words,
    `R/rho <v> above 0.25` and `MA not certified`. The trace-closure Warning of every solve
    becomes one information line per run (the gap between two closures, information, not a
    verdict). `palace.json` `SurfaceResponse.Coverage.Flagged` no longer lists a cell whose
    only flag is `NO-RANGE`, and the end-of-run `NO-RANGE` count prints only when a range
    file is loaded. `"Coverage": {"ValidatedRegime"}` is a parked developer option (refused
    by name without `"Developer": true`; no range file is shipped). No value, solve,
    marking or library changes. SchemaVer 8-0-0
  - THE CHECK OF A CONFIGURATION AGAINST ITS RESPONSE LIBRARY RUNS ONCE, AT LOAD, before
    the mesh is read, in the solve form and in `"PreflightOnly"`, and prints ONE message
    with every mismatch (quantity, the entry by its `Index`, `Type` and attributes, both
    values, the relative mismatch, its estimated effect, the tolerance, the action): a
    designer's chip learns a permittivity slip in seconds where it learned one entry per
    run after the mesh load (33 minutes on a designers' flip chip). ONE RULE: a mismatch is
    tolerated while its estimated effect on a corrected value stays within 1 % (factors
    times the relative mismatch: 1 for a layer's permittivity and thickness, the substrate's
    measured on two designs' windows; per target interface its layer mismatches add and the
    largest substrate one is added; reference.md gives the factors and where they were
    measured); inside, one warning; beyond, the run stops, with NO OVERRIDE
    (`"Developer": true` does not relax it).
    CHANGED ACTIONS: a layer `Thickness` or `Permittivity` that differs from the library's
    `InterfaceLayers` by less than the rule allows (1 %; the designers' 11.47 against
    a library's 11.45: 0.175 %) WARNS where it stopped the run; a substrate permittivity next
    to the corrected metal beyond the rule (about 0.5 % alone) STOPS the run after the mesh
    where it only warned, and at load when no configured material other than vacuum could
    pass beside the layer mismatches. NEW: `Fabrication.LengthUnit`, when the library states
    it (`um`, `µm`, `μm`, `micron`, `nm`, `mm`, `cm`, `m`), must name `Model.L0`.
    Unchanged in action, now at load and all listed: `MatchingRadius` against the largest
    `EdgeDistances`, a missing `InterfaceLayers` entry, the library file and `Version`, a
    parked metal law without the developer switch, `TargetInterfaces`.
    Every configured material that would stop the run after the mesh if it touched the
    corrected metal is NAMED at load (with its value and mismatch), and the load verdict then
    reads "TOLERATED unless" it touches the metal. An effect of exactly 1 % is within the
    rule. THE AIR MUST BE EXACTLY VACUUM (Palace's own test, as for the vacuum side of `MA`
    and `SA` interfaces, whose raw values are otherwise zero): with an `MA` or `SA` target
    and no vacuum material the run stops at load; a material nearer to 1 than to the
    library's substrate that is not vacuum (air entered as 1.0006) is never judged as a
    substrate slip: named at load, it stops the run after the mesh if it touches the metal.
    With any corrected target (`MA` included: a metal edge with no substrate is not what the
    library was built for) and no material of the substrate class (not vacuum and nearer to
    the library's substrate than to 1) the run stops at load, naming each material
    that is neither with its value: a substrate of 3.9 against a library made for 11.45, or
    one entered as vacuum, stops in seconds (reference.md: the decision table at load).
    A tolerated mismatch is named in the outputs: the flags word `library mismatch tolerated` on every affected row of `surface-participation.csv` and
    `SurfaceResponse.LibraryMismatch` in `palace.json`; without a mismatch nothing is added.
  - The length-weighted census (the PREFLIGHT line, `SurfaceResponse.Resolution.Preflight`:
    `EdgeLength`, `LengthM90`, the shares, the coarse-cells flag and the word) weighs each
    patch's length and knots by its model weight (verified positive): an interpolated
    feature's two or three members counted its edge two or three times (a designers' flip
    chip: 3.927 m of census length against 2.776
    m once per edge). Patch counts are unchanged; the energy-weighted census already
    weighed them.
  - ONE resolution verdict: the PREFLIGHT word (the log's PREFLIGHT line before each
    level's solve, `Status` of `SurfaceResponse.Resolution.Preflight` in `palace.json`, the
    `resolution` of a `"PreflightOnly"` run's `surface-participation.csv`) is the
    per-iteration table's word, by its rule and its verdict data, on the census weighted
    by length (edge patches by cell length, spatial patches as edge patches with as many
    knots), with no convergence evaluated: `MARGINAL`, `UNDER-RESOLVED`, `NO-BARS` or
    `NOT-EVALUATED`, never `RESOLVED`, as the table prints it at a first iteration. The
    share-rule words of the census bars are retired from the PREFLIGHT (`EdgeStatus` and
    `SpatialStatus` are removed; the shares stay as numbers); `LengthM90` (``h/R``) is
    added; a `"PreflightOnly"` row reads the bare word (no `PREFLIGHT ` label). The words
    of a solved level and the energy-weighted census are unchanged. The word can read
    better than the table's later word (a short, coarse, energy-carrying stretch such as a
    junction lead); hence a FLAG: where a cell of a row's census is coarser than twice
    the bar, the PREFLIGHT line gives the share of the length beyond it and the largest
    ``h/R``, `palace.json` `Preflight[]` records `UnderResolvedLengthShare` and `LargestM`
    on every row, and the `"PreflightOnly"` rows of that `Type` add the flags word
    `cells beyond twice the bar: may read UNDER-RESOLVED after the solve`. Without the
    flag no cell is beyond twice the bar (``h/R`` read by the census's bin centre, as the
    table's ``m_{90}``), so the table's energy-weighted ``m_{90}`` cannot exceed it on
    that mesh whatever the field: a pre-solve `MARGINAL` without the flag cannot become
    `UNDER-RESOLVED` after the solve.
  - `surface-participation.csv` gains one column, `S corner+cluster share (fraction)`: the
    share of a row's within-R corrected energy served by corner and spatial-cluster models
    (whose trace map onto 16-knot rings is an open item); `flags` adds
    `S >= 5% on corner and cluster models` at 5 % or more and `S over 10% on cluster models`
    above 10 % on cluster models. No value changes.
  - Purely inductive `Impedance` sheets are served by the perfect-conductor library models
    by default: `"ImpedanceMetalRule"` of `Solver.SurfaceResponseCorrection` defaults to
    `"Apply"` (was `"Off"`; the electrostatic block's default stays `"Off"`). A sheet with
    `Rs` = `Cs` = 0 and a configured `Ls` at or below `"ImpedanceMetalMaxInductance"`
    (default 9e-14 H/sq, the validated bound; validated at the designers' 4.5e-14 H/sq in
    driven and eigenmode runs, cracked and uncracked) is corrected; the run prints this
    regime in one line (it replaces the line "experimental and not validated") and
    `surface-participation.csv` flags the rows `Ls sheet served as PEC`. On a run with an
    `Impedance` boundary, a sheet above the bound, with `Rs` or `Cs`, or a metal of another
    law is not served by the perfect-conductor models: its values stay raw unless the
    library holds a model of its own law, one line names it with the bound, its rows are
    flagged `metal law not served: raw there`, and `"UnmatchedPolicy"` decides as before. A
    run without an `Impedance` boundary runs exactly as `"Off"`, also under an explicit
    `"Apply"` (which printed its lines and wrote its record there before): a metal of
    another law there (a `Conductivity` metal alone) is named only by the law-caused
    warning, as under `"Off"`. `"Record"`,
    `"ImpedanceMetalFilmOverrides"` and an `"ImpedanceMetalMaxInductance"` above 9e-14
    H/sq need `"Developer": true`. The
    bound compares the configured `Ls` (H/sq as written), inclusive: a sheet at exactly the
    bound was refused for about 6 % of characteristic lengths (the scaled `Ls` times its
    scale landed one unit in the last place above it), and the recorded served `Ls` values
    are the configured numbers. SchemaVer 7-0-0
  - Synchronized with *Palace* 0.20.0. The surface-response outputs follow its time-averaged
    energy convention ([PR 937](https://github.com/awslabs/palace/pull/937)): an energy of a
    complex (frequency-domain) field is the time average, one half of the Hermitian form,
    applied once by the owner of each energy, as for *Palace*'s own energies. What changes
    for a designer, by solver:
      + eigenmode and boundary mode: every energy in J of the surface-response outputs
        halves, exactly when the field is unchanged, as `domain-E.csv`: `E_norm`, `E_surf`
        and `W within-R raw` of `surface-Q-corrected.csv`; `E_out`, `E_ann` and the local
        energies of `surface-Q-edge.csv` and `surface-Q-edge-local.csv`; the "(J)" columns
        of the `surface-response-*-energy.csv` tables (model, model-raw, patch, exclusion,
        uncovered and its by-key and by-feature forms, loop-remainder, conductor-consistency,
        domain-boundary, spatial-cluster-priority), of the near-edge energy budget, of the
        coverage, resolution, loop-fit (the loop weight and the withheld domain correction),
        magnetic and impedance-candidate tables; and the "(J)" values of `SurfaceResponse` in
        `palace.json`.
      + port-driven (lumped or wave ports): the same energies do not change, up to the solver
        tolerance (the ports deliver 1 W of time-averaged incident power, so the fields scale
        by sqrt(2) and the energies by 2 x 1/2); the linear quantities scale by sqrt(2): the
        device traces and loop scales (V), the Ls-law gate normalizations and leg values (V),
        and the conductor-consistency `Amplitude` and `Normalization`.
      + driven by a surface current or through Floquet ports: the energies halve.
      + electrostatic: nothing changes.
      + the participations (`p_surf`, `p_out`, ...), Q, shares, fractions, correction
        factors, residuals, statuses and verdicts do not change (ratios of like energies),
        nor does the fabrication-process library (static responses, stored unchanged).
  - A frequency-dependent material `"Permittivity"` (*Palace* 0.20.0) is refused with the
    surface-response correction and with the interface keys it builds on
    (`"EdgeDistances"`, `"LocalizeEdgeEnergy"` and `"SaveLocalEdgeEnergy"`,
    `"FluxRecovery"`, the ownership keys and `"Region"`), for every problem type, with a
    message naming the key.
  - The response-corrected adaptive driven sweep (the parked self-consistent path, with
    `"Developer": true` and `"CorrectionMode"` `"SelfConsistent"` or `"Both"`) follows
    *Palace*'s certification of the PROM and its sample refinement (its samples are solved
    from zero, as *Palace*'s, not from the previous corrected field: values move within the
    corrected solve's `"SolveTol"`); every output file of the
    surface-response product is written through *Palace*'s checked writers (a failed
    write aborts).

#### New Features

  - Surface-response-corrected runs write `surface-participation.csv`: per adaptive level
    the raw and corrected participations with their verdicts and a convergence indicator.
    The optional `"Verdict"` member of the `Resolution.BarsFile` file replaces the built-in
    verdict data. SchemaVer 6-1-1
  - A `"PreflightOnly"` run of a surface-response-corrected electrostatic, eigenmode or
    driven configuration evaluates the resolution census of the first mesh without a solve
    (no system matrix is assembled and nothing is solved; a wave port's boundary-mode
    matrices are built with the space operator): the knots and PREFLIGHT lines,
    `SurfaceResponse.Resolution` of `palace.json` and, with `"PatchEnergy"`,
    `surface-response-resolution-patches.csv`, as a corrected run writes them before its
    first solve, and one row set of `surface-participation.csv` whose resolution reads
    the table's word before a solve and whose values are `nan` (not computed); a
    preflight that fails closed by its default `GeometryFailurePolicy` prints that
    summary and then builds the census, as a solve-form run does, before its exit
    status 3. Where a solve-form run stops before its census
    (`"UnmatchedPolicy": "Error"` with missing requirements, a library its reader refuses, a
    response model whose fabricated domain matrix its matrix step refuses, a library file
    its model loop refuses: a basis-point, domain-matrix or surface-matrix file that cannot
    be read, holds a cell the table parser refuses (text, a binary file, a cell beyond the
    header's columns), does not fit its basis points, or lacks a mapped coupon's response
    (a header-only surface file), blend sources without the same coupon interfaces, or
    interface mappings without surface
    matrices, an explicit `"GeometryFailurePolicy": "FailClosed"` with a fallback) no
    census is built,
    one line says why (for library files every refused model once, with every refused
    file's message, the first being the solve form's; also kept as
    `SurfaceResponse.ResolutionCensusNotEvaluated` of
    `palace.json`), and the exit status is the one the preflight had before. The model
    files are read by non-aborting forms of their readers, which the solve form's readers
    call and stop with the same message; a library CSV file is read in one pass by the
    table parser's own rules without stopping (`TryParseTable`, the same table bit for
    bit), since a build without MFEM exceptions (the containers) cannot catch the parser's
    stop. A solve form stopped by a non-numeric library cell now prints a message that
    names the file (`Response matrix file "<path>": Could not parse CSV table from
    string! ...`, where it printed `Could not parse CSV table from string!`); its exit
    status is unchanged. Other corrected solve runs are unchanged. The census
    reuses the manifest pass's identification: on a designers' chip (one node of 192
    cores) it was printed about 114 minutes after the start (measured
    before the matrix-step prediction, which adds one eigendecomposition per model on the
    root; not re-measured), against 72 to 76 minutes for a solve-form run to reach the same
    census.
  - Added the configuration keys of the surface-response correction and of the features
    below to the schema, all optional (a configuration accepted before is accepted with the
    same meaning): `Solver.SurfaceResponseCorrection`; `Solver.Electrostatic`
    `"ResponseMatrix"`, `"AggregateResponseMatrix"`, `"ResponseArchive"` and
    `"ResponseCorrection"`; the edge and ownership keys of
    `Boundaries.Postprocessing.Dielectric[]` (`"EdgeAttributes"`, `"AutomaticEdges"`,
    `"EdgeExcludeAttributes"`, `"EdgeExcludeSegments"`, `"EdgeExcludeSegmentTolerance"`,
    `"EdgeExcludeSegmentSagittas"`, `"EdgeDistances"`, `"EdgeDistanceSmoothing"`,
    `"LocalizeEdgeEnergy"`, `"SaveLocalEdgeEnergy"`, `"OwnershipDataFile"`,
    `"OwnershipGroup"`, `"OwnershipSlot"`, `"OwnershipQuadratureOrder"`,
    `"EdgeFrameNormal"`, `"EdgeRefinement"`, `"FluxRecovery"`, `"Region"`);
    `Model.Axisymmetric`; `Model.Refinement.MaxSizePredicted` and `"MinElementSize"`;
    `Boundaries.PrescribedPotential`; `"Length"` of a lumped port and of its
    `"Elements"`; and the impedance-metal rule's `"ImpedanceMetalRule"` (`"Off"`,
    `"Record"`, `"Apply"`), `"ImpedanceMetalFilmOverrides"` and
    `"ImpedanceMetalMaxInductance"` of `Solver.SurfaceResponseCorrection` (EXPERIMENTAL);
    and `"ConductorConsistency": {"Exclusion"}` (`"PerExcitation"`, `"Sticky"`) of the
    surface-response correction, whose default `"PerExcitation"` CHANGES the meaning of a
    configuration whose run has a conductor-consistency exclusion before a later excitation
    (that excitation's correction no longer omits the patch unless its own test fails; set
    `"Sticky"` for the former rule). SchemaVer 6-1-0
  - Added the edge-resolution census of the surface-response correction (lane
    rescensus-product; outputs only, ON by default whenever the correction is on:
    `"Resolution": {"Enabled", "BarsFile"}`): one added point location of every applied
    patch's contour knots per AMR level (the construction's own locator, no fallback
    search); per patch the lower median of the owning elements' longest edge over the
    matching radius (a knot on a face shared by several elements: the largest of their
    longest edges; on a given mesh independent of the rank count);
    `surface-response-resolution.csv` (per cell, interface `Type` and plane: the
    energy-weighted shares of RESOLVED / MARGINAL / UNDER-RESOLVED /
    NOT-EVALUATED patches, the `m10` / `m50` / `m90` quantiles, a status on INTERIM bars
    valid at `R` = 1.9 µm only, NO-BARS otherwise), the `palace.json` block
    `SurfaceResponse.Resolution`, a PREFLIGHT line before each solve, a banner line per
    `Type`, and with `PatchEnergy` `surface-response-resolution-patches.csv`. No corrected
    value or coverage status changes. RESOLVED is necessary, not sufficient.
  - Added the coverage outputs of the surface-response correction, on every corrected run of
    both solve families: the
    closed partition of the raw within-`R` energy `W` by the correction's ownership, one
    quadrature pass per target interface; `surface-response-coverage.csv` (per cell and
    interface `Type`: `W`, the uncovered / missing / method-gap / pending-model /
    no-status-model shares of `W`, the covered correction factor, its validated range, an
    INTERIM coverage status and a regime status, the untargeted raw energy of the `Type` and
    the loop-remainder terms), `surface-response-uncovered-energy-by-key.csv`
    (per Missing-class key) and `terminal-C-corrected-coverage.csv`; the coverage columns
    APPENDED to `surface-Q-corrected.csv` (both families) and
    `surface-response-confidence.csv`; `SurfaceResponse.Coverage` in `palace.json` and an
    end-of-run banner per cell that is not `OK`; the config block `Coverage` (`Policy`
    `"Warn"` only, `MaxUncoveredShare` 0.05 and `BuildTarget` 0.02 fixed interim bars, an
    optional `ValidatedRegime` sidecar checked against the library's SHA-256). Every existing column keeps its
    position, name and value. The geometry cache is version 25 (the uncovered features' key
    records; an older cache is refused by name).
    The coverage table also prints the signed method-gap share before its clamp at 0; a
    cell whose partition serves more than the correction removes reads `INCONSISTENT`
    (never `OK`); a served claim end one round-off short of its edge segment's end is
    snapped to it in the coverage partition only (every corrected value unchanged), and
    `SurfaceResponse.Coverage.ClaimEndCensus` counts such ends (served and kept) by
    direction: round-off short with the exact vertex energy concerned, round-off beyond
    (the end past the vertex: nothing at stake), near; a
    `ValidatedRegime` range requires `N` >= 1 and a `Provenance`; a served model without a
    qualification record counts as no-status; the log prints, per interface `Type`, the
    Missing share of the target perimeter LENGTH before the solve. A configured target
    interface on which no loaded model has a fabricated surface is a coverage interface
    (its within-`R` raw energy uncovered in its `Type`'s `W`) with the status
    `TARGET-NOT-SERVED`, a startup warning and `SurfaceResponse.Coverage.
    ConfiguredTargetsNotServed`; a run with no coverage cell writes `NoCellsReason` and
    `"Complete": true` only with that reason. The coverage block names the AMR `Level` that
    wrote it of `Levels` (`Refinement.MaxIts` + 1) and gets `"Final": true` only after the
    AMR loop has ended, so the block of a run killed during a later level is not read as
    final. A cell whose coverage is `OK` but whose regime is `NO-RANGE` is flagged per cell
    (banner and `SurfaceResponse.Coverage.Flagged`) as an `OK` + `OUT` cell is, with the
    words "computed, NOT VALIDATED in this regime". A `ValidatedRegime` file with a
    `Driven` range is refused by name (the class lumps every driven excitation).

  - Added `Coverage.Policy` `"MarkNotApplied"` to the surface-response correction (both solve
    families; `"Warn"` stays the default and prints every value as before): a corrected value
    built from a cell whose coverage status is `PARTIAL`, `NOT-ACHIEVABLE`, `INCONSISTENT`,
    `TARGET-NOT-SERVED` or `NOT-MEASURED` prints `nan` and its value moves to an appended
    `<column> partial (uncovered part raw)` column; the corrected normalizations, frequency
    shifts, self-consistent frequencies and terminal matrices built from such a cell are
    marked the same way (palace.json nulls + `FixedTracePartial` / `SelfConsistentPartial`;
    `MarkedCells` / `MarkedMeasurements` / `Flagged[].Marked`; banner lines); every corrected
    table states each cell's mark explicitly in integer `coverage mark ...` columns appended
    last (`1` = marked; the headers are listed in the reference). No computed value changes.

  - Added `Diagnostics.UncoveredFeatureEnergy` to the Maxwell `SurfaceResponseCorrection`
    (off by default; refused by the electrostatic driver): the raw within-`R` energy of the
    uncovered requirements (`surface-response-uncovered-energy.csv`) split per type by
    identification FEATURE into `surface-response-uncovered-feature-energy.csv` (the label's
    portions, length and bounding box, the energy per target interface; overlapping features
    of one type share a `+`-joined label), one quadrature pass per interface; per mode, type
    and interface the rows sum to the uncovered table (lane ctx002-miss, decision 1567 (b)).

  - Added `Solver.Electrostatic.ResponseArchive` (`Mode` `"None"` / `"Archive"` /
    `"Reduce"` / `"Estimate"`, `Directory`, `BlockSize`, `SourceTiming`,
    `RecycleInitialGuess`, `EstimateRequest`): the archived coupon-response generation
    modes as configuration keys, so that a run launched with a configuration file only (a
    container job) can select them. The
    `PALACE_RESPONSE_ARCHIVE_ONLY` / `_REDUCE_ONLY` / `_ESTIMATE_ONLY` / `_ARCHIVE_DIR` /
    `_BLOCK_SIZE` / `_SOURCE_TIMING` / `_RECYCLE_INITIAL_GUESS` / `_ESTIMATE_REQUEST`
    variables stay as aliases with their meaning unchanged; a key and its variable that are
    both set must agree. Results are unchanged.
  - Added `Diagnostics.Directory` and `Diagnostics.ConductorConsistencyRecords` to the
    surface-response correction (`Solver.Electrostatic.ResponseCorrection` and
    `Solver.SurfaceResponseCorrection`; decision 1685): with a directory name set, the
    per-patch tables `surface-response-patches.csv`, `surface-response-patch-energy.csv`,
    `surface-response-traces.csv` and `surface-response-patch-loops.csv` are written into
    that sub-directory of `Problem.Output` (the same names and bytes; each AMR
    `iterationNN/` gets its own), and `"ConductorConsistencyRecords": "Directory"` moves
    the `Records` array of `SurfaceResponse.Diagnostics.ConductorConsistency` out of
    `palace.json` into `surface-response-conductor-consistency-records.json` there (the
    same array text; `palace.json` lacks only that member). Defaults: absent and
    `"Metadata"`, every output as before.

  - Added the magnetic (inductive) defect of finite-thickness metal to the Maxwell
    self-consistent surface-response correction for the 2D translational model families
    (mx-sc-magnetic DESIGN-v2 option (C); USER decision 1221 D2-D8; decisions 1227 / 1229):
    OPTIONAL Version-3 library leaves `MagneticFabricatedMatrix` / `MagneticThinMatrix` /
    `MagneticSource` (the model's electrostatic coupon pair re-run at ``eps = 1``, the exact
    2D duality ``S_m = Q_dom,vac``; assembled by `examples/cpw2d/attach_magnetic_leaves.py`)
    with the top-level `MagneticRules`; the K-side defect ``dK = P_m^T W dS_m P_m`` on the
    flux trace of ``curl E`` in every corrected operator (uniform and adaptive driven, the
    self-consistent eigenproblem, its shift-invert and nonlinear refine); the refusal rules
    r1-r4 by name (r4 two-level: FINDING / REFUSAL); the first-order fixed-trace frequency
    shift as the corrected frequency of record (`fixed-trace frequency shift`,
    `f fixed-trace corrected (GHz)`, `magnetic term applied`, `Re{f self-consistent
    (magnetic)} (GHz)` in `surface-response-confidence.csv`), `surface-response-magnetic.csv`,
    the preflight manifest's `Library.Magnetic` / `Summary.MagneticElectricOnly` and the
    coverage report's `MagneticElectricOnly` row; the config keys
    `SurfaceResponseCorrection.Magnetic` (default true) and `AllowElectricOnlyFallback`
    (default false). A library without the leaves is unchanged, bitwise. The spatial
    families stay electric-only (named); the kinetic (London) term is out of scope
    (OI-MAGNETIC-KINETIC-LONDON-TERM). Review fix-ups (decisions 1235 / 1237): the magnetic
    half of the shift is withheld on the electric loop-fit fallback patches (`magnetic
    fallback patches`, `withheld magnetic domain correction (J)`), `magnetic confidence pass`
    gains the loop-response failure fraction (`magnetic loop-response fraction above limit`),
    the magnetic records are written whenever a leaf was read, the manifest keys are
    Maxwell-only, `magnetic term applied` means a non-zero ``dK``, and
    `attach_magnetic_leaves.py` writes relocatable (relative) leaf paths.
    Phase 2 (decision 1242): the r3 same-mesh rule is TWO-TIER: the optional
    `MagneticRules.SameMeshToleranceWithMeshIdentity` applies to a twin whose
    `MagneticSource.MeshOfRecordSha256Match` records both coupon meshes identical to the
    model's meshes of record (absent = `SameMeshTolerance` everywhere).
    Decision 1250 (Phase 2, the (ii) resonator): the fixed-trace shift follows the SAME
    coverage rule as the electric fixed-trace correction (printed on partial coverage, the
    uncovered requirements kept raw): it is no longer withheld by the closure-independent
    confidence gate (whose coverage criterion >= 1 - 1e-10 made it NaN on every device
    with ends); the magnetic loop gate stays; `fixed-trace shift coverage` prints the
    matched energy share the shift rests on. Decision 1251 (INTERIM): the optional library
    rule `MagneticRules.ConductorConsistencyTolerance` sets the MAGNETIC conductor-consistency
    gate's tolerance (absent = the electric gate's 0.02, bitwise as before); the probe form is
    OI-MAGNETIC-DUAL-GATE-PROBE-DISCRETISATION (Phase 3).
    Phase 3 (a) + (b) (decisions 1284 / 1292 / 1302 / 1309 / 1323 / 1334): the magnetic
    conductor-consistency gate reads the LIFTED probe form (its legs lifted off the sheet
    and set back from the singular edges into their piece; the library leaves
    `MagneticRules.ProbeLiftOverR` / `ProbeSetbackOverR` / `ProbeSetbackWidthFraction` /
    `ProbeSheetNudgeOverR` / `ProbeBandHeights`; FLAT = 0 / 0 = the electric legs, the
    comparison form of record); `MagneticRules.ConductorConsistencyTolerance` may be
    per-order (`{"3": t3, "4": t4, ...}`, the solve's order applies, an absent order is
    refused); the magnetic TERM's chain, loop fit, dual gate, fixed-trace shift and ``dK``
    read the conductor anchors set back `MagneticRules.StateAnchorSetbackOverR` x R into
    their piece as twin lines (R-MAG-17; absent = 0 = the legacy anchors, bitwise), an
    anchor within `MagneticRules.StateAnchorEdgeToleranceOverR` x R (default 0.1) of an edge
    snapped onto it (one rule for the term's state anchors and the probe form's anchor end;
    the tolerance below half the narrowest model gap when either setback is in effect; the
    snapped counts recorded as `Correction.MagneticStateAnchorSnapped`); the DIAGNOSTIC
    CONTROLS `ConductorConsistency.ControlFlipProbePiece` and `ControlStateAnchorSetbackOverR`
    (`RecordOnly` only; never of record); the per-model / per-patch magnetic energies
    (`magnetic domain correction (J)` / `magnetic thin energy (J)` in
    `surface-response-magnetic.csv` and, with `Diagnostics.PatchEnergy`,
    `surface-response-patch-energy.csv`) and `Diagnostics.WithinRVolumeEnergy` (the raw
    field's within-R magnetic / electric energies and volume per patch; a per-model sum
    over-counts a multi-edge model by its edge count, OI-WITHIN-R-MULTI-EDGE-OVERCOUNT); the
    records `Correction.MagneticProbeForm` / `MagneticProbeLines` / `MagneticGateTolerance` /
    `MagneticStateAnchorSetbackOverR` / `MagneticStateAnchorPatches`; the attach tool's
    matching flags (`--probe-lift`, `--probe-setback`, `--probe-width-fraction`,
    `--probe-nudge`, `--probe-band-heights`, `--state-anchor-setback` (default 0.5, the
    form of record; 0 = the legacy form), `--state-anchor-edge-tolerance`). Named open items:
    OI-MAGNETIC-DUAL-GATE-3D-ENCLOSED-FLUX (the dual gate's deviation on a device is the 3D
    field's enclosed flux at junctions / ends, excluded by the calibrated per-order
    tolerance) and OI-MAGNETIC-SC-GMRES-STAGNATION (the self-consistent GMRES with ``dK``
    stagnates on a device; the fixed-trace shift of record is unaffected).
  - Added `"GeometryFailurePolicy"` (`"RecordRaw"` | `"FailClosed"`) to the electrostatic
    `ResponseCorrection` and the Maxwell `SurfaceResponseCorrection` blocks (readiness
    hardening, USER decision 1000 (b)): the former geometry-reachable abort sites of the
    identification / placement (no matched feature, every coupon outside the mesh, no metal
    perimeter, overlapping coupon supports without an owner, a refused placement, an alias
    key mismatch) become recorded raw fallbacks on a solve run (`SurfaceResponse.
    Diagnostics.GeometryFallbacks`, `surface-response-geometry-fallbacks.csv`, a LOW
    CONFIDENCE banner with the raw share per interface type) and a deferred abort after the
    manifest on the preflight; nothing changes when no site fires. Readiness hardening-2
    adds the classes `UnlocatedPatch` (a coupon's metal-edge reference off the mesh beyond
    the boundary tolerance `1e-3 R`; a reference within it is read as located) and
    `InterfaceMappingConflict` (a matched model whose interface slots cannot bind the
    feature's distinct interface target maps; the feature is read as Missing, never mapped
    through slot 0), makes the same-class coupon-support overlap test oriented (the
    frame-aligned boxes of the placed basis points by the separating-axis test; the
    cluster-over-vertex priority keeps the mesh-aligned test), and makes the discovery
    (`discover_surface_response_requirements.py`) record the Missing requirements it can
    give no placeholder (`UnboxableRequirements`: an Unboxable key; `UnbuildableRequirements`:
    a key the coupon generator refuses by name) in the manifest and the closure history with
    a `DISCOVERY FALLBACKS` listing instead of dying. The stage-F fail-safe adds the class
    `LongitudinalAxisOutOfTolerance` (a pair / stack quadrature cell whose placed frame's
    AxisW lies beyond |cos| 0.95 of the sample's own segment, the former fail-closed check of
    `LongitudinalCellOffsets`: the cell is excluded from the correction and kept raw with
    the DomainBoundary portions, one entry per affected feature; the tolerance unchanged).
  - New optional `config["Model"]["Refinement"]["MinElementSize"]`, an element-size floor
    for AMR marking (a length in mesh units; default `0` = off, bitwise the previous
    marking). With a floor, only elements whose volume-based size `h` (the `h` of the mesh
    statistics in the log) is at least twice the floor are eligible for marking, the Dörfler
    threshold is computed over the eligible elements only, and the error estimate and the
    `"Tol"` test keep every element. Each AMR iteration logs the floor, the ineligible count
    with its share of the error and the number of elements below the floor created by the
    conformity closure; AMR stops cleanly when no element is eligible.
  - Added `"CorrectionMode"` for electrostatic and Maxwell surface-response correction,
    with `"PostprocessOnly"`, `"SelfConsistent"`, and backward-compatible `"Both"`
    execution modes.
  - Documented the caveat OI-MAXWELL-SC-ELECTRIC-ONLY in the Maxwell response-correction
    reference: the self-consistent driven / eigenmode correction adds only the electric
    (electrostatic-coupon) defect, no magnetic / inductive defect of the thick metal, so its
    corrected frequencies, self-consistent participations and corrected fields are not
    validated (the ordinary port and S-parameter outputs are the raw thin-metal result; the
    fixed-trace participations of record are much less affected); a fix is under design.
  - Accelerated aggregated electrostatic coupon response matrices with batched surface
    Gram assembly and one collective domain-matrix reduction. Added probe-only and
    concurrent-coupon controls to the process-library planner, plus exhaustive virtual
    closure discovery so a partial library cannot hide nested geometry requirements.
    Multi-slot spatial coupons now partition MA/MS surfaces by interface slot and
    conductor before summing repeated target-interface contributions. Electrostatic runs
    also report additive per-model response energies and model metadata.
  - Added `"PatchEnergy"` to the electrostatic `"ResponseCorrection"` (default `false`):
    `surface-response-patch-energy.csv` with the domain correction and fabricated
    interface energies of every applied patch (model, feature, longitudinal cell, weight,
    origin) per source and evaluation, the per-patch rows of
    `surface-response-model-energy.csv`. Added an optional `"Region"` to interface
    dielectric postprocessing entries: a quadrature-level spatial filter (axis-aligned box
    and / or the translational cells of given segments: along-range between the end cuts,
    in-plane transverse distance), so that one interface can be partitioned into entries
    on the same attributes without re-tagging mesh faces. Without either option every
    output is unchanged.
  - Corner-class placement rules of the automatic surface-response correction: at a matched
    corner whose arms are not perpendicular the second arm's straight cells begin where
    that arm exits the corner coupon's matching square, at `R / max(|cos t|, |sin t|)`
    from the vertex (t the corner angle; unchanged at 90 degrees), recorded under
    `Diagnostics.CornerArmTrim`; and a feature without a library model keeps its raw
    within-`R` surface energy in the corrected interface energies, reported per feature
    type in `surface-response-uncovered-energy.csv` (written only when something is
    uncovered); the part of an uncovered portion inside a matched spatial cluster's
    support box is removed at placement (recorded under
    `Diagnostics.Uncovered.ClippedBySpatialSupport` and in the CSV), and trimmed corners
    whose coupon is not applied are counted under
    `Diagnostics.CornerArmTrim.ExcludedCoupons`. Response-geometry cache version 9.
  - The translational `SurfaceMortar` projects the device trace onto the coupon's own hat
    basis: the constrained metal-band knots where the fabricated metal meets the matching
    contour (listed as `ZeroTraceIndices`, or inserted at library load from the model's
    topology and `Fabrication.MetalThickness` for a library that publishes the free knots
    only; a library whose knots disagree with its constraints is refused) are vertices of
    the basis without a coefficient, so the adjacent hats no longer span the band. Recorded
    under `Diagnostics.ConsistentMortar`; libraries without a derivable band and the
    `Collocated` coupling are unchanged. Response-geometry cache version 10.
  - Rounded corners in the automatic surface-response correction: the corner coupon
    planner's refined trace basis rule (`all-rings-follow-metal`) now applies to rounded
    corners as to sharp ones (the box basis does not depend on `CornerRadius`; the legacy
    rule it replaces failed the held-out self-check on every 90-degree corner), the
    library load checks rounded refined coupons against the rule at their angle, and a
    matched rounded corner's arm cells (or an unmatched neighbour's uncovered portion)
    begin where the arm exits the coupon's matching square — its identification claim ends
    `r / tan(t / 2)` beyond it, a stretch that was modelled by nothing — recorded under
    `Diagnostics.CornerArmExtension` (a stretch with nothing to extend is counted in
    `UnhostedLength` with a warning; a virtual corner on a natural truncation plane exits
    at the mirror arm trim's `s_half`). Sharp corners are unchanged. Response-geometry
    cache version 16 (the stacks own-cell pre-image, the corner-arm extension records and
    the foreign-context lane's `Chain: false` context pieces; a version-15 cache of either
    earlier schema is refused).
  - The surface-response requirements manifest writes `KeyText` beside `Hash` on every
    Signature-bearing record (the version-2 `Requirements[]` records and the
    `Identification.Features[]` entries): the exact text the key is the SHA-256 of (the
    signature with its `Type`, nlohmann's compact dump), so a consumer verifies
    `sha256(KeyText) == Hash` and `parse(KeyText) == Signature` without re-serialising doubles
    (nlohmann's Grisu2 writes some doubles with 17 digits where Python's repr is shorter; the
    Python key checks of the coupon tooling use `KeyText` when present and keep the legacy
    float hash for manifests without it).
  - The mirror-formed cluster REQUIREMENT CONTRACT of the automatic surface-response
    correction (`SURFACE-RESPONSE-IDENTIFICATION.md` (e); impl-B5 CONTRACT.md v3): an Unmerged
    mirror-formed SpatialEdgeCluster / CurvedEdge configuration across a natural truncation
    plane is emitted as a merged feature with its contract (the canonical Frame with its
    Chirality, the real / image split of the serialised portions mapped exactly by geometry,
    the lengths over R; refused by name otherwise) and as a Missing / Exact requirement
    carrying `MirrorFormedContract`; a Signature-keyed model of its key is placed on the
    configuration's REAL half only (weight = the real length fraction, the model's edges on
    the real portions, the image half evaluated by even extension within the coupon's reach,
    the touched real cells owned by the coupon and no longer DomainBoundary raw); a model
    carrying the consumer's `MirrorFormed` ENTRY stamp must agree with the contract (fail
    closed). The patch dry run `surface-response-patches.csv` gains a `DomainWeight` column.
    Response-geometry cache version 17 (the mirror-formed placement of every spatial patch;
    a version-16 cache is refused).
  - The conductor-consistency gate's probes of a spatial cluster patch snap onto the nearest
    device arc chord when within the identified arc's recorded largest chord sagitta (plus
    the arc-fit tolerance); straight edges are unmoved; the snaps are recorded as
    `ProbeSnaps` per patch and excitation and travel with the placed patch; the tolerance
    is unchanged. A placed spatial model whose trace data is not readable at the placement
    (a discovery placeholder: its basis-points file is never written) is placed without
    snaps (`ProbeSnapsComputed` false; the geometry-only preflight and the discovery
    closure proceed) and refused by a solve on every rank alike. Response-geometry cache
    version 18 (the `ProbeSnaps` and `ProbeSnapsComputed` of every spatial patch; a
    version-17 cache is refused). A probe of a mirror-formed coupon's image half, which
    the solve evaluates at its reflection into the domain, is snapped at that reflected
    position (the displacement carried back through the same mirror planes); every other
    probe is read as before. Response-geometry cache version 19 (a version-18 cache, whose
    image-side probes are un-snapped, is refused).
  - Added the METAL-GRADED corner trace basis rule (the concave-MS remedy, supervisor
    decisions 590 / 629 O-9; a decide-step branch, not adopted by this entry): the record key
    `TraceBasis.MetalKnotGrading` grades the metal-interior knots from each metal-arm
    crossing into the metal arc exactly as `FreeKnotGrading` grades the free knots (scaled by
    `min(1, MetalArc / FreeKnotGradingReferenceFreeArcOverR)` on a short metal arc; the same
    reference constant); the generator's `--trace-basis all-rings-follow-metal-metal-graded`
    (MetalInteriorKnots 9, RingSize 20, 220 knots). An absent key is the recorded
    equal-fraction rule, so every published model reads unchanged; the load-time file check
    refuses a graded coupon whose record lacks the key. The planner's default rule is
    unchanged. The identification `Conventions.CornerTraceBasisRule` text is a dated history
    of the known rules asserting no per-model fact.
  - Spatial-cluster contract rule T5 (excluded-corner containment) of the automatic
    surface-response correction: a free corner whose claim window lies partly inside a
    cluster's contract box and partly outside every cluster box (the corner coupon the
    spatial-cluster priority drops, leaving its outside arms uncorrected) widens that box by
    `0.25 R` steps until the corner's `R`-neighbourhood is contained, recorded under
    `SpatialSupport.CornerContainment` (the stage-1 growth history merged into the re-read
    record, `Growth.Reread`); a grown cluster without a library model falls back to the
    ungrown key's model (`CornerContainment.Fallback`, reason `NoMatchedModel`) and the
    grown signature is listed as a Missing-class `FallbackPlaced` requirement
    (`Summary.FallbackPlaced`; none for an `Unboxable` grown signature). The electrostatic
    solver writes the diagnostic-only companion table
    `surface-response-spatial-cluster-priority-clipped-correction.csv` (the dropped corner
    coupons' fabricated / raw / clipped-correction energies; every existing output
    unchanged). Response-geometry cache version 20 (a version-19 cache, which may hold the
    ungrown placement of a cluster this rule grows, is refused).
  - Continuation ownership of the automatic surface-response correction: a translational
    stretch that ENTERS a spatial cluster's contract box without continuing one of its
    claims is now judged per cell on the cell's inside interval against the cluster's
    `Chain: true` context pieces (previously only a stretch lying wholly inside the box was):
    a context-edge cell inside the box is owned (weight 0) and a context-edge cell
    straddling a box face is re-cut at the face exactly as a continuing cell (the inside
    part owned, the outside part kept with the weight scaled by kept / cell). Recorded and
    left untouched: a cell the box lies strictly inside of and every `Chain: false` cell of
    such a stretch. New diagnostics leaf
    `TranslationalStretchesPartiallyInsideSpatialSupport` (written only when it has
    records), `PartialStretch: true` on the `ContinuationOwnership.OwnedCells` entries the
    rule touches and a `ContinuationOwnership.PartialStretch` counter block (only when it
    acted). A placement rule: keys, requirements and the response-geometry cache (version
    20) are unchanged.
  - Self-consistent corrected eigenmode runs write the number of corrected linear solves
    that did not reach `config["Solver"]["Linear"]["Tol"]` within `"MaxIts"` to
    `palace.json` as `SurfaceResponse.CorrectedSolver.Eigenmode.NonConvergedLinearSolves`
    (the total over the whole solve, beside `LinearSolves` and `LinearIterations`; 0 in
    every run whose solves converged), and print one warning at the end of the log when it
    is above 0. The example configurations with an AMR size cap now state
    `"MaxSizePredicted"` beside `"MaxSize"` (checked by
    `python3 examples/test_max_size_predicted.py`), and the transmon refinement example no
    longer pre-refines edge tubes by default (`"EdgeRefinement"` stays an explicit opt-in).

#### Interface Changes

  - The surface-response correction is post-processing only by default (the v1 default
    path): `"CorrectionMode"` defaults to `"PostprocessOnly"` (was `"Both"`, which also
    ran the self-consistent corrected solves), and the options set aside from the default
    path need the new developer switch `"Developer": true` of the response-correction block
    (refused by name without it): `"CorrectionMode"` `"SelfConsistent"` / `"Both"`.
    `"SolveTol"` and `"TranslationalDomainCorrection"` configure that solve only: they are
    read with it and, in a block without it, not read and named in one warning
    (`keys ignored`). With the switch the run prints one
    line `DEVELOPER OPTIONS: ...` and `palace.json` records the list as
    `SurfaceResponse.Developer`. Without the switch the parked columns are not written
    (their values are still computed): the fixed-flux and self-consistent columns of
    `surface-Q-corrected.csv` (Maxwell 8 of 33; electrostatic 6 per interface and 2
    totals), the `corrected` matrices of `terminal-*-corrected.csv`, the self-consistent
    columns of `surface-response-confidence.csv` and
    `SurfaceResponse.TerminalCapacitance.SelfConsistent`. The legacy options need the
    switch too: `"PatchConstruction": "Legacy"`,
    `"ConductorConsistency": {"Exclusion": "Sticky"}`, and
    `config["Model"]["CrackInternalBoundaryElements"] = true` in a
    correction run where it would crack an interior perfect conductor (PEC, Ground,
    Terminal, PrescribedPotential; refused at mesh load, naming the per-sheet `"Crack"`);
    the key stays accepted where the interior metal it cracks is impedance (Ls) sheets.
    The device-run writers of the examples default to `--crack-policy default` (the v1
    form; `--developer` writes the records' `legacy` crack with `"Developer": true`), and
    the preflight writer's `legacy` default writes `"Developer": true` with it. The
    alternatives need the switch: `"TraceCoupling": "SurfaceMortar"` (its quadrature key
    `"MortarOversampling"` is not read without it, and named in the `keys ignored`
    warning) and the hand-placed form (a block without `"Library"`). The
    library's magnetic (kinetic) leaves are read only with the switch (by default with
    it; without it the run is that of `"Magnetic": false`, byte-identical to a library
    without leaves); `"Magnetic": true`, `"AllowElectricOnlyFallback": true` and
    `"ConductorConsistency": {"ControlStateAnchorSetbackOverR"}` need it. So do the
    two-dimensional boundary-mode correction, the conductivity or rational-impedance
    coupons of a library where the device carries that law (refused when the library is
    read; listed as `DEVELOPER OPTIONS (library)` with the switch; a coupon of a law the
    device does not carry is never matched and not refused), `"Crack": true` on a
    `Conductivity` or `RationalImpedance` boundary of a correction run, and
    `"ConductorConsistency": {"ControlFlipProbePiece"}`. And the conductor-consistency
    records leave `palace.json` by default: `"Diagnostics": {"ConductorConsistencyRecords"}`
    defaults to `"Directory"` (was `"Metadata"`), which no longer requires
    `Diagnostics.Directory` (without it the file is written to `conductor-consistency/` of
    the output); `palace.json` keeps the count and the gate summary, and the readers of the
    examples (`spatial_qualification.py`'s gate, `coverage_report.py`) read the Records from
    that file where `palace.json` lacks them. Nothing is removed.
    SchemaVer 4-6-0

  - The config writers of the examples write the v1 form by default: the recorded settings
    less the parked options (`developer_form.py`); `--developer` writes the recorded
    developer form (the recorded parked settings, the default of record
    `"CorrectionMode": "Both"` stated where they relied on it, and `"Developer": true`), and
    each writer prints which form it wrote. No bar or threshold changes: values from the v1
    form are not comparable with records made before this cut until the benchmark's
    collocated twins are read.

  - A non-finite numeric cell of a CSV table written through the table formatter or the
    CsvNumber helper is a bare token (`nan` for any NaN, where the directly written tables
    printed `+nan`; `+inf` / `-inf`, and `inf` / `-inf` in an integer column) instead of
    being padded to the column width, so that a reader that takes only a bare token as
    missing (pandas' default `read_csv`) reads the column as numbers; every finite, text and
    header cell is unchanged (the padded `NULL` of a missing cell stays), and the writers
    that print unsigned, unpadded numbers are as before (they can print `inf` or `-nan`,
    also bare).

  - The Maxwell `surface-Q-corrected.csv` and `surface-response-confidence.csv` column
    `matched edge fraction` is renamed `matched edge length fraction` (a length figure, not
    an energy coverage; values unchanged).

  - With the new optional `config["Model"]["Refinement"]["MaxSizePredicted"]` set to `true`,
    AMR checks `config["Model"]["Refinement"]["MaxSize"]` against the predicted size
    of each refined mesh (the last solved number of degrees of freedom scaled by the
    element-count growth of the refinement step) and stops before the solve when the
    prediction exceeds it; the last solved mesh is the final result and `palace.json`
    records the stop under `"AdaptiveMeshRefinement"`/`"MaxSizeStop"`. Previously the
    check used only the last solved size, so the final solve could exceed `"MaxSize"` by
    the growth of one Dörfler step (observed 2-3x); that remains the default rule (off by default, as upstream).

#### Performance Improvements

  - The surface-response geometry passes at chip scale, every output unchanged (lane
    v1prefast; decision 2403 (3)): the corner-arm trim and extension read the patches near
    each corner arm (an AABB tree over the own-edge cells) instead of every patch; the
    identification's T5 corner-containment read skips a corner site whose inside length
    cannot reach the quantum before clipping it against every cluster box, and
    `ClusterSupport` reads the excluded and the truncation segments from two lists built
    once; the preflight's knife-edge census runs its two identifications at R (1 -/+ 0.01)
    on ranks 1 and 2 beside the root's own (on one rank as before). On the designers' flip
    chip (1.6 M patches, 192 ranks) the census comes after 36.5 minutes in `PreflightOnly`
    (114.0 before) and 32.8 minutes in the solve form (74.4), the identification of
    `PreflightOnly` taking 2.3 minutes (54 before); on the CTX-002 chip the `PreflightOnly`
    census comes after 20.4 minutes (61.7). Of what remains on the flip chip, the mesh load
    is 15 minutes.
  - The surface-response correction's per-evaluation memory no longer grows on every rank
    with the GLOBAL patch and point counts (lane fix-response-memory; USER decision 1439,
    decisions 1446 / 1447; results bitwise): the gathered per-patch records of every
    response evaluation (patch contributions, traces, loop shares, path residuals, within-R
    volume energies)
    and the conductor-consistency gate's `Records` live on the root only in the drivers
    (`SurfaceResponseOperator::PatchRecordResidency::ROOT_ONLY`; the default stays
    `REPLICATED` for a caller that reads them on every rank: OI-RESIDENCY-DEFAULT-REPLICATED);
    the per-patch CSV tables are written per excitation (`TableWithCSVFile::
    AppendRowsAndClear`, the former bytes); the raw-claim maps hold a rank's local patches
    only and an excluded patch's claim is fetched from its owner; the placed-patch
    configuration is kept local between AMR constructions and re-assembled exactly before
    the next; the domain-boundary containment test runs in candidate blocks of at most
    `PALACE_RESPONSE_CONTAINMENT_BLOCK_POINTS` points (default 2^20; results-neutral) and
    the raw-claim partition's co-located-interval scan uses an insertion-ordered feature
    index (the same first match). Still replicated by design: the response SETUP's
    placed-patch list (the identification's `config.patches` and the placement's copy of it,
    O(global patches) on every rank until the construction ends) and the placement's other
    setup transients, the model library (the metadata of every model and the matrices of
    every model a rank's patches reference), the derived / file matrix caches (released
    after the construction when no AMR construction follows), the
    within-R box list of a measurement and the PatchEnergy model-raw-claim partition
    (OI-CONDUCTOR-CONSISTENCY-RECORDS-PALACE-JSON-SIZE stays open: the file's size is
    unchanged).
  - `palace.json` is streamed to the file by the JSON serializer (no in-memory text of the
    whole record) and every metadata write is checked: a short write (a full disk, an I/O
    error) aborts instead of leaving an incomplete file (OI-METADATA-WRITE-UNCHECKED-ENOSPC).
    In the drivers the response operator's final metadata write moves its diagnostics
    record instead of copying it, and the root spills the conductor-consistency `Records`
    of each excitation / mode to a scratch file in the output directory
    (`.palace-conductor-consistency-records.spill`, removed after the write) that is
    spliced into `palace.json` byte for byte; the domain-boundary containment cache's
    signature is kept on the root only; a rank frees the matrices of the models none of its
    patches references; under glibc the allocator is trimmed after the identification and
    after the construction; the PatchEnergy claims are held once per rank (the drivers take
    them from the operator into the partition's list); without AMR (`Refinement.MaxIts` 0)
    the per-process response-matrix caches are released after the construction (lane
    fix-response-memory round 2; decisions 1475 / 1484; results bitwise).
  - The point locator of the surface-response correction (the candidate exchange of the
    Maxwell contour quadrature and of the electrostatic trace points, its rank-local gslib
    fallback) and the Maxwell contour-line exchange run in rounds of at most 4,096 local
    points per rank, so that a rank whose routing boxes cover most of the device (a rank
    owning large elements) receives at most ~0.5 GB of queries per round instead of every
    point of the device; every int MPI count of these exchanges and of the owner plan is
    checked in 64 bit before use (an abort instead of a wrapped count). Results are bitwise
    unchanged (every query of a point travels in one round; the answers, the owner rule and
    the output order are unchanged). In the Maxwell construction the placed-patch lists
    (the shared geometry's list, the placement's list, now moved into the applied list, and
    the applied list), the matrices of the models no local patch references and, without
    AMR, the response-matrix caches are released before the line and point exchanges
    instead of after the construction. `PALACE_RESPONSE_EXCHANGE_SUMMARY=1` prints the
    per-rank receive counts of the candidate exchange and the most-loaded ranks' routing
    boxes (log only); `PALACE_RESPONSE_EXCHANGE_BLOCK_QUERIES` is a test-only override of
    the round size (lane fix-response-memory round 3; decision 1501; results bitwise).

#### Bug Fixes

  - Fixed the libCEED element restrictions of H(curl) spaces of order 2 and higher on meshes
    with pyramids, which left out the DOF transformation of the pyramids' triangular faces.
    Every libCEED operator with such a space was wrong there: the system operators of
    eigenmode, driven, transient and magnetostatic simulations, the discrete gradient and
    curl used in postprocessing (the electric field of electrostatic simulations and every
    surface and interface quantity computed from it), and the p-multigrid prolongation.
    Meshes of tetrahedra, prisms and hexahedra are not affected.
  - The per-type portion tables and the loop-fit tables of the Maxwell correction write
    their header with their first rows whatever print that is (a table first written at a
    later mode or frequency, such as the conductor-consistency table, was appended without
    a header).

  - `Mpi::BroadcastLarge` now bounds every single `MPI_Bcast` it issues IN BYTES, by the
    compile-time constant `mpi::BROADCAST_LARGE_CHUNK_BYTES` = 2^28 bytes (the element count
    per call = that / `sizeof(T)`, at least 1; an explicit chunk is capped at it), instead of
    INT_MAX elements. On Open MPI 5.0.8 + libfabric-aws 2.6.0 the broadcast of the mesh's
    serialized text in chunks of up to INT_MAX `MPI_CHAR` elements (~2^31 bytes) segfaults in
    the transport (one broadcast of 2^28 `MPI_DOUBLE` elements = 2^31 bytes passes; the
    failing size is not located), so the restart from a nonconforming mesh whose serialized
    text exceeds 2^31 bytes crashed in the mesh distribution
    (OI-NC-LOAD-BYTESTRING-ABOVE-2GIB-SEGFAULT). Results are bitwise unchanged (the same
    bytes are delivered in the same order).
  - The conductor-consistency gate of the surface-response correction no longer aborts
    ("excluded patch ... without a raw claim record") on the Maxwell path at more than one
    rank when an excluded patch is held by a rank other than its raw claim's owner (the
    applied position vs the placed index; a regression of the round-1 memory fix): the
    excluded patches are gathered and each owner serves its raw claims (lane fix-gate-claim;
    decision 1533; results bitwise for every run that completed before).
  - The Maxwell surface-response patch-energy CSV no longer rebuilds a response's per-patch
    maximum loop residuals for every row on the root (rows x loop records: a chip-scale
    eigenmode postprocessing stalled there, OI-PATCH-ENERGY-CSV-QUADRATIC-ROOT-WRITER), and
    the conductor-consistency gate no longer scans the whole patch table per excluded patch;
    both are indexed once (decision 1544; the same bytes are written).
  - The surface-response construction holds less replicated data per rank at its setup
    peak (lane setupmem-impl, decision 1546): a rank other than the root builds the model
    matrices only of the models its own patches reference, after the placement (the root
    still builds every model: the records and the log are unchanged), and the placement
    moves the configured patch list instead of copying it (the original values of the
    fields it writes are kept in a side table for the readers of the original list);
    two log-only memory-census instants are added inside the patch build. Results are
    bitwise unchanged (the same functions of the same inputs; only where the originals
    are read from changes). A model-matrix failure on some ranks only (the root, or a
    referencing rank) is exchanged after the model loop and after the referenced build,
    so that every rank fails instead of waiting in the next collective (decision 1553),
    and a placement that throws leaves the shared AMR geometry's patch list as a completed
    placement does (decision 1554).

## [0.20.0] - 2026-11-23

#### New Features

  - Magnetostatic simulations now model domains with a `LondonDepth`: the London term is
    added to the curl-curl operator and the kinetic energy of the supercurrent to the
    extracted inductance, so conductors of any thickness and width can be modeled as
    volumes. SchemaVer 4-4-1
    [PR 1111](https://github.com/awslabs/palace/pull/1111).
  - `validate-config` accepts several configuration files and validates them in one Julia
    session.
    [PR 1180](https://github.com/awslabs/palace/pull/1180).
  - Added `Boundaries.LumpedPort[].ExcitationAmplitude` and `ExcitationPhase` (degrees),
    a complex weight a = A e^{iφ} on the unit-power incident wave of each lumped port, so
    ports sharing an `Excitation` index can be driven with different amplitudes and phases.
    S-parameters of single-port excitations are normalized by the weight; with a phase that
    is not a multiple of 180 degrees `port-V.csv` and `port-I.csv` report the real and
    imaginary parts of the incident voltage and current. A phase that is not a multiple of
    180 degrees is an error for transient simulations. SchemaVer 4-3-0
    [PR 1169](https://github.com/awslabs/palace/pull/1169).
  - Added scalar frequency-dependent `Permittivity` objects, with a `HighFrequencyLimit`
    and additive Drude, Debye, Lorentz, PoleResidue, and native Djordjevic-Sarkar terms.
    Frequency-domain operators, ports, postprocessing, energy/error estimates, adaptive
    sweeps, and circuit synthesis use the material model at the solution frequency.
    Compatible static limits are supported; transient problems are not. With the default
    solver settings, eigenmode problems whose only frequency-dependent terms are Drude,
    Debye, Lorentz, or PoleResidue terms are solved as linear eigenvalue problems, with
    SLEPc and ARPACK; with other frequency-dependent terms, the initial guesses of the
    nonlinear eigensolver take these pole terms exactly. General PoleResidue models remain
    usable with gain warnings; a PoleResidue term with a nonzero pole specifies no stored
    electric energy, so the electric energy outputs of a material with one are omitted.
    SchemaVer 5-1-0
    [PR 877](https://github.com/awslabs/palace/pull/877).
  - Far-field extraction (`Boundaries.Postprocessing.FarField`) on interior surfaces, so
    that far fields can be computed when the outer boundary is not a valid integration
    surface, or away from the truncation. The fields of each face are taken from the side
    facing the center of the bounding box of the surface. SchemaVer 4-2-3
    [PR 1110](https://github.com/awslabs/palace/pull/1110).

#### Interface Changes

  - **Breaking:** Frequency-domain energies and powers are time averages of the peak
    phasors: energies carry a factor 1/4 and powers 1/2. Lumped, wave and Floquet port
    excitations carry 1 W of time-averaged incident power; a 50 Ω lumped port has 10 V and
    0.2 A peak. Impressed current sources are unchanged: a surface current source drives
    1/√Z₀ A ≈ 51.5 mA peak (the `I_inc` of `surface-I.csv`), and a current dipole its
    `"Moment"` in A·m. SchemaVer 6-0-0
    [PR 937](https://github.com/awslabs/palace/pull/937).

    To compare with results from 0.19 or earlier, for the same configuration:

    | Solver / excitation                 | Fields, V/I, probes, `rE`, electric/magnetic flux | Energies and power flux |
    |:----------------------------------- |:------------------------------------------------- |:----------------------- |
    | Driven lumped or wave port          | ×√2                                               | same                    |
    | Driven Floquet or impressed current | same                                              | ×1/2                    |
    | Eigenmode or boundary mode          | same                                              | ×1/2                    |
    | Transient or static                 | same                                              | same                    |

    S-parameters, impedances, eigenfrequencies, Q, κ, EPR, participation ratios and error
    indicators are not rescaled. The wave port mode field `E0_<i>` of the boundary output is
    the 1 W mode, ×√2 in eigenmode simulations too. In a mixed-port sweep, use the row for
    each excitation; a drive that combines current and port sources has no single factor.
    Write results to a new output directory rather than appending to restart CSVs from an
    earlier version: the headers and units match, but the values do not. The factors that
    remain against COMSOL and HFSS (time domain port drives, current sources, the HFSS peak
    electric energy, eigenmode scales, COMSOL's quality factor) are listed under
    [Comparing with other tools](https://awslabs.github.io/palace/dev/reference/#Lumped-ports-and-wave-ports).

  - With adaptive mesh refinement on nonconforming meshes, the SLEPc Krylov-Schur eigenvalue
    solves of later adaptation iterations can start from the eigenvectors of the previous
    one, transferred to the refined mesh, instead of the start vector of
    `config["Solver"]["Eigenmode"]["StartVector"]`, and then need fewer Krylov restarts. The
    eigenvalues of adaptive simulations change at the level of the solver tolerances, and
    fewer converged eigenpairs beyond the requested ones can be reported. SchemaVer 4-4-2
    [PR 1203](https://github.com/awslabs/palace/pull/1203).

  - Rebalancing a nonconforming mesh during adaptive refinement now keeps each root of its
    refinement tree on one process, when the load stays within half of
    `Model.Refinement.MaximumImbalance`, so that saving the adapted mesh needs no
    redistribution; the partition, and so round-off, can differ from before.
    [PR 1214](https://github.com/awslabs/palace/pull/1214).

  - When `Solver.Linear.Type` is `"Default"` (or omitted) for frequency domain problems,
    and for the wave port and boundary mode eigenvalue solves, *Palace* now prefers
    STRUMPACK over SuperLU_DIST, then MUMPS, as the sparse direct solver when built with
    several of them. Set `"Type": "SuperLU"` to keep using SuperLU_DIST. SchemaVer 4-2-0
    [PR 1148](https://github.com/awslabs/palace/pull/1148).

  - The `mode-Z.csv` output of boundary-mode solves now reports the complex characteristic
    impedances, `Re{Z_PV}` and `Im{Z_PV}` with `Z_PV = |V|² / 2P*` as for wave ports, and
    `Re{Z_VI}` and `Im{Z_VI}` with `Z_VI = V / I` (oriented with the mode power) instead of
    `|V| / |I|`. The per-unit-length inductance and capacitance are now
    `L = Re{n_eff Z} / c0` and `C = Re{n_eff / Z} / c0`, which are unchanged for lossless
    modes.
    [PR 1119](https://github.com/awslabs/palace/pull/1119).

  - Circuit synthesis (`"AdaptiveCircuitSynthesis"`) now rejects a nonzero Floquet wave
    vector (`"FloquetWaveVector"`), as it already did Floquet ports, instead of synthesizing
    L⁻¹, R⁻¹, and C matrices which are not symmetric: the Floquet terms make the system
    matrices Hermitian rather than complex symmetric. Such configurations must be updated.
    [PR 1112](https://github.com/awslabs/palace/pull/1112).

  - `"AdaptiveGSOrthogonalization": "MGS"` now orthogonalizes the reduced-order basis of
    adaptive driven sweeps with `"CGS2"`, which is as stable with fewer reductions. The
    minimal rational interpolation of the error estimate always uses `"CGS2"`.
    SchemaVer 4-3-1
    [PR 1187](https://github.com/awslabs/palace/pull/1187).

  - Added a warning for interior `"Conductivity"` boundaries without `"Crack"`, such as the
    surface of a conductor meshed as a dielectric volume. Since
    [PR 994](https://github.com/awslabs/palace/pull/994) these are uncracked sheets with the
    admittance of both faces, which lowers the conductor loss compared to previous versions:
    set `"Crack": true` to restore the previous results, or `"Crack": false` to silence the
    warning.
    [PR 1168](https://github.com/awslabs/palace/pull/1168).

  - Uniform lumped ports and surface current sources now take their width from a
    cross-section across the port, and their length as their area over their width (in 2D,
    the width counts the edges crossing it). This fixes setup aborting for ports with curved
    ends, folded ports, and cylindrical ports. Ports which are not strips along
    `"Direction"`, such as L-shaped ports, and 2D ports with `"Direction"` normal to their
    edge, now abort.
    [PR 1194](https://github.com/awslabs/palace/pull/1194).

  - `config["Solver"]["Linear"]["InitialGuess"]` now defaults to `true` for adaptive
    frequency sweeps too. There, the HDM solve of the sample which would complete the
    convergence memory starts from the PROM solution, to certify the PROM: if the linear
    solver converges, to the tighter of `"Tol"` and `"AdaptiveTol"`, and changes it by less
    than `"AdaptiveTol"`, the sweep stops, and this solution is added to the basis as the
    last sample (unless it is already in the span of the basis, for example when the PROM
    solution already meets this tolerance); otherwise, the sample is solved again from zero
    (reusing the operators and the preconditioner), as are all other samples. With GMRES and
    left preconditioning (the defaults), a certification which fails stops as soon as the
    change exceeds `"AdaptiveTol"`. With an explicit `"InitialGuess": true`, adaptive sweeps
    used to start each HDM solve from the solution at the previous sample. When the PROM is
    certified, the last HDM solve of the excitation needs fewer linear solver iterations,
    and the results of adaptive sweeps change by about the accuracy of the PROM. Set
    `"InitialGuess": false` to restore the previous default behavior. SchemaVer 5-0-0
    [PR 1126](https://github.com/awslabs/palace/pull/1126).

#### Performance Improvements

  - Added a hypre patch, a prototype for
    [hypre-space/hypre#1654](https://github.com/hypre-space/hypre/issues/1654), with which
    AMS gathers the coarse levels of its auxiliary space BoomerAMG hierarchies onto fewer
    processes during the solve and runs the additive Pi-space corrections concurrently.
    `config["Solver"]["Linear"]["AMSAgglomerationTarget"]` sets the target number of rows
    per process (`0` disables it). Results agree up to the summation order. SchemaVer 5-2-0
    [PR 1154](https://github.com/awslabs/palace/pull/1154).
  - When AMS is the coarse solver of geometric multigrid for a singular operator (see
    [`"AMSSingularOperator"`](https://awslabs.github.io/palace/dev/config/reference/#config-solver-linear-amssingularoperator)),
    each AMS iteration now applies 8 to 32 smoothing sweeps before and after its auxiliary
    space corrections, instead of `"MGSmoothIts"`: the more, the larger the finer multigrid
    levels are next to the coarsest one. This lowers the number of linear solver iterations.
    [PR 1159](https://github.com/awslabs/palace/pull/1159).
  - The AMS solver, used with a fixed number of iterations, now has a tolerance of zero and
    flags its zero initial guess. With a hypre patch, it then computes no residual norms,
    each a fine-level matrix-vector product and a reduction over all processes, and the
    first relaxation of its auxiliary space solves skips the exchange of their zero initial
    guess. Results are unchanged.
    [PR 1238](https://github.com/awslabs/palace/pull/1238).
  - With `Model.Refinement.SaveAdaptMesh`, a nonconforming adapted mesh is saved without
    gathering it on the root: every process prints its part of the file, and the root only
    writes the compressed parts. The saved file holds the same mesh and refinement tree as
    before, numbered differently (the roots are listed process by process). Meshes the
    parallel save does not support are still gathered, with a warning.
    [PR 1214](https://github.com/awslabs/palace/pull/1214).
  - A restart from a nonconforming adapted mesh which *Palace* saved reads the mesh in
    parallel: every process reads a slice of the file (the root decompresses a compressed
    file and sends each process its slice), the coarse mesh is rebuilt from the roots of the
    refinement tree, and each process refines its roots and reads their nodes, instead of
    every process reading the whole mesh.
    [PR 1215](https://github.com/awslabs/palace/pull/1215).
  - Added a STRUMPACK patch for
    [pghysels/STRUMPACK#174](https://github.com/pghysels/STRUMPACK/pull/174), which *Palace*
    always enables: the STRUMPACK solvers store their fill-reducing permutation once per
    node, in an MPI shared memory window, instead of on every process, for the sequential
    column orderings (all but ParMETIS and PT-Scotch) without compression. If the window
    cannot be allocated, or with `OMPI_MCA_osc=^sm,ucx` for Open MPI, every process keeps
    its copy. Results are unchanged.
    [PR 1231](https://github.com/awslabs/palace/pull/1231).
  - Complex operators on CPUs which pack matching pairs of real and imaginary libCEED terms
    now also pair the H1 diffusion terms and pack the real terms without a partner, such as
    those of the lossless materials of partially lossy models, so that each of them acts on
    the real and imaginary parts in one application. This speeds up every multigrid level
    above the coarsest, and its H1 auxiliary operator, in driven and eigenmode simulations
    with losses. Operators without any pair, as in lossless problems, keep their original
    action, which needs less memory. Results change at most at the level of round-off.
    [PR 1173](https://github.com/awslabs/palace/pull/1173).
  - With `config["Solver"]["Linear"]["ComplexCoarseSolve"]`, SuperLU_DIST (with or without
    `"SuperLU3DCommunicator"`), STRUMPACK (with or without compression), and MUMPS now
    factor the complex-valued coarse-level matrix in complex arithmetic instead of the
    equivalent real-valued system of twice the size, which reduces the memory of the factors
    by about 40%. This changes the behavior of an explicit `"ComplexCoarseSolve": true`: the
    factorization differs, so results change at the level of the linear solver tolerance,
    and by more with STRUMPACK compression, which now compresses the factors of the
    complex-valued matrix. With `"DropSmallEntries"`, they drop the small entries of the
    complex-valued matrix, by the rule which is applied to the real-valued system, and with
    `"PCMatReal"` the coarse-level matrix is real-valued and is factored as before. When the
    imaginary part of the complex-valued matrix is exactly zero, as for a model without
    losses, these solvers factor its real part in real arithmetic, which halves the memory
    of the factors again. The log says which matrix is factored, and whether in complex
    arithmetic or as a real-valued system. With the AMS and BoomerAMG linear solvers, which
    cannot be given the real-valued system of twice the size (AMS aborted, and simulations
    with BoomerAMG did not converge), `"ComplexCoarseSolve": true` is now disabled with a
    warning. Wave ports keep it for their mode solver, which uses a sparse direct solver.
    SchemaVer 4-4-0
    [PR 1193](https://github.com/awslabs/palace/pull/1193).
  - The linear solves of driven simulations without adaptive frequency sweep now start from
    the minimal-residual combination of the solutions at up to four previous frequencies,
    and GMRES and FGMRES scale nonzero initial guesses to minimize the initial residual
    (except when certifying the PROM solution of an adaptive frequency sweep). This reduces
    the iterations of frequency sweeps. SchemaVer 5-0-1
    [PR 1229](https://github.com/awslabs/palace/pull/1229).
  - Each excitation of a uniform driven sweep now starts its first linear solve from zero
    instead of from the last solution of the previous excitation, which saves a few
    iterations for every excitation after the first.
    [PR 1219](https://github.com/awslabs/palace/pull/1219).
  - Added a STRUMPACK patch for
    [pghysels/STRUMPACK#171](https://github.com/pghysels/STRUMPACK/pull/171), which lowers
    the memory of the distributed STRUMPACK solver by 12 bytes per coarse matrix row on
    each process: the row-to-front map and the column scaling were stored for the whole
    matrix on every process. Results are unchanged.
    [PR 1198](https://github.com/awslabs/palace/pull/1198).
  - STRUMPACK is now built with a PBLAS block size of 128 instead of 32 in the superbuild
    and in Spack builds of STRUMPACK without SLATE, the default of the Spack package since
    [PR 1223](https://github.com/awslabs/palace/pull/1223). Spack builds of STRUMPACK with
    SLATE keep the block size of 512 of SLATE. The block size is set with a STRUMPACK
    patch for [pghysels/STRUMPACK#169](https://github.com/pghysels/STRUMPACK/pull/169).
    The patch for [pghysels/STRUMPACK#165](https://github.com/pghysels/STRUMPACK/pull/165)
    ([PR 1170](https://github.com/awslabs/palace/pull/1170)) is still applied, but its
    lower-latency distributed triangular solves for a few right-hand sides only handle
    block sizes up to 32 and are thus no longer used. With a block size of 128, the time
    of the STRUMPACK coarse solver is 4-17% lower than with these solves and a block size
    of 32, and its memory is lower.
    [PR 1197](https://github.com/awslabs/palace/pull/1197).
  - Used two levels of aggressive coarsening in the BoomerAMG coarse solver of the
    multigrid preconditioner of the divergence-free projection, instead of one. Its
    V-cycle is latency bound at high process counts, and the shorter hierarchy makes it
    faster.
    [PR 1149](https://github.com/awslabs/palace/pull/1149).
  - Uniform frequency sweeps of driven simulations release the linear solver and system
    matrices after the last solve, before the final field output and error estimation.
    This lowers the peak memory when the final postprocessing sets it.
    [PR 1087](https://github.com/awslabs/palace/pull/1087).
  - Added a hypre patch, based on
    [hypre-space/hypre#1652](https://github.com/hypre-space/hypre/pull/1652), for
    multicomponent vectors in hybrid Gauss-Seidel relaxation and Gaussian elimination, and
    the BoomerAMG preconditioner now applies one V-cycle to several right-hand sides at
    once, such as the real and imaginary parts of a complex vector, with one halo exchange
    per level for all of them. This speeds up the divergence-free projection at high
    process counts. Results are unchanged.
    [PR 1145](https://github.com/awslabs/palace/pull/1145).
  - Added a hypre patch for
    [hypre-space/hypre#1651](https://github.com/hypre-space/hypre/pull/1651), which skips
    the halo exchange of a zero initial guess in the hybrid Gauss-Seidel relaxation of
    BoomerAMG and AMS, and marked the zero initial guess of the BoomerAMG preconditioner as
    zero. This saves one halo exchange per level and V-cycle, which is latency bound on the
    coarse levels at high process counts. Results are unchanged.
    [PR 1139](https://github.com/awslabs/palace/pull/1139).
  - Added a STRUMPACK patch for
    [pghysels/STRUMPACK#163](https://github.com/pghysels/STRUMPACK/pull/163), which lowers
    the communication cost of its distributed solve phase: point-to-point exchanges of only
    the nonempty messages instead of all-to-all collectives, cached right-hand side
    redistribution patterns and extend-add index maps, and the row interchanges of the
    forward solve applied in one exchange. This speeds up the repeated solves of STRUMPACK
    as the coarse solver at high process counts, with unchanged results.
    [PR 1143](https://github.com/awslabs/palace/pull/1143).
  - SuperLU_DIST now uses relaxed supernodes of up to 200 columns instead of 30. This gives
    fewer, larger blocks and many fewer messages in its distributed triangular solves,
    which speeds up the repeated solves of SuperLU_DIST as the coarse solver at high
    process counts. Results change at the level of round-off. It uses a patch for
    [mfem/mfem#5731](https://github.com/mfem/mfem/pull/5731), which adds the setter of the
    relaxed supernode size to `mfem::SuperLUSolver`.
    [PR 1166](https://github.com/awslabs/palace/pull/1166).
  - Added SuperLU_DIST patches for
    [xiaoyeli/superlu_dist#250](https://github.com/xiaoyeli/superlu_dist/pull/250),
    [xiaoyeli/superlu_dist#251](https://github.com/xiaoyeli/superlu_dist/pull/251),
    [xiaoyeli/superlu_dist#252](https://github.com/xiaoyeli/superlu_dist/pull/252) and
    [xiaoyeli/superlu_dist#253](https://github.com/xiaoyeli/superlu_dist/pull/253), which
    lower the communication cost of SuperLU_DIST as the coarse solver at high process
    counts: point-to-point exchanges of only the nonempty messages instead of dense
    all-to-all collectives in the symbolic factorization and in the final permutation of
    the solution, a static schedule reused by refactorizations with the same sparsity
    pattern, cached row indices for the redistributions of the right-hand side, which on
    several nodes go through one process per node, the matrix values for MC64 gathered on
    one process only, and fewer barriers in the triangular solves. Results are unchanged.
    [PR 1162](https://github.com/awslabs/palace/pull/1162).
  - Skipped the eigenvalue estimate of the Jacobi preconditioner of the error estimator's
    flux recovery, on which the preconditioned conjugate gradient iterates do not depend.
    This speeds up the construction of the error estimator, which is repeated at every
    adaptive mesh refinement iteration.
    [PR 1134](https://github.com/awslabs/palace/pull/1134).
  - Solved the second divergence-free projection in each shift-and-invert application of
    the linearized quadratic eigenvalue problem (eigenmode problems with damping) to a
    relative tolerance of 10⁻³. Its right-hand side is already at round-off, which the
    projection previously reduced to machine precision.
    [PR 1138](https://github.com/awslabs/palace/pull/1138).
  - Operator norm estimates, such as the maximum eigenvalue estimate of each Chebyshev
    smoother, now use the Lanczos method with the same tolerance in all builds, instead of
    SLEPc or the power iteration. Each Chebyshev smoother also starts its estimate from the
    dominant vector of its previous operator, and estimates real operators in real
    arithmetic. This makes the multigrid preconditioner setup several times faster at
    each frequency of a driven sweep.
    [PR 1140](https://github.com/awslabs/palace/pull/1140).
  - Fused the residual update, Jacobi preconditioner and inner product of each conjugate
    gradient iteration into a single pass over memory on CPUs with one thread per process,
    and computed inner products on CPUs with several partial sums, complex ones in a single
    pass over memory. Each conjugate gradient iteration also updates the solution in the
    same pass as the next search direction, with any preconditioner and on CPUs and GPUs.
    This speeds up the error estimator's flux recovery and the GMRES orthogonalization of
    driven and eigenmode simulations.
    [PR 1146](https://github.com/awslabs/palace/pull/1146).
  - Expanded the reduced-order solution at each frequency of an adaptive sweep into the
    full space in one pass over memory on CPUs, reading each basis vector once.
    [PR 1165](https://github.com/awslabs/palace/pull/1165).
  - Sped up the offline phase of adaptive driven sweeps and lowered its memory use. Each
    Gram-Schmidt pass of a new basis vector, and of GMRES with `"CGS"` or `"CGS2"`
    orthogonalization, computes all inner products in one pass over the basis with one
    reduction, the projections onto the basis are summed in one reduction per update, and
    the error estimate no longer stores the full-size solutions.
    [PR 1187](https://github.com/awslabs/palace/pull/1187).
  - The lumped and wave port outputs of adaptive driven sweeps are evaluated without
    expanding the reduced-order solution at each frequency.
    [PR 1188](https://github.com/awslabs/palace/pull/1188).
  - The wave and Floquet port terms of the reduced-order model of adaptive driven sweeps
    are evaluated from the projections of the port modes, without assembling and projecting
    full-size vectors at each frequency, and the reduced system is assembled and solved once
    per frequency for all excitations (once per frequency and excitation with field output).
    A frequency-dependent boundary condition without a reduced-order term now stops the
    simulation instead of falling back to assembling it at each frequency.
    [PR 1189](https://github.com/awslabs/palace/pull/1189).
  - The reduced system of adaptive driven sweeps is solved with partial-pivoting LU instead
    of full-pivoting QR, falling back to the QR when the LU solution is not finite or
    inaccurate.
    [PR 1190](https://github.com/awslabs/palace/pull/1190).
  - The Floquet port DtN correction and the wave-port modal correction now use a single MPI
    reduction per operator application instead of four per rank-one term, and only touch
    the port degrees of freedom instead of full-length vectors.
    [PR 1089](https://github.com/awslabs/palace/pull/1089).
  - The Floquet port DtN correction and the wave-port modal correction now pack their
    rank-one terms per port, so that each term only touches the degrees of freedom of its
    own port, and compute their dot products on GPUs with block reductions instead of one
    thread per term.
    [PR 1155](https://github.com/awslabs/palace/pull/1155).
  - GPU builds now share one Umpire device memory pool between MFEM and hypre, created on
    each rank's GPU, and use standard host memory instead of the slow Umpire host allocator.
    hypre no longer reserves its own 4 GiB device pool per rank, which reduces peak GPU
    memory use substantially for small and medium problems (5.2 GiB to 1.25 GiB for the
    spheres example on one NVIDIA L4) and speeds up setup and postprocessing.
    [PR 1085](https://github.com/awslabs/palace/pull/1085).
  - Sped up the per-port wave-port PROM of adaptive driven sweeps for port modes whose shape
    changes with frequency, as for lossy lines, where every exact port solve adds a basis
    vector. Its projections and residual Gram matrix are now built with dense matrix
    products instead of one inner product per entry, a reduced solve checks the residuals of
    only the modes nearest to the target, and the reduced model is only used with a basis of
    at most `2√N` vectors for `N` port unknowns, beyond which exact port solves are cheaper.
    [PR 1119](https://github.com/awslabs/palace/pull/1119).
  - Driven simulations apply the system matrix through the finest level of the multigrid
    preconditioner and no longer assemble it separately for uniform sweeps, reducing CPU
    time by 1–5% and peak memory by up to 6%.
    [PR 980](https://github.com/awslabs/palace/pull/980).

#### Bug Fixes

  - Fixed the boundary-mode `Sn` output, which lacked the 1/μ₀ that converts B to H.
    [PR 937](https://github.com/awslabs/palace/pull/937).
  - Fixed S-parameter CSVs being written for driven excitations that are not a single
    power-normalized port, such as an impressed current. Their port overlaps are not
    S-parameters, and a fresh run removes the S-parameter CSVs of an earlier run.
    [PR 937](https://github.com/awslabs/palace/pull/937).
  - Fixed sparse direct solves with `config["Solver"]["Linear"]["ComplexCoarseSolve"]`
    which could crash, hang, or return NaN when a process owns no unknowns of the solve: in
    `BoundaryMode` and wave port mode solves, and on the coarsest multigrid level with the
    solvers which factor the real-valued system of twice the size.
    [PR 1206](https://github.com/awslabs/palace/pull/1206).
  - Fixed slow and failing linear solves with the second-order absorbing boundary condition
    (`"Absorbing"` with `"Order": 2`) without `"ComplexCoarseSolve"`: the sign of its term
    is now flipped in the real-valued coarse solve, as it already was for AMS and BoomerAMG,
    which now flip it only there.
    [PR 1252](https://github.com/awslabs/palace/pull/1252).
  - A failed write of an output file, for example on a full disk, aborts the run instead
    of leaving the file empty or truncated, and a failed rewrite keeps the previous file.
    [PR 1256](https://github.com/awslabs/palace/pull/1256).
  - Fixed boundary terms that could take the anisotropic material tensor of another
    boundary. The scalar coefficient of an `"Impedance"`, `"Conductivity"` or
    `"RationalImpedance"` boundary or a lumped port could be replaced by the tensor of an
    `"Absorbing"` boundary, wave port or Floquet port on an anisotropic material when the
    first diagonal entry of the tensor was equal to the scalar to a relative 1e-9. In one
    `"BoundaryMode"` example this changed the propagation constant by 6%.
    [PR 1209](https://github.com/awslabs/palace/pull/1209).
  - The `"Problem"`/`"Iteration"` field of `palace.json` is now updated only after the AMR
    iteration's output is in its final location. With `SaveAdaptIterations`, it used to be
    written before the output was moved into `iterationXX`, so a reader polling it could
    find a partially filled `iterationXX` directory.
    [PR 1235](https://github.com/awslabs/palace/pull/1235).
  - The `Initialization` time now includes MPI initialization and the parsing of the
    configuration file, which were counted in `Total` but in no row of the report (up to
    tens of seconds on many nodes). MPI initialization, which includes waiting for all
    processes to start, has its own row, `MPI Startup`.
    [PR 1242](https://github.com/awslabs/palace/pull/1242).
  - Fixed the error estimates being all zero when MFEM uses a device but libCEED a CPU
    backend, as with `"Device": "Debug"` or a GPU run with a CPU `"Backend"`, which kept
    adaptive mesh refinement from refining. Default GPU runs were not affected.
    [PR 1241](https://github.com/awslabs/palace/pull/1241).
  - Fixed magnetostatic solves that did not converge, and gave wrong inductances, when the
    PEC boundaries allow a constant vector potential: no PEC at all, or PEC only on planes
    normal to one axis, such as planar films with a natural far field. For singular
    operators, AMS now relaxes the coarsest level of its auxiliary space solves instead of
    using Gaussian elimination.
    [PR 1244](https://github.com/awslabs/palace/pull/1244).
  - Fixed the pivoting of distributed BLR fronts in STRUMPACK with a patch from upstream.
    With BLR compression (`"STRUMPACKCompressionType"` `"BLR"`, `"BLR-HODLR"` or
    `"ZFP-BLR-HODLR"`) on two or more MPI processes, the factorization of an indefinite
    matrix was inaccurate whatever the tolerance, so the linear solves needed more
    iterations and could return wrong solutions with `"ComplexCoarseSolve"`.
    [PR 1222](https://github.com/awslabs/palace/pull/1222).
  - Fixed adaptive frequency sweeps which could keep sampling up to `"AdaptiveMaxSamples"`,
    with spurious resonances in the band, when `"AdaptiveTol"` is close to the accuracy of
    the HDM solutions. *Palace* now refines the samples which the PROM does not reproduce to
    `"AdaptiveTol"`, no longer aborts on a sample in the span of the basis, and warns when
    the sampling reaches `"AdaptiveMaxSamples"` before converging.
    [PR 1218](https://github.com/awslabs/palace/pull/1218).
  - Fixed the flux recovery in the error estimators to use the projection which minimizes
    the error it is measured in (ε⁻¹-weighted for the electric flux, μ-weighted for the
    magnetic field). The unweighted projection overestimated the error in the
    low-permittivity elements next to a dielectric interface, by up to an order of magnitude
    next to a sapphire substrate, and AMR spent refinement there. The error indicator
    decreases for models with material interfaces, and their AMR meshes change.
    [PR 1072](https://github.com/awslabs/palace/pull/1072).
  - Scaling a nonconforming mesh without nodes (for example with
    `Model.RemoveCurvature`) now also scales the top-level vertex coordinates kept by its
    refinement tree, so that its saved adapted mesh and restarts from it have the right
    length unit.
    [PR 1214](https://github.com/awslabs/palace/pull/1214).
  - Fixed wave port sweeps with a degenerate pair of port modes, where the polarization
    selected within the pair could jump between frequencies: the fields of the sweep were
    then not those of one excitation, and an adaptive sweep could excite another
    polarization than it was trained with. Every port mode solve now starts from the same
    vector, so that its result, and with it the polarization, depends on the frequency
    only, also in `BoundaryMode` simulations. The index of a mode above an exactly
    degenerate pair can still depend on roundoff.
    [PR 1208](https://github.com/awslabs/palace/pull/1208).
  - Fixed driven frequency sweeps aborting with "... not found in sample frequencies!" when
    overlapping `"Samples"` ranges compute the same frequency with different rounding.
    [PR 1192](https://github.com/awslabs/palace/pull/1192).
  - Fixed the domain magnetic energies of adaptive driven sweeps evaluated from the
    reduced-order model, which lost up to 10⁻² relative accuracy at low frequency, where
    the curls of the basis vectors nearly cancel, and were only checked against the full
    evaluation for the first excitation. All reduced outputs are now checked against the
    full evaluation at the first frequency of every excitation (the first one without field
    output), which is the one written to the output files, reverting to the full
    postprocessing with a warning if they differ by more than a tenth of `"AdaptiveTol"`
    (or about 10⁻⁹ for tighter tolerances).
    [PR 1188](https://github.com/awslabs/palace/pull/1188).
  - Fixed the coupling of a lumped port in the circuit of adaptive circuit synthesis
    (`rom-coupled-G.csv`, `rom-coupled-H.csv`, and `rom-coupled-S.csv`). The port was
    coupled through a unit selector on its terminal row, but the basis vector of that row,
    the projection of the port mode normalized in the W norm, does not have a unit voltage,
    and the other basis vectors carry small voltages. It is now coupled through its voltage
    on every basis vector. The S-parameters of the circuit of a shorted coaxial line driven
    by a lumped port match those of the reduced-order model to 1e-6 instead of 2e-2, and
    those of the `coaxial_lumped_wave_synth` regression case to 7e-13 instead of about
    6e-4. These files change for every circuit synthesis with lumped ports.
    [PR 877](https://github.com/awslabs/palace/pull/877).
  - Fixed circuit synthesis with wave ports aborting, without writing `port-S.csv` or the
    `rom-*` files, when a wave port mode did not converge at the complex frequency of a
    synthesized eigenvalue estimate. The HDM errors of such estimates are now written as
    NaN, with a warning, and the other outputs are written as before. For this, an ARPACK
    eigenvalue solve which reaches its iteration limit now returns the eigenpairs which have
    converged, as SLEPc does, instead of aborting.
    [PR 877](https://github.com/awslabs/palace/pull/877).
  - Fixed errors of up to 1e-4 in the S-parameters of a synthesized circuit with wave ports
    in a lossy or frequency-dependent material, relative to the reduced-order model the
    circuit is made from, while the reported fit residuals were below `AdaptiveTol`. The two
    parts of a wave port fit were realized on different sets of port directions, and the
    rounding-level imaginary part of a real wave port mode was kept as a circuit node.
    The circuits of existing inputs change: of a wave port with poles in its fit also in
    vacuum (the S-parameters of the circuit of the `adapter_driven_synth` regression case
    by 5e-8), and of a wave port with a real mode, which loses its `waveport_<n>_im` node
    and one auxiliary direction for each pole of its fit (in the
    `cylinder_driven_wave_tm_lossy_synth` regression case, also its exported reference
    impedance by 5e-4).
    [PR 877](https://github.com/awslabs/palace/pull/877).
  - Changed the fits of circuit synthesis for surface conductivity boundaries, and for
    rational impedance boundaries which are not split exactly, to bound their error relative
    to the fitted response at each frequency instead of relative to its largest value over
    the band, which did not bound the error where the response is small. The polynomial
    part of a fit with poles, of a surface conductivity or a rational impedance boundary
    (split exactly or not), is now realized on the directions of its poles only, as for the
    frequency-dependent materials: on the full mass, a direction below the rank cutoff kept
    the polynomial part alone, which for a series RLC sheet below its resonance left the
    term in the circuit 3e-5 from that of the reduced model where the fit residual was
    1e-13. The circuits of existing inputs with these boundaries change, including the
    support points of their AAA fits: the S-parameters of the circuit of the
    `coaxial_wave_conductivity_synth` regression case differ from those of its reduced-order
    model by 5e-6 instead of 9e-6.
    [PR 877](https://github.com/awslabs/palace/pull/877).
  - Added a warning when eigenvalue estimates of a synthesized circuit
    (`rom-eigenvalues.csv`) have Im{f} < 0. These roots grow in time, which their Q does not
    show. They come from the poles of the AAA fits of frequency-dependent boundaries, such
    as wave ports in lossy or frequency-dependent media, which are not constrained to the
    stable half-plane: the circuit can match the reduced-order model on the band and not be
    stable, for example in a transient circuit simulation.
    [PR 877](https://github.com/awslabs/palace/pull/877).
  - Fixed the quasi-Newton nonlinear eigenvalue solver stalling below
    `config["Solver"]["Eigenmode"]["PreconditionerLagTol"]`, where it kept a stale
    preconditioner and could stop without the requested modes. Below that residual, the
    preconditioner is now updated when the tolerance would otherwise not be reached within
    the next `"PreconditionerLag"` iterations.
    [PR 877](https://github.com/awslabs/palace/pull/877).
  - Added a STRUMPACK patch for
    [pghysels/STRUMPACK#164](https://github.com/pghysels/STRUMPACK/pull/164), which fixes
    the solves after the matrix values change when STRUMPACK compression is enabled
    (`Solver.Linear.STRUMPACKCompressionType` other than `"None"`) and the reordering is
    reused (`Solver.Linear.ReorderingReuse`, the default). With the STRUMPACK patch of
    [PR 1143](https://github.com/awslabs/palace/pull/1143), these solves could abort with
    an MPI error at high process counts, or give a poor preconditioner which needs more
    linear solver iterations.
    [PR 1164](https://github.com/awslabs/palace/pull/1164).
  - Fixed SuperLU_DIST solves with `config["Solver"]["Linear"]["SuperLU3DCommunicator"]`,
    which gave wrong solutions and crashed at exit. SuperLU_DIST patches for
    [xiaoyeli/superlu_dist#254](https://github.com/xiaoyeli/superlu_dist/pull/254) and
    [xiaoyeli/superlu_dist#255](https://github.com/xiaoyeli/superlu_dist/pull/255) fix the
    solutions with several right-hand sides and uneven rows, and the refactorizations with
    the same sparsity pattern, on a 3D process grid. An MFEM patch for
    [mfem/mfem#5732](https://github.com/mfem/mfem/pull/5732) frees the factors on all the
    process layers. ParMETIS ordering, which SuperLU_DIST does not support with a 3D process
    grid and with which it hung, now falls back to METIS with a warning, and a patch for
    [xiaoyeli/superlu_dist#256](https://github.com/xiaoyeli/superlu_dist/pull/256) makes
    SuperLU_DIST return an error for it instead of hanging.
    [PR 1167](https://github.com/awslabs/palace/pull/1167).
  - Fixed STRUMPACK sparse direct solves which could crash, fail with `ZERO_PIVOT`, or hang
    when a process holds very few rows of the matrix, for example a coarse multigrid level
    at a high MPI rank count. The sparsity pattern of the matrix passed to STRUMPACK is no
    longer declared symmetric, since hypre's parallel matrix products can store explicit
    zeros without their transposes.
    [PR 1137](https://github.com/awslabs/palace/pull/1137).
  - Added an MFEM patch for [mfem/mfem#5717](https://github.com/mfem/mfem/pull/5717),
    which fixes hangs and memory corruption on nonconforming meshes with triangular faces
    when the face neighbors of some processes have only quadrilateral faces.
    [PR 1124](https://github.com/awslabs/palace/pull/1124).
  - Added an MFEM patch for [mfem/mfem#5711](https://github.com/mfem/mfem/pull/5711),
    which fixes current dipole sources (`config["Domains"]["CurrentDipole"]`) with
    `config["Solver"]["Order"]` of 2 and above on meshes with tetrahedra, prisms or
    pyramids. Their amplitude, sign and direction were wrong, depending on the element
    containing the dipole.
    [PR 1121](https://github.com/awslabs/palace/pull/1121).
  - Fixed the Floquet port MaxOrder cap reading past the boundary element Jacobian, which
    gave a garbage element size and could randomly drop Floquet orders from the solve and
    the `port-floquet-S.csv` output.
    [PR 1091](https://github.com/awslabs/palace/pull/1091).
  - Fixed the ParaView and grid function output of the out-of-plane magnetic flux density
    `B` of 2D simulations and boundary modes, which was scaled by `L0 / Lc`, the mesh length
    unit over the reference length, along with the surface current `J_s` (and the Poynting
    vector `S` of 2D simulations), while the magnetic energy density `U_m` was scaled by
    `(L0 / Lc)^2`. This put them off by orders of magnitude for meshes in micrometers with
    the default `Lc`. The CSV postprocessing was not affected.
    [PR 1098](https://github.com/awslabs/palace/pull/1098).
  - Added an MFEM patch for [mfem/mfem#5721](https://github.com/mfem/mfem/pull/5721),
    which makes the refinement and rebalancing of nonconforming meshes in parallel
    independent of the arrival order of MPI messages. The numbering of adapted meshes, the
    meshes written with `"SaveAdaptMesh"`, and the round-off of the results computed on
    them no longer depend on the timing of MPI messages.
    [PR 1129](https://github.com/awslabs/palace/pull/1129).
  - Added an MFEM patch for [mfem/mfem#5686](https://github.com/mfem/mfem/pull/5686), which
    fixes nondeterministic, non-converging multigrid solves of multi-GPU runs with GPU-aware
    MPI ([#1061](https://github.com/awslabs/palace/issues/1061)).
    [PR 1226](https://github.com/awslabs/palace/pull/1226).
  - Fixed a segmentation fault reading 2D meshes with boundary elements which are not an
    edge of any element, for example a conductor face tagged on an outer boundary of the
    mesh, which the conductor (a hole in the mesh) covers. These boundary elements, which
    no domain element is attached to, are now removed with a warning.
    [PR 1100](https://github.com/awslabs/palace/pull/1100).
  - Fixed driven-sweep restarts with more than one excitation aborting with
    `Writing data table at incorrect index`, because the frequencies reloaded from the CSV
    tables (rounded to their printed precision) were compared exactly with the recomputed
    ones.
    [PR 1088](https://github.com/awslabs/palace/pull/1088).
  - Fixed the S-parameters of numeric wave ports with lossy port modes, for example with
    `"Conductivity"` or `"Impedance"` boundaries on the port boundary, or lossy dielectrics.
    The wave port boundary condition and excitation now use the complex propagation constant
    of the port mode, including its attenuation, instead of its real part. The S-parameters
    of lossy port modes are now the amplitudes of their pseudo-waves, from the unconjugated
    reaction of the field with the port mode, instead of the conjugated (power) projection
    onto it, which does not separate the incident and reflected waves of a lossy mode: this
    changes their values, but not those of lossless propagating modes. Below cutoff, the
    propagation constant is now the root with `Im{kₙ} < 0`, which decays away from the port,
    also for the modes of `"BoundaryMode"` simulations: for a lossless cross-section, the
    sign of its imaginary part followed the roundoff of the mode solve. The characteristic
    impedance `Z_PV` of a lossy port mode is now complex, `|V|² / P*` with the complex power
    `P` of the mode, consistent with the total-field impedance in `port-Z.csv`. Circuit
    synthesis now also fits the complex propagation constant of a lossy port mode, so that
    the synthesized port loads carry the attenuation of the port lines.
    [PR 1119](https://github.com/awslabs/palace/pull/1119).
  - Fixed a crash of circuit synthesis (`"AdaptiveCircuitSynthesis"`) when several wave
    ports share a boundary, for example inactive wave ports extracting the S-parameters of
    a second mode at the boundaries of the active ports.
    [PR 1119](https://github.com/awslabs/palace/pull/1119).
  - Fixed the characteristic impedance `Z_PV` of lossy modes in 2D boundary-mode solves, for
    example with `"Conductivity"` or `"Impedance"` boundaries on the conductors: the term of
    the mode power with the longitudinal electric field was conjugated wrongly, which added
    a spurious imaginary part to the power. This made `Z_PV` too low, by a factor of 3.2 for
    a strongly lossy coplanar-strip line at 1 GHz, and the power normalization of lossy
    modes wrong. Results are unchanged for lossless modes.
    [PR 1119](https://github.com/awslabs/palace/pull/1119).
  - Fixed the second-order absorbing boundary condition (`"Absorbing"` with `"Order": 2`)
    with Floquet periodic boundary conditions (`"FloquetWaveVector"`). Previously the
    absorbing boundary reduced to first order for an obliquely incident TE plane wave.
    For TM plane waves, the second-order condition is still the first-order one.
    [PR 1112](https://github.com/awslabs/palace/pull/1112).
  - Fixed the reduced-order model of adaptive frequency sweeps with a Floquet wave vector
    which does not depend on frequency (no `"FloquetReferenceFrequency"`), which projected
    the stiffness matrix as if it were symmetric.
    [PR 1112](https://github.com/awslabs/palace/pull/1112).
  - Fixed SLP and TOAR (`PEPLinear: false`) eigenmode solves aborting on a GPU.
    [PR 1068](https://github.com/awslabs/palace/pull/1068).
  - Fixed TOAR (`PEPLinear: false`) missing the eigenmodes near the target with the
    default `"Scaling": true`: the polynomial operators passed to SLEPc were not scaled,
    while the target, the preconditioner, and the returned eigenvalues were. TOAR can still
    return spurious eigenvalues near the target, with large backward errors; the default
    `"PEPLinear": true` is not affected.
    [PR 1068](https://github.com/awslabs/palace/pull/1068).
  - Fixed the preconditioners of the linear solves of `"BoundaryMode"` simulations with
    p-multigrid (`"MGMaxLevels"` greater than 1) or with `"PCMatReal"`, which left out loss
    terms of the system matrix, such as the conductivity of conductors meshed as volumes or
    the resistance of `"Impedance"` boundaries. Where these terms dominate the system
    matrix, the linear solves of the eigenvalue solver did not converge, and the modes were
    wrong. MUMPS now compresses the factorizations of boundary mode and wave port mode
    solves only with `"STRUMPACKCompressionType": "BLR"`. SchemaVer 4-2-2
    [PR 1178](https://github.com/awslabs/palace/pull/1178).
  - Fixed underestimated eigenmode resistive lumped-port loss when the field does not match
    the port mode, causing too-small coupling-rate magnitudes and too-large external quality
    factors.
    [PR 1086](https://github.com/awslabs/palace/pull/1086).
  - Fixed the grid function output (`Problem.OutputFormats.GridFunction`) of driven
    simulations with several excitations, where each excitation overwrote the saved fields
    of the previous one, so that only the fields of the last excitation were kept. As for
    the ParaView output, the fields of excitation `i` are now written to the subfolder
    `gridfunction/driven/excitation_i/`. The mesh, rank, and error indicator files are
    shared by all excitations and stay in `gridfunction/driven/`. The output of driven
    simulations with a single excitation is unchanged.
    [PR 1220](https://github.com/awslabs/palace/pull/1220).
  - Fixed SLEPc eigenvalue solves on more than one process that occasionally hung, failed
    (for example with "The inner product is not well defined"), or returned wrong modes
    when the LAPACK library did not return identical results on all processes, such as MKL
    on some Intel processors. All processes now use the first process's solution of the
    small dense problems of SLEPc's solvers.
    [PR 1251](https://github.com/awslabs/palace/pull/1251).
  - Fixed probe postprocessing aborting with "invalid host pointer access" in Debug builds
    on GPUs, and with MFEM's debug device.
    [PR 1042](https://github.com/awslabs/palace/pull/1042).
  - Added an MFEM patch for [mfem/mfem#5780](https://github.com/mfem/mfem/pull/5780),
    which fixes the unit tests on MFEM's debug device that could segfault at random in the
    AMS setup, as seen on 2 processes: with hypre built without GPU support, the debug
    device could protect a heap page of hypre with unrelated data on it.
    [PR 1266](https://github.com/awslabs/palace/pull/1266).

#### Documentation

  - The `microstrip` regression cases use the mesh of the new dispersive microstrip example,
    which also has the side walls and the ends of the line as boundaries: the inductance of
    `microstrip_pec` increases by 0.5%, from 5.641 to 5.669 pH, and that of
    `microstrip_superconductor` by 0.05%.
    [PR 877](https://github.com/awslabs/palace/pull/877).
  - Updated the configuration schema descriptions of `"AdaptiveCircuitSynthesis"`, which
    listed frequency-dependent boundary conditions it supports as unsupported, and of the
    absorbing boundary `"Order"`, and documented the limitations of the second-order
    absorbing boundary condition: it only improves on the first-order one for TE fields, and
    it is not passive for tangential wavenumbers above √2 times the wavenumber in the
    material. SchemaVer 4-2-1
    [PR 1112](https://github.com/awslabs/palace/pull/1112).

#### Build system

  - Fixed the unit-test build against Catch2 3.15.2 and later, which failed intermittently
    after linking: `catch_discover_tests` lists the tests through a file that every MPI
    process wrote at once. Only the root process writes it now.
    [PR 1248](https://github.com/awslabs/palace/pull/1248).
  - Fixed MPI runs of several unit tests at once (`mpirun -np 2 palace-unit-tests`) with
    Catch2 3.9.0 and later, which run the tests in random order by default: the processes
    used different seeds, ran the tests in different orders, and their collectives did not
    match. All processes now use the seed of the root.
    [PR 1248](https://github.com/awslabs/palace/pull/1248).
  - With MPICH and the MPI implementations derived from it, the superbuild patches the BLACS
    of ScaLAPACK to receive point-to-point messages with derived datatypes packed and unpack
    them locally. This works around an MPICH bug for processes on the same node
    ([pmodels/mpich#8169](https://github.com/pmodels/mpich/issues/8169), fixed on the main
    branch of MPICH by [pmodels/mpich#8172](https://github.com/pmodels/mpich/pull/8172),
    which no release up to 5.0.2 and 5.1.0b1 has), seen on macOS arm64, which corrupts or
    aborts the complex-valued distributed factorizations of STRUMPACK and MUMPS: a packed
    message received in several shared memory cells cannot be unpacked when a cell boundary
    splits a complex-valued element. Builds with other MPI implementations are unchanged.
    The ScaLAPACK of a Spack build is not patched.
    [PR 1193](https://github.com/awslabs/palace/pull/1193).
  - The superbuild now builds MUMPS with complex double precision (`zmumps`), which
    *Palace* now requires with `PALACE_WITH_MUMPS`, also in the Spack package
    (`mumps@5.3:+double+complex`).
    [PR 1193](https://github.com/awslabs/palace/pull/1193).
  - The superbuild now builds ParMETIS, which *Palace* links directly, whenever it is
    enabled. Previously it was built only with SuperLU_DIST, STRUMPACK, or MUMPS, so
    configuring with cuDSS as the only sparse direct solver failed.
    [PR 1133](https://github.com/awslabs/palace/pull/1133).
  - The superbuild's sub-makes (the unit tests, GSLIB, LIBXSMM, libCEED, SLEPc, and MFEM's
    test step) now inherit the make jobserver, so GNU make 4.3 and later build them in
    parallel instead of serially. The superbuild also skips Eigen's BLAS, LAPACK, and demos,
    which *Palace* does not use, and clones PETSc and SLEPc shallowly.
    [PR 1180](https://github.com/awslabs/palace/pull/1180).
  - The Spack package now builds STRUMPACK without SLATE by default, as the superbuild does.
    [PR 1223](https://github.com/awslabs/palace/pull/1223).
  - ParMETIS is now optional and disabled by default: enable it with
    `-DPALACE_WITH_PARMETIS=ON` (Spack variant `+parmetis`). Otherwise, `"ParMETIS"` for
    `config["Solver"]["Linear"]["ColumnOrdering"]` falls back to METIS with a warning.
    SchemaVer 4-4-3
    [PR 1195](https://github.com/awslabs/palace/pull/1195).
  - Building the unit tests now requires MFEM with exceptions (`PALACE_MFEM_USE_EXCEPTIONS`
    in the superbuild, `mfem+exceptions` in Spack). The tests of error handling expect
    `MFEM_VERIFY` and `MFEM_ABORT` to throw, and without exceptions the first of them ended
    the whole test run. Without them, `make palace-tests` in the superbuild now fails right
    away instead of after building all the dependencies.
    [PR 1240](https://github.com/awslabs/palace/pull/1240).
  - The superbuild builds fmt 12.1.0 instead of 12.2.0, which does not compile with nvcc,
    and builds Umpire with it instead of letting Umpire install its own copy of fmt over it.
    [PR 1228](https://github.com/awslabs/palace/pull/1228).
  - Removed the `PALACE_WITH_GPU_AWARE_MPI` option: *Palace* now asks MPI at runtime whether
    it is GPU aware.
    [PR 1226](https://github.com/awslabs/palace/pull/1226).
  - Spack builds with `+cudss` now allow cuDSS 0.8, which the builtin `mfem` package capped
    at 0.7.1, with a patch for
    [spack/spack-packages#7693](https://github.com/spack/spack-packages/pull/7693), which
    also applies [mfem/mfem#5758](https://github.com/mfem/mfem/pull/5758).
    [PR 1211](https://github.com/awslabs/palace/pull/1211).
  - With cuDSS, the superbuild now patches MFEM with
    [mfem/mfem#5758](https://github.com/mfem/mfem/pull/5758), as the Spack build does. With
    cuDSS 0.8, the parallel cuDSS solver of MFEM 4.10 did not set the communicator for
    device buffers, and `"Type": "cuDSS"` in `config["Solver"]["Linear"]` crashed in the
    first factorization.
    [PR 1250](https://github.com/awslabs/palace/pull/1250).
  - The unit tests now pass with MFEM's debug device, and the long tests run them on it, on
    1 and 2 processes, to catch host accesses to data whose valid copy is on the device. The
    installed test tree now also discovers each case on the number of processes it runs on,
    which adds the `[Parallel]`-only regression cases to the Spack test jobs.
    [PR 1241](https://github.com/awslabs/palace/pull/1241).
  - The local Spack repository (`spack_repo/local`) now has an `mfem` package which derives
    from the builtin one, so *Palace* builds from it without linking the builtin MFEM recipe
    into it by hand. Remove such a link (`spack_repo/local/packages/mfem`) before updating.
    [PR 1151](https://github.com/awslabs/palace/pull/1151).

## [0.19.0] - 2026-11-02

#### New Features

  - Added total (geometric + kinetic) inductance extraction for thin-film superconductors. A
    `Boundaries.Superconductor` sheet carries the one-sided London kinetic sheet inductance
    `L_ksq = mu0 * lambda * coth(d/lambda)` (via `PenetrationDepth`/`Thickness`, or directly
    via `KineticInductance`); a `FluxLoop` over such a film imposes a prescribed fluxoid and
    the total inductance is read from the field energy.
    [PR 929](https://github.com/awslabs/palace/pull/929).
  - Added `Solver.Linear.LondonPCShift`, a preconditioner-only gauge shift that keeps the
    London magnetostatic solve SPD-solvable by AMS. Magnetostatic simulations with
    `Superconductor` or `FluxLoop` boundaries also enable
    `Solver.Linear.MGAuxiliarySmoother` by default. SchemaVer 3-0-1
    [PR 929](https://github.com/awslabs/palace/pull/929),
    [PR 1051](https://github.com/awslabs/palace/pull/1051).
  - Introduced the `"AbsTol"` option for linear solvers, defaulting to 0.0, to allow using
    an absolute tolerance when defining convergence. SchemaVer 2-1-0
    [PR 734](https://github.com/awslabs/palace/pull/734).
  - Added Method of Manufactured Solutions (MMS) verification for the electrostatic solver,
    confirming optimal `L2` convergence of the potential and recovered field on affine
    hexahedra and curved (second-order) tetrahedra, with polynomial exactness checked on
    affine hexahedra and tetrahedra, and covering non-homogeneous Dirichlet and Neumann
    boundaries and anisotropic permittivity. The manufactured source terms are derived
    either analytically or, optionally, via
    [Enzyme-AD](https://github.com/EnzymeAD/Enzyme), an LLVM-level automatic
    differentiation tool, enabled with a new `PALACE_WITH_ENZYME=ON` build option or the
    `+enzyme` Spack variant (used only to build the unit tests).
    [PR 901](https://github.com/awslabs/palace/pull/901).
  - Added `Boundaries.FluxLoop[].Excitation`, the port excitation convention for flux
    loops: loops sharing a positive excitation index are driven together in one
    magnetostatic solve, each at its own `FluxAmounts`, and `0` holds a loop at zero
    fluxoid. With grouped excitations `terminal-Phi.csv` reports the stored energy and loop
    fluxoids of each excitation instead of the inductance matrix, and the ParaView fields
    are saved per excitation. The default, each loop its own excitation, is unchanged.
    SchemaVer 3-1-0
    [PR 979](https://github.com/awslabs/palace/pull/979).
  - Error estimation for adaptive mesh refinement now uses broken finite element spaces,
    allowing discontinuous fields across interior boundaries without cracking the mesh. Mesh
    cracking is now a per-boundary physical modeling choice (see Interface Changes).
    [PR 994](https://github.com/awslabs/palace/pull/994).
  - Added `Boundaries.Superconductor[].TwoSided`, which models a London film as two coupled
    coincident faces (the exact two-port response of a slab of thickness `d`), capturing
    screening on both faces when `d` is comparable to `lambda`. With the same field on both
    faces, the film reduces to `L_ksq = (mu0 * lambda / 2) * coth(d/(2 * lambda))`. The
    magnetostatic solution space is broken across the film, so the mesh is not cracked.
    Requires the `PenetrationDepth`/`Thickness` form, a 3D mesh, and an interior film.
    SchemaVer 4-1-0
    [PR 992](https://github.com/awslabs/palace/pull/992).

#### Interface Changes

  - Renamed the `FluxLoop` keys `FluxLoopPEC` to `FilmAttributes` and `Regularization` to
    `PecPenetrationDepth`, and limited each `FluxLoop` to a single hole (the schema still
    accepts arrays for `HoleAttributes`/`FluxAmounts`). Existing `FluxLoop` configurations
    must be updated. SchemaVer 2-0-0
    [PR 929](https://github.com/awslabs/palace/pull/929).
  - Added `Solver.Eigenmode.NonlinearMaxIts` (default 100), the maximum number of
    quasi-Newton iterations of each attempt to converge an eigenpair in the hybrid nonlinear
    eigenvalue solver. `Solver.Eigenmode.MaxIts` no longer limits these iterations, and
    only applies to the linear or quadratic eigenvalue solver computing the initial guesses
    of the hybrid solver. SchemaVer 2-2-0
    [PR 1036](https://github.com/awslabs/palace/pull/1036).
  - Removed the `"NSamples"` option of the wave port `"VoltagePath"` and of the
    `"Impedance"` and `"Voltage"` boundary postprocessing paths. Coordinate-path line
    integrals are now exact for the discrete field, so the number of samples no longer has
    any effect, and configurations specifying it must be updated. SchemaVer 3-0-0
    [PR 1049](https://github.com/awslabs/palace/pull/1049).
  - Added the `"Crack"` option to `"Impedance"`, `"RationalImpedance"`, and `"Conductivity"`
    boundaries, which selects the physical model of an interior sheet: a thin film coupling
    the fields on its two sides (`false`, the mesh is not modified) or a conductor much
    thicker than its penetration depth with two independent faces, each with the given
    surface impedance as on an exterior boundary (`true`, the mesh is cracked along the
    boundary). The default is `false`.
    [PR 994](https://github.com/awslabs/palace/pull/994).
  - The mesh is now only cracked by default along interior boundaries with boundary
    conditions which apply to each side separately (PMC or zero charge, and absorbing
    boundaries), and along impedance, rational impedance, and conductivity boundaries with
    `"Crack": true`. Interior PEC, ground, wave port PEC, impedance, rational impedance,
    conductivity, and surface current boundaries are no longer cracked by default. For a
    conductor much thicker than its skin depth with different currents on its two faces,
    this lowers the conductor loss compared to the cracked mesh: use `"Crack": true` to keep
    the previous results. Cracked impedance and rational impedance boundaries no longer
    split the given surface impedance between their two faces: use `"Crack": true` with
    twice the impedance (twice `"Rs"` and `"Ls"`, half `"Cs"`) to keep the previous results.
    Specifying `config["Model"]["CrackInternalBoundaryElements"]` is deprecated: `true`
    restores the previous behavior of cracking all interior boundaries with boundary
    conditions except for lumped ports, with split impedances, and `false` disables cracking
    entirely, which is an error for boundaries with `"Crack": true`. SchemaVer 4-0-0
    [PR 994](https://github.com/awslabs/palace/pull/994).
  - Results change for interior boundaries which are no longer cracked by default. Cracking
    required a local bisection of the elements next to the boundary, which is no longer
    done: on the same mesh, cracked and uncracked PEC sheets give the same results, but the
    coarser mesh changes results within the discretization error, for example by about 3.4%
    for the metal-air and metal-substrate participation ratios of the `cpw` regression
    cases, and by up to 2.8% for their external quality factors. Cracking an interior
    surface current excitation decoupled the fields on the two sides of the source sheet:
    the inductance of the `microstrip` regression case, for example, decreases from 7.43 pH
    to 5.64 pH with PEC conductors (below the parallel plate estimate of 6.28 pH, which
    neglects the fringing fields), and by 12% with superconducting conductors. Interior
    conductivity boundaries which are not cracked, also with
    `config["Model"]["CrackInternalBoundaryElements"]` set to `false`, now model the
    conductor surfaces on both faces of the sheet, with an admittance `2/Z` instead of
    `1/Z`.
    [PR 994](https://github.com/awslabs/palace/pull/994).
  - Interior boundaries which are to be cracked (with `"Crack": true`, or with boundary
    conditions which apply to each side separately) are an error where the mesh cannot be
    cracked along them: on a mesh with hexahedra, prisms, or pyramids refined by
    `config["Model"]["Refinement"]["Boxes"]` or `["Spheres"]`, which makes it nonconformal,
    and on a nonconformal input mesh, for example from a previous adaptive simulation,
    which is not already cracked along them. The mesh was previously left uncracked.
    [PR 994](https://github.com/awslabs/palace/pull/994).
  - Surface flux postprocessing with `"TwoSided": false` on interior boundaries which are no
    longer cracked gives the average of the fluxes on the two sides, instead of the sum of
    the fluxes of the two faces of the cracked mesh.
    [PR 994](https://github.com/awslabs/palace/pull/994).

#### Performance Improvements

  - Reduced matrix-free material and boundary operator work by omitting elements outside
    shared nonzero-coefficient support unions, including compatible complex operator
    packing.
    [PR 987](https://github.com/awslabs/palace/pull/987).
  - Accelerated adaptive online sweeps by replacing 2D wave-port EVP solves with a per-port
    PROM trained on the offline modes. Also reusing reduced operators and wave-port state
    across excitations, evaluating the wave-port modal correction and excitation of the
    PROM from port-space mode pairings instead of per-frequency assembly of mode vectors on
    the full mesh, evaluating default domain energies in reduced coordinates, and evaluating
    port quantities without unnecessary magnetic-field reconstruction when possible.
    [PR 909](https://github.com/awslabs/palace/pull/909).
  - Restricted the assembly of wave-port mode vectors (excitation and modal correction) to
    the port boundary elements instead of looping over all boundary elements of the mesh,
    and assembled the per-frequency wave-port mode forms (S-parameter projection and modal
    reactions) in a single sweep over the port elements.
    [PR 909](https://github.com/awslabs/palace/pull/909).
  - CSV output tables are now appended row by row instead of rewritten at every measurement,
    so long transient and single-excitation driven runs write O(N) rather than O(N²) bytes.
    [PR 942](https://github.com/awslabs/palace/pull/942).

#### Bug Fixes

  - Fixed the sign of the second-order absorbing boundary condition term (`"Absorbing"` with
    `"Order": 2`): it now matches β = μ⁻¹c₀/(2iω) in the documentation. Before, the extra
    term made the boundary reflect more than the first-order condition at oblique incidence.
    [PR 1070](https://github.com/awslabs/palace/pull/1070).
  - Fixed cracking of interior boundaries in meshes with quadrilateral faces, where a
    cracked face with three vertices on the crack front (at a corner of the crack) merged
    with its twin and broke parallel mesh construction.
    [PR 992](https://github.com/awslabs/palace/pull/992).
  - Fixed numeric wave ports on superconducting lines selecting an evanescent mode: the
    kinetic inductance makes the line modes slower than bulk waves, above the eigenvalue
    target of the port mode solve, which only sought modes below it. The port modes are
    now sought on both sides of the target. Added `Boundaries.WavePort[].Target`, the
    effective index near which the port modes are sought, to select a mode far from the
    default target. Wave ports warn when the selected mode is evanescent. SchemaVer 2-3-0
    [PR 1067](https://github.com/awslabs/palace/pull/1067).
  - Fixed the 2D magnetic field energy (`E_mag`, `p_mag`) of driven, eigenmode, and
    transient simulations and of boundary modes (the out-of-plane part), and the 2D error
    estimates of all problem types but electrostatics, which also drive adaptive mesh
    refinement. Both used the out-of-plane magnetic flux density scaled by the element
    Jacobian determinant, which could put `E_mag` off by orders of magnitude.
    [PR 1069](https://github.com/awslabs/palace/pull/1069).
  - Fixed the sign of `Im{Z}` in `port-Z.csv` for wave ports with a `"VoltagePath"`: the
    total-field impedance was computed as `|V|² / P` instead of `|V|² / conj(P)`, giving
    `conj(Z)`. `Re{Z}` and `Z_PV` are unchanged.
    [PR 1066](https://github.com/awslabs/palace/pull/1066).
  - Fixed the final ParaView save (mesh partition and error indicator) overwriting a saved
    field when sources are written out of index order, for example a magnetostatic
    `SurfaceCurrent` index above a `FluxLoop` index.
    [PR 979](https://github.com/awslabs/palace/pull/979).
  - Fixed boundary coefficient terms being added to attributes outside their boundary when
    attributes with equal properties shared one material entry (a term stamped per attribute
    or per port element was counted once per attribute on all of them). This affects every
    coefficient combining several boundary terms: the boundary mode (wave port) system
    matrices and the driven extra-system matrix (surface conductivity, rational impedance,
    wave ports), whose results change for configurations with several such attributes, and
    the 3D preconditioner, which now matches the intended matrix.
    [PR 962](https://github.com/awslabs/palace/pull/962).
  - Fixed the flux error estimators on meshes with several element geometry types, where
    every estimate after the first (and the imaginary part of complex fields) reused stale
    field data for all geometry types but the first.
    [PR 962](https://github.com/awslabs/palace/pull/962).
  - Fixed `"BoundaryMode"` simulations with p-multigrid preconditioning
    (`"MGMaxLevels"` > 1) aborting in the Chebyshev smoother setup. The Chebyshev and Jacobi
    smoother spectral estimates now accept any finite, nonzero diagonal, scaling by its
    absolute value (unchanged for positive diagonals).
    [PR 977](https://github.com/awslabs/palace/pull/977).
  - Fixed a bug with first-order absorbing boundary condition on the out-of-plane field
    component in `"BoundaryMode"` simulations and wave port mode solves, which sometimes
    used incorrect boundary material properties.
    [PR 977](https://github.com/awslabs/palace/pull/977).
  - Fixed `GitTag` in `palace.json` reporting `"UNKNOWN"` for builds without usable Git
    metadata. The version can now be supplied with the `PALACE_GIT_COMMIT_ID` CMake option;
    the Spack package sets it from the resolved commit (`v<version>` for releases, from
    0.19.0 on), and the containers embed the `git describe` of the built commit.
    [PR 974](https://github.com/awslabs/palace/pull/974).
  - Corrected the `FluxLoop` hole inductance, which was low in 0.18 (e.g. the
    `circular_hole` PEC limit moves from 2.00 to 2.85 pH): 0.18 pinned the film's tangential
    field to a regularized 2D surface-curl solution, whereas the film is now driven by a
    curl-free cut generator and relaxes freely. Re-run existing `FluxLoop` cases; the
    regression references were re-baselined.
    [PR 929](https://github.com/awslabs/palace/pull/929).
  - Fixed CG, GMRES, and FGMRES returning NaN when the initial residual is exactly zero
    (e.g. a zero right-hand side) and no absolute tolerance is set.
    [PR 941](https://github.com/awslabs/palace/pull/941).
  - Fixed driven-sweep restarts (`"Restart"` > 1) aborting with
    `Column ... not found in table`, because reloaded CSV tables kept stale column-name
    indices.
    [PR 944](https://github.com/awslabs/palace/pull/944).
  - Fixed the ordering of the initial guesses of the hybrid nonlinear eigenvalue solver
    when a linear eigenpair has a residual equal to the low-error threshold, for example
    when several residuals are exactly zero: the ordering was undefined, and could try most
    of the low-error guesses after the high-error ones.
    [PR 999](https://github.com/awslabs/palace/pull/999).
  - Fixed the sign of the London penetration depth term for the out-of-plane field component
    in `"BoundaryMode"` simulations and wave ports.
    [PR 1025](https://github.com/awslabs/palace/pull/1025).
  - Fixed issues affecting non-cracked 2D nonconformal meshes, and the wave ports of
    non-cracked 3D nonconformal meshes, where the edges on the unrefined side were not
    constrained. This patches MFEM with
    [mfem/mfem#5531](https://github.com/mfem/mfem/pull/5531).
    [PR 1023](https://github.com/awslabs/palace/pull/1023).
  - Fixed boundary mode simulations extracting their cross-section from a nonconformal 3D
    mesh (such as a mesh saved by an adaptive mesh refinement run), which lost all boundary
    conditions of the cross-section.
    [PR 1023](https://github.com/awslabs/palace/pull/1023).
  - Fixed a memory leak in the AMS preconditioner at order 1, including the lowest-order
    level of p-multigrid: the interpolation matrices that hypre builds from the vertex
    coordinates were never freed, so memory grew at every operator update, such as each
    frequency of a driven sweep (about 1.9 GiB per frequency for 3.9M unknowns on GPU).
    [PR 1038](https://github.com/awslabs/palace/pull/1038).
  - Fixed integer range expansion in configuration files rewriting brackets inside strings,
    so a path such as `"mesh[1-3].msh"` is no longer read as `"mesh[1,2,3].msh"`, in both
    *Palace* and `scripts/validate-config`.
    [PR 970](https://github.com/awslabs/palace/pull/970).
  - Fixed the electric surface flux (`surface-F.csv`) and boundary surface charge (`Q_s`)
    of complex fields in lossy dielectrics, which used the real permittivity only; they now
    use the complex permittivity `ε(1 − i tanδ)` of the material model, consistent with the
    solver. Also corrected the documented `Q_s` units to C/m².
    [PR 953](https://github.com/awslabs/palace/pull/953).
  - Fixed a heap overflow in 2D surface flux postprocessing with an explicit `"Center"`.
    [PR 953](https://github.com/awslabs/palace/pull/953).
  - Fixed a segmentation fault of the hybrid nonlinear eigenvalue solver when the
    quasi-Newton refinement converges none of the eigenpairs, which is now reported as an
    error since fewer eigenmodes than requested were found.
    [PR 1036](https://github.com/awslabs/palace/pull/1036).
  - Fixed numeric wave-port and boundary-mode mode ordering on lossy cross-sections (e.g.
    with absorbing or finite-conductivity boundaries) that could select a strongly
    evanescent mode instead of the intended propagating mode. Exact and reduced (adaptive
    sweep) mode solves now rank modes by the complex distance of the propagation constant
    from the target.
    [PR 1019](https://github.com/awslabs/palace/pull/1019).
  - Fixed crashes in the libCEED setup on nonconforming domain submeshes and in the flux
    error estimators on processes without elements of a submesh.
    [PR 1022](https://github.com/awslabs/palace/pull/1022).
  - Fixed silent misparsing of Nastran coordinates with implicit or Fortran `D` exponents
    (e.g. `-7.-1` was read as `-7.0` instead of `-0.7`), and read blank or
    trailing-trimmed fixed-width fields as zero instead of aborting.
    [PR 951](https://github.com/awslabs/palace/pull/951).
  - Fixed reading Nastran meshes exported by Gmsh, which have no `BEGIN BULK` line. Free
    field (comma-separated) cards that omit trailing fields no longer reuse the previous
    field (`CTRIA3,1,1,1,2,3` became a degenerate 6-node triangle), free field and
    blank-field continuation lines are read, and the final line no longer needs a newline.
    `GRID` points with a nonzero coordinate system (`CP`, or a `GRDSET` default) now abort,
    since coordinate systems are not applied.
    [PR 951](https://github.com/awslabs/palace/pull/951).
  - Fixed driven frequency sampling edge cases: `"NSample": 1` now samples `"MinFreq"`
    instead of producing a NaN frequency, and a `"FreqStep"` sweep with `"MaxFreq"` below
    `"MinFreq"`, a nonpositive step, or too many steps for an `int` now aborts instead of
    sampling only `"MinFreq"`, crashing, or sweeping below `"MinFreq"`. `validate-config`
    now also rejects circuit synthesis without a positive `"AdaptiveTol"`, as *Palace*
    already did, and the schema rejects 0 GHz samples, where the driven solve is singular.
    SchemaVer 2-2-1
    [PR 952](https://github.com/awslabs/palace/pull/952).
  - Fixed adaptive frequency sweeps hanging when all samples were nearly the same
    frequency, and circuit synthesis with a single distinct sample, which now aborts at
    startup instead of failing in the sweep.
    [PR 952](https://github.com/awslabs/palace/pull/952).
  - Fixed out-of-bounds writes for `"Impedance"` and `"Conductivity"` boundary attributes
    outside the range of the mesh boundary attributes (for example `0` or a typo). They are
    now ignored with a warning, like other unknown boundary attributes.
    [PR 1002](https://github.com/awslabs/palace/pull/1002).
  - Fixed coordinate-path line integrals (wave port `"VoltagePath"` and the `"Impedance"`
    and `"Voltage"` postprocessing paths in `"BoundaryMode"` simulations) losing accuracy
    on locally refined meshes: a few tenths of a percent in `Z_PV` on the CPW examples, and
    up to about 1% after several levels of adaptive refinement towards conductor edges,
    where `Z_PV` then changed erratically from one level to the next. They are now exact for
    the discrete field on affine elements, and the wave port voltage no longer needs GSLIB
    or a point search over the full 3D mesh.
    [PR 1049](https://github.com/awslabs/palace/pull/1049).
  - Fixed wave port modes depending on the mesh partitioning: for some partitionings, the
    finite element spaces of the wave port cross-section were inconsistent between
    processes, and the port eigenvalue problem gave spurious modes. This patches MFEM with
    [mfem/mfem#5639](https://github.com/mfem/mfem/pull/5639), which orders the shared
    entities of `ParSubMesh` consistently across processes.
    [PR 1056](https://github.com/awslabs/palace/pull/1056).
  - Fixed wave port boundary conditions after nonconforming adaptive mesh refinement in
    parallel: the boundary edges of the port cross-section could miss the attribute of the
    adjacent boundary when that boundary was owned by another process, so its boundary
    condition was not applied on the port. Likewise, an uncracked internal boundary
    crossing the port (now the default; before with
    `"CrackInternalBoundaryElements": false` or a nonconformal input mesh) could miss its
    boundary condition on the port in parallel.
    [PR 1056](https://github.com/awslabs/palace/pull/1056).
  - Fixed driven-sweep restarts (`"Restart"` > 1) with `"FarField"` postprocessing, which
    always aborted while reloading the CSV tables because the far-field table was never
    read back. It now reloads and validates like the other tables, matching its
    one-row-per-direction layout.
    [PR 1075](https://github.com/awslabs/palace/pull/1075).
  - Fixed interface dielectric energies on uncracked interior metal sheets (now the default;
    before with `config["Model"]["CrackInternalBoundaryElements"]` disabled), which used
    the energy of the average of the fields of both sides for type `"Default"`, and a single
    face for types `"MA"` and `"MS"` on sheets with the same material on both sides, such as
    air bridges. They now give the same results as on a cracked mesh.
    [PR 994](https://github.com/awslabs/palace/pull/994).
  - Fixed silently wrong results of multi-GPU simulations with GPU-aware MPI, such as wrong
    S-parameters and slow or failed GMRES convergence, caused by a race in the device
    matrix-vector product of hypre. The superbuild now patches hypre v3.2.0 with
    [hypre-space/hypre#1629](https://github.com/hypre-space/hypre/pull/1629); for Spack
    builds, the patch is in
    [spack/spack-packages#7198](https://github.com/spack/spack-packages/pull/7198).
    [PR 1062](https://github.com/awslabs/palace/pull/1062).

#### Documentation

  - The developer testing notes now describe `TestMeshPath` for the unit test meshes and
    `PALACE_TEST_DATA_DIR` instead of the removed `MFEM_DATA_PATH`, and
    `test/data/mesh/README.md` records the source and MFEM version of each test mesh.
    [PR 1073](https://github.com/awslabs/palace/pull/1073).
  - Fixed the sign in the denominators of the finite-thickness surface impedance formula in
    the reference docs, which now match the implementation.
    [PR 1074](https://github.com/awslabs/palace/pull/1074).
  - Documented in the developer notes that schema validation tightenings which only
    reject configurations that never produced valid solver output are bug fixes
    (ADDITION), not breaking (MODEL) changes.
    [PR 952](https://github.com/awslabs/palace/pull/952).
  - Restored the paragraph on the outputs of the lumped and wave port case of the coaxial
    example, which was cut off mid-equation.
    [PR 1018](https://github.com/awslabs/palace/pull/1018).
  - The installation guide and quick start now install Spack from the release series
    *Palace* is tested with (`releases/v1.2`) instead of Spack `develop`.
    [PR 990](https://github.com/awslabs/palace/pull/990).

#### Build system

  - *Palace* now requires a C++20 compiler, including for CUDA and HIP device code, so
    CUDA builds need CUDA 12.0. The Spack recipe keeps
    `cxxstd=17` available for 0.16 to 0.18.
    [PR 654](https://github.com/awslabs/palace/pull/654).
  - GPU Spack builds from 0.19 no longer cap Umpire at 2025.12; MFEM is built as C++20 to
    match. 0.16 to 0.18 keep Umpire 2025.12 as C++17.
    [PR 654](https://github.com/awslabs/palace/pull/654).
  - Bumped the MFEM dependency to v4.10 and dropped the backport patches merged
    upstream (PRs 4983, 5246, 5415, 5124, and the Gmsh reader rewrite), keeping
    only the still-unmerged patches (PRs 3847, 5353, 5494, 5502).
    [PR 918](https://github.com/awslabs/palace/pull/918).
  - Bumped the MUMPS dependency of the superbuild to 5.9.1. The superbuild now uses the
    ParMETIS and METIS it builds instead of letting MUMPS download its own.
    [PR 1026](https://github.com/awslabs/palace/pull/1026).
  - A `BLA_VENDOR` passed to CMake now takes precedence over the BLAS/LAPACK vendor
    detected from the environment, and AOCL is located through CMake's `FindBLAS` and
    `FindLAPACK`. *Palace* now requires CMake 3.27 or later.
    [PR 1026](https://github.com/awslabs/palace/pull/1026).
  - `palace --version` now reports the version of every linked dependency, and the build
    identifier set with the `PALACE_BUILD_ID` CMake option (the Spack package sets it to
    its DAG hash).
    [PR 1028](https://github.com/awslabs/palace/pull/1028).
  - The superbuild dependencies are now pinned to release tags instead of commit SHAs and
    were bumped to their latest releases. The Spack package now uses a libCEED release
    instead of a development commit.
    [PR 1028](https://github.com/awslabs/palace/pull/1028).
  - The `cudss` variant of the Spack package now requires *Palace* 0.19 or later, which
    builds against MFEM 4.10: the MFEM `cudss` variant does not exist for the MFEM 4.9 of
    0.16 to 0.18.1, so `palace@0.18+cudss` could not be concretized.
    [PR 1053](https://github.com/awslabs/palace/pull/1053).
  - Added an MFEM patch for [mfem/mfem#5698](https://github.com/mfem/mfem/pull/5698)
    (`NCMesh::GetRootTransforms`), used to identify the interior boundary sides of
    nonconforming meshes.
    [PR 994](https://github.com/awslabs/palace/pull/994).
  - Added an MFEM patch for [mfem/mfem#5610](https://github.com/mfem/mfem/pull/5610),
    which keeps host hypre matrix products from reading a device-resident row cache.
    [PR 992](https://github.com/awslabs/palace/pull/992).
  - Bumped the superbuild dependencies to their latest releases: libCEED v1.0.1,
    SuperLU_DIST v9.3.1, json-schema-validator 2.4.1, fmt 12.2.0, Umpire v2026.07.2, and
    Enzyme v0.0.303. ButterflyPACK stays at v4.1.0, as STRUMPACK 8.1.0 does not link
    against ButterflyPACK 5.0.0.
    [PR 1114](https://github.com/awslabs/palace/pull/1114).
  - Added an MFEM patch for [mfem/mfem#5714](https://github.com/mfem/mfem/pull/5714),
    which restores reading Gmsh meshes whose boundary elements have a different order than
    the volume elements, such as linear tetrahedra with second-order boundary triangles, as
    produced by Gmsh when optimizing a second-order mesh with Netgen. The Gmsh reader of
    MFEM 4.10 rejected them with "Variable order Gmsh meshes are not supported", whereas
    *Palace* 0.18 read them as linear meshes.
    [PR 1122](https://github.com/awslabs/palace/pull/1122).

## [0.18.1] - 2026-09-21

#### New Features

  - Added machine-readable schema compatibility data to the installed schema files and
    generated the developer documentation table from the same source
    [PR 881](https://github.com/awslabs/palace/pull/881).
  - Improve BoundaryMode linear solver convergence when lossy boundary conditions are present
    by including complex terms in the real-valued preconditioner, or using an exact
    complex-valued preconditioner when `"ComplexCoarseSolve"` is true. Numeric wave ports
    can override this choice per port with `"ComplexCoarseSolve"`, without changing the full
    3D system preconditioner. SchemaVer 1-7-0
    [PR 921](https://github.com/awslabs/palace/pull/921).

#### Interface Changes

  - The adapted mesh saved by `config["Model"]["Refinement"]["SaveAdaptMesh"]` is now
    gzip-compressed and written with a `.meshgz` extension when *Palace* is built with `zlib`
    support (otherwise it is written uncompressed as `.mesh`, as before). This mesh can be
    used directly in Palace. [PR 923](https://github.com/awslabs/palace/pull/923).

#### Performance Improvements

  - Reduced repeated work in complex operators [PR
    932](https://github.com/awslabs/palace/pull/932).
  - Reduced repeated work in PROM construction [PR
    930](https://github.com/awslabs/palace/pull/930).
  - Adaptive driven sweeps now buffer CSV output tables in memory during the online phase
    and rewrite them at most once every 10 s (and at the end of the sweep), instead of
    rewriting the growing tables at every output frequency.
    [PR 938](https://github.com/awslabs/palace/pull/938).

#### Bug Fixes

  - Fixed issues with nonconformal AMR corrupting the wave port mesh and VoltagePath line
    integral in 3D simulations. [PR 919](https://github.com/awslabs/palace/pull/919).
  - Fixed issues with reading from MFEM meshes (`.mesh`) deadlocking simulations
    or leading to incorrect results. [PR 927](https://github.com/awslabs/palace/pull/927).
  - Fixed numeric wave ports ignoring their configured `"MaxIts"` and `"KSPTol"`, so the
    port-mode linear solve inherited the tolerance of the full 3D linear solver instead.
    [PR 921](https://github.com/awslabs/palace/pull/921).
  - Fixed cuDSS builds with MPI implementations other than Open MPI. The cuDSS MPI
    communication layer is now built against the selected MPI, by the superbuild from the
    source in the cuDSS archive and by Spack via `cudss+mpi`, instead of loading NVIDIA's
    Open MPI binary. [PR 826](https://github.com/awslabs/palace/pull/826).
  - Fixed a reference-counter overflow in nonconformal meshes with high-valence vertices,
    which could assert or crash during adaptive mesh refinement.
    [PR 934](https://github.com/awslabs/palace/pull/934).
  - Fixed the units of the far-field output `farfield-rE.csv`, which was off by a factor
    `Lc / Z₀` and so depended on the characteristic length `Lc`.
    [PR 936](https://github.com/awslabs/palace/pull/936).
  - Fixed the quadrature data offset in the H(div) mass build QFunction, which left the cached
    curl-curl block of 2D-in-3D H(div) operators unwritten and partially overwrote the mass
    block. [PR 931](https://github.com/awslabs/palace/pull/931).
  - Fixed non-unitary wave-port S-parameters for modes with a longitudinal electric field
    (TM/hybrid modes), where a lossless shorted guide could report `|S11| > 1`. A modal
    correction to the port Robin operator restores unitarity and reciprocity across driven,
    adaptive, eigenmode, and circuit-synthesis solves; TEM/TE modes are unchanged.
    [PR 886](https://github.com/awslabs/palace/pull/886).
  - Fixed the SLP nonlinear eigensolver to fall back to the Hybrid method when SLEPc is
    unavailable or an ARPACK backend is selected (SLP is realized only through SLEPc), resolved
    before the Hybrid interpolation is built, and to handle a null damping matrix and a missing
    nonlinear term, preventing crashes on otherwise valid problems.
    [PR 886](https://github.com/awslabs/palace/pull/886).
  - Reported synthesized eigenvalue estimates with the same `Q = |ω|/(2|Im ω|)` convention as
    the eigenmode postprocessor, and report HDM absolute and backward residuals for every root
    instead of filtering by fixed quality-factor or coordinate-norm cutoffs.
    [PR 886](https://github.com/awslabs/palace/pull/886).

#### Build system

  - Improved PETSc and SLEPc CMake configuration diagnostics to include output from failed
    compile and runtime probes [PR 880](https://github.com/awslabs/palace/pull/880).
  - Pinned the libCEED commit in the Spack recipe per Palace release instead of building an
    unpinned libCEED `develop`. [PR 940](https://github.com/awslabs/palace/pull/940).

## [0.18.0] - 2026-09-09

#### New Features

  - Improved configuration schema validation errors to list valid enum values through nested
    schema alternatives and suggest the canonical capitalization for case-only mismatches
    [PR 890](https://github.com/awslabs/palace/pull/890).
  - Added libCEED-based GPU evaluation for domain and boundary visualization fields in
    ParaView and grid function output, preserving native primary-field spaces and
    discontinuous element-local derived fields
    [PR 825](https://github.com/awslabs/palace/pull/825).
  - Added libCEED-based GPU evaluation for boundary postprocessing functionals, avoiding
    host coefficient projection for energy, flux, current, port, and far-field reductions
    [PR 824](https://github.com/awslabs/palace/pull/824).
  - Added a `RationalImpedance` boundary condition: a surface (Robin) impedance boundary
    whose per-square impedance is an arbitrary rational function of frequency,
    `Zs(s) = N(s)/D(s)` with `s = iω`, given by numerator and denominator polynomial
    coefficients (its zeros and poles), generalizing the parallel-RLC `Impedance` boundary
    to any passive lumped-network response. SchemaVer 1-1-0
    [PR 770](https://github.com/awslabs/palace/pull/770).
  - Added support for the [NVIDIA cuDSS](https://docs.nvidia.com/cuda/cudss/index.html)
    high-performance CUDA sparse direct solver, selectable via
    `config["Solver"]["Linear"]["Type"] = "cuDSS"` for the KSP and wave port solvers.
    Enable at build time with `PALACE_WITH_CUDSS=ON` (requires `PALACE_WITH_CUDA=ON`).
    SchemaVer 1-3-0 [PR 717](https://github.com/awslabs/palace/pull/717).
  - Added control over how *inactive* surface current ports are treated during a
    magnetostatic inductance sweep, via `config["Solver"]["Magnetostatic"]["InactivePorts"]`
    (`"Open"`, the default, or `"Short"`). `"Open"` applies a natural (zero surface current)
    boundary condition so no induced current flows across an inactive port, while `"Short"`
    treats it as a perfect conductor (PEC) allowing induced screening currents. The mode can
    be overridden per port with `config["Boundaries"]["SurfaceCurrent"][...]["InactiveMode"]`.
    SchemaVer 1-5-0 [PR 831](https://github.com/awslabs/palace/pull/831).
  - Added flux-loop excitation support for magnetostatic simulations, enabling
    flux-trapping analysis with surface-curl postprocessing, flux-only
    inductance extraction, and corresponding examples.
    [PR 461](https://github.com/awslabs/palace/pull/461).
  - Added support for magnetostatic inductance extraction with combined `SurfaceCurrent` and
    `FluxLoop` excitations, which was previously rejected. The self-inductance blocks come
    from the usual cross-energies, while the current-flux mutual inductance is measured from
    the magnetic flux linked through each current port's aperture during the flux excitations,
    since the two excitation types are energy-orthogonal and their coupling does not appear in
    the cross-energy matrix. Every surface-current element in a mixed simulation specifies
    an oriented `"Aperture"` containing its spanning-surface `"Attributes"` and Cartesian
    `"Direction"`; multielement ports combine the element linkages using the same weights as
    their parallel current excitations. Mixed extraction requires all current ports to be
    open when inactive and always reports one complete bare inductance matrix. SchemaVer
    1-6-0 [PR 861](https://github.com/awslabs/palace/pull/861).
  - Improve hybrid nonlinear eigensolver convergence when using 2nd-order absorbing and
    rational impedance boundary conditions. [PR 792](https://github.com/awslabs/palace/pull/792).
  - Enable circuit synthesis with frequency-dependent boundary conditions, including
    wave ports with a new optional `IncludeInSynthesis` wave-port flag. SchemaVer 1-4-0
    [PR 811](https://github.com/awslabs/palace/pull/811).
  - Added a `"Problem"`/`"Iteration"` field to the `palace.json` metadata file, recording the
    1-based index of the most recently completed adaptive mesh refinement (AMR) iteration
    (`1` for the initial solve, matching the `iterationXX` archive subdirectory).
    [PR 906](https://github.com/awslabs/palace/pull/906).
  - Emit a `SavedAdaptedMesh` block to `palace.json` when
    `config["Model"]["Refinement"]["SaveAdaptMesh"]` is enabled and adaptation was performed,
    recording the true (conforming) topological entity counts of the adapted mesh broken down
    by element and face geometry (vertices, edges, faces per face geometry, and cells per
    element geometry), together with the sorted domain and boundary attributes. The block is
    written for every adapted mesh: the top-level `palace.json` describes the final mesh, and
    (with `SaveAdaptIterations` also enabled) each `iterationX/palace.json` describes that
    iteration's mesh, so each saved mesh is paired with a matching block.
    [PR 911](https://github.com/awslabs/palace/pull/911).

#### Interface Changes

  - Eigenmode simulations with frequency-dependent boundary conditions (e.g. waveports,
    second-order absorbing, surface conductivity) now evaluate the boundary conditions at
    the true complex frequency/eigenmode. Previously, the boundary conditions were
    evaluated using only the real component of the eigenmode.
    [PR 778](https://github.com/awslabs/palace/pull/778).
  - Enabled the `RationalImpedance` boundary condition in eigenmode simulations, handled
    by the nonlinear eigenvalue solver and evaluated at the complex eigenfrequency.
    SchemaVer 1-3-1 [PR 778](https://github.com/awslabs/palace/pull/778).
  - Excited lumped ports may now carry reactance (`"L"` and/or `"C"`, or the
    surface-parameter equivalents `"Ls"`/`"Cs"`, including purely reactive ports with
    `"R": 0` / `"Rs": 0`). The reactance enters the system matrix as a physical
    termination; the incident drive is normalized to the port `"R"` when positive, or to
    an internal unit reference impedance for a purely reactive port. Previously, excited
    lumped ports were required to be purely resistive. As part of this change, *passive*
    purely reactive lumped ports (`R = 0` with `L` and/or `C`) now report a finite
    S-parameter row referenced to the internal unit impedance, where previously their
    `port-S.csv` row was identically zero; similarly, the incident columns `inc<idx>` of
    `port-V.csv` and `port-I.csv` at a purely reactive excited port now report finite
    unit-reference values, where previously both were zero. These are bookkeeping
    quantities referenced to the internal unit impedance, not physical incident waves.
    [PR 841](https://github.com/awslabs/palace/pull/841).

#### Performance Improvements

  - Omit material terms with mathematically exact-zero coefficients from fine-level partial
    assembly while retaining the configured coarse sparse structure for symbolic reuse.
    [PR 876](https://github.com/awslabs/palace/pull/876).
  - Reused screened magnetostatic operators and AMS preconditioners across
    open-port excitations with the same inactive shorted-port set, avoiding
    repeated assembly and setup during multi-port inductance sweeps.
    [PR 879](https://github.com/awslabs/palace/pull/879).

#### Bug Fixes

  - Fixed a distributed NLEPS fallback crash when a rank has no essential true
    degrees of freedom and the nonlinear eigensolver exhausts its initial
    guesses. The solve now reports the unconverged result rather than
    dereferencing a null DoF list.
    [PR 922](https://github.com/awslabs/palace/pull/922).
  - Fixed `linalg::MatrixSqrt`/`MatrixPow` returning incorrect results for fully
    anisotropic material tensors (all three off-diagonal entries nonzero, e.g. a crystal
    with generically-rotated `MaterialAxes`). The closed-form 3x3 eigen-decomposition used
    an incorrect invariant and eigenvector expressions that break for degenerate spectra
    (such as rotated uniaxial tensors); this branch now uses MFEM's robust symmetric
    eigensolver. Affected `√(μ⁻¹ε)`, `(με)^(-1/2)`, `ε·√(I + tanδ·tanδᵀ)` in the material
    operator and the flux error-estimator weights.
    [PR 912](https://github.com/awslabs/palace/pull/912).
  - Fixed adaptive iteration output archiving overwriting earlier meshes and made its
    filesystem updates more robust. [PR 892](https://github.com/awslabs/palace/pull/892).
  - Fixed ParaView output for multiple driven excitations deleting fields from earlier
    excitations. [PR 892](https://github.com/awslabs/palace/pull/892).
  - Fixed a false-positive in the domain material-coverage check that could abort numeric
    wave-port simulations under nonconformal AMR on multiple MPI ranks. The check now
    validates only domains owned by a local volume element, ignoring ghost/neighbor
    bookkeeping attributes. [PR 888](https://github.com/awslabs/palace/pull/888).
  - Reject all-zero numerator or denominator polynomial coefficients in
    `config["Boundaries"]["RationalImpedance"]` during schema validation, matching the
    existing checks in the configuration parser and providing users with earlier feedback.
    SchemaVer 1-5-1 [PR 874](https://github.com/awslabs/palace/pull/874).
  - Fixed validation of domain material coverage so postprocessing and retained mesh
    attributes without a corresponding `config["Domains"]["Materials"]` entry are rejected
    instead of silently assigning zero material coefficients to retained volumes. [PR
    840](https://github.com/awslabs/palace/pull/840).
  - Fixed smoother spectral estimates to pass a Hermitian similarity operator to the HEP
    solver [PR 837](https://github.com/awslabs/palace/pull/837).
  - Fixed saving output to non-shared filesystems [PR 813](https://github.com/awslabs/palace/pull/813).
  - Fixed S-parameter post-processing for mixed Floquet + lumped/wave port configurations.
    Previously, `MeasureSParameter()` skipped processing when Floquet ports coexisted with
    lumped or wave ports, leaving `port-S.csv` with incorrect values. The fix applies a power
    normalization factor, restoring S-matrix unitarity and reciprocity for the mixed
    case [PR 804](https://github.com/awslabs/palace/pull/804).
  - Fixed a bug where a mesh saved with `"SaveAdaptMesh"` during nonconformal AMR on a
    periodic problem could not be used as the input mesh for a subsequent solve. Loading now
    detects the already-periodic mesh and skips the periodicity construction.
    [PR 810](https://github.com/awslabs/palace/pull/810).
  - Fixed an overflow in the mesh broadcast path used when loading nonconformal meshes by
    introducing a chunked MPI broadcast call handle arbitrarily large meshes.
    [PR 810](https://github.com/awslabs/palace/pull/810).
  - Fixed a bug in the 2D mode eigensolver that sometimes led to crashes when waveports
    were used in conjunction with non-zero conductivity materials.
    [PR 778](https://github.com/awslabs/palace/pull/778).
  - Fixed a bug in the 2D mode eigensolver for lossy cross-sections, where a complex
    permittivity term was missing. [PR 778](https://github.com/awslabs/palace/pull/778).
  - Fixed numeric wave port and automatic boundary-mode targeting for anisotropic
    materials, which could select an evanescent mode instead of the mode with the largest
    propagation constant. [PR 836](https://github.com/awslabs/palace/pull/836).
  - Fixed intermittent MPI deadlocks in output-directory unit tests when ranks observed
    filesystem changes at different times. [PR
    839](https://github.com/awslabs/palace/pull/839).

#### Documentation

  - Added a circuit-extraction guide covering the reduced-order circuit
    workflow, including the treatment of nonlinear frequency-dependent
    boundary terms.
    [PR 812](https://github.com/awslabs/palace/pull/812).
  - Corrected the configuration validation script name in the docs
    (`validate-config`). [PR 891](https://github.com/awslabs/palace/pull/891).
  - Clarified how to interpret the elapsed time report: indented rows are exclusive
    sub-categories, parent rows are remainders, and Min/Max/Avg are over MPI ranks.
    [PR 885](https://github.com/awslabs/palace/pull/885).

## [0.17.0] - 2026-06-28

#### New Features

  - Added a new feature guide for adaptive frequency sweeps in driven simulations, with
    CPW lumped-port and transmon examples, plus reference documentation for the adaptive
    PROM/MRI algorithm and practical guidance on tolerances, convergence, and validation
    [PR 702](https://github.com/awslabs/palace/pull/702).
  - Added a
    [SchemaVer](https://docs.snowplow.io/docs/pipeline-components-and-non-pipeline-components/iglu/common-architecture/schemaver/)
    version to the root configuration JSON schema via the standard `"$id"` field
    (`"urn:palace:schema:1-0-0"`). The schema contract is now versioned independently of
    the *Palace* release so downstream tools can reason about configuration
    compatibility. See `docs/src/developer/notes.md` for details. [PR
    765](https://github.com/awslabs/palace/pull/765).
  - Added `"IncludeInSynthesis"` boolean flag to lumped port configuration (default
    `true`). When adaptive driven circuit synthesis is enabled
    (`config["Solver"]["Driven"]["AdaptiveCircuitSynthesis"]`), this flag controls
    whether the port contributes a port-mode basis vector to the reduced-order model.
    The boundary condition is still enforced for `"IncludeInSynthesis": false` ports —
    only the row/column in the synthesized circuit matrices is omitted. Useful for
    keeping passive terminations in the simulation for correct physics without paying
    their basis-vector cost. Excited ports must always be included; the parser rejects
    configurations otherwise. [PR 741](https://github.com/awslabs/palace/pull/741).
  - Added `BoundaryMode` solver for mode analysis on 2D meshes, supporting both
    standalone 2D meshes and 2D boundary submesh extracted from a 3D mesh [PR
    657](https://github.com/awslabs/palace/pull/657).
  - Added support for 2D meshes in all simulation types [PR
    657](https://github.com/awslabs/palace/pull/657).
  - Changed waveport formulation to enable support for impedance boundary conditions [PR
    657](https://github.com/awslabs/palace/pull/657).
  - Enabled voltage and impedance postprocessing for waveport and BoundaryMode by
    specifying a line from ground to conductor, either through mesh coordinates or
    boundary attribute in `config["Boundaries"]["Postprocessing"]["Impedance"]` [PR
    657](https://github.com/awslabs/palace/pull/657).
  - Adaptive mesh refinement is now supported for `BoundaryMode` simulations on a 2D
    submesh extracted from a 3D input mesh. [PR
    727](https://github.com/awslabs/palace/pull/727).
  - Palace now writes a fully-resolved copy of the run configuration to
    `<output_dir>/<input>_resolved.json` at startup, with every implicit default and
    `"Default"` sentinel filled in (options delegated to an external library, such as
    `"ColumnOrdering"` and `"PCSide"`, remain `"Default"`). The sidecar passes schema
    validation and can be used to re-run the same simulation deterministically; it can
    also be printed without running via `--dry-run` [PR
    719](https://github.com/awslabs/palace/pull/719).
  - Enable S-parameters calculations when both lumped and wave ports are present [PR
    743](https://github.com/awslabs/palace/pull/743).
  - Added Floquet port boundary conditions for frequency-domain driven simulations on
    periodic structures, with support for oblique incidence, multi-order diffraction,
    TE/TM/circular polarization, and adaptive frequency sweeps. [PR
    693](https://github.com/awslabs/palace/pull/693).
  - Upgraded documentation of Palace JSON config. Added script to generate config page
    automatically from fields in the JSON schema. Enables linking to individual tags and
    better html visuals. [PR 716](https://github.com/awslabs/palace/pull/716)
  - Palace now prints mesh statistics after conformal adaptive refinement, making mesh
    quality changes easier to track during AMR runs.
    [PR 750](https://github.com/awslabs/palace/pull/750).

#### Interface Changes

  - The regression test suite has moved from the Julia harness in
    `test/examples/runtests.jl` to Catch2 cases registered in the
    `palace-unit-tests` binary. Run with
    `mpirun -np $NUM_PROC_TEST palace-unit-tests "[Regression]"` or via
    CTest as `ctest -R "^regression-"`. Regression fixtures now live with
    the rest of the installed test data under `test/data/regression/`:
    configs and meshes in `input/`, reference CSVs in `ref/`. Shared example
    inputs may be symlinked in the source tree, but installed test data is
    self-contained.
    Re-baselining is now `palace -np N test/data/regression/input/<case>/<config>.json`
    followed by copying the live `postpro/` CSVs into the matching
    `test/data/regression/ref/<case>/<subdir>/` directory.
    [PR 735](https://github.com/awslabs/palace/pull/735).
  - Box and sphere region refinement (specified under
    `config["Model"]["Refinement"]["Boxes"]` and
    `config["Model"]["Refinement"]["Spheres"]`) is now applied to the mesh before
    partitioning rather than after, and is no longer reflected as a distinct level of
    the geometric multigrid mesh hierarchy used when
    `config["Solver"]["Linear"]["MGUseMesh"]` is enabled. The finest solve mesh is
    unchanged. [PR 727](https://github.com/awslabs/palace/pull/727).
  - Fixed a bug where `farfield-rE.csv` did not distinguish per-excitation rows in
    driven simulations with multiple lumped-port excitations. The file now includes an
    `exc` column so each (frequency, angle, excitation) row is uniquely identifiable.
    Existing single-excitation output gains the column too (a constant value), shifting
    downstream column indices by one. [PR
    746](https://github.com/awslabs/palace/pull/746).
  - `Solver.Eigenmode.MaxSize` now resolves to a concrete value at parse time, using
    SLEPc's default `max(2 * N, N + 15)` for both SLEPc and ARPACK (previously the unset
    default was backend-dependent). `Solver.Eigenmode.MaxIts` now defaults to
    `1'000'000` [PR 719](https://github.com/awslabs/palace/pull/719).

#### Performance Improvements

  - Cached the GSLIB point locator used for wave-port voltage-path line integrals,
    avoiding repeated setup work during adaptive sweeps with wave ports while
    preserving results.
    [PR 768](https://github.com/awslabs/palace/pull/768).

#### Bug Fixes

  - Fixed the driven solver ParaView output writing the nondimensional angular frequency
    to the `Time` field instead of the physical frequency in GHz. [PR
    754](https://github.com/awslabs/palace/pull/754).
  - Fixed a contradiction between the configuration schema and runtime behavior where
    omitting `config["Problem"]["Output"]` (allowed by the schema) aborted at startup
    with "Invalid output directory, got empty string". `Output` now defaults to
    `"postpro"` when omitted [PR 753](https://github.com/awslabs/palace/pull/753).
  - Fixed impedance boundary mesh-cracking scaling to apply per-attribute rather than
    per-block, allowing mixed cracked and uncracked attributes within a single impedance
    boundary definition. [PR 748](https://github.com/awslabs/palace/pull/748).
  - Revert change of adaptive PROM from interpolation back to projection. [PR
    733](https://github.com/awslabs/palace/pull/733).
  - Fixed nondeterministic current-dipole assembly when the source lies on a shared mesh
    entity by distributing the underlying MFEM delta source over all containing
    elements. [PR 742](https://github.com/awslabs/palace/pull/742).
  - Fixed a bug where nonconformal AMR would lead to erroneous results when waveports
    are present. Part of [PR 657](https://github.com/awslabs/palace/pull/657).
  - For `BoundaryMode` simulations on a 3D mesh, the non-dimensionalization reference
    length `Lc` is now computed from the extracted 2D solve mesh rather than the loaded
    3D mesh, so reported scales match the geometry actually being solved. [PR
    727](https://github.com/awslabs/palace/pull/727).
  - Fixed a bug where `BoundaryMode` simulations on a 3D mesh could abort during
    non-dimensionalization in parallel runs. The 2D submesh extraction was computing its
    tangent-plane normal against a stale attribute list, collapsing the normal and the
    reference length `Lc` to zero. Pre-existing from [PR
    657](https://github.com/awslabs/palace/pull/657). [PR
    727](https://github.com/awslabs/palace/pull/727).
  - Fixed a bug causing incorrect S-parameters when lumped port pairs have different R
    values [PR 743](https://github.com/awslabs/palace/pull/743).
  - Fixed a bug where mesh cracking refinement (`RefineCrackElements`) would fail on
    periodic meshes [PR 777](https://github.com/awslabs/palace/pull/777).
  - Fixed a bug where the magnetic energy field calculation in 2D simulations used the
    incorrect permeability vector components [PR
    782](https://github.com/awslabs/palace/pull/782).

#### Documentation

  - Corrected second-order absorbing boundary condition documentation and examples to
    match the code behavior.
    [PR 738](https://github.com/awslabs/palace/pull/738).

#### Build system

  - Fixed the ButterflyPACK superbuild, updated ButterflyPACK to v4.1.0, and enabled
    ButterflyPACK coverage in CI.
    [PR 667](https://github.com/awslabs/palace/pull/667).

## [0.16.1] - 2026-04-24

#### New Features

  - Added a new table reporting memory consumption for various stages in the
    simulation. [PR 708](https://github.com/awslabs/palace/pull/708)

  - Expanded JSON schema validation to cover required fields, mutual exclusion constraints
    (e.g., `PEC`/`Ground`, `PMC`/`ZeroCharge`), array type validation, and numeric bounds.
    Many runtime checks are now caught earlier at configuration parsing time with clearer
    error messages [PR 635](https://github.com/awslabs/palace/pull/635).

#### Bug Fixes

  - Reduced memory usage when `MaxIts` for GMRES is larger than the number of
    required iterations. [PR 715](https://github.com/awslabs/palace/pull/715)
  - Improved IO performance for simulations with Adaptive Mesh Refinement. Now,
    files from previous iterations are moved instead of being copied. The latest
    output is always available at the top level of the output directory (as
    symlinks to the latest `iterationXX` folder). [PR
    673](https://github.com/awslabs/palace/pull/673)
  - Fixed port excitation parsing to correctly handle signed integer values
    [PR 635](https://github.com/awslabs/palace/pull/635).
  - Fixed linear solver stalling in the nonlinear eigenvalue solver by adding an absolute
    tolerance and adapting the relative tolerance to the outer residual
    [PR 694](https://github.com/awslabs/palace/pull/694).
  - Fixed divergence of the Quasi-Newton nonlinear eigensolver by adding an
    Armijo line search to damp overshooting Newton steps
    [PR 713](https://github.com/awslabs/palace/pull/713).
  - Replaced boundary mass matrix CG solve workaround for lumped port excitation
    projection with direct interpolation via MFEM's
    `ProjectBdrCoefficientTangent`. An MFEM patch fixes the underlying face DOF
    orientation bug for Nedelec elements at order >= 2 in parallel, removing the
    need for the expensive mass solve
    [PR 684](https://github.com/awslabs/palace/pull/684).
  - Fixed bug when using a `SurfaceCurrent` in a driven simulation without a pre-existing excitation
    [PR 661](https://github.com/awslabs/palace/pull/661).
  - Added warning when STRUMPACK is selected for magnetostatic problems, which
    will fail due to the singular curl-curl operator
    [PR 680](https://github.com/awslabs/palace/pull/680).
  - Embedded JSON schemas are now regenerated at build time instead of configure time,
    so schema edits no longer require a full CMake reconfigure
    [PR 679](https://github.com/awslabs/palace/pull/679).
  - Backport positive-weight simplex quadrature rules from MFEM
    ([mfem/mfem#5246](https://github.com/mfem/mfem/pull/5246))
    [PR 706](https://github.com/awslabs/palace/pull/706).

## [0.16.0] - 2026-03-05

#### New Features

  - Added `PALACE_BUILD_WITH_SANITIZERS` CMake option and `asan` Spack variant
    for AddressSanitizer and UndefinedBehaviorSanitizer support. [PR
    642](https://github.com/awslabs/palace/pull/642)
  - Added reporting for memory used by Palace. Peak values are also stored in
    `postpro/palace.json`. [PR 629](https://github.com/awslabs/palace/pull/629)
  - Added JSON schema validation for configuration files at runtime, with embedded schemas
    and detailed error messages including valid enum options
    [PR 613](https://github.com/awslabs/palace/pull/613).
  - Added support for electric current dipole source excitations. Current dipoles can be
    configured as an array in `config["Domains"]["CurrentDipole"]` with specifications for
    `"Index"`, `"Moment"`, `"Center"`, and `"Direction"`. Consult the
    [documentation](https://awslabs.github.io/palace/dev/examples/antenna/#alternative-short-dipole-using-the-electric-current-dipole-operator)
    for additional information. [PR 556](https://github.com/awslabs/palace/pull/556).
  - Added an option to specify the number of AMS iterations in `config["Solver"]["Linear"]["AMSMaxIts"]`,
    with the default value changed from 1 to the solution order, and changed some of the AMS options.
    [PR 625](https://github.com/awslabs/palace/pull/625).
  - Added support for circuit synthesis from the ROM in adaptive driven simulations
    (`config["Solver"]["Driven"]["AdaptiveCircuitSynthesis"]`). Modifies ROM to add `LumpedPort`
    fields to make circuit. Prints out the projected matrices of the driven problem as well as the
    basis orthogonalization matrix. Extended tests on aspects related to ROM and LumpedPortOp. [PR 326](https://github.com/awslabs/palace/pull/326).
  - Upgraded internal orthogonalization routines to optionally allow taking non-identity weights in
    inner product. Part of [PR 326](https://github.com/awslabs/palace/pull/326).

#### Bug Fixes

  - Fixed a bug in the computation of the waveport maximum propagation constant in parallel simulations
    [PR 580](https://github.com/awslabs/palace/pull/580).
  - Fixed a bug where `BdrFieldCoefficient::Eval` used a comma operator instead of addition
    when computing the average of field values on interior boundaries, causing incorrect
    results for interior trace postprocessing.
    [PR 621](https://github.com/awslabs/palace/pull/621).
  - Fixed a bug in complex diagonal operator multiplication where the imaginary part used an
    incorrect term in `AddMult` and `AddMultHermitianTranspose`.
    [PR 633](https://github.com/awslabs/palace/pull/633).

#### Interface Changes

  - `config["Boundaries"]["Terminal"]` must now be used for electrostatic simulations to
    compute capacitance matrices, and `config["Boundaries"]["SurfaceCurrent"]` for
    magnetostatic simulations to compute inductance matrices. Previously, terminals and
    surface current sources could also be specified under `"LumpedPort"`.
    [PR 624](https://github.com/awslabs/palace/pull/624).
  - `config["Solver"]["Eigenmode"]["Target"]` is now always required for eigenmode
    simulations. Previously, it could be omitted if a `"Driven"` section was also present.
    [PR 624](https://github.com/awslabs/palace/pull/624).

## [0.15.0] - 2025-12-2

#### New Features

  - Added support for nonlinear eigenvalue problems arising from frequency-dependent
    boundary conditions. Two nonlinear eigensolvers are now available and can be specified
    by setting the `config["Solver"]["Eigenmode"]["NonlinearType"]` option to `"Hybrid"`
    (default) or `"SLP"`. The nonlinear eigensolver will automatically be used if
    frequency-dependent boundary conditions are used. [PR
    467](https://github.com/awslabs/palace/pull/467).
  - Added support for extraction of electric fields in the radiative zone.
    Consult the
    [documentation](https://awslabs.github.io/palace/dev/features/farfield) for
    additional information. [PR
    449](https://github.com/awslabs/palace/pull/449).
  - Added support for more granular tests and for computing the testing coverage.
    Consult the
    [documentation](https://awslabs.github.io/palace/dev/developer/testing/) for
    additional information. [PR
    398](https://github.com/awslabs/palace/pull/398), [PR
    480](https://github.com/awslabs/palace/pull/480).
  - Added option to export field solutions as MFEM grid functions to visualize with
    [GLVis](https://glvis.org). Output formats can be specified in `config["Problem"]["OutputFormats"]`.
    [PR 518](https://github.com/awslabs/palace/pull/518).
  - Added an option to drop small entries (below machine epsilon) from the matrix used in the sparse
    direct solver. This can be specified with `config["Solver"]["Linear"]["DropSmallEntries"]`.
    [PR 476](https://github.com/awslabs/palace/pull/476).

#### Interface Changes

  - `config["Boundaries"]["Periodic"]` is now a dictionary where all periodic boundary
    pairs, built into the mesh file or not, should be specified in
    `config["Boundaries"]["Periodic"]["BoundaryPairs"]` and a single global Floquet wave
    vector can be specified in `config["Boundaries"]["Periodic"]["FloquetWaveVector"]`. [PR
    471](https://github.com/awslabs/palace/pull/471).
  - Changed Paraview and GLVis output fields from nondimensional to SI units
    [PR 532](https://github.com/awslabs/palace/pull/532).
  - The name of the unit test executable is now `palace-unit-tests` (instead of `unit-tests`)
    [PR 549](https://github.com/awslabs/palace/pull/549).

#### Bug Fixes

  - Changed wave port eigenproblem shift and sorting to fix an issue with the mode ordering.
    The first mode now has the largest propagation constant, closest to the TEM limit, and
    subsequent modes are ordered by decreasing propagation constant. [PR
    448](https://github.com/awslabs/palace/pull/448).
  - Fixed an issue where Gmsh meshes with built-in periodicity (specified in the mesh file)
    were failing. [PR
    471](https://github.com/awslabs/palace/pull/471).
  - Fixed bug where a mesh from a previous nonconformal adaptation could not be loaded to
    use in a non-amr simulation. [PR
    497](https://github.com/awslabs/palace/pull/497).
  - Fixed bug where `"CrackInternalBoundaryElements"` would result in incorrect results for
    some lumped port boundary conditions. [PR
    505](https://github.com/awslabs/palace/pull/505).
  - Changed the sign of the Floquet phase factor from $\exp(+ik\cdot x)$ to $\exp(-ik\cdot
    x)$. [PR 510](https://github.com/awslabs/palace/pull/510) and [PR
    526](https://github.com/awslabs/palace/pull/526).
  - Update EM constants to CODATA Recommended Values of the Fundamental Physical Constants
    2022 [PR 525](https://github.com/awslabs/palace/pull/525).
  - Scaled Rs/Ls/Cs of impedance boundary conditions affected by mesh cracking, fixing bug
    where `"CrackInternalBoundaryElements"` would lead to incorrect results. [PR
    544](https://github.com/awslabs/palace/pull/544).
  - Fixed Paraview/MFEM output for transient simulations. [PR 561](https://github.com/awslabs/palace/pull/561).
  - Fixed units of energy. Results were incorrectly reported in joules, while the correct
    unit was nanojoules. Results are now correctly in joules. [PR
    541](https://github.com/awslabs/palace/pull/541).
  - Fixed a bug in transient solver implicit formulation affecting the magnetic flux density fields computed with
    `"GeneralizedAlpha"` or `"RungeKutta"` transient solver types (`config["Solver"]["Transient"]["Type"]`)
    [PR 568](https://github.com/awslabs/palace/pull/568).

#### Build system

  - When building on a machine with a GPU, `CMAKE_CUDA_ARCHITECTURES` has now a
    default value that is automatically determined to be the consistent with the
    GPU (when no GPU is detected, `CMAKE_CUDA_ARCHITECTURES` is set to 70). [PR
    530](https://github.com/awslabs/palace/pull/530).
  - The spack Palace build now supports compilation of unit tests with Spack
    tests. [PR 549](https://github.com/awslabs/palace/pull/549).
  - MFEM is no longer unconditionally compiled alongside of Palace and can now
    be built separately. This is used in the Palace spack build, which now
    explicitly depends on MFEM. [PR
    550](https://github.com/awslabs/palace/pull/550).

## [0.14.0] - 2025-08-20

  - Added `--version` command line flag for displaying Palace version information.
  - Fixed a small regression bug for boundary postprocessing when specifying
    `"Side": "LargerRefractiveIndex"`, introduced as part of v0.13.0.
  - Added an improvement to numeric wave ports to avoid targeting evanescent modes at
    higher operating frequencies. Also finite conductivity boundaries
    (`config["Boundaries"]["Conductivity"]`) are automatically marked as PEC for the wave
    port mode solve (previously these were marked as PMC unless specified under
    `"WavePortPEC"`).
  - Fixed a bug in divergence-free projection for problems without essential or mixed
    boundary conditions.
  - Added `"MakeSimplex"` and `"MakeHexahedral"` mesh options to convert an input mesh to
    all tetrahedra or all hexahedra. Also adds `"SerialUniformLevels"` option to
    `config["Model"]["Refinement"]` for testing or debugging.
  - Added `config["Model"]["CrackInternalBoundaryElements"]` which will separate or "crack" the mesh
    along all internal boundaries. This improves the performance of error estimation and AMR
    as the recovered smooth fields do not enforce additional erroneous continuity at
    internal boundaries. This will change the default behaviour in the case of internal
    impedance boundary conditions, and can be disabled by setting this option to false.
  - Added support for exact periodic boundary conditions, these can be specified as part of
    the mesh file (where supported by the format) or by specification of `"DonorAttributes"`
    and `"ReceiverAttributes"` in `config["Boundaries"]["Periodic"]` which will attempt to
    match the mesh on the boundaries specified by the donor and receiver attributes. This is
    only possible the mesh on the donor and receiver match exactly, non-matching meshes are
    not supported.
  - Exposed configuration linear solver and eigen solver options for the wave port
    subproblem. These can now be specified as part of the `config["Boundaries"]["WavePort"]`
    configuration. The defaults align with the previously hardcoded values.
  - Nonconformal adaptation is now supported for WavePort boundary conditions. This was
    achieved through a patch applied to MFEM to support `mfem::ParSubMesh` on external
    nonconformal surface subdomains.
  - Added adaptive time-stepping capability for transient simulations. The new ODE integrators
    rely on the SUNDIALS library and can be specified by setting the
    `config["Solver"]["Transient"]["Type"]` option to `"CVODE"` or `"ARKODE"`.
  - Added an option to use the complex-valued system matrix for the coarse level solve (sparse
    direct solve) instead of the real-valued approximation. This can be specified with
    `config["Solver"]["Linear"]["ComplexCoarseSolve"]`.
  - Fix bug in London equation implementation where a curl-curl term was added to the
    stiffness operator instead of a mass term.
  - Added support for Floquet periodic boundary conditions with phase-delay constraints.
    The Floquet wave vector can be specified along with periodic boundaries in the
    `config["Boundaries"]["Periodic"]` configuration.
  - Fixed bug where the default `config["Model"]["CrackDisplacementFactor"]` would cause
    lumped port discovery to fail by reducing default from `1e-3` to `1e-12`.
  - Fixed bug in Nastran mesh reader where carriage returns (`\r`) in the mesh file could
    cause a failure to read the mesh.
  - Changed post-processing, so that appropriate measurements are always written to disk as CSV
    files. This is a breaking change, since users are now required to specify a valid output
    folder in `config["Problem"]["Output"]`. Previously an empty string `""` would suppress file
    printing. ParaView printing is still controlled and suppressed by the `"Save"` or `"SaveStep"`
    options in the `config["Solver"][...]`.
  - Changed measurement and printing during post-processing substantially. This change is not
    user-facing, but will enable future multi-excitation support. All measurements are now
    performed as part of the `PostOperator` class, which is templated on solver type. The public
    interface of the `PostOperator` has been simplified: measurements can no longer be called
    individually, but are grouped in a `MeasurePrintAll` function. The actual printing of
    measurements has moved out of individual solvers and is orchestrated by the `PostOperatorCSV`
    class. Other helper classes include a small `Table` class (for data storage and formatting)
    and a light wrapper `TableWithCSVFile` (with file interaction).
  - Changed unit conversion interface: this was moved out of `IOData` into a separate `Units` class.
  - Changed csv output to print more digits.
  - Changed boundary indexing to enforce that the value of `Index` is unique between `LumpedPort`,
    `WavePort`, `SurfaceCurrent` and `Terminal`. This is a breaking change — e.g. currently a lumped
    port and a wave port could share an index.
  - Added multi-excitation support. `Excitation` in `LumpedPort` and `WavePort` boundary conditions
    can now be specified either as bools or excitation integers. Allow driven simulations to have
    multiple excitations indices in same configuration file to be simulated consecutively during the
    same Palace run. The measurement output of multiple excitations is printed to the same csv files
    with distinguish columns being post-indexed by the excitation index.
  - Changed `RomOperator` to split out rational interpolation component into the separate
    `MinimalRationalInterpolation` class to support multiple excitations. The MRI is unique to each
    excitation, but the PROM is shared between excitations.
  - Added `PortExcitations` to manage excitation pattern and print it to json metadata.
  - Update spack configuration to address `0.12` API breakage between libxsmm and libCEED.
  - Support spack vended libCEED and GSLIB builds.
  - Fix bug where `MakeSimplex` would fail for higher order meshes.
  - Fix bug in STRUMPACK build where dependency SCALAPACK library install directory could
    change based on build environment.
  - Introduce `"Samples"` mechanism for driven simulation, which allows for specifying a
    range of sample frequencies by combining `"Linear"`, `"Log"` and `"Point"`
    specifications. These are joined together effectively allowing for variable resolution
    sampling (for examples a coarse sampling supported by a fine sampling in a subset
    region). This is in addition to the prior existing `"MinFreq"`, `"MaxFreq"`,
    `"FreqStep"` interface which is now implicitly converted to a `"Linear"` sample scheme
    internally.
  - Introduce `"Save"` keyword for `"Driven"` simulation type, which will save specific
    frequency choices as opposed to the existing (and maintained `"SaveStep"` based on a
    regular sampling). Only frequencies contained within the set of `"Samples"` are
    supported, as no interpolation is performed.
  - Fix bugs in post processing where `E_cap`, `E_ind`, and `mode_port_kappa` and other
    dependent quantities were not dimensionalized correctly.
  - Refactor `PostOperator` usage of `Measurement` to be entirely non-dimensional, until
    `PostOperatorCSV` which dimensionalizes all measurements before writing to file. Reduces
    the risk of mixed unit bugs throughout `PostOperator`.
  - Fix bug in using `"MakeSimplex"` which would cause undefined behaviour for higher order
    meshes.
  - Fix bug when combining OpenMP and GPU builds in reduction operations over `Vector`.
  - Fix race condition that would affect OpenMP parallelism with periodic boundaries (exact
    and Floquet).
  - Fix race condition in `mfem::DenseTensor::operator()` with OpenMP, due to class member
    variable access.
  - Fix race condition in `DofToQuad` methods within mfem.
  - Normalize eigenmodes so their mean phase is a positive real number.
  - Added `scnlib` as a dependency to Palace.
  - Added parsing ability of existing csv files into the `Table` class. Fix `"Restart"`
    behaviour with the multi-excitation feature of driven solver.

## [0.13.0] - 2024-05-20

  - Changed default value of `config["Solver"]["PartialAssemblyOrder"]` in order to activate
    operator partial assembly by default for all operators in all simulation types.
  - Changed the normalization of computed eigenmodes for consistency across different domain
    decompositions. Eigenvectors are now normalized with respect to the mass matrix for unit
    domain electric field energy.
  - Added documentation for various timer categories and improved timing breakdown of
    various sections of a simulation.
  - Changed mesh files for the cavity and CPW examples, including prism, hexahedral, and
    tetrahedral meshes for the cylindrical cavity and correcting the wave port dimensions
    for the coplanar wave guide.
  - Fixed a few bugs and issues in the implementation of numeric wave ports for driven
    simulations.
  - Added GPU support for *Palace* via its dependencies, and added the
    `config["Solver"]["Device"]` and `config["Solver"]["Backend"]` options for runtime
    configuration of the MFEM device (`"CPU"` or `"GPU"`) and libCEED backend, with suitable
    defaults for users.
  - Added a new section to the documentation on
    [Parallelism and GPU support](https://awslabs.github.io/palace/dev/guide/parallelism/).
  - Removed use of `mfem::SparseMatrix` and replaced with HYPRE's `hypre_CSRMatrix` when
    needed for full assembly, wrapped as `palace::hypre::HypreCSRMatrix`.
  - Added `"Active"` configuration file parameter for lumped and wave port boundaries to
    disable the associated boundary condition and only use the surface for postprocessing.
  - Changed the smooth flux space for the electrostatic error estimator to fix performance
    on problems with material interfaces.
  - Fixed error estimation bug affecting time-dependent simulation types (driven, transient,
    eigenmode) where the recovery of the electric flux density also needs to be taken into
    account in addition to the magnetic field.
  - Fixed a bug related to mesh cleaning for unspecified domains and mesh partitioning.
  - Change computation of domain energy postprocessing for electrostatic and magnetostatic
    simulation types in order to improve performance.
  - Fixed a bug when computing the energy associated with lumped elements with more than
    one nonzero R, L, or C. This also affects the inductive EPR for lumped inductors with
    and associated parallel capacitance.
  - Fixed a bug for coaxial lumped ports which led to incorrect extraction of the geometric
    parameters, especially when coarsely-meshed or non-axis-aligned.
  - Added boundary postprocessing functionality for surface flux including electric,
    magnetic, and power given by the Poynting vector. This results in some breaking changes
    to the configuration file interface, see
    `config["Boundaries"]["Postprocessing"]["SurfaceFlux"]` and
    `config["Boundaries"]["Postprocessing"]["Dielectric"]`. In addition, related
    configuration file keyword changes to for consistency were made to
    `config["Domains"]["Postprocessing"]["Probe"]` and
    `config["Model"]["Refinement"]["Boxes"]`.
  - Fixed a bug in MFEM for nonconformal AMR meshes with internal boundaries affecting
    non-homogeneous Dirichlet boundary conditions for electrostatic simulations (see
    [#236](https://github.com/awslabs/palace/pull/236)).

## [0.12.0] - 2023-12-21

  - Added support for operator partial assembly for high-order finite element spaces based
    on libCEED for mixed and non-tensor product element meshes. This option is disabled by
    default, but can be activated using `config["Solver"]["PartialAssemblyOrder"]` set to
    some number less than or equal to `"Order"`.
  - Added flux-based error estimation, reported in `error-estimate.csv`. This computes the
    difference between the numerical gradient (electrostatics) or curl (otherwise) of the
    solution, and a smoother approximation obtained through a global mass matrix inversion.
    The results are reported in `error-estimates.csv` within the `"Output"` folder.
  - Added Adaptive Mesh Refinement (AMR), specified in the `config["Model"]["Refinement"]`,
    for all problem types aside from transient. To enable AMR, a user must specify
    `"MaxIts"`, while all other options have reasonable defaults. Nonconformal (all mesh
    types) and conformal (simplex meshes) refinement are supported.
  - Added support for non-axis-aligned lumped ports and current sources. Key words `"X"`,
    `"Y"`, `"Z"` and `"R"`, with optional prefix `"+"` or `"-"` still work, but now
    directions can be specified as vectors with 3 components. Users will be warned, and
    ultimately errored, if the specified directions do not agree with axis directions
    discovered from the geometry.
  - Added output of lumped port voltage and current for eigenmode simulations.
  - Added dimensionalized output for energies, voltages, currents, and field values based on
    a choice of the characteristic magnetic field strength used for nondimensionalization.
  - Added output of electric and magnetic field energies and participation ratios in regions
    of the domain, specified with `config["Domains"]["Postprocessing"]["Energy"]` and
    written to `domain-E.csv`. This replaces
    `config["Domains"]["Postprocessing"]["Dielectric"]` and `domain-Q.csv`.
  - Added improved `Timer` and `BlockTimer` classes with more timing categories for
    reporting simulation runtime.
  - Changed implementation of complex-valued linear algebra to use new `ComplexVector` and
    `ComplexOperator` types, which are based on the underlying `mfem::Vector` and
    `mfem::Operator` classes, instead of PETSc. PETSc is now fully optional and only
    required when SLEPc eigenvalue solver support is requested. Krylov solvers for real- and
    complex-valued linear systems are implemented via the built-in `IterativeSolver`
    classes.
  - Changed implementation of PROMs for adaptive fast frequency sweep to use the Eigen
    library for sequential dense linear algebra.
  - Changed implementation of numeric wave ports to use MFEM's `SubMesh` functionality. As
    of [#3379](https://github.com/mfem/mfem/pull/3379) in MFEM, this has full ND and RT
    basis support. For now, support for nonconforming mesh boundaries is limited.
  - Added build dependencies on [libCEED](https://github.com/CEED/libCEED), including
    [LIBXSMM](https://github.com/libxsmm/libxsmm) and [MAGMA](https://icl.utk.edu/magma/)
    to support CPU- and GPU-based operator partial assembly.
  - Added unit test framework for all integrators based on
    [Catch2](https://github.com/catchorg/Catch2), which also includes some automated
    benchmarking capabilities for operator assembly and application.
  - Added improved OpenMP support in `palace` wrapper script and CI tests.
  - Added Apptainer/Singularity container build definition for *Palace*.
  - Fixed bugs related to thread-safety for OpenMP builds and parallel tetrahedral meshes in
    the upstream MFEM library.

## [0.11.2] - 2023-07-14

  - Fixed a regression bug affecting meshes which have domain elements which are not
    assigned material properties in the configuration file.
  - Changed layout and names of `palace/` source directory for better organization.
  - Added many updates to build system: Removed use of Git submodules to download
    dependencies relying instead directly on CMake's ExternalProject, patch GSLIB dependency
    for shared library builds, add CI tests with ARPACK-NG instead of SLEPc, update all
    dependency versions including MFEM to incorporate bug fixes and improvements. This
    affects the Spack package recipe, so a new recipe is distributed as part of *Palace* in
    in `spack/` which will keep compatibility with `main` while changes are upstreamed to
    the built-in Spack repository.
  - Added a basic Makefile with some useful targets for development.

## [0.11.1] - 2023-05-03

  - Fixed a bug for interface dielectric loss postprocessing including when using
    `--dry-run`.
  - Fixed a regression bug affecting second-order absorbing boundary conditions.
  - Fixed a bug when the number of processors was large enough that some subdomains own no
    true dofs.
  - Fixed memory bug in destruction of SuperLU solver.
  - Fixed some typos in the documentation.
  - Fixed a bug in the cavity convergence study example for the coarsest hexahedral mesh
    case, as well as some other minor Julia fixes in the `examples/` and `scripts/`
    directories.
  - Added updates to superbuild including better ARPACK support, though still experimental
    compared to SLEPc.
  - Added updated submodules for superbuild.

## [0.11.0] - 2023-01-26

  - Initial public release on GitHub.
  - Added support for anisotropic materials through the use of `"MaterialAxes"` under
    `"Materials"`. These material axes allow for specifying symmetric material property
    tensors through a sum of weighted outer products of normal vectors. This can be used in
    conjunction with scalar material properties or multiple different anisotropic
    properties, though all properties are required to use the same set of basis vectors.
    Demonstration given for the coplanar waveguide example which now utilizes a sapphire
    substrate.
  - Added to default postprocessing outputs: Bulk and interface dielectric loss writes
    energy-participation ratios in addition to quality factors associated with lossy
    regions, IO coupling $\kappa$ in addition to quality factor for eigenmode
    simulations, lumped element energy for transient simulations similar to eigenmode and
    frequency domain driven simulations.
  - Changed configuration file syntax to simplify support for multielement lumped ports,
    where each index for a multielement port should now use the `"Elements"` object to
    describe the attributes for each port element.
  - Changed configuration file keywords to better reflect meaning:
    `config["Boundaries"]["Current"]` => `config["Boundaries"]["SurfaceCurrent"]`, and
    `"UseGMG"`, `"UsePCShifted"`, `"MGCycleIts"`, and `"MGSmoothIts"` under
    `config["Solver"]["Linear"]`.
  - Changed geometric multigrid implementation to generalize across problem types (added
    for electrostatics) and use Jacobi-preconditioned Chebyshev smoothing with optional
    auxiliary space smoothing at each multigrid level as well. The auxiliary space matrices
    are constructed directly in H1 now to allow for matrix-free/partial-assembly support
    eventually. Geometric multigrid is turned on by default.
  - Added structured simulation metadata output in the form of a JSON file `palace.json`
    written to the directory specified by `config["Problem"]["Output"]`.
  - Added JSON Schema files for the configuration file format as well as a helper script
    to check a provided configuration file against the schema.
  - Added optional interface to GSLIB library which enables `"Probe"` functionality to
    sample the computed electric and magnetic fields at points in space.
  - Added preliminary support for builds using Spack. This includes an option in the build
    system for user-supplied dependencies as an alternative to the superbuild
    configuration.
  - Added updated submodules for superbuild, fixing a bug in SuperLU_DIST solver causing
    communication hangs for certain numbers of processors.

## [0.10.0] - 2022-10-04

  - Added interfaces to ARPACK and FEAST eigenvalue solvers.
  - Added divergence-free projection for eliminating spurious DC modes for eigenmode
    solves.
  - Added option for visualizing fields on mesh boundaries, output alongside the full 3D
    solution fields for ParaVeiew visualization.
  - Added real and imaginary fields output for complex-valued phasors, and electric and
    magnetic energy density fields.
  - Added convergence study for the cavity example application, and example Julia code for
    automated mesh generation, running of the solver, and postprocessing.
  - Fixed bugs in mesh preprocessing related to precision of nodal coordinates when
    distributing a serial mesh to each process.
  - Added option for `"Driven"` and `"Transient"` simulations to accelerate postprocessing
    by only considering port boundaries.
  - Added `-dry-run`/`--dry-run` command line option to check configuration file for errors
    and exit.
  - Changed dependency build process to rely less on PETSc build system, no longer give
    option to build BLAS/LAPACK libraries from source, added detection for OpenBLAS and AMD
    Optimizing CPU Libraries (AOCL) blis/libflame.
  - Fixed issues when compiling with Intel or Clang compilers.
  - Added CI test builds to test various build options and added easier default build
    options in CMake configuration.
  - Added `clang-format` and `JuliaFormatter` configurations for automated C++, Julia, and
    Markdown code formatting.

## [0.9.0] - 2022-07-11

  - Added region-based mesh refinement with box or sphere refinement regions.
  - Added new sparse direct solver interfaces to SuperLU_DIST and STRUMPACK libraries, to
    complement the existing MUMPS interface.
  - Changed configuration file keywords for linear solver and preconditioner parameters to
    make more options available to the user and clarify them.
  - Added check that all external boundary surface attributes must have an associated
    boundary condition specified in the configuration file (adds new `"PMC"` and
    `"ZeroCharge"` boundaries for natural boundary conditions).
  - Added advanced configuration options for GMRES/FGMRES/SLEPc eigenvalue solver
    orthogonalization, making CGS2 (classical GS with iterative refinement) the default.

## [0.8.0] - 2022-06-06

  - Added new comprehensive example applications and detailed tutorials in the
    documentation.
  - Added CI pipeline for build and regression testing.
  - Fixed bugs for shared library builds, which are now the default for dependency builds.
  - Added complex lumped port voltage and current output for computing impedance and
    admittance parameters from frequency domain driven simulations.
  - Fixed a bug in mesh parsers for COMSOL and Nastran formats for high-order hexahedra,
    prism, and pyramid element types.
  - Changed adaptive frequency sweep implementation to introduce optimizations which, in
    particular, speed up driven simulations with wave ports.
  - Fixed bug related to domain dielectric postprocessing and surface dielectric loss
    postprocessing now automatically detects correct side for substrate-air interfaces
    (always evaluates on vacuum side).

## [0.7.0] - 2022-03-25

  - Added adaptive fast frequency sweep implementation for frequency domain driven
    simulations, activated with `config["Solver"]["Driven"]["AdaptiveTol"] > 0.0`.
  - Changed ASCII file postprocessing to write out .csv files instead of generic whitespace
    delimited .txt for easier parsing.
  - Changed field output postprocessing for visualization to save the mesh in the original
    mesh length units (rather than nondimensionalized ones).
  - Added electric and magnetic field energy output to files as function of mode, frequency,
    time, or solution step, in `domain-E.csv`.
  - Added option to specify circuit R/L/C parameters for lumped port boundaries as an
    alternative to the surface Rs/Ls/Cs ones.
  - Added improved and expanded documentation.
  - Added new MFEM features, including support for ND and RT spaces on wedge or prism
    elements and bug fixes.

## [0.6.0] - 2022-02-17

  - Changed implementation of boundary conditions and excitations to do a better job
    separating spatial discretization and algebraic operators. This is also in preparation
    for future additions of boundary conditions which are not quadratic in the frequency
    such as wave ports, as well as linear algebra improvements.
  - Added a `config["Boundaries"]["Current"]` configuration file parameter to specify a
    surface current excitation, and removed the lumped port `"Termination"` option as it is
    no longer relevant. Magnetostatics simulations use this `"Current"` excitation rather
    than `"Port"` now.
  - Changed lumped port and surface current excitations to now use a unit circuit voltage
    and circuit current excitation, rather than unit electric field or surface current
    density.
  - Added numeric wave port boundaries for frequency domain driven simulations. To
    differentiate, lumped port boundaries are now specified using the keyword
    `"LumpedPort"` instead of just `"Port"` in the configuration file.
  - Fixed farfield boundaries to account for non-vacuum materials at the boundary and added
    second-order option for absorbing boundary condition (ABC).
  - Added finite conductivity surface impedance boundary condition, with optional thickness
    correction.
  - Added dielectric loss calculation specializations for metal-air, metal-substrate, and
    substrate-air interfaces in accordance with
    [this paper](https://aip.scitation.org/doi/10.1063/1.3637047).
  - Changed implementation of frequency domain linear algebra to eliminate the use of PETSc
    matrix types, relying instead on wrappers of the real and imaginary Hypre matrices
    making up a complex operator.
  - Added geometric multigrid AMS solvers (in two variants) using h- or p-coarsening, which
    employs either a specialized smoother for Nedelec space coarsening or geometric
    coarsening for the auxiliary space solves.
  - Changed build process to include PETSc as a submodule and configure as part the standard
    CMake build process.

## [0.5.0] - 2021-11-08

  - Changed postprocessing implementation which leads to some configuration file changes,
    namely `config["Domains"]["Postprocessing"]` and
    `config["Boundaries"]["Postprocessing"]`.
  - Changed some filenames for postprocessed quantities to be more intuitive.
  - Added dielectric (bulk and interface) loss calculations for all simulation types.
    Correctly handles two-sided internal surfaces using `"Side"` specification.
  - Added mutual capacitance matrix extraction to electrostatic solver.
  - Added built-in mesh support for COMSOL (`.mphbin` and `.mphtxt`) and Nastran (`.nas` and
    `.bdf`) mesh formats.
  - Removed Gmsh dependency for mesh format conversion.
  - Added some improvements to the default build flags and dependency build options.

## [0.4.0] - 2021-09-27

  - Added `"Capacitance"` and `"Inductance"` options to extract surface charges and fluxes
    in frequency and time domain simulation types.
  - Changed `"Transient"` solver implementation to use real-valued scalars (user must use a
    domain electrical conductivity as opposed to a loss tangent for loss modeling), for a
    roughly 2x performance improvement.
  - Added new MFEM features, including upgrades to high-order Nedelec element spaces on
    tetrahedra which no longer require global dof reordering.
  - Changed system matrix assembly to use Hypre P^T A P within MFEM and keep the full dense
    element matrices when the coefficient is nonzero (`skip_zeros = 0` always). The
    stiffness and mass matrices thus always share the same sparsity pattern now.

## [0.3.0] - 2021-09-02

  - Added electrostatic solver for capacitance matrix extraction.
  - Added `"Termination"` option for ports to allow for open or short circuiting but still
    enabling port postprocessing.
  - Changed timer to separate out time spent for disk I/O during initialization and
    postprocessing.
  - Added options to job wrappers for new instance types and subnet ID option.
  - Added changes from Hypre and MFEM dependencies to include improvements and bug fixes
    from those libraries.

## [0.2.0] - 2021-06-25

  - Added a proper [changelog](./CHANGELOG.md) for the project.
  - Added time domain solver with various excitation options to accompany existing frequency
    domain driven and eigenmode simulation types.
  - Added wrapper scripts in [`bin/`](./bin) for launching parallel simulations either
    interactively or as batch jobs using PBS.
  - Changed JSON configuration file keywords for frequency and time intervals.
  - Added domain conductivity models for normal and superconducting metals.
  - Fixed a bug for coaxial lumped port postprocessing.
  - Added visualization of surface currents on domain boundaries to accompany existing
    $\bm {E}$- and $\bm{B}$-field visualization.
  - Changed default length scale for nondimensionalization.
  - Changed default filenames for port-related postprocessed quantities.
  - Fixed many smaller issues and added various other improvements.

## [0.1.0] - 2020-11-03

  - Initial development release.
