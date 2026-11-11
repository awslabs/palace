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

#### New Features

  - Added `Model.LocalReorder`, which numbers the elements of each process along a Hilbert
    curve through their centroids once the mesh is partitioned (one curve through the whole
    mesh, so that the order does not depend on the partition). The edges, faces and degrees
    of freedom, numbered in the order of the elements, become local in space, which reduces
    the cost of each linear solver iteration, mostly in the smoothers of geometric
    multigrid. Unlike `Model.ReorderElements`, it has no serial cost and keeps the partition
    and, with nonconforming adaptive mesh refinement, the refined meshes: the root elements
    of the refinement tree are numbered along the curve with their vertex orders, which
    define how they are refined, so that the refined and rebalanced meshes are ordered as
    well, and the rebalancing numbers the elements in the order of the input mesh, which
    keeps its partitions. A probe on the boundary between elements takes the value in the
    element first in the order of the input mesh among those containing it. Results change
    by round-off: the faces take their orientation from the element first in the new order,
    which changes the basis of the face degrees of freedom but not the global phase of
    eigenvectors (see the eigenvector phase convention under Interface Changes). Conforming
    meshes refined after their distribution (conformal adaptive or parallel uniform
    refinement, where the marking of tetrahedra for refinement depends on the local edge
    order), conforming meshes with curved nodes, and nonconforming input meshes (such as
    saved adapted meshes) keep their element order. [PR TBD].
  - Added `config["Solver"]["Driven"]["AdaptiveIndicatorSamples"]`: with adaptive mesh
    refinement, an adaptive frequency sweep can compute the error indicator from its
    reduced-order model at this number of fixed frequencies of the sweep, uniformly spaced,
    instead of from its full-order samples. The samples follow the greedy search for the
    largest reduced-order model error, whose choices change with round-off, and so did the
    marked elements. With fixed frequencies, the marking no longer depends on the sample
    set, which also makes `"AdaptivePoleExclusion"` safe with refinement. The fixed
    frequencies under-weight the resonances narrower than their spacing, which the
    refinement then does not converge, and the indicator norm compared with
    `config["Model"]["Refinement"]["Tol"]` becomes their root mean square. The default `0`
    uses the full-order samples. [PR TBD].
  - Added `config["Solver"]["Driven"]["AdaptiveIndicatorMaxQ"]`: with
    `"AdaptiveIndicatorSamples"`, the error indicator for adaptive mesh refinement is also
    computed at the resonances of the reduced-order model in the band of the sweep with a
    quality factor below this value, which uniformly spaced frequencies only weigh when one
    falls within a few linewidths. The resonances are the poles whose eigenvectors satisfy
    the full-order equations to `1e-2` relative to the size of their terms, and each weighs
    as much as the fixed frequencies together. Keep it well below the quality factors of
    weakly coupled modes, which the reduced-order model resolves only when the sampling
    visits them. The default `0.0` adds none. [PR TBD].
  - Added `Solver.Linear.AMSSmoothIts`, the number of smoothing sweeps on the Nédélec space
    in each AMS iteration when AMS is the coarse solver of geometric multigrid. It defaults
    to 8 when `Solver.Linear.AMSSingularOperator` is enabled (by default, magnetostatic
    simulations without superconductor boundaries), which reduces the number of linear
    solver iterations, and to `Solver.Linear.MGSmoothIts` otherwise, as before.
    [PR 1159](https://github.com/awslabs/palace/pull/1159).
  - Added `config["Solver"]["Driven"]["CoarseReuseTol"]`: frequency sweeps can reuse the
    coarse-level sparse direct factorization of the multigrid preconditioner for frequencies
    within this relative distance of the frequency where it was computed, instead of
    factoring the coarse level for every frequency. Solves still converge to the linear
    solver tolerance; a solve which needs 50% more iterations than with new factorizations
    factors the coarse level again and continues. If the option is omitted, uniform sweeps
    reuse factorizations within `0.05` and adaptive sweeps do not (`0.0`): their samples are
    rarely that close, and a different coarse solve can change which samples the greedy
    selection picks. An explicit value applies to both sweep types, and `0.0` disables
    reuse. Preconditioners other than geometric multigrid are not affected. [PR TBD].
  - Added `config["Solver"]["Linear"]["MGAuxiliarySmoothOrder"]`, the polynomial order of
    the Chebyshev smoother in the auxiliary (gradient) space of the distributive relaxation
    smoother of geometric multigrid. If omitted, it is the order of the smoother of the
    primary space, as before, except with `Solver.Linear.ComplexCoarseSolve`, where it is
    `max(2 p, 4)` for solution order `p`: there, the lower order needs no more linear solver
    iterations and makes each iteration cheaper. [PR TBD].
  - The mesh information printed for the initial and every refined mesh includes the numbers
    of tetrahedra whose minimum dihedral angle is below 5, 2 and 1 degrees, and the smallest
    angle, if any: such slivers can slow down the convergence of the iterative linear
    solvers (see `Solver.Linear.MGSliverAngle`). [PR TBD].

#### Interface Changes

  - Reworked the performance report in `palace.json` (new `SchemaVersion` 2). Every
    `iterationX/palace.json` of an adaptive run now describes that iteration alone, as if it
    were a separate simulation (times, counts and peak memory restart at each iteration, and
    the work between iterations is charged to the iteration using the new mesh), and a
    `Cumulative` block holds the totals over all iterations so far. `ElapsedTime` keys are
    now the timer names (for example `KSP_COARSE_SOLVE` instead of `CoarseSolve`, which made
    `Setup`, `Solve` and `Construction` ambiguous); `Durations` holds the average over MPI
    processes instead of the root process's value, with `Min`, `Max` and `MaxRank` maps next
    to it. `LinearSolver/TotalSolves` and `TotalIts` now count all main linear solves of the
    iteration rather than those of the last solver object, and other solvers (for example
    the divergence-free projection and the error estimator) are reported under
    `LinearSolver/<Role>`. New rows `Startup`, `Coarse Setup` (coarse solver setup and
    factorization, previously part of `Coarse Solve`, or of `Setup` for MUMPS) and
    `Adapted Mesh Output` (previously part of `Rebalancing`) appear in the log and in
    `palace.json`, together with inclusive tallies (multigrid levels, smoothers, Krylov
    operator and orthogonalization, halo exchanges, reductions and waits), solver counters,
    resolved solver choices and memory attribution. The log tables keep their layout. An
    opt-in check of the coarse solve residual is enabled with
    `config["Problem"]["Verbose"] >= 2` or the `PALACE_PERF_VERIFY` environment variable.
  - When `Solver.Linear.Type` is `"Default"` (or omitted) for frequency domain problems,
    and for the wave port and boundary mode eigenvalue solves, *Palace* now prefers
    STRUMPACK over SuperLU_DIST, then MUMPS, as the sparse direct solver when built with
    several of them. Set `"Type": "SuperLU"` to keep using SuperLU_DIST. SchemaVer 4-2-0
    [PR 1148](https://github.com/awslabs/palace/pull/1148).
  - The global phase of the eigenvectors of eigenmode simulations is fixed by the projection
    of the electric field on a smooth reference field defined on the bounding box of the
    mesh (`∫ w ⋅ E dΩ` real and positive, see the reference documentation) instead of the
    mean of the degrees of freedom, which depended on the orientations of the edges and
    faces of the mesh, on the order of the elements and on the number of processes. The
    eigenvectors of lossless modes are real, with a definite sign. Relative to previous
    versions, this changes the phases of port voltages and currents, the phases of fluxes
    and field outputs, and the signs of the external coupling rates (`κ_ext`) and of signed
    energy participation ratios of eigenmodes: all gauge choices (Q, `|κ_ext|`, energies and
    participation magnitudes are unchanged), which now no longer depend on the numbering of
    the mesh or on the number of processes, except within degenerate eigenspaces and when
    the mean of the degrees of freedom is used as a fallback (with a warning, when the
    projections on both reference fields are negligible). [PR TBD].
  - `config["Model"]["LocalReorder"]` defaults to true: the elements of each process are
    numbered along a Hilbert curve once the mesh is partitioned (see `Model.LocalReorder`
    under New Features), which reduces the cost of each linear solver iteration with the
    partition, the elements and their refinement unchanged. Results change by round-off,
    except for probes at points on the boundary between elements (a vertex, an edge or a
    face), which take the value in the element first in the order of the input mesh: this
    can be another element than before, so that the probed values of fields which are
    discontinuous there (such as the normal electric field at a material interface) can
    change by O(1). Set `"LocalReorder": false` for the previous numbering and probes. [PR
    TBD].
  - Geometric multigrid solves exactly on the patches of degenerate tetrahedra of the
    complex-valued operators of driven and eigenmode simulations (new
    `Solver.Linear.MGSliverAngle`, default `-1`: tetrahedra with a minimum dihedral angle
    below 2 degrees, and none for the real-valued operators of electrostatic, magnetostatic
    and transient simulations; `0` disables it). Slivers have basis functions close to
    linearly dependent: the Chebyshev smoothers, with their diagonal scaling, cannot correct
    their local modes, and these can dominate the number of linear solver iterations. On
    every level above the coarsest, the true dofs of the mesh entities incident to the
    vertices of each sliver form a patch, also across processes (new
    `Solver.Linear.MGSliverLayers`: `1`, or `0` for its own entities; by default `-1`, the
    former unless their dense solves are expensive compared with the smoothing of the level,
    else none), and an exact solve with the dense block of the level operator on each patch
    precedes the coarse-grid correction, with the residual updated. A level whose patches
    would hold more than 1% of its true dofs keeps no patch solves, with a warning. The
    patch solves are not available with OpenMP or device backends, and a warning reports
    meshes with slivers whose patch solves such a backend disables. Meshes without slivers
    are not affected; on meshes with slivers, eigenmode simulations can need far fewer
    linear solver iterations. [PR TBD].
  - The `"Default"` option of `config["Solver"]["Linear"]["ColumnOrdering"]`, also its
    default, now uses the parallel PT-Scotch ordering (below) for `"Type": "STRUMPACK"`
    built with PT-Scotch for matrices of at least 50000 rows on more than one process, and
    otherwise (fewer than 50000 rows or one process), sequential METIS, as before. Results
    with STRUMPACK change at the solver tolerance: use `"ColumnOrdering": "METIS"` to keep
    the previous ordering. `"Type": "STRUMPACK-MP"` is not affected: as before, it ignores
    `"ColumnOrdering"` and the STRUMPACK compression options (it always orders with METIS
    and does not compress). [PR TBD].
  - `config["Solver"]["Linear"]["EstimatorTol"]`, the relative tolerance of the flux
    projection solves of the error estimator, defaults to the smaller of `1e-4` and 0.01
    times the refinement tolerance `config["Model"]["Refinement"]["Tol"]`, but not less than
    `1e-10`, instead of `1e-6`, for all simulation types: `1e-4` for the default refinement
    tolerance of `1e-2` (an explicit value always applies). The flux projection solves need
    far fewer iterations, while adaptive mesh refinement marks the same elements, up to
    elements at the marking threshold, as long as the relative error indicator is well above
    the tolerance. Very accurate solutions, whose indicator norm is within about an order of
    magnitude of `"EstimatorTol"`, need a tighter tolerance (for example `1e-8` or `1e-10`):
    the error of the projections adds about `"EstimatorTol"` to the indicator. Set
    `"EstimatorTol": 1.0e-6` for the previous behavior. [PR TBD].
  - The `"ParMETIS"` option of `config["Solver"]["Linear"]["ColumnOrdering"]` is now a
    deprecated alias for `"METIS"`, since *Palace* no longer depends on ParMETIS.
    SchemaVer 5-0-0 [PR 1130](https://github.com/awslabs/palace/pull/1130).
  - The flux-recovery error estimator now only runs by default when adaptive mesh refinement
    uses it (`Model.Refinement.MaxIts` > 0, except for transient simulations). Simulations
    without AMR no longer write `error-indicators.csv` or the `Indicator` field, and skip
    the estimator setup and the per-solve flux projections, which took a large fraction of
    the run time and memory of some simulations. Set the new
    `Solver.Linear.AlwaysComputeErrorEstimation` to `true` to compute the error indicator
    without AMR. [PR 1128](https://github.com/awslabs/palace/pull/1128).
  - The error estimator also no longer runs after the last refinement allowed by
    `Model.Refinement.MaxIts`, whose error indicator marked nothing: the final iteration of
    such a simulation does not write `error-indicators.csv` or the `Indicator` field, and
    the AMR summary reports its indicator norm as not computed, followed by the last
    computed one (of the mesh before the last refinement). Whether the final mesh satisfies
    `Model.Refinement.Tol` is then not checked; `Solver.Linear.AlwaysComputeErrorEstimation`
    computes the indicator of the final mesh and writes its `error-indicators.csv` as
    before. [PR TBD].
  - Adaptive frequency sweeps now start the HDM solve at each sampled frequency from the
    PROM solution there (with `"InitialGuess": true`, they used to start from the solution
    at the previous sample), and `config["Solver"]["Linear"]["InitialGuess"]` now defaults
    to `true` for all driven simulations. This reduces the linear solver iterations of the
    offline phase. The results of adaptive sweeps change, and so can the number and the
    frequencies of the samples: with an `"AdaptiveTol"` tighter than the accuracy of the
    linear solves, a sweep can now converge instead of sampling up to
    `"AdaptiveMaxSamples"`. Set `"InitialGuess": false` to restore the previous behavior.
    [PR 1126](https://github.com/awslabs/palace/pull/1126).
  - Driven simulations apply the system matrix through the finest level of the multigrid
    preconditioner and no longer assemble it separately for uniform sweeps, reducing CPU
    time by 1–5% and peak memory by up to 6%.
    [PR 980](https://github.com/awslabs/palace/pull/980).
  - GMRES and FGMRES now stop when the residual norm decreased by less than 0.1% over the
    last 100 iterations, instead of iterating to `config["Solver"]["Linear"]["MaxIts"]`: the
    solve has then reached the accuracy attainable for the preconditioned system (for
    example near high-Q resonances of driven simulations, or with a tolerance near the
    round-off floor), and further iterations would not change the solution. The new
    `config["Solver"]["Linear"]["StagnationIts"]` (`0` disables the test) and
    `"StagnationTol"` control the test. [PR TBD].
  - Without `config["Solver"]["Linear"]["MGSmoothOrder"]`, eigenmode simulations of solution
    order 2 use Chebyshev smoothers of order 8 instead of 4 on the multigrid levels whose
    shifted operator is nearly positive semi-definite, as estimated at setup by a few
    Lanczos iterations on the Jacobi-scaled real part. Higher orders reduce the linear
    solver iterations of nearly definite problems but not of strongly indefinite ones, which
    keep order 4. The multigrid preconditioner keeps the raised orders unless a short GMRES
    probe at setup (with the whole multigrid cycle from a random vector, alternating the two
    orders and stopped once the prediction is clear, after at most 20 iterations of each, or
    40 when the prediction is then close to the threshold) predicts that they need more than
    25% more iterations than order 4: with Floquet periodic and absorbing boundaries, the
    real-valued coarse solve approximates the complex operator poorly, and order 8 would
    increase the linear solver iterations instead. [PR TBD].
  - Without `config["Solver"]["Linear"]["ComplexCoarseSolve"]`, the STRUMPACK coarse-level
    solver of the multigrid preconditioner of eigenmode simulations can factorize the
    complex-valued coarse matrix instead of its real-valued approximation (the sum of the
    real and imaginary parts). At setup, if 10 iterations of GMRES on the coarse level
    preconditioned by the approximation do not reduce the residual below `1e-4`, for example
    with Floquet periodic or absorbing boundaries, the complex-valued matrix is factorized
    too, and kept if short GMRES probes with the whole multigrid cycle, as for the Chebyshev
    orders, predict that it needs at least 1.5 times fewer iterations than the approximation
    at the raised and at the base Chebyshev orders: this can reduce the linear solver
    iterations several times. The complex-valued factorization keeps about three times the
    memory of the real-valued one and peaks higher, more so with the later factorizations of
    nonlinear eigenvalue solves, and its setup takes several times longer. It is not tried
    if it is estimated not to fit, with a reserve for the later allocations, in the memory
    available to the process (the available memory of the node or the headroom of its memory
    control group), or if that is unknown; limits enforced without control groups and shared
    nodes are not seen. Its setup prints its estimate first and stops the simulation on a
    failure, with a hint to set `"ComplexCoarseSolve": false`. The choice is reused for the
    later shifts and adaptive mesh refinement iterations while the coarse level changes by
    at most a factor of 1.5. With the complex-valued coarse solve, the auxiliary space
    smoothers use the order `max(2 p, 4)`, as with an explicit `"ComplexCoarseSolve": true`.
    This applies on Linux to STRUMPACK on the host with fourth-kind Chebyshev smoothers;
    SuperLU and MUMPS coarse solvers, other simulation types and wave port mode solvers keep
    the real-valued approximation unless `"ComplexCoarseSolve": true` is given. SchemaVer
    6-0-0 [PR TBD].
  - `config["Solver"]["Linear"]["GSOrthogonalization"]` now defaults to `"CGS2"` instead of
    `"MGS"`, and with `"CGS2"` GMRES and FGMRES reorthogonalize each Krylov basis vector one
    iteration late (classical Gram-Schmidt with delayed reorthogonalization), in the same
    two passes over the basis as the first projection of the next vector and with one global
    reduction per iteration, instead of the `k + 1` reductions and `2 k` passes over basis
    vectors of MGS against `k` basis vectors. The iteration counts do not change and the
    basis stays orthogonal to round-off, while the orthogonalization of long unrestarted
    solves is much faster. The SLEPc eigenvalue solvers also use CGS2 by default. The CGS
    variants of GMRES, FGMRES and the reduced-order model basis make one pass over the basis
    for all inner products and one for all updates. [PR TBD].
  - With adaptive mesh refinement, the eigenvalue solve of each adaptation iteration can
    start from the eigenvectors of the previous one, transferred to the refined and
    rebalanced mesh (new `Solver.Eigenmode.WarmStart`, default `true`), which needs fewer
    linear solves. The warm start is skipped when the previous solve converged within its
    first Krylov cycle, which it cannot shorten, and the first Krylov cycle of a
    warm-started solve uses 1e-2 times the linear solver tolerance (not below 1e-12): the
    eigenvalues of the eigenpairs converging in that cycle are set by its first linear
    solves. [PR TBD].
  - The tolerance of the linear solves of the SLEPc Krylov-Schur solver can be relaxed as
    the requested eigenpairs converge (new `Solver.Eigenmode.LinearSolverRelaxation`,
    default `0`, off), bounded by the separation of each requested eigenvalue from the
    others and by its distance into the transformed spectrum. A solve is repeated without
    relaxation if it leaves fewer converged eigenpairs than requested, accepted a requested
    eigenpair whose Ritz value did not settle between restarts, or leaves one with a
    backward error above twice the tolerance times `max(1, |λ - σ| / |σ|)` for the shift σ
    (about the factor by which the residuals of the linear solves enter it); a relaxed
    linear solve which converges falsely is repeated with the strict tolerance. With it,
    `eig.csv` and the other outputs list only the `N` requested eigenpairs, instead of all
    converged ones. [PR TBD].
  - For problems with damping, the eigenvalue search excludes strongly damped modes with
    quality factor below `Solver.Eigenmode.MinQ` (new, default `0.54`) and growing modes
    with `|Q|` < 10, which could stall the convergence of the others. With
    `Problem.Verbose` >= 1, the log reports how many Ritz values of the last iteration lie
    outside of the search region, and how many of them converged. The default subspace
    dimension `Solver.Eigenmode.MaxSize` of the SLEPc solver is `max(3 N, N + 28)`, instead
    of `max(2 N, N + 15)`, for `N` ≥ 6 when `MinQ` is enabled. [PR TBD].
  - The new `Solver.Eigenmode.IntermediateTol` and
    `Solver.Eigenmode.IntermediateLinearSolverTol` (default: the final tolerances) set
    looser eigenvalue and linear solver tolerances for the adaptation iterations before the
    last one, whose eigenvectors only drive the mesh refinement. The computed eigenvalues
    change within the solver tolerances. [PR TBD].
  - With adaptive mesh refinement, the shift of the shift-and-invert spectral transformation
    of each adaptation iteration after the first moves from the target into the band of the
    requested eigenvalues of the previous one (new `Solver.Eigenmode.AdaptiveShift`, default
    `true`), 0.45 of the way to the farthest of them, clear of them and of clusters of close
    ones: the modes below the target, which shift-and-invert at the target makes the most
    dominant ones, are then less dominant than every requested mode, which are found with
    fewer restarts, and the errors of the linear solves enter the eigenvalues far from the
    target much less amplified. The search, and the tolerance as a relative error of the
    eigenvalues, do not change. Adaptation iterations after the first need fewer linear
    solves (the shift stays at the target after a solve which converged within two restarts
    there). A solve whose shift is above the middle of the target and the farthest requested
    eigenvalue (modes entering the requested set lowered it) is repeated once with the shift
    from the eigenvalues found, else (or with fewer converged eigenpairs than requested)
    with the shift at the target. More converged eigenpairs beyond `N` can be reported. [PR
    TBD].
  - The linear solves after the first one of each eigenvalue solve can be deflated with the
    harmonic Ritz vectors of the preconditioned operator closest to zero, from the Krylov
    subspace of the first solve (new `Solver.Eigenmode.LinearSolverDeflation`, the number of
    vectors, default `0`, off): the later solves use a two-level (A-DEF1) deflated
    preconditioner, which removes the eigenvalues near zero of the multigrid-preconditioned
    shifted operator that slow down GMRES. With `-1`, each eigenvalue solve deflates the
    harmonic Ritz values of magnitude below 0.25 (at most 64, possibly none) if GMRES for
    the projected operator predicts a gain at the measured cost of the deflation, solves
    without deflation verify the gain (or stop the deflation), and while it pays the vectors
    of the values below 0.5 and then 1, from the same first solve, are tried when the model
    predicts more gain and kept if they pay more. The deflation basis is stored in single
    precision on CPUs. Deflation is opt-in: it can shift the eigenvalues of modes far from
    the shift by more than the tolerance, which their backward errors do not show. [PR TBD].

#### Performance Improvements

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
    driven and eigenmode simulations. [PR 1146](https://github.com/awslabs/palace/pull/1146).
  - The Floquet port DtN correction and the wave-port modal correction now use a single MPI
    reduction per operator application instead of four per rank-one term, and only touch
    the port degrees of freedom instead of full-length vectors.
  - The Krylov-Schur eigenvalue solver of eigenmode simulations tests convergence after
    every step of the Arnoldi factorization of a cycle, instead of only after the full
    subspace of `Solver.Eigenmode.MaxSize` vectors, and ends each cycle after the first as
    soon as the requested eigenpairs converge and their Ritz values have settled to the
    tolerance since the previous check, for subspaces of up to 128 vectors (SLEPc patch
    `extern/patch/slepc/patch_krylovschur_check_interval.diff`,
    `EPSKrylovSchurSetCheckInterval`; builds against an unpatched SLEPc keep the
    end-of-cycle test). The last cycle of an eigenvalue solve of several cycles, as those of
    the adaptation iterations, extended its subspace to `MaxSize` vectors, one linear solve
    each, however early its eigenpairs converged; they now meet the same tolerance with
    fewer linear solves. A check costs a dense eigensolve of the small projected problem,
    which outweighs the linear solves it saves with larger subspaces: they keep the
    end-of-cycle test. A check does not end the first cycle, nor a cycle while two of the
    converged eigenvalues are within a relative distance of 1e-6 of each other: the other
    copies of a multiple eigenvalue (as with symmetries) enter the Krylov subspace only
    through rounding errors, after more steps, so such cycles keep the end-of-cycle test.
    The wave port mode solvers keep the end-of-cycle test. The problem setup guide describes
    degenerate eigenvalues and the PETSc option `-eps_krylovschur_check_interval 0`, which
    keeps the end-of-cycle test. [PR TBD].
  - The HDM solves of adaptive frequency sweeps start from zero instead of the PROM solution
    when the true residual of the PROM solution is at least half the norm of the right-hand
    side, the bound of the true residual check of GMRES: on multiscale models with many
    lumped ports, such initial guesses have residuals dominated by stiff rows which the
    preconditioned residual does not see, so that the solve from the guess fails the check
    and is solved again from zero. The check costs two inner products; the solves from the
    other PROM solutions are unchanged. [PR TBD].
  - The serial preprocessing of large meshes on the processes holding the serial mesh is
    faster, for the same preprocessed mesh: Gmsh 2.2 meshes (ASCII or binary, first-order
    elements, no periodicity) are read by a memory-mapped parser instead of MFEM's reader,
    the first-order nodes of a linear mesh are added once its topology is final, a mesh read
    without nodes is finalized once, with the refinement marking, the interior boundary
    cracking searches the element connectivity through a per-element array of face
    neighbors, and the root process extracts and writes the parts of the mesh distribution
    with one thread per process of its node, while the other processes wait without
    spinning. Added an MFEM patch with hash-indexed `DSTable` and `STable3D` and the
    element-to-element table built by counting, for the edge, face and element tables of
    meshes with tens of millions of elements. The environment variable `PALACE_GMSH_FAST=0`
    reads every Gmsh mesh with MFEM, and `PALACE_SETUP_THREADS` sets the number of threads
    of the root process for the mesh distribution. [PR TBD].
  - The serial mesh preprocessing of the processes holding the serial mesh (the root, or one
    process per node) uses one thread per process of the node for the connectivity search
    and the edge marking of the interior boundary cracking, and for the edge lengths and the
    marking of the tetrahedra for refinement, while the other processes wait without
    spinning. The preprocessed mesh is unchanged. [PR TBD].
  - Nonconforming adaptive mesh refinement no longer builds the whole serial mesh on every
    process: for a conforming input mesh, one process per node sends the root elements of
    the nonconforming mesh, from which every process builds the parallel mesh, and the
    curved nodes of its own elements; for a nonconforming input mesh (a restart from an
    adapted mesh), the refinement tree is replayed in parallel and the mesh is rebalanced to
    the serial partitioning. Added an MFEM patch which gives ghost indices only to the edges
    and faces of the ghost elements of nonconforming meshes. This lowers the memory of the
    mesh distribution and of restarts from large adapted meshes several times, for the same
    parallel mesh. [PR TBD].
  - With glibc, the memory freed after the mesh distribution and before each adaptive solve
    is returned to the system (`malloc_trim`), which lowers the resident memory of the
    processes. [PR TBD].
  - Nonconforming meshes are now rebalanced with a PT-Scotch partition of the element
    adjacency graph, instead of MFEM's split of the refinement tree leaves into equal
    chunks, which fragmented the subdomains over the AMR iterations. Processes have fewer
    neighbors, more of them on the same node, and exchange less data after a rebalance.
    [PR 1130](https://github.com/awslabs/palace/pull/1130).
  - Adapted meshes saved with `config["Model"]["Refinement"]["SaveAdaptMesh"]` are
    compressed with the fastest zlib level, for slightly larger files, and written in a
    background thread when a CPU of the root process's node is idle. [PR TBD].
  - The root process writes an adapted nonconforming simplex mesh with linear nodes from the
    gathered refinement tree, instead of building the serial parallel mesh and its nodal
    space, for the same file. With nonconforming adaptive mesh refinement, the processes
    which read the serial mesh, one per node, read it at the same time instead of after the
    root. [PR TBD].
  - Faster writing of the adapted mesh saved with
    `config["Model"]["Refinement"]["SaveAdaptMesh"]` when no CPU of the root process's node
    is idle, as with one MPI process per core: threads on the CPUs of the other processes of
    the node, which wait without spinning meanwhile, format and compress the mesh file in
    independent chunks that form one gzip stream, and the gathering of the nonconforming
    mesh on the root process no longer prunes it. The decompressed file is unchanged. [PR
    TBD].
  - Added a STRUMPACK patch which lowers the memory of each process of the distributed
    sparse direct solver: the fill-reducing permutation is kept once per node in an MPI-3
    shared memory window, with the sequential and the parallel PT-Scotch orderings, and no
    arrays of the global matrix dimension are kept on every process. This lowers the memory
    per process of large coarse solves at high process counts; the shared memory windows
    need a `/dev/shm` large enough for the permutation, and the environment variable
    `STRUMPACK_SHARED_PERM=0` keeps it on every process instead (`STRUMPACK_MEMORY_VERIFY`
    checks the map of the rows to the fronts). [PR TBD].
  - Added an MFEM patch with which the parallel conforming interpolation of a finite element
    space on a nonconforming mesh stores rows only for the ghost DOFs it references, instead
    of for every ghost DOF of the coarse mesh replicated on every process. This removes a
    large transient of the memory of every finite element space build or update on adapted
    meshes. [PR TBD].
  - On CPUs, every volume and boundary operator assembles its quadrature data, and data
    which is constant on each element up to the quadrature weight (affine elements with
    elementwise constant coefficients) is stored once per element instead of once per
    quadrature point, which lowers the memory of the partially assembled operators. [PR
    TBD].
  - Eigenmode simulations release the eigenvalue solver, the linear solver, the
    preconditioner, the divergence-free projection and the system matrices before
    postprocessing and error estimation, and construct the error estimator after the
    eigenvalue solve, so that postprocessing and error estimation no longer add to the peak
    memory of the solve. [PR TBD].
  - Complex operators also apply the libCEED operators of their real or imaginary part
    without a partner in the other part, and the H1 diffusion operators, in one packed pass
    for both parts, so that the fine and auxiliary space operators of the multigrid
    preconditioners of driven and eigenmode simulations use the packed fast path. [PR TBD].
  - On CPUs, small libCEED operators on simplices (single-component H1 and boundary
    operators, and real square operators) are applied through dense element matrices instead
    of partial assembly when this moves fewer bytes. The element matrices are assembled on a
    copy of the operator, and again after in-place updates of its inputs. [PR TBD].
  - On affine simplices with elementwise constant coefficients, H(curl) and H(div) operators
    use an exact quadrature rule with fewer points, the dimension of the polynomial space,
    with the basis represented in an orthonormal polynomial basis, and the geometry factors
    of affine elements are computed once per element. The exactness is verified at setup;
    otherwise the default rule is used. [PR TBD].
  - Products with HypreParMatrix prolongations on the host apply their identity rows as
    copies, their transpose through a gathered copy, and the real and imaginary parts of
    complex vectors in one pass with one halo exchange, with the same results as hypre, and
    complex operators apply the prolongation of conforming spaces to their real and
    imaginary parts with one halo exchange. [PR TBD].
  - The Chebyshev smoothers of the geometric multigrid preconditioner update the iterate in
    the same pass over memory as the next search direction. Results are unchanged. [PR TBD].
  - The element-assembled terms of CPU operators are applied by a kernel which processes
    four rows of the element matrices at a time, and their element matrices are assembled
    over chunks of elements, which bounds the temporary memory of the assembly. [PR TBD].
  - On CPUs, the Chebyshev smoothers of the geometric multigrid preconditioner of frequency
    domain and transient simulations apply the operator of each multigrid level to their
    search directions through single-precision copies of its element matrices, with the
    orientations of the Nédélec restrictions folded in and the real and imaginary parts
    applied in one pass, and the auxiliary space operators of the distributive relaxation
    smoother store their element matrices in single precision. The arithmetic stays in
    double precision, and the residuals of the smoothed iterates, the coarsest level and the
    system matrix are unchanged, so the preconditioner only changes at the level of
    single-precision rounding. Imaginary volume terms below that level relative to the real
    part, for example a small dielectric loss tangent, are left out of these applications.
    With unchanged linear solver iterations, this shortens the multigrid V-cycle where the
    smoothers are bound by their operator applications (orders 3 and 4). The element
    matrices take additional memory and setup time (about 5 kB and 0.2 ms per tetrahedron at
    order 3, 16 kB and 0.5 ms at order 4). Setting the new
    `config["Solver"]["Linear"]["MGSmoothSinglePrecision"]` (default `true`) to `false`
    disables it. [PR TBD].
  - With GMRES and FGMRES, the distributive relaxation smoother of the geometric multigrid
    preconditioner no longer computes the residual of its iterate between the Chebyshev
    smoother of the primary space and the auxiliary space correction: the auxiliary space
    correction of the pre-smoother uses the residual before the last Chebyshev step, and the
    first Chebyshev step of the post-smoother the residual before the auxiliary space
    correction. This saves two of the four double-precision operator applications per
    multigrid level and V-cycle, which makes the V-cycle faster (also without
    single-precision smoothing), while the linear solver iterations change little. As the
    preconditioner is no longer symmetric, CG and MINRES instead update these residuals
    through the single-precision operator applications of the smoothers, for the same
    smoother up to rounding. Since single-precision rounding of the curl-curl element
    matrices does not preserve their null space, each multigrid level measures the relative
    error of the single-precision operator on the range of the discrete gradient when its
    operators are set, and keeps exact residuals where it exceeds 0.1, as on meshes with
    elements much smaller than the wavelength. Driven and transient simulations also measure
    the error on the auxiliary space, scaled by the inverse diagonal of its operator, which
    large terms of the operator on gradients (such as the inductance of impedance
    boundaries, or lumped ports) do not hide: where it exceeds 100, as at low frequencies,
    the residual projected on the auxiliary space is computed exactly. The errors and the
    residuals used are reported per level in the `Chebyshev` entries of the linear solver in
    `palace.json`. The new `config["Solver"]["Linear"]["MGAuxiliaryResidual"]` (`"Lagged"`,
    `"Updated"`, or `"Exact"` for the previous behavior; `"Default"` picks by Krylov solver)
    selects them. [PR TBD].
  - The Chebyshev and Jacobi smoothers use the exact diagonal of the partially assembled
    operators of all multigrid levels. Previously, on adapted (nonconforming) meshes, the
    diagonal was approximated with the absolute values of the constraint weights of the
    hanging dofs, which overweights their masters and leaves out the couplings between the
    constrained dofs of each element, and for Nédélec spaces of order p > 1 on tetrahedra
    and prisms, libCEED approximated the diagonal of the face dofs with the absolute values
    of their orientation transforms. The exact diagonal is assembled from the libCEED
    assembled QFunction, forming only the element matrix entries it needs, with the
    communication of the transpose of the prolongation. This reduces the number of linear
    solver iterations of eigenmode simulations on adapted meshes, most at orders 3 and 4,
    and slightly on conforming meshes, for a diagonal assembly of a small fraction of the
    solve time. [PR TBD].
  - The discrete interpolators of the multigrid smoothers and the transfer operators between
    multigrid levels apply the restriction to the true dofs of their range (a selection of
    local dofs) as copies of contiguous runs instead of a sparse matrix product, with the
    same results. On adapted meshes, the restriction and its transpose were a large part of
    the local work of these operators. [PR TBD].
  - Square operators with eliminated essential boundary conditions apply the prolongation of
    their space to the input with zeros at the essential true dofs without forming a copy of
    it, and their `AddMult` adds the restriction of the result to the output directly
    instead of through a temporary vector, with the same results. On adapted meshes, the
    products with the prolongation and its transpose, with these copies, took a large part
    of the operator applications of the multigrid preconditioner. [PR TBD].
  - The prolongation of a nonconforming space stores its constraint coefficients as 16-bit
    indices into a table of their distinct values (at most a few thousand on adapted meshes)
    instead of one double per entry, which reduces the memory streamed by the hanging-dof
    rows of its products and of the products with its transpose, with the same results. [PR
    TBD].
  - Operators whose local action consists of element matrices, such as the single-precision
    approximations applied to the search directions of the Chebyshev smoothers and the
    element-assembled auxiliary-space operators of the distributive relaxation smoothers,
    gather the element coefficients directly from the true dofs and from the values of the
    other local dofs of the prolongation (hanging, ghost and essential dofs), and add their
    results the same way, instead of forming local vectors around the element products, with
    the same results. With the three previous changes, the multigrid preconditioner with
    Nédélec elements of order 3 and 4 is substantially faster on adapted meshes, with the
    same results. [PR TBD].
  - Added libCEED patches for faster partially assembled operators on CPUs: the blocked CPU
    backends, among them the default `/cpu/self/xsmm/blocked`, apply operators with at least
    512 elements to blocks of 32 elements instead of 8 (smaller operators keep blocks of 8),
    the LIBXSMM backend keeps the kernels of each basis contraction instead of looking them
    up at every call, the bases of simplices contract all components of their vector values
    at the quadrature points (gradients, curls, H(curl) and H(div) interpolation) in one
    matrix product, and the restrictions of H(curl) and H(div) spaces skip the terms of
    their orientation transforms that are zero on every element. On Graviton4 with one MPI
    process per core, partially assembled operators on tetrahedra apply 10% to 20% faster at
    orders 3 and 4 and 25% to 35% faster at order 2. Results change at round-off level: the
    larger element blocks change the order in which the contributions of the elements
    sharing a degree of freedom are summed. [PR TBD].
  - Frequency sweeps (uniform and adaptive) no longer assemble the exact diagonals of the
    multigrid levels of the preconditioner at every frequency: from the second frequency on,
    they are combined from the diagonals of the stiffness, damping and mass terms, assembled
    once, as the diagonal is linear in the coefficients of the terms. The diagonals agree to
    round-off. Problems with frequency-dependent boundary conditions or Floquet periodicity
    assemble them at every frequency, as before. [PR TBD].
  - With single-precision multigrid smoothing
    (`config["Solver"]["Linear"]["MGSmoothSinglePrecision"]`, the default), frequency sweeps
    likewise form the single-precision element matrices of the smoothers of the levels above
    the coarsest from those of the stiffness, damping and mass terms, assembled once with
    their diagonals, instead of assembling them at every frequency. The element matrices
    then carry the single-precision rounding of the terms instead of that of their sum,
    which changes the iterates of the linear solver at the level of that rounding. The
    element matrices of the terms take about twice the memory of those of the smoothers. [PR
    TBD].
  - Faster STRUMPACK solve phase, which is applied once per Krylov iteration when STRUMPACK
    is the coarse or the only solver: the elimination tree is mapped to processes by factor
    memory instead of factorization flops, which balances the solve; added STRUMPACK patches
    that store the inverses of the separator blocks of the distributed fronts so their
    triangular solves become matrix products (selective inversion, only for blocks with a
    1-norm condition number below 9e10: the nearly singular blocks of driven problems at low
    frequency keep the backward stable triangular solves), stream the local factors once for
    a few right-hand sides, and use a PBLAS block size of 128 in the superbuild. Results
    change at round-off level. The environment variables `STRUMPACK_SI_MAX_COND` and
    `STRUMPACK_SI_VERBOSE` set the condition number limit of the inversion and report its
    decisions. [PR TBD].
  - Added `config["Solver"]["Linear"]["STRUMPACKSelectiveInversion"]` for the selective
    inversion above, which adds to the time of each factorization and saves time in each
    solve: by default it is used for eigenmode, boundary mode, and transient simulations,
    whose factorizations serve many solves, and not for driven, electrostatic, and
    magnetostatic simulations, where a factorization serves the Krylov solve of one
    frequency or of the sources (driven sweeps were slightly faster without it, also when
    the coarse factorization is reused across frequencies). [PR TBD].
  - Added a STRUMPACK patch that balances its solve phase further: the factor memory mapping
    of the elimination tree also counts a fixed cost per front, so that a subtree of very
    many small fronts (for example a nearly disconnected region of the mesh) is not left to
    one or a few processes, and where a few processes share a subtree they are split between
    the two child subtrees so that the largest work per process is minimal, instead of by
    rounding the proportional share. The solve also no longer traverses the local subtrees
    to size its workspace in every solve. Results change at round-off level. The environment
    variables `STRUMPACK_PM_FRONT_COST`, `STRUMPACK_PM_DP` and `STRUMPACK_PM_MINMAX` tune
    the mapping. [PR TBD].
  - STRUMPACK (`"Type": "STRUMPACK"`) built with PT-Scotch now computes its default column
    ordering (the default `config["Solver"]["Linear"]["ColumnOrdering"]`) in parallel,
    instead of with METIS on a single process, whose time grows with the size of the matrix
    while the other processes wait and which gathers the whole graph on that process.
    PT-Scotch computes the distributed separators on at most 256 processes (the environment
    variable `STRUMPACK_PTSCOTCH_NPROCS` sets their number) and METIS orders the subgraph of
    each process. The fill-in matches that of METIS, and ordering a coarse problem with tens
    of millions of unknowns on thousands of processes takes seconds instead of minutes, once
    per factorization (each adaptive mesh refinement iteration). The parallel ordering is
    used from 50000 rows on (see Interface Changes): small matrices on few processes, such
    as those of the wave port mode problems, which are factorized at each exact port solve,
    are ordered faster by METIS. [PR TBD].
  - Added a hypre patch, a prototype for
    [hypre-space/hypre#1654](https://github.com/hypre-space/hypre/issues/1654), with which
    AMS gathers the coarse levels of its auxiliary space BoomerAMG hierarchies onto fewer
    processes during the solve and runs the additive Pi-space corrections concurrently.
    These levels have few rows per process and are latency bound at high process counts. The
    AMS preconditioner also flags a zero initial guess as zero. Results agree up to the
    summation order. [PR 1154](https://github.com/awslabs/palace/pull/1154).
  - Used two levels of aggressive coarsening in the BoomerAMG coarse solver of the multigrid
    preconditioner of the divergence-free projection, instead of one. Its V-cycle is latency
    bound at high process counts, and the shorter hierarchy makes it faster.
    [PR 1149](https://github.com/awslabs/palace/pull/1149).
  - Added SuperLU_DIST patches for
    [xiaoyeli/superlu_dist#250](https://github.com/xiaoyeli/superlu_dist/pull/250),
    [xiaoyeli/superlu_dist#251](https://github.com/xiaoyeli/superlu_dist/pull/251),
    [xiaoyeli/superlu_dist#252](https://github.com/xiaoyeli/superlu_dist/pull/252) and
    [xiaoyeli/superlu_dist#253](https://github.com/xiaoyeli/superlu_dist/pull/253), which
    lower the communication cost of SuperLU_DIST as the coarse solver at high process
    counts: point-to-point exchanges of only the nonempty messages instead of dense
    all-to-all collectives in the symbolic factorization and in the final permutation of the
    solution, a static schedule reused by refactorizations with the same sparsity pattern,
    cached row indices for the redistributions of the right-hand side, which on several
    nodes go through one process per node, the matrix values for MC64 gathered on one
    process only, and fewer barriers in the triangular solves. Results are unchanged.
    [PR 1162](https://github.com/awslabs/palace/pull/1162).
  - Solved the divergence-free projection of the result of each shift-and-invert linear
    solve of the eigenmode solver to a relative tolerance of 10⁻³ (or
    `config["Solver"]["Linear"]["DivFreeTol"]`, if larger), as the second projection of the
    linearized quadratic problem, unless the target is far below the eigenfrequencies (see
    the bug fix for spurious eigenpairs at such targets). For a projected basis the
    right-hand side of the solve is divergence-free, so the divergence of its result only
    comes from the error of the linear solve, which the projection previously reduced to
    machine precision: on fine geometric features this took several times the iterations of
    the second projection. `"DivFreeTol"` remains the tolerance for the starting vectors.
    [PR TBD].
  - Added an entity-block Jacobi preconditioner for the flux recovery of the error
    estimator: the true degrees of freedom of each edge, face, and element interior form a
    dense block of the mass matrix, which corrects the near linear dependence of the
    high-order basis functions of an entity on anisotropic elements. It replaces Jacobi in
    CPU builds without OpenMP, with `config["Solver"]["Linear"]["EstimatorMG"]` false and a
    partially assembled mass matrix, for H(div) and H(curl) recovery spaces of order
    (`config["Solver"]["Order"]`) two and three on meshes made of tetrahedra, prisms, or
    both, or of triangles, and of order four on meshes made only of tetrahedra or only of
    triangles. The conjugate gradient iterations of the recovery drop by a factor of 1.4
    to 2.9 on isotropic meshes and of 2 to 5.6 on strongly anisotropic tetrahedra, also on
    nonconforming meshes. They stop on the residual norm of Jacobi preconditioning, so that
    the recovered flux at a given `EstimatorTol` is about as accurate as with Jacobi: its
    error in the mass norm is within about a factor of two of that with Jacobi. The other
    cases (GPUs, OpenMP, hexahedra, pyramids, quadrilaterals, other orders, H1 and L2
    recovery spaces) keep Jacobi. [PR TBD].
  - Applied the entity-block Jacobi preconditioner of the error estimator's flux recovery
    through the upper triangles of its symmetric inverse blocks, which halves the memory
    traffic of its applications: on fourth-order tetrahedra, the preconditioner of a complex
    recovery is applied 1.3 to 1.5 times faster, with the same iterations. The
    preconditioner is the same, up to rounding. [PR TBD].
  - Started each flux recovery of the error estimator after the first one on a mesh (for the
    next frequency sample, terminal or mode) from the closest combination of the previous
    recovered fluxes in the mass norm, so that the initial error is never larger than from
    zero: the last 20 fluxes in driven simulations, whose span approximates the new flux as
    the frequency samples of an adaptive sweep span the solution, and the last one
    otherwise. This reduces the conjugate gradient iterations of the recovery in frequency
    sweeps, with the same stopping criterion. [PR TBD].
  - Adaptive frequency sweeps can skip, in the search for the largest error of the
    reduced-order model, its resonances narrower than the new
    `config["Solver"]["Driven"]["AdaptivePoleExclusion"]` (relative half-width, default `0`,
    off; for example `1e-4`), once their poles have converged to `"AdaptiveTol"` times this
    width. Within a few linewidths of a weakly coupled mode with a quality factor of `1e6`
    to `1e11`, the error stays of order one whatever the number of samples, while the error
    indicator is largest at every resolved pole: sweeps kept sampling such modes, and each
    of these samples restarted the count of `"AdaptiveConvergenceMemory"`, which made the
    number of samples erratic. With the exclusion, sweeps over bands with such modes need
    fewer samples, for the same S-parameters. The excluded resonances are listed in the log,
    and frequencies of the sweep within them are still searched. With adaptive mesh
    refinement, the excluded samples are also missing from the error indicator, which then
    changes the refinement path and can end on a coarser mesh, so the exclusion is opt-in:
    it suits adaptive sweeps without refinement. [PR TBD].

#### Bug Fixes

  - Fixed STRUMPACK sparse direct solves which could crash, fail with `ZERO_PIVOT`, or hang
    when a process holds very few rows of the matrix, for example a coarse multigrid level
    at a high MPI rank count. The sparsity pattern of the matrix passed to STRUMPACK is no
    longer declared symmetric, since hypre's parallel matrix products can store explicit
    zeros without their transposes. [PR 1137](https://github.com/awslabs/palace/pull/1137).
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
  - Adaptive frequency sweeps no longer abort with
    `Linearly dependent vector added to PROM basis` when the real or imaginary part of an
    HDM solution started from the PROM solution is already in the span of the reduced basis
    to round-off, for example the small imaginary part of a nearly lossless problem: only
    the independent part is added to the basis, with a warning. A sample which adds no
    direction still aborts. [PR TBD].
  - Adaptive frequency sweeps no longer stop with a false convergence of GMRES at
    frequencies close to a resonance when the solution is accurate: the left-preconditioned
    residual can weigh the components of the residual by many orders of magnitude, and an
    initial guess such as the PROM solution of an HDM sample leaves residual components
    which the preconditioned residual does not see. A solve from a nonzero initial guess
    whose true residual check fails is now solved again from zero. If the check still fails,
    the solution is continued with right preconditioning, which minimizes the true residual:
    if the continuation satisfies the check within a relative distance of 10⁻⁴ of the
    solution, the solution is kept, else (an unstable preconditioner moves it by about its
    norm) the false convergence still stops the simulation. Solves which pass the check are
    unchanged. [PR TBD].
  - GMRES with left preconditioning (the default) and a multigrid preconditioner (geometric
    multigrid, AMS or BoomerAMG) now checks the true residual of a solve which converged on
    the preconditioned residual: the simulation stops with an error if the true residual is
    not reduced below 0.5 (or 10⁷ times the tolerance) of the norm of the right-hand side,
    instead of continuing with an invalid solution (false convergence, as with an unstable
    preconditioner). The check costs one operator application per linear solve. [PR TBD].
  - Fixed wrong results of simulations on nonconforming meshes (adaptive mesh refinement, or
    restarts from an adapted mesh) of order 4 and above on tetrahedra (3 and above on
    hexahedra, 2 and above for the Nédélec space), when a face shared by two processes has a
    rotated relative orientation, which only meshes with inverted or tangled elements have:
    the conforming interpolation of its DOFs was permuted on one of the processes. Whether a
    mesh partitioning cuts such a face is a matter of chance; in eigenmode simulations,
    GMRES then reported convergence in a few iterations and the eigenpairs were invalid
    (with large backward errors). An MFEM patch fixes it. [PR TBD].
  - The SLEPc Krylov-Schur eigenvalue solver (linear problems, and quadratic ones with
    `config["Solver"]["Eigenmode"]["PEPLinear"]`) no longer accepts an eigenpair whose Ritz
    value just appeared: the error estimate cannot see the error of the linear solves, which
    could leave eigenvalues several times `config["Solver"]["Eigenmode"]["Tol"]` off, for
    example after the Ritz values outside of the `"MinQ"` search region were purged, or in
    warm-started adaptation iterations which converged in their first Krylov cycle. An
    eigenpair now converges only once the nearest Ritz value of the previous restart is
    within a relative distance of `Tol`, or of up to 10⁻³ and a tenth of the distance to the
    second nearest one, except in the first Krylov cycle of a solve without warm start; one
    which does not settle within two more restarts converges on its error estimate. This
    takes up to one more restart per solve, and fewer eigenpairs beyond the requested `N`
    converge, so `eig.csv` can list fewer of them (the requested ones are unaffected). [PR
    TBD].
  - With `config["Solver"]["Linear"]["GSOrthogonalization"]` `"CGS2"` (the default), the
    SLEPc eigenvalue solvers flagged a breakdown of the Krylov basis whenever the second
    Gram-Schmidt pass reduced the norm of a new basis vector at all, instead of below 1/√2
    of its norm after the first pass. Near an invariant subspace, as after a warm start of
    an adaptive mesh refinement iteration, the solve then stopped with `DIVERGED_BREAKDOWN`
    and the simulation with fewer converged eigenpairs than requested. [PR TBD].
  - Fixed spurious eigenpairs, converged with small backward errors, and gradient components
    in the eigenvectors of eigenmode simulations with a target far below the
    eigenfrequencies (for example, a target of 0.5 GHz for modes at 15 GHz without damping,
    or of 0.015 GHz for modes at 12 to 16 GHz with damping). Gradients are eigenvectors of
    the shift-and-invert operator which it amplifies r = |λ - σ| / |σ| times more than an
    eigenvector with eigenvalue λ, for the shift σ, and the divergence-free projections of
    the linear solve results to a relative tolerance τ leave about τ of their gradient
    components, which grow in the Krylov basis as τ r approaches 1. These projections used τ
    = 10⁻³ (the second one of the quadratic problem since PR 1138, all of them since the
    change above). The eigenvalue solver now keeps τ r <= 0.1 for the requested eigenvalues,
    checked after each restart for the wanted Ritz values and before the solve for the
    eigenvalues of the previous adaptation iteration: otherwise the solve stops at the end
    of its Krylov cycle, without further linear solves, and is repeated with τ = 0.01 / r,
    also below `config["Solver"]["Linear"]["DivFreeTol"]`, at most twice, and a warning
    follows a solve whose last repeat still has τ r > 0.1. With damping, r is about the
    ratio of the eigenfrequencies to the target, without damping its square, so τ = 10⁻³ is
    kept for targets down to about 100 times below the eigenfrequencies with damping and 10
    times without. [PR TBD].
  - Added a STRUMPACK patch for
    [pghysels/STRUMPACK#164](https://github.com/pghysels/STRUMPACK/pull/164), which fixes
    the solves after the matrix values change when STRUMPACK compression is enabled
    (`Solver.Linear.STRUMPACKCompressionType` other than `"None"`) and the reordering is
    reused (`Solver.Linear.ReorderingReuse`, the default). With the STRUMPACK patch of
    [PR 1143](https://github.com/awslabs/palace/pull/1143), these solves could abort with an
    MPI error at high process counts, or give a poor preconditioner which needs more linear
    solver iterations. [PR 1164](https://github.com/awslabs/palace/pull/1164).

#### Build system

  - The superbuild now always builds the dependencies which *Palace* links directly.
    Previously, configuring with cuDSS as the only sparse direct solver failed.
    [PR 1133](https://github.com/awslabs/palace/pull/1133).
  - *Palace* no longer depends on ParMETIS, and depends on PT-Scotch instead. The superbuild
    builds PT-Scotch 7.0.18, and builds SuperLU_DIST, STRUMPACK, and MUMPS with METIS alone.
    SuperLU_DIST is patched to keep its METIS column ordering, also its default ordering,
    without ParMETIS. From version 0.20, the Spack package depends on `scotch` instead of
    `parmetis`, and leaves the use of ParMETIS by the sparse direct solvers to their own
    variants. *Palace* links ParMETIS only for sparse direct solvers built with it, given
    `PARMETIS_DIR`. The containers are built without ParMETIS.
    [PR 1130](https://github.com/awslabs/palace/pull/1130).
  - The superbuild patches the PT-Scotch ordering of STRUMPACK (above), and from
    version 0.20 the Spack package depends on `strumpack+scotch` with the same
    patch. *Palace* no longer requests `MPI_THREAD_MULTIPLE` for STRUMPACK with PT-Scotch,
    which does not need it, only for STRUMPACK with SLATE. Unless `SCOTCH_PTHREAD_NUMBER` or
    `SCOTCH_DETERMINISTIC` are set, *Palace* sets them so that the orderings run on one
    thread per process and are reproducible, also with a Scotch built with threads. [PR
    TBD].
  - The Spack package of `palace@develop` applies the dependency patches of the superbuild,
    in its order, and requires the versions of the dependencies which these patches are made
    against where they are not upstream, so that a later version cannot build without them:
    hypre 3.2.0 or 3.2.1, MFEM 4.10.0, libCEED 1.0.1, which a local libceed package adds,
    and SLEPc 3.26.0. [PR TBD].

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
