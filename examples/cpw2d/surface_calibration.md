# Thin-metal surface-participation calibration

This study compares an infinitely thin CPW model with a matched
fabrication-resolved model containing 100 nm metal, 50 nm substrate overetch,
80 degree sidewalls, and 10 nm corner radii.

For each physical edge segment `e`, fit the thin-metal bulk and surface annuli
at the smallest radii to

```text
p_bulk,ann,e(R) / R = A_e + B_e R
p_surf,ann,e(R)     = a_e + b_e R.
```

`A_e` is the local thin-sheet singular amplitude. The surface-annulus slope
gives the smooth contribution inside the matching radius,
`p_smooth,in,e(R) = max(0, b_e R)`. For interface `j`, a matched thin and
fabrication-resolved 2D pair defines the additive process coefficients

```text
A = sum_e A_e
C_in,j(R)  = (p_fab,in,j(R) - p_thin,smooth,in,j(R)) / A
C_out,j(R) = (p_fab,out,j(R) - p_thin,out,j(R)) / A.
```

The corrected target participation is

```text
p_corrected,j(R) =
    p_thin,out,j(R) + p_thin,smooth,in,j(R)
    + C_in,j(R) sum_e A_e + C_out,j(R) sum_e A_e.
```

The mesh-sensitive raw inside energy is used only to obtain the complementary
outside integral, `p_out = p_total - p_in`. It is never blended into the
corrected inside term. Refinement-driven growth of the unresolved raw
singular core therefore cannot leak back into the corrected result. `C_out`
is an additive fabrication correction, not a scale applied to all outside
energy.

The coefficients are calibrated on the 10/6 um CPW and can be applied
unchanged to other geometries made with the same process. Choose `R` where the
fits and coefficient tables vary slowly, subject to both:

```text
fabrication feature size << R
2 R < distance to the next geometric feature
```

The matching radius is a numerical regularization parameter and does not have
to be shared by SA, MS, and MA. Select and freeze one radius window per
interface using independent 2D calibration and validation geometries. The
calibration and target simulations should use the same FEM order and comparable
matching-tube resolution. Never tune the radii against the
fabrication-resolved result being used as a blind validation target.

Generate and run the study from this directory:

```text
cd mesh
julia generate_surface_calibration_meshes.jl
cd ..
../../build/bin/palace -np 4 cpw2d_surface_calibration_thin.json
../../build/bin/palace -np 4 cpw2d_surface_calibration_fabricated.json
../../build/bin/palace -np 4 cpw2d_surface_validation_thin.json
../../build/bin/palace -np 4 cpw2d_surface_validation_fabricated.json
python3 local_edge_correction.py calibrate \
  --thin postpro/surface_calibration_thin \
  --fabricated postpro/surface_calibration_fabricated \
  --output postpro/local-edge-process.csv
python3 local_edge_correction.py apply \
  --calibration postpro/local-edge-process.csv \
  --input postpro/surface_validation_thin \
  --reference postpro/surface_validation_fabricated \
  --output postpro/surface-validation-local-corrected.csv
```

The first command writes the process coefficients. The second writes each
radius result and the selected multi-radius summary, including the four
additive terms, fit diagnostics, AMR spreads, and reference error.

## Coupled response for nearby edges

The self-consistent response-matrix correction uses an isolated one-edge
coupon only while its matching patch is disjoint from every other patch. For
two edges separated by width `w` and matching radius `R`, use a coupled
two-edge coupon when `w < 2R`. At `w = 2R` the contours touch but do not
overlap, so either representation is geometrically admissible.

`mesh/mesh_edge_pair_coupon.jl` creates matched thin and fabricated coupons. By
default the two edges are equipotential, as for a ground-plane cutout. For
example, for a 2 um cutout and `R = 2 um`:

```text
mkdir -p /tmp/edge-pair
julia mesh/mesh_edge_pair_coupon.jl thin 2.0 /tmp/edge-pair/edge_pair_thin.msh
julia mesh/mesh_edge_pair_coupon.jl fabricated 2.0 \
  /tmp/edge-pair/edge_pair_fabricated.msh
python3 generate_edge_pair_response.py \
  --cutout-width 2.0 --radius 2.0 --metal-thickness 0.1 \
  --output /tmp/edge-pair
../../build/bin/palace -np 4 /tmp/edge-pair/edge_pair_thin.json
../../build/bin/palace -np 4 /tmp/edge-pair/edge_pair_fabricated.json
```

The generator phase-aligns the periodic contour basis and inserts both the
lower thin-sheet junction and upper fabricated-metal junction at every
conductor cut. `--metal-thickness` must equal `t_metal` in the fabricated mesh.
This gives the thin and fabricated coupons the same trace space. A basis which
interpolates across a conductor cut imposes an artificial voltage discontinuity
and can create a large, basis-dependent MA response.

For edges belonging to different conductors, append `different` to both mesh
commands and add `--different-conductors` to the Python command. The generator
then removes all knots on both conductor cuts from the free contour basis.
Every free trace is zero at all four cut endpoints. One additional
conductor-voltage trace is zero on conductor A, one on conductor B, and holds
conductor B at one volt through `TerminalAttributes`. The emitted version-2
process library records the two open dielectric contour paths for Maxwell
reconstruction and can be used directly by Palace's automatic matcher.

For the opposite local topology, where two outward-facing edges bound a narrow
physical metal strip, append `strip` to both mesh commands and add `--strip` to
the Python command:

```text
julia mesh/mesh_edge_pair_coupon.jl thin 2.0 \
  /tmp/edge-pair/strip_thin.msh strip
julia mesh/mesh_edge_pair_coupon.jl fabricated 2.0 \
  /tmp/edge-pair/strip_fabricated.msh strip
python3 generate_edge_pair_response.py \
  --cutout-width 2.0 --radius 2.0 --metal-thickness 0.1 --strip \
  --thin-mesh /tmp/edge-pair/strip_thin.msh \
  --fabricated-mesh /tmp/edge-pair/strip_fabricated.msh \
  --output /tmp/edge-pair/strip
```

The strip does not intersect the matching contour, so all contour knots remain
free and the strip center is the single PEC reference. Palace identifies this
topology from the two local gap directions; the two sides may belong to
different perimeter loops even though the metal between them is physically
connected.

For a paired topology, Palace first uses an entry within
`SeparationTolerance`. Otherwise it interpolates between the nearest lower and
upper separations with the same conductor-state representation. Each coupon is
evaluated on its own matching contour and the two responses are combined with
linear weights; Palace does not extrapolate. The separation-bracket span
divided by `R` enters the 3D Maxwell confidence diagnostic, so the process
library should sample rapidly changing small gaps more densely.

As an intentionally coarse check, the independently calibrated 1 and 3 um
same-conductor coupons were interpolated to the held-out 2 um coupon. Across
four smooth trace excitations, most SA/MS/MA fixed-trace and fixed-flux errors
were `0.6--12.2%`, while the domain-defect error was `9.6--19.3%`. One SA
relative error was much larger because its held-out energy was nearly zero.
This 2 um bracket has span/R equal to one and therefore fails the current
confidence threshold. The result supports interpolation as a library-coverage
mechanism, but not treating a coarse separation grid as converged calibration
data.

A held-out `SameConductorStrip` check used independently calibrated 1 and 3 um
strip coupons to interpolate the 2 um claw-shield strips in the coarse,
first-order single-transmon Maxwell model. The automatic matcher interpolated
154 longitudinal intervals with the same geometric coverage as a direct 2 um
strip coupon. Across both computed modes and all SA/MS/MA interfaces, the
interpolated corrected participation differed from the direct-coupon result by
`0.19--0.70%` for fixed trace and `0.19--0.78%` for fixed flux. This is good
evidence that strip-width interpolation works for that process and field
family, but the `span/R = 1` confidence warning remains intentionally
conservative because the same bracket was less accurate for gap coupons.
Repeating the direct 2 um strip calibration with 48 rather than 24 periodic
hat functions changed the same transmon participations by at most `0.18%` for
fixed trace and `0.25%` for fixed flux. The 24-function strip basis is therefore
adequate for these comparisons.

The generated matrices can also be configured as an explicit
`ResponseCorrection.Models` entry with one patch at the midpoint of the pair. A
separate one-edge model can still be reused for other nonoverlapping edges in the
same simulation. Every patch contributes to every target in its model's
`Interfaces` mapping.

### Building a process-library sweep

`build_surface_process_library.py` automates the straight-edge coupon workflow
for one fabrication process. It generates the isolated-edge response plus
matched thin and fabricated meshes for every paired-edge sample, writes and
runs the Palace response-matrix configurations, and packages the result with
any existing corner libraries into one portable directory. For example, from
`examples/cpw2d`:

```text
python3 build_surface_process_library.py \
  --output /tmp/al-100nm-process \
  --name al-100nm-overetch-50nm-R2um \
  --palace ../../build/bin/palace --ranks 4 \
  --matching-radius 2.0 \
  --metal-thickness 0.1 --overetch 0.05 \
  --sidewall-angle 80 --top-radius 0.01 --bottom-radius 0.01 \
  --substrate-permittivity 11.47 \
  --sa-thickness 0.002 --sa-permittivity 4.0 \
  --ms-thickness 0.002 --ms-permittivity 11.47 \
  --ma-thickness 0.002 --ma-permittivity 10.0 \
  --separations 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 \
  --base-library /path/to/corner/process-library.json
```

The final file is
`/tmp/al-100nm-process/library/process-library.json`. The build manifest and
individual coupon inputs remain under the output directory so an interrupted
build can be resumed with the same command. A coupon is reused only when its
recorded fabrication, mesh, and response-basis specification is unchanged and
all four response matrices exist. Use `--prepare-only` to generate the meshes
and Palace configurations without running the solves.

The packaged library uses schema version 3 and records the SA, MS, and MA layer
thicknesses and permittivities used by the coupon configurations. Palace checks
these values against every selected runtime dielectric interface before applying
the library. This prevents a response calibrated with one interface-layer model
from being silently reused with another.

The substrate and interface-layer options drive both the generated Palace coupon
configs and the version-3 process metadata. Changing any of them invalidates the
resumable coupon specification, so matrices generated with the old values are not
reused.

If a base library already contains an `IsolatedEdge` model, the builder reuses
it and does not generate a duplicate. `--skip-isolated` disables isolated-edge
generation explicitly, for example when preparing only an additional paired
sweep.

`--separations` applies one sweep to `SameConductorGap`,
`DifferentConductorGap`, and `SameConductorStrip`. The topology-specific
`--same-conductor-gaps`, `--different-conductor-gaps`, and
`--same-conductor-strips` options override that common sweep; pass a
topology-specific option with no values to omit that topology.

### Coverage preflight and coupon planning

Run Palace with `--surface-response-preflight` to classify a device without
assembling or solving its field equations. A direct preflight verifies coverage
against an existing library. For an empty or partial seed, first discover the
library-independent canonical closure in one user-visible operation:

```text
python3 discover_surface_response_requirements.py \
  /path/to/device-config.json \
  --output /tmp/process-requirements \
  --palace ../../build/bin/palace
```

The helper performs only cheap geometry preflights with virtual model
descriptors; it never meshes or solves a coupon. This prevents adding a corner
or spatial model from changing the residual partition and exposing new coupon
requirements after expensive generation. Its final
`surface-response-requirements.json` records exact, interpolated, and missing
coverage against the original source library. Convert those entries into a
fabrication-aware work plan with:

```text
python3 prepare_surface_response_coupons.py \
  /tmp/process-requirements/surface-response-requirements.json \
  --output /tmp/process-coupon-plan
```

The input process library may be a metadata-only version-3 seed with an empty
`Models` array. This bootstraps a new fabrication process without a hand-built
coupon. Palace accepts an empty library only for geometry preflight; an ordinary
correction run still requires at least one qualified model.

For a practical 3D classification fixture, prepare the checked-in coarse transmon
without changing its solve configuration:

```text
python3 ../transmon/prepare_surface_response_preflight.py \
  --output /tmp/transmon-surface-preflight.json
palace --surface-response-preflight /tmp/transmon-surface-preflight.json
```

The helper resolves the existing mesh path and attaches the metadata seed at
`../transmon/benchmark/transmon_surface_process_seed.json`, but does not run an eigenmode
solve or create fabrication-resolved device geometry.

The plan routes isolated and paired edges, parallel multi-edge clusters, 90
degree convex and concave corners, and exact spatial edge clusters to candidate
canonical coupon builders. It retains an explicit endpoint or junction requirement
only when the manifest provides the required classified closed mask. Unsupported
corner angles remain explicitly identified in `coupon-plan.json`.

For a PEC process library, generate, qualify, cache, and merge all missing
coupons with:

```text
python3 prepare_surface_response_coupons.py \
  /path/to/surface-response-requirements.json \
  --output /tmp/process-coupon-plan \
  --execute \
  --palace ../../build/bin/palace \
  --orders 2 3
```

Automatic execution requires version-3 fabrication metadata in microns,
including metal thickness, overetch, sidewall angle, top and bottom rounding,
substrate permittivity, and SA/MS/MA layer properties. The selected fine mesh
size must place at least two elements across the metal thickness and overetch.
Rounded curves instead use a minimum geometric sampling requirement, so a
small rounding radius does not force that size throughout the volume.

A coupon library is independent of the FEM order of the device simulation. Do
not build a complete dense library at every candidate order. First run only the
order- and mesh-convergence probes:

```text
python3 prepare_surface_response_coupons.py \
  /path/to/surface-response-requirements.json \
  --output /tmp/process-probe-study \
  --execute --probe-study-only \
  --palace ../../build/bin/palace \
  --orders 2 3 4
```

After selecting the lowest converged coupon order, build the full response
library once at that order.

Independent corner and spatial families can run concurrently. For example, on
a 192-core node, `--ranks 64 --coupon-jobs 3` permits up to three 64-rank
coupon jobs. The product of these values should not exceed the available MPI
ranks, and aggregate memory must fit on the node. Output libraries and
qualification reports retain deterministic plan order even when jobs finish in
a different order.

Each cache key includes the complete coupon geometry, fabrication process,
mesh and FEM settings, and generator source fingerprints. A six-field
thin/fabricated probe must pass the p-order
convergence limits before the full trace response matrices are computed.
Parallel-edge clusters, corners, and spatial coupons must also pass a
fixed-final-order mesh-resolution study. Their default coarse-to-fine factors
are `2 1`; use `--cluster-h-factors`, `--corner-h-factors`, or
`--spatial-h-factors` to add levels. Curved fabrication details in these
generated coupons use at least quadratic mesh geometry. The full matrices must
then pass definiteness and independent held-out trace tests. For aggregated
response matrices, Palace evaluates all basis traces once at surface quadrature
points and assembles the dense energy matrix as batched Gram products. This
replaces one complete surface-postprocessing traversal per basis pair while
preserving the localized non-aggregated fallback and CSV schema. Only qualified
entries are merged into `library/process-library.json`;
`qualification-manifest.json` records every source report. `--force`
invalidates completed solves as well as the qualification result.

Candidate spatial coupons are deliberately fail-closed. For a manifold planar
metal sheet, the perimeter has degree two everywhere: a trace cap or T/cross
branch is represented by corners, smooth chains, and possibly an exact spatial
cluster, not by a graph `Endpoint` or `Junction`. Degree-one and degree-three
perimeter vertices do not arise from an ordinary manifold sheet. Automatic
generation never invents their closure: it requires an exact closed mask and
rejects point-touch or otherwise nonmanifold masks. Exact spatial-cluster
requirements describe local edge segments and their inward metal sides. If
those half-strip descriptions overlap
for different conductors, they do not uniquely reconstruct the device's local
plan-view mask. Palace preflight therefore extracts the connected metal sheets,
clips their surface facets to the interaction neighborhood, and attaches them to
the spatial signature automatically. The mesher uses those facets as an exact
plan-view mask; an older edge-only signature still fails closed when its inferred
conductors overlap. The cache hashes a canonical boundary of the facet union, not
the source triangulation, so equivalent mesh refinement does not create a new
coupon family. The planner also aggregates requirements with the same canonical
boundary before execution. The generated spatial model stores this `PlanViewBoundary`.
During correction, Palace clips rank-local metal facets and exchanges only the
small neighborhood polygons. A masked model must match the same canonical
boundary; legacy edge-only models remain available only for unambiguous
nonoverlapping conductor strips.

The fabricated spatial mesher lofts the complete classified mask. Every physical
boundary receives the configured sidewall taper, top rounding, and substrate
overetch; artificial `ContinuationSegments` remain vertical and receive no
overetch. Polygonal samples of a circular mask chain are reconstructed as exact
OCC arcs, and trench generation uses a continuous offset-mask shell rather than
independent edge strips. Spatial coupon preparation promotes geometric mesh
order to at least two, optimizes the curved elements, and rejects a mesh with any
nonpositive scaled Jacobian. It independently gates FEM-order convergence and local
mesh-size convergence. `--spatial-h-factors` gives coarse-to-fine multipliers on
`--spatial-lc-fine`; its default `2 1` compares both thin and fabricated probes at
the maximum requested FEM order and reuses the finest completed FEM-order solve.
Probe and held-out qualification reject a coupon when either resolved response is
not stable.

For the rounded-endpoint development coupon at p=3, reducing `lc_fine` from
`0.025 um` to `0.0125 um` changed the fabricated domain matrix by `0.067%`, MA
by `13.14%`, MS by `2.39%` in Frobenius norm (`16.01%` for the worst probe
energy), and SA by `8.23%`. The geometry is valid, but this coupon does not pass
the production convergence limits and needs another h-refinement level. Fresh
quadratic meshes for the rounded endpoint, a transmon-derived L mask, and a
two-component mask had final minimum scaled Jacobians `0.1013`, `0.1007`, and
`0.1002`, respectively. These mesh checks do not substitute for response
convergence.

Execution continues across independent requirements after a candidate fails.
Qualified entries are cached and merged into a partial process library, while
`coupon-plan.json` and `qualification-manifest.json` retain every failure. The
command returns a nonzero status and `Execution.Complete` remains false until the
entire preflight manifest is covered.

Automatic finite-impedance coupon generation is not yet qualified. The planner
preserves the complete boundary-law signature, but `--execute` currently
rejects generated non-PEC coupons instead of silently substituting a PEC
calibration.

Palace never extrapolates a paired response. Each requested topology should
therefore include the smallest separation supported by the fabrication design
rules and an endpoint at `2R`. Features wider than `2R` have disjoint matching
tubes and use the isolated-edge response. The spacing between samples controls
both interpolation error and the `max library distance` confidence diagnostic.
The example's `0.5 um` spacing gives a normalized interpolation span of `0.25`
for `R = 2 um`, but the appropriate grid still needs an independent held-out
coupon convergence check for each fabrication process.

A flip-chip CPW sweep using 100 nm metal, 50 nm overetch, 80 degree sidewalls,
10 nm rounding, and `R = 2 um` tested cutout widths
`1, 2, 3, 4, 6, 10, 20, 50 um`. Widths below 4 um used coupled coupons and the
others used independent one-edge coupons. Against fabrication-resolved
references, the maximum absolute corrected errors over the sweep were:

| Interface | L1    | L2    |
|:--------- | -----:| -----:|
| SA        | 5.41% | 0.30% |
| MS        | 2.19% | 0.32% |
| MA        | 9.11% | 4.24% |

The corresponding raw thin-metal errors reached 169% for L1 MS and 61% for MA.
The weakest width-1 L1 MS thick reference required nine AMR iterations; using
only three iterations overstated its correction error by more than a factor of
four.

## Applying the calibration to 3D geometries

The calibrated coefficients depend on the local fabrication cross-section, not
on the global CPW layout. In 3D, `surface-Q-edge.csv` integrates the outside and
annular energies over the complete metal perimeter. The local field amplitude
may vary along that perimeter, but the same cross-sectional coefficients can be
applied to each local contribution and therefore to the integrated energies.

Apply the process table to any thin-metal Palace output with matching
interface properties, edge radii, fabrication stack, and FEM order:

```text
python3 local_edge_correction.py apply \
  --calibration postpro/local-edge-process.csv \
  --input ../transmon/postpro/transmon_surface_amr \
  --output ../transmon/postpro/transmon_surface_amr/surface-Q-local-corrected.csv
```

Repeated `--radius R_UM` options select an explicit averaging window. Omitting
them lets the script select the smallest AMR-stable window and retains the
complete radius sweep. A 3D mesh still needs enough resolution in each matching annulus
`R <= d < 2R`; a zero or abruptly changing annular energy indicates that the
selected radius is under-resolved. Corners, edge junctions, and nearby features
also violate the locally extruded cross-section assumption, so the radius sweep
should show a stable interval between the fabrication scale and the next
geometric feature. Components with multiple fabrication stacks should use
separate dielectric entries and edge-attribute sets for each stack. For
truncated or extruded models, use `EdgeExcludeAttributes` to omit artificial
cut-plane boundaries from the physical metal perimeter.

`"EdgeRefinement"` enforces matching-tube resolution geometrically:

```json
"EdgeRefinement": {
  "Radius": 0.2,
  "ElementsPerRadius": 3,
  "OuterRadiusFactor": 2.0,
  "CoreIndicatorWeight": 0.0
}
```

Before the first solve, Palace refines every element intersecting the `2R`
tube until its diameter is at most `R/3`. During later AMR, elements wholly
inside `d<R` receive zero field-error weight because this core is replaced by
the calibrated model. Elements crossing `R` keep full weight. This prevents
high-energy edges or layers from monopolizing refinement while weak edges
remain too coarse. The field solution and raw postprocessing are unchanged;
only mesh selection and the AMR stopping norm use the weighted indicator.
`Radius` must also be one of the interface's `EdgeDistances`.

## Local edge diagnostics

For geometry-aware corrections, add `"LocalizeEdgeEnergy": true` to an
interface dielectric entry which already specifies `EdgeAttributes` and
`EdgeDistances`. Palace then writes `surface-Q-edge-local.csv`. Each row assigns
the annular energy to the nearest physical perimeter segment and records its
one-based edge index, endpoints, and length. Each row contains both the
interface energy assigned to the edge's nearest-segment region (`p_total`), the
energy inside the matching radius (`p_in`), the interface annular energy
(`p_ann`), and electric energy in the volume tube with the same radius
(`p_bulk_ann`). Summing all local `p_total` rows for a fixed state, excitation,
interface, and radius recovers `p_surf`. Summing `p_in` recovers
`p_surf - p_out`, and summing `p_ann` recovers the corresponding value in
`surface-Q-edge.csv`.

For automatically extracted 3D edges, `p_vertex_in` is the part of `p_in` within
an along-edge distance `R` of a corner, endpoint, or junction, and
`p_bulk_vertex_ann` is the corresponding part of `p_bulk_ann`. The correction scripts
use these directly when reporting surface- and bulk-weighted unmodeled vertex fractions.
They continue to accept older CSV files and estimate those fractions from endpoint
counts and segment lengths when the exact columns are absent.

`"EdgeDistanceSmoothing": f` optionally replaces the hard cutoffs at `R` and
`2R` with cubic transitions whose relative half-width is `f`. The local and
aggregate tables use the same weights, and the sum of local `p_in` values still
recovers `p_surf - p_out`. A modest value such as `0.2` can reduce AMR
oscillation caused by quadrature points moving across a hard cutoff. The
default is `0.0`. In the CPW studies this reduction was modest compared with
the benefit of selecting an AMR-stable multi-radius window. Smoothing is
therefore optional; it does not replace choosing a matching region which is
resolved by the target mesh.

For an ideal thin-sheet edge, `|E|^2` scales as `1/r`. The volume-tube energy
therefore scales as `R`, so `p_bulk_ann / R` is a common local singular-amplitude
proxy. Unlike the interface annulus, it samples both sides of the sheet and
does not depend on selecting SA, MS, or MA as the proxy trace. Palace reuses
this volume integral when multiple dielectric entries have the same perimeter
and radius list.

For a two-sided, polarization-aware diagnostic, also set
`"EdgeFrameNormal"` from the substrate or process side toward vacuum. The
localized bulk tube is then split into `top` and `bottom` sides and local
`n`, `m`, and `l` polarizations. The experimental polarized reducer uses:

```text
MA normal       <- top_n
MS normal       <- bottom_n
SA normal       <- top_n
SA tangential   <- top_m
```

The `top_l` channel is electric field parallel to a 3D edge. It is absent from
the 2D calibration and is reported as an unmodeled fraction when applying the
calibration in 3D. A large longitudinal fraction or a target channel mixture
which differs strongly from the 2D calibration identifies a geometry where a
locally extruded 2D correction is not trustworthy. The reducer reports
`descriptor_mismatch` as the total-variation distance between the normalized
six-channel calibration and target descriptors. It is zero for identical
mixtures and approaches one for disjoint mixtures, and directly reduces the
reported model confidence.

Create and apply the experimental polarized calibration with:

```text
python3 local_edge_polarized_correction.py calibrate \
  --thin postpro/surface_calibration_thin \
  --fabricated postpro/surface_calibration_fabricated \
  --output postpro/local-edge-polarized-process.csv
python3 local_edge_polarized_correction.py apply \
  --calibration postpro/local-edge-polarized-process.csv \
  --input postpro/surface_validation_thin \
  --reference postpro/surface_validation_fabricated \
  --output postpro/surface-validation-polarized-corrected.csv
```

This decomposition does not make an inapplicable 2D edge model valid. It
reduces cross-coupling between physically different field components and
provides a sharper confidence diagnostic when a nearby conductor changes the
local field. The scalar and polarized models should both be retained during
validation until the polarized transfer has been tested on genuinely 3D
corners and junctions.

The segment ordering is deterministic for a fixed mesh, but AMR can subdivide a
geometric edge and therefore change numeric segment indices. Coordinates should
be used when matching individual segments across meshes or AMR states. For
automatic 3D extraction, `component` and `chain` provide the physical topology:
a chain joins locally straight segments and stops at a corner, endpoint, or
junction. In 2D, physical edges are points and therefore have zero segment
length. In 3D, `p_ann / L_edge` is the annular participation per unit edge
length.

The correction reducers pool automatic 3D segments with the same component
and chain before fitting. This makes the fitted response insensitive to AMR
subdivision of a physical edge and gives short segments enough samples to
identify a singular amplitude. Manual selections and legacy tables continue
to fit each numeric edge independently. Because the available process coupon
is a straight-edge model, the reducers also report
`unmodeled_vertex_length_fraction`: the fraction of total automatic chain
length within one matching radius of a corner, endpoint, or junction. The
reducers also weight this local fraction by each segment's assigned surface
energy and bulk-annulus energy. `unmodeled_vertex_fraction` is the largest of
the geometric, surface-weighted, and bulk-weighted estimates and reduces model
confidence. It does not supply a corner correction; a nonzero value means the
reported participation still applies a straight-edge response in a
geometrically unsupported neighborhood.

Validate the partition and create a compact line-density table with:

```text
python3 local_edge_diagnostics.py \
  --input ../cpw3d_surface/postpro/surface_validation_thin_L50
```

The reduced table also reports the fitted singular amplitude, smooth surface
density, separate bulk and surface fit residuals, singular fraction, and model
confidence for every physical edge segment.
Inspect their spatial distribution before trusting an aggregate correction:
isolated low-confidence segments identify corners, junctions, under-resolved
short segments, or nearby features which violate the locally straight edge
model. Legacy local tables without `p_in` can still produce this fit report,
but cannot validate the inside-energy partition.

## Geometry-aware local correction

The local bulk proxy distinguishes a true thin-sheet edge singularity from a
large smooth field caused by nearby conductors. The surface-annulus fit
separately estimates the smooth interface energy:

```text
p_bulk_ann,e(R) / R = A_e + B_e R.
p_surf_ann,e(R)     = a_e + b_e R.
```

Here `A_e` estimates the singular edge amplitude, while the term proportional
to `R` in the bulk fit is the leading smooth-volume contamination. The local
singular fraction is
`f_e(R) = clamp(A_e / (p_bulk_ann,e(R) / R), 0, 1)` and is diagnostic only.
The local fit confidence depends on the bulk and surface residuals, not on
`f_e`: a smooth-dominated edge can still have a well-resolved smooth
decomposition.

For each fabrication interface, the 2D calibration stores additive
inside-singular and outside-delta coefficients. Application replaces the raw
inside term with `max(0,b_e R) + C_in A_e` and adds `C_out A_e` to the resolved
outside term. A failed fit therefore produces a large uncertainty instead of
falling back to the divergent raw inside energy. The aggregate confidence is
correction-weighted. It combines fit confidence with singular identifiability
`f_e * q_fit` and the fraction of the reported participation contributed by
the calibrated singular terms. Thus a smooth-dominated region remains
high-confidence when its answer is mostly resolved outside/smooth energy, but
is marked low-confidence when an ill-conditioned singular intercept controls
the correction. The output records `fit_confidence`,
`singular_identifiability`, and `modeled_fraction` separately.

Create one process calibration from matched thin and fabrication-resolved 2D
simulations:

```text
python3 local_edge_correction.py calibrate \
  --thin postpro/surface_calibration_thin \
  --fabricated postpro/surface_calibration_fabricated \
  --output postpro/local-edge-process.csv
```

Apply it to a thin-metal target:

```text
python3 local_edge_correction.py apply \
  --calibration postpro/local-edge-process.csv \
  --input ../transmon/postpro/transmon_surface_amr \
  --output ../transmon/postpro/transmon_surface_amr/surface-Q-local-corrected.csv
```

Target interface indices which differ from the 2D calibration use mappings
such as `--interface-map 4=1`. Repeat the option for every nonidentity mapping.
The calibration records the FEM order and number of radii used to fit the edge
amplitude. Application rejects mismatches because the unresolved finite-element
trace, especially for MS, can depend appreciably on both choices. It also
records and validates the interface type, physical layer thickness, and
permittivity, as well as whether H(div) normal-flux recovery was enabled, so an
SA coefficient cannot silently be applied to an incompatible target interface
or field representation.

By default, the script averages the corrected participation over three
neighboring radii. A correct local replacement should be independent of the
artificial matching radius, so `radius_spread` is reported as a model
uncertainty instead of being amplified by extrapolation to `R=0`. If Palace
AMR iteration directories are available, the script selects the smallest
radius window whose mean changes by less than 2% over the final three AMR
states. If no window meets that threshold, it uses the least-varying window
and reports the residual `amr_spread`. Explicit repeated `--radius R_UM`
options override automatic selection. Large `radius_spread`, `amr_spread`, or
fit residual values are warnings that the result needs more mesh resolution
or that the local edge model is not applicable. The summary also reports how
many perimeter segments receive nonzero bulk-tube quadrature samples at the
least-resolved selected radius. A fit RMS larger than
`--fit-residual-scale` is explicitly reported as a poor local edge fit; such a
result is diagnostic, not a corrected
participation suitable for validation. When both calibration simulations save
AMR iterations, the process table also records the spread of the inside
coefficient and outside delta over its final three states. Application warns
when either calibration quantity varies by more than 5%. The reported
`relative_uncertainty` is the largest of the fit residual, half the radius
window range, the target AMR spread, the two calibration AMR spreads, and the
correction-confidence deficit.
It is an observable consistency bound, not a statistical confidence interval.

Validation should report aggregate and per-layer SA/MS/MA errors against
fabrication-resolved references, correction spread over the final AMR states,
radius-window spread, and the separate bulk/surface fit residuals. Extruded
CPW, flip-chip CPW, and genuinely three-dimensional corners test distinct
assumptions and should remain separate validation cases.
