<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Matched mesh-growth pilot results — 2026-09-10

## Contract

This study targets the mesh, not quadrature changes. All production comparisons
below use the same dedicated `soca-green-job` c7g.16xlarge, **64 MPI ranks**, one
frozen executable, p5, PCG + p-multigrid/BoomerAMG, relative tolerance 1e-8, zero
initial guesses, and identical three-source inputs. Sources are a conductor state
and affine x/y matching traces. They are exactly representable on both boundary
mesh families. Whole-family interface integration is unchanged.

Geometry: a closed strip, 100 nm metal, 50 nm overetch, 90-degree walls, zero
rounding, substrate permittivity 11.47, and coupon radius 0.5 um. The prism control
uses explicit full-gap etching to match the CAD candidate; the original finite-
collar default and retained libraries are unchanged. Material volumes and physical
interface areas agree at about 1e-13 relative.

**Scope:** complete 3-by-3 pilot matrices, not the full production nodal-trace basis
or an absolute continuum-accuracy certification. No new library is installed.

## The direct acceleration result

| Case | p5 H1 DOFs | PCG iterations, sources 1/2/3 | Worker + reducer |
|---|---:|---:|---:|
| Thin prism | 15,803,902 | 45 / 92 / 91 | 515.11 s |
| Thin baseline tet | 319,617 | 7 / 7 / 7 | 22.25 s |
| Thin refined-surface tet, ~2 nm surface altitude | 1,293,172 | 8 / 10 / 10 | 30.33 s |
| Fabricated prism | 28,468,757 | 43 / 87 / 86 | **>655.05 s** |
| Fabricated baseline tet | 750,847 | 8 / 8 / 8 | 25.77 s |
| Fabricated finest-surface tet, ~0.51 nm altitude | 10,136,797 | 10 / 11 / 11 | 129.83 s |

The fabricated prism worker hit its 600-second limit after two completed sources.
Only the missing source was re-solved, with the same essential-boundary contract;
then the complete reference matrices were reduced. Its accuracy data are complete,
but its uninterrupted runtime was not measured. **655.05 s is a lower bound**, not
a reconstructed exact time. The two source archives already completed were retained.
The earlier interrupted 32-rank run is not used in these speed comparisons.

The practically interesting asymmetric pilot pair is:

- thin ~2 nm surface tet: **17.0x** faster than the thin prism, about **12.2x fewer
  H1 DOFs**, and worst domain-energy error about **0.0250%** over the entire pilot span;
- fabricated ~0.51 nm surface tet: **>5.0x** faster than the fabricated prism,
  about **2.8x fewer H1 DOFs**, with the accuracy caveats below.

Together these two pilot generations took **160.15 s**, versus a prism lower bound
of **1170.16 s**: **>7.3x** for this small pilot basis. This is not an extrapolated
full-library wall-clock claim.

## Growth and accuracy tradeoff

The first three tet profiles have an exactly identical input surface mesh. Their
only changes are volume growth/minimum/maximum settings. The later profiles refine
near-metal surface resolution while retaining aggressive growth into the bulk.

| Fabricated tet profile | Elements | H1 DOFs | Time | MA matrix difference | MS | SA |
|---|---:|---:|---:|---:|---:|---:|
| Baseline | 33,382 | 750,847 | 25.77 s | 7.412% | 0.439% | 1.826% |
| Faster bulk growth | 31,124 | 703,837 | 26.09 s | 9.057% | 0.441% | 1.930% |
| Protected volume edges | 44,005 | 971,617 | 27.87 s | 7.131% | 0.404% | 1.582% |
| Stronger near-metal surface/volume | 122,092 | 2,688,982 | 42.44 s | 4.821% | 0.203% | 0.837% |
| ~2.07 nm surface altitude | 154,144 | 3,394,267 | 52.69 s | 2.784% | 0.111% | 0.388% |
| ~1.03 nm surface altitude | 256,013 | 5,704,837 | 76.41 s | 1.240% | 0.044% | 0.072% |
| ~0.51 nm surface altitude | 450,223 | 10,136,797 | 129.83 s | **0.344%** | **0.0034%** | **0.159%** |

Differences are full symmetric-matrix Frobenius differences against the matched
p5 prism control. Compact upper-triangular CSV matrices are reconstructed before
norms or spectral comparisons. MA was underestimated by the coarser tet meshes.

**Fewer elements alone was not a win:** faster bulk growth reduced the baseline
tet DOFs by roughly 6–7%, but increased iterations and slightly increased runtime.
It also worsened MA accuracy. Stronger surface resolution was more effective than
volume refinement alone for this coupon.

The finest fabricated mesh has worst relative energy errors over arbitrary linear
combinations of the three pilot sources of approximately:

- domain: **0.00367%**;
- MA: **0.709%**;
- MS: **0.00635%**;
- SA: **0.186%**.

Thus it is within 1% for this pilot, but **does not satisfy a universal 0.5% MA gate**.
The 0.344% MA matrix norm must not be misreported as the worst-direction error.
The prism reference's own continuum uncertainty remains separate.

## Evidence of targeted 3D coarsening

For the finest fabricated tet mesh, elements with centroids more than 0.2 um from
all process planes fell from **51,816 prisms to 2,448 tets**, while covering roughly
half the volume in both meshes. This centroid-based statistic is not an exact
geometric subdomain measure, but it shows that fine XY topology is not propagated
through far z regions.

The finest fabricated mesh has almost as many *elements* as the prism control,
but fewer p5 unknowns and much better element conditioning: maximum kappa about
**28.2 versus 167.2**, with 10–11 rather than 43–87 PCG iterations. Its peak
process-tree RSS was about **30.3 GiB versus 62.7 GiB** for the prism control.

Surface minimum and volume minimum can now differ: the finest interface uses a
0.5 nm target, while the volume minimum remains 2 nm and grows rapidly away from
the features. This is a mesh policy, not a lowered solution order or solver tolerance.

## Controls and remaining work

- Frozen surface fingerprints match before/after volume generation, optimization,
  and serialization. Independent repeated baseline meshes were byte-identical.
- Source-timing instrumentation on/off produced byte-identical archives.
- Matched mesh/process/excitation hashes and executable/rank/host identity are
  checked by the benchmark tools.
- The efficient asymmetric thin/fabricated pair's domain-defect difference is
  approximately 0.088% relative to the prism defect, or 0.0081% relative to the thin
  domain matrix. This does not prove stability of the assembled device correction.
- Thin raw edge-inclusive SPR is not used as an accuracy gate.
- Full production trace excitations, larger/more complex coupons, reference h/p
  convergence, and final-device observables still require qualification.

Evidence: `/tmp/coupon-volume-growth-20260910/matched-final-analysis.json`
and the `results64-*` directories locally; the corresponding
`coupon_volume_growth_*_20260910` directories on the dedicated host. No existing
reference library or benchmark mesh was overwritten.
