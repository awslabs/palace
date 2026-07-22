# Three-dimensional corner response coupon

This directory generates a local electrostatic response model for a 90 degree
convex corner. The coupon frame is the frame used by Palace's automatic corner
placement:

- the corner is at the origin;
- the two incident metal-edge arms point along positive `x` and `y`; and
- positive `z` points from the substrate into vacuum.

The fabricated prototype uses 100 nm metal, 50 nm substrate overetch, and 80
degree sidewalls. The thin and fabricated geometries can use either a sharp
plan-view corner or a circular fillet with a specified mask radius.

## Generate the prototype

From `examples/cpw3d_surface`:

```sh
mkdir -p corner_coupon/mesh
julia corner_coupon/mesh_corner_coupon.jl \
  thin corner_coupon/mesh/corner_thin.msh 0.5
julia corner_coupon/mesh_corner_coupon.jl \
  fabricated corner_coupon/mesh/corner_fabricated.msh 0.5

python3 corner_coupon/generate_corner_response.py \
  --output corner_coupon/calibration \
  --thin-mesh corner_coupon/mesh/corner_thin.msh \
  --fabricated-mesh corner_coupon/mesh/corner_fabricated.msh \
  --corner-radius 0.5 \
  --ring-size 8 \
  --order 1

../../build/bin/palace --serial corner_coupon/calibration/thin.json
../../build/bin/palace --serial corner_coupon/calibration/fabricated.json

python3 corner_coupon/finalize_corner_response.py \
  corner_coupon/calibration

../../build/bin/palace --serial \
  corner_coupon/calibration/heldout-thin.json
../../build/bin/palace --serial \
  corner_coupon/calibration/heldout-fabricated.json

python3 corner_coupon/finalize_corner_response.py \
  corner_coupon/calibration
```

The mesh generator imports `Gmsh.jl`. Use an instantiated Julia environment
containing that package. Omit the final mesh-generator argument and
`--corner-radius` to generate a sharp corner.

For the default process, the matching surface has nine closed square contours with
eight trace knots per contour. The resulting 72 coefficients include every corner
of the box matching surface and contours at the trench bottom, thin-metal plane,
and metal top. The coupon mesh resolves the intersections of these process planes
with the matching surface. Pass
`--metal-thickness` and `--overetch-depth` when those process dimensions differ
from the defaults. Each prescribed
potential is stored as a triangulated `x,y,z,V,triangle` surface trace.

`finalize_corner_response.py` checks the domain and surface matrices, aggregates
the per-edge response localized to the union of the physical-edge radius-`R`
tubes into compact matrices, and evaluates a smooth
held-out boundary excitation defined on a finer matching-surface
triangulation.

For a rounded convex corner, the generated library uses `(rho, rho, 0)` as the
local PEC reference. The virtual sharp corner at the origin is outside the
rounded conductor and must not be used to gauge the coupon trace.

At runtime Palace joins connected curved perimeter facets, estimates
`rho = arc length / total turn`, and matches both corner angle and radius. A
matched corner replaces the fillet and the remaining distance `R - rho` along
both straight arms. If a sub-`R` fillet has no compatible radius-aware model,
Palace retains straight-edge response through that neighborhood and prints a
warning.

## Calibration cost and checks

For `R = 2 um` and `rho = 0.5 um`, the first-order local meshes contain 174,115
thin and 214,641 fabricated tetrahedra. On a local Apple M3 Pro, the two
72-source calibrations took 584 and 656 seconds when run concurrently. Peak
memory was 2.5 and 2.8 GB. There are no AMR iterations, and BoomerAMG is used
for the electrostatic solves.

For the current meshes, the held-out excitation gave:

| Coupon | Domain energy error | Maximum SA/MS/MA energy error |
|:--|--:|--:|
| Thin | `2.13e-4` | `4.23e-5` |
| Fabricated | `2.13e-4` | `2.43e-5` |

These errors validate the spatial trace representation and response-matrix
assembly for one smooth excitation. They do not validate the physical corner
correction in a device.

## Rectangular-island validation

The scripts in `../corner_validation` compare a corrected zero-thickness
rounded rectangular metal island with an independently meshed
fabrication-resolved island. The process library combines the `rho = 0.5 um`
corner entry with the existing straight-edge entry. For a second-order device
solve, Palace detected all four fillets automatically. The participation-ratio
errors relative to the fabricated reference were:

| Interface | Raw | Fixed trace | Fixed flux | Self-consistent |
|:--|--:|--:|--:|--:|
| SA | `-1.99%` | `-6.09%` | `-1.19%` | `+2.24%` |
| MS | `-2.10%` | `-6.13%` | `+2.97%` | `-10.97%` |
| MA | `-53.03%` | `+1.72%` | `-13.09%` | `-17.16%` |

The fixed-trace/fixed-flux closure spread is 9.0--12.6%. Maxwell correction is
postprocessing-only, so those two columns, together with the reported closure
spread, are the relevant current estimate and confidence diagnostic.

## Remaining validation

Before adding these matrices to a production process library:

1. Add any missing cross-sectional process rounding and fabrication details.
2. Repeat at higher FEM order and with refined fabricated meshes.
3. Check matching-radius and trace-basis convergence.
4. Sweep plan-view corner radii and validate radius interpolation or matching.
5. Generate and validate the corresponding concave-corner coupon.
6. Combine the corner entry with the process's isolated- and paired-edge
   entries.
7. Validate additional corner angles and nearby-edge configurations.

The generated `process-library.json` contains only the convex-corner entry. It
is therefore not by itself a complete library for a practical device.
