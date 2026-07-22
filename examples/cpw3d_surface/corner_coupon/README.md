# Three-dimensional corner response coupon

This directory generates a local electrostatic response model for a 90 degree
convex corner. The coupon frame is the frame used by Palace's automatic corner
placement:

- the corner is at the origin;
- the two incident metal-edge arms point along positive `x` and `y`; and
- positive `z` points from the substrate into vacuum.

The fabricated prototype uses 100 nm metal, 50 nm substrate overetch, and 80
degree sidewalls. The current geometry has sharp process corners; fabrication
rounding has not yet been added.

## Generate the prototype

From `examples/cpw3d_surface`:

```sh
mkdir -p corner_coupon/mesh
julia corner_coupon/mesh_corner_coupon.jl \
  thin corner_coupon/mesh/corner_thin.msh
julia corner_coupon/mesh_corner_coupon.jl \
  fabricated corner_coupon/mesh/corner_fabricated.msh

python3 corner_coupon/generate_corner_response.py \
  --output corner_coupon/calibration \
  --thin-mesh corner_coupon/mesh/corner_thin.msh \
  --fabricated-mesh corner_coupon/mesh/corner_fabricated.msh \
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
containing that package.

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

## Calibration cost and checks

The first-order local meshes contain 235,780 thin and 207,361 fabricated
tetrahedra. On a local Apple M3 Pro, the two 72-source calibrations took 728 and
631 seconds when run concurrently. Peak memory was about 3.5 GB for each run.
There are no AMR iterations, and BoomerAMG is used for the electrostatic solves.

For the current meshes, the held-out excitation gave:

| Coupon | Domain energy error | Maximum SA/MS/MA energy error |
|:--|--:|--:|
| Thin | `2.11e-4` | `2.05e-5` |
| Fabricated | `2.05e-4` | `9.90e-6` |

These errors validate the spatial trace representation and response-matrix
assembly for one smooth excitation. They do not validate the physical corner
correction in a device.

## Rectangular-island validation

The scripts in `../corner_validation` compare a corrected zero-thickness
rectangular metal island with an independently meshed fabrication-resolved
island. The process library combines this corner entry with the existing
straight-edge entry. For a second-order device solve, Palace detected all four
corners and replaced the complete radius-`R` neighborhood on both incident
arms. The participation-ratio errors relative to the fabricated reference were:

| Interface | Raw | Fixed trace | Fixed flux | Self-consistent |
|:--|--:|--:|--:|--:|
| SA | `-3.32%` | `-5.26%` | `+0.05%` | `+6.32%` |
| MS | `-5.28%` | `-5.25%` | `+3.72%` | `+1.77%` |
| MA | `-53.88%` | `+3.31%` | `-10.88%` | `-15.09%` |

The fixed-trace/fixed-flux closure spread is 9.7--12.8%. Maxwell correction is
postprocessing-only, so those two columns, together with the reported closure
spread, are the relevant current estimate and confidence diagnostic.

## Remaining validation

Before adding these matrices to a production process library:

1. Add the actual process rounding and any other missing fabrication details.
2. Repeat at higher FEM order and with refined fabricated meshes.
3. Check matching-radius and trace-basis convergence.
4. Generate and validate the corresponding concave-corner coupon.
5. Combine the corner entry with the process's isolated- and paired-edge
   entries.
6. Validate additional corner angles and nearby-edge configurations.

The generated `process-library.json` contains only the convex-corner entry. It
is therefore not by itself a complete library for a practical device.
