# Extruded CPW surface-participation validation

This study checks the thin-metal surface correction in a three-dimensional
mesh against an independently fabrication-resolved reference. The CPW
cross-section is extruded with a structured mesh, so the converged 3D
electrostatic result should also reproduce the corresponding 2D result.

[`corner_coupon/`](corner_coupon/README.md) contains the separate prototype
workflow for generating a spatial three-dimensional convex-corner response
matrix.

The geometry is the held-out 20/12 um CPW from the 2D validation study. The
fabrication stack is 100 nm metal, 50 nm substrate overetch, 80 degree
sidewalls, and 10 nm corner radii. Correction coefficients must remain fixed
at the values calibrated on the separate 10/6 um 2D CPW.

Generate the meshes from the repository root:

```text
julia --project=examples examples/cpw3d_surface/mesh/generate_validation_meshes.jl
```

The generator writes thin-metal extrusions at 50 and 200 um, plus a
fabrication-resolved 50 um extrusion. The thin configurations exclude the
front and back boundary attributes from the physical metal-edge model. The
50 and 200 um results should therefore agree without a length extrapolation;
their comparison verifies that artificial extrusion-end edges were removed.
The fabrication-resolved participation is already length-independent for the
translationally invariant electrostatic solution.

Run the 50 um pair and apply the frozen 2D calibration:

```text
build/bin/palace -np 4 examples/cpw3d_surface/cpw3d_surface_validation_thin.json
build/bin/palace -np 4 examples/cpw3d_surface/cpw3d_surface_validation_fabricated.json
python3 examples/cpw2d/apply_surface_calibration.py \
  --calibration examples/cpw2d/postpro/surface-calibration.csv \
  --input examples/cpw3d_surface/postpro/surface_validation_thin_L50 \
  --radius SA=0.2 \
  --radius MS=0.2 \
  --radius MA=0.2 \
  --output examples/cpw3d_surface/postpro/surface_validation_thin_L50/surface-Q-corrected.csv
```

Compare the corrected thin result against both the 3D fabricated result and
`examples/cpw2d/postpro/surface_validation_fabricated/surface-Q.csv`.

The generated base thin mesh is intended to exercise the 3D diagnostics, not
to provide a mesh-converged thin reference. With order 3 and no AMR, its
edge-decomposed MS and MA energies agree with the first 2D AMR level rather
than the final 2D result. At `R = 0.2 um`, the raw thin SA/MS/MA errors relative
to the 3D fabrication-resolved result are approximately
`+17.5%/+18.3%/-42.2%`; the corrected errors are
`+2.6%/+16.8%/+14.8%`. Refine the thin mesh until the corrected result at each
preselected matching radius is stable. Radius selection may differ by
interface, but it must be frozen using an independent 2D cross-geometry
validation at the same FEM order. Do not tune any radius using this 3D
fabrication-resolved reference.

The extrusion checks are independent of this cross-section refinement:

- The 3D fabrication-resolved SA/MS/MA values agree with the converged 2D
  values within `0.12%/0.84%/1.50%`.
- The thin 50 and 200 um extrusions agree within `0.0041%` for all raw,
  edge-decomposed, and corrected participation values. This verifies that
  `EdgeExcludeAttributes` removes the artificial extrusion-end edges.

Both committed configurations use `Refinement.MaxIts = 0`, so their error
indicators are diagnostic only. They use a relaxed `EstimatorTol` to avoid
oversolving an estimator that does not drive AMR. Restore a strict estimator
tolerance before enabling refinement.

When electrostatic response correction is enabled, `surface-Q-corrected.csv`
reports two postprocessing-only closures evaluated on the unchanged raw
thin-metal potential in addition to the historical raw result and the
self-consistent corrected solve. The fixed-trace closure preserves the
potential on the coupon contour; the fixed-flux closure transforms that trace
to preserve the coupon flux. Their per-interface difference is reported as
`trace closure spread`.

On the order-3, 20 um test extrusion, errors relative to the converged 2D
fabricated SA/MS/MA values were `+0.11%/-2.51%/+6.98%` for postprocessed
fixed-trace, `+6.60%/+1.06%/+0.66%` for postprocessed fixed-flux, and
`+1.53%/-0.82%/+4.35%` for the self-consistent corrected solve. The last set
exactly reproduces the result obtained before the postprocessing-only columns
were added. The maximum trace-closure spread was 8.39%, so neither
postprocessing-only closure passes the 5% confidence threshold.

## Driven Maxwell validation

The compact wave-port case exercises the postprocessing-only Maxwell path
against a fabrication-resolved driven reference. It uses the same held-out
20/12 um CPW and fabrication stack and 50 and 200 um extrusions. Structured
layers preserve a 25 um axial spacing along the propagation direction. The
50 um line keeps the artificial wave-port endpoint neighborhoods below the
current 10% confidence limit; the 200 um line reduces their fraction from
0.08 to 0.02 so a length sweep can distinguish endpoint error from
straight-edge response error. Both lines remain electrically short at 5 GHz.

Generate the first-order thin and fabricated meshes:

```text
julia --project=examples \
  examples/cpw3d_surface/mesh/generate_maxwell_validation_meshes.jl
```

Prepare raw thin, corrected thin, and fabricated-reference configurations.
`--library` must name the independently generated process library for the
100 nm metal, 50 nm overetch, 80 degree sidewalls, 10 nm rounding, and 2 um
matching radius:

```text
python3 examples/cpw3d_surface/prepare_maxwell_validation.py \
  --library /path/to/process-library-3d.json \
  --output /tmp/cpw3d-maxwell-validation \
  --order 2 \
  --length 50
```

Use `--length 200` to prepare the endpoint-sensitivity pair.

Run the three generated configurations and compare them:

```text
build/bin/palace -np 2 /tmp/cpw3d-maxwell-validation/thin_raw.json
build/bin/palace -np 2 /tmp/cpw3d-maxwell-validation/thin_corrected.json
build/bin/palace -np 2 /tmp/cpw3d-maxwell-validation/fabricated_reference.json
python3 examples/cpw3d_surface/compare_maxwell_validation.py \
  --baseline /tmp/cpw3d-maxwell-validation/thin_raw \
  --corrected /tmp/cpw3d-maxwell-validation/thin_corrected \
  --fabricated /tmp/cpw3d-maxwell-validation/fabricated_reference
```

The comparison first verifies that enabling correction did not alter any raw
output, then reports normalization, surface-energy numerator, and participation
errors separately. `p_surf corrected` retains the original fixed-trace
(Dirichlet) response. The additional `p_surf corrected fixed-flux` column uses
the local fixed-flux (Neumann) response. Neither closure is selected as
universally correct: their difference measures how strongly the unresolved
fabricated field depends on a closure assumption that a postprocessing-only
method cannot recover from the thin-metal field. The per-interface spread is
reported in `surface-Q-corrected.csv`; its maximum is reported in
`surface-response-confidence.csv` and fails confidence above 5%.

The fabrication-resolved mesh assigns separate physical groups to trace and
ground MS and MA surfaces. On freshly generated 200 um meshes at order 1 with
no AMR, fixed-trace SA/MS/MA participation errors relative to the split
fabricated reference are `+3.37%/+3.24%/+22.63%`; fixed-flux errors are
`+9.92%/+7.03%/+22.83%`, and self-consistent errors are
`+4.64%/+4.71%/+22.29%`. Raw errors are
`-20.11%/-35.79%/-64.74%`. The maximum closure spread is `8.02%`, so the
confidence gate correctly rejects the postprocessing-only values. Geometric
coverage is complete apart from an endpoint-neighborhood fraction of 0.02.

On the 50 um split meshes at order 2, fixed-trace errors are
`+3.38%/-1.93%/+7.65%`, fixed-flux errors are
`+10.12%/+1.68%/+7.80%`, and self-consistent errors are
`+4.71%/-0.50%/+7.13%`. Raw errors are
`-10.48%/-21.15%/-59.85%`. The maximum closure spread is `8.26%`, and the
shorter line has an endpoint-neighborhood fraction of 0.08. The order-2
Maxwell MA error agrees with the order-2 2D electrostatic error below, so
Maxwell postprocessing itself does not add a separate MA discrepancy.

The large order-1 MA error is primarily a discretization effect, not an MS/MA
grouping artifact. A compact 2D electrostatic p-refinement study on the same
cross-section and response library gives fixed-trace SA/MS/MA errors:

| Order | SA | MS | MA |
| ---: | ---: | ---: | ---: |
| 1 | `+2.69%` | `+3.56%` | `+22.76%` |
| 2 | `+1.09%` | `-0.70%` | `+7.71%` |
| 3 | `-0.78%` | `+0.41%` | `+6.27%` |
| 4 | `-0.32%` | `-0.84%` | `+4.66%` |
| 5 | `-1.19%` | `-0.45%` | `+4.14%` |
| 6 | `-0.77%` | `-1.18%` | `+3.31%` |

The corrected MA value is stable to less than 0.5% between consecutive orders
above order 2, while the fabrication-resolved reference is still converging.
The reference changes by 0.55% from order 5 to 6, so the remaining asymptotic
MA bias is smaller than the order-6 error and is likely a few percent. Over the
same sweep, raw MA remains 40-50% low and raw SA/MS drift upward.

Two independent resolution checks do not explain the residual. Increasing the
isolated-edge trace basis from 48 to 96 functions changes fixed-trace SA/MS/MA
by `-0.05%/-0.03%/+0.78%`. Refining the entire 2 um matching tube to about
0.17 um element size changes the order-4 result from
`-0.32%/-0.84%/+4.66%` to `+0.06%/+0.32%/+5.66%`. The production mesh is
therefore adequate at the matching contour, and refining the full tube by
default would add substantial cost without removing the MA bias.

This supports the response correction but leaves a small MA model discrepancy
to investigate. A third-order or higher 3D Maxwell run using the split
fabricated mesh remains necessary before setting a converged Maxwell accuracy
gate. Do not compare this compact Maxwell case directly to the large-domain 2D
electrostatic result: its ground, substrate, and vacuum extents are
intentionally much smaller.

A remote runner must make the process-library JSON and its matrix files visible
at the exact path stored in the generated configuration; staging only the
Palace config and mesh is insufficient.

## Eigenmode Maxwell validation

The same 200 um extrusion also provides a compact standing-wave eigenmode
check. The front and back faces are left as PMC boundaries instead of being
modeled as metal, and their attributes are excluded from the automatic edge
selection. This avoids introducing artificial PEC-cap edges with a fabrication
mapping unrelated to the CPW metal. Prepare the configurations with:

```text
python3 examples/cpw3d_surface/prepare_maxwell_validation.py \
  --library /path/to/process-library-3d.json \
  --output /tmp/cpw3d-eigenmode-validation \
  --problem eigenmode \
  --eigenmode-target 100 \
  --order 1 \
  --length 200
```

Run the same three generated configurations, then compare the first common
mode:

```text
python3 examples/cpw3d_surface/compare_maxwell_validation.py \
  --baseline /tmp/cpw3d-eigenmode-validation/thin_raw \
  --corrected /tmp/cpw3d-eigenmode-validation/thin_corrected \
  --fabricated /tmp/cpw3d-eigenmode-validation/fabricated_reference \
  --mode 1
```

With the 2 um isolated-edge library and freshly generated split first-order
meshes, the first thin and fabricated quasi-TEM modes are 290.96 and
292.83 GHz. Raw SA/MS/MA participation errors are
`-24.91%/-35.36%/-66.47%`; fixed-trace correction gives
`-2.95%/+0.45%/+20.42%`, and fixed-flux correction gives
`+3.36%/+4.13%/+20.54%`. The maximum closure spread is 8.02%, so Palace
correctly fails the confidence gate. This case verifies eigenmode extraction,
automatic PEC edge matching, and corrected-output plumbing, but the order-1
MA result is not an accuracy pass. Eigenmode correction remains
postprocessing-only and therefore has no self-consistent columns.
