# Extruded CPW surface-participation validation

This study checks the thin-metal surface correction in a three-dimensional
mesh against an independently fabrication-resolved reference. The CPW
cross-section is extruded with a structured mesh, so the converged 3D
electrostatic result should also reproduce the corresponding 2D result.

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

On the generated mesh at order 2 with no AMR, the fixed-trace SA/MS/MA
participation errors relative to the same-order fabricated run were
`+3.37%/-1.93%/+14.90%`, while the fixed-flux errors were
`+10.13%/+1.66%/+7.64%`. Raw errors were
`-10.48%/-21.15%/-59.85%`. The maximum closure spread was `8.27%`, so the
closure diagnostic correctly rejects this result even though the matched
fraction is one, the loop residual is below 0.05, and the artificial endpoint
neighborhood fraction is 0.08.

At order 4, the 50 um corrected SA/MS/MA errors were
`+0.06%/-0.78%/+12.15%`, compared with raw errors of
`-0.37%/-5.67%/-53.13%`. Increasing the line length to 200 um while preserving
the 25 um axial mesh spacing gave corrected errors of
`-0.62%/-0.76%/+11.75%`. The corrected values at the two lengths agree within
`0.01%` for every interface, and the 200 um confidence diagnostics report a
matched fraction of one, loop residual of `4.8e-5`, and endpoint-neighborhood
fraction of 0.02. The persistent MA excess is therefore not caused by the
artificial endpoints or insufficient axial resolution. These order-4 runs
predate the fixed-flux diagnostic. The order-2 closure dependence shows that
the remaining error is specifically a postprocessing-only trace-closure
limitation, rather than evidence that the coupon response matrix is inaccurate
for a prescribed trace. Repeat the order-4 corrected run with the current
executable before using MA as an accuracy gate.

The order-3 to order-4 changes in the 50 um fabricated SA/MS/MA references were
`+0.99%/+1.45%/+2.37%`; account for that remaining reference-discretization
uncertainty when interpreting the errors above. Do not compare this compact
Maxwell case directly to the large-domain 2D electrostatic result: its ground,
substrate, and vacuum extents are intentionally much smaller.

A remote runner must make the process-library JSON and its matrix files visible
at the exact path stored in the generated configuration; staging only the
Palace config and mesh is insufficient.
