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
20/12 um CPW and fabrication stack, a 50 um extrusion, and two structured
layers along the propagation direction. The length keeps the artificial
wave-port endpoint neighborhoods below the current 10% confidence limit while
remaining electrically short at 5 GHz.

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
  --order 2
```

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
errors separately. On the generated mesh at order 2 with no AMR, the corrected
SA/MS/MA participation errors relative to the same-order fabricated run were
`+3.37%/-1.93%/+14.90%`, compared with raw errors of
`-10.48%/-21.15%/-59.85%`. The confidence diagnostics passed: the matched
fraction was one, the loop residual was below 0.05, and the artificial endpoint
neighborhood fraction was 0.08.

The order-2 fabricated MA result remains about 11% below the converged 2D
cross-section value. It is therefore not yet a converged three-dimensional
reference. Repeat the corrected and fabricated runs at order 3 before using
this case as an accuracy gate. A remote runner must make the process-library
JSON and its matrix files visible at the exact path stored in the generated
configuration; staging only the Palace config and mesh is insufficient.
