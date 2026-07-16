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
  --output examples/cpw3d_surface/postpro/surface_validation_thin_L50/surface-Q-corrected.csv
```

Compare the corrected thin result against both the 3D fabricated result and
`examples/cpw2d/postpro/surface_validation_fabricated/surface-Q.csv`.
