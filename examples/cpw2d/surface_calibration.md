# Thin-metal surface-participation calibration

This study compares an infinitely thin CPW model with a matched
fabrication-resolved model containing 100 nm metal, 50 nm substrate overetch,
80 degree sidewalls, and 10 nm corner radii.

For interface `j` and matching radius `R`, define

```text
p_fab,inside,j(R) = p_fab,raw,j - p_fab,out,j(R)
S_j(R) = p_fab,out,j(R) / p_thin,out,j(R)
C_j(R) = p_fab,inside,j(R) / p_thin,ann,SA(R)
p_corrected,j(R) =
    S_j(R) p_thin,out,j(R) + C_j(R) p_thin,ann,SA(R)
```

The thin-model SA annulus is used as a shared local edge-amplitude proxy for
SA, MS, and MA. The MS and MA surface-field traces themselves converge too
irregularly under AMR to be good amplitude proxies.

The two coefficients are calibrated on the 10/6 um CPW and applied unchanged
to a separate 20/12 um CPW. Choose `R` where both coefficient tables vary
slowly, subject to both:

```text
fabrication feature size << R
2 R < distance to the next geometric feature
```

For this fabrication stack, the cross-geometry validation favors
approximately 0.2-0.3 um. Larger radii are retained in the sweep to show the
breakdown as the annulus begins to sample the macroscopic CPW geometry.

At `R = 0.2 um`, applying the 10/6 um coefficients unchanged to the 20/12 um
CPW changes the SA/MS/MA errors from `+63%/+40%/-37%` to approximately
`-1.2%/-3.2%/-3.7%`. Over the complete thin-model AMR history, the
corresponding participation spreads change from `63%/57%/25%` to
`0.8%/19%/16%`. Restricting the comparison to the last four refinements gives
corrected spreads below 6.3% for all interfaces.

Generate and run the study from this directory:

```text
cd mesh
julia generate_surface_calibration_meshes.jl
cd ..
../../build/bin/palace -np 4 cpw2d_surface_calibration_thin.json
../../build/bin/palace -np 4 cpw2d_surface_calibration_fabricated.json
../../build/bin/palace -np 4 cpw2d_surface_validation_thin.json
../../build/bin/palace -np 4 cpw2d_surface_validation_fabricated.json
python3 surface_calibration.py
```

The reducer writes the final-radius sweep to
`postpro/surface-calibration.csv`, the thin-metal AMR history to
`postpro/surface-calibration-convergence.csv`, and the blind transfer result to
`postpro/surface-calibration-validation.csv`.
