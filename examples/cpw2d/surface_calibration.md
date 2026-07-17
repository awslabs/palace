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

The matching radius is a numerical regularization parameter and does not have
to be shared by SA, MS, and MA. Select and freeze one radius per interface
using independent 2D calibration and validation geometries. The calibration
and target simulations should use the same FEM order and comparable AMR
convergence because the finite-resolution surface traces, especially MS, can
otherwise shift the fitted coefficients. Never tune the radii against the 3D
fabrication-resolved result being used as a validation target.

For the order-2 study in this directory, the cross-geometry validation favors
approximately 0.2-0.3 um for all three interfaces. Larger radii are retained
in the sweep to show the breakdown as the annulus begins to sample the
macroscopic CPW geometry. An order-3 calibration can favor different frozen
radii for the three interfaces because the MS surface trace converges more
slowly than the shared SA annulus proxy.

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

## Applying the calibration to 3D geometries

The calibrated coefficients depend on the local fabrication cross-section, not
on the global CPW layout. In 3D, `surface-Q-edge.csv` integrates the outside and
annular energies over the complete metal perimeter. The local field amplitude
may vary along that perimeter, but the same cross-sectional coefficients can be
applied to each local contribution and therefore to the integrated energies.

Apply the coefficient table to any thin-metal Palace output with matching
interface indices, dielectric properties, edge radii, and fabrication stack:

```text
python3 apply_surface_calibration.py \
  --calibration postpro/surface-calibration.csv \
  --input ../transmon/postpro/transmon_surface_coarse \
  --radius SA=0.2 \
  --radius MS=0.2 \
  --radius MA=0.2 \
  --output ../transmon/postpro/transmon_surface_coarse/surface-Q-corrected.csv
```

With `--radius`, the result contains one corrected participation ratio and
quality factor per interface. Omitting the option retains the complete radius
sweep. A 3D mesh still needs enough resolution in each matching annulus
`R <= d < 2R`; a zero or abruptly changing annular energy indicates that the
selected radius is under-resolved. Corners, edge junctions, and nearby features
also violate the locally extruded cross-section assumption, so the radius sweep
should show a stable interval between the fabrication scale and the next
geometric feature. Components with multiple fabrication stacks should use
separate dielectric entries and edge-attribute sets for each stack. For
truncated or extruded models, use `EdgeExcludeAttributes` to omit artificial
cut-plane boundaries from the physical metal perimeter.
