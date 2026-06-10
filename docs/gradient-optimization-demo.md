# Gradient-Based Shape Optimization with Palace

**Date**: 2026-02-23 | **Branch**: `hughcars/derivative-propagation`

## Result

Optimized a cylindrical PEC cavity radius to match a target eigenfrequency of 3.000 GHz. **Converged in 3 iterations** using Newton steps with analytic mesh sensitivity gradients:

| Iter | R (cm)   | f (GHz)  | Residual   | df/dR (GHz/cm) |
|:---- |:-------- |:-------- |:---------- |:-------------- |
| 0    | 2.740000 | 2.975631 | -24.4 MHz  | -1.0860        |
| 1    | 2.717561 | 3.000201 | +0.20 MHz  | -1.0950        |
| 2    | 2.717744 | 2.999998 | -0.002 MHz | -1.0949        |

The gradient `df/dR * R / f = -1.000000` matches the analytical TM₀₁₀ formula to 6 digits.

## How It Works

The sensitivity pipeline differentiates Palace's CurlCurlMass eigenvalue integrand with respect to the element Jacobian at each quadrature point, using **forward-mode automatic differentiation** via a scalar dual-number type. The chain rule contraction with mesh shape function gradients is handled natively by libCEED's `CEED_EVAL_GRAD` transpose. No finite differences are used — the gradient is exact.

```
Eigensolve → Eigenvector → QFunction AD → df/dJ at quad pts → B^T contraction → df/dx per node → chain with dx/dR → df/dR
```

Each optimization iteration runs one Palace eigensolve plus one CeedOperator pass (same cost as a single K·u evaluation). The mesh is morphed by scaling (x,y) coordinates proportionally to R.

## Technical Details

  - **Differentiation backend**: Dual-number forward-mode AD (9 passes per quad point). Enzyme reverse-mode AD also supported as a compile-time option (`-DPALACE_WITH_ENZYME=ON`).
  - **Integrand**: Matches `f_apply_hdivmass_33` exactly — H(div) Piola stiffness (`J/detJ`) and H(curl) Piola mass (`adj(J)^T/detJ`). Verified against FD at 1e-8 in standalone spike test.
  - **Eigenvector phase rotation**: The eigenvalue solver returns a complex phase; an automatic rotation maximizes the real part before sensitivity evaluation.
  - **Nondimensionalization**: Palace's internal `Lc/L0` mesh scaling is handled in the conversion factor.
  - **Scope**: Lossless linear eigenvalue problem, 3D, single process. Quadratic EVP and MPI support are follow-ups.

## Files

| File                                        | Role                                                                   |
|:------------------------------------------- |:---------------------------------------------------------------------- |
| `palace/fem/qfunctions/33/meshsens_33_qf.h` | Dual-number type, templated integrands, Enzyme/dual-number dispatch    |
| `palace/fem/integ/meshsensitivity.cpp`      | Custom CeedOperator: active mesh nodes, passive eigenvector            |
| `palace/models/meshsensitivity.cpp`         | Model layer: coefficient construction, unit conversion, phase rotation |
| `examples/cylinder/optimize_radius.py`      | Optimization loop: Palace → gradient → mesh morph → repeat             |
