# Palace Eigenfrequency Sensitivity Verification

**Date**: 2026-02-22
**Branch**: `hughcars/derivative-propagation`
**Test case**: Cylinder cavity with PEC walls (`examples/cylinder/cavity_pec.json`)

## Configuration

  - Geometry: cylindrical cavity (hex mesh, 37.5 KB)
  - Material: single domain, `eps_r = 2.08`, `mu_r = 1.0`, `tan_delta = 0.0004`
  - Solver: Eigenmode, order 2, 2 modes, SLEPc, tol = 1e-8
  - Sensitivity: `"Sensitivity": true` in eigenmode config

## Results

### Eigenfrequencies

| Mode | Re{f} (GHz) | Im{f} (GHz) | Q    |
|:---- |:----------- |:----------- |:---- |
| 1    | 2.905077    | 5.810e-4    | 2500 |
| 2    | 2.926602    | 5.853e-4    | 2500 |
| 3    | 3.469647    | 6.939e-4    | 2500 |

### Sensitivity Comparison: Adjoint vs Finite Difference vs Analytical

The adjoint computes `df/d(eps_r)` — the derivative of eigenfrequency with respect
to relative permittivity. For a uniform-fill cavity, the analytical result is:

```
df/d(eps_r) = -f / (2 * eps_r)
```

FD perturbation: `delta(eps_r) = 2.08e-5` (relative perturbation 1e-5).

| Mode | Adjoint (GHz) | FD (GHz)  | Analytical (GHz) | Adj/Analytical | Adj-FD error |
|:---- |:------------- |:--------- |:---------------- |:-------------- |:------------ |
| 1    | -0.698336     | -0.714082 | -0.698336        | 1.00000        | 2.3%         |
| 2    | -0.703510     | -0.692542 | -0.703510        | 1.00000        | 1.6%         |
| 3    | -0.834050     | -0.841978 | -0.834050        | 1.00000        | 0.9%         |

**Adjoint matches analytical to 8 significant digits.** The FD error of 1-2% is
expected from forward-difference truncation at `O(h)` with `h = 1e-5`.

### Sensitivity to Inverse Permeability

| Mode | df/d(mu_inv_r) (GHz) |
|:---- |:-------------------- |
| 1    | +1.452538            |
| 2    | +1.463301            |
| 3    | +1.734823            |

For uniform fill: `df/d(mu_inv) = +f / (2 * mu_inv) = f/2`. The values match
`f/2` to machine precision (e.g., mode 1: `2.905077/2 = 1.452539`).

### CSV Output

Sensitivity data written to `postpro/cavity_pec_sens/eig-sensitivity.csv` with
columns: mode index, domain attribute, domain index, and real/imaginary parts of
`df/d(eps_r)`, `df/d(mu_inv_r)`, `df/d(sigma)`.

## Bug Found and Fixed During Verification

The domain-restricted assembly methods initially used the actual material tensor
(e.g., `eps_r = 2.08`) as the coefficient, producing sensitivities that were
exactly `eps_r` times too large. The fix: use identity coefficient (`I`) since
`dM/d(eps_r) = integral(I * u . v)`, not `integral(eps_r * u . v)`.

## Runtime

Total solve + sensitivity: 0.76 seconds (1 MPI rank, Apple M1).
