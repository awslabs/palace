# Domain decomposition preconditioning

Palace provides a one-level restricted additive Schwarz (RAS) preconditioner backed by
hypre's RAS-ILU(k) implementation. The partition of the assembled parallel finite element
matrix defines one subdomain per MPI rank. Neighboring off-process unknowns extend the
local ILU solves, and the correction is restricted to locally owned unknowns.

RAS is selected through the linear solver configuration:

```json
{
  "Solver": {
    "Linear": {
      "Type": "RAS",
      "KSPType": "GMRES",
      "RASFillLevel": 1
    }
  }
}
```

`RASFillLevel` controls the fill level ``k`` in each ILU(k) subdomain solve. Larger values
may improve convergence at the cost of additional setup time and memory. On a serial or
globally block-diagonal matrix, Palace uses block ILU because there is no subdomain
interface.

## Limitations

- RAS is nonsymmetric and must be used with GMRES or FGMRES. Palace selects GMRES when
  `KSPType` is omitted and rejects incompatible Krylov methods.
- This is a one-level method without a global coarse correction. Its iteration count is
  therefore not expected to remain constant as the number of subdomains increases.
- RAS currently supports CPU execution only. Palace builds GPU hypre without the unified
  memory required by hypre's RAS-ILU implementation.
- AMS remains the recommended choice for magnetostatic curl-curl systems. RAS is **not**
  suitable there: when the curl-curl operator has a gradient nullspace, ILU-based RAS
  stalls and the Krylov solver fails to converge, producing incorrect results. Use RAS for
  well-posed positive-definite systems such as electrostatics.
