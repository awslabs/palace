# Domain decomposition preconditioning

Palace provides two domain decomposition preconditioners: a one-level restricted additive
Schwarz (RAS) method and a two-level balancing domain decomposition by constraints (BDDC)
method. Both use the existing mesh partition, so each MPI rank owns one subdomain and no
repartitioning is performed.

## Restricted additive Schwarz (RAS)

RAS is backed by hypre's RAS-ILU(k) implementation. The partition of the assembled parallel
finite element matrix defines one subdomain per MPI rank. Neighboring off-process unknowns
extend the local ILU solves, and the correction is restricted to locally owned unknowns.

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

## Balancing domain decomposition by constraints (BDDC)

BDDC is a two-level method: in addition to subdomain-local corrections it solves a small
coarse problem built from primal constraints on the subdomain interface. That coarse
correction is what keeps the iteration count nearly independent of the number of
subdomains, which a one-level method such as RAS cannot achieve.

BDDC is provided by PETSc's `PCBDDC` and is selected with:

```json
{
  "Solver": {
    "Linear": {
      "Type": "BDDC",
      "KSPType": "CG"
    }
  }
}
```

Palace supplies PETSc with the unassembled subdomain matrix, the local-to-global dof
mapping, the essential boundary dofs, and, for H(curl) problems, the discrete gradient so
that BDDC can account for the curl-curl kernel.

BDDC is symmetric, so it is used with CG.

### BDDC limitations

- Requires a Palace build with PETSc/SLEPc support.
- Requires a conforming mesh. The mapping from local to global dofs must have a single
  ±1 entry per local dof, which is not the case for nonconforming (hanging) dofs.
- Supports real-valued problems only, so electrostatic and magnetostatic simulations.
  Driven and eigenmode simulations are rejected.
- CPU execution only.
- Palace's PETSc is built with complex scalars because SLEPc requires them, so a
  real-valued BDDC solve carries the cost of complex arithmetic.
- The magnetostatic curl-curl operator is singular when it has a gradient nullspace. BDDC's
  subdomain solves are then singular as well, and a regularized preconditioner matrix is
  required. Prefer AMS for such cases until that support is added.
