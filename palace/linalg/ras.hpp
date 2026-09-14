// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_RAS_HPP
#define PALACE_LINALG_RAS_HPP

#include <mfem.hpp>

namespace palace
{

//
// A wrapper for hypre's one-level restricted additive Schwarz (RAS) preconditioner with
// ILU(k) subdomain solves. The parallel matrix partition defines one subdomain per MPI
// rank and hypre extends each local solve with the neighboring off-process unknowns.
//
// RAS is a nonsymmetric, one-level domain decomposition preconditioner, so it requires
// GMRES or FGMRES and, without a coarse correction, its iteration count grows with the
// number of subdomains. CPU only: Palace builds GPU hypre without the unified memory
// required by RAS-ILU.
//
class RasSolver : public mfem::HypreILU
{
public:
  RasSolver(int fill_level = 1, int print = 0);

  void SetOperator(const mfem::Operator &op) override;
};

}  // namespace palace

#endif  // PALACE_LINALG_RAS_HPP
