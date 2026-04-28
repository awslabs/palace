// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_MODE_ASSEMBLY_HPP
#define PALACE_MODELS_MODE_ASSEMBLY_HPP

#include <complex>
#include <memory>
#include <tuple>
#include <mfem.hpp>

// Frequency-independent and frequency-dependent assembly of the 2D boundary-mode GEP
// blocks, shared by BoundaryModeOperator (2D domain path) and ModeEigenSolver's bare
// ctor (3D wave port submesh path).

namespace palace
{

class ComplexVector;
class FarfieldBoundaryOperator;
class FiniteElementSpace;
class MaterialOperator;
class SurfaceConductivityOperator;
class SurfaceImpedanceOperator;

namespace mode_assembly
{

using ComplexHypreParMatrix = std::tuple<std::unique_ptr<mfem::HypreParMatrix>,
                                         std::unique_ptr<mfem::HypreParMatrix>>;

// Atn = -(mu^{-1} grad_t u, v). ND / H1 gradient coupling, real-only.
ComplexHypreParMatrix AssembleAtn(const FiniteElementSpace &nd_fespace,
                                  const FiniteElementSpace &h1_fespace,
                                  const MaterialOperator &mat_op);

// Btt = (mu^{-1} u, v). ND mass, real-only (positive).
ComplexHypreParMatrix AssembleBtt(const FiniteElementSpace &nd_fespace,
                                  const MaterialOperator &mat_op);

// Att = mu_cc^{-1} curl-curl  -  omega^2 eps mass  -  sigma (mu^{-1} mass) + BC-t
// (impedance / absorbing / conductivity). Frequency- and shift-dependent.
ComplexHypreParMatrix
AssembleAtt(const FiniteElementSpace &nd_fespace, const MaterialOperator &mat_op,
            const mfem::Vector *normal, SurfaceImpedanceOperator &surf_z_op,
            FarfieldBoundaryOperator &farfield_op,
            SurfaceConductivityOperator &surf_sigma_op, double omega, double sigma);

// Ann = -(mu^{-1} grad u, grad v) + omega^2 (eps u, v) + BC-n. Frequency-dependent.
// farfield_op and surf_sigma_op contribute impedance / loss terms on the H1 block.
ComplexHypreParMatrix AssembleAnn(const FiniteElementSpace &h1_fespace,
                                  const MaterialOperator &mat_op,
                                  const mfem::Vector *normal,
                                  SurfaceImpedanceOperator &surf_z_op,
                                  FarfieldBoundaryOperator &farfield_op,
                                  SurfaceConductivityOperator &surf_sigma_op, double omega);

// Alias the ND and H1 halves of a pre-loaded eigenvector e0 = [e_t_tilde; e_n_tilde] as
// et / en, and apply the Vardapetyan–Demkowicz back-transform en := ẽn / (i·kn) so en
// holds the physical En. Pure scalar op, no MPI.
void ApplyVDBackTransform(ComplexVector &e0, std::complex<double> kn, int nd_size,
                          int h1_size, ComplexVector &et, ComplexVector &en);

}  // namespace mode_assembly

}  // namespace palace

#endif  // PALACE_MODELS_MODE_ASSEMBLY_HPP
