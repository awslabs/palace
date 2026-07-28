// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_BOUNDARY_MODE_OPERATOR_HPP
#define PALACE_MODELS_BOUNDARY_MODE_OPERATOR_HPP

#include <complex>
#include <memory>
#include <tuple>
#include <vector>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "fem/singularassembly.hpp"
#include "fem/singularfeatures.hpp"
#include "fem/singularsystem.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "utils/iodata.hpp"

namespace palace
{

class FarfieldBoundaryOperator;
class SurfaceConductivityOperator;
class SurfaceRationalImpedanceOperator;
class SurfaceImpedanceOperator;

//
// Top-level operator for 2D boundary mode analysis, analogous to SpaceOperator for 3D
// problems. Owns FE spaces and boundary operators. The driver is expected to fold any
// 3D-parent frame information into mesh + iodata before construction.
//
class BoundaryModeOperator
{
public:
  struct SingularCoefficientNorms
  {
    double nd_gradient = 0.0;
    double nd_rotational = 0.0;
    double h1_gradient = 0.0;
  };

  struct SingularFieldEnergies
  {
    double electric_transverse = 0.0;
    double electric_normal = 0.0;
    double magnetic_transverse = 0.0;
    double magnetic_normal = 0.0;
  };

  // Caller owns mat_op; it must outlive this operator.
  BoundaryModeOperator(
      const IoData &iodata, const std::vector<std::unique_ptr<Mesh>> &mesh,
      const MaterialOperator &mat_op,
      const fem::singular::TriangleFeatureTopology *singular_features = nullptr,
      const std::vector<fem::singular::GlobalVertexId> *source_vertex_ids = nullptr);

  // Access FE spaces.
  FiniteElementSpace &GetNDSpace() { return nd_fespaces.GetFinestFESpace(); }
  FiniteElementSpace &GetH1Space() { return h1_fespaces.GetFinestFESpace(); }
  FiniteElementSpace &GetRTSpace() { return rt_fespaces.GetFinestFESpace(); }
  FiniteElementSpace &GetCurlSpace() { return *l2_curl_fespace; }
  const FiniteElementSpace &GetNDSpace() const { return nd_fespaces.GetFinestFESpace(); }
  const FiniteElementSpace &GetH1Space() const { return h1_fespaces.GetFinestFESpace(); }
  const FiniteElementSpace &GetRTSpace() const { return rt_fespaces.GetFinestFESpace(); }
  const FiniteElementSpace &GetCurlSpace() const { return *l2_curl_fespace; }

  // Access space hierarchies (for error estimation and multigrid).
  FiniteElementSpaceHierarchy &GetNDSpaceHierarchy() { return nd_fespaces; }
  FiniteElementSpaceHierarchy &GetH1SpaceHierarchy() { return h1_fespaces; }
  FiniteElementSpaceHierarchy &GetH1AuxSpaceHierarchy() { return h1_aux_fespaces; }
  FiniteElementSpaceHierarchy &GetRTSpaceHierarchy() { return rt_fespaces; }
  const FiniteElementSpaceHierarchy &GetNDSpaceHierarchy() const { return nd_fespaces; }
  const FiniteElementSpaceHierarchy &GetH1SpaceHierarchy() const { return h1_fespaces; }
  const FiniteElementSpaceHierarchy &GetH1AuxSpaceHierarchy() const
  {
    return h1_aux_fespaces;
  }
  const FiniteElementSpaceHierarchy &GetRTSpaceHierarchy() const { return rt_fespaces; }

  // Access per-level essential BC true DOF lists for the block system.
  std::vector<mfem::Array<int>> &GetNDDbcTDofLists() { return nd_dbc_tdof_lists; }
  std::vector<mfem::Array<int>> &GetH1DbcTDofLists() { return h1_dbc_tdof_lists; }
  std::vector<mfem::Array<int>> &GetH1AuxDbcTDofLists() { return h1_aux_dbc_tdof_lists; }

  // Access material and boundary operators.
  const MaterialOperator &GetMaterialOp() const { return mat_op; }
  SurfaceImpedanceOperator &GetSurfZOp() { return *surf_z_op; }
  FarfieldBoundaryOperator &GetFarfieldOp() { return *farfield_op; }
  SurfaceConductivityOperator &GetSurfSigmaOp() { return *surf_sigma_op; }
  SurfaceRationalImpedanceOperator &GetSurfRZOp() { return *surf_rz_op; }

  // Frequency-independent block matrices, assembled in the constructor.
  //   Atn = -(mu^{-1} grad_t u, v)           (ND / H1 gradient coupling)
  //   Btn = -Atn^T                            (negative transpose coupling)
  //   Btt = (mu^{-1} u, v)                    (ND mass, positive)
  const mfem::HypreParMatrix *GetAtnr() const { return Atnr.get(); }
  const mfem::HypreParMatrix *GetAtni() const { return Atni.get(); }
  const mfem::HypreParMatrix *GetBtnr() const { return Btnr.get(); }
  const mfem::HypreParMatrix *GetBtni() const { return Btni.get(); }
  const mfem::HypreParMatrix *GetBtt() const { return Bttr.get(); }
  const mfem::HypreParMatrix *GetCombinedGradient() const
  {
    return singular_gradient.get();
  }

  using ComplexHypreParMatrix = std::tuple<std::unique_ptr<mfem::HypreParMatrix>,
                                           std::unique_ptr<mfem::HypreParMatrix>>;

  // Frequency-dependent block assembly on the finest FE spaces, rebuilt each solve.
  //   Att = mu_cc^{-1}(curl_t u, curl_t v) - omega^2 (eps u, v) - sigma (mu^{-1} u, v)
  //         + BC-t (impedance, absorbing, conductivity)
  //   Ann = -(mu^{-1} grad u, grad v) + omega^2 (eps u, v) + BC-n
  ComplexHypreParMatrix AssembleAtt(std::complex<double> omega, double sigma) const;
  ComplexHypreParMatrix AssembleAnn(std::complex<double> omega) const;

  // Alias the ND and H1 halves of a pre-loaded eigenvector e0 = [e_t_tilde; e_n_tilde]
  // as et / en, then apply the Vardapetyan–Demkowicz back-transform en := ẽn / (i·kn)
  // so en holds the physical En.
  void ApplyVDBackTransform(ComplexVector &e0, std::complex<double> kn, ComplexVector &et,
                            ComplexVector &en) const;

  // Poynting power P = (1/2) conj(kn)/omega · etᴴ Btt et + i/(2·omega) · etᴴ Atn En for
  // physical (et, En). Used for mode power normalization and impedance postprocessing.
  std::complex<double> ComputePoyntingPower(double omega, std::complex<double> kn,
                                            const ComplexVector &et,
                                            const ComplexVector &en) const;
  SingularCoefficientNorms ComputeSingularCoefficientNorms(const ComplexVector &et,
                                                           const ComplexVector &en) const;
  SingularFieldEnergies ComputeSingularFieldEnergies(double omega, std::complex<double> kn,
                                                     const ComplexVector &et,
                                                     const ComplexVector &en) const;

  const fem::singular::TriangleDofTopology &GetSingularDofTopology() const
  {
    MFEM_VERIFY(singular_dofs, "BoundaryMode operator has no singular DOF topology!");
    return *singular_dofs;
  }
  const fem::singular::ParallelDofNumbering &GetSingularDofNumbering() const
  {
    MFEM_VERIFY(singular_numbering,
                "BoundaryMode operator has no singular parallel numbering!");
    return *singular_numbering;
  }

  // Access mesh.
  Mesh &GetMesh() { return *solve_mesh; }
  const Mesh &GetMesh() const { return *solve_mesh; }
  MPI_Comm GetComm() const { return solve_mesh->GetComm(); }

  // Access solver order.
  int GetSolverOrder() const { return solver_order; }

  // True vector sizes.
  int GetNDTrueVSize() const
  {
    return nd_fespaces.GetFinestFESpace().GetTrueVSize() +
           (singular_numbering ? static_cast<int>(singular_numbering->nd.owned_size) : 0);
  }
  int GetH1TrueVSize() const
  {
    return h1_fespaces.GetFinestFESpace().GetTrueVSize() +
           (singular_numbering ? static_cast<int>(singular_numbering->h1.owned_size) : 0);
  }
  HYPRE_BigInt GetNDGlobalTrueVSize() const
  {
    return nd_fespaces.GetFinestFESpace().GlobalTrueVSize() +
           (singular_numbering ? singular_numbering->nd.global_size : 0);
  }
  HYPRE_BigInt GetH1GlobalTrueVSize() const
  {
    return h1_fespaces.GetFinestFESpace().GlobalTrueVSize() +
           (singular_numbering ? singular_numbering->h1.global_size : 0);
  }
  bool HasSingularEnrichment() const { return singular_features != nullptr; }
  const mfem::Array<int> &GetCombinedDbcTDofList() const { return combined_dbc_tdof_list; }

private:
  const IoData &iodata;
  int solver_order;

  // Non-owning pointer to the caller-owned solve mesh.
  Mesh *solve_mesh;

  // FE collections and space hierarchies.
  std::vector<std::unique_ptr<mfem::ND_FECollection>> nd_fecs;
  std::vector<std::unique_ptr<mfem::H1_FECollection>> h1_fecs;
  std::vector<std::unique_ptr<mfem::H1_FECollection>> h1_aux_fecs;
  std::vector<std::unique_ptr<mfem::RT_FECollection>> rt_fecs;
  FiniteElementSpaceHierarchy nd_fespaces, h1_fespaces, h1_aux_fespaces, rt_fespaces;
  std::vector<mfem::Array<int>> nd_dbc_tdof_lists, h1_dbc_tdof_lists, h1_aux_dbc_tdof_lists;

  // L2 curl space for B-field in 2D.
  std::unique_ptr<mfem::L2_FECollection> l2_curl_fec;
  std::unique_ptr<FiniteElementSpace> l2_curl_fespace;

  // Optional triangular singular enrichment. Combined true vectors use
  // process-local ordering [standard, owned enrichment] independently in the
  // ND and H1 blocks.
  const fem::singular::TriangleFeatureTopology *singular_features;
  const std::vector<fem::singular::GlobalVertexId> *source_vertex_ids;
  std::unique_ptr<fem::singular::TriangleDofTopology> singular_dofs;
  std::unique_ptr<fem::singular::ParallelDofNumbering> singular_numbering;
  fem::singular::ParallelSparseEnrichmentMatrices singular_mu_matrices;
  fem::singular::ParallelSparseEnrichmentMatrices singular_epsilon_matrices;
  fem::singular::ParallelSparseEnrichmentMatrices singular_epsilon_imag_matrices;
  std::unique_ptr<mfem::HypreParMatrix> singular_gradient;
  std::unique_ptr<mfem::HypreParMatrix> singular_epsilon_nd_mass;
  std::unique_ptr<mfem::HypreParMatrix> singular_epsilon_h1_mass;
  std::unique_ptr<mfem::HypreParMatrix> singular_mu_nd_curl_curl;
  mfem::Array<int> singular_nd_essential_true_dofs;
  mfem::Array<int> singular_h1_essential_true_dofs;
  mfem::Array<int> combined_dbc_tdof_list;

  // Material and boundary operators. mat_op is caller-owned; the rest are constructed
  // from mat_op + the solve mesh at construction time.
  const MaterialOperator &mat_op;
  std::unique_ptr<SurfaceImpedanceOperator> surf_z_op;
  std::unique_ptr<FarfieldBoundaryOperator> farfield_op;
  std::unique_ptr<SurfaceConductivityOperator> surf_sigma_op;
  std::unique_ptr<SurfaceRationalImpedanceOperator> surf_rz_op;

  // Frequency-independent block matrices (assembled in the constructor).
  std::unique_ptr<mfem::HypreParMatrix> Atnr, Atni, Btnr, Btni, Bttr;

  // DBC attributes.
  mfem::Array<int> dbc_bcs;

  // Setup helper.
  void SetUpFESpaces(const std::vector<std::unique_ptr<Mesh>> &mesh);
  void SetUpSingularEnrichment();
};

}  // namespace palace

#endif  // PALACE_MODELS_BOUNDARY_MODE_OPERATOR_HPP
