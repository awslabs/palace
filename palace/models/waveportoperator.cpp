// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "waveportoperator.hpp"
#include "fem/bilinearform.hpp"
#include "linalg/amg.hpp"
#include "linalg/ams.hpp"
#include "linalg/arpack.hpp"
#include "linalg/blockprecond.hpp"
#include "linalg/gmg.hpp"
#include "linalg/hypre.hpp"
#include "linalg/iterative.hpp"
#include "linalg/mumps.hpp"
#include "linalg/rap.hpp"
#include "linalg/solver.hpp"

#include <fmt/ranges.h>
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "fem/interpolator.hpp"
#include "linalg/amg.hpp"
#include "linalg/ams.hpp"
#include "linalg/arpack.hpp"
#include "linalg/blockprecond.hpp"
#include "linalg/gmg.hpp"
#include "linalg/iterative.hpp"
#include "linalg/mumps.hpp"
#include "linalg/slepc.hpp"
#include "linalg/strumpack.hpp"
#include "linalg/superlu.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/surfaceconductivityoperator.hpp"
#include "models/surfaceimpedanceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

namespace
{

void GetEssentialTrueDofs(mfem::ParGridFunction &E0t, mfem::ParGridFunction &E0n,
                          mfem::ParGridFunction &port_E0t, mfem::ParGridFunction &port_E0n,
                          mfem::ParTransferMap &port_nd_transfer,
                          mfem::ParTransferMap &port_h1_transfer,
                          const mfem::Array<int> &dbc_attr,
                          mfem::Array<int> &port_nd_dbc_tdof_list,
                          mfem::Array<int> &port_h1_dbc_tdof_list)
{
  auto &nd_fespace = *E0t.ParFESpace();
  auto &h1_fespace = *E0n.ParFESpace();
  auto &port_nd_fespace = *port_E0t.ParFESpace();
  auto &port_h1_fespace = *port_E0n.ParFESpace();
  const auto &mesh = *nd_fespace.GetParMesh();

  mfem::Array<int> dbc_marker, nd_dbc_tdof_list, h1_dbc_tdof_list;
  mesh::AttrToMarker(mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0, dbc_attr,
                     dbc_marker);
  nd_fespace.GetEssentialTrueDofs(dbc_marker, nd_dbc_tdof_list);
  h1_fespace.GetEssentialTrueDofs(dbc_marker, h1_dbc_tdof_list);

  Vector tE0t(nd_fespace.GetTrueVSize()), tE0n(h1_fespace.GetTrueVSize());
  tE0t.UseDevice(true);
  tE0n.UseDevice(true);
  tE0t = 0.0;
  tE0n = 0.0;
  linalg::SetSubVector(tE0t, nd_dbc_tdof_list, 1.0);
  linalg::SetSubVector(tE0n, h1_dbc_tdof_list, 1.0);
  E0t.SetFromTrueDofs(tE0t);
  E0n.SetFromTrueDofs(tE0n);
  port_nd_transfer.Transfer(E0t, port_E0t);
  port_h1_transfer.Transfer(E0n, port_E0n);

  Vector port_tE0t(port_nd_fespace.GetTrueVSize()),
      port_tE0n(port_h1_fespace.GetTrueVSize());
  port_tE0t.UseDevice(true);
  port_tE0n.UseDevice(true);
  port_E0t.ParallelProject(port_tE0t);
  port_E0n.ParallelProject(port_tE0n);
  {
    const auto *h_port_tE0t = port_tE0t.HostRead();
    const auto *h_port_tE0n = port_tE0n.HostRead();
    for (int i = 0; i < port_tE0t.Size(); i++)
    {
      if (h_port_tE0t[i] != 0.0)
      {
        port_nd_dbc_tdof_list.Append(i);
      }
    }
    for (int i = 0; i < port_tE0n.Size(); i++)
    {
      if (h_port_tE0n[i] != 0.0)
      {
        port_h1_dbc_tdof_list.Append(i);
      }
    }
  }
}

void GetInitialSpace(const mfem::ParFiniteElementSpace &nd_fespace,
                     const mfem::ParFiniteElementSpace &h1_fespace,
                     const mfem::Array<int> &dbc_tdof_list, ComplexVector &v)
{
  // Initial space which satisfies Dirichlet BCs.
  const int nd_size = nd_fespace.GetTrueVSize(), h1_size = h1_fespace.GetTrueVSize();
  v.SetSize(nd_size + h1_size);
  v.UseDevice(true);
  v = std::complex<double>(1.0, 0.0);
  // linalg::SetRandomReal(nd_fespace.GetComm(), v);
  linalg::SetSubVector(v, nd_size, nd_size + h1_size, 0.0);
  linalg::SetSubVector(v, dbc_tdof_list, 0.0);
}

void Normalize(const GridFunction &S0t, GridFunction &E0t, GridFunction &E0n,
               mfem::LinearForm &sr, mfem::LinearForm &si)
{
  // Normalize grid functions to a chosen polarization direction and unit power, |E x H⋆| ⋅
  // n, integrated over the port surface (+n is the outward mesh normal). The n x H
  // coefficients are updated implicitly as the only store references to the Et, En grid
  // functions. We choose a (rather arbitrary) phase constraint to at least make results for
  // the same port consistent between frequencies/meshes.

  // |E x H⋆| ⋅ n = |E ⋅ (-n x H⋆)|. This also updates the n x H coefficients depending on
  // Et, En. Update linear forms for postprocessing too.
  std::complex<double> dot[2] = {
      {sr * S0t.Real(), si * S0t.Real()},
      {-(sr * E0t.Real()) - (si * E0t.Imag()), -(sr * E0t.Imag()) + (si * E0t.Real())}};
  Mpi::GlobalSum(2, dot, S0t.ParFESpace()->GetComm());
  auto scale = std::abs(dot[0]) / (dot[0] * std::sqrt(std::abs(dot[1])));
  ComplexVector::AXPBY(scale, E0t.Real(), E0t.Imag(), 0.0, E0t.Real(), E0t.Imag());
  ComplexVector::AXPBY(scale, E0n.Real(), E0n.Imag(), 0.0, E0n.Real(), E0n.Imag());
  ComplexVector::AXPBY(scale, sr, si, 0.0, sr, si);

  // This parallel communication is not required since wave port boundaries are true one-
  // sided boundaries.
  // E0t.Real().ExchangeFaceNbrData();  // Ready for parallel comm on shared faces for n x H
  // E0t.Imag().ExchangeFaceNbrData();  // coefficients evaluation
  // E0n.Real().ExchangeFaceNbrData();
  // E0n.Imag().ExchangeFaceNbrData();
}

// Helper for BdrSubmeshEVectorCoefficient and BdrSubmeshHVectorCoefficient.
enum class ValueType
{
  REAL,
  IMAG
};

// Return as a vector coefficient the boundary mode electric field.
template <ValueType Type>
class BdrSubmeshEVectorCoefficient : public mfem::VectorCoefficient
{
private:
  const GridFunction &Et, &En;
  const mfem::ParSubMesh &submesh;
  const std::unordered_map<int, int> &submesh_parent_elems;
  mfem::IsoparametricTransformation T_loc;
  const double scaling;

public:
  BdrSubmeshEVectorCoefficient(const GridFunction &Et, const GridFunction &En,
                               const mfem::ParSubMesh &submesh,
                               const std::unordered_map<int, int> &submesh_parent_elems,
                               double scaling = 1.0)
    : mfem::VectorCoefficient(Et.Real().VectorDim()), Et(Et), En(En), submesh(submesh),
      submesh_parent_elems(submesh_parent_elems), scaling(scaling)
  {
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    // Always do the GridFunction evaluation in the submesh.
    mfem::ElementTransformation *T_submesh = nullptr;
    if (T.mesh == submesh.GetParent())
    {
      MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
                  "BdrSubmeshEVectorCoefficient requires ElementType::BDR_ELEMENT when not "
                  "used on a SubMesh!");
      auto it = submesh_parent_elems.find(T.ElementNo);
      if (it == submesh_parent_elems.end())
      {
        // Just return zero for a parent boundary element not in the submesh.
        V.SetSize(vdim);
        V = 0.0;
        return;
      }
      else
      {
        submesh.GetElementTransformation(it->second, &T_loc);
        T_loc.SetIntPoint(&ip);
        T_submesh = &T_loc;
      }
    }
    else if (T.mesh == &submesh)
    {
      MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::ELEMENT,
                  "BdrSubmeshEVectorCoefficient requires ElementType::ELEMENT when used on "
                  "a SubMesh!");
      T_submesh = &T;
    }
    else
    {
      MFEM_ABORT("Invalid mesh for BdrSubmeshEVectorCoefficient!");
    }

    // Compute Eₜ + n ⋅ Eₙ . The normal returned by GetNormal points out of the
    // computational domain, so we reverse it (direction of propagation is into the domain).
    double normal_data[3];
    mfem::Vector normal(normal_data, vdim);
    BdrGridFunctionCoefficient::GetNormal(*T_submesh, normal);
    if constexpr (Type == ValueType::REAL)
    {
      Et.Real().GetVectorValue(*T_submesh, ip, V);
      auto Vn = En.Real().GetValue(*T_submesh, ip);
      V.Add(-Vn, normal);
    }
    else
    {
      Et.Imag().GetVectorValue(*T_submesh, ip, V);
      auto Vn = En.Imag().GetValue(*T_submesh, ip);
      V.Add(-Vn, normal);
    }
    V *= scaling;
  }
};

// Computes boundary mode n x H, where +n is the outward mesh normal: n x H =
// -1/(iωμ) (ikₙ Eₜ + ∇ₜ Eₙ), using the tangential and normal electric field component grid
// functions evaluated on the (single-sided) boundary element.
template <ValueType Type>
class BdrSubmeshHVectorCoefficient : public mfem::VectorCoefficient
{
private:
  const GridFunction &Et, &En;
  const MaterialOperator &mat_op;
  const mfem::ParSubMesh &submesh;
  const std::unordered_map<int, int> &submesh_parent_elems;
  mfem::IsoparametricTransformation T_loc;
  std::complex<double> kn;
  double omega;

public:
  BdrSubmeshHVectorCoefficient(const GridFunction &Et, const GridFunction &En,
                               const MaterialOperator &mat_op,
                               const mfem::ParSubMesh &submesh,
                               const std::unordered_map<int, int> &submesh_parent_elems,
                               std::complex<double> kn, double omega)
    : mfem::VectorCoefficient(Et.Real().VectorDim()), Et(Et), En(En), mat_op(mat_op),
      submesh(submesh), submesh_parent_elems(submesh_parent_elems), kn(kn), omega(omega)
  {
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    // Always do the GridFunction evaluation in the submesh.
    mfem::ElementTransformation *T_submesh = nullptr;
    if (T.mesh == submesh.GetParent())
    {
      MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
                  "BdrSubmeshHVectorCoefficient requires ElementType::BDR_ELEMENT when not "
                  "used on a SubMesh!");
      auto it = submesh_parent_elems.find(T.ElementNo);
      if (it == submesh_parent_elems.end())
      {
        // Just return zero for a parent boundary element not in the submesh.
        V.SetSize(vdim);
        V = 0.0;
        return;
      }
      else
      {
        submesh.GetElementTransformation(it->second, &T_loc);
        T_loc.SetIntPoint(&ip);
        T_submesh = &T_loc;
      }
    }
    else if (T.mesh == &submesh)
    {
      MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::ELEMENT,
                  "BdrSubmeshHVectorCoefficient requires ElementType::ELEMENT when used on "
                  "a SubMesh!");
      T_submesh = &T;
    }
    else
    {
      MFEM_ABORT("Invalid mesh for BdrSubmeshHVectorCoefficient!");
    }

    // Get the attribute in the neighboring domain element of the parent mesh.
    int attr = [&T, this]()
    {
      int i = -1, o, iel1, iel2;
      if (T.mesh == submesh.GetParent())
      {
        MFEM_ASSERT(
            T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
            "BdrSubmeshHVectorCoefficient requires ElementType::BDR_ELEMENT when not "
            "used on a SubMesh!");
        T.mesh->GetBdrElementFace(T.ElementNo, &i, &o);
      }
      else if (T.mesh == &submesh)
      {
        MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::ELEMENT,
                    "BdrSubmeshHVectorCoefficient requires ElementType::ELEMENT when used "
                    "on a SubMesh!");
        submesh.GetParent()->GetBdrElementFace(submesh.GetParentElementIDMap()[T.ElementNo],
                                               &i, &o);
      }
      else
      {
        MFEM_ABORT("Invalid mesh for BdrSubmeshHVectorCoefficient!");
      }
      submesh.GetParent()->GetFaceElements(i, &iel1, &iel2);
      return submesh.GetParent()->GetAttribute(iel1);
    }();

    // Compute Re/Im{-1/i (ikₙ Eₜ + ∇ₜ Eₙ)} (t-gradient evaluated in boundary element).
    double U_data[3];
    mfem::Vector U(U_data, vdim);
    if constexpr (Type == ValueType::REAL)
    {
      Et.Real().GetVectorValue(*T_submesh, ip, U);
      U *= -kn.real();

      double dU_data[3];
      mfem::Vector dU(dU_data, vdim);
      En.Imag().GetGradient(*T_submesh, dU);
      U -= dU;
    }
    else
    {
      Et.Imag().GetVectorValue(*T_submesh, ip, U);
      U *= -kn.real();

      double dU_data[3];
      mfem::Vector dU(dU_data, vdim);
      En.Real().GetGradient(*T_submesh, dU);
      U += dU;
    }

    // Scale by 1/(ωμ) with μ evaluated in the neighboring element.
    V.SetSize(U.Size());
    mat_op.GetInvPermeability(attr).Mult(U, V);
    V *= (1.0 / omega);
  }
};

}  // namespace

WavePortData::WavePortData(const config::WavePortData &data,
                           const config::BoundaryData &boundaries,
                           const config::DomainData &domains, ProblemType problem_type,
                           const config::LinearSolverData &linear, const Units &units,
                           const MaterialOperator &mat_op,
                           mfem::ParFiniteElementSpace &nd_fespace,
                           mfem::ParFiniteElementSpace &h1_fespace,
                           const mfem::Array<int> &dbc_attr)
  : mat_op(mat_op), excitation(data.excitation), active(data.active)
{
  mode_idx = data.mode_idx;
  d_offset = data.d_offset;
  kn0 = 0.0;
  omega0 = 0.0;

  // Construct the SubMesh.
  MFEM_VERIFY(!data.attributes.empty(), "Wave port boundary found with no attributes!");
  const auto &mesh = *nd_fespace.GetParMesh();
  attr_list.Append(data.attributes.data(), data.attributes.size());
  auto port_submesh_ptr = std::make_unique<mfem::ParSubMesh>(
      mfem::ParSubMesh::CreateFromBoundary(mesh, attr_list));

  // Add internal boundary elements for edges where the port face intersects other boundary
  // faces (PEC, impedance, conductivity, absorbing). ParSubMesh::CreateFromBoundary only
  // creates boundary elements at the geometric boundary of the selected face region, but
  // internal intersections need boundary elements for the 2D eigenvalue problem BCs.
  {
    std::vector<int> internal_bdr_attrs;
    for (auto a : boundaries.pec.attributes)
    {
      internal_bdr_attrs.push_back(a);
    }
    for (auto a : boundaries.auxpec.attributes)
    {
      internal_bdr_attrs.push_back(a);
    }
    for (const auto &d : boundaries.impedance)
    {
      for (auto a : d.attributes)
      {
        internal_bdr_attrs.push_back(a);
      }
    }
    for (const auto &d : boundaries.conductivity)
    {
      for (auto a : d.attributes)
      {
        internal_bdr_attrs.push_back(a);
      }
    }
    for (auto a : boundaries.farfield.attributes)
    {
      internal_bdr_attrs.push_back(a);
    }
    mesh::AddSubMeshInternalBoundaryElements(*port_submesh_ptr, attr_list,
                                             internal_bdr_attrs);
  }

  port_mesh = std::make_unique<Mesh>(std::move(port_submesh_ptr));
  port_normal = mesh::GetSurfaceNormal(*port_mesh);

  port_nd_fec = std::make_unique<mfem::ND_FECollection>(nd_fespace.GetMaxElementOrder(),
                                                        port_mesh->Dimension());
  port_h1_fec = std::make_unique<mfem::H1_FECollection>(h1_fespace.GetMaxElementOrder(),
                                                        port_mesh->Dimension());
  port_nd_fespace = std::make_unique<FiniteElementSpace>(*port_mesh, port_nd_fec.get());
  port_h1_fespace = std::make_unique<FiniteElementSpace>(*port_mesh, port_h1_fec.get());

  GridFunction E0t(nd_fespace), E0n(h1_fespace);
  port_E0t = std::make_unique<GridFunction>(*port_nd_fespace, true);
  port_E0n = std::make_unique<GridFunction>(*port_h1_fespace, true);
  port_E = std::make_unique<GridFunction>(*port_nd_fespace, true);

  port_nd_transfer = std::make_unique<mfem::ParTransferMap>(
      mfem::ParSubMesh::CreateTransferMap(E0t.Real(), port_E0t->Real()));
  port_h1_transfer = std::make_unique<mfem::ParTransferMap>(
      mfem::ParSubMesh::CreateTransferMap(E0n.Real(), port_E0n->Real()));

  // Remap submesh attributes so that domain elements get the adjacent volume element
  // attribute (matching material definitions) and boundary edges get the adjacent boundary
  // face attribute (matching BC definitions). This must happen AFTER CreateTransferMap
  // (which depends on the original parent-child attribute relationship). Then rebuild the
  // CEED attribute maps from the remapped submesh.
  {
    auto &port_submesh = static_cast<mfem::ParSubMesh &>(port_mesh->Get());
    mesh::RemapSubMeshAttributes(port_submesh);
    mesh::RemapSubMeshBdrAttributes(port_submesh, attr_list);
  }
  port_mesh->RebuildCeedAttributes();

  // Construct submesh-specific MaterialOperator (uses remapped submesh for CEED attribute
  // mapping) and boundary condition operators. The surface operators are constructed with
  // the parent 3D mesh for bdr_attributes validation (modifying a ParSubMesh's
  // bdr_attributes corrupts MFEM internal state), but use port_mat_op for CEED boundary
  // attribute lookups which correctly reference the remapped submesh.
  port_mat_op = std::make_unique<MaterialOperator>(domains.materials, boundaries.periodic,
                                                   problem_type, *port_mesh);
  port_surf_z_op = std::make_unique<SurfaceImpedanceOperator>(
      boundaries.impedance, boundaries.cracked_attributes, units, *port_mat_op, mesh);
  port_farfield_op = std::make_unique<FarfieldBoundaryOperator>(
      boundaries.farfield, problem_type, *port_mat_op, mesh);
  port_surf_sigma_op = std::make_unique<SurfaceConductivityOperator>(
      boundaries.conductivity, problem_type, units, *port_mat_op, mesh);

  // Construct mapping from parent (boundary) element indices to submesh (domain)
  // elements.
  {
    const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
    const mfem::Array<int> &parent_elems = port_submesh.GetParentElementIDMap();
    for (int i = 0; i < parent_elems.Size(); i++)
    {
      submesh_parent_elems[parent_elems[i]] = i;
    }
  }

  // Extract Dirichlet BC true dofs for the port FE spaces.
  {
    mfem::Array<int> port_nd_dbc_tdof_list, port_h1_dbc_tdof_list;
    GetEssentialTrueDofs(E0t.Real(), E0n.Real(), port_E0t->Real(), port_E0n->Real(),
                         *port_nd_transfer, *port_h1_transfer, dbc_attr,
                         port_nd_dbc_tdof_list, port_h1_dbc_tdof_list);
    int nd_tdof_offset = port_nd_fespace->GetTrueVSize();
    port_dbc_tdof_list.Reserve(port_nd_dbc_tdof_list.Size() + port_h1_dbc_tdof_list.Size());
    for (auto tdof : port_nd_dbc_tdof_list)
    {
      port_dbc_tdof_list.Append(tdof);
    }
    for (auto tdof : port_h1_dbc_tdof_list)
    {
      port_dbc_tdof_list.Append(tdof + nd_tdof_offset);
    }
  }

  // Create vector for initial space for eigenvalue solves and eigenmode solution.
  GetInitialSpace(*port_nd_fespace, *port_h1_fespace, port_dbc_tdof_list, v0);
  e0.SetSize(port_nd_fespace->GetTrueVSize() + port_h1_fespace->GetTrueVSize());
  e0.UseDevice(true);

  // Set the shift-and-invert target as the maximum propagation constant.
  double c_min = mat_op.GetLightSpeedMax().Min();
  Mpi::GlobalMin(1, &c_min, nd_fespace.GetComm());
  MFEM_VERIFY(c_min > 0.0 && c_min < mfem::infinity(),
              "Invalid material speed of light detected in WavePortOperator!");
  mu_eps_max = 1.0 / (c_min * c_min) * 1.1;  // Add a safety factor for maximum
                                             // propagation constant possible

  // Configure a communicator for the processes which have elements for this port.
  MPI_Comm comm = nd_fespace.GetComm();
  int color = (port_nd_fespace->GetVSize() > 0 || port_h1_fespace->GetVSize() > 0)
                  ? 0
                  : MPI_UNDEFINED;
  MPI_Comm_split(comm, color, Mpi::Rank(comm), &port_comm);
  MFEM_VERIFY((color == 0 && port_comm != MPI_COMM_NULL) ||
                  (color == MPI_UNDEFINED && port_comm == MPI_COMM_NULL),
              "Unexpected error splitting communicator for wave port boundaries!");
  port_root = (color == MPI_UNDEFINED) ? Mpi::Size(comm) : Mpi::Rank(comm);
  Mpi::GlobalMin(1, &port_root, comm);
  MFEM_VERIFY(port_root < Mpi::Size(comm), "No root process found for port!");

  // Configure the boundary mode solver. Matrix assembly is MPI-collective on the FE space
  // communicator (all processes), so the config + construction must happen on all
  // processes. The solver_comm (port_comm) restricts solver setup to port processes only.
  {
    mode_solver = std::make_unique<ModeEigenSolver>(
        *port_mat_op, &port_normal, port_surf_z_op.get(), port_farfield_op.get(),
        port_surf_sigma_op.get(), *port_nd_fespace, *port_h1_fespace, port_dbc_tdof_list,
        mode_idx, data.max_size, data.eig_tol, EigenvalueSolver::WhichType::LARGEST_REAL,
        linear, data.eigen_solver, data.verbose, port_comm);
  }

  // Configure port mode sign convention: 1ᵀ Re{-n x H} >= 0 on the "upper-right quadrant"
  // of the wave port boundary, in order to deal with symmetry effectively.
  {
    Vector bbmin, bbmax;
    mesh::GetAxisAlignedBoundingBox(*port_mesh, bbmin, bbmax);
    const int dim = port_mesh->SpaceDimension();

    double la = 0.0, lb = 0.0;
    int da = -1, db = -1;
    for (int d = 0; d < dim; d++)
    {
      double diff = bbmax(d) - bbmin(d);
      if (diff > la)
      {
        lb = la;
        la = diff;
        db = da;
        da = d;
      }
      else if (diff > lb)
      {
        lb = diff;
        db = d;
      }
    }
    MFEM_VERIFY(da >= 0 && db >= 0 && da != db,
                "Unexpected wave port geometry for normalization!");
    double ca = 0.5 * (bbmax[da] + bbmin[da]), cb = 0.5 * (bbmax[db] + bbmin[db]);

    auto TDirection = [da, db, ca, cb, dim](const Vector &x, Vector &f)
    {
      MFEM_ASSERT(x.Size() == dim,
                  "Invalid dimension mismatch for wave port mode normalization!");
      f.SetSize(dim);
      if (x[da] >= ca && x[db] >= cb)
      {
        f = 1.0;
      }
      else
      {
        f = 0.0;
      }
    };
    mfem::VectorFunctionCoefficient tfunc(dim, TDirection);
    port_S0t = std::make_unique<GridFunction>(*port_nd_fespace);
    port_S0t->Real().ProjectCoefficient(tfunc);
  }

  // Store voltage path coordinates if provided.
  // Coordinates are already nondimensionalized in IoData::NondimensionalizeInputs.
  if (data.voltage_path.size() >= 2)
  {
    has_voltage_coords = true;
    for (const auto &pt : data.voltage_path)
    {
      mfem::Vector p(pt.size());
      for (int d = 0; d < static_cast<int>(pt.size()); d++)
      {
        p(d) = pt[d];
      }
      voltage_path.push_back(std::move(p));
    }
    voltage_n_samples = data.n_samples;
  }
}

WavePortData::WavePortData(
    const MaterialOperator &mat_op_ref, const mfem::Vector *normal_,
    SurfaceImpedanceOperator *surf_z_op_, FarfieldBoundaryOperator *farfield_op_,
    SurfaceConductivityOperator *surf_sigma_op_, FiniteElementSpace &nd_fespace,
    FiniteElementSpace &h1_fespace, const mfem::Array<int> &dbc_tdof_list, int num_modes_,
    int num_vec_, double eig_tol_, EigenvalueSolver::WhichType which_eig_,
    const config::LinearSolverData &linear_, EigenSolverBackend eigen_backend_,
    int verbose_, MPI_Comm solver_comm, const ModeEigenSolverMultigridConfig *mg_config)
  : mat_op(mat_op_ref), excitation(0), active(true), is_2d_direct(true)
{
  mode_idx = 0;
  d_offset = 0.0;
  kn0 = 0.0;
  omega0 = 0.0;

  mode_solver = std::make_unique<ModeEigenSolver>(
      mat_op_ref, normal_, surf_z_op_, farfield_op_, surf_sigma_op_, nd_fespace, h1_fespace,
      dbc_tdof_list, num_modes_, num_vec_, eig_tol_, which_eig_, linear_, eigen_backend_,
      verbose_, solver_comm, mg_config);
}

WavePortData::~WavePortData()
{
  // Free the solver before the communicator on which it is based.
  mode_solver.reset();
  if (!is_2d_direct && port_comm != MPI_COMM_NULL)
  {
    MPI_Comm_free(&port_comm);
  }
}

void WavePortData::Initialize(double omega)
{
  if (omega == omega0)
  {
    return;
  }

  // Solve the generalized eigenvalue problem for the desired wave port mode using the
  // ModeEigenSolver. Frequency-independent matrices were assembled in the constructor.
  const double sigma = -omega * omega * mu_eps_max;
  std::complex<double> lambda;
  {
    bool has_solver = (port_comm != MPI_COMM_NULL);
    auto result = mode_solver->Solve(omega, sigma, has_solver, has_solver ? &v0 : nullptr);
    if (has_solver)
    {
      MFEM_VERIFY(result.num_converged >= mode_idx,
                  "Wave port eigensolver did not converge!");
      lambda = mode_solver->GetEigenvalue(mode_idx - 1);
    }
  }
  Mpi::Broadcast(1, &lambda, port_root, port_mesh->GetComm());

  // Extract the eigenmode solution and postprocess. The extracted eigenvalue is λ =
  // 1 / (-kₙ² - σ).
  kn0 = std::sqrt(-sigma - 1.0 / lambda);
  omega0 = omega;

  // Separate the computed field out into eₜ and eₙ and and transform back to true
  // electric field variables: Eₜ = eₜ and Eₙ = eₙ / ikₙ.
  {
    if (port_comm != MPI_COMM_NULL)
    {
      mode_solver->GetEigenvector(mode_idx - 1, e0);
      linalg::NormalizePhase(port_comm, e0);
    }
    else
    {
      MFEM_ASSERT(e0.Size() == 0,
                  "Unexpected non-empty port FE space in wave port boundary mode solve!");
    }
    e0.Real().Read();  // Ensure memory is allocated on device before aliasing
    e0.Imag().Read();
    Vector e0tr(e0.Real(), 0, port_nd_fespace->GetTrueVSize());
    Vector e0nr(e0.Real(), port_nd_fespace->GetTrueVSize(),
                port_h1_fespace->GetTrueVSize());
    Vector e0ti(e0.Imag(), 0, port_nd_fespace->GetTrueVSize());
    Vector e0ni(e0.Imag(), port_nd_fespace->GetTrueVSize(),
                port_h1_fespace->GetTrueVSize());
    e0tr.UseDevice(true);
    e0nr.UseDevice(true);
    e0ti.UseDevice(true);
    e0ni.UseDevice(true);
    ComplexVector::AXPBY(1.0 / (1i * kn0), e0nr, e0ni, 0.0, e0nr, e0ni);
    port_E0t->Real().SetFromTrueDofs(e0tr);  // Parallel distribute
    port_E0t->Imag().SetFromTrueDofs(e0ti);
    port_E0n->Real().SetFromTrueDofs(e0nr);
    port_E0n->Imag().SetFromTrueDofs(e0ni);
  }

  // Configure the linear forms for computing S-parameters (projection of the field onto the
  // port mode). Normalize the mode for a chosen polarization direction and unit power,
  // |E x H⋆| ⋅ n, integrated over the port surface (+n is the outward mesh normal).
  {
    const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
    BdrSubmeshHVectorCoefficient<ValueType::REAL> port_nxH0r_func(
        *port_E0t, *port_E0n, mat_op, port_submesh, submesh_parent_elems, kn0, omega0);
    BdrSubmeshHVectorCoefficient<ValueType::IMAG> port_nxH0i_func(
        *port_E0t, *port_E0n, mat_op, port_submesh, submesh_parent_elems, kn0, omega0);
    {
      port_sr = std::make_unique<mfem::LinearForm>(&port_nd_fespace->Get());
      port_sr->AddDomainIntegrator(new VectorFEDomainLFIntegrator(port_nxH0r_func));
      port_sr->UseFastAssembly(false);
      port_sr->UseDevice(false);
      port_sr->Assemble();
      port_sr->UseDevice(true);
    }
    {
      port_si = std::make_unique<mfem::LinearForm>(&port_nd_fespace->Get());
      port_si->AddDomainIntegrator(new VectorFEDomainLFIntegrator(port_nxH0i_func));
      port_si->UseFastAssembly(false);
      port_si->UseDevice(false);
      port_si->Assemble();
      port_si->UseDevice(true);
    }
    Normalize(*port_S0t, *port_E0t, *port_E0n, *port_sr, *port_si);
  }
}

std::unique_ptr<mfem::VectorCoefficient>
WavePortData::GetModeExcitationCoefficientReal() const
{
  const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
  return std::make_unique<
      RestrictedVectorCoefficient<BdrSubmeshHVectorCoefficient<ValueType::REAL>>>(
      attr_list, *port_E0t, *port_E0n, mat_op, port_submesh, submesh_parent_elems, kn0,
      omega0);
}

std::unique_ptr<mfem::VectorCoefficient>
WavePortData::GetModeExcitationCoefficientImag() const
{
  const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
  return std::make_unique<
      RestrictedVectorCoefficient<BdrSubmeshHVectorCoefficient<ValueType::IMAG>>>(
      attr_list, *port_E0t, *port_E0n, mat_op, port_submesh, submesh_parent_elems, kn0,
      omega0);
}

std::unique_ptr<mfem::VectorCoefficient>
WavePortData::GetModeFieldCoefficientReal(double scaling) const
{
  const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
  return std::make_unique<
      RestrictedVectorCoefficient<BdrSubmeshEVectorCoefficient<ValueType::REAL>>>(
      attr_list, *port_E0t, *port_E0n, port_submesh, submesh_parent_elems, scaling);
}

std::unique_ptr<mfem::VectorCoefficient>
WavePortData::GetModeFieldCoefficientImag(double scaling) const
{
  const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
  return std::make_unique<
      RestrictedVectorCoefficient<BdrSubmeshEVectorCoefficient<ValueType::IMAG>>>(
      attr_list, *port_E0t, *port_E0n, port_submesh, submesh_parent_elems, scaling);
}

double WavePortData::GetExcitationPower() const
{
  // The computed port modes are normalized such that the power integrated over the port is
  // 1: ∫ (E_inc x H_inc⋆) ⋅ n dS = 1.
  return HasExcitation() ? 1.0 : 0.0;
}

std::complex<double> WavePortData::GetPower(GridFunction &E, GridFunction &B) const
{
  // Compute port power, (E x H) ⋅ n = E ⋅ (-n x H), integrated over the port surface using
  // the computed E and H = μ⁻¹ B fields, where +n is the outward mesh normal. The
  // BdrSurfaceCurrentVectorCoefficient computes -n x H for the outward normal. The linear
  // form is reconstructed from scratch each time due to changing H.
  MFEM_VERIFY(E.HasImag() && B.HasImag(),
              "Wave ports expect complex-valued E and B fields in port power "
              "calculation!");
  auto &nd_fespace = *E.ParFESpace();
  const auto &mesh = *nd_fespace.GetParMesh();
  BdrSurfaceCurrentVectorCoefficient nxHr_func(B.Real(), mat_op);
  BdrSurfaceCurrentVectorCoefficient nxHi_func(B.Imag(), mat_op);
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker = mesh::AttrToMarker(bdr_attr_max, attr_list);
  std::complex<double> dot;
  {
    mfem::LinearForm pr(&nd_fespace);
    pr.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(nxHr_func), attr_marker);
    pr.UseFastAssembly(false);
    pr.UseDevice(false);
    pr.Assemble();
    pr.UseDevice(true);
    dot = -(pr * E.Real()) - 1i * (pr * E.Imag());
  }
  {
    mfem::LinearForm pi(&nd_fespace);
    pi.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(nxHi_func), attr_marker);
    pi.UseFastAssembly(false);
    pi.UseDevice(false);
    pi.Assemble();
    pi.UseDevice(true);
    dot += -(pi * E.Imag()) + 1i * (pi * E.Real());
  }
  Mpi::GlobalSum(1, &dot, nd_fespace.GetComm());
  return dot;
}

std::complex<double> WavePortData::GetSParameter(GridFunction &E) const
{
  // Compute port S-parameter, or the projection of the field onto the port mode:
  // (E x H_inc⋆) ⋅ n = E ⋅ (-n x H_inc⋆), integrated over the port surface.
  MFEM_VERIFY(E.HasImag(),
              "Wave ports expect complex-valued E and B fields in port S-parameter "
              "calculation!");
  port_nd_transfer->Transfer(E.Real(), port_E->Real());
  port_nd_transfer->Transfer(E.Imag(), port_E->Imag());
  std::complex<double> dot(-((*port_sr) * port_E->Real()) - ((*port_si) * port_E->Imag()),
                           -((*port_sr) * port_E->Imag()) + ((*port_si) * port_E->Real()));
  Mpi::GlobalSum(1, &dot, port_nd_fespace->GetComm());
  return dot;
}

std::complex<double> WavePortData::GetVoltage(GridFunction &E) const
{
  // Compute voltage V = ∫ E · dl along the configured path segments.
  // Uses GSLIB interpolation on the 3D parent mesh E field.
  if (!has_voltage_coords)
  {
    return 0.0;
  }
  MFEM_VERIFY(E.HasImag(),
              "Wave ports expect complex-valued E field in port voltage calculation!");
  std::complex<double> V(0.0, 0.0);
  for (std::size_t k = 0; k + 1 < voltage_path.size(); k++)
  {
    V.real(V.real() + fem::ComputeLineIntegral(voltage_path[k], voltage_path[k + 1],
                                               E.Real(), voltage_n_samples));
    V.imag(V.imag() + fem::ComputeLineIntegral(voltage_path[k], voltage_path[k + 1],
                                               E.Imag(), voltage_n_samples));
  }
  return V;
}

std::complex<double> WavePortData::GetExcitationVoltage() const
{
  // TODO: The port mode field (port_E0t) lives on the 2D port submesh, and GSLIB cannot
  // find 3D points on a 2D submesh. To implement this, the port mode field must first be
  // transferred back to the 3D parent mesh before calling ComputeLineIntegral.
  Mpi::Warning("GetExcitationVoltage is not yet implemented for wave port boundaries!\n");
  return 0.0;
}

std::complex<double> WavePortData::GetCharacteristicImpedance() const
{
  // TODO: Same limitation as GetExcitationVoltage — the port mode field lives on the 2D
  // submesh. Requires transfer to parent mesh before GSLIB interpolation can work.
  Mpi::Warning(
      "GetCharacteristicImpedance is not yet implemented for wave port boundaries!\n");
  return 0.0;
}

WavePortOperator::WavePortOperator(const config::BoundaryData &boundaries,
                                   const config::DomainData &domains,
                                   const config::SolverData &solver,
                                   ProblemType problem_type, const Units &units,
                                   const MaterialOperator &mat_op,
                                   mfem::ParFiniteElementSpace &nd_fespace,
                                   mfem::ParFiniteElementSpace &h1_fespace)
  : suppress_output(false), fc(units.Dimensionalize<Units::ValueType::FREQUENCY>(1.0)),
    kc(1.0 / units.Dimensionalize<Units::ValueType::LENGTH>(1.0))
{
  MFEM_VERIFY(nd_fespace.GetParMesh() == h1_fespace.GetParMesh(),
              "Mesh mismatch in WavePortOperator FE spaces!");
  SetUpBoundaryProperties(boundaries, domains, solver, problem_type, units, mat_op,
                          nd_fespace, h1_fespace);
  PrintBoundaryInfo(units, *nd_fespace.GetParMesh());
}

WavePortOperator::WavePortOperator(const IoData &iodata, const MaterialOperator &mat_op,
                                   mfem::ParFiniteElementSpace &nd_fespace,
                                   mfem::ParFiniteElementSpace &h1_fespace)
  : WavePortOperator(iodata.boundaries, iodata.domains, iodata.solver, iodata.problem.type,
                     iodata.units, mat_op, nd_fespace, h1_fespace)
{
}

void WavePortOperator::SetUpBoundaryProperties(const config::BoundaryData &boundaries,
                                               const config::DomainData &domains,
                                               const config::SolverData &solver,
                                               ProblemType problem_type, const Units &units,
                                               const MaterialOperator &mat_op,
                                               mfem::ParFiniteElementSpace &nd_fespace,
                                               mfem::ParFiniteElementSpace &h1_fespace)
{

  // Check that wave port boundary attributes have been specified correctly.
  const auto &mesh = *nd_fespace.GetParMesh();
  MFEM_VERIFY(boundaries.waveport.empty() || mesh.Dimension() == 3,
              "Wave port boundaries are only available for 3D simulations!");
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  if (!boundaries.waveport.empty())
  {
    mfem::Array<int> bdr_attr_marker(bdr_attr_max), port_marker(bdr_attr_max);
    bdr_attr_marker = 0;
    port_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    for (const auto &[idx, data] : boundaries.waveport)
    {
      for (auto attr : data.attributes)
      {
        MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
                    "Port boundary attribute tags must be non-negative and correspond to "
                    "boundaries in the mesh!");
        MFEM_VERIFY(bdr_attr_marker[attr - 1],
                    "Unknown port boundary attribute " << attr << "!");
        MFEM_VERIFY(!data.active || !port_marker[attr - 1],
                    "Boundary attribute is assigned to more than one wave port!");
        port_marker[attr - 1] = 1;
      }
    }
  }

  // List of all boundaries which will be marked as essential (Dirichlet) for the purposes
  // of computing wave port modes: PEC and AuxPEC surfaces. Other wave ports are also marked
  // as Dirichlet in case two wave ports are touching and share one or more edges.
  // Impedance, conductivity, and absorbing BCs are handled through their respective
  // boundary operators in the eigenvalue problem.
  mfem::Array<int> dbc_bcs;
  dbc_bcs.Reserve(static_cast<int>(boundaries.pec.attributes.size() +
                                   boundaries.auxpec.attributes.size()));
  for (auto attr : boundaries.pec.attributes)
  {
    if (attr <= 0 || attr > bdr_attr_max)
    {
      continue;
    }
    dbc_bcs.Append(attr);
  }
  for (auto attr : boundaries.auxpec.attributes)
  {
    if (attr <= 0 || attr > bdr_attr_max)
    {
      continue;
    }
    dbc_bcs.Append(attr);
  }
  // If user accidentally specifies a surface as both "PEC" and "WavePortPEC", this is fine
  // so allow for duplicates in the attribute list.
  dbc_bcs.Sort();
  dbc_bcs.Unique();

  // Set up wave port data structures.
  for (const auto &[idx, data] : boundaries.waveport)
  {
    mfem::Array<int> port_dbc_bcs(dbc_bcs);
    for (const auto &[other_idx, other_data] : boundaries.waveport)
    {
      if (other_idx == idx || !other_data.active)
      {
        continue;
      }
      for (auto attr : other_data.attributes)
      {
        if (std::binary_search(data.attributes.begin(), data.attributes.end(), attr))
        {
          continue;
        }
        port_dbc_bcs.Append(attr);
      }
    }
    port_dbc_bcs.Sort();
    port_dbc_bcs.Unique();
    ports.try_emplace(idx, data, boundaries, domains, problem_type, solver.linear, units,
                      mat_op, nd_fespace, h1_fespace, port_dbc_bcs);
  }
  MFEM_VERIFY(
      ports.empty() || problem_type == ProblemType::DRIVEN ||
          problem_type == ProblemType::EIGENMODE,
      "Wave port boundaries are only available for frequency domain driven simulations!");
}

void WavePortOperator::PrintBoundaryInfo(const Units &units, const mfem::ParMesh &mesh)
{
  if (ports.empty())
  {
    return;
  }
  fmt::memory_buffer buffer{};
  auto out = fmt::appender{buffer};

  // Print out BC info for all active port attributes.
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    for (auto attr : data.GetAttrList())
    {
      fmt::format_to(
          out, " {:d}: Index = {:d}, mode = {:d}, d = {:.3e} m,  n = ({:+.1f})\n", attr,
          idx, data.mode_idx, units.Dimensionalize<Units::ValueType::LENGTH>(data.d_offset),
          fmt::join(data.port_normal, ","));
    }
  }
  if (buffer.size() > 0)
  {
    Mpi::Print("\nConfiguring Robin impedance BC for wave ports at attributes:\n");
    Mpi::Print("{}", fmt::to_string(buffer));
    buffer.clear();
  }

  // Print some information for excited wave ports.
  for (const auto &[idx, data] : ports)
  {
    if (!data.HasExcitation())
    {
      continue;
    }
    for (auto attr : data.GetAttrList())
    {
      fmt::format_to(out, " {:d}: Index = {:d}\n", attr, idx);
    }
  }
  if (buffer.size() > 0)
  {
    Mpi::Print("\nConfiguring wave port excitation source term at attributes:\n");
    Mpi::Print("{}", fmt::to_string(buffer));
  }
}

const WavePortData &WavePortOperator::GetPort(int idx) const
{
  auto it = ports.find(idx);
  MFEM_VERIFY(it != ports.end(), "Unknown wave port index requested!");
  return it->second;
}

mfem::Array<int> WavePortOperator::GetAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    attr_list.Append(data.GetAttrList());
  }
  return attr_list;
}

void WavePortOperator::Initialize(double omega)
{
  bool init = false, first = true;
  for (const auto &[idx, data] : ports)
  {
    init = init || (data.omega0 != omega);
    first = first && (data.omega0 == 0.0);
  }
  if (!init)
  {
    return;
  }
  BlockTimer bt(Timer::WAVE_PORT);
  if (!suppress_output)
  {
    Mpi::Print(
        "\nCalculating boundary modes at wave ports for ω/2π = {:.3e} GHz ({:.3e})\n",
        omega * fc / (2.0 * M_PI), omega);
  }
  for (auto &[idx, data] : ports)
  {
    data.Initialize(omega);
    if (!suppress_output)
    {
      if (first)
      {
        Mpi::Print(" Number of global unknowns for port {:d}:\n"
                   "  H1: {:d}, ND: {:d}\n",
                   idx, data.GlobalTrueH1Size(), data.GlobalTrueNDSize());
      }
      Mpi::Print(" Port {:d}, mode {:d}: kₙ = {:.3e}{:+.3e}i m⁻¹\n", idx, data.mode_idx,
                 data.kn0.real() * kc, data.kn0.imag() * kc);
    }
  }
}

void WavePortOperator::AddExtraSystemBdrCoefficients(double omega,
                                                     MaterialPropertyCoefficient &fbr,
                                                     MaterialPropertyCoefficient &fbi)
{
  // Add wave port boundaries to the bilinear form. This looks a lot like the lumped port
  // boundary, except the iω / Z_s coefficient goes to ikₙ / μ where kₙ is specific to the
  // port mode at the given operating frequency (note only the real part of the propagation
  // constant contributes).
  Initialize(omega);
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    const MaterialOperator &mat_op = data.mat_op;
    MaterialPropertyCoefficient muinv_func(mat_op.GetBdrAttributeToMaterial(),
                                           mat_op.GetInvPermeability());
    muinv_func.RestrictCoefficient(mat_op.GetCeedBdrAttributes(data.GetAttrList()));
    // fbr.AddCoefficient(muinv_func.GetAttributeToMaterial(),
    //                    muinv_func.GetMaterialProperties(),
    //                    -data.kn0.imag());
    fbi.AddCoefficient(muinv_func.GetAttributeToMaterial(),
                       muinv_func.GetMaterialProperties(), data.kn0.real());
  }
}

void WavePortOperator::AddExcitationBdrCoefficients(int excitation_idx, double omega,
                                                    SumVectorCoefficient &fbr,
                                                    SumVectorCoefficient &fbi)
{
  // Re/Im{-U_inc} = Re/Im{+2 (-iω) n x H_inc}, which is a function of E_inc as computed by
  // the modal solution (stored as a grid function and coefficient during initialization).
  Initialize(omega);
  for (const auto &[idx, data] : ports)
  {
    if (data.excitation != excitation_idx)
    {
      continue;
    }
    fbr.AddCoefficient(data.GetModeExcitationCoefficientImag(), 2.0 * omega);
    fbi.AddCoefficient(data.GetModeExcitationCoefficientReal(), -2.0 * omega);
  }
}

// === ModeEigenSolver implementation (merged from modeeigensolver.cpp) ===

namespace
{

constexpr bool skip_zeros = false;

}  // namespace

ModeEigenSolver::ModeEigenSolver(
    const MaterialOperator &mat_op, const mfem::Vector *normal,
    SurfaceImpedanceOperator *surf_z_op, FarfieldBoundaryOperator *farfield_op,
    SurfaceConductivityOperator *surf_sigma_op, const FiniteElementSpace &nd_fespace,
    const FiniteElementSpace &h1_fespace, const mfem::Array<int> &dbc_tdof_list,
    int num_modes, int num_vec, double eig_tol, EigenvalueSolver::WhichType which_eig,
    const config::LinearSolverData &linear, EigenSolverBackend eigen_backend, int verbose,
    MPI_Comm solver_comm, const ModeEigenSolverMultigridConfig *mg_config)
  : num_modes(num_modes), num_vec(num_vec), eig_tol(eig_tol), which_eig(which_eig),
    linear(linear), eigen_backend(eigen_backend), verbose(verbose), mat_op(mat_op),
    normal(normal), surf_z_op(surf_z_op), farfield_op(farfield_op),
    surf_sigma_op(surf_sigma_op), nd_fespace(nd_fespace), h1_fespace(h1_fespace),
    dbc_tdof_list(dbc_tdof_list), mg_config(mg_config)
{
  nd_size = nd_fespace.GetTrueVSize();
  h1_size = h1_fespace.GetTrueVSize();

  // Assemble frequency-independent matrices. These use the FE space communicator
  // (all processes that share the mesh), NOT solver_comm.
  //
  // Atn: gradient coupling -(mu^{-1} grad_t u, v).
  std::tie(Atnr, Atni) = AssembleAtn();

  // Btn: NEGATED transpose of Atn. Atn = -(mu^{-1} grad ·, ·), so Atn^T = -(mu^{-1} ·,
  // grad·). The physical Btn from Eq 2 is +(mu^{-1} ·, grad·) (positive), so Btn = -Atn^T.
  Btnr.reset(Atnr->Transpose());
  *Btnr *= -1.0;
  if (Atni)
  {
    Btni.reset(Atni->Transpose());
    *Btni *= -1.0;
  }

  // Assemble Btt (negative mu^{-1} ND mass) and build the block B matrix.
  {
    auto [bttr, btti] = AssembleBtt();
    Bttr = std::make_unique<mfem::HypreParMatrix>(*bttr);  // Keep a copy for external use

    // Construct the zero diagonal block for the H1 portion of B.
    Vector d(h1_size);
    d.UseDevice(false);  // SparseMatrix constructor uses Vector on host
    d = 0.0;
    mfem::SparseMatrix diag(d);
    auto Dnn = std::make_unique<mfem::HypreParMatrix>(
        h1_fespace.Get().GetComm(), h1_fespace.Get().GlobalTrueVSize(),
        h1_fespace.Get().GetTrueDofOffsets(), &diag);

    auto [Br, Bi] =
        BuildSystemMatrixB(bttr.get(), btti.get(), Btnr.get(), Btni.get(), Dnn.get());
    opB = std::make_unique<ComplexWrapperOperator>(std::move(Br), std::move(Bi));
  }

  // Configure linear and eigenvalue solvers.
  // - BoundaryMode: solver_comm == MPI_COMM_NULL -> use FE space communicator (all procs).
  // - WavePort on port procs: solver_comm is a valid subset communicator -> use it.
  // - WavePort on non-port procs: solver_comm == MPI_COMM_NULL but FE space has no local
  //   DOFs -> skip solver setup (ksp and eigen remain null).
  const bool use_mg =
      mg_config && mg_config->nd_fespaces && mg_config->nd_fespaces->GetNumLevels() > 1;
  if (solver_comm != MPI_COMM_NULL)
  {
    if (use_mg)
    {
      SetUpMultigridLinearSolver(solver_comm);
    }
    else
    {
      SetUpLinearSolver(solver_comm);
    }
    SetUpEigenSolver(solver_comm);
  }
  else if (nd_size > 0)
  {
    // Standalone mode (BoundaryMode): all procs have DOFs, use FE space comm.
    if (use_mg)
    {
      SetUpMultigridLinearSolver(nd_fespace.GetComm());
    }
    else
    {
      SetUpLinearSolver(nd_fespace.GetComm());
    }
    SetUpEigenSolver(nd_fespace.GetComm());
  }
  // else: non-port process with empty FE space -- no solvers needed.
}

ModeEigenSolver::~ModeEigenSolver()
{
  // Free the solvers before any communicator they depend on is freed externally.
  ksp.reset();
  eigen.reset();
}

void ModeEigenSolver::AssembleFrequencyDependent(double omega, double sigma)
{
  auto [Attr, Atti] = AssembleAtt(omega, sigma);
  auto [Annr_local, Anni_local] = AssembleAnn(omega);

  // Compute the shifted (2,1) block: -sigma * Btn_r. For sigma = -kn_target^2, this
  // gives +kn_target^2 * Btn_r. Btn is real-only (no imaginary part).
  std::unique_ptr<mfem::HypreParMatrix> shifted_Btnr;
  if (Btnr && std::abs(sigma) > 0.0)
  {
    shifted_Btnr = std::make_unique<mfem::HypreParMatrix>(*Btnr);
    *shifted_Btnr *= -sigma;
  }

  auto [Ar, Ai] =
      BuildSystemMatrixA(Attr.get(), Atti.get(), Atnr.get(), Atni.get(), Annr_local.get(),
                         Anni_local.get(), shifted_Btnr.get());
  opA = std::make_unique<ComplexWrapperOperator>(std::move(Ar), std::move(Ai));
}

ModeEigenSolver::SolveResult ModeEigenSolver::Solve(double omega, double sigma,
                                                    bool has_solver,
                                                    const ComplexVector *initial_space)
{
  // Assemble frequency-dependent matrices (MPI collective on FE space communicator).
  AssembleFrequencyDependent(omega, sigma);

  // Only processes with solvers participate in the eigenvalue solve. For BoundaryMode
  // (has_solver=true on all processes) all processes solve. For WavePort, only port
  // processes solve while non-port processes skip.
  if (!has_solver)
  {
    return {0, sigma};
  }

  MFEM_VERIFY(ksp && eigen, "ModeEigenSolver::Solve called on process without solvers!");

  if (block_pc_ptr)
  {
    // Multigrid path: assemble preconditioner operators at all levels and set on the
    // block-diagonal preconditioner. The outer Krylov solver gets the monolithic opA.
    att_mg_op = AssembleAttPreconditioner(omega, sigma);
    ann_mg_op = AssembleAnnPreconditioner(omega);
    block_pc_ptr->SetBlockOperators(*att_mg_op, *ann_mg_op);

    // Set the off-diagonal operator -sigma*Btn for block lower-triangular preconditioning.
    // This captures the shift-and-invert coupling that dominates the off-diagonal.
    if (Btnr && std::abs(sigma) > 0.0)
    {
      auto sBtnr = std::make_unique<mfem::HypreParMatrix>(*Btnr);
      *sBtnr *= -sigma;
      shifted_Btn_op = std::make_unique<ComplexWrapperOperator>(std::move(sBtnr), nullptr);
      block_pc_ptr->SetOffDiagonalOperator(shifted_Btn_op.get());
    }
    else
    {
      block_pc_ptr->SetOffDiagonalOperator(nullptr);
    }

    ksp->SetOperators(*opA, *opA);  // opA passed twice; pc uses block_pc_ptr
  }
  else
  {
    // Sparse direct path: precondition with real part of the full block system.
    ComplexWrapperOperator opP(opA->Real(), nullptr);
    ksp->SetOperators(*opA, opP);
  }
  eigen->SetOperators(*opB, *opA, EigenvalueSolver::ScaleType::NONE);

  if (initial_space)
  {
    eigen->SetInitialSpace(*initial_space);
  }

  int num_conv = eigen->Solve();

  // Build a permutation sorted by proximity to the shift target so that mode ordering is
  // consistent across eigensolver backends (ARPACK vs SLEPc sort eigenvalues differently).
  // The shift sigma = -kn_target^2, so kn_target = sqrt(-sigma). Sorting by ascending
  // |Re{kn} - kn_target| puts the mode closest to the target first.
  const double kn_target = std::sqrt(-sigma);
  mode_perm.resize(num_conv);
  std::iota(mode_perm.begin(), mode_perm.end(), 0);
  std::sort(mode_perm.begin(), mode_perm.end(),
            [this, sigma, kn_target](int a, int b)
            {
              auto kn_a = std::sqrt(-sigma - 1.0 / eigen->GetEigenvalue(a));
              auto kn_b = std::sqrt(-sigma - 1.0 / eigen->GetEigenvalue(b));
              return std::abs(kn_a.real() - kn_target) < std::abs(kn_b.real() - kn_target);
            });

  return {num_conv, sigma};
}

std::complex<double> ModeEigenSolver::GetEigenvalue(int i) const
{
  return eigen->GetEigenvalue(mode_perm[i]);
}

double ModeEigenSolver::GetError(int i, EigenvalueSolver::ErrorType type) const
{
  return eigen->GetError(mode_perm[i], type);
}

void ModeEigenSolver::GetEigenvector(int i, ComplexVector &x) const
{
  eigen->GetEigenvector(mode_perm[i], x);
}

// Stiffness matrix (shifted): Att = (mu_cc^{-1} curl_t x u, curl_t x v)
//                                   - omega^2 (eps u, v)
//                                   - sigma (mu^{-1} u, v)
//                                   + BC-t impedance/absorbing/conductivity.
ModeEigenSolver::ComplexHypreParMatrix ModeEigenSolver::AssembleAtt(double omega,
                                                                    double sigma) const
{
  MaterialPropertyCoefficient muinv_cc_func(mat_op.GetAttributeToMaterial(),
                                            normal ? mat_op.GetInvPermeability()
                                                   : mat_op.GetCurlCurlInvPermeability());
  if (normal)
  {
    muinv_cc_func.NormalProjectedCoefficient(*normal);
  }

  MaterialPropertyCoefficient eps_shifted_func(
      mat_op.GetAttributeToMaterial(), mat_op.GetPermittivityReal(), -omega * omega);
  eps_shifted_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                  mat_op.GetInvPermeability(), -sigma);
  // London superconductor contribution: +(1/lambda_L^2)(et, ft).
  if (mat_op.HasLondonDepth())
  {
    eps_shifted_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                    mat_op.GetInvLondonDepth(), 1.0);
  }

  // Assemble Att real part: domain + boundary contributions in a single BilinearForm,
  // following the same pattern as SpaceOperator. Boundary impedance/absorbing terms use
  // Palace's ceed-based BilinearForm with VectorFEMassIntegrator.
  int max_bdr_attr = mat_op.MaxCeedBdrAttribute();
  MaterialPropertyCoefficient fbr(max_bdr_attr), fbi(max_bdr_attr);

  // Add boundary contributions from impedance, absorbing, and conductivity operators.
  if (surf_z_op)
  {
    surf_z_op->AddStiffnessBdrCoefficients(1.0, fbr);        // Ls -> real
    surf_z_op->AddDampingBdrCoefficients(omega, fbi);        // Rs -> imag
    surf_z_op->AddMassBdrCoefficients(-omega * omega, fbr);  // Cs -> real
  }
  if (farfield_op)
  {
    farfield_op->AddDampingBdrCoefficients(omega, fbi);  // 1st-order ABC
  }
  if (surf_sigma_op)
  {
    surf_sigma_op->AddExtraSystemBdrCoefficients(omega, fbr, fbi);
  }

  BilinearForm att(nd_fespace);
  att.AddDomainIntegrator<CurlCurlMassIntegrator>(muinv_cc_func, eps_shifted_func);
  if (!fbr.empty())
  {
    att.AddBoundaryIntegrator<VectorFEMassIntegrator>(fbr);
  }

  // Assemble domain + boundary contributions via Palace's BilinearForm (libCEED).
  auto Attr_assembled =
      ParOperator(att.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble();
  // Assemble imaginary part: domain (loss tangent, conductivity) + boundary (Rs, ABC).
  std::unique_ptr<mfem::HypreParMatrix> Atti_assembled;
  {
    bool has_imag = mat_op.HasLossTangent() || mat_op.HasConductivity() || !fbi.empty();
    if (has_imag)
    {
      // Coefficients must outlive the BilinearForm (integrators store raw pointers).
      int n_attr = mat_op.GetAttributeToMaterial().Size();
      MaterialPropertyCoefficient negepstandelta_func(n_attr);
      MaterialPropertyCoefficient fi_domain(n_attr);
      if (mat_op.HasLossTangent())
      {
        negepstandelta_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                           mat_op.GetPermittivityImag(), -omega * omega);
      }
      if (mat_op.HasConductivity())
      {
        fi_domain.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetConductivity(),
                                 omega);
      }
      BilinearForm atti(nd_fespace);
      if (!negepstandelta_func.empty())
      {
        atti.AddDomainIntegrator<VectorFEMassIntegrator>(negepstandelta_func);
      }
      if (!fi_domain.empty())
      {
        atti.AddDomainIntegrator<VectorFEMassIntegrator>(fi_domain);
      }
      if (!fbi.empty())
      {
        atti.AddBoundaryIntegrator<VectorFEMassIntegrator>(fbi);
      }
      Atti_assembled =
          ParOperator(atti.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble();
    }
  }

  return {std::move(Attr_assembled), std::move(Atti_assembled)};
}

// Coupling matrix: Atn = -(mu^{-1} grad_t u, v).
ModeEigenSolver::ComplexHypreParMatrix ModeEigenSolver::AssembleAtn() const
{
  MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                         mat_op.GetInvPermeability(), -1.0);
  BilinearForm atn(h1_fespace, nd_fespace);
  atn.AddDomainIntegrator<MixedVectorGradientIntegrator>(muinv_func);
  return {ParOperator(atn.FullAssemble(skip_zeros), h1_fespace, nd_fespace, false)
              .StealParallelAssemble(),
          nullptr};
}

// H1 block: Ann = -(mu^{-1} grad u, grad v) + omega^2(eps u, v) + BC-n impedance.
//
// This replaces ModeEigenSolver's Ann (which is just -eps mass) with the full normal
// curl-curl equation from Equation 2. The stiffness (diffusion) is NEGATIVE (IBP of
// Laplacian), the mass is POSITIVE (+omega^2 eps), and impedance boundary terms use
// NEGATED signs relative to the tangential Att terms (from the double cross product
// in the impedance BC: [n x (n x E)]_z = -En).
ModeEigenSolver::ComplexHypreParMatrix ModeEigenSolver::AssembleAnn(double omega) const
{
  // H1 stiffness: -(mu^{-1} grad u, grad v) -- negative sign from IBP of Laplacian.
  MaterialPropertyCoefficient neg_muinv_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetInvPermeability(), -1.0);
  if (normal)
  {
    neg_muinv_func.NormalProjectedCoefficient(*normal);
  }

  // H1 mass: +omega^2(eps u, v) -- positive sign in the normal equation.
  MaterialPropertyCoefficient poseps_h1_func(mat_op.GetAttributeToMaterial(),
                                             normal ? mat_op.GetPermittivityReal()
                                                    : mat_op.GetPermittivityScalar(),
                                             omega * omega);
  if (normal)
  {
    poseps_h1_func.NormalProjectedCoefficient(*normal);
  }

  // London superconductor contribution: +(1/lambda_L^2)(en, fn).
  // In the normal curl-curl equation, the London term is +(1/lambda_L^2) En, which
  // enters as a positive H1 mass (same sign as the omega^2 eps mass).
  if (mat_op.HasLondonDepth())
  {
    if (!normal)
    {
      // 2D mode analysis: use pre-computed scalar out-of-plane London depth.
      poseps_h1_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                    mat_op.GetInvLondonDepthScalar());
    }
    else if (normal)
    {
      // 3D wave port: use normal projection of the full tensor (handled by caller via
      // NormalProjectedCoefficient on the assembled form).
      const auto &ild = mat_op.GetInvLondonDepth();
      mfem::DenseTensor ild_scalar(1, 1, ild.SizeK());
      for (int k = 0; k < ild.SizeK(); k++)
      {
        ild_scalar(0, 0, k) = ild(0, 0, k);
      }
      poseps_h1_func.AddCoefficient(mat_op.GetAttributeToMaterial(), ild_scalar);
      poseps_h1_func.NormalProjectedCoefficient(*normal);
    }
  }

  // Boundary impedance for en: -(iw/Zs)(en, fn)_gamma -- NEGATIVE sign from the
  // double cross product in the impedance BC: [n x (n x E)]_z = -En, giving
  // -(iw/Zs)(En, fn) on the boundary. Note: opposite sign from et boundary term.
  // This follows the same pattern as the quadratic formulation's M0(2,2) assembly.
  int max_bdr_attr = mat_op.MaxCeedBdrAttribute();
  MaterialPropertyCoefficient nn_fbr(max_bdr_attr), nn_fbi(max_bdr_attr);
  if (surf_z_op)
  {
    surf_z_op->AddStiffnessBdrCoefficients(-1.0, nn_fbr);
    surf_z_op->AddDampingBdrCoefficients(-omega, nn_fbi);
    surf_z_op->AddMassBdrCoefficients(omega * omega, nn_fbr);
  }
  if (farfield_op && farfield_op->GetAttrList().Size() > 0)
  {
    // The farfield operator's AddDampingBdrCoefficients adds a tensor-valued inverse
    // impedance (for VectorFEMassIntegrator on ND). For the scalar H1 MassIntegrator,
    // we need the scalar inverse impedance sqrt(eps/mu). Extract the (0,0) component
    // of the tensor for each boundary material and add as scalar coefficients.
    const auto &farfield_attrs = farfield_op->GetAttrList();
    const auto &inv_z = mat_op.GetInvImpedance();
    const auto &bdr_attr_to_mat = mat_op.GetBdrAttributeToMaterial();
    for (auto attr : farfield_attrs)
    {
      int mat_idx =
          (attr > 0 && attr <= bdr_attr_to_mat.Size()) ? bdr_attr_to_mat[attr - 1] : -1;
      double inv_z0_scalar = (mat_idx >= 0) ? inv_z(0, 0, mat_idx) : 1.0;
      auto ceed_attrs = mat_op.GetCeedBdrAttributes(attr);
      if (ceed_attrs.Size() > 0)
      {
        nn_fbi.AddMaterialProperty(ceed_attrs, inv_z0_scalar, -omega);
      }
    }
  }
  if (surf_sigma_op)
  {
    // For conductivity: negate the coefficients (opposite sign from et boundary term).
    MaterialPropertyCoefficient cond_r(max_bdr_attr), cond_i(max_bdr_attr);
    surf_sigma_op->AddExtraSystemBdrCoefficients(omega, cond_r, cond_i);
    if (!cond_r.empty())
    {
      cond_r *= -1.0;
      nn_fbr.AddCoefficient(cond_r.GetAttributeToMaterial(),
                            cond_r.GetMaterialProperties());
    }
    if (!cond_i.empty())
    {
      cond_i *= -1.0;
      nn_fbi.AddCoefficient(cond_i.GetAttributeToMaterial(),
                            cond_i.GetMaterialProperties());
    }
  }

  // Assemble Ann real: H1 stiffness (negative Laplacian) + mass (positive) + boundary.
  // Use DiffusionMassIntegrator to combine stiffness and mass in a single form.
  BilinearForm annr(h1_fespace);
  annr.AddDomainIntegrator<DiffusionMassIntegrator>(neg_muinv_func, poseps_h1_func);
  if (!nn_fbr.empty())
  {
    annr.AddBoundaryIntegrator<MassIntegrator>(nn_fbr);
  }
  auto Annr_assembled =
      ParOperator(annr.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble();

  // Assemble Ann imaginary: loss tangent + boundary impedance.
  std::unique_ptr<mfem::HypreParMatrix> Anni_assembled;
  {
    bool has_imag = mat_op.HasLossTangent() || !nn_fbi.empty();
    if (has_imag)
    {
      // Coefficient must outlive the BilinearForm (integrators store raw pointers).
      int n_attr = mat_op.GetAttributeToMaterial().Size();
      MaterialPropertyCoefficient posepsi_h1_func(n_attr);
      if (mat_op.HasLossTangent())
      {
        if (normal)
        {
          // 3D wave port: project the full tensor using the surface normal.
          posepsi_h1_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                         mat_op.GetPermittivityImag(), omega * omega);
          posepsi_h1_func.NormalProjectedCoefficient(*normal);
        }
        else if (!normal)
        {
          // 2D mode analysis: use pre-computed scalar imaginary permittivity.
          posepsi_h1_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                         mat_op.GetPermittivityImagScalar(), omega * omega);
        }
      }
      BilinearForm anni(h1_fespace);
      if (!posepsi_h1_func.empty())
      {
        anni.AddDomainIntegrator<MassIntegrator>(posepsi_h1_func);
      }
      if (!nn_fbi.empty())
      {
        anni.AddBoundaryIntegrator<MassIntegrator>(nn_fbi);
      }
      Anni_assembled =
          ParOperator(anni.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble();
    }
  }

  return {std::move(Annr_assembled), std::move(Anni_assembled)};
}

// Mass matrix: Btt = (mu^{-1} u, v).
ModeEigenSolver::ComplexHypreParMatrix ModeEigenSolver::AssembleBtt() const
{
  MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                         mat_op.GetInvPermeability());
  BilinearForm btt(nd_fespace);
  btt.AddDomainIntegrator<VectorFEMassIntegrator>(muinv_func);
  return {ParOperator(btt.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble(),
          nullptr};
}

ModeEigenSolver::ComplexHypreParMatrix ModeEigenSolver::BuildSystemMatrixA(
    const mfem::HypreParMatrix *Attr, const mfem::HypreParMatrix *Atti,
    const mfem::HypreParMatrix *Atnr, const mfem::HypreParMatrix *Atni,
    const mfem::HypreParMatrix *Annr, const mfem::HypreParMatrix *Anni,
    const mfem::HypreParMatrix *shifted_Btnr) const
{
  // Construct the 2x2 block matrices for the eigenvalue problem A e = lambda B e.
  // The (1,0) block is -sigma * Btn from the shift-and-invert transformation.
  // Without shift (sigma=0), this block is zero (upper block-triangular).
  mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
  blocks(0, 0) = Attr;
  blocks(0, 1) = Atnr;
  blocks(1, 0) = shifted_Btnr;  // -sigma * Btn (nullptr when sigma=0)
  blocks(1, 1) = Annr;
  std::unique_ptr<mfem::HypreParMatrix> Ar(mfem::HypreParMatrixFromBlocks(blocks));

  std::unique_ptr<mfem::HypreParMatrix> Ai;
  if (Atti || Atni || Anni)
  {
    // HypreParMatrixFromBlocks requires at least one non-null block per row and column
    // to determine sizes. Since (1,0) is always null (shifted Btn is real-only), add
    // zero diagonal placeholders when an entire block row or column would be null.
    std::unique_ptr<mfem::HypreParMatrix> Dtt_zero, Dnn_zero;
    if (!Atti && !Atni)
    {
      Vector d(nd_size);
      d.UseDevice(false);
      d = 0.0;
      mfem::SparseMatrix diag(d);
      Dtt_zero = std::make_unique<mfem::HypreParMatrix>(
          nd_fespace.Get().GetComm(), nd_fespace.Get().GlobalTrueVSize(),
          nd_fespace.Get().GetTrueDofOffsets(), &diag);
    }
    if (!Anni)
    {
      Vector d(h1_size);
      d.UseDevice(false);
      d = 0.0;
      mfem::SparseMatrix diag(d);
      Dnn_zero = std::make_unique<mfem::HypreParMatrix>(
          h1_fespace.Get().GetComm(), h1_fespace.Get().GlobalTrueVSize(),
          h1_fespace.Get().GetTrueDofOffsets(), &diag);
    }
    blocks(0, 0) = Atti ? Atti : Dtt_zero.get();
    blocks(0, 1) = Atni;
    blocks(1, 0) = nullptr;  // Shifted Btn is real-only
    blocks(1, 1) = Anni ? Anni : Dnn_zero.get();
    Ai.reset(mfem::HypreParMatrixFromBlocks(blocks));
  }

  // Eliminate boundary true dofs constrained by Dirichlet BCs.
  Ar->EliminateBC(dbc_tdof_list, Operator::DIAG_ONE);
  if (Ai)
  {
    Ai->EliminateBC(dbc_tdof_list, Operator::DIAG_ZERO);
  }

  return {std::move(Ar), std::move(Ai)};
}

ModeEigenSolver::ComplexHypreParMatrix ModeEigenSolver::BuildSystemMatrixB(
    const mfem::HypreParMatrix *Bttr, const mfem::HypreParMatrix *Btti,
    const mfem::HypreParMatrix *Btnr, const mfem::HypreParMatrix *Btni,
    const mfem::HypreParMatrix *Dnn) const
{
  // Construct the 2x2 block matrices for the eigenvalue problem A e = lambda B e.
  // B = [Btt, 0; Btn, Dnn] where Btn = Atn^T and Dnn is the zero diagonal.
  mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
  blocks(0, 0) = Bttr;
  blocks(0, 1) = nullptr;
  blocks(1, 0) = Btnr;
  blocks(1, 1) = Dnn;
  std::unique_ptr<mfem::HypreParMatrix> Br(mfem::HypreParMatrixFromBlocks(blocks));

  std::unique_ptr<mfem::HypreParMatrix> Bi;
  if (Btti || Btni)
  {
    // NOTE: Currently unreachable (Btt, Btn are real for real permeability). If complex
    // permeability is added, zero placeholder blocks would be needed here too (same as
    // the imaginary A block above) to prevent HypreParMatrixFromBlocks sizing errors.
    blocks(0, 0) = Btti;
    blocks(0, 1) = nullptr;
    blocks(1, 0) = Btni;
    blocks(1, 1) = nullptr;
    Bi.reset(mfem::HypreParMatrixFromBlocks(blocks));
  }

  // Eliminate boundary true dofs constrained by Dirichlet BCs.
  Br->EliminateBC(dbc_tdof_list, Operator::DIAG_ZERO);
  if (Bi)
  {
    Bi->EliminateBC(dbc_tdof_list, Operator::DIAG_ZERO);
  }

  return {std::move(Br), std::move(Bi)};
}

void ModeEigenSolver::SetUpLinearSolver(MPI_Comm comm)
{
  // GMRES iterative solver preconditioned with a sparse direct solver.
  auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(comm, verbose);
  gmres->SetInitialGuess(false);
  gmres->SetRelTol(linear.tol);
  gmres->SetMaxIter(linear.max_it);
  gmres->SetRestartDim(linear.max_size);
  gmres->EnableTimer();

  // Select sparse direct solver type. Multigrid (AMS, BoomerAMG) is not applicable for
  // the combined ND+H1 block system, so fall back to a sparse direct solver.
  LinearSolver pc_type = linear.type;
  if (pc_type == LinearSolver::DEFAULT || pc_type == LinearSolver::AMS ||
      pc_type == LinearSolver::BOOMER_AMG)
  {
#if defined(MFEM_USE_SUPERLU)
    pc_type = LinearSolver::SUPERLU;
#elif defined(MFEM_USE_STRUMPACK)
    pc_type = LinearSolver::STRUMPACK;
#elif defined(MFEM_USE_MUMPS)
    pc_type = LinearSolver::MUMPS;
#else
    MFEM_ABORT("ModeEigenSolver requires building with SuperLU_DIST, STRUMPACK, "
               "or MUMPS!");
#endif
  }
  else if (pc_type == LinearSolver::SUPERLU)
  {
#if !defined(MFEM_USE_SUPERLU)
    MFEM_ABORT("Solver was not built with SuperLU_DIST support, please choose a "
               "different solver!");
#endif
  }
  else if (pc_type == LinearSolver::STRUMPACK || pc_type == LinearSolver::STRUMPACK_MP)
  {
#if !defined(MFEM_USE_STRUMPACK)
    MFEM_ABORT("Solver was not built with STRUMPACK support, please choose a "
               "different solver!");
#endif
  }
  else if (pc_type == LinearSolver::MUMPS)
  {
#if !defined(MFEM_USE_MUMPS)
    MFEM_ABORT("Solver was not built with MUMPS support, please choose a "
               "different solver!");
#endif
  }

  auto pc = std::make_unique<MfemWrapperSolver<ComplexOperator>>(
      [&]() -> std::unique_ptr<mfem::Solver>
      {
        if (pc_type == LinearSolver::SUPERLU)
        {
#if defined(MFEM_USE_SUPERLU)
          return std::make_unique<SuperLUSolver>(comm, linear.sym_factorization,
                                                 linear.superlu_3d, true, verbose - 1);
#endif
        }
        else if (pc_type == LinearSolver::STRUMPACK ||
                 pc_type == LinearSolver::STRUMPACK_MP)
        {
#if defined(MFEM_USE_STRUMPACK)
          return std::make_unique<StrumpackSolver>(
              comm, linear.sym_factorization, linear.strumpack_compression_type,
              linear.strumpack_lr_tol, linear.strumpack_butterfly_l,
              linear.strumpack_lossy_precision, true, verbose - 1);
#endif
        }
        else if (pc_type == LinearSolver::MUMPS)
        {
#if defined(MFEM_USE_MUMPS)
          return std::make_unique<MumpsSolver>(comm, MatrixSymmetry::UNSYMMETRIC,
                                               linear.sym_factorization,
                                               linear.strumpack_lr_tol, true, verbose - 1);
#endif
        }
        MFEM_ABORT("Unsupported linear solver type for boundary mode solver!");
        return {};
      }());
  pc->SetSaveAssembled(false);
  pc->SetDropSmallEntries(false);
  ksp = std::make_unique<ComplexKspSolver>(std::move(gmres), std::move(pc));
}

void ModeEigenSolver::SetUpMultigridLinearSolver(MPI_Comm comm)
{
  MFEM_VERIFY(mg_config && mg_config->nd_fespaces && mg_config->h1_fespaces &&
                  mg_config->h1_aux_fespaces,
              "Multigrid linear solver requires ND, H1, and H1 auxiliary space "
              "hierarchies!");
  const int print = verbose - 1;

  // Determine coarse solver types for the ND and H1 blocks from the config. The default
  // type is resolved in IoData (sparse direct for frequency-domain problems, AMS
  // otherwise). For the H1 block, AMS is not applicable — use BoomerAMG instead.
  LinearSolver nd_pc_type = linear.type;
  LinearSolver h1_pc_type = linear.type;

  if (nd_pc_type == LinearSolver::BOOMER_AMG)
  {
    Mpi::Warning(comm,
                 "BoomerAMG is not well-suited for the Nedelec system matrix, consider "
                 "using another solver.\n");
  }
  if (h1_pc_type == LinearSolver::AMS)
  {
    h1_pc_type = LinearSolver::BOOMER_AMG;
    Mpi::Print(" Multigrid coarse solve: AMS for ND block, BoomerAMG for H1 block\n");
  }

  // Helper to create a coarse solver from the type.
  auto MakeCoarseSolver =
      [&](LinearSolver type) -> std::unique_ptr<MfemWrapperSolver<ComplexOperator>>
  {
    switch (type)
    {
      case LinearSolver::AMS:
        MFEM_VERIFY(mg_config->h1_aux_fespaces,
                    "AMS coarse solver requires auxiliary H1 space!");
        return std::make_unique<MfemWrapperSolver<ComplexOperator>>(
            std::make_unique<HypreAmsSolver>(
                mg_config->nd_fespaces->GetFESpaceAtLevel(0),
                mg_config->h1_aux_fespaces->GetFESpaceAtLevel(0), linear.ams_max_it,
                linear.mg_smooth_it, linear.ams_vector_interp, linear.ams_singular_op,
                linear.amg_agg_coarsen, print),
            true, linear.complex_coarse_solve, linear.drop_small_entries,
            linear.reorder_reuse);
      case LinearSolver::BOOMER_AMG:
        return std::make_unique<MfemWrapperSolver<ComplexOperator>>(
            std::make_unique<BoomerAmgSolver>(1, linear.mg_smooth_it,
                                              linear.amg_agg_coarsen, print),
            true, linear.complex_coarse_solve, linear.drop_small_entries,
            linear.reorder_reuse);
      case LinearSolver::SUPERLU:
#if defined(MFEM_USE_SUPERLU)
        return std::make_unique<MfemWrapperSolver<ComplexOperator>>(
            std::make_unique<SuperLUSolver>(comm, linear.sym_factorization,
                                            linear.superlu_3d, true, print),
            false, linear.complex_coarse_solve, linear.drop_small_entries,
            linear.reorder_reuse);
#else
        MFEM_ABORT("Solver was not built with SuperLU_DIST support!");
        return {};
#endif
      case LinearSolver::STRUMPACK:
      case LinearSolver::STRUMPACK_MP:
#if defined(MFEM_USE_STRUMPACK)
        return std::make_unique<MfemWrapperSolver<ComplexOperator>>(
            std::make_unique<StrumpackSolver>(
                comm, linear.sym_factorization, linear.strumpack_compression_type,
                linear.strumpack_lr_tol, linear.strumpack_butterfly_l,
                linear.strumpack_lossy_precision, true, print),
            false, linear.complex_coarse_solve, linear.drop_small_entries,
            linear.reorder_reuse);
#else
        MFEM_ABORT("Solver was not built with STRUMPACK support!");
        return {};
#endif
      case LinearSolver::MUMPS:
#if defined(MFEM_USE_MUMPS)
        return std::make_unique<MfemWrapperSolver<ComplexOperator>>(
            std::make_unique<MumpsSolver>(comm, MatrixSymmetry::UNSYMMETRIC,
                                          linear.sym_factorization, linear.strumpack_lr_tol,
                                          true, print),
            false, linear.complex_coarse_solve, linear.drop_small_entries,
            linear.reorder_reuse);
#else
        MFEM_ABORT("Solver was not built with MUMPS support!");
        return {};
#endif
      default:
        MFEM_ABORT("Unsupported coarse solver type for multigrid boundary mode solver!");
        return {};
    }
  };

  // ND block: p-multigrid with Hiptmair distributive relaxation smoothing.
  const auto nd_P = mg_config->nd_fespaces->GetProlongationOperators();
  const auto nd_G =
      mg_config->nd_fespaces->GetDiscreteInterpolators(*mg_config->h1_aux_fespaces);
  auto nd_gmg = std::make_unique<GeometricMultigridSolver<ComplexOperator>>(
      comm, MakeCoarseSolver(nd_pc_type), nd_P, &nd_G, linear.mg_cycle_it,
      linear.mg_smooth_it, linear.mg_smooth_order, linear.mg_smooth_sf_max,
      linear.mg_smooth_sf_min, linear.mg_smooth_cheby_4th);
  nd_gmg->EnableTimer();

  // H1 block: p-multigrid with Chebyshev smoothing.
  const auto h1_P = mg_config->h1_fespaces->GetProlongationOperators();
  auto h1_gmg = std::make_unique<GeometricMultigridSolver<ComplexOperator>>(
      comm, MakeCoarseSolver(h1_pc_type), h1_P, nullptr, linear.mg_cycle_it,
      linear.mg_smooth_it, linear.mg_smooth_order, linear.mg_smooth_sf_max,
      linear.mg_smooth_sf_min, linear.mg_smooth_cheby_4th);
  h1_gmg->EnableTimer();

  // Combine into block-diagonal preconditioner.
  auto block_pc = std::make_unique<BlockDiagonalPreconditioner<ComplexOperator>>(
      nd_size, std::move(nd_gmg), std::move(h1_gmg));
  block_pc_ptr = block_pc.get();

  // Outer Krylov solver — use the user-configured type from Solver.Linear.KSPType.
  std::unique_ptr<IterativeSolver<ComplexOperator>> krylov;
  switch (linear.krylov_solver)
  {
    case KrylovSolver::CG:
      krylov = std::make_unique<CgSolver<ComplexOperator>>(comm, verbose);
      break;
    case KrylovSolver::FGMRES:
      {
        auto fgmres = std::make_unique<FgmresSolver<ComplexOperator>>(comm, verbose);
        fgmres->SetRestartDim(linear.max_size);
        krylov = std::move(fgmres);
      }
      break;
    case KrylovSolver::GMRES:
    case KrylovSolver::DEFAULT:
    default:
      {
        auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(comm, verbose);
        gmres->SetRestartDim(linear.max_size);
        krylov = std::move(gmres);
      }
      break;
  }
  krylov->SetInitialGuess(linear.initial_guess);
  krylov->SetRelTol(linear.tol);
  krylov->SetMaxIter(linear.max_it);
  krylov->EnableTimer();

  ksp = std::make_unique<ComplexKspSolver>(std::move(krylov), std::move(block_pc));
}

std::unique_ptr<ComplexMultigridOperator>
ModeEigenSolver::AssembleAttPreconditioner(double omega, double sigma) const
{
  MFEM_VERIFY(mg_config && mg_config->nd_fespaces && mg_config->h1_aux_fespaces,
              "AssembleAttPreconditioner requires ND and H1 auxiliary hierarchies!");
  const auto n_levels = mg_config->nd_fespaces->GetNumLevels();
  auto B = std::make_unique<ComplexMultigridOperator>(n_levels);

  // Material coefficients (same at all levels — indexed by element attribute, not by p).
  MaterialPropertyCoefficient muinv_cc_func(mat_op.GetAttributeToMaterial(),
                                            normal ? mat_op.GetInvPermeability()
                                                   : mat_op.GetCurlCurlInvPermeability());
  if (normal)
  {
    muinv_cc_func.NormalProjectedCoefficient(*normal);
  }

  // Preconditioner mass coefficient: -omega^2 eps - sigma/mu (+ London). When
  // pc_mat_shifted is enabled (config), use absolute values to ensure well-conditioning
  // near the shift target (where -omega^2 eps + kn^2/mu ≈ 0). This follows
  // SpaceOperator's pc_mat_shifted pattern.
  const bool shifted = linear.pc_mat_shifted;
  const double mass_coeff = shifted ? std::abs(omega * omega) : (-omega * omega);
  const double shift_coeff = shifted ? std::abs(sigma) : (-sigma);

  MaterialPropertyCoefficient eps_shifted_pc(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetPermittivityReal(), mass_coeff);
  eps_shifted_pc.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                mat_op.GetInvPermeability(), shift_coeff);
  if (mat_op.HasLondonDepth())
  {
    eps_shifted_pc.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                  mat_op.GetInvLondonDepth(), 1.0);
  }

  // Boundary coefficients for preconditioner (real part only, matching SpaceOperator).
  int max_bdr_attr = mat_op.MaxCeedBdrAttribute();
  MaterialPropertyCoefficient fbr(max_bdr_attr);
  if (surf_z_op)
  {
    surf_z_op->AddStiffnessBdrCoefficients(1.0, fbr);
    surf_z_op->AddMassBdrCoefficients(shifted ? std::abs(omega * omega) : (-omega * omega),
                                      fbr);
  }
  if (farfield_op)
  {
    farfield_op->AddDampingBdrCoefficients(omega, fbr);
  }
  if (surf_sigma_op)
  {
    surf_sigma_op->AddExtraSystemBdrCoefficients(omega, fbr, fbr);
  }

  // Assemble ND operators at all levels using the hierarchy-aware assembly pattern
  // (matching SpaceOperator::AssembleOperators). Creates a BilinearForm on the finest
  // space and uses Assemble(fespaces) to produce operators at all levels.
  constexpr bool assemble_q_data = false;
  {
    BilinearForm att(mg_config->nd_fespaces->GetFinestFESpace());
    att.AddDomainIntegrator<CurlCurlMassIntegrator>(muinv_cc_func, eps_shifted_pc);
    if (!fbr.empty())
    {
      att.AddBoundaryIntegrator<VectorFEMassIntegrator>(fbr);
    }
    auto att_ops = att.Assemble(*mg_config->nd_fespaces, skip_zeros);

    // Auxiliary H1 operators for Hiptmair smoothing (matching SpaceOperator::
    // AssembleAuxOperators — uses DiffusionIntegrator with the mass coefficient).
    BilinearForm att_aux(mg_config->h1_aux_fespaces->GetFinestFESpace());
    att_aux.AddDomainIntegrator<DiffusionIntegrator>(eps_shifted_pc);
    if (!fbr.empty())
    {
      att_aux.AddBoundaryIntegrator<DiffusionIntegrator>(fbr);
    }
    auto att_aux_ops = att_aux.Assemble(*mg_config->h1_aux_fespaces, skip_zeros);

    for (std::size_t l = 0; l < n_levels; l++)
    {
      const auto &nd_fespace_l = mg_config->nd_fespaces->GetFESpaceAtLevel(l);
      auto B_l = std::make_unique<ComplexParOperator>(std::move(att_ops[l]), nullptr,
                                                      nd_fespace_l);
      if (mg_config->nd_dbc_tdof_lists && l < mg_config->nd_dbc_tdof_lists->size())
      {
        B_l->SetEssentialTrueDofs((*mg_config->nd_dbc_tdof_lists)[l],
                                  Operator::DiagonalPolicy::DIAG_ONE);
      }
      B->AddOperator(std::move(B_l));

      const auto &h1_aux_l = mg_config->h1_aux_fespaces->GetFESpaceAtLevel(l);
      auto B_aux_l = std::make_unique<ComplexParOperator>(std::move(att_aux_ops[l]),
                                                          nullptr, h1_aux_l);
      if (mg_config->h1_aux_dbc_tdof_lists && l < mg_config->h1_aux_dbc_tdof_lists->size())
      {
        B_aux_l->SetEssentialTrueDofs((*mg_config->h1_aux_dbc_tdof_lists)[l],
                                      Operator::DiagonalPolicy::DIAG_ONE);
      }
      B->AddAuxiliaryOperator(std::move(B_aux_l));
    }
  }

  return B;
}

std::unique_ptr<ComplexMultigridOperator>
ModeEigenSolver::AssembleAnnPreconditioner(double omega) const
{
  MFEM_VERIFY(mg_config && mg_config->h1_fespaces,
              "AssembleAnnPreconditioner requires H1 hierarchy!");
  const auto n_levels = mg_config->h1_fespaces->GetNumLevels();
  auto B = std::make_unique<ComplexMultigridOperator>(n_levels);

  // Material coefficients matching AssembleAnn (real part only). The negative diffusion
  // sign from the IBP is preserved — this matches the actual operator.
  MaterialPropertyCoefficient neg_muinv_func(mat_op.GetAttributeToMaterial(),
                                             mat_op.GetInvPermeability(), -1.0);
  if (normal)
  {
    neg_muinv_func.NormalProjectedCoefficient(*normal);
  }

  MaterialPropertyCoefficient poseps_h1_func(mat_op.GetAttributeToMaterial(),
                                             normal ? mat_op.GetPermittivityReal()
                                                    : mat_op.GetPermittivityScalar(),
                                             omega * omega);
  if (normal)
  {
    poseps_h1_func.NormalProjectedCoefficient(*normal);
  }
  if (mat_op.HasLondonDepth())
  {
    if (!normal)
    {
      poseps_h1_func.AddCoefficient(mat_op.GetAttributeToMaterial(),
                                    mat_op.GetInvLondonDepthScalar());
    }
    else if (normal)
    {
      const auto &ild = mat_op.GetInvLondonDepth();
      mfem::DenseTensor ild_scalar(1, 1, ild.SizeK());
      for (int k = 0; k < ild.SizeK(); k++)
      {
        ild_scalar(0, 0, k) = ild(0, 0, k);
      }
      poseps_h1_func.AddCoefficient(mat_op.GetAttributeToMaterial(), ild_scalar);
      poseps_h1_func.NormalProjectedCoefficient(*normal);
    }
  }

  // Boundary coefficients (real part only).
  int max_bdr_attr = mat_op.MaxCeedBdrAttribute();
  MaterialPropertyCoefficient nn_fbr(max_bdr_attr);
  if (surf_z_op)
  {
    surf_z_op->AddStiffnessBdrCoefficients(-1.0, nn_fbr);
    surf_z_op->AddMassBdrCoefficients(omega * omega, nn_fbr);
  }
  if (surf_sigma_op)
  {
    MaterialPropertyCoefficient cond_r(max_bdr_attr);
    surf_sigma_op->AddExtraSystemBdrCoefficients(omega, cond_r, cond_r);
    if (!cond_r.empty())
    {
      cond_r *= -1.0;
      nn_fbr.AddCoefficient(cond_r.GetAttributeToMaterial(),
                            cond_r.GetMaterialProperties());
    }
  }

  // Assemble H1 operators at all levels using hierarchy-aware assembly.
  {
    BilinearForm ann(mg_config->h1_fespaces->GetFinestFESpace());
    ann.AddDomainIntegrator<DiffusionMassIntegrator>(neg_muinv_func, poseps_h1_func);
    if (!nn_fbr.empty())
    {
      ann.AddBoundaryIntegrator<MassIntegrator>(nn_fbr);
    }
    auto ann_ops = ann.Assemble(*mg_config->h1_fespaces, skip_zeros);

    for (std::size_t l = 0; l < n_levels; l++)
    {
      const auto &h1_fespace_l = mg_config->h1_fespaces->GetFESpaceAtLevel(l);
      auto B_l = std::make_unique<ComplexParOperator>(std::move(ann_ops[l]), nullptr,
                                                      h1_fespace_l);
      if (mg_config->h1_dbc_tdof_lists && l < mg_config->h1_dbc_tdof_lists->size())
      {
        B_l->SetEssentialTrueDofs((*mg_config->h1_dbc_tdof_lists)[l],
                                  Operator::DiagonalPolicy::DIAG_ONE);
      }
      B->AddOperator(std::move(B_l));
    }
  }

  return B;
}

void ModeEigenSolver::SetUpEigenSolver(MPI_Comm comm)
{
  constexpr int print = 0;
  EigenSolverBackend type = eigen_backend;

  if (type == EigenSolverBackend::SLEPC)
  {
#if !defined(PALACE_WITH_SLEPC)
    MFEM_ABORT("Solver was not built with SLEPc support, please choose a "
               "different solver!");
#endif
  }
  else if (type == EigenSolverBackend::ARPACK)
  {
#if !defined(PALACE_WITH_ARPACK)
    MFEM_ABORT("Solver was not built with ARPACK support, please choose a "
               "different solver!");
#endif
  }
  else  // Default choice
  {
#if defined(PALACE_WITH_SLEPC)
    type = EigenSolverBackend::SLEPC;
#elif defined(PALACE_WITH_ARPACK)
    type = EigenSolverBackend::ARPACK;
#else
#error "ModeEigenSolver requires building with ARPACK or SLEPc!"
#endif
  }

  if (type == EigenSolverBackend::ARPACK)
  {
#if defined(PALACE_WITH_ARPACK)
    auto arpack = std::make_unique<arpack::ArpackEPSSolver>(comm, print);
    arpack->SetNumModes(num_modes, num_vec);
    arpack->SetTol(eig_tol);
    arpack->SetWhichEigenpairs(which_eig);
    arpack->SetLinearSolver(*ksp);
    eigen = std::move(arpack);
#endif
  }
  else  // EigenSolverBackend::SLEPC
  {
#if defined(PALACE_WITH_SLEPC)
    auto slepc = std::make_unique<slepc::SlepcEPSSolver>(comm, print);
    slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
    slepc->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GEN_NON_HERMITIAN);
    slepc->SetNumModes(num_modes, num_vec);
    slepc->SetTol(eig_tol);
    slepc->SetWhichEigenpairs(which_eig);
    slepc->SetLinearSolver(*ksp);
    eigen = std::move(slepc);
#endif
  }
}

}  // namespace palace
