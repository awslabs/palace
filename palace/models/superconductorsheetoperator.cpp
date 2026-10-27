// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "superconductorsheetoperator.hpp"

#include <cmath>
#include <map>
#include <set>
#include "fem/bilinearform.hpp"
#include "fem/brokenspace.hpp"
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/integrator.hpp"
#include "linalg/rap.hpp"
#include "models/boundaryattributes.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

SuperconductorSheetOperator::SuperconductorSheetOperator(
    const std::vector<config::SuperconductorData> &superconductor,
    const std::unordered_set<int> &cracked_attributes, const Units &units,
    const MaterialOperator &mat_op, const mfem::ParMesh &mesh)
  : mat_op(mat_op)
{
  SetUpBoundaryProperties(superconductor, cracked_attributes, mesh);
  PrintBoundaryInfo(units, mesh);
}

SuperconductorSheetOperator::SuperconductorSheetOperator(const IoData &iodata,
                                                         const MaterialOperator &mat_op,
                                                         const mfem::ParMesh &mesh)
  : SuperconductorSheetOperator(iodata.boundaries.superconductor,
                                iodata.boundaries.cracked_attributes, iodata.units, mat_op,
                                mesh)
{
}

void SuperconductorSheetOperator::SetUpBoundaryProperties(
    const std::vector<config::SuperconductorData> &superconductor,
    const std::unordered_set<int> &cracked_attributes, const mfem::ParMesh &mesh)
{
  // Check that superconductor sheet boundary attributes have been specified correctly.
  MeshBoundaryAttributes mesh_attrs(mesh);
  if (!superconductor.empty())
  {
    CheckBoundaryAttributes(mesh_attrs, superconductor, "superconductor sheet");
  }

  // Kinetic sheet inductance L_ksq [H/sq]: supplied directly, or from (lambda, d) via the
  // finite-thickness London sheet L_ksq = lambda*coth(d/lambda) (nondimensional, mu0
  // absorbed). Reduces to the thin-film Pearl limit lambda^2/d for d << lambda and
  // saturates at lambda for d >> lambda, capturing the through-thickness current profile
  // without resolving the film thickness.
  boundaries.reserve(superconductor.size());
  bool two_sided = false;
  for (const auto &data : superconductor)
  {
    // Two-sided sheet: the diagonal self-term inductance of each side is λ·tanh(d/λ), so
    // 1/Ls = coth(d/λ)/λ, and the cross term is assembled with it in
    // AssembleBrokenSheetMatrix. Otherwise the single-sheet L_ksq applies.
    const double Ls =
        data.two_sided
            ? data.lambda_L * std::tanh(data.thickness / data.lambda_L)
            : ((data.Ls > 0.0) ? data.Ls
                               : KineticSheetInductance(data.lambda_L, data.thickness));
    MFEM_VERIFY(Ls > 0.0,
                "Superconductor sheet has non-positive kinetic sheet inductance!");
    auto &bdr = boundaries.emplace_back();
    bdr.Ls = Ls;
    bdr.lambda_L = data.lambda_L;
    bdr.thickness = data.thickness;
    bdr.two_sided = data.two_sided;
    two_sided = two_sided || data.two_sided;
    bdr.attr_list.Reserve(static_cast<int>(data.attributes.size()));
    for (auto attr : data.attributes)
    {
      if (!mesh_attrs.Contains(attr))
      {
        continue;  // Can just ignore if wrong
      }
      bdr.attr_list.Append(attr);
      // Per-attribute scaling to account for increased area when using mesh cracking.
      bdr.attr_scaling[attr] = cracked_attributes.contains(attr) ? 2.0 : 1.0;
    }
  }
  if (!two_sided)
  {
    return;
  }
  MFEM_VERIFY(mesh.Dimension() == 3,
              "Superconductor boundary \"TwoSided\" requires a 3D mesh!");

  // The two sides of a two-sided film are the elements on either side of its faces, so
  // the film must be interior to the mesh.
  int num_exterior = 0;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    const int attr = mesh.GetBdrAttribute(be);
    for (const auto &bdr : boundaries)
    {
      if (bdr.two_sided && bdr.attr_list.Find(attr) >= 0)
      {
        num_exterior +=
            mesh.GetFaceInformation(mesh.GetBdrElementFaceIndex(be)).IsBoundary();
      }
    }
  }
  Mpi::GlobalSum(1, &num_exterior, mesh.GetComm());
  MFEM_VERIFY(num_exterior == 0, "Superconductor boundary \"TwoSided\" has "
                                     << num_exterior
                                     << " face(s) on the exterior boundary of the mesh; "
                                        "the film must be interior!");
}

void SuperconductorSheetOperator::PrintBoundaryInfo(const Units &units,
                                                    const mfem::ParMesh &mesh)
{
  if (boundaries.empty())
  {
    return;
  }

  fmt::memory_buffer buffer{};
  auto out = fmt::appender{buffer};
  using VT = Units::ValueType;

  fmt::format_to(out, "\nConfiguring superconductor sheet BC at attributes:\n");
  for (const auto &bdr : boundaries)
  {
    for (auto attr : bdr.attr_list)
    {
      if (bdr.two_sided)
      {
        // Each side's self-term, and the single sheet ½λ·coth(d/2λ) the pair reduces to
        // when both sides carry the same field.
        fmt::format_to(
            out,
            " {:d}: L_ksq = {:.3e} H/sq per side, {:.3e} H/sq common mode, n = ({:+.1f})\n",
            attr, units.Dimensionalize<VT::INDUCTANCE>(bdr.Ls),
            units.Dimensionalize<VT::INDUCTANCE>(
                0.5 * KineticSheetInductance(bdr.lambda_L, 0.5 * bdr.thickness)),
            fmt::join(mesh::GetSurfaceNormal(mesh, attr), ","));
        continue;
      }
      fmt::format_to(out, " {:d}: L_ksq = {:.3e} H/sq, n = ({:+.1f})\n", attr,
                     units.Dimensionalize<VT::INDUCTANCE>(bdr.Ls),
                     fmt::join(mesh::GetSurfaceNormal(mesh, attr), ","));
    }
  }
  Mpi::Print("{}", fmt::to_string(buffer));
}

double SuperconductorSheetOperator::KineticSheetInductance(double lambda_L,
                                                           double thickness)
{
  // Finite-thickness London sheet inductance L_ksq = lambda * coth(d/lambda), with
  // coth(x) = 1/tanh(x). Reduces to lambda^2/d for d << lambda; saturates at lambda for
  // d >> lambda.
  return lambda_L / std::tanh(thickness / lambda_L);
}

mfem::Array<int> SuperconductorSheetOperator::GetAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &bdr : boundaries)
  {
    attr_list.Append(bdr.attr_list);
  }
  return attr_list;
}

void SuperconductorSheetOperator::AddStiffnessBdrCoefficients(
    double coeff, MaterialPropertyCoefficient &fb) const
{
  // Kinetic sheet inductance boundaries: add (coeff / L_ksq) as a tangential surface mass.
  for (const auto &bdr : boundaries)
  {
    if (bdr.two_sided)
    {
      continue;
    }
    for (auto attr : bdr.attr_list)
    {
      const double s = bdr.attr_scaling.at(attr);
      fb.AddMaterialProperty(mat_op.GetCeedBdrAttributes(attr), coeff / (bdr.Ls * s));
    }
  }
}

std::vector<int> SuperconductorSheetOperator::GetTwoSidedAttributes() const
{
  std::vector<int> attr_list;
  for (const auto &bdr : boundaries)
  {
    if (bdr.two_sided)
    {
      attr_list.insert(attr_list.end(), bdr.attr_list.begin(), bdr.attr_list.end());
    }
  }
  return attr_list;
}

namespace
{

// Boundary attribute of each mesh face carrying one of attrs (0 otherwise), also for faces
// whose boundary element is on another process or is a slave of a nonconforming master,
// propagated through the true DOFs of the lowest-order H(div) space.
std::vector<int> GetFaceAttributes(const mfem::ParMesh &cmesh, const std::set<int> &attrs)
{
  auto &mesh = const_cast<mfem::ParMesh &>(cmesh);
  mfem::RT_FECollection fec(0, mesh.Dimension());
  mfem::ParFiniteElementSpace fespace(&mesh, &fec);
  mfem::Vector num(fespace.GetVSize()), den(fespace.GetVSize());
  num = 0.0;
  den = 0.0;
  mfem::Array<int> dofs;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    const int attr = mesh.GetBdrAttribute(be);
    if (attrs.contains(attr))
    {
      fespace.GetFaceDofs(mesh.GetBdrElementFaceIndex(be), dofs);
      const int d = (dofs[0] >= 0) ? dofs[0] : -1 - dofs[0];
      num(d) = attr;
      den(d) = 1.0;
    }
  }
  mfem::HypreParMatrix P(*fespace.Dof_TrueDof_Matrix());
  P.HostReadWrite();
  {
    hypre_ParCSRMatrix *hP = P;
    for (auto *b : {hypre_ParCSRMatrixDiag(hP), hypre_ParCSRMatrixOffd(hP)})
    {
      auto *data = hypre_CSRMatrixData(b);
      for (int k = 0; data && k < hypre_CSRMatrixNumNonzeros(b); k++)
      {
        data[k] = std::abs(data[k]);
      }
    }
  }
  mfem::Vector t(P.Width());
  P.MultTranspose(num, t);
  P.Mult(t, num);
  P.MultTranspose(den, t);
  P.Mult(t, den);
  std::vector<int> face_attr(mesh.GetNumFaces(), 0);
  for (int f = 0; f < mesh.GetNumFaces(); f++)
  {
    fespace.GetFaceDofs(f, dofs);
    const int d = (dofs[0] >= 0) ? dofs[0] : -1 - dofs[0];
    face_attr[f] = (den(d) > 0.0) ? static_cast<int>(std::lround(num(d) / den(d))) : 0;
  }
  return face_attr;
}

// Add the tangential surface mass ∫ coeff u_t·v_t over the local face i of element e, with
// the shape functions of the element, to elmat (in the reference element DOFs).
void AddElementFaceMass(const mfem::ParMesh &mesh, int e, int i, double coeff,
                        const mfem::FiniteElement &fe, mfem::IsoparametricTransformation &T,
                        mfem::DenseMatrix &elmat)
{
  const mfem::Element &el = *mesh.GetElement(e);
  const int nfv = el.GetNFaceVertices(i);
  const int *fv = el.GetFaceVertices(i);
  const auto face_geom = (nfv == 3) ? mfem::Geometry::TRIANGLE : mfem::Geometry::SQUARE;
  const mfem::IntegrationRule &verts = *mfem::Geometries.GetVertices(el.GetGeometryType());
  const mfem::IntegrationRule &ir =
      mfem::IntRules.Get(face_geom, fem::DefaultIntegrationOrder::Get(mesh, face_geom));
  double V[4][3];
  for (int k = 0; k < nfv; k++)
  {
    const auto &v = verts.IntPoint(fv[k]);
    V[k][0] = v.x;
    V[k][1] = v.y;
    V[k][2] = v.z;
  }
  const int ndof = fe.GetDof();
  mfem::DenseMatrix vshape(ndof, 3);
  mfem::Vector un(ndof);
  for (int q = 0; q < ir.GetNPoints(); q++)
  {
    // Point on the face in the reference element, and the tangent vectors of the face
    // parametrization.
    const mfem::IntegrationPoint &ipf = ir.IntPoint(q);
    const double x = ipf.x, y = ipf.y;
    double X[3], dx[3], dy[3];
    for (int c = 0; c < 3; c++)
    {
      if (nfv == 3)
      {
        X[c] = V[0][c] + x * (V[1][c] - V[0][c]) + y * (V[2][c] - V[0][c]);
        dx[c] = V[1][c] - V[0][c];
        dy[c] = V[2][c] - V[0][c];
      }
      else
      {
        X[c] = (1.0 - x) * (1.0 - y) * V[0][c] + x * (1.0 - y) * V[1][c] + x * y * V[2][c] +
               (1.0 - x) * y * V[3][c];
        dx[c] = (1.0 - y) * (V[1][c] - V[0][c]) + y * (V[2][c] - V[3][c]);
        dy[c] = (1.0 - x) * (V[3][c] - V[0][c]) + x * (V[2][c] - V[1][c]);
      }
    }
    mfem::IntegrationPoint ip;
    ip.Set3(X[0], X[1], X[2]);
    T.SetIntPoint(&ip);
    const mfem::DenseMatrix &J = T.Jacobian();
    double t1[3], t2[3], n[3];
    for (int r = 0; r < 3; r++)
    {
      t1[r] = J(r, 0) * dx[0] + J(r, 1) * dx[1] + J(r, 2) * dx[2];
      t2[r] = J(r, 0) * dy[0] + J(r, 1) * dy[1] + J(r, 2) * dy[2];
    }
    n[0] = t1[1] * t2[2] - t1[2] * t2[1];
    n[1] = t1[2] * t2[0] - t1[0] * t2[2];
    n[2] = t1[0] * t2[1] - t1[1] * t2[0];
    const double dA = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    for (auto &c : n)
    {
      c /= dA;
    }
    fe.CalcVShape(T, vshape);
    for (int a = 0; a < ndof; a++)
    {
      un(a) = vshape(a, 0) * n[0] + vshape(a, 1) * n[1] + vshape(a, 2) * n[2];
    }
    const double w = coeff * ipf.weight * dA;
    for (int b = 0; b < ndof; b++)
    {
      for (int a = 0; a < ndof; a++)
      {
        elmat(a, b) += w * (vshape(a, 0) * vshape(b, 0) + vshape(a, 1) * vshape(b, 1) +
                            vshape(a, 2) * vshape(b, 2) - un(a) * un(b));
      }
    }
  }
}

}  // namespace

std::unique_ptr<mfem::HypreParMatrix>
SuperconductorSheetOperator::AssembleBrokenSheetMatrix(
    const FiniteElementSpace &nd_fespace, const FiniteElementSpace &broken_nd_fespace,
    const CrackSides &sides) const
{
  MFEM_VERIFY(broken_nd_fespace.IsBroken() && &broken_nd_fespace.Get() == &nd_fespace.Get(),
              "Two-sided superconductor sheets need a view of the ND space broken across "
              "the films!");
  const auto &mesh = nd_fespace.GetParMesh();
  MPI_Comm comm = nd_fespace.GetComm();

  // Coefficients a + b (each side) and b (sum of the sides) per film attribute.
  std::map<int, std::pair<double, double>> coeffs;
  std::set<int> attrs;
  for (const auto &bdr : boundaries)
  {
    if (!bdr.two_sided)
    {
      continue;
    }
    const double D = bdr.thickness / bdr.lambda_L;
    const double a = 1.0 / (bdr.lambda_L * std::tanh(D)),
                 b = 1.0 / (bdr.lambda_L * std::sinh(D));
    for (auto attr : bdr.attr_list)
    {
      coeffs[attr] = {a + b, b};
      attrs.insert(attr);
    }
  }

  // The DOFs of the two sides of a film face are versions of one another: an entity at a
  // junction of films would have more than two.
  {
    int num_junction = 0;
    for (int t = 0; t < nd_fespace.GetTrueVSize(); t++)
    {
      num_junction += (broken_nd_fespace.GetNumVersions(t) > 2);
    }
    Mpi::GlobalSum(1, &num_junction, comm);
    MFEM_VERIFY(num_junction == 0,
                "Superconductor boundaries \"TwoSided\" cannot meet each other ("
                    << num_junction << " DOF(s) with more than two sides)!");
  }

  // Each side: (a + b) ∫ u·v over the film faces of each element, on the L-vector of the
  // broken space, lifted with its prolongation, together with the one-sided sheets.
  std::unique_ptr<mfem::HypreParMatrix> C;
  {
    const auto face_attr = GetFaceAttributes(mesh, attrs);
    auto C_l = std::make_unique<mfem::SparseMatrix>(broken_nd_fespace.GetVSize());
    mfem::Array<int> faces, orient, dofs;
    mfem::DofTransformation dof_trans;
    mfem::IsoparametricTransformation T;
    mfem::DenseMatrix elmat;
    int num_unknown = 0;
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      if (!sides.split[e])
      {
        continue;
      }
      const auto geom = mesh.GetElementGeometry(e);
      const int face_bit0 = mfem::Geometry::NumVerts[geom] + mfem::Geometry::NumEdges[geom];
      mesh.GetElementFaces(e, faces, orient);
      const mfem::FiniteElement &fe = *nd_fespace.Get().GetFE(e);
      elmat.SetSize(fe.GetDof());
      elmat = 0.0;
      bool any = false;
      for (int i = 0; i < faces.Size(); i++)
      {
        if (!((sides.split[e] >> (face_bit0 + i)) & 1u))
        {
          continue;
        }
        auto it = coeffs.find(face_attr[faces[i]]);
        if (it == coeffs.end())
        {
          num_unknown++;
          continue;
        }
        mesh.GetElementTransformation(e, &T);
        AddElementFaceMass(mesh, e, i, it->second.first, fe, T, elmat);
        any = true;
      }
      if (any)
      {
        broken_nd_fespace.GetElementDofs(e, dofs, dof_trans);
        if (!dof_trans.IsIdentity())
        {
          dof_trans.TransformDual(elmat);
        }
        C_l->AddSubMatrix(dofs, dofs, elmat, 0);
      }
    }
    Mpi::GlobalSum(1, &num_unknown, comm);
    MFEM_VERIFY(num_unknown == 0, num_unknown << " split face(s) of a broken space with no "
                                                 "two-sided superconductor attribute!");

    // One-sided sheets: each boundary element is integrated by its adjacent element with
    // the DOFs of its side, so that a sheet standing on one side of a two-sided film reads
    // the versions of that side (both sides read a sheet away from the films).
    std::map<int, double> one_sided;
    for (const auto &bdr : boundaries)
    {
      if (!bdr.two_sided)
      {
        for (auto attr : bdr.attr_list)
        {
          one_sided[attr] = 1.0 / (bdr.Ls * bdr.attr_scaling.at(attr));
        }
      }
    }
    for (int be = 0; be < mesh.GetNBE() && !one_sided.empty(); be++)
    {
      auto it = one_sided.find(mesh.GetBdrAttribute(be));
      if (it == one_sided.end())
      {
        continue;
      }
      int e, info;
      mesh.GetBdrElementAdjacentElement(be, e, info);
      const mfem::FiniteElement &fe = *nd_fespace.Get().GetFE(e);
      elmat.SetSize(fe.GetDof());
      elmat = 0.0;
      mesh.GetElementTransformation(e, &T);
      AddElementFaceMass(mesh, e, info / 64, it->second, fe, T, elmat);
      broken_nd_fespace.GetElementDofs(e, dofs, dof_trans);
      if (!dof_trans.IsIdentity())
      {
        dof_trans.TransformDual(elmat);
      }
      C_l->AddSubMatrix(dofs, dofs, elmat, 0);
    }
    auto C_h = linalg::BlockDiagonalParMatrix(
        comm, broken_nd_fespace.GlobalVSize(), broken_nd_fespace.GlobalVSize(),
        broken_nd_fespace.GetDofOffsets(), broken_nd_fespace.GetDofOffsets(),
        std::move(C_l));
    C = linalg::RAP(*C_h, *broken_nd_fespace.Dof_TrueDof_Matrix());
  }

  // Sum of the sides: −b ∫ (u⁺ + u⁻)·(v⁺ + v⁻) = −b Sᵀ M S, with M the film boundary mass
  // of the unbroken space.
  {
    MaterialPropertyCoefficient fb(mat_op.MaxCeedBdrAttribute());
    for (const auto &[attr, ab] : coeffs)
    {
      fb.AddMaterialProperty(mat_op.GetCeedBdrAttributes(attr), ab.second);
    }
    BilinearForm m(nd_fespace);
    if (!fb.empty())  // A process may have no film boundary elements.
    {
      m.AddBoundaryIntegrator<VectorFEMassIntegrator>(fb);
    }
    auto M = ParOperator(m.FullAssemble(false), nd_fespace).StealParallelAssemble();
    auto S = broken_nd_fespace.BuildVersionMap(FiniteElementSpace::VersionMap::SIDES);
    auto X = linalg::RAP(*M, *S);
    C.reset(mfem::Add(1.0, *C, -1.0, *X));
  }
  return C;
}

}  // namespace palace
