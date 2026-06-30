// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "interpolator.hpp"

#include <deque>
#include <limits>
#include <memory>

#include <algorithm>
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

#include "fem/libceed/basis.hpp"
#include "fem/libceed/ceed.hpp"
#include "fem/libceed/functional.hpp"
#include "fem/libceed/restriction.hpp"
#include "fem/output_functionals.hpp"
#include "utils/diagnostic.hpp"

PalacePragmaDiagnosticPush
PalacePragmaDiagnosticDisableUnused

#include "fem/qfunctions/33/eval_33_qf.h"

PalacePragmaDiagnosticPop

namespace palace
{

namespace
{

constexpr auto GSLIB_BB_TOL = 0.01;  // MFEM defaults, slightly reduced bounding box
constexpr auto GSLIB_NEWTON_TOL = 1.0e-12;

}  // namespace

#if defined(MFEM_USE_GSLIB)

// Evaluates fields at the (fixed, located) probe points through libCEED: one small
// operator per locally owned point (volume element restriction and basis tabulated at
// the point, geometry on the fly from mesh node gradients). Results for all points are
// reduced over all processes (unowned and not-found points contribute zero, matching
// the GSLIB default interpolation value).
class CeedProbeEvaluator
{
private:
  std::vector<fem::CeedGroupOperator> groups;
  std::deque<std::unique_ptr<mfem::IntegrationRule>> point_irs;  // Tabulation lifetime
  mutable Vector field_staging, local_out;
  MPI_Comm comm;
  bool valid = true;

public:
  CeedProbeEvaluator(const mfem::FindPointsGSLIB &op,
                     const mfem::ParFiniteElementSpace &fespace)
    : comm(fespace.GetComm())
  {
    const auto &mesh = *fespace.GetParMesh();
    const auto map_type = fespace.FEColl()->GetMapType(mesh.Dimension());
    if (mesh.Dimension() != 3 || mesh.SpaceDimension() != 3 ||
        (map_type != mfem::FiniteElement::H_CURL && map_type != mfem::FiniteElement::H_DIV))
    {
      valid = false;
      return;
    }
    const mfem::FiniteElementSpace &mesh_fespace = *mesh.GetNodes()->FESpace();
    const int npts = op.GetCode().Size();
    const int rank = Mpi::Rank(comm);
    field_staging.SetSize(fespace.GetVSize());
    field_staging.UseDevice(true);
    field_staging = 0.0;
    local_out.SetSize(npts * 3);
    local_out.UseDevice(true);

    // Claim each point on exactly one process: for points on element borders the GSLIB
    // ownership may not be consistent across the (all) searching processes, so resolve
    // ties globally with a minimum reduction.
    std::vector<int> owner(npts);
    for (int i = 0; i < npts; i++)
    {
      owner[i] =
          (op.GetCode()[i] != 2 && op.GetProc()[i] == static_cast<unsigned int>(rank))
              ? rank
              : std::numeric_limits<int>::max();
    }
    Mpi::GlobalMin(npts, owner.data(), comm);

    Ceed ceed = ceed::internal::GetCeedObjects()[0];
    for (int i = 0; i < npts; i++)
    {
      if (owner[i] != rank)
      {
        continue;  // Not found (default value 0) or owned by another process
      }
      const int elem = op.GetElem()[i];
      const auto geom = mesh.GetElementGeometry(elem);
      const int dim = mesh.Dimension();

      // FindPointsGSLIB::GetReferencePosition returns mfem reference coordinates of
      // the original mesh element (point-major).
      auto ir = std::make_unique<mfem::IntegrationRule>(1);
      mfem::IntegrationPoint &ip = ir->IntPoint(0);
      ip.Set3(op.GetReferencePosition()[i * dim + 0],
              op.GetReferencePosition()[i * dim + 1],
              op.GetReferencePosition()[i * dim + 2]);
      ip.weight = 1.0;
      const std::vector<int> indices = {elem};

      std::vector<ceed::CeedFunctionalFieldInput> inputs;
      std::vector<std::pair<std::string, int>> field_sources;
      {
        CeedElemRestriction mesh_restr = FiniteElementSpace::BuildCeedElemRestriction(
            mesh_fespace, ceed, geom, indices, /*is_interp*/ true);
        const mfem::FiniteElement *mesh_fe =
            mesh_fespace.FEColl()->FiniteElementForGeometry(geom);
        CeedBasis mesh_basis;
        ceed::InitBasisAtPoints(*mesh_fe, *ir, mesh_fespace.GetVDim(), ceed, &mesh_basis);
        CeedVector mesh_nodes_vec;
        ceed::InitCeedVector(*mesh.GetNodes(), ceed, &mesh_nodes_vec);
        inputs.push_back(
            {"x", mesh_nodes_vec, mesh_restr, mesh_basis, ceed::EvalMode::Grad});
        CeedElemRestriction field_restr;
        ceed::InitRestriction(fespace, indices, false, /*is_interp*/ true, false, ceed,
                              &field_restr);
        const mfem::FiniteElement *fe = fespace.FEColl()->FiniteElementForGeometry(geom);
        CeedBasis field_basis;
        ceed::InitBasisAtPoints(*fe, *ir, fespace.GetVDim(), ceed, &field_basis);
        CeedVector field_vec;
        ceed::InitCeedVector(field_staging, ceed, &field_vec);
        inputs.push_back(
            {"u_1", field_vec, field_restr, field_basis, ceed::EvalMode::Interp});
        field_sources.emplace_back("u_1", 0);

        // Output: 3 components for point i (byVDIM layout).
        CeedElemRestriction out_restr;
        const CeedInt offset[1] = {3 * i};
        PalaceCeedCall(ceed, CeedElemRestrictionCreate(ceed, 1, 1, 3, 1, (CeedSize)npts * 3,
                                                       CEED_MEM_HOST, CEED_COPY_VALUES,
                                                       offset, &out_restr));

        ceed::CeedQFunctionInfo info;
        if (map_type == mfem::FiniteElement::H_CURL)
        {
          info.apply_qf = f_eval_probe_hcurl_33;
          info.apply_qf_path = PalaceQFunctionRelativePath(f_eval_probe_hcurl_33_loc);
        }
        else
        {
          info.apply_qf = f_eval_probe_hdiv_33;
          info.apply_qf_path = PalaceQFunctionRelativePath(f_eval_probe_hdiv_33_loc);
        }
        CeedOperator point_op;
        ceed::AssembleCeedPointEvaluator(info, nullptr, 0, ceed, inputs, 3, out_restr,
                                         &point_op);
        groups.push_back({ceed, point_op, std::move(field_sources)});

        PalaceCeedCall(ceed, CeedVectorDestroy(&mesh_nodes_vec));
        PalaceCeedCall(ceed, CeedVectorDestroy(&field_vec));
        PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&mesh_restr));
        PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&field_restr));
        PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&out_restr));
        PalaceCeedCall(ceed, CeedBasisDestroy(&mesh_basis));
        PalaceCeedCall(ceed, CeedBasisDestroy(&field_basis));
      }
      point_irs.push_back(std::move(ir));
    }
  }

  ~CeedProbeEvaluator()
  {
    for (auto &group : groups)
    {
      PalaceCeedCall(group.ceed, CeedOperatorDestroy(&group.op));
      if (group.out_vec)
      {
        PalaceCeedCall(group.ceed, CeedVectorDestroy(&group.out_vec));
      }
    }
  }

  bool IsValid() const { return valid; }

  // Evaluate the field at all probe points (byVDIM ordering). Collective.
  std::vector<double> Eval(const mfem::ParGridFunction &U) const
  {
    local_out = 0.0;
    fem::ApplyAddGroupOperators(groups, {&U}, local_out);
    std::vector<double> vals(local_out.Size());
    const double *d = local_out.HostRead();
    std::copy(d, d + local_out.Size(), vals.begin());
    Mpi::GlobalSum(static_cast<int>(vals.size()), vals.data(), comm);
    return vals;
  }
};

#endif

InterpolationOperator::InterpolationOperator(const std::map<int, config::ProbeData> &probe,
                                             const Units &units,
                                             FiniteElementSpace &nd_space)
#if defined(MFEM_USE_GSLIB)
  : op(nd_space.GetParMesh().GetComm()), v_dim_fes(nd_space.Get().GetVectorDim())
{
  auto &mesh = nd_space.GetParMesh();
  // Set up probes interpolation. All processes search for all points.
  if (probe.empty())
  {
    return;
  }
  const int dim = mesh.SpaceDimension();
  MFEM_VERIFY(
      mesh.Dimension() == dim,
      "Probe postprocessing functionality requires mesh dimension == space dimension!");
  const int npts = static_cast<int>(probe.size());
  mfem::Vector xyz(npts * dim);
  op_idx.resize(npts);
  int i = 0;
  for (const auto &[idx, data] : probe)
  {
    for (int d = 0; d < dim; d++)
    {
      // Use default ordering byNODES.
      xyz(d * npts + i) = data.center[d];
    }
    op_idx[i++] = idx;
  }
  op.Setup(mesh, GSLIB_BB_TOL, GSLIB_NEWTON_TOL, npts);
  op.FindPoints(xyz, mfem::Ordering::byNODES);
  op.SetDefaultInterpolationValue(0.0);
  i = 0;
  for (const auto &[idx, data] : probe)
  {
    if (op.GetCode()[i++] == 2)
    {
      Mpi::Warning(
          "Probe {:d} at ({:.3e}) m could not be found!\n Using default value 0.0!\n", idx,
          fmt::join(units.Dimensionalize<Units::ValueType::LENGTH>(data.center), ", "));
    }
  }
}
#else
{
  MFEM_CONTRACT_VAR(GSLIB_BB_TOL);
  MFEM_CONTRACT_VAR(GSLIB_NEWTON_TOL);
  MFEM_VERIFY(probe.empty(), "InterpolationOperator class requires MFEM_USE_GSLIB!");
}
#endif

InterpolationOperator::InterpolationOperator(const IoData &iodata,
                                             FiniteElementSpace &nd_space)
  : InterpolationOperator(iodata.domains.postpro.probe, iodata.units, nd_space)
{
}

InterpolationOperator::~InterpolationOperator() = default;

std::vector<double> InterpolationOperator::ProbeField(const mfem::ParGridFunction &U)
{
#if defined(MFEM_USE_GSLIB)
  const int npts = op.GetCode().Size();

  // Use the libCEED evaluation path when supported and enough probe points are present
  // to amortize operator assembly/JIT. Tiny probe sets are cheaper through the existing
  // GSLIB interpolation path, especially in driven postprocessing where only a few
  // fixed monitor points are evaluated once or twice.
  constexpr int ceed_probe_min_points = 8;
  if (SurfaceFunctional::Enabled() && npts >= ceed_probe_min_points)
  {
    auto &eval = ceed_probes[U.FESpace()];
    if (!eval)
    {
      eval = std::make_unique<CeedProbeEvaluator>(op, *U.ParFESpace());
    }
    if (eval->IsValid())
    {
      return eval->Eval(U);
    }
  }

  // Interpolated vector values are returned from GSLIB interpolator with the same ordering
  // as the source grid function, which we transform to byVDIM for output.
  const int vdim = U.VectorDim();
  std::vector<double> vals(npts * vdim);
  if (U.FESpace()->GetOrdering() == mfem::Ordering::byVDIM)
  {
    mfem::Vector v(vals.data(), npts * vdim);
    op.Interpolate(U, v);
  }
  else
  {
    mfem::Vector v(npts * vdim);
    op.Interpolate(U, v);
    for (int d = 0; d < vdim; d++)
    {
      for (int i = 0; i < npts; i++)
      {
        vals[i * vdim + d] = v(d * npts + i);
      }
    }
  }
  return vals;
#else
  MFEM_ABORT("InterpolationOperator class requires MFEM_USE_GSLIB!");
  return {};
#endif
}

std::vector<std::complex<double>> InterpolationOperator::ProbeField(const GridFunction &U)
{
  std::vector<double> vr = ProbeField(U.Real());
  if (U.HasImag())
  {
    std::vector<double> vi = ProbeField(U.Imag());
    std::vector<std::complex<double>> vals(vr.size());
    std::transform(vr.begin(), vr.end(), vi.begin(), vals.begin(),
                   [](double xr, double xi) { return std::complex<double>(xr, xi); });
    return vals;
  }
  else
  {
    return {vr.begin(), vr.end()};
  }
}

namespace fem
{

void InterpolateFunction(const mfem::GridFunction &U, mfem::GridFunction &V)
{
#if defined(MFEM_USE_GSLIB)
  // Generate list of points where the grid function will be evaluated. If the grid function
  // to interpolate is an H1 space of the same order as the mesh nodes, we can use the
  // mesh node points directly. Otherwise, for a different basis order or type, we generate
  // the interpolation points in the physical space manually.
  auto &dest_mesh = *V.FESpace()->GetMesh();
  MFEM_VERIFY(dest_mesh.GetNodes(), "Destination mesh has no nodal FE space!");
  const int dim = dest_mesh.SpaceDimension();
  mfem::Vector xyz;
  mfem::Ordering::Type ordering;
  const auto *dest_fec_h1 =
      dynamic_cast<const mfem::H1_FECollection *>(V.FESpace()->FEColl());
  const auto *dest_nodes_h1 = dynamic_cast<const mfem::H1_FECollection *>(
      dest_mesh.GetNodes()->FESpace()->FEColl());
  int dest_fespace_order = V.FESpace()->GetMaxElementOrder();
  int dest_nodes_order = dest_mesh.GetNodes()->FESpace()->GetMaxElementOrder();
  dest_mesh.GetNodes()->HostRead();
  if (dest_fec_h1 && dest_nodes_h1 && dest_fespace_order == dest_nodes_order)
  {
    xyz.MakeRef(*dest_mesh.GetNodes(), 0, dest_mesh.GetNodes()->Size());
    ordering = dest_mesh.GetNodes()->FESpace()->GetOrdering();
  }
  else
  {
    int npts = 0, offset = 0;
    for (int i = 0; i < dest_mesh.GetNE(); i++)
    {
      npts += V.FESpace()->GetFE(i)->GetNodes().GetNPoints();
    }
    xyz.SetSize(npts * dim);
    mfem::DenseMatrix pointmat;
    for (int i = 0; i < dest_mesh.GetNE(); i++)
    {
      const mfem::FiniteElement &fe = *V.FESpace()->GetFE(i);
      mfem::ElementTransformation &T = *dest_mesh.GetElementTransformation(i);
      T.Transform(fe.GetNodes(), pointmat);
      for (int j = 0; j < pointmat.Width(); j++)
      {
        for (int d = 0; d < dim; d++)
        {
          // Use default ordering byNODES.
          xyz(d * npts + offset + j) = pointmat(d, j);
        }
      }
      offset += pointmat.Width();
    }
    ordering = mfem::Ordering::byNODES;
  }
  const int npts = xyz.Size() / dim;

  // Set up the interpolator.
  auto &src_mesh = *U.FESpace()->GetMesh();
  MFEM_VERIFY(src_mesh.GetNodes(), "Source mesh has no nodal FE space!");
  auto *src_pmesh = dynamic_cast<mfem::ParMesh *>(&src_mesh);
  MPI_Comm comm = (src_pmesh) ? src_pmesh->GetComm() : MPI_COMM_SELF;
  mfem::FindPointsGSLIB op(comm);
  op.Setup(src_mesh, GSLIB_BB_TOL, GSLIB_NEWTON_TOL, npts);

  // Perform the interpolation and fill the target GridFunction (see MFEM's field-interp
  // miniapp).
  const int vdim = U.VectorDim();
  mfem::Vector vals(npts * vdim);
  op.SetDefaultInterpolationValue(0.0);
  op.SetL2AvgType(mfem::FindPointsGSLIB::NONE);
  op.Interpolate(xyz, U, vals, ordering);
  const auto *dest_fec_l2 =
      dynamic_cast<const mfem::L2_FECollection *>(V.FESpace()->FEColl());
  if (dest_fec_h1 || dest_fec_l2)
  {
    if (dest_fec_h1 && dest_fespace_order != dest_nodes_order)
    {
      // H1 with order != mesh order needs to handle duplicated interpolation points.
      mfem::Vector elem_vals;
      mfem::Array<int> vdofs;
      int offset = 0;
      for (int i = 0; i < dest_mesh.GetNE(); i++)
      {
        const mfem::FiniteElement &fe = *V.FESpace()->GetFE(i);
        const int elem_npts = fe.GetNodes().GetNPoints();
        elem_vals.SetSize(elem_npts * vdim);
        for (int d = 0; d < vdim; d++)
        {
          for (int j = 0; j < elem_npts; j++)
          {
            // Arrange element values byNODES to align with GetElementVDofs.
            int idx = (U.FESpace()->GetOrdering() == mfem::Ordering::byNODES)
                          ? d * npts + offset + j
                          : (offset + j) * vdim + d;
            elem_vals(d * elem_npts + j) = vals(idx);
          }
        }
        const auto *dof_trans = V.FESpace()->GetElementVDofs(i, vdofs);
        if (dof_trans)
        {
          dof_trans->TransformPrimal(elem_vals);
        }
        V.SetSubVector(vdofs, elem_vals);
        offset += elem_npts;
      }
    }
    else
    {
      // Otherwise, H1 and L2 copy interpolated values to vdofs.
      MFEM_ASSERT(V.Size() == vals.Size(),
                  "Unexpected size mismatch for interpolated values and grid function!");
      V = vals;
    }
  }
  else
  {
    // H(div) or H(curl) use ProjectFromNodes.
    mfem::Vector elem_vals, v;
    mfem::Array<int> vdofs;
    int offset = 0;
    for (int i = 0; i < dest_mesh.GetNE(); i++)
    {
      const mfem::FiniteElement &fe = *V.FESpace()->GetFE(i);
      mfem::ElementTransformation &T = *dest_mesh.GetElementTransformation(i);
      const int elem_npts = fe.GetNodes().GetNPoints();
      elem_vals.SetSize(elem_npts * vdim);
      for (int d = 0; d < vdim; d++)
      {
        for (int j = 0; j < elem_npts; j++)
        {
          // Arrange element values byVDIM for ProjectFromNodes.
          int idx = (U.FESpace()->GetOrdering() == mfem::Ordering::byNODES)
                        ? d * npts + offset + j
                        : (offset + j) * vdim + d;
          elem_vals(j * vdim + d) = vals(idx);
        }
      }
      const auto *dof_trans = V.FESpace()->GetElementVDofs(i, vdofs);
      v.SetSize(vdofs.Size());
      fe.ProjectFromNodes(elem_vals, T, v);
      if (dof_trans)
      {
        dof_trans->TransformPrimal(v);
      }
      V.SetSubVector(vdofs, v);
      offset += elem_npts;
    }
  }
#else
  MFEM_ABORT("InterpolateFunction requires MFEM_USE_GSLIB!");
#endif
}

void InterpolateFunction(const mfem::Vector &xyz, const mfem::GridFunction &U,
                         mfem::Vector &vals, mfem::Ordering::Type ordering)
{
#if defined(MFEM_USE_GSLIB)
  // Set up a single-use interpolator, then delegate to the reusable-op overload. Callers
  // that interpolate repeatedly on a fixed mesh (e.g. wave-port voltage-path line
  // integrals during a frequency sweep) should instead build the FindPointsGSLIB once via
  // SetupInterpolator and call the overload below — the geometric Setup is the expensive
  // step and depends only on the mesh, not the field values.
  auto &src_mesh = *U.FESpace()->GetMesh();
  MFEM_VERIFY(src_mesh.GetNodes(), "Source mesh has no nodal FE space!");
  const int dim = src_mesh.SpaceDimension();
  const int npts = xyz.Size() / dim;
  auto *src_pmesh = dynamic_cast<mfem::ParMesh *>(&src_mesh);
  MPI_Comm comm = (src_pmesh) ? src_pmesh->GetComm() : MPI_COMM_SELF;
  mfem::FindPointsGSLIB op(comm);
  op.Setup(src_mesh, GSLIB_BB_TOL, GSLIB_NEWTON_TOL, npts);
  InterpolateFunction(op, xyz, U, vals, ordering);
#else
  MFEM_ABORT("InterpolateFunction requires MFEM_USE_GSLIB!");
#endif
}

void SetupInterpolator(mfem::FindPointsGSLIB &op, mfem::Mesh &mesh)
{
#if defined(MFEM_USE_GSLIB)
  op.Setup(mesh, GSLIB_BB_TOL, GSLIB_NEWTON_TOL);
#else
  MFEM_ABORT("SetupInterpolator requires MFEM_USE_GSLIB!");
#endif
}

void InterpolateFunction(mfem::FindPointsGSLIB &op, const mfem::Vector &xyz,
                         const mfem::GridFunction &U, mfem::Vector &vals,
                         mfem::Ordering::Type ordering)
{
#if defined(MFEM_USE_GSLIB)
  // Perform the interpolation using a pre-Setup point locator (op.Setup already called on
  // the source mesh). Setup is geometric — it depends only on the mesh — so a fixed-mesh
  // caller reuses op across field values without paying the O(num_elements) hash build
  // each call. The ordering of the returned values matches the source grid function.
  const int dim = U.FESpace()->GetMesh()->SpaceDimension();
  const int npts = xyz.Size() / dim;
  const int vdim = U.VectorDim();
  MFEM_VERIFY(vals.Size() == npts * vdim, "Incorrect size for interpolated values vector!");
  op.SetDefaultInterpolationValue(0.0);
  op.SetL2AvgType(mfem::FindPointsGSLIB::NONE);
  op.Interpolate(xyz, U, vals, ordering);
#else
  MFEM_ABORT("InterpolateFunction requires MFEM_USE_GSLIB!");
#endif
}

#if defined(MFEM_USE_GSLIB)
namespace
{

// Shared quadrature + dot-product core for the line integral. `interp` fills `vals`
// (byNODES ordering, length npts*vdim) with the field sampled at the quadrature points;
// it is the only part that touches GSLIB, so callers can supply either a fresh single-use
// locator or a cached pre-Setup one.
template <typename InterpFn>
double LineIntegralCore(const mfem::Vector &p1, const mfem::Vector &p2,
                        const mfem::ParGridFunction &field, int quad_order, InterpFn interp)
{
  const int dim = p1.Size();
  MFEM_VERIFY(p2.Size() == dim, "ComputeLineIntegral: p1 and p2 must have same dimension!");

  // Direction vector dl = p2 - p1. The full integral is V = ∫₀¹ F(r(t)) · dl dt.
  mfem::Vector dl(dim);
  for (int d = 0; d < dim; d++)
  {
    dl(d) = p2(d) - p1(d);
  }
  MFEM_VERIFY(dl.Norml2() > 0.0, "ComputeLineIntegral: p1 and p2 must be distinct points!");

  // Gauss-Legendre quadrature on the segment [0, 1]. The caller is responsible for
  // choosing a high enough order to resolve the field variation along the line.
  const mfem::IntegrationRule &ir = mfem::IntRules.Get(mfem::Geometry::SEGMENT, quad_order);
  const int npts = ir.GetNPoints();

  // Generate physical coordinates along the line: r(t_i) = p1 + t_i * dl.
  // Stored in byNODES ordering: [x0,x1,...,xN, y0,y1,...,yN, (z0,...,zN)].
  mfem::Vector xyz(npts * dim);
  for (int i = 0; i < npts; i++)
  {
    double t = ir.IntPoint(i).x;
    for (int d = 0; d < dim; d++)
    {
      xyz(d * npts + i) = p1(d) + t * dl(d);
    }
  }

  // Interpolate the vector field at the quadrature points.
  const int vdim = field.VectorDim();
  mfem::Vector vals(npts * vdim);
  interp(xyz, vals);

  // Compute the dot product F · dl at each quadrature point and sum. Only the first
  // min(vdim, dim) components contribute; higher components (if vdim > dim) are ignored
  // since dl has length dim.
  const int ndot = std::min(vdim, dim);
  double result = 0.0;
  for (int i = 0; i < npts; i++)
  {
    double dot = 0.0;
    for (int d = 0; d < ndot; d++)
    {
      // vals is always in byNODES ordering from InterpolateFunction.
      double val = vals(d * npts + i);
      dot += val * dl(d);
    }
    result += ir.IntPoint(i).weight * dot;
  }
  return result;
}

}  // namespace
#endif

double ComputeLineIntegral(const mfem::Vector &p1, const mfem::Vector &p2,
                           const mfem::ParGridFunction &field, int quad_order)
{
#if defined(MFEM_USE_GSLIB)
  return LineIntegralCore(
      p1, p2, field, quad_order, [&field](const mfem::Vector &xyz, mfem::Vector &vals)
      { InterpolateFunction(xyz, field, vals, mfem::Ordering::byNODES); });
#else
  MFEM_ABORT("ComputeLineIntegral requires MFEM_USE_GSLIB!");
  return 0.0;
#endif
}

double ComputeLineIntegral(mfem::FindPointsGSLIB &op, const mfem::Vector &p1,
                           const mfem::Vector &p2, const mfem::ParGridFunction &field,
                           int quad_order)
{
#if defined(MFEM_USE_GSLIB)
  // As ComputeLineIntegral, but reuses a pre-Setup point locator (op.Setup already called
  // on the field's mesh). Avoids rebuilding the GSLIB spatial hash on every call — the
  // dominant cost when integrating repeatedly on a fixed mesh.
  return LineIntegralCore(
      p1, p2, field, quad_order, [&op, &field](const mfem::Vector &xyz, mfem::Vector &vals)
      { InterpolateFunction(op, xyz, field, vals, mfem::Ordering::byNODES); });
#else
  MFEM_ABORT("ComputeLineIntegral requires MFEM_USE_GSLIB!");
  return 0.0;
#endif
}

}  // namespace fem

}  // namespace palace
