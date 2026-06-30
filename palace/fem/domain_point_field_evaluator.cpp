// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "domain_point_field_evaluator.hpp"

#include <algorithm>
#include <map>
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/libceed/basis.hpp"
#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/functional.hpp"
#include "fem/libceed/integrator.hpp"
#include "fem/libceed/restriction.hpp"
#include "fem/mesh.hpp"
#include "models/materialoperator.hpp"
#include "utils/diagnostic.hpp"

PalacePragmaDiagnosticPush
PalacePragmaDiagnosticDisableUnused

#include "fem/qfunctions/22/eval_22_qf.h"
#include "fem/qfunctions/33/eval_33_qf.h"

PalacePragmaDiagnosticPop

namespace palace
{

namespace
{

// Holds libCEED object references created during operator assembly for destruction once
// the assembled operator owns them.
struct CeedAssemblyScratch
{
  Ceed ceed;
  std::vector<CeedVector> vecs;
  std::vector<CeedElemRestriction> restrs;
  std::vector<CeedBasis> bases;

  CeedAssemblyScratch(Ceed ceed) : ceed(ceed) {}
  CeedAssemblyScratch(const CeedAssemblyScratch &) = delete;
  CeedAssemblyScratch &operator=(const CeedAssemblyScratch &) = delete;

  ~CeedAssemblyScratch()
  {
    for (auto &v : vecs)
    {
      PalaceCeedCall(ceed, CeedVectorDestroy(&v));
    }
    for (auto &r : restrs)
    {
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&r));
    }
    for (auto &b : bases)
    {
      PalaceCeedCall(ceed, CeedBasisDestroy(&b));
    }
  }
};

}  // namespace

DomainPointFieldEvaluator::DomainPointFieldEvaluator(
    Kind kind, const Mesh &mesh, const MaterialOperator &mat_op,
    const mfem::ParFiniteElementSpace *nd_fespace,
    const mfem::ParFiniteElementSpace *rt_fespace,
    const mfem::ParFiniteElementSpace &target_fespace, double scaling)
  : kind(kind), nd_fespace(nd_fespace), rt_fespace(rt_fespace)
{
  const bool needs_nd =
      kind == Kind::FIELD_E || kind == Kind::ENERGY_E || kind == Kind::POYNTING;
  const bool needs_rt =
      kind == Kind::FIELD_B || kind == Kind::ENERGY_M || kind == Kind::POYNTING;
  MFEM_VERIFY((!needs_nd || nd_fespace) && (!needs_rt || rt_fespace),
              "Missing finite element space for domain field evaluator!");
  Assemble(mesh, mat_op, target_fespace, scaling);
}

DomainPointFieldEvaluator::~DomainPointFieldEvaluator()
{
  fem::DestroyGroupOperators(groups);
  fem::DestroyGroupOperators(buffer_groups);
}

void DomainPointFieldEvaluator::Assemble(const Mesh &mesh, const MaterialOperator &mat_op,
                                         const mfem::ParFiniteElementSpace &target_fespace,
                                         double scaling)
{
  const int dim = mesh.Dimension();
  const int sdim = mesh.SpaceDimension();
  if (!((dim == 2 && sdim == 2) || (dim == 3 && sdim == 3)))
  {
    valid = false;
    return;
  }
  const bool scalar_magnetic_field = (dim == 2);
  const mfem::ParMesh &pmesh = mesh.Get();
  const mfem::FiniteElementSpace &mesh_fespace = *pmesh.GetNodes()->FESpace();

  // Point-major buffer layout matching VTK tuple order and MFEM's domain VTU point
  // order (element order, then refined/nodal points within each element). The output
  // restriction scatters each point tuple contiguously so CeedParaViewDataCollection can
  // write the payload without repacking.
  buffer_num_comp = target_fespace.GetVDim();
  buffer_bases.assign(pmesh.GetNE(), -1);
  int buffer_points = 0;
  const int vtu_lod = target_fespace.GetMaxElementOrder();
  for (int e = 0; e < pmesh.GetNE(); e++)
  {
    const auto *RefG =
        mfem::GlobGeometryRefiner.Refine(pmesh.GetElementBaseGeometry(e), vtu_lod, 1);
    buffer_bases[e] = buffer_points;
    buffer_points += RefG->RefPts.GetNPoints();
  }
  buffer_size = buffer_points * buffer_num_comp;

  // Group the elements by geometry type.
  std::map<mfem::Geometry::Type, std::vector<int>> geom_elems;
  for (int e = 0; e < pmesh.GetNE(); e++)
  {
    geom_elems[pmesh.GetElementGeometry(e)].push_back(e);
  }

  // QFunction context: derived quantities need a scaling factor followed by the material
  // property table. Plain field probes only need geometry and field values.
  const bool field_value_kind = (kind == Kind::FIELD_E || kind == Kind::FIELD_B);
  std::vector<CeedIntScalar> ctx;
  if (!field_value_kind)
  {
    ctx.resize(1);
    ctx[0].second = scaling;
    const auto &coeff = (kind == Kind::ENERGY_E)
                            ? mat_op.GetPermittivityReal()
                            : (scalar_magnetic_field ? mat_op.GetCurlCurlInvPermeability()
                                                     : mat_op.GetInvPermeability());
    MaterialPropertyCoefficient coeff_func(mat_op.GetAttributeToMaterial(), coeff);
    auto mat_ctx = ceed::PopulateCoefficientContext(
        (scalar_magnetic_field && kind != Kind::ENERGY_E) ? 1 : dim, &coeff_func);
    ctx.insert(ctx.end(), mat_ctx.begin(), mat_ctx.end());
  }

  field_staging.SetSize(std::max(nd_fespace ? nd_fespace->GetVSize() : 0,
                                 rt_fespace ? rt_fespace->GetVSize() : 0));
  field_staging.UseDevice(true);
  field_staging = 0.0;

  Ceed ceed = ceed::internal::GetCeedObjects()[0];
  for (const auto &geom_indices : geom_elems)
  {
    const mfem::Geometry::Type geom = geom_indices.first;
    const auto &indices = geom_indices.second;
    CeedAssemblyScratch scratch(ceed);

    // Evaluation points are the nodal points of the (interpolatory) target space.
    const mfem::FiniteElement *target_fe =
        target_fespace.FEColl()->FiniteElementForGeometry(geom);
    MFEM_VERIFY(target_fe, "Unable to get target finite element for field evaluator!");
    const mfem::IntegrationRule &nodes_ir = target_fe->GetNodes();
    const int num_pts = nodes_ir.GetNPoints();

    // Element attributes (libCEED local, for material lookup) and mesh node gradients
    // (for on the fly geometry evaluation at the points).
    std::vector<ceed::CeedFunctionalFieldInput> inputs;
    std::vector<std::pair<std::string, int>> field_sources;
    if (!field_value_kind)
    {
      auto &elem_attr = elem_attrs.emplace_back(indices.size());
      const auto &loc_attr = mesh.GetCeedAttributes();
      for (std::size_t k = 0; k < indices.size(); k++)
      {
        elem_attr[k] = loc_attr.at(pmesh.GetAttribute(indices[k]));
      }
      CeedElemRestriction attr_restr;
      CeedBasis attr_basis;
      CeedVector attr_vec;
      PalaceCeedCall(
          ceed, CeedElemRestrictionCreateStrided(ceed, indices.size(), 1, 1, indices.size(),
                                                 CEED_STRIDES_BACKEND, &attr_restr));
      {
        // Note: ceed::GetCeedTopology(CEED_TOPOLOGY_LINE) == 1.
        mfem::Vector Bt(num_pts), Gt(num_pts), qX(num_pts), qW(num_pts);
        Bt = 1.0;
        Gt = 0.0;
        qX = 0.0;
        qW = 0.0;
        PalaceCeedCall(ceed, CeedBasisCreateH1(ceed, CEED_TOPOLOGY_LINE, 1, 1, num_pts,
                                               Bt.GetData(), Gt.GetData(), qX.GetData(),
                                               qW.GetData(), &attr_basis));
      }
      ceed::InitCeedVector(elem_attr, ceed, &attr_vec);
      inputs.push_back({"attr", attr_vec, attr_restr, attr_basis, ceed::EvalMode::Interp});
      scratch.vecs.push_back(attr_vec);
      scratch.restrs.push_back(attr_restr);
      scratch.bases.push_back(attr_basis);
    }
    {
      CeedElemRestriction mesh_restr = FiniteElementSpace::BuildCeedElemRestriction(
          mesh_fespace, ceed, geom, indices, /*is_interp*/ true);
      const mfem::FiniteElement *mesh_fe =
          mesh_fespace.FEColl()->FiniteElementForGeometry(geom);
      CeedBasis mesh_basis;
      ceed::InitBasisAtPoints(*mesh_fe, nodes_ir, mesh_fespace.GetVDim(), ceed,
                              &mesh_basis);
      CeedVector mesh_nodes_vec;
      ceed::InitCeedVector(*mesh_fespace.GetMesh()->GetNodes(), ceed, &mesh_nodes_vec);
      inputs.push_back({"x", mesh_nodes_vec, mesh_restr, mesh_basis, ceed::EvalMode::Grad});
      scratch.vecs.push_back(mesh_nodes_vec);
      scratch.restrs.push_back(mesh_restr);
      scratch.bases.push_back(mesh_basis);
    }

    // Field inputs evaluated at the nodal points.
    auto AddFieldInput =
        [&](const std::string &name, int source, const mfem::ParFiniteElementSpace &fespace)
    {
      CeedElemRestriction restr;
      CeedBasis basis;
      CeedVector vec;
      ceed::InitRestriction(fespace, indices, false, /*is_interp*/ true, false, ceed,
                            &restr);
      const mfem::FiniteElement *fe = fespace.FEColl()->FiniteElementForGeometry(geom);
      ceed::InitBasisAtPoints(*fe, nodes_ir, fespace.GetVDim(), ceed, &basis);
      ceed::InitCeedVector(field_staging, ceed, &vec);
      inputs.push_back({name, vec, restr, basis, ceed::EvalMode::Interp});
      field_sources.emplace_back(name, source);
      scratch.vecs.push_back(vec);
      scratch.restrs.push_back(restr);
      scratch.bases.push_back(basis);
    };
    if (kind == Kind::FIELD_E || kind == Kind::ENERGY_E || kind == Kind::POYNTING)
    {
      AddFieldInput("u_1", 0, *nd_fespace);
    }
    if (kind == Kind::FIELD_B || kind == Kind::ENERGY_M || kind == Kind::POYNTING)
    {
      AddFieldInput(kind == Kind::POYNTING ? "u_2" : "u_1", 1, *rt_fespace);
    }

    // Output restriction scattering nodal values into the target grid function.
    CeedElemRestriction out_restr;
    ceed::InitRestriction(target_fespace, indices, false, /*is_interp*/ true, false, ceed,
                          &out_restr);
    scratch.restrs.push_back(out_restr);

    ceed::CeedQFunctionInfo info;
    switch (kind)
    {
      case Kind::FIELD_E:
        info.apply_qf = (dim == 2) ? f_eval_probe_hcurl_22 : f_eval_probe_hcurl_33;
        info.apply_qf_path = (dim == 2)
                                 ? PalaceQFunctionRelativePath(f_eval_probe_hcurl_22_loc)
                                 : PalaceQFunctionRelativePath(f_eval_probe_hcurl_33_loc);
        break;
      case Kind::FIELD_B:
        info.apply_qf = (dim == 2) ? f_eval_probe_hdiv_22 : f_eval_probe_hdiv_33;
        info.apply_qf_path = (dim == 2)
                                 ? PalaceQFunctionRelativePath(f_eval_probe_hdiv_22_loc)
                                 : PalaceQFunctionRelativePath(f_eval_probe_hdiv_33_loc);
        break;
      case Kind::ENERGY_E:
        info.apply_qf = (dim == 2) ? f_eval_energy_e_22 : f_eval_energy_e_33;
        info.apply_qf_path = (dim == 2)
                                 ? PalaceQFunctionRelativePath(f_eval_energy_e_22_loc)
                                 : PalaceQFunctionRelativePath(f_eval_energy_e_33_loc);
        break;
      case Kind::ENERGY_M:
        info.apply_qf = (dim == 2) ? f_eval_energy_m_22 : f_eval_energy_m_33;
        info.apply_qf_path = (dim == 2)
                                 ? PalaceQFunctionRelativePath(f_eval_energy_m_22_loc)
                                 : PalaceQFunctionRelativePath(f_eval_energy_m_33_loc);
        break;
      case Kind::POYNTING:
        info.apply_qf = (dim == 2) ? f_eval_poynting_22 : f_eval_poynting_33;
        info.apply_qf_path = (dim == 2)
                                 ? PalaceQFunctionRelativePath(f_eval_poynting_22_loc)
                                 : PalaceQFunctionRelativePath(f_eval_poynting_33_loc);
        break;
    }

    CeedOperator op;
    ceed::AssembleCeedPointEvaluator(info, ctx.empty() ? nullptr : ctx.data(),
                                     ctx.size() * sizeof(CeedIntScalar), ceed, inputs,
                                     target_fespace.GetVDim(), out_restr, &op);
    groups.push_back({ceed, op, std::move(field_sources)});

    // Assemble a second operator for the direct VTU point-buffer path. This evaluates at
    // MFEM's refined VTU points (not necessarily the same order as the target L2 nodes)
    // and scatters into point-major tuples consumed by CeedParaViewDataCollection.
    const mfem::IntegrationRule &vtu_ir =
        mfem::GlobGeometryRefiner.Refine(geom, vtu_lod, 1)->RefPts;
    const int num_vtu_pts = vtu_ir.GetNPoints();
    std::vector<ceed::CeedFunctionalFieldInput> buffer_inputs;
    std::vector<std::pair<std::string, int>> buffer_field_sources;
    if (!field_value_kind)
    {
      auto &elem_attr = elem_attrs.emplace_back(indices.size());
      const auto &loc_attr = mesh.GetCeedAttributes();
      for (std::size_t k = 0; k < indices.size(); k++)
      {
        elem_attr[k] = loc_attr.at(pmesh.GetAttribute(indices[k]));
      }
      CeedElemRestriction attr_restr;
      CeedBasis attr_basis;
      CeedVector attr_vec;
      PalaceCeedCall(
          ceed, CeedElemRestrictionCreateStrided(ceed, indices.size(), 1, 1, indices.size(),
                                                 CEED_STRIDES_BACKEND, &attr_restr));
      {
        mfem::Vector Bt(num_vtu_pts), Gt(num_vtu_pts), qX(num_vtu_pts), qW(num_vtu_pts);
        Bt = 1.0;
        Gt = 0.0;
        qX = 0.0;
        qW = 0.0;
        PalaceCeedCall(ceed, CeedBasisCreateH1(ceed, CEED_TOPOLOGY_LINE, 1, 1, num_vtu_pts,
                                               Bt.GetData(), Gt.GetData(), qX.GetData(),
                                               qW.GetData(), &attr_basis));
      }
      ceed::InitCeedVector(elem_attr, ceed, &attr_vec);
      buffer_inputs.push_back(
          {"attr", attr_vec, attr_restr, attr_basis, ceed::EvalMode::Interp});
      scratch.vecs.push_back(attr_vec);
      scratch.restrs.push_back(attr_restr);
      scratch.bases.push_back(attr_basis);
    }
    {
      CeedElemRestriction mesh_restr = FiniteElementSpace::BuildCeedElemRestriction(
          mesh_fespace, ceed, geom, indices, /*is_interp*/ true);
      const mfem::FiniteElement *mesh_fe =
          mesh_fespace.FEColl()->FiniteElementForGeometry(geom);
      CeedBasis mesh_basis;
      ceed::InitBasisAtPoints(*mesh_fe, vtu_ir, mesh_fespace.GetVDim(), ceed, &mesh_basis);
      CeedVector mesh_nodes_vec;
      ceed::InitCeedVector(*mesh_fespace.GetMesh()->GetNodes(), ceed, &mesh_nodes_vec);
      buffer_inputs.push_back(
          {"x", mesh_nodes_vec, mesh_restr, mesh_basis, ceed::EvalMode::Grad});
      scratch.vecs.push_back(mesh_nodes_vec);
      scratch.restrs.push_back(mesh_restr);
      scratch.bases.push_back(mesh_basis);
    }
    auto AddBufferFieldInput =
        [&](const std::string &name, int source, const mfem::ParFiniteElementSpace &fespace)
    {
      CeedElemRestriction restr;
      CeedBasis basis;
      CeedVector vec;
      ceed::InitRestriction(fespace, indices, false, /*is_interp*/ true, false, ceed,
                            &restr);
      const mfem::FiniteElement *fe = fespace.FEColl()->FiniteElementForGeometry(geom);
      ceed::InitBasisAtPoints(*fe, vtu_ir, fespace.GetVDim(), ceed, &basis);
      ceed::InitCeedVector(field_staging, ceed, &vec);
      buffer_inputs.push_back({name, vec, restr, basis, ceed::EvalMode::Interp});
      buffer_field_sources.emplace_back(name, source);
      scratch.vecs.push_back(vec);
      scratch.restrs.push_back(restr);
      scratch.bases.push_back(basis);
    };
    if (kind == Kind::FIELD_E || kind == Kind::ENERGY_E || kind == Kind::POYNTING)
    {
      AddBufferFieldInput("u_1", 0, *nd_fespace);
    }
    if (kind == Kind::FIELD_B || kind == Kind::ENERGY_M || kind == Kind::POYNTING)
    {
      AddBufferFieldInput(kind == Kind::POYNTING ? "u_2" : "u_1", 1, *rt_fespace);
    }
    CeedElemRestriction buffer_out_restr;
    {
      std::vector<CeedInt> offsets(indices.size() * num_vtu_pts);
      for (std::size_t e = 0; e < indices.size(); e++)
      {
        const int point_base = buffer_bases[indices[e]];
        for (int j = 0; j < num_vtu_pts; j++)
        {
          offsets[e * num_vtu_pts + j] = (point_base + j) * buffer_num_comp;
        }
      }
      PalaceCeedCall(ceed, CeedElemRestrictionCreate(
                               ceed, static_cast<CeedInt>(indices.size()), num_vtu_pts,
                               buffer_num_comp, 1, (CeedSize)buffer_size, CEED_MEM_HOST,
                               CEED_COPY_VALUES, offsets.data(), &buffer_out_restr));
    }
    scratch.restrs.push_back(buffer_out_restr);
    CeedOperator buffer_op;
    ceed::AssembleCeedPointEvaluator(
        info, ctx.empty() ? nullptr : ctx.data(), ctx.size() * sizeof(CeedIntScalar), ceed,
        buffer_inputs, buffer_num_comp, buffer_out_restr, &buffer_op);
    buffer_groups.push_back({ceed, buffer_op, std::move(buffer_field_sources)});
  }
}

void DomainPointFieldEvaluator::Eval(const GridFunction *E, const GridFunction *B,
                                     Vector &out) const
{
  MFEM_VERIFY(valid, "Eval called on an invalid (unassembled) DomainPointFieldEvaluator!");
  MFEM_VERIFY((E || kind == Kind::FIELD_B || kind == Kind::ENERGY_M) &&
                  (B || kind == Kind::FIELD_E || kind == Kind::ENERGY_E),
              "Missing field grid function for domain field evaluator!");
  out = 0.0;
  fem::ApplyAddGroupOperators(groups, {E ? &E->Real() : nullptr, B ? &B->Real() : nullptr},
                              out);
  if (E ? E->HasImag() : B->HasImag())
  {
    fem::ApplyAddGroupOperators(groups,
                                {E ? &E->Imag() : nullptr, B ? &B->Imag() : nullptr}, out);
  }
}

void DomainPointFieldEvaluator::EvalBuffer(const Vector &u, Vector &buffer) const
{
  MFEM_VERIFY(valid,
              "EvalBuffer called on an invalid (unassembled) DomainPointFieldEvaluator!");
  MFEM_VERIFY(kind == Kind::FIELD_E || kind == Kind::FIELD_B,
              "Vector EvalBuffer is only valid for linear field evaluators!");
  MFEM_ASSERT(buffer.Size() == buffer_size,
              "Invalid buffer size for DomainPointFieldEvaluator::EvalBuffer!");
  buffer = 0.0;
  fem::ApplyAddGroupOperators(buffer_groups,
                              {kind == Kind::FIELD_E ? &u : nullptr,
                               kind == Kind::FIELD_B ? &u : nullptr, nullptr, nullptr},
                              buffer);
}

void DomainPointFieldEvaluator::EvalBuffer(const GridFunction *E, const GridFunction *B,
                                           Vector &buffer) const
{
  MFEM_VERIFY(valid,
              "EvalBuffer called on an invalid (unassembled) DomainPointFieldEvaluator!");
  MFEM_VERIFY((E || kind == Kind::FIELD_B || kind == Kind::ENERGY_M) &&
                  (B || kind == Kind::FIELD_E || kind == Kind::ENERGY_E),
              "Missing field grid function for domain field evaluator!");
  MFEM_ASSERT(buffer.Size() == buffer_size,
              "Invalid buffer size for DomainPointFieldEvaluator::EvalBuffer!");
  buffer = 0.0;
  fem::ApplyAddGroupOperators(buffer_groups,
                              {E ? &E->Real() : nullptr, B ? &B->Real() : nullptr}, buffer);
  if (E ? E->HasImag() : B->HasImag())
  {
    fem::ApplyAddGroupOperators(
        buffer_groups, {E ? &E->Imag() : nullptr, B ? &B->Imag() : nullptr}, buffer);
  }
}

}  // namespace palace
