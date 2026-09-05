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
    const mfem::ParFiniteElementSpace &target_fespace, double scaling,
    bool build_gridfunction, bool build_buffer)
  : kind(kind), nd_fespace(nd_fespace), rt_fespace(rt_fespace)
{
  const bool needs_nd = kind == Kind::FIELD_E || kind == Kind::FIELD_H1 ||
                        kind == Kind::ENERGY_E || kind == Kind::POYNTING ||
                        kind == Kind::MODE_SN;
  const bool needs_rt = kind == Kind::FIELD_B || kind == Kind::ENERGY_M ||
                        kind == Kind::POYNTING || kind == Kind::MODE_SN;
  MFEM_VERIFY((!needs_nd || nd_fespace) && (!needs_rt || rt_fespace),
              "Missing finite element space for domain field evaluator!");
  MFEM_VERIFY(build_gridfunction || build_buffer,
              "Domain point evaluator has no requested output path!");
  MFEM_VERIFY((mesh.Dimension() == 2 && mesh.SpaceDimension() == 2) ||
                  (mesh.Dimension() == 3 && mesh.SpaceDimension() == 3),
              "Domain point field evaluation requires a 2D or 3D volume mesh!");
  Assemble(mesh, mat_op, target_fespace, scaling, build_gridfunction, build_buffer);
}

DomainPointFieldEvaluator::~DomainPointFieldEvaluator()
{
  fem::DestroyGroupOperators(groups);
  fem::DestroyGroupOperators(buffer_groups);
}

void DomainPointFieldEvaluator::Assemble(const Mesh &mesh, const MaterialOperator &mat_op,
                                         const mfem::ParFiniteElementSpace &target_fespace,
                                         double scaling, bool build_gridfunction,
                                         bool build_buffer)
{
  const int dim = mesh.Dimension();
  const bool scalar_magnetic_field = (dim == 2);
  const mfem::ParMesh &pmesh = mesh.Get();
  const mfem::FiniteElementSpace &mesh_fespace = *pmesh.GetNodes()->FESpace();

  const int vtu_lod = target_fespace.GetMaxElementOrder();
  if (build_buffer)
  {
    buffer_num_comp = target_fespace.GetVDim();
    buffer_bases.assign(pmesh.GetNE(), -1);
    int num_points = 0;
    for (int e = 0; e < pmesh.GetNE(); e++)
    {
      buffer_bases[e] = num_points;
      num_points +=
          mfem::GlobGeometryRefiner.Refine(pmesh.GetElementBaseGeometry(e), vtu_lod, 1)
              ->RefPts.GetNPoints();
    }
    buffer_size = num_points * buffer_num_comp;
  }

  // Group the elements by geometry type.
  std::map<mfem::Geometry::Type, std::vector<int>> geom_elems;
  for (int e = 0; e < pmesh.GetNE(); e++)
  {
    geom_elems[pmesh.GetElementGeometry(e)].push_back(e);
  }

  const bool field_value_kind =
      kind == Kind::FIELD_E || kind == Kind::FIELD_B || kind == Kind::FIELD_H1;
  std::vector<CeedIntScalar> ctx;
  if (!field_value_kind)
  {
    ctx.resize(1);
    ctx[0].second = scaling;
    const auto &coeff = (kind == Kind::ENERGY_E)
                            ? mat_op.GetPermittivityReal()
                            : ((scalar_magnetic_field || kind == Kind::MODE_SN)
                                   ? mat_op.GetCurlCurlInvPermeability()
                                   : mat_op.GetInvPermeability());
    MaterialPropertyCoefficient coeff_func(mat_op.GetAttributeToMaterial(), coeff);
    auto mat_ctx = ceed::PopulateCoefficientContext(
        ((scalar_magnetic_field && kind != Kind::ENERGY_E) || kind == Kind::MODE_SN) ? 1
                                                                                     : dim,
        &coeff_func);
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
        info.apply_qf =
            scalar_magnetic_field ? f_eval_probe_integral_22 : f_eval_probe_hdiv_33;
        info.apply_qf_path = scalar_magnetic_field
                                 ? PalaceQFunctionRelativePath(f_eval_probe_integral_22_loc)
                                 : PalaceQFunctionRelativePath(f_eval_probe_hdiv_33_loc);
        break;
      case Kind::FIELD_H1:
        info.apply_qf = (dim == 2) ? f_eval_probe_l2_22 : f_eval_probe_l2_33;
        info.apply_qf_path = (dim == 2)
                                 ? PalaceQFunctionRelativePath(f_eval_probe_l2_22_loc)
                                 : PalaceQFunctionRelativePath(f_eval_probe_l2_33_loc);
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
      case Kind::MODE_SN:
        MFEM_VERIFY(dim == 2, "Boundary-mode Sn output requires a 2D mesh!");
        info.apply_qf = f_eval_mode_sn_22;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_eval_mode_sn_22_loc);
        break;
    }

    auto AssembleAtPoints = [&](const mfem::IntegrationRule &ir,
                                CeedElemRestriction out_restr,
                                std::vector<fem::CeedGroupOperator> &operators)
    {
      const int num_pts = ir.GetNPoints();
      std::vector<ceed::CeedFunctionalFieldInput> inputs;
      std::vector<std::pair<std::string, int>> field_sources;
      if (!field_value_kind)
      {
        CeedElemRestriction attr_restr;
        CeedBasis attr_basis;
        CeedVector attr_vec;
        PalaceCeedCall(ceed, CeedElemRestrictionCreateStrided(
                                 ceed, indices.size(), 1, 1, indices.size(),
                                 CEED_STRIDES_BACKEND, &attr_restr));
        mfem::Vector Bt(num_pts), Gt(num_pts), qX(num_pts), qW(num_pts);
        Bt = 1.0;
        Gt = 0.0;
        qX = 0.0;
        qW = 0.0;
        PalaceCeedCall(ceed, CeedBasisCreateH1(ceed, CEED_TOPOLOGY_LINE, 1, 1, num_pts,
                                               Bt.GetData(), Gt.GetData(), qX.GetData(),
                                               qW.GetData(), &attr_basis));
        ceed::InitCeedVector(elem_attrs.back(), ceed, &attr_vec);
        inputs.push_back(
            {"attr", attr_vec, attr_restr, attr_basis, ceed::EvalMode::Interp});
        scratch.vecs.push_back(attr_vec);
        scratch.restrs.push_back(attr_restr);
        scratch.bases.push_back(attr_basis);
      }

      CeedElemRestriction mesh_restr = FiniteElementSpace::BuildCeedElemRestriction(
          mesh_fespace, ceed, geom, indices, /*is_interp*/ true);
      const mfem::FiniteElement *mesh_fe =
          mesh_fespace.FEColl()->FiniteElementForGeometry(geom);
      CeedBasis mesh_basis;
      ceed::InitCachedBasisFromRule(*mesh_fe, ir, mesh_fespace.GetVDim(), ceed,
                                    &mesh_basis);
      CeedVector mesh_nodes_vec;
      ceed::InitCeedVector(*mesh_fespace.GetMesh()->GetNodes(), ceed, &mesh_nodes_vec);
      inputs.push_back({"x", mesh_nodes_vec, mesh_restr, mesh_basis, ceed::EvalMode::Grad});
      scratch.vecs.push_back(mesh_nodes_vec);
      scratch.restrs.push_back(mesh_restr);
      scratch.bases.push_back(mesh_basis);

      auto AddFieldInput = [&](const std::string &name, int source,
                               const mfem::ParFiniteElementSpace &fespace)
      {
        CeedElemRestriction restr;
        CeedBasis basis;
        CeedVector vec;
        ceed::InitRestriction(fespace, indices, false, /*is_interp*/ true, false, ceed,
                              &restr);
        const mfem::FiniteElement *fe = fespace.FEColl()->FiniteElementForGeometry(geom);
        ceed::InitCachedBasisFromRule(*fe, ir, fespace.GetVDim(), ceed, &basis);
        ceed::InitCeedVector(field_staging, ceed, &vec);
        inputs.push_back({name, vec, restr, basis, ceed::EvalMode::Interp});
        field_sources.emplace_back(name, source);
        scratch.vecs.push_back(vec);
        scratch.restrs.push_back(restr);
        scratch.bases.push_back(basis);
      };
      if (kind == Kind::FIELD_E || kind == Kind::FIELD_H1 || kind == Kind::ENERGY_E ||
          kind == Kind::POYNTING || kind == Kind::MODE_SN)
      {
        AddFieldInput("u_1", 0, *nd_fespace);
      }
      if (kind == Kind::FIELD_B || kind == Kind::ENERGY_M || kind == Kind::POYNTING ||
          kind == Kind::MODE_SN)
      {
        AddFieldInput((kind == Kind::POYNTING || kind == Kind::MODE_SN) ? "u_2" : "u_1", 1,
                      *rt_fespace);
      }

      CeedOperator op;
      ceed::AssembleCeedPointEvaluator(info, ctx.empty() ? nullptr : ctx.data(),
                                       ctx.size() * sizeof(CeedIntScalar), ceed, inputs,
                                       target_fespace.GetVDim(), out_restr, &op);
      operators.push_back({ceed, op, std::move(field_sources)});
      operators.back().mesh_nodes = pmesh.GetNodes();
      operators.back().mesh_node_fields = {"grad_x"};
      fem::CacheGroupOperatorFieldVectors(operators.back());
    };

    if (!field_value_kind)
    {
      auto &elem_attr = elem_attrs.emplace_back(indices.size());
      const auto &loc_attr = mesh.GetCeedAttributes();
      for (std::size_t k = 0; k < indices.size(); k++)
      {
        elem_attr[k] = loc_attr.at(pmesh.GetAttribute(indices[k]));
      }
    }
    if (build_gridfunction)
    {
      const mfem::FiniteElement *target_fe =
          target_fespace.FEColl()->FiniteElementForGeometry(geom);
      MFEM_VERIFY(target_fe, "Unable to get target finite element for field evaluator!");
      CeedElemRestriction out_restr;
      ceed::InitRestriction(target_fespace, indices, false, /*is_interp*/ true, false, ceed,
                            &out_restr);
      scratch.restrs.push_back(out_restr);
      AssembleAtPoints(target_fe->GetNodes(), out_restr, groups);
    }
    if (build_buffer)
    {
      const mfem::IntegrationRule &ir =
          mfem::GlobGeometryRefiner.Refine(geom, vtu_lod, 1)->RefPts;
      const int num_pts = ir.GetNPoints();
      std::vector<CeedInt> offsets(indices.size() * num_pts);
      for (std::size_t e = 0; e < indices.size(); e++)
      {
        const int point_base = buffer_bases[indices[e]];
        for (int j = 0; j < num_pts; j++)
        {
          offsets[e * num_pts + j] = (point_base + j) * buffer_num_comp;
        }
      }
      CeedElemRestriction out_restr;
      PalaceCeedCall(
          ceed, CeedElemRestrictionCreate(ceed, static_cast<CeedInt>(indices.size()),
                                          num_pts, buffer_num_comp, 1,
                                          static_cast<CeedSize>(buffer_size), CEED_MEM_HOST,
                                          CEED_COPY_VALUES, offsets.data(), &out_restr));
      scratch.restrs.push_back(out_restr);
      AssembleAtPoints(ir, out_restr, buffer_groups);
    }
  }

  // Passive field vectors are re-pointed to caller-owned data on every apply.
  fem::DetachGroupOperatorFieldVectors(groups);
  fem::DetachGroupOperatorFieldVectors(buffer_groups);
  buffer_bases.clear();
  buffer_bases.shrink_to_fit();
  field_staging.Destroy();
  MFEM_ASSERT(field_staging.Capacity() == 0,
              "Domain evaluator staging allocation was not released!");
}

void DomainPointFieldEvaluator::Eval(const GridFunction *E, const GridFunction *B,
                                     Vector &out) const
{
  const bool needs_e =
      kind == Kind::ENERGY_E || kind == Kind::POYNTING || kind == Kind::MODE_SN;
  const bool needs_b =
      kind == Kind::ENERGY_M || kind == Kind::POYNTING || kind == Kind::MODE_SN;
  MFEM_VERIFY((!needs_e || E) && (!needs_b || B),
              "Missing field grid function for domain field evaluator!");
  if (E && B)
  {
    MFEM_VERIFY(E->HasImag() == B->HasImag(),
                "Mismatched real and complex fields in domain field evaluator!");
  }

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
  MFEM_VERIFY(kind == Kind::FIELD_E || kind == Kind::FIELD_B || kind == Kind::FIELD_H1,
              "Vector EvalBuffer is only valid for primary field evaluators!");
  MFEM_VERIFY(buffer.Size() == buffer_size, "Invalid domain point buffer size!");
  buffer = 0.0;
  fem::ApplyAddGroupOperators(
      buffer_groups,
      {(kind == Kind::FIELD_E || kind == Kind::FIELD_H1) ? &u : nullptr,
       kind == Kind::FIELD_B ? &u : nullptr},
      buffer);
}

void DomainPointFieldEvaluator::EvalBuffer(const GridFunction *E, const GridFunction *B,
                                           Vector &buffer) const
{
  const bool needs_e =
      kind == Kind::ENERGY_E || kind == Kind::POYNTING || kind == Kind::MODE_SN;
  const bool needs_b =
      kind == Kind::ENERGY_M || kind == Kind::POYNTING || kind == Kind::MODE_SN;
  MFEM_VERIFY((!needs_e || E) && (!needs_b || B),
              "Missing field grid function for domain field evaluator!");
  MFEM_VERIFY(buffer.Size() == buffer_size, "Invalid domain point buffer size!");
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
