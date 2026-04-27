// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "spaceoperator.hpp"

#include <set>
#include <type_traits>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "fem/libceed/ceed.hpp"  // for <ceed.h> before coeff_qf.h
#include "fem/mesh.hpp"
#include "fem/multigrid.hpp"
#include "fem/qfunctions/coeff/coeff_qf.h"
#include "linalg/hypre.hpp"
#include "linalg/rap.hpp"
#include "models/pml.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

using namespace std::complex_literals;

SpaceOperator::SpaceOperator(const config::SolverData &solver,
                             const config::DomainData &domains,
                             const config::BoundaryData &boundaries,
                             ProblemType problem_type, const Units &units,
                             const std::vector<std::unique_ptr<Mesh>> &mesh)
  : pc_mat_real(solver.linear.pc_mat_real), pc_mat_shifted(solver.linear.pc_mat_shifted),
    print_hdr(true), print_prec_hdr(true),
    dbc_attr(SetUpBoundaryProperties(boundaries.pec, *mesh.back())),
    nd_fecs(fem::ConstructFECollections<mfem::ND_FECollection>(
        solver.order, mesh.back()->Dimension(), solver.linear.mg_max_levels,
        solver.linear.mg_coarsening, false)),
    h1_fecs(fem::ConstructFECollections<mfem::H1_FECollection>(
        solver.order, mesh.back()->Dimension(), solver.linear.mg_max_levels,
        solver.linear.mg_coarsening, false)),
    rt_fecs(fem::ConstructFECollections<mfem::RT_FECollection>(
        solver.order - 1, mesh.back()->Dimension(),
        solver.linear.estimator_mg ? solver.linear.mg_max_levels : 1,
        solver.linear.mg_coarsening, false)),
    nd_fespaces(fem::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
        solver.linear.mg_max_levels, mesh, nd_fecs, &dbc_attr, &nd_dbc_tdof_lists)),
    h1_fespaces(fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
        solver.linear.mg_max_levels, mesh, h1_fecs, &dbc_attr, &h1_dbc_tdof_lists)),
    rt_fespaces(fem::ConstructFiniteElementSpaceHierarchy<mfem::RT_FECollection>(
        solver.linear.estimator_mg ? solver.linear.mg_max_levels : 1, mesh, rt_fecs)),
    mat_op(domains.materials, boundaries.periodic, problem_type, *mesh.back()),
    current_dipole_op(domains.current_dipole, units, *mesh.back()),
    farfield_op(boundaries.farfield, problem_type, mat_op, *mesh.back()),
    surf_sigma_op(boundaries.conductivity, problem_type, units, mat_op, *mesh.back()),
    surf_z_op(boundaries.impedance, boundaries.cracked_attributes, units, mat_op,
              *mesh.back()),
    lumped_port_op(boundaries.lumpedport, units, mat_op, *mesh.back()),
    wave_port_op(boundaries, solver, problem_type, units, mat_op, GetNDSpace(),
                 GetH1Space()),
    surf_j_op(boundaries.current, *mesh.back()),
    port_excitation_helper(lumped_port_op, wave_port_op, surf_j_op, current_dipole_op)
{
  // Check Excitations.
  if (problem_type == ProblemType::DRIVEN)
  {
    MFEM_VERIFY(!port_excitation_helper.Empty(),
                "Driven problems must specify at least one excitation!");
  }
  else if (problem_type == ProblemType::EIGENMODE)
  {
    MFEM_VERIFY(port_excitation_helper.Empty(),
                "Eigenmode problems must not specify any excitation!");
  }
  else if (problem_type == ProblemType::TRANSIENT)
  {
    MFEM_VERIFY(
        port_excitation_helper.Size() == 1,
        "Transient problems currently only support a single excitation per simulation!");
  }
  else
  {
    MFEM_ABORT("Internal Error: Solver type incompatible with SpaceOperator.");
  }

  // Finalize setup.
  CheckBoundaryProperties();

  // Print essential BC information.
  if (dbc_attr.Size())
  {
    Mpi::Print("\nConfiguring Dirichlet PEC BC at attributes:\n");
    utils::PrettyPrint(dbc_attr);
  }
}

SpaceOperator::SpaceOperator(const IoData &iodata,
                             const std::vector<std::unique_ptr<Mesh>> &mesh)
  : SpaceOperator(iodata.solver, iodata.domains, iodata.boundaries, iodata.problem.type,
                  iodata.units, mesh)
{
}

mfem::Array<int> SpaceOperator::SetUpBoundaryProperties(const config::PecBoundaryData &pec,
                                                        const mfem::ParMesh &mesh)
{
  // Check that boundary attributes have been specified correctly.
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> bdr_attr_marker;
  if (!pec.empty())
  {
    bdr_attr_marker.SetSize(bdr_attr_max);
    bdr_attr_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    std::set<int> bdr_warn_list;
    for (auto attr : pec.attributes)
    {
      // MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
      //             "PEC boundary attribute tags must be non-negative and correspond to "
      //             "attributes in the mesh!");
      // MFEM_VERIFY(bdr_attr_marker[attr - 1],
      //             "Unknown PEC boundary attribute " << attr << "!");
      if (attr <= 0 || attr > bdr_attr_max || !bdr_attr_marker[attr - 1])
      {
        bdr_warn_list.insert(attr);
      }
    }
    if (!bdr_warn_list.empty())
    {
      Mpi::Print("\n");
      Mpi::Warning("Unknown PEC boundary attributes!\nSolver will just ignore them!");
      utils::PrettyPrint(bdr_warn_list, "Boundary attribute list:");
      Mpi::Print("\n");
    }
  }

  // Mark selected boundary attributes from the mesh as essential (Dirichlet).
  mfem::Array<int> dbc_bcs;
  dbc_bcs.Reserve(static_cast<int>(pec.attributes.size()));
  for (auto attr : pec.attributes)
  {
    if (attr <= 0 || attr > bdr_attr_max || !bdr_attr_marker[attr - 1])
    {
      continue;  // Can just ignore if wrong
    }
    dbc_bcs.Append(attr);
  }
  return dbc_bcs;
}

void SpaceOperator::CheckBoundaryProperties()
{
  // Mark selected boundary attributes from the mesh as having some Dirichlet, Neumann, or
  // mixed BC applied.
  const mfem::ParMesh &mesh = GetMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  const auto dbc_marker = mesh::AttrToMarker(bdr_attr_max, dbc_attr);
  const auto farfield_marker = mesh::AttrToMarker(bdr_attr_max, farfield_op.GetAttrList());
  const auto surf_sigma_marker =
      mesh::AttrToMarker(bdr_attr_max, surf_sigma_op.GetAttrList());
  const auto surf_z_Rs_marker = mesh::AttrToMarker(bdr_attr_max, surf_z_op.GetRsAttrList());
  const auto surf_z_Ls_marker = mesh::AttrToMarker(bdr_attr_max, surf_z_op.GetLsAttrList());
  const auto lumped_port_Rs_marker =
      mesh::AttrToMarker(bdr_attr_max, lumped_port_op.GetRsAttrList());
  const auto lumped_port_Ls_marker =
      mesh::AttrToMarker(bdr_attr_max, lumped_port_op.GetLsAttrList());
  const auto wave_port_marker =
      mesh::AttrToMarker(bdr_attr_max, wave_port_op.GetAttrList());
  mfem::Array<int> aux_bdr_marker(dbc_marker.Size());
  for (int i = 0; i < dbc_marker.Size(); i++)
  {
    aux_bdr_marker[i] =
        (dbc_marker[i] || farfield_marker[i] || surf_sigma_marker[i] ||
         surf_z_Rs_marker[i] || surf_z_Ls_marker[i] || lumped_port_Rs_marker[i] ||
         lumped_port_Ls_marker[i] || wave_port_marker[i]);
    if (aux_bdr_marker[i])
    {
      aux_bdr_attr.Append(i + 1);
    }
  }
  // aux_bdr_marker = 1;  // Mark all boundaries (including material interfaces
  //                      // added during mesh preprocessing)
  //                      // As tested, this does not eliminate all DC modes!
  for (std::size_t l = 0; l < GetH1Spaces().GetNumLevels(); l++)
  {
    GetH1Spaces().GetFESpaceAtLevel(l).Get().GetEssentialTrueDofs(
        aux_bdr_marker, aux_bdr_tdof_lists.emplace_back());
  }

  // A final check that no boundary attribute is assigned multiple boundary conditions.
  const auto surf_z_marker = mesh::AttrToMarker(bdr_attr_max, surf_z_op.GetAttrList());
  const auto lumped_port_marker =
      mesh::AttrToMarker(bdr_attr_max, lumped_port_op.GetAttrList());
  const auto surf_j_marker = mesh::AttrToMarker(bdr_attr_max, surf_j_op.GetAttrList());
  for (int i = 0; i < dbc_marker.Size(); i++)
  {
    MFEM_VERIFY(dbc_marker[i] + farfield_marker[i] + surf_sigma_marker[i] +
                        surf_z_marker[i] + lumped_port_marker[i] + wave_port_marker[i] +
                        surf_j_marker[i] <=
                    1,
                "Boundary attributes should not be specified with multiple BC!");
  }
}

namespace
{

void PrintHeader(const mfem::ParFiniteElementSpace &h1_fespace,
                 const mfem::ParFiniteElementSpace &nd_fespace,
                 const mfem::ParFiniteElementSpace &rt_fespace, bool &print_hdr)
{
  if (print_hdr)
  {
    Mpi::Print("\nAssembling system matrices, number of global unknowns:\n"
               " H1 (p = {:d}): {:d}, ND (p = {:d}): {:d}, RT (p = {:d}): {:d}\n Operator "
               "assembly level: {}\n",
               h1_fespace.GetMaxElementOrder(), h1_fespace.GlobalTrueVSize(),
               nd_fespace.GetMaxElementOrder(), nd_fespace.GlobalTrueVSize(),
               rt_fespace.GetMaxElementOrder(), rt_fespace.GlobalTrueVSize(),
               (nd_fespace.GetMaxElementOrder() >= BilinearForm::pa_order_threshold)
                   ? "Partial"
                   : "Full");

    const auto &mesh = *nd_fespace.GetParMesh();
    const auto geom_types = mesh::CheckElements(mesh).GetGeomTypes();
    Mpi::Print(" Mesh geometries:\n");
    for (auto geom : geom_types)
    {
      const auto *fe = nd_fespace.FEColl()->FiniteElementForGeometry(geom);
      MFEM_VERIFY(fe, "MFEM does not support ND spaces on geometry = "
                          << mfem::Geometry::Name[geom] << "!");
      const int q_order = fem::DefaultIntegrationOrder::Get(mesh, geom);
      Mpi::Print("  {}: P = {:d}, Q = {:d} (quadrature order = {:d}){}\n",
                 mfem::Geometry::Name[geom], fe->GetDof(),
                 mfem::IntRules.Get(geom, q_order).GetNPoints(), q_order,
                 (geom == geom_types.back()) ? "" : ",");
    }
  }
  print_hdr = false;
}

// One PML integrator to add to a BilinearForm. Callers build these (typically one for
// curl-curl and one for mass, for each of Re/Im branches) and pass a list to
// AddIntegrators. The caller owns the lifetime of the ctx buffer; the BilinearForm
// takes a copy at libCEED operator build time via CEED_COPY_VALUES.
enum class PMLIntegKind : char
{
  CurlCurl,
  VectorFEMass
};
struct PMLIntegrator
{
  const void *ctx = nullptr;
  std::size_t ctx_size = 0;
  PMLIntegKind kind = PMLIntegKind::CurlCurl;
  bool imag_part = false;  // true ⇒ use Im(μ̃⁻¹) or Im(ε̃)
};

void AddIntegrators(BilinearForm &a, const MaterialPropertyCoefficient *df,
                    const MaterialPropertyCoefficient *f,
                    const MaterialPropertyCoefficient *dfb,
                    const MaterialPropertyCoefficient *fb,
                    const MaterialPropertyCoefficient *fp, bool assemble_q_data = false,
                    const std::vector<PMLIntegrator> *pml = nullptr)
{
  if (df && !df->empty() && f && !f->empty())
  {
    a.AddDomainIntegrator<CurlCurlMassIntegrator>(*df, *f);
  }
  else
  {
    if (df && !df->empty())
    {
      a.AddDomainIntegrator<CurlCurlIntegrator>(*df);
    }
    if (f && !f->empty())
    {
      a.AddDomainIntegrator<VectorFEMassIntegrator>(*f);
    }
  }
  if (dfb && !dfb->empty() && fb && !fb->empty())
  {
    a.AddBoundaryIntegrator<CurlCurlMassIntegrator>(*dfb, *fb);
  }
  else
  {
    if (dfb && !dfb->empty())
    {
      a.AddBoundaryIntegrator<CurlCurlIntegrator>(*dfb);
    }
    if (fb && !fb->empty())
    {
      a.AddBoundaryIntegrator<VectorFEMassIntegrator>(*fb);
    }
  }
  if (fp && !fp->empty())
  {
    a.AddDomainIntegrator<MixedVectorWeakCurlIntegrator>(*fp);
    a.AddDomainIntegrator<MixedVectorCurlIntegrator>(*fp, true);
  }
  if (pml)
  {
    for (const auto &p : *pml)
    {
      if (!p.ctx)
      {
        continue;
      }
      if (p.kind == PMLIntegKind::CurlCurl)
      {
        a.AddDomainIntegrator<CurlCurlPMLIntegrator>(p.ctx, p.ctx_size, p.imag_part);
      }
      else
      {
        a.AddDomainIntegrator<VectorFEMassPMLIntegrator>(p.ctx, p.ctx_size, p.imag_part);
      }
    }
  }
  if (assemble_q_data)
  {
    a.AssembleQuadratureData();
  }
}

void AddAuxIntegrators(BilinearForm &a, const MaterialPropertyCoefficient *f,
                       const MaterialPropertyCoefficient *fb, bool assemble_q_data = false)
{
  if (f && !f->empty())
  {
    a.AddDomainIntegrator<DiffusionIntegrator>(*f);
  }
  if (fb && !fb->empty())
  {
    a.AddBoundaryIntegrator<DiffusionIntegrator>(*fb);
  }
  if (assemble_q_data)
  {
    a.AssembleQuadratureData();
  }
}

auto AssembleOperator(const FiniteElementSpace &fespace,
                      const MaterialPropertyCoefficient *df,
                      const MaterialPropertyCoefficient *f,
                      const MaterialPropertyCoefficient *dfb,
                      const MaterialPropertyCoefficient *fb,
                      const MaterialPropertyCoefficient *fp, bool skip_zeros = false,
                      bool assemble_q_data = false,
                      const std::vector<PMLIntegrator> *pml = nullptr)
{
  BilinearForm a(fespace);
  AddIntegrators(a, df, f, dfb, fb, fp, assemble_q_data, pml);
  return a.Assemble(skip_zeros);
}

auto AssembleOperators(const FiniteElementSpaceHierarchy &fespaces,
                       const MaterialPropertyCoefficient *df,
                       const MaterialPropertyCoefficient *f,
                       const MaterialPropertyCoefficient *dfb,
                       const MaterialPropertyCoefficient *fb,
                       const MaterialPropertyCoefficient *fp, bool skip_zeros = false,
                       bool assemble_q_data = false, std::size_t l0 = 0,
                       const std::vector<PMLIntegrator> *pml = nullptr)
{
  BilinearForm a(fespaces.GetFinestFESpace());
  AddIntegrators(a, df, f, dfb, fb, fp, assemble_q_data, pml);
  return a.Assemble(fespaces, skip_zeros, l0);
}

auto AssembleAuxOperators(const FiniteElementSpaceHierarchy &fespaces,
                          const MaterialPropertyCoefficient *f,
                          const MaterialPropertyCoefficient *fb, bool skip_zeros = false,
                          bool assemble_q_data = false, std::size_t l0 = 0)
{
  BilinearForm a(fespaces.GetFinestFESpace());
  AddAuxIntegrators(a, f, fb, assemble_q_data);
  return a.Assemble(fespaces, skip_zeros, l0);
}

// Short-lived holder for a (possibly-scaled, possibly-filtered) copy of the
// MaterialOperator's packed PML QFunction context buffer. Used to feed
// CurlCurlPMLIntegrator / VectorFEMassPMLIntegrator registrations: the context is copied
// at libCEED operator build time, so different scales (e.g. a0.real() vs −a0.imag() in
// the preconditioner) need distinct buffers.
struct PMLContextBuffer
{
  std::vector<CeedIntScalar> buf;
  bool valid = false;
};

// Which PML regions to include. The MaterialOperator packs all PML profiles into one
// context; we filter at SpaceOperator level so that:
//   - GetStiffnessMatrix/GetMassMatrix include only FIXED/CFS profiles (their ω₀ is
//     baked into the context at setup).
//   - GetExtraSystemMatrix(ω) includes only FREQUENCY_DEPENDENT profiles, and refreshes
//     their ω field to the live solve frequency before assembly.
// Filtering works by zeroing the attr→profile entries for excluded profiles: the
// QFunction falls through to coeff = 0 for those attributes.
enum class PMLFilter : char
{
  All,
  FixedCFS,            // profiles with formulation == FIXED or CFS
  FrequencyDependent,  // profiles with formulation == FREQUENCY_DEPENDENT
};

PMLContextBuffer MakePMLContext(const MaterialOperator &mat_op, double scale,
                                PMLFilter filter = PMLFilter::All)
{
  PMLContextBuffer out;
  if (!mat_op.HasPML() || scale == 0.0)
  {
    return out;
  }
  if (filter == PMLFilter::FrequencyDependent && !mat_op.HasFrequencyDependentPML())
  {
    return out;
  }
  const CeedIntScalar *src = mat_op.GetPMLContextData();
  const std::size_t n = mat_op.GetPMLContextSize() / sizeof(CeedIntScalar);
  out.buf.assign(src, src + n);
  pml::SetPMLContextScale(out.buf.data(), scale);

  // Apply filter: null out attr→profile entries whose profile doesn't match.
  if (filter != PMLFilter::All)
  {
    const auto &profiles = mat_op.GetPMLProfiles();
    const auto &attr_to_profile = mat_op.GetPMLAttrToProfile();
    const int num_attr = static_cast<int>(attr_to_profile.size());
    auto keep_profile = [&](int pidx)
    {
      if (pidx < 0 || pidx >= static_cast<int>(profiles.size()))
      {
        return false;
      }
      const auto f = profiles[pidx].formulation;
      const bool is_fd = (f == PMLStretchFormulation::FREQUENCY_DEPENDENT);
      return (filter == PMLFilter::FrequencyDependent) ? is_fd : !is_fd;
    };
    bool any = false;
    for (int k = 0; k < num_attr; k++)
    {
      const int pidx = attr_to_profile[k];
      if (keep_profile(pidx))
      {
        any = true;
      }
      else
      {
        // ctx[2 + k] is the attr→profile entry (see PackProfileContextAll).
        out.buf[2 + k].first = -1;
      }
    }
    if (!any)
    {
      return PMLContextBuffer{};
    }
  }
  out.valid = true;
  return out;
}

// Convenience: build a PMLIntegrator list for the curl-curl and/or mass parts at the
// given real/imag branch, using a single scaled context. Returns empty when there's no
// PML or when scale==0. The context buffer must outlive the BilinearForm::Assemble()
// call; callers typically hold a local PMLContextBuffer and pass `spec.buf.data()`.
std::vector<PMLIntegrator> MakePMLIntegrators(const PMLContextBuffer &ctx_buf,
                                              bool include_curlcurl, bool include_mass,
                                              bool imag_part)
{
  std::vector<PMLIntegrator> out;
  if (!ctx_buf.valid)
  {
    return out;
  }
  const std::size_t bytes = ctx_buf.buf.size() * sizeof(CeedIntScalar);
  if (include_curlcurl)
  {
    out.push_back({ctx_buf.buf.data(), bytes, PMLIntegKind::CurlCurl, imag_part});
  }
  if (include_mass)
  {
    out.push_back({ctx_buf.buf.data(), bytes, PMLIntegKind::VectorFEMass, imag_part});
  }
  return out;
}

}  // namespace

template <typename OperType>
std::unique_ptr<OperType>
SpaceOperator::GetStiffnessMatrix(Operator::DiagonalPolicy diag_policy)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient df(mat_op.MaxCeedAttribute()), f(mat_op.MaxCeedAttribute()),
      fb(mat_op.MaxCeedBdrAttribute()), fc(mat_op.MaxCeedAttribute());
  AddStiffnessCoefficients(1.0, df, f);
  AddStiffnessBdrCoefficients(1.0, fb);
  AddRealPeriodicCoefficients(1.0, f);
  AddImagPeriodicCoefficients(1.0, fc);

  // PML integrators. Real branch: +1 · Re(μ̃⁻¹) curl-curl. Imag branch: +1 · Im(μ̃⁻¹)
  // curl-curl. Only FIXED/CFS profiles contribute here (their reference ω₀ is baked
  // into the context at setup). FREQUENCY_DEPENDENT profiles contribute only through
  // GetExtraSystemMatrix(ω) where ω is the live solve frequency.
  auto pml_ctx_re = MakePMLContext(mat_op, 1.0, PMLFilter::FixedCFS);
  auto pml_ctx_im = MakePMLContext(mat_op, 1.0, PMLFilter::FixedCFS);
  auto pml_re = MakePMLIntegrators(pml_ctx_re, /*curl=*/true, /*mass=*/false, false);
  auto pml_im = MakePMLIntegrators(pml_ctx_im, /*curl=*/true, /*mass=*/false, true);

  int empty[2] = {(df.empty() && f.empty() && fb.empty() && pml_re.empty()),
                  (fc.empty() && pml_im.empty())};
  Mpi::GlobalMin(2, empty, GetComm());
  if (empty[0] && empty[1])
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  std::unique_ptr<Operator> kr, ki;
  if (!empty[0])
  {
    kr = AssembleOperator(GetNDSpace(), &df, &f, nullptr, &fb, nullptr, skip_zeros, false,
                          &pml_re);
  }
  if (!empty[1])
  {
    // Imaginary branch: periodic imaginary coefficient (fc) plus PML imag μ̃⁻¹ curl-curl.
    ki = AssembleOperator(GetNDSpace(), nullptr, nullptr, nullptr, nullptr, &fc, skip_zeros,
                          false, &pml_im);
  }
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    auto K =
        std::make_unique<ComplexParOperator>(std::move(kr), std::move(ki), GetNDSpace());
    K->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return K;
  }
  else
  {
    MFEM_VERIFY(!ki, "Unexpected imaginary part in GetStiffnessMatrix<Operator>!");
    auto K = std::make_unique<ParOperator>(std::move(kr), GetNDSpace());
    K->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return K;
  }
}

template <typename OperType>
std::unique_ptr<OperType>
SpaceOperator::GetDampingMatrix(Operator::DiagonalPolicy diag_policy)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient f(mat_op.MaxCeedAttribute()),
      fb(mat_op.MaxCeedBdrAttribute());
  AddDampingCoefficients(1.0, f);
  AddDampingBdrCoefficients(1.0, fb);
  int empty = (f.empty() && fb.empty());
  Mpi::GlobalMin(1, &empty, GetComm());
  if (empty)
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  auto c = AssembleOperator(GetNDSpace(), nullptr, &f, nullptr, &fb, nullptr, skip_zeros);
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    auto C = std::make_unique<ComplexParOperator>(std::move(c), nullptr, GetNDSpace());
    C->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return C;
  }
  else
  {
    auto C = std::make_unique<ParOperator>(std::move(c), GetNDSpace());
    C->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return C;
  }
}

template <typename OperType>
std::unique_ptr<OperType> SpaceOperator::GetMassMatrix(Operator::DiagonalPolicy diag_policy)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient fr(mat_op.MaxCeedAttribute()), fi(mat_op.MaxCeedAttribute()),
      fbr(mat_op.MaxCeedBdrAttribute()), fbi(mat_op.MaxCeedBdrAttribute());
  AddRealMassCoefficients(1.0, fr);
  AddRealMassBdrCoefficients(1.0, fbr);
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    AddImagMassCoefficients(1.0, fi);
  }
  // PML ε̃ contributions for FIXED/CFS profiles. The frequency-dependent ones go into
  // GetExtraSystemMatrix(ω).
  PMLContextBuffer pml_ctx_re, pml_ctx_im;
  std::vector<PMLIntegrator> pml_re, pml_im;
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    pml_ctx_re = MakePMLContext(mat_op, 1.0, PMLFilter::FixedCFS);
    pml_ctx_im = MakePMLContext(mat_op, 1.0, PMLFilter::FixedCFS);
    pml_re = MakePMLIntegrators(pml_ctx_re, /*curl=*/false, /*mass=*/true, false);
    pml_im = MakePMLIntegrators(pml_ctx_im, /*curl=*/false, /*mass=*/true, true);
  }
  int empty[2] = {(fr.empty() && fbr.empty() && pml_re.empty()),
                  (fi.empty() && fbi.empty() && pml_im.empty())};
  Mpi::GlobalMin(2, empty, GetComm());
  if (empty[0] && empty[1])
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  std::unique_ptr<Operator> mr, mi;
  if (!empty[0])
  {
    mr = AssembleOperator(GetNDSpace(), nullptr, &fr, nullptr, &fbr, nullptr, skip_zeros,
                          false, &pml_re);
  }
  if (!empty[1])
  {
    mi = AssembleOperator(GetNDSpace(), nullptr, &fi, nullptr, &fbi, nullptr, skip_zeros,
                          false, &pml_im);
  }
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    auto M =
        std::make_unique<ComplexParOperator>(std::move(mr), std::move(mi), GetNDSpace());
    M->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M;
  }
  else
  {
    auto M = std::make_unique<ParOperator>(std::move(mr), GetNDSpace());
    M->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M;
  }
}

template <typename OperType>
std::unique_ptr<OperType>
SpaceOperator::GetExtraSystemMatrix(double omega, Operator::DiagonalPolicy diag_policy)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient dfbr(mat_op.MaxCeedBdrAttribute()),
      dfbi(mat_op.MaxCeedBdrAttribute()), fbr(mat_op.MaxCeedBdrAttribute()),
      fbi(mat_op.MaxCeedBdrAttribute());
  AddExtraSystemBdrCoefficients(omega, dfbr, dfbi, fbr, fbi);

  // FREQUENCY_DEPENDENT PML contributions. Refresh the live ω in the (shared) FD
  // regions of the master context, then build branch-local, scale-local copies. The
  // system composition is A = K + iω C − ω² M + A2. For FD regions, μ̃⁻¹ lives in A2
  // with prefactor +1 (curl-curl); ε̃ lives in A2 with prefactor −ω² (mass).
  mat_op.RefreshPMLContextFrequency(omega);
  PMLContextBuffer pml_ctx_curl_re, pml_ctx_mass_re, pml_ctx_curl_im, pml_ctx_mass_im;
  std::vector<PMLIntegrator> pml_re, pml_im;
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    pml_ctx_curl_re = MakePMLContext(mat_op, 1.0, PMLFilter::FrequencyDependent);
    pml_ctx_mass_re = MakePMLContext(mat_op, -omega * omega, PMLFilter::FrequencyDependent);
    pml_ctx_curl_im = MakePMLContext(mat_op, 1.0, PMLFilter::FrequencyDependent);
    pml_ctx_mass_im = MakePMLContext(mat_op, -omega * omega, PMLFilter::FrequencyDependent);
    auto re_curl = MakePMLIntegrators(pml_ctx_curl_re, true, false, false);
    auto re_mass = MakePMLIntegrators(pml_ctx_mass_re, false, true, false);
    auto im_curl = MakePMLIntegrators(pml_ctx_curl_im, true, false, true);
    auto im_mass = MakePMLIntegrators(pml_ctx_mass_im, false, true, true);
    pml_re.insert(pml_re.end(), re_curl.begin(), re_curl.end());
    pml_re.insert(pml_re.end(), re_mass.begin(), re_mass.end());
    pml_im.insert(pml_im.end(), im_curl.begin(), im_curl.end());
    pml_im.insert(pml_im.end(), im_mass.begin(), im_mass.end());
  }

  int empty[2] = {(dfbr.empty() && fbr.empty() && pml_re.empty()),
                  (dfbi.empty() && fbi.empty() && pml_im.empty())};
  Mpi::GlobalMin(2, empty, GetComm());
  if (empty[0] && empty[1])
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  std::unique_ptr<Operator> ar, ai;
  if (!empty[0])
  {
    ar = AssembleOperator(GetNDSpace(), nullptr, nullptr, &dfbr, &fbr, nullptr, skip_zeros,
                          false, &pml_re);
  }
  if (!empty[1])
  {
    ai = AssembleOperator(GetNDSpace(), nullptr, nullptr, &dfbi, &fbi, nullptr, skip_zeros,
                          false, &pml_im);
  }
  if constexpr (std::is_same<OperType, ComplexOperator>::value)
  {
    auto A =
        std::make_unique<ComplexParOperator>(std::move(ar), std::move(ai), GetNDSpace());
    A->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return A;
  }
  else
  {
    MFEM_VERIFY(!ai, "Unexpected imaginary part in GetExtraSystemMatrix<Operator>!");
    auto A = std::make_unique<ParOperator>(std::move(ar), GetNDSpace());
    A->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return A;
  }
}

template <typename OperType, typename ScalarType>
std::unique_ptr<OperType>
SpaceOperator::GetSystemMatrix(ScalarType a0, ScalarType a1, ScalarType a2,
                               const OperType *K, const OperType *C, const OperType *M,
                               const OperType *A2)
{
  return BuildParSumOperator({a0, a1, a2, ScalarType{1}}, {K, C, M, A2});
}

std::unique_ptr<Operator> SpaceOperator::GetInnerProductMatrix(double a0, double a2,
                                                               const ComplexOperator *K,
                                                               const ComplexOperator *M)
{
  const auto *PtAP_K = (K) ? dynamic_cast<const ComplexParOperator *>(K) : nullptr;
  const auto *PtAP_M = (M) ? dynamic_cast<const ComplexParOperator *>(M) : nullptr;
  return BuildParSumOperator(
      {a0, a2}, {PtAP_K ? PtAP_K->Real() : nullptr, PtAP_M ? PtAP_M->Real() : nullptr});
}

namespace
{

template <typename OperType>
auto BuildLevelParOperator(std::unique_ptr<Operator> &&br, std::unique_ptr<Operator> &&bi,
                           const FiniteElementSpace &fespace);

template <>
auto BuildLevelParOperator<Operator>(std::unique_ptr<Operator> &&br,
                                     std::unique_ptr<Operator> &&bi,
                                     const FiniteElementSpace &fespace)
{
  MFEM_VERIFY(
      !bi,
      "Should not be constructing a real-valued ParOperator with non-zero imaginary part!");
  return std::make_unique<ParOperator>(std::move(br), fespace);
}

template <>
auto BuildLevelParOperator<ComplexOperator>(std::unique_ptr<Operator> &&br,
                                            std::unique_ptr<Operator> &&bi,
                                            const FiniteElementSpace &fespace)
{
  return std::make_unique<ComplexParOperator>(std::move(br), std::move(bi), fespace);
}

}  // namespace

void SpaceOperator::AssemblePreconditioner(
    std::complex<double> a0, std::complex<double> a1, std::complex<double> a2, double a3,
    std::vector<std::unique_ptr<Operator>> &br_vec,
    std::vector<std::unique_ptr<Operator>> &br_aux_vec,
    std::vector<std::unique_ptr<Operator>> &bi_vec,
    std::vector<std::unique_ptr<Operator>> &bi_aux_vec)
{
  constexpr bool skip_zeros = false, assemble_q_data = false;
  MaterialPropertyCoefficient dfr(mat_op.MaxCeedAttribute()),
      dfi(mat_op.MaxCeedAttribute()), fr(mat_op.MaxCeedAttribute()),
      fi(mat_op.MaxCeedAttribute()), dfbr(mat_op.MaxCeedBdrAttribute()),
      dfbi(mat_op.MaxCeedBdrAttribute()), fbr(mat_op.MaxCeedBdrAttribute()),
      fbi(mat_op.MaxCeedBdrAttribute()), fpi(mat_op.MaxCeedAttribute()),
      fpr(mat_op.MaxCeedAttribute());
  AddStiffnessCoefficients(a0.real(), dfr, fr);
  AddStiffnessCoefficients(a0.imag(), dfi, fi);
  AddStiffnessBdrCoefficients(a0.real(), fbr);
  AddStiffnessBdrCoefficients(a0.imag(), fbi);
  AddDampingCoefficients(a1.real(), fr);
  AddDampingCoefficients(a1.imag(), fi);
  AddDampingBdrCoefficients(a1.real(), fbr);
  AddDampingBdrCoefficients(a1.imag(), fbi);
  AddRealMassCoefficients(pc_mat_shifted ? std::abs(a2.real()) : a2.real(), fr);
  AddRealMassCoefficients(a2.imag(), fi);
  AddRealMassBdrCoefficients(pc_mat_shifted ? std::abs(a2.real()) : a2.real(), fbr);
  AddRealMassBdrCoefficients(a2.imag(), fbi);
  AddImagMassCoefficients(a2.real(), fi);
  AddImagMassCoefficients(-a2.imag(), fr);
  AddExtraSystemBdrCoefficients(a3, dfbr, dfbi, fbr, fbi);
  AddRealPeriodicCoefficients(a0.real(), fr);
  AddRealPeriodicCoefficients(a0.imag(), fi);
  AddImagPeriodicCoefficients(a0.real(), fpi);
  AddImagPeriodicCoefficients(-a0.imag(), fpr);

  // PML integrators for the complex-valued preconditioner. The preconditioner mirrors
  // the system matrix A = a0 K + a1 C + a2 M + A2(a3). For PML μ̃⁻¹ = μr + i μi and
  // ε̃ = εr + i εi:
  //   a0·μ̃⁻¹ · (u,v)curl has real part a0r·μr − a0i·μi and imag part a0r·μi + a0i·μr
  //   a2·ε̃ · (u,v)      has real part a2r·εr − a2i·εi and imag part a2r·εi + a2i·εr
  // For the shifted variant, only the real part of a2 gets |·|; the imag-branch cross
  // term (−a2i·εr in the real branch) is shaped by the real part of ε̃ so we match the
  // mass helpers above (which gate real-part scaling on pc_mat_shifted). FIXED/CFS
  // profiles contribute at all frequencies here; FREQUENCY_DEPENDENT ones use ω = a3.
  mat_op.RefreshPMLContextFrequency(a3);
  const double a2r = pc_mat_shifted ? std::abs(a2.real()) : a2.real();
  const double a2i = a2.imag();
  PMLContextBuffer fc_re_mu_re, fc_re_mu_im, fc_im_mu_re, fc_im_mu_im;
  PMLContextBuffer fc_re_ep_re, fc_re_ep_im, fc_im_ep_re, fc_im_ep_im;
  PMLContextBuffer fd_re_mu_re, fd_re_mu_im, fd_im_mu_re, fd_im_mu_im;
  PMLContextBuffer fd_re_ep_re, fd_re_ep_im, fd_im_ep_re, fd_im_ep_im;
  fc_re_mu_re = MakePMLContext(mat_op, a0.real(), PMLFilter::FixedCFS);
  fc_re_mu_im = MakePMLContext(mat_op, -a0.imag(), PMLFilter::FixedCFS);
  fc_im_mu_re = MakePMLContext(mat_op, a0.imag(), PMLFilter::FixedCFS);
  fc_im_mu_im = MakePMLContext(mat_op, a0.real(), PMLFilter::FixedCFS);
  fc_re_ep_re = MakePMLContext(mat_op, a2r, PMLFilter::FixedCFS);
  fc_re_ep_im = MakePMLContext(mat_op, -a2i, PMLFilter::FixedCFS);
  fc_im_ep_re = MakePMLContext(mat_op, a2i, PMLFilter::FixedCFS);
  fc_im_ep_im = MakePMLContext(mat_op, a2r, PMLFilter::FixedCFS);
  fd_re_mu_re = MakePMLContext(mat_op, 1.0, PMLFilter::FrequencyDependent);
  fd_im_mu_im = MakePMLContext(mat_op, 1.0, PMLFilter::FrequencyDependent);
  fd_re_ep_re = MakePMLContext(mat_op, -a3 * a3, PMLFilter::FrequencyDependent);
  fd_im_ep_im = MakePMLContext(mat_op, -a3 * a3, PMLFilter::FrequencyDependent);
  // Real branch integrators: curl-curl (Re(μ̃⁻¹) with +a0r, Im(μ̃⁻¹) with −a0i) + mass
  // (Re(ε̃) with a2r, Im(ε̃) with −a2i) + FD cross-terms at ω=a3.
  std::vector<PMLIntegrator> pml_re, pml_im;
  auto append = [](std::vector<PMLIntegrator> &dst, std::vector<PMLIntegrator> &&src)
  {
    for (auto &v : src)
    {
      dst.push_back(v);
    }
  };
  append(pml_re, MakePMLIntegrators(fc_re_mu_re, true, false, false));
  append(pml_re, MakePMLIntegrators(fc_re_mu_im, true, false, true));
  append(pml_re, MakePMLIntegrators(fc_re_ep_re, false, true, false));
  append(pml_re, MakePMLIntegrators(fc_re_ep_im, false, true, true));
  append(pml_re, MakePMLIntegrators(fd_re_mu_re, true, false, false));
  append(pml_re, MakePMLIntegrators(fd_re_ep_re, false, true, false));
  // Imag branch integrators.
  append(pml_im, MakePMLIntegrators(fc_im_mu_re, true, false, false));
  append(pml_im, MakePMLIntegrators(fc_im_mu_im, true, false, true));
  append(pml_im, MakePMLIntegrators(fc_im_ep_re, false, true, false));
  append(pml_im, MakePMLIntegrators(fc_im_ep_im, false, true, true));
  append(pml_im, MakePMLIntegrators(fd_im_mu_im, true, false, true));
  append(pml_im, MakePMLIntegrators(fd_im_ep_im, false, true, true));

  int empty[2] = {(dfr.empty() && fr.empty() && dfbr.empty() && fbr.empty() &&
                   fpr.empty() && pml_re.empty()),
                  (dfi.empty() && fi.empty() && dfbi.empty() && fbi.empty() &&
                   fpi.empty() && pml_im.empty())};
  Mpi::GlobalMin(2, empty, GetComm());
  if (!empty[0])
  {
    br_vec = AssembleOperators(GetNDSpaces(), &dfr, &fr, &dfbr, &fbr, &fpr, skip_zeros,
                               assemble_q_data, /*l0=*/0, &pml_re);
    br_aux_vec =
        AssembleAuxOperators(GetH1Spaces(), &fr, &fbr, skip_zeros, assemble_q_data);
  }
  if (!empty[1])
  {
    bi_vec = AssembleOperators(GetNDSpaces(), &dfi, &fi, &dfbi, &fbi, &fpi, skip_zeros,
                               assemble_q_data, /*l0=*/0, &pml_im);
    bi_aux_vec =
        AssembleAuxOperators(GetH1Spaces(), &fi, &fbi, skip_zeros, assemble_q_data);
  }
}

void SpaceOperator::AssemblePreconditioner(
    std::complex<double> a0, std::complex<double> a1, std::complex<double> a2, double a3,
    std::vector<std::unique_ptr<Operator>> &br_vec,
    std::vector<std::unique_ptr<Operator>> &br_aux_vec)
{
  constexpr bool skip_zeros = false, assemble_q_data = false;
  MaterialPropertyCoefficient dfr(mat_op.MaxCeedAttribute()), fr(mat_op.MaxCeedAttribute()),
      dfbr(mat_op.MaxCeedBdrAttribute()), fbr(mat_op.MaxCeedBdrAttribute());
  AddStiffnessCoefficients(a0.real(), dfr, fr);
  AddStiffnessBdrCoefficients(a0.real(), fbr);
  AddDampingCoefficients(a1.imag(), fr);
  AddDampingBdrCoefficients(a1.imag(), fbr);
  AddAbsMassCoefficients(pc_mat_shifted ? std::abs(a2.real()) : a2.real(), fr);
  AddRealMassBdrCoefficients(pc_mat_shifted ? std::abs(a2.real()) : a2.real(), fbr);
  AddExtraSystemBdrCoefficients(a3, dfbr, dfbr, fbr, fbr);
  AddRealPeriodicCoefficients(a0.real(), fr);

  // Real-valued preconditioner approximation for PML: use Re(μ̃⁻¹) for the stiffness
  // piece and Re(ε̃) for the mass piece, with the same sign/shifted rules as the bulk
  // material. FIXED/CFS regions use their reference ω₀ (baked in at setup);
  // FREQUENCY_DEPENDENT ones pick up ω = a3. Im(μ̃⁻¹) / Im(ε̃) are dropped — this is an
  // approximation and the surrounding Krylov solver absorbs the leftover.
  mat_op.RefreshPMLContextFrequency(a3);
  const double a2r = pc_mat_shifted ? std::abs(a2.real()) : a2.real();
  auto fc_mu = MakePMLContext(mat_op, a0.real(), PMLFilter::FixedCFS);
  auto fc_ep = MakePMLContext(mat_op, a2r, PMLFilter::FixedCFS);
  auto fd_mu = MakePMLContext(mat_op, 1.0, PMLFilter::FrequencyDependent);
  auto fd_ep = MakePMLContext(mat_op, -a3 * a3, PMLFilter::FrequencyDependent);
  std::vector<PMLIntegrator> pml;
  auto append = [](std::vector<PMLIntegrator> &dst, std::vector<PMLIntegrator> &&src)
  {
    for (auto &v : src)
    {
      dst.push_back(v);
    }
  };
  append(pml, MakePMLIntegrators(fc_mu, true, false, false));
  append(pml, MakePMLIntegrators(fc_ep, false, true, false));
  append(pml, MakePMLIntegrators(fd_mu, true, false, false));
  append(pml, MakePMLIntegrators(fd_ep, false, true, false));

  int empty = (dfr.empty() && fr.empty() && dfbr.empty() && fbr.empty() && pml.empty());
  Mpi::GlobalMin(1, &empty, GetComm());
  if (!empty)
  {
    br_vec = AssembleOperators(GetNDSpaces(), &dfr, &fr, &dfbr, &fbr, nullptr, skip_zeros,
                               assemble_q_data, /*l0=*/0, &pml);
    br_aux_vec =
        AssembleAuxOperators(GetH1Spaces(), &fr, &fbr, skip_zeros, assemble_q_data);
  }
}

void SpaceOperator::AssemblePreconditioner(
    double a0, double a1, double a2, double a3,
    std::vector<std::unique_ptr<Operator>> &br_vec,
    std::vector<std::unique_ptr<Operator>> &br_aux_vec)
{
  constexpr bool skip_zeros = false, assemble_q_data = false;
  MaterialPropertyCoefficient dfr(mat_op.MaxCeedAttribute()), fr(mat_op.MaxCeedAttribute()),
      dfbr(mat_op.MaxCeedBdrAttribute()), fbr(mat_op.MaxCeedBdrAttribute());
  AddStiffnessCoefficients(a0, dfr, fr);
  AddStiffnessBdrCoefficients(a0, fbr);
  AddDampingCoefficients(a1, fr);
  AddDampingBdrCoefficients(a1, fbr);
  AddAbsMassCoefficients(pc_mat_shifted ? std::abs(a2) : a2, fr);
  AddRealMassBdrCoefficients(pc_mat_shifted ? std::abs(a2) : a2, fbr);
  AddExtraSystemBdrCoefficients(a3, dfbr, dfbr, fbr, fbr);
  AddRealPeriodicCoefficients(a0, fr);

  // Real-valued PML preconditioner: Re(μ̃⁻¹) curl-curl and Re(ε̃) mass. Same filtering
  // split by formulation as the complex variant. Im parts are dropped by design.
  mat_op.RefreshPMLContextFrequency(a3);
  const double a2r = pc_mat_shifted ? std::abs(a2) : a2;
  auto fc_mu = MakePMLContext(mat_op, a0, PMLFilter::FixedCFS);
  auto fc_ep = MakePMLContext(mat_op, a2r, PMLFilter::FixedCFS);
  auto fd_mu = MakePMLContext(mat_op, 1.0, PMLFilter::FrequencyDependent);
  auto fd_ep = MakePMLContext(mat_op, -a3 * a3, PMLFilter::FrequencyDependent);
  std::vector<PMLIntegrator> pml;
  auto append = [](std::vector<PMLIntegrator> &dst, std::vector<PMLIntegrator> &&src)
  {
    for (auto &v : src)
    {
      dst.push_back(v);
    }
  };
  append(pml, MakePMLIntegrators(fc_mu, true, false, false));
  append(pml, MakePMLIntegrators(fc_ep, false, true, false));
  append(pml, MakePMLIntegrators(fd_mu, true, false, false));
  append(pml, MakePMLIntegrators(fd_ep, false, true, false));

  int empty = (dfr.empty() && fr.empty() && dfbr.empty() && fbr.empty() && pml.empty());
  Mpi::GlobalMin(1, &empty, GetComm());
  if (!empty)
  {
    br_vec = AssembleOperators(GetNDSpaces(), &dfr, &fr, &dfbr, &fbr, nullptr, skip_zeros,
                               assemble_q_data, /*l0=*/0, &pml);
    br_aux_vec =
        AssembleAuxOperators(GetH1Spaces(), &fr, &fbr, skip_zeros, assemble_q_data);
  }
}

template <typename OperType, typename ScalarType>
std::unique_ptr<OperType> SpaceOperator::GetPreconditionerMatrix(ScalarType a0,
                                                                 ScalarType a1,
                                                                 ScalarType a2, double a3)
{
  // When partially assembled, the coarse operators can reuse the fine operator quadrature
  // data if the spaces correspond to the same mesh. When appropriate, we build the
  // preconditioner on all levels based on the actual complex-valued system matrix. The
  // coarse operator is always fully assembled.
  if (print_prec_hdr)
  {
    Mpi::Print("\nAssembling multigrid hierarchy:\n");
  }
  MFEM_VERIFY(GetH1Spaces().GetNumLevels() == GetNDSpaces().GetNumLevels(),
              "Multigrid hierarchy mismatch for auxiliary space preconditioning!");

  const auto n_levels = GetNDSpaces().GetNumLevels();
  std::vector<std::unique_ptr<Operator>> br_vec(n_levels), bi_vec(n_levels),
      br_aux_vec(n_levels), bi_aux_vec(n_levels);
  if (std::is_same<OperType, ComplexOperator>::value && !pc_mat_real)
  {
    AssemblePreconditioner(a0, a1, a2, a3, br_vec, br_aux_vec, bi_vec, bi_aux_vec);
  }
  else
  {
    AssemblePreconditioner(a0, a1, a2, a3, br_vec, br_aux_vec);
  }

  auto B = std::make_unique<BaseMultigridOperator<OperType>>(n_levels);
  for (bool aux : {false, true})
  {
    for (std::size_t l = 0; l < n_levels; l++)
    {
      const auto &fespace_l =
          aux ? GetH1Spaces().GetFESpaceAtLevel(l) : GetNDSpaces().GetFESpaceAtLevel(l);
      const auto &dbc_tdof_lists_l = aux ? h1_dbc_tdof_lists[l] : nd_dbc_tdof_lists[l];
      auto &br_l = aux ? br_aux_vec[l] : br_vec[l];
      auto &bi_l = aux ? bi_aux_vec[l] : bi_vec[l];
      if (print_prec_hdr)
      {
        Mpi::Print(" Level {:d}{} (p = {:d}): {:d} unknowns", l, aux ? " (auxiliary)" : "",
                   fespace_l.GetMaxElementOrder(), fespace_l.GlobalTrueVSize());
        const auto *b_spm = dynamic_cast<const hypre::HypreCSRMatrix *>(br_l.get());
        if (!b_spm)
        {
          b_spm = dynamic_cast<const hypre::HypreCSRMatrix *>(bi_l.get());
        }
        if (b_spm)
        {
          HYPRE_BigInt nnz = b_spm->NNZ();
          Mpi::GlobalSum(1, &nnz, fespace_l.GetComm());
          Mpi::Print(", {:d} NNZ\n", nnz);
        }
        else
        {
          Mpi::Print("\n");
        }
      }
      auto B_l =
          BuildLevelParOperator<OperType>(std::move(br_l), std::move(bi_l), fespace_l);
      B_l->SetEssentialTrueDofs(dbc_tdof_lists_l, Operator::DiagonalPolicy::DIAG_ONE);
      if (aux)
      {
        B->AddAuxiliaryOperator(std::move(B_l));
      }
      else
      {
        B->AddOperator(std::move(B_l));
      }
    }
  }

  print_prec_hdr = false;
  return B;
}

void SpaceOperator::AddStiffnessCoefficients(double coeff, MaterialPropertyCoefficient &df,
                                             MaterialPropertyCoefficient &f)
{
  // Contribution from material permeability. MaterialOperator zeros out entries for PML
  // attributes, so this line contributes nothing there; PML contributions come from a
  // separate CurlCurlPMLIntegrator queued at BilinearForm-registration time.
  df.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetInvPermeability(), coeff);

  // Contribution for London superconductors.
  if (mat_op.HasLondonDepth())
  {
    f.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetInvLondonDepth(), coeff);
  }
}

void SpaceOperator::AddStiffnessBdrCoefficients(double coeff,
                                                MaterialPropertyCoefficient &fb)
{
  // Robin BC contributions due to surface impedance and lumped ports (inductance).
  surf_z_op.AddStiffnessBdrCoefficients(coeff, fb);
  lumped_port_op.AddStiffnessBdrCoefficients(coeff, fb);
}

void SpaceOperator::AddDampingCoefficients(double coeff, MaterialPropertyCoefficient &f)
{
  // Contribution for domain conductivity.
  if (mat_op.HasConductivity())
  {
    f.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetConductivity(), coeff);
  }
}

void SpaceOperator::AddDampingBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb)
{
  // Robin BC contributions due to surface impedance, lumped ports, and absorbing
  // boundaries (resistance).
  farfield_op.AddDampingBdrCoefficients(coeff, fb);
  surf_z_op.AddDampingBdrCoefficients(coeff, fb);
  lumped_port_op.AddDampingBdrCoefficients(coeff, fb);
}

void SpaceOperator::AddRealMassCoefficients(double coeff, MaterialPropertyCoefficient &f)
{
  // MaterialOperator zeros bulk ε for PML attributes; PML ε̃ is applied through a
  // VectorFEMassPMLIntegrator queued separately at BilinearForm-registration time.
  f.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetPermittivityReal(), coeff);
}

void SpaceOperator::AddRealMassBdrCoefficients(double coeff,
                                               MaterialPropertyCoefficient &fb)
{
  // Robin BC contributions due to surface impedance and lumped ports (capacitance).
  surf_z_op.AddMassBdrCoefficients(coeff, fb);
  lumped_port_op.AddMassBdrCoefficients(coeff, fb);
}

void SpaceOperator::AddImagMassCoefficients(double coeff, MaterialPropertyCoefficient &f)
{
  // Contribution for loss tangent: ε -> ε * (1 - i tan(δ)).
  if (mat_op.HasLossTangent())
  {
    f.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetPermittivityImag(), coeff);
  }
  // PML imag ε̃ is applied through a separate VectorFEMassPMLIntegrator.
}

void SpaceOperator::AddAbsMassCoefficients(double coeff, MaterialPropertyCoefficient &f)
{
  f.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetPermittivityAbs(), coeff);
}

void SpaceOperator::AddExtraSystemBdrCoefficients(double omega,
                                                  MaterialPropertyCoefficient &dfbr,
                                                  MaterialPropertyCoefficient &dfbi,
                                                  MaterialPropertyCoefficient &fbr,
                                                  MaterialPropertyCoefficient &fbi)
{
  // Contribution for second-order farfield boundaries and finite conductivity boundaries.
  farfield_op.AddExtraSystemBdrCoefficients(omega, dfbr, dfbi);
  surf_sigma_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);

  // Contribution for numeric wave ports.
  wave_port_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);
}

void SpaceOperator::AddRealPeriodicCoefficients(double coeff,
                                                MaterialPropertyCoefficient &f)
{
  // Floquet periodicity contributions.
  if (mat_op.HasWaveVector())
  {
    f.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetFloquetMass(), coeff);
  }
}

void SpaceOperator::AddImagPeriodicCoefficients(double coeff,
                                                MaterialPropertyCoefficient &f)
{
  // Floquet periodicity contributions.
  if (mat_op.HasWaveVector())
  {
    f.AddCoefficient(mat_op.GetAttributeToMaterial(), mat_op.GetFloquetCurl(), coeff);
  }
}

bool SpaceOperator::GetExcitationVector(int excitation_idx, Vector &RHS)
{
  // Time domain excitation vector.
  RHS.SetSize(GetNDSpace().GetTrueVSize());
  RHS.UseDevice(true);
  RHS = 0.0;
  bool nnz = AddExcitationVector1Internal(excitation_idx, RHS);
  linalg::SetSubVector(RHS, nd_dbc_tdof_lists.back(), 0.0);
  return nnz;
}

bool SpaceOperator::GetExcitationVector(int excitation_idx, double omega,
                                        ComplexVector &RHS)
{
  // Frequency domain excitation vector: RHS = iω RHS1 + RHS2(ω).
  RHS.SetSize(GetNDSpace().GetTrueVSize());
  RHS.UseDevice(true);
  RHS = 0.0;
  bool nnz1 = AddExcitationVector1Internal(excitation_idx, RHS.Real());
  RHS *= 1i * omega;
  bool nnz2 = AddExcitationVector2Internal(excitation_idx, omega, RHS);
  linalg::SetSubVector(RHS, nd_dbc_tdof_lists.back(), 0.0);
  return nnz1 || nnz2;
}

void SpaceOperator::GetLumpedPortExcitationVectorPrimaryEt(int port_idx,
                                                           ComplexVector &Et_primary,
                                                           bool zero_metal)
{
  const auto &data = GetLumpedPortOp().GetPort(port_idx);

  SumVectorCoefficient fb(GetMesh().SpaceDimension());
  mfem::Array<int> attr_list;
  for (const auto &elem : data.elems)
  {
    const double Rs = 1.0 * data.GetToSquare(*elem);
    const double Einc = std::sqrt(
        Rs / (elem->GetGeometryWidth() * elem->GetGeometryLength() * data.elems.size()));
    fb.AddCoefficient(elem->GetModeCoefficient(Einc));
    attr_list.Append(elem->GetAttrList());
  }

  Et_primary.SetSize(GetNDSpace().GetTrueVSize());
  Et_primary.UseDevice(true);
  Et_primary = 0.0;

  const auto &mesh = GetNDSpace().GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker;
  mesh::AttrToMarker(bdr_attr_max, attr_list, attr_marker);

  GridFunction rhs(GetNDSpace());
  rhs = 0.0;
  rhs.Real().ProjectBdrCoefficientTangent(fb, attr_marker);
  GetNDSpace().GetRestrictionMatrix()->Mult(rhs.Real(), Et_primary.Real());

  if (zero_metal)
  {
    linalg::SetSubVector(Et_primary.Real(), GetNDDbcTDofLists().back(), 0.0);
  }
}

void SpaceOperator::GetLumpedPortExcitationVectorPrimaryHtcn(int port_idx,
                                                             ComplexVector &Htcn_primary,
                                                             bool zero_metal)
{
  const auto &data = lumped_port_op.GetPort(port_idx);

  SumVectorCoefficient fb(GetMesh().SpaceDimension());
  mfem::Array<int> attr_list;
  for (const auto &elem : data.elems)
  {
    const double Rs = 1.0 * data.GetToSquare(*elem);
    const double Hinc = 1.0 / std::sqrt(Rs * elem->GetGeometryWidth() *
                                        elem->GetGeometryLength() * data.elems.size());
    fb.AddCoefficient(elem->GetModeCoefficient(Hinc));
    attr_list.Append(elem->GetAttrList());
  }

  Htcn_primary.SetSize(GetNDSpace().GetTrueVSize());
  Htcn_primary.UseDevice(true);
  Htcn_primary = 0.0;

  const auto &mesh = GetNDSpace().GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker;
  mesh::AttrToMarker(bdr_attr_max, attr_list, attr_marker);

  GridFunction rhs(GetNDSpace());
  rhs = 0.0;
  rhs.Real().ProjectBdrCoefficientTangent(fb, attr_marker);
  GetNDSpace().GetRestrictionMatrix()->Mult(rhs.Real(), Htcn_primary.Real());

  if (zero_metal)
  {
    linalg::SetSubVector(Htcn_primary.Real(), GetNDDbcTDofLists().back(), 0.0);
  }
}

bool SpaceOperator::GetExcitationVector1(int excitation_idx, ComplexVector &RHS1)
{
  // Assemble the frequency domain excitation term with linear frequency dependence
  // (coefficient iω, see GetExcitationVector above, is accounted for later).
  RHS1.SetSize(GetNDSpace().GetTrueVSize());
  RHS1.UseDevice(true);
  RHS1 = 0.0;
  bool nnz1 = AddExcitationVector1Internal(excitation_idx, RHS1.Real());
  linalg::SetSubVector(RHS1.Real(), nd_dbc_tdof_lists.back(), 0.0);
  return nnz1;
}

bool SpaceOperator::GetExcitationVector2(int excitation_idx, double omega,
                                         ComplexVector &RHS2)
{
  RHS2.SetSize(GetNDSpace().GetTrueVSize());
  RHS2.UseDevice(true);
  RHS2 = 0.0;
  bool nnz2 = AddExcitationVector2Internal(excitation_idx, omega, RHS2);
  linalg::SetSubVector(RHS2, nd_dbc_tdof_lists.back(), 0.0);
  return nnz2;
}

bool SpaceOperator::AddExcitationVector1Internal(int excitation_idx, Vector &RHS1)
{
  // Assemble the time domain excitation -g'(t) J or frequency domain excitation -iω J.
  // The g'(t) or iω factors are not accounted for here, they are accounted for in the time
  // integration or frequency sweep later.
  MFEM_VERIFY(RHS1.Size() == GetNDSpace().GetTrueVSize(),
              "Invalid T-vector size for AddExcitationVector1Internal!");

  // Boundary sources
  SumVectorCoefficient fb(GetMesh().SpaceDimension());
  lumped_port_op.AddExcitationBdrCoefficients(excitation_idx, fb);
  surf_j_op.AddExcitationBdrCoefficients(fb);  // No excitation_idx — currently in all

  // Domain sources (current dipoles) - use integrator-based approach
  bool has_current_dipoles = !current_dipole_op.Empty();

  int empty[2] = {(fb.empty()), (!has_current_dipoles)};
  Mpi::GlobalMin(2, empty, GetComm());
  if (empty[0] && empty[1])
  {
    return false;
  }

  mfem::LinearForm rhs1(&GetNDSpace().Get());

  // Add boundary integrators
  if (!empty[0])
  {
    rhs1.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fb));
  }

  // Add domain integrators for current dipoles
  if (!empty[1])
  {
    current_dipole_op.AddExcitationDomainIntegrators(rhs1);
  }

  rhs1.UseFastAssembly(false);
  rhs1.UseDevice(false);
  rhs1.Assemble();
  rhs1.UseDevice(true);
  GetNDSpace().GetProlongationMatrix()->AddMultTranspose(rhs1, RHS1);
  return true;
}

bool SpaceOperator::AddExcitationVector2Internal(int excitation_idx, double omega,
                                                 ComplexVector &RHS2)
{
  // Assemble the contribution of wave ports to the frequency domain excitation term at
  // the specified frequency.
  MFEM_VERIFY(RHS2.Size() == GetNDSpace().GetTrueVSize(),
              "Invalid T-vector size for AddExcitationVector2Internal!");
  SumVectorCoefficient fbr(GetMesh().SpaceDimension()), fbi(GetMesh().SpaceDimension());
  wave_port_op.AddExcitationBdrCoefficients(excitation_idx, omega, fbr, fbi);
  int empty = (fbr.empty() && fbi.empty());
  Mpi::GlobalMin(1, &empty, GetComm());
  if (empty)
  {
    return false;
  }
  {
    mfem::LinearForm rhs2(&GetNDSpace().Get());
    rhs2.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbr));
    rhs2.UseFastAssembly(false);
    rhs2.UseDevice(false);
    rhs2.Assemble();
    rhs2.UseDevice(true);
    GetNDSpace().GetProlongationMatrix()->AddMultTranspose(rhs2, RHS2.Real());
  }
  {
    mfem::LinearForm rhs2(&GetNDSpace().Get());
    rhs2.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbi));
    rhs2.UseFastAssembly(false);
    rhs2.UseDevice(false);
    rhs2.Assemble();
    rhs2.UseDevice(true);
    GetNDSpace().GetProlongationMatrix()->AddMultTranspose(rhs2, RHS2.Imag());
  }
  return true;
}

void SpaceOperator::GetConstantInitialVector(ComplexVector &v)
{
  v.SetSize(GetNDSpace().GetTrueVSize());
  v.UseDevice(true);
  v = 1.0;
  linalg::SetSubVector(v.Real(), nd_dbc_tdof_lists.back(), 0.0);
}

void SpaceOperator::GetRandomInitialVector(ComplexVector &v)
{
  v.SetSize(GetNDSpace().GetTrueVSize());
  v.UseDevice(true);
  linalg::SetRandom(GetNDSpace().GetComm(), v);
  linalg::SetSubVector(v, nd_dbc_tdof_lists.back(), 0.0);
}

template std::unique_ptr<Operator>
    SpaceOperator::GetStiffnessMatrix(Operator::DiagonalPolicy);
template std::unique_ptr<ComplexOperator>
    SpaceOperator::GetStiffnessMatrix(Operator::DiagonalPolicy);

template std::unique_ptr<Operator>
    SpaceOperator::GetDampingMatrix(Operator::DiagonalPolicy);
template std::unique_ptr<ComplexOperator>
    SpaceOperator::GetDampingMatrix(Operator::DiagonalPolicy);

template std::unique_ptr<Operator> SpaceOperator::GetMassMatrix(Operator::DiagonalPolicy);
template std::unique_ptr<ComplexOperator>
    SpaceOperator::GetMassMatrix(Operator::DiagonalPolicy);

template std::unique_ptr<Operator>
SpaceOperator::GetExtraSystemMatrix(double, Operator::DiagonalPolicy);
template std::unique_ptr<ComplexOperator>
SpaceOperator::GetExtraSystemMatrix(double, Operator::DiagonalPolicy);

template std::unique_ptr<Operator>
SpaceOperator::GetSystemMatrix<Operator, double>(double, double, double, const Operator *,
                                                 const Operator *, const Operator *,
                                                 const Operator *);
template std::unique_ptr<ComplexOperator>
SpaceOperator::GetSystemMatrix<ComplexOperator, std::complex<double>>(
    std::complex<double>, std::complex<double>, std::complex<double>,
    const ComplexOperator *, const ComplexOperator *, const ComplexOperator *,
    const ComplexOperator *);

template std::unique_ptr<Operator>
SpaceOperator::GetPreconditionerMatrix<Operator, double>(double, double, double, double);
template std::unique_ptr<ComplexOperator>
SpaceOperator::GetPreconditionerMatrix<ComplexOperator, std::complex<double>>(
    std::complex<double>, std::complex<double>, std::complex<double>, double);

}  // namespace palace
