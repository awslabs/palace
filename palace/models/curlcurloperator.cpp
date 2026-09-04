// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "curlcurloperator.hpp"

#include <algorithm>
#include <set>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/multigrid.hpp"
#include "linalg/hypre.hpp"
#include "linalg/rap.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

namespace
{

// ±1 aligning the boundary-face normal with the flux loop direction. A hole aperture is an
// *internal* surface, whose face orientations the mesh need not make consistent (EnsureNCMesh
// scrambles them, and some generators emit them mixed), which would otherwise make c measure
// only part of the aperture. Matches ComputeFluxThroughSurface's DIRECTION_BASED convention.
class FluxOrientationCoefficient : public mfem::Coefficient
{
private:
  const mfem::Vector dir;
  mutable mfem::Vector nor;

public:
  FluxOrientationCoefficient(const std::vector<double> &d, int sdim)
    : dir(const_cast<double *>(d.data()), sdim), nor(sdim)
  {
    MFEM_VERIFY(static_cast<int>(d.size()) >= sdim,
                "London flux loop direction must have " << sdim << " components!");
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    T.SetIntPoint(&ip);
    mfem::CalcOrtho(T.Jacobian(), nor);
    return (nor * dir < 0.0) ? -1.0 : 1.0;
  }
};

}  // namespace

CurlCurlOperator::CurlCurlOperator(const config::BoundaryData &boundaries,
                                   const config::SolverData &solver,
                                   const std::vector<config::MaterialData> &materials,
                                   ProblemType problem_type, const Units &units,
                                   const std::vector<std::unique_ptr<Mesh>> &mesh)
  : print_hdr(true),
    dbc_attr(SetUpBoundaryProperties(boundaries.pec, boundaries.fluxloop,
                                     boundaries.superconductor, *mesh.back())),
    nd_fecs(fem::ConstructFECollections<mfem::ND_FECollection>(
        solver.order, mesh.back()->Dimension(), solver.linear.mg_max_levels,
        solver.linear.mg_coarsening, false)),
    h1_fecs(fem::ConstructFECollections<mfem::H1_FECollection>(
        solver.order, mesh.back()->Dimension(), solver.linear.mg_max_levels,
        solver.linear.mg_coarsening, false)),
    rt_fec(std::make_unique<mfem::RT_FECollection>(solver.order - 1,
                                                   mesh.back()->Dimension())),
    nd_fespaces(fem::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
        solver.linear.mg_max_levels, mesh, nd_fecs, &dbc_attr, &dbc_tdof_lists)),
    h1_fespaces(fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
        solver.linear.mg_max_levels, mesh, h1_fecs)),
    rt_fespace(*mesh.back(), rt_fec.get()),
    mat_op(materials, boundaries.periodic, problem_type, *mesh.back()),
    surf_j_op(boundaries.current, *mesh.back()),
    sc_sheet_op(boundaries.superconductor, boundaries.cracked_attributes, units, mat_op,
                *mesh.back()),
    london_pc_shift_(solver.linear.london_pc_shift)
{
  // In 2D, curl maps H(curl) → L2 (scalar), so add an L2 space for B = curl A.
  if (mesh.back()->Dimension() == 2)
  {
    l2_curl_fec = std::make_unique<mfem::L2_FECollection>(
        solver.order - 1, mesh.back()->Dimension(), mfem::BasisType::GaussLegendre,
        mfem::FiniteElement::INTEGRAL);
    l2_curl_fespace = std::make_unique<FiniteElementSpace>(*mesh.back(), l2_curl_fec.get());
  }

  // Finalize setup.
  CheckBoundaryProperties();

  // Print essential BC information.
  if (dbc_attr.Size())
  {
    Mpi::Print("\nConfiguring Dirichlet BC at attributes:\n");
    utils::PrettyPrint(dbc_attr);
  }
}

CurlCurlOperator::CurlCurlOperator(const IoData &iodata,
                                   const std::vector<std::unique_ptr<Mesh>> &mesh)
  : CurlCurlOperator(iodata.boundaries, iodata.solver, iodata.domains.materials,
                     iodata.problem.type, iodata.units, mesh)
{
  surf_flux_op = SurfaceFluxOperator(iodata);

  // Identify London flux films (FluxLoopPEC ∩ Superconductor). For these the whole film is
  // a free unknown governed by the sheet term; the flux excitation uses the shifted London
  // penalty (RHS = M_sheet·a_h in GetFluxExcitationVector), not the whole-film PEC clamp.
  for (const auto &[idx, data] : surf_flux_op)
  {
    for (int attr : data.fluxloop_pec)
    {
      if (london_flux_film_attr_.count(attr))
      {
        london_flux_loops_.insert(idx);
        break;
      }
    }
  }

  // Build the fluxoid circulation functional c for each London flux loop (needs
  // surf_flux_op, set above, and the discrete curl matrix).
  SetUpLondonFluxConstraints();
}

mfem::Array<int> CurlCurlOperator::SetUpBoundaryProperties(
    const config::PecBoundaryData &pec, const std::map<int, config::FluxLoopData> &fluxloop,
    const std::vector<config::SuperconductorData> &superconductor,
    const mfem::ParMesh &mesh)
{
  // Check that boundary attributes have been specified correctly.
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> bdr_attr_marker(bdr_attr_max);
  bdr_attr_marker = 0;
  for (auto attr : mesh.bdr_attributes)
  {
    bdr_attr_marker[attr - 1] = 1;
  }
  if (!pec.empty())
  {
    std::set<int> bdr_warn_list;
    for (auto attr : pec.attributes)
    {
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
      continue;
    }
    dbc_bcs.Append(attr);
  }
  // Superconductor-sheet attributes (finite-λ films). A FluxLoopPEC film that is also a
  // superconductor sheet is a "London flux film": its whole surface is a free unknown
  // governed by the sheet term, so it is NOT marked essential here. Its hole fluxoid is an
  // integral constraint enforced by the range-space two-solve.
  std::set<int> sc_attrs;
  for (const auto &data : superconductor)
  {
    for (auto attr : data.attributes)
    {
      sc_attrs.insert(attr);
    }
  }

  // Add flux loop boundary attributes as essential boundaries, except London flux films.
  std::set<int> flux_attrs;
  for (const auto &[idx, data] : fluxloop)
  {
    for (auto attr : data.fluxloop_pec)
    {
      MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max && bdr_attr_marker[attr - 1],
                  "Unknown FluxLoopPEC boundary attribute "
                      << attr << " for FluxLoop index " << idx << "!");
      if (sc_attrs.count(attr))
      {
        // London flux film: keep the interior free (sheet term), do not clamp Dirichlet.
        london_flux_film_attr_.insert(attr);
      }
      else
      {
        flux_attrs.insert(attr);
      }
    }
    for (auto attr : data.hole_attributes)
    {
      MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max && bdr_attr_marker[attr - 1],
                  "Unknown FluxLoop HoleAttributes boundary attribute "
                      << attr << " for FluxLoop index " << idx << "!");
    }
  }
  for (auto attr : flux_attrs)
  {
    dbc_bcs.Append(attr);  // Each attribute added only once
  }
  return dbc_bcs;
}

void CurlCurlOperator::SetUpLondonFluxConstraints()
{
  if (london_flux_loops_.empty())
  {
    return;
  }
  // The fluxoid circulation functional requires B = curl A in an H(div) (RT) space so that
  // ∫_hole B·n is a well-defined face functional. In 2D curl maps to a scalar L2 space and
  // this construction does not apply.
  MFEM_VERIFY(!l2_curl_fespace,
              "London flux films (FluxLoopPEC ∩ Superconductor) require a 3D mesh!");

  const auto &pmesh = static_cast<const mfem::ParMesh &>(GetMesh());
  const Operator &curl = GetCurlMatrix();
  auto &curl_fespace = GetCurlSpace();
  for (int idx : london_flux_loops_)
  {
    const auto &data = surf_flux_op.GetSource(idx);

    // One generator, one α, one uniform hole functional c ⇒ a single scalar fluxoid per
    // loop. A loop listing several holes can only pin their sum, so model each hole as its
    // own loop.
    MFEM_VERIFY(
        data.hole_attributes.size() == 1,
        "FluxLoop index "
            << idx << " lists " << data.hole_attributes.size()
            << " holes, but the cut-cohomology two-solve supports exactly one hole per "
               "flux loop. Model each hole as a separate FluxLoop excitation.");

    // Fluxoid functional c = Curlᵀ·f, where f is the hole-cap flux functional on the RT
    // (curl) space: f·B = ∫_hole B·n (assembled via VectorFEBoundaryFluxLFIntegrator over
    // the hole attributes). Then cᵀA = fᵀ(Curl A) = ∫_hole (curl A)·n = ∮_∂hole A·dl
    // (Stokes), exact for any ND order. We zero c on the PEC essential set so the
    // fluxoid-mode solve A_h = K⁻¹c has A_h = 0 there under DIAG_ONE elimination.
    mfem::Array<int> hole_attr(const_cast<int *>(data.hole_attributes.data()),
                               static_cast<int>(data.hole_attributes.size()));
    auto hole_marker = mesh::BdrAttrToMarker(pmesh, hole_attr, true);

    FluxOrientationCoefficient sgn(data.direction, pmesh.SpaceDimension());
    mfem::LinearForm lf(&curl_fespace.Get());
    lf.AddBoundaryIntegrator(new mfem::VectorFEBoundaryFluxLFIntegrator(sgn), hole_marker);
    lf.UseFastAssembly(false);
    lf.UseDevice(false);
    lf.Assemble();
    lf.UseDevice(true);

    Vector f(curl_fespace.GetTrueVSize());
    f.UseDevice(true);
    f = 0.0;
    curl_fespace.GetProlongationMatrix()->AddMultTranspose(lf, f);

    Vector c(GetNDSpace().GetTrueVSize());
    c.UseDevice(true);
    c = 0.0;
    curl.MultTranspose(f, c);
    linalg::SetSubVector(c, dbc_tdof_lists.back(), 0.0);
    london_flux_constraint_[idx] = std::move(c);
  }
}

void CurlCurlOperator::CheckBoundaryProperties()
{
  // A final check that no boundary attribute is assigned multiple boundary conditions.
  const mfem::ParMesh &mesh = GetMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  const auto dbc_marker = mesh::AttrToMarker(bdr_attr_max, dbc_attr);
  const auto surf_j_marker = mesh::AttrToMarker(bdr_attr_max, surf_j_op.GetAttrList());
  for (int i = 0; i < dbc_marker.Size(); i++)
  {
    MFEM_VERIFY(dbc_marker[i] + surf_j_marker[i] <= 1,
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
    const auto geom_types = mesh::CheckElements(mesh).GetGeomTypes(mesh.Dimension());
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

    Mpi::Print("\nAssembling multigrid hierarchy:\n");
  }
}

}  // namespace

std::unique_ptr<Operator> CurlCurlOperator::GetStiffnessMatrix()
{
  // When partially assembled, the coarse operators can reuse the fine operator quadrature
  // data if the spaces correspond to the same mesh.
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);

  constexpr bool skip_zeros = false;
  MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                         mat_op.GetCurlCurlInvPermeability());
  BilinearForm k(GetNDSpace());
  k.AddDomainIntegrator<CurlCurlIntegrator>(muinv_func);
  // Add the thin-film superconductor kinetic sheet term (1/L_ksq) A_t · v_t as a tangential
  // surface mass on the ND space.
  MaterialPropertyCoefficient fbr(mat_op.MaxCeedBdrAttribute());
  sc_sheet_op.AddStiffnessBdrCoefficients(1.0, fbr);
  if (!fbr.empty())
  {
    k.AddBoundaryIntegrator<VectorFEMassIntegrator>(fbr);
  }
  // k.AssembleQuadratureData();
  auto k_vec = k.Assemble(GetNDSpaces(), skip_zeros);
  auto K = std::make_unique<MultigridOperator>(GetNDSpaces().GetNumLevels());
  for (std::size_t l = 0; l < GetNDSpaces().GetNumLevels(); l++)
  {
    const auto &nd_fespace_l = GetNDSpaces().GetFESpaceAtLevel(l);
    if (print_hdr)
    {
      Mpi::Print(" Level {:d} (p = {:d}): {:d} unknowns", l,
                 nd_fespace_l.GetMaxElementOrder(), nd_fespace_l.GlobalTrueVSize());
      if (const auto *k_spm = dynamic_cast<const hypre::HypreCSRMatrix *>(k_vec[l].get()))
      {
        HYPRE_BigInt nnz = k_spm->NNZ();
        Mpi::GlobalSum(1, &nnz, nd_fespace_l.GetComm());
        Mpi::Print(", {:d} NNZ\n", nnz);
      }
      else
      {
        Mpi::Print("\n");
      }
    }
    auto K_l = std::make_unique<ParOperator>(std::move(k_vec[l]), nd_fespace_l);
    K_l->SetEssentialTrueDofs(dbc_tdof_lists[l], Operator::DiagonalPolicy::DIAG_ONE);
    K->AddOperator(std::move(K_l));
  }

  print_hdr = false;
  return K;
}

std::unique_ptr<Operator> CurlCurlOperator::GetPreconditionerMatrix()
{
  // Preconditioner-only gauge shift for the London magnetostatic path. Without a London
  // flux film (or with a non-positive shift) there is no residual gradient null to lift, so
  // we return nullptr and the caller preconditions with the stiffness matrix itself.
  if (london_flux_loops_.empty() || !(london_pc_shift_ > 0.0))
  {
    return nullptr;
  }

  // Assemble K̃_pc = A_curlcurl + (1/L_ksq) M_sheet + london_pc_shift_ · (1/µ) ∫|A|². The
  // first two terms match GetStiffnessMatrix (the preconditioner tracks the true operator);
  // the volume mass is a reluctance-weighted regularization that lifts the film-constant
  // gradient mode so GᵀK̃_pcG is SPD for AMS. It touches the preconditioner only.
  constexpr bool skip_zeros = false;
  MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                         mat_op.GetCurlCurlInvPermeability());
  BilinearForm k(GetNDSpace());
  k.AddDomainIntegrator<CurlCurlIntegrator>(muinv_func);
  MaterialPropertyCoefficient mass_func(mat_op.GetAttributeToMaterial(),
                                        mat_op.GetInvPermeability(), london_pc_shift_);
  k.AddDomainIntegrator<VectorFEMassIntegrator>(mass_func);
  MaterialPropertyCoefficient fbr(mat_op.MaxCeedBdrAttribute());
  sc_sheet_op.AddStiffnessBdrCoefficients(1.0, fbr);
  if (!fbr.empty())
  {
    k.AddBoundaryIntegrator<VectorFEMassIntegrator>(fbr);
  }
  auto k_vec = k.Assemble(GetNDSpaces(), skip_zeros);
  auto K = std::make_unique<MultigridOperator>(GetNDSpaces().GetNumLevels());
  for (std::size_t l = 0; l < GetNDSpaces().GetNumLevels(); l++)
  {
    const auto &nd_fespace_l = GetNDSpaces().GetFESpaceAtLevel(l);
    auto K_l = std::make_unique<ParOperator>(std::move(k_vec[l]), nd_fespace_l);
    K_l->SetEssentialTrueDofs(dbc_tdof_lists[l], Operator::DiagonalPolicy::DIAG_ONE);
    K->AddOperator(std::move(K_l));
  }
  return K;
}

const Operator &
CurlCurlOperator::GetScreenedStiffnessMatrix(const mfem::Array<int> &extra_dbc_attr)
{
  // Key the cache on the sorted, deduplicated shorted-attribute set so that different call
  // orderings of the same ports map to one entry. Steps sharing a key are solved
  // back-to-back, so the newest entry is the only reuse candidate.
  std::set<int> unique_attr(extra_dbc_attr.begin(), extra_dbc_attr.end());
  std::vector<int> key(unique_attr.begin(), unique_attr.end());
  if (!screened_stiffness_cache.empty() && screened_stiffness_cache.back().key == key)
  {
    return *screened_stiffness_cache.back().K;
  }

  // Assemble the stiffness matrix with the base plus extra essential attributes. Size the
  // marker from the global max attribute (BdrAttrToMarker) so a rank whose partition lacks
  // the shorted-port boundary does not build an undersized, corrupt marker.
  mfem::Array<int> merged_attr(dbc_attr);
  merged_attr.Append(extra_dbc_attr);

  const auto &pmesh = static_cast<const mfem::ParMesh &>(GetMesh());
  auto merged_marker = mesh::BdrAttrToMarker(pmesh, merged_attr, true);

  constexpr bool skip_zeros = false;
  MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                         mat_op.GetCurlCurlInvPermeability());
  BilinearForm k(GetNDSpace());
  k.AddDomainIntegrator<CurlCurlIntegrator>(muinv_func);
  // Include the thin-film superconductor kinetic sheet term so shorted-port assemblies
  // carry the same kinetic inductance as the base stiffness matrix.
  MaterialPropertyCoefficient fbr(mat_op.MaxCeedBdrAttribute());
  sc_sheet_op.AddStiffnessBdrCoefficients(1.0, fbr);
  if (!fbr.empty())
  {
    k.AddBoundaryIntegrator<VectorFEMassIntegrator>(fbr);
  }
  auto k_vec = k.Assemble(GetNDSpaces(), skip_zeros);

  // Add the new entry before binding essential DOFs, since SetEssentialTrueDofs keeps a
  // shallow reference to the tdof list and that list must stay put. Drop the older entry
  // once two are cached; the one left behind is untouched, so its operator stays valid.
  if (screened_stiffness_cache.size() == 2)
  {
    screened_stiffness_cache.pop_front();
  }
  auto &entry = screened_stiffness_cache.emplace_back();
  entry.key = std::move(key);
  entry.dbc_tdof_lists.assign(GetNDSpaces().GetNumLevels(), mfem::Array<int>());
  auto K = std::make_unique<MultigridOperator>(GetNDSpaces().GetNumLevels());
  for (std::size_t l = 0; l < GetNDSpaces().GetNumLevels(); l++)
  {
    const auto &nd_fespace_l = GetNDSpaces().GetFESpaceAtLevel(l);
    nd_fespace_l.Get().GetEssentialTrueDofs(merged_marker, entry.dbc_tdof_lists[l]);
    // London flux films add no essential DOFs (the whole film is free, governed by the
    // sheet term; the hole fluxoid is an integral constraint), so the screened essential
    // set here is exactly the PEC + shorted-port set — no perimeter augmentation.
    auto K_l = std::make_unique<ParOperator>(std::move(k_vec[l]), nd_fespace_l);
    K_l->SetEssentialTrueDofs(entry.dbc_tdof_lists[l], Operator::DiagonalPolicy::DIAG_ONE);
    K->AddOperator(std::move(K_l));
  }
  entry.K = std::move(K);
  return *entry.K;
}

void CurlCurlOperator::ZeroEssentialTrueDofs(const mfem::Array<int> &extra_dbc_attr,
                                             Vector &v) const
{
  // Zero the excitation on the merged essential set (base Dirichlet plus shorted-port
  // attributes) so DIAG_ONE elimination sees a zero RHS on every constrained true DOF,
  // including edges an active port shares with a shorted one.
  mfem::Array<int> merged_attr(dbc_attr);
  merged_attr.Append(extra_dbc_attr);

  const auto &pmesh = static_cast<const mfem::ParMesh &>(GetMesh());
  auto merged_marker = mesh::BdrAttrToMarker(pmesh, merged_attr, true);

  mfem::Array<int> tdof_list;
  GetNDSpace().Get().GetEssentialTrueDofs(merged_marker, tdof_list);
  linalg::SetSubVector(v, tdof_list, 0.0);
}

void CurlCurlOperator::GetCurrentExcitationVector(int idx, Vector &RHS)
{
  // Assemble the surface current excitation +J. The SurfaceCurrentOperator assembles -J
  // (meant for time or frequency domain Maxwell discretization, so we multiply by -1 to
  // retrieve +J).
  SumVectorCoefficient fb(GetMesh().SpaceDimension());
  surf_j_op.AddExcitationBdrCoefficients(idx, fb);
  RHS.SetSize(GetNDSpace().GetTrueVSize());
  RHS.UseDevice(true);
  RHS = 0.0;
  int empty = (fb.empty());
  Mpi::GlobalMin(1, &empty, GetComm());
  if (empty)
  {
    return;
  }
  mfem::LinearForm rhs(&GetNDSpace().Get());
  rhs.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fb));
  rhs.UseFastAssembly(false);
  rhs.UseDevice(false);
  rhs.Assemble();
  rhs.UseDevice(true);
  GetNDSpace().GetProlongationMatrix()->AddMultTranspose(rhs, RHS, -1.0);
  linalg::SetSubVector(RHS, dbc_tdof_lists.back(), 0.0);
}

template <ProblemType T>
void CurlCurlOperator::GetFluxExcitationVector(int idx, Vector &RHS,
                                               PostOperator<T> &post_op)
{
  GetFluxExcitationVector(idx, RHS, post_op, nullptr);
}

template <ProblemType T>
void CurlCurlOperator::GetFluxExcitationVector(int idx, Vector &RHS,
                                               PostOperator<T> &post_op,
                                               Vector *boundary_values)
{
  RHS.SetSize(GetNDSpace().GetTrueVSize());
  RHS.UseDevice(true);
  RHS = 0.0;

  // Solve the surface problem for a_h, the film field carrying the hole fluxoid ∮a_h·dl =
  // Φ. Pure-PEC uses the uniform-perimeter Dirichlet clamp; London uses the curl-free cut
  // cohomology generator, whose gradient part the 3D solve absorbs (no spurious penalty
  // energy).
  const bool london = IsLondonFluxLoop(idx);
  Vector flux_solution;
  if (boundary_values)
  {
    SolveSurfaceCurlProblem(idx, post_op, *boundary_values, london);
    flux_solution = *boundary_values;  // Use the pre-allocated result
  }
  else
  {
    flux_solution = SolveSurfaceCurlProblem(idx, post_op, london);
  }

  // Early exit if boundary values are zero
  double boundary_norm = linalg::Norml2(GetComm(), flux_solution);
  if (boundary_norm < 1e-12)
  {
    return;
  }

  if (london)
  {
    // London flux film: the shifted penalty ½∫_Σ (1/L_ksq)|A_t − a_h|² has Euler-Lagrange
    // equation [K_cc + M_sheet] A = M_sheet·a_h, so with the sheet-inclusive base operator
    // K on the LHS the RHS is M_sheet·a_h. The film is absent from dbc_tdof_lists (only PEC
    // DOFs are pinned), so its interior relaxes freely; as λ→0 the penalty drives A_t → a_h
    // (geometric limit) and at finite λ the film penetrates and stores kinetic energy.
    if (!M_sheet_)
    {
      MaterialPropertyCoefficient fbr(mat_op.MaxCeedBdrAttribute());
      sc_sheet_op.AddStiffnessBdrCoefficients(1.0, fbr);
      BilinearForm m(GetNDSpace());
      if (!fbr.empty())
      {
        m.AddBoundaryIntegrator<VectorFEMassIntegrator>(fbr);
      }
      auto m_mat = m.Assemble(GetNDSpaces(), false);
      M_sheet_ = std::make_unique<ParOperator>(std::move(m_mat.back()), GetNDSpace());
    }
    // Normalize the generator so its exact fluxoid cᵀa_h = ∮a_h·dl equals Φ. The driver
    // still enforces the constraint on the full solution via α; this just keeps the drive
    // a_h normalized.
    {
      double c_ah = MeasureLondonHoleFlux(idx, flux_solution);
      double phi = surf_flux_op.GetSource(idx).GetExcitationFlux();
      if (std::abs(c_ah) > 1.0e-30)
      {
        double s = phi / c_ah;
        flux_solution *= s;
        if (boundary_values)
        {
          *boundary_values = flux_solution;
        }
      }
    }
    M_sheet_->Mult(flux_solution, RHS);
    linalg::SetSubVector(RHS, dbc_tdof_lists.back(), 0.0);
    return;
  }

  // Pure-PEC flux loop (unchanged): clamp the whole film Dirichlet.
  // Lift with the curl-curl-only operator and prescribe the whole-film boundary
  // distribution.
  if (!K_orig_)
  {
    MaterialPropertyCoefficient muinv_func(mat_op.GetAttributeToMaterial(),
                                           mat_op.GetInvPermeability());
    BilinearForm k(GetNDSpace());
    k.AddDomainIntegrator<CurlCurlIntegrator>(muinv_func);
    auto k_mat = k.Assemble(GetNDSpaces(), false);
    K_orig_ = std::make_unique<ParOperator>(std::move(k_mat.back()), GetNDSpace());
  }

  // Compute RHS = -K × boundary_values for boundary-interior coupling
  K_orig_->Mult(flux_solution, RHS);
  RHS *= -1.0;

  // Set boundary DOF entries to the prescribed values
  linalg::SetSubVector(RHS, dbc_tdof_lists.back(), flux_solution);
}

template <ProblemType T>
Vector CurlCurlOperator::SolveSurfaceCurlProblem(int flux_loop_idx,
                                                 PostOperator<T> &post_op,
                                                 bool harmonic_generator) const
{
  // Validate flux loop index exists
  MFEM_VERIFY(surf_flux_op.Size() > 0, "No flux loops configured!");
  surf_flux_op.GetSource(flux_loop_idx);  // Will throw if not found

  Vector result;
  surf_flux_op.SolveSurfaceCurlProblem(flux_loop_idx, GetMesh(), GetNDSpace(), post_op,
                                       result, harmonic_generator);
  return result;
}

template <ProblemType T>
void CurlCurlOperator::SolveSurfaceCurlProblem(int flux_loop_idx, PostOperator<T> &post_op,
                                               Vector &result,
                                               bool harmonic_generator) const
{
  // Validate flux loop index exists
  MFEM_VERIFY(surf_flux_op.Size() > 0, "No flux loops configured!");
  surf_flux_op.GetSource(flux_loop_idx);  // Will throw if not found

  surf_flux_op.SolveSurfaceCurlProblem(flux_loop_idx, GetMesh(), GetNDSpace(), post_op,
                                       result, harmonic_generator);
}

// Explicit template instantiations for PostOperator<ProblemType::MAGNETOSTATIC>
template void CurlCurlOperator::GetFluxExcitationVector<ProblemType::MAGNETOSTATIC>(
    int idx, Vector &RHS, PostOperator<ProblemType::MAGNETOSTATIC> &post_op);

template void CurlCurlOperator::GetFluxExcitationVector<ProblemType::MAGNETOSTATIC>(
    int idx, Vector &RHS, PostOperator<ProblemType::MAGNETOSTATIC> &post_op,
    Vector *boundary_values);

template Vector CurlCurlOperator::SolveSurfaceCurlProblem<ProblemType::MAGNETOSTATIC>(
    int flux_loop_idx, PostOperator<ProblemType::MAGNETOSTATIC> &post_op,
    bool harmonic_generator) const;

template void CurlCurlOperator::SolveSurfaceCurlProblem<ProblemType::MAGNETOSTATIC>(
    int flux_loop_idx, PostOperator<ProblemType::MAGNETOSTATIC> &post_op, Vector &result,
    bool harmonic_generator) const;

}  // namespace palace
