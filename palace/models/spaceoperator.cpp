// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "spaceoperator.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <numbers>
#include <optional>
#include <set>
#include <type_traits>
#include <utility>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "fem/libceed/ceed.hpp"
#include "fem/mesh.hpp"
#include "fem/multigrid.hpp"
#include "fem/qfunctions/coeff/coeff_qf.h"
#include "linalg/hypre.hpp"
#include "linalg/iterative.hpp"
#include "linalg/jacobi.hpp"
#include "linalg/ksp.hpp"
#include "linalg/rap.hpp"
#include "models/environmentoperator.hpp"
#include "models/floquetportoperator.hpp"
#include "models/pml.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"
#include "utils/units.hpp"

namespace palace
{

using namespace std::complex_literals;

SpaceOperator::SpaceOperator(const config::SolverData &solver,
                             const config::DomainData &domains,
                             const config::BoundaryData &boundaries,
                             ProblemType problem_type, const Units &units,
                             const std::vector<std::unique_ptr<Mesh>> &mesh)
  : pc_mat_real(solver.linear.pc_mat_real), pc_mat_shifted(solver.linear.pc_mat_shifted),
    pc_coarse_real(UseRealCoarseApproximation(solver.linear)),
    pml_subdomain(solver.linear.pml_subdomain_solver == PMLSubdomainSolver::DIRECT),
    pml_smoother_threshold(solver.linear.pml_smoother_threshold), print_hdr(true),
    print_prec_hdr(true), dbc_attr(SetUpBoundaryProperties(boundaries.pec, *mesh.back())),
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
    mat_op(domains.materials, boundaries.periodic, problem_type, *mesh.back(), domains.pml,
           solver.order),
    current_dipole_op(domains.current_dipole, units, *mesh.back()),
    surface_terms(boundaries, problem_type, units, mat_op, *mesh.back()),
    lumped_port_op(boundaries.lumpedport, units, mat_op, *mesh.back()),
    wave_port_op(boundaries, domains, solver, problem_type, units, mat_op, GetNDSpace(),
                 GetH1Space()),
    floquet_port_op(boundaries.floquetport, boundaries.periodic, problem_type, units,
                    mat_op, GetNDSpace().Get()),
    surf_j_op(boundaries.current, *mesh.back()),
    port_excitation_helper(lumped_port_op, wave_port_op, floquet_port_op, surf_j_op,
                           current_dipole_op)
{
  // The PML is only available in 3D. MaterialOperator treats PML regions as regular
  // materials on 2D meshes (as is required for the wave port cross-section submeshes), and
  // for other simulation types (with a warning).
  MFEM_VERIFY(
      (problem_type != ProblemType::DRIVEN && problem_type != ProblemType::EIGENMODE) ||
          (mesh.back()->Dimension() == 3 && mesh.back()->SpaceDimension() == 3) ||
          domains.pml.empty(),
      "PML regions are only supported for 3D simulations!");

  // Configure the two shared libCEED element sets before any geometry data is built.
  SetUpCeedElementSets(mesh);
  SetUpFrequencyDependentPermittivityMassOperators();

  // The second-order absorbing boundary condition is not passive for large tangential
  // Floquet wave vectors.
  surface_terms.GetFarfieldOp().CheckFloquetPassivity(solver, problem_type, units,
                                                      *mesh.back());

  // In 2D, curl maps H(curl) → L2 (scalar), so we need an L2 FE space for B = curl E.
  // Must use INTEGRAL map type so the discrete interpolator recognizes this as the curl
  // target space.
  if (mesh.back()->Dimension() == 2)
  {
    const int l2_order = solver.order - 1;
    const int dim = mesh.back()->Dimension();
    l2_curl_fecs.push_back(std::make_unique<mfem::L2_FECollection>(
        l2_order, dim, mfem::BasisType::GaussLegendre, mfem::FiniteElement::INTEGRAL));
    l2_curl_fespaces = std::make_unique<FiniteElementSpaceHierarchy>(
        fem::ConstructFiniteElementSpaceHierarchy<mfem::L2_FECollection>(1, mesh,
                                                                         l2_curl_fecs));
  }

  // Finalize setup.
  CheckBoundaryProperties();

  // Print essential BC information.
  if (dbc_attr.Size())
  {
    Mpi::Print("\nConfiguring Dirichlet PEC BC at attributes:\n");
    utils::PrettyPrint(dbc_attr);
  }

  // The PML integrators are only assembled on the elements of the PML regions, and evaluate
  // the material tensors at the physical coordinates of their quadrature points.
  if (mat_op.HasPML())
  {
    PrintPML(units);

    // Boundary conditions on the boundaries of the PML regions are not transformed by the
    // PML stretch. A port there would be computed without it (a wave port which crosses a
    // PML region solves its mode on the unstretched cross-section), so ports must not touch
    // the PML regions, as for config["Boundaries"]["PML"]. The other boundary conditions
    // there use the background material properties of the PML regions, with a warning,
    // except an absorbing boundary, which can terminate the layer: on its back, or across
    // its ends on an open side of an extruded box, where the stretch is one at the
    // interface with the physical region and the layer absorbs the field beyond it.
    const auto pml_attributes = mat_op.GetPMLAttributes();
    const std::set<int> pml_attr(pml_attributes.begin(), pml_attributes.end());
    const mfem::ParMesh &pmesh = GetMesh();
    int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
    Mpi::GlobalMax(1, &bdr_attr_max, GetComm());
    auto PMLBoundaryAttributes =
        [&](const std::vector<mfem::Array<int>> &attr_lists, const std::set<int> &skip = {})
    {
      mfem::Array<int> bdr_marker(bdr_attr_max);
      bdr_marker = 0;
      for (const auto &attr_list : attr_lists)
      {
        for (auto attr : attr_list)
        {
          const auto it = mesh.back()->GetCeedBdrAttributes().find(attr);
          if (attr > bdr_attr_max || skip.contains(attr) ||
              it == mesh.back()->GetCeedBdrAttributes().end())
          {
            continue;
          }
          for (const auto &[nbr_attr, ceed_attr] : it->second)
          {
            if (pml_attr.contains(nbr_attr))
            {
              bdr_marker[attr - 1] = 1;
            }
          }
        }
      }
      Mpi::GlobalMax(bdr_attr_max, bdr_marker.HostReadWrite(), GetComm());
      std::vector<int> pml_bdr_attr;
      for (int i = 0; i < bdr_attr_max; i++)
      {
        if (bdr_marker[i])
        {
          pml_bdr_attr.push_back(i + 1);
        }
      }
      return pml_bdr_attr;
    };
    const auto port_attr =
        PMLBoundaryAttributes({lumped_port_op.GetAttrList(), wave_port_op.GetAttrList(),
                               floquet_port_op.GetAttrList()});
    MFEM_VERIFY(port_attr.empty(),
                "Ports on boundary attributes adjacent to PML regions are not supported, "
                "since they are computed without the PML stretch ("
                    << fmt::format("{}", fmt::join(port_attr, ", ")) << ")!");
    const auto &farfield_attr = surface_terms.GetFarfieldOp().GetAttrList();
    const auto pml_bdr_attr =
        PMLBoundaryAttributes({surface_terms.GetAttrList(), surf_j_op.GetAttrList()},
                              std::set<int>(farfield_attr.begin(), farfield_attr.end()));
    if (!pml_bdr_attr.empty())
    {
      Mpi::Warning(GetComm(),
                   "Boundary conditions on boundary attributes adjacent to PML regions "
                   "use the background material properties of the PML regions and are "
                   "not transformed by the PML stretch ({})!\n",
                   fmt::join(pml_bdr_attr, ", "));
    }
  }
}

void SpaceOperator::PrintPML(const Units &units) const
{
  constexpr std::array<const char *, 6> face_name = {"-x", "+x", "-y", "+y", "-z", "+z"};
  const double length_scale = units.GetMeshLengthRelativeScale();
  for (const auto &layer : mat_op.GetPMLLayers())
  {
    const auto &p = layer.GetStretch();
    Mpi::Print("\nConfiguring PML regions at attributes:\n");
    utils::PrettyPrint(layer.GetAttributes());
    if (p.frequency_dependent)
    {
      Mpi::Print(" Frequency-dependent stretch");
    }
    else
    {
      Mpi::Print(" Static stretch at f₀ = {:.3e} GHz",
                 units.Dimensionalize<Units::ValueType::FREQUENCY>(
                     p.reference_frequency / (2.0 * std::numbers::pi)));
    }
    Mpi::Print(", conductivity profile σ = σ_max r^{:d}\n", p.order);
    const auto sigma = [&units](double s)
    { return units.Dimensionalize<Units::ValueType::CONDUCTIVITY>(s); };
    for (int f = 0; f < 6; f++)
    {
      if (p.geometry.thickness[f] <= 0.0)
      {
        continue;
      }
      Mpi::Print("  {}: thickness = {:.3e}, element size = {:.3e}, ",
                 p.extruded ? "extruded" : face_name[f],
                 p.geometry.thickness[f] * length_scale, p.element_size[f] * length_scale);
      Mpi::Print("σ_max = {:.3e} S/m, theoretical reflection = {:.1e}\n",
                 sigma(p.sigma_max[f]), pml::TheoreticalReflection(p, f));
    }
  }
}

SpaceOperator::SpaceOperator(const IoData &iodata,
                             const std::vector<std::unique_ptr<Mesh>> &mesh)
  : SpaceOperator(iodata.solver, iodata.domains, iodata.boundaries, iodata.problem.type,
                  iodata.units, mesh)
{
  // Validate excitations after wave port setup is complete.
  CheckExcitations(iodata.problem.type);
}

void SpaceOperator::SetUpCeedElementSets(const std::vector<std::unique_ptr<Mesh>> &mesh)
{
  // The active domain set is the lossy volume set, the union of the exact nonzero supports
  // of the loss-type MaterialOperator coefficients, each extracted on its own (summing the
  // tensors first could cancel entries of opposite sign and drop an attribute). This
  // MaterialOperator has no imaginary permeability property; add it here if one is
  // introduced.
  MaterialPropertyCoefficient loss_epsilon(mat_op.GetAttributeToMaterial(),
                                           mat_op.GetPermittivityImag());
  MaterialPropertyCoefficient loss_conductivity(mat_op.GetAttributeToMaterial(),
                                                mat_op.GetConductivity());
  std::vector<int> active_attr;
  loss_epsilon.AddNonzeroAttributes(active_attr);
  loss_conductivity.AddNonzeroAttributes(active_attr);
  // A frequency-dependent material is lossy too, and its unit mass, energy and flux forms
  // are supported on its attributes only. Its exact zero poles are in the conductivity
  // above, its other terms are not in any loss coefficient.
  for (std::size_t m = 0; m < mat_op.NumFrequencyDependentPermittivityMaterials(); m++)
  {
    if (mat_op.HasFrequencyDependentPermittivityModel(m))
    {
      for (auto attr :
           mat_op.GetCeedAttributes(mat_op.GetFrequencyDependentPermittivityAttributes(m)))
      {
        active_attr.push_back(attr);
      }
    }
  }
  std::sort(active_attr.begin(), active_attr.end());
  active_attr.erase(std::unique(active_attr.begin(), active_attr.end()), active_attr.end());

  // Without loss, the active domain set of a problem with PML regions is the elements
  // outside of them: the bulk material terms vanish on the PML regions (see
  // MaterialOperator::GetBulkAttributeToMaterial), whose terms are separate integrators on
  // their own elements, so the matrix-free operators skip the PML elements for them.
  int has_lossy = !active_attr.empty();
  Mpi::GlobalMax(1, &has_lossy, mat_op.GetMesh().GetComm());
  if (!has_lossy && mat_op.HasPML())
  {
    for (const auto &[attr, local_attr] : mat_op.GetMesh().GetCeedAttributes())
    {
      if (!mat_op.IsPMLCeedAttribute(local_attr))
      {
        active_attr.push_back(local_attr);
      }
    }
    std::sort(active_attr.begin(), active_attr.end());
  }

  // All boundary bilinear-form coefficients are supported on this global attribute union.
  std::vector<int> active_bdr_global;
  auto AddAttributes = [&active_bdr_global](const mfem::Array<int> &attr)
  { active_bdr_global.insert(active_bdr_global.end(), attr.begin(), attr.end()); };
  AddAttributes(surface_terms.GetAttrList());
  AddAttributes(lumped_port_op.GetAttrList());
  AddAttributes(wave_port_op.GetAttrList());
  AddAttributes(floquet_port_op.GetAttrList());
  std::sort(active_bdr_global.begin(), active_bdr_global.end());
  active_bdr_global.erase(std::unique(active_bdr_global.begin(), active_bdr_global.end()),
                          active_bdr_global.end());

  // Convert the fine-mesh process-local active set back to global attributes. Then map both
  // global unions directly into the local numbering of every h-level mesh. The p-levels
  // share the finest mesh and ordering.
  std::vector<int> active_global;
  for (const auto &[attr, local_attr] : mat_op.GetMesh().GetCeedAttributes())
  {
    if (std::binary_search(active_attr.begin(), active_attr.end(), local_attr))
    {
      active_global.push_back(attr);
    }
  }
  for (const auto &mesh_l : mesh)
  {
    auto active_attr_l = mesh_l->GetCeedAttributes(active_global);
    auto active_bdr_attr_l = mesh_l->GetCeedBdrAttributes(active_bdr_global);
    mesh_l->SetCeedActiveAttributes(
        std::vector<int>(active_attr_l.begin(), active_attr_l.end()),
        std::vector<int>(active_bdr_attr_l.begin(), active_bdr_attr_l.end()));
  }
}

void SpaceOperator::SetUpFrequencyDependentPermittivityMassOperators()
{
  frequency_dependent_permittivity_mass.resize(
      mat_op.NumFrequencyDependentPermittivityMaterials());
  for (std::size_t i = 0; i < frequency_dependent_permittivity_mass.size(); i++)
  {
    if (!mat_op.HasFrequencyDependentPermittivityA2(i) ||
        !mat_op.HasFrequencyDependentPermittivitySupport(i))
    {
      continue;
    }
    MaterialPropertyCoefficient f(mat_op.MaxCeedAttribute());
    f.AddMaterialProperty(
        mat_op.GetCeedAttributes(mat_op.GetFrequencyDependentPermittivityAttributes(i)),
        1.0);

    // Every rank builds the same cache entry. Ranks without this support assemble an empty,
    // dimensionally valid local operator, while the outer ParOperator supplies the common
    // parallel restriction/prolongation. B_m is only ever applied (in the A2 sums, and
    // borrowed by the eigenvalue, ROM and domain postprocessing operators), so a partially
    // assembled B_m omits the elements outside the lossy element set, which contains the
    // material, and a fully assembled one (order below PartialAssemblyOrder) stores no
    // explicit zeros.
    constexpr bool skip_zeros = true;
    BilinearForm a(GetNDSpace());
    a.SkipZeroCoefficientElements();
    if (!f.empty())
    {
      a.AddDomainIntegrator<VectorFEMassIntegrator>(f);
    }
    auto b = a.Assemble(skip_zeros);
    frequency_dependent_permittivity_mass[i] =
        std::make_unique<ParOperator>(std::move(b), GetNDSpace());
  }
}

std::unique_ptr<Operator> SpaceOperator::BuildFrequencyDependentPermittivityA2Operator(
    std::unique_ptr<Operator> &&base, const std::vector<double> &coeff) const
{
  MFEM_VERIFY(coeff.size() == frequency_dependent_permittivity_mass.size(),
              "Invalid dispersive material coefficient array!");
  auto sum = std::make_unique<SumOperator>(GetNDSpace().GetVSize());
  if (base)
  {
    sum->AddOperator(std::move(base));
  }
  for (std::size_t i = 0; i < coeff.size(); i++)
  {
    if (coeff[i] != 0.0 && frequency_dependent_permittivity_mass[i])
    {
      // B_m is owned by SpaceOperator and intentionally borrowed by each A2 sum.
      sum->AddOperator(frequency_dependent_permittivity_mass[i]->LocalOperator(), coeff[i]);
    }
  }
  return sum;
}

std::unique_ptr<ComplexOperator> SpaceOperator::GetFrequencyDependentPermittivityMassMatrix(
    std::size_t material_idx, Operator::DiagonalPolicy diag_policy) const
{
  // The cache entry exists on every rank exactly when the material has a nonlinear term
  // with global support (SetUpFrequencyDependentPermittivityMassOperators), so the null
  // contract is rank-uniform.
  const auto &B = frequency_dependent_permittivity_mass.at(material_idx);
  if (!B)
  {
    return {};
  }
  auto B_op =
      std::make_unique<ComplexParOperator>(nullptr, &B->LocalOperator(), GetNDSpace());
  B_op->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
  return B_op;
}

void SpaceOperator::AddFrequencyDependentPermittivityA2Coefficient(
    std::size_t material_idx, double coeff, MaterialPropertyCoefficient &f) const
{
  if (coeff != 0.0 && mat_op.HasFrequencyDependentPermittivitySupport(material_idx))
  {
    f.AddMaterialProperty(
        mat_op.GetCeedAttributes(
            mat_op.GetFrequencyDependentPermittivityAttributes(material_idx)),
        coeff);
  }
}

void SpaceOperator::GetFrequencyDependentPermittivityA2Coefficients(
    std::complex<double> omega, std::vector<double> &real, std::vector<double> &imag,
    bool include_poles) const
{
  real.resize(mat_op.NumFrequencyDependentPermittivityMaterials());
  imag.resize(mat_op.NumFrequencyDependentPermittivityMaterials());
  const std::complex<double> s = 1i * omega;
  for (std::size_t i = 0; i < real.size(); i++)
  {
    // A material without mesh support contributes nothing and is not evaluated, since it
    // may be singular at this frequency. Without the pole terms, only the
    // Djordjevic-Sarkar terms remain.
    const std::complex<double> value =
        !mat_op.HasFrequencyDependentPermittivitySupport(i) ? 0.0
        : include_poles ? mat_op.EvaluateFrequencyDependentPermittivityA2(i, s)
                        : mat_op.EvaluateDjordjevicSarkarPermittivityA2(i, s);
    real[i] = value.real();
    imag[i] = value.imag();
  }
}

void SpaceOperator::AddFrequencyDependentPermittivityA2Coefficients(
    std::complex<double> omega, MaterialPropertyCoefficient &fr,
    MaterialPropertyCoefficient &fi) const
{
  const std::complex<double> s = 1i * omega;
  for (std::size_t i = 0; i < mat_op.NumFrequencyDependentPermittivityMaterials(); i++)
  {
    // The real part is mass-like (s²χ), so it follows the same shift as the ε∞ mass; a
    // strongly dispersive material (Re χ > ε∞) would otherwise make the shifted
    // preconditioner indefinite. The imaginary part is loss-like and is left as is.
    if (!mat_op.HasFrequencyDependentPermittivitySupport(i))
    {
      continue;
    }
    const std::complex<double> g = mat_op.EvaluateFrequencyDependentPermittivityA2(i, s);
    AddFrequencyDependentPermittivityA2Coefficient(
        i, pc_mat_shifted ? std::abs(g.real()) : g.real(), fr);
    AddFrequencyDependentPermittivityA2Coefficient(i, g.imag(), fi);
  }
}

void SpaceOperator::AddFrequencyDependentPermittivityA2Coefficients(
    std::complex<double> omega, MaterialPropertyCoefficient &f) const
{
  const std::complex<double> s = 1i * omega;
  for (std::size_t i = 0; i < mat_op.NumFrequencyDependentPermittivityMaterials(); i++)
  {
    // Real preconditioners accumulate the real and imaginary slots into the same form, as
    // for other A2 terms, with the mass-like real part shifted like the ε∞ mass.
    if (!mat_op.HasFrequencyDependentPermittivitySupport(i))
    {
      continue;
    }
    const std::complex<double> g = mat_op.EvaluateFrequencyDependentPermittivityA2(i, s);
    AddFrequencyDependentPermittivityA2Coefficient(
        i, (pc_mat_shifted ? std::abs(g.real()) : g.real()) + g.imag(), f);
  }
}

void SpaceOperator::CheckExcitations(ProblemType problem_type) const
{
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
  const auto &farfield_op = surface_terms.GetFarfieldOp();
  const auto &surf_sigma_op = surface_terms.GetSurfaceConductivityOp();
  const auto &surf_z_op = surface_terms.GetSurfaceImpedanceOp();
  const auto &surf_rz_op = surface_terms.GetRationalImpedanceOp();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  const auto dbc_marker = mesh::AttrToMarker(bdr_attr_max, dbc_attr);
  const auto farfield_marker = mesh::AttrToMarker(bdr_attr_max, farfield_op.GetAttrList());
  const auto surf_sigma_marker =
      mesh::AttrToMarker(bdr_attr_max, surf_sigma_op.GetAttrList());
  const auto surf_z_Rs_marker = mesh::AttrToMarker(bdr_attr_max, surf_z_op.GetRsAttrList());
  const auto surf_z_Ls_marker = mesh::AttrToMarker(bdr_attr_max, surf_z_op.GetLsAttrList());
  const auto surf_rz_marker = mesh::AttrToMarker(bdr_attr_max, surf_rz_op.GetAttrList());
  const auto lumped_port_Rs_marker =
      mesh::AttrToMarker(bdr_attr_max, lumped_port_op.GetRsAttrList());
  const auto lumped_port_Ls_marker =
      mesh::AttrToMarker(bdr_attr_max, lumped_port_op.GetLsAttrList());
  const auto wave_port_marker =
      mesh::AttrToMarker(bdr_attr_max, wave_port_op.GetAttrList());
  const auto floquet_port_marker =
      mesh::AttrToMarker(bdr_attr_max, floquet_port_op.GetAttrList());
  mfem::Array<int> aux_bdr_marker(dbc_marker.Size());
  for (int i = 0; i < dbc_marker.Size(); i++)
  {
    aux_bdr_marker[i] = (dbc_marker[i] || farfield_marker[i] || surf_sigma_marker[i] ||
                         surf_z_Rs_marker[i] || surf_rz_marker[i] || surf_z_Ls_marker[i] ||
                         lumped_port_Rs_marker[i] || lumped_port_Ls_marker[i] ||
                         wave_port_marker[i] || floquet_port_marker[i]);
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
                        surf_z_marker[i] + surf_rz_marker[i] + lumped_port_marker[i] +
                        wave_port_marker[i] + floquet_port_marker[i] + surf_j_marker[i] <=
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
  }
  print_hdr = false;
}

// A PML integrator to add to a BilinearForm along with its QFunction context (see
// fem/qfunctions/coeff/pml_qf.h). The libCEED operators copy the context when they are
// built, so it only needs to outlive the call to BilinearForm::Assemble.
enum class PMLIntegKind : char
{
  CURL_CURL,         // c μ̃⁻¹ (curl u, curl v)
  MASS,              // c ε̃ (u, v)
  CURL_CURL_MASS,    // c_μ μ̃⁻¹ (curl u, curl v) + c_ε ε̃ (u, v)
  FLOQUET_MASS,      // c [k ×]ᵀ μ̃⁻¹ [k ×] (u, v)
  FLOQUET_CROSS,     // c ([k ×]ᵀ μ̃⁻¹ curl u, v) - c (μ̃⁻¹ [k ×] u, curl v)
  DIFFUSION,         // c ε̃ (∇u, ∇v), auxiliary space
  FLOQUET_DIFFUSION  // c [k ×]ᵀ μ̃⁻¹ [k ×] (∇u, ∇v), auxiliary space
};
struct PMLIntegrator
{
  std::vector<CeedIntScalar> ctx;
  PMLIntegKind kind;
  std::vector<int> attributes;  // Domain attributes of the PML regions
};

// The PML integrators store their quadrature data whatever the assemble_q_data of the other
// integrators of the form, including in the system matrices and on devices (see
// PMLIntegratorBase).
void AddPMLIntegrators(BilinearForm &a, const std::vector<PMLIntegrator> *pml, bool aux)
{
  if (!pml)
  {
    return;
  }
  for (const auto &p : *pml)
  {
    const void *ctx = p.ctx.data();
    const std::size_t bytes = p.ctx.size() * sizeof(CeedIntScalar);
    MFEM_ASSERT(aux == (p.kind == PMLIntegKind::DIFFUSION ||
                        p.kind == PMLIntegKind::FLOQUET_DIFFUSION),
                "Invalid PML integrator for the " << (aux ? "auxiliary" : "primary")
                                                  << " space operator!");
    switch (p.kind)
    {
      case PMLIntegKind::CURL_CURL:
        a.AddDomainIntegratorOnAttributes<CurlCurlPMLIntegrator>(p.attributes, ctx, bytes);
        break;
      case PMLIntegKind::MASS:
        a.AddDomainIntegratorOnAttributes<VectorFEMassPMLIntegrator>(p.attributes, ctx,
                                                                     bytes);
        break;
      case PMLIntegKind::CURL_CURL_MASS:
        a.AddDomainIntegratorOnAttributes<CurlCurlMassPMLIntegrator>(p.attributes, ctx,
                                                                     bytes);
        break;
      case PMLIntegKind::FLOQUET_MASS:
        a.AddDomainIntegratorOnAttributes<FloquetMassPMLIntegrator>(p.attributes, ctx,
                                                                    bytes);
        break;
      case PMLIntegKind::FLOQUET_CROSS:
        a.AddDomainIntegratorOnAttributes<FloquetCrossPMLIntegrator>(p.attributes, ctx,
                                                                     bytes);
        break;
      case PMLIntegKind::DIFFUSION:
        a.AddDomainIntegratorOnAttributes<DiffusionPMLIntegrator>(p.attributes, ctx, bytes);
        break;
      case PMLIntegKind::FLOQUET_DIFFUSION:
        a.AddDomainIntegratorOnAttributes<FloquetDiffusionPMLIntegrator>(p.attributes, ctx,
                                                                         bytes);
        break;
    }
  }
}

template <typename... CoeffType>
bool AreExactlyZero(const CoeffType &...coeff)
{
  return (... && coeff.IsExactlyZero());
}

// Add every configured (non-empty) term, including terms whose current value is exactly
// zero. This storage-based selection is used only for coarse preconditioner assembly, which
// retains the configured sparse structure for symbolic reuse; normal operator assembly uses
// AddIntegrators, which drops exactly-zero coefficients first.
void AddConfiguredIntegrators(
    BilinearForm &a, const MaterialPropertyCoefficient *df,
    const MaterialPropertyCoefficient *f, const MaterialPropertyCoefficient *dfb,
    const MaterialPropertyCoefficient *fb, const MaterialPropertyCoefficient *fp,
    const MaterialPropertyCoefficient *frb, const MaterialPropertyCoefficient *fpb,
    bool assemble_q_data = false, const std::vector<PMLIntegrator> *pml = nullptr)
{
  const bool has_df = df && !df->empty();
  const bool has_f = f && !f->empty();
  if (has_df && has_f)
  {
    a.AddDomainIntegrator<CurlCurlMassIntegrator>(*df, *f);
  }
  else
  {
    if (has_df)
    {
      a.AddDomainIntegrator<CurlCurlIntegrator>(*df);
    }
    if (has_f)
    {
      a.AddDomainIntegrator<VectorFEMassIntegrator>(*f);
    }
  }

  const bool has_dfb = dfb && !dfb->empty();
  const bool has_fb = fb && !fb->empty();
  if (has_dfb && has_fb)
  {
    a.AddBoundaryIntegrator<CurlCurlMassIntegrator>(*dfb, *fb);
  }
  else
  {
    if (has_dfb)
    {
      a.AddBoundaryIntegrator<CurlCurlIntegrator>(*dfb);
    }
    if (has_fb)
    {
      a.AddBoundaryIntegrator<VectorFEMassIntegrator>(*fb);
    }
  }
  if (fp && !fp->empty())
  {
    a.AddDomainIntegrator<MixedVectorWeakCurlIntegrator>(*fp);
    a.AddDomainIntegrator<MixedVectorCurlIntegrator>(*fp, true);
  }

  // Floquet terms of the second-order absorbing boundary condition.
  if (frb && !frb->empty())
  {
    a.AddBoundaryIntegrator<RotatedVectorFEMassIntegrator>(*frb);
  }
  if (fpb && !fpb->empty())
  {
    a.AddBoundaryIntegrator<MixedVectorWeakCurlIntegrator>(*fpb);
    a.AddBoundaryIntegrator<MixedVectorCurlIntegrator>(*fpb, true);
  }
  AddPMLIntegrators(a, pml, false);
  if (assemble_q_data)
  {
    a.AssembleQuadratureData();
  }
}

// Coefficients which are exactly zero on every attribute are dropped before the integrator
// type is chosen, so that a combined integrator degrades to its nonzero block (a mass
// integrator instead of a curl-curl + mass integrator with a zero curl-curl coefficient).
// Element subsets then handle coefficients which vanish on only part of the mesh.
const MaterialPropertyCoefficient *
ActiveCoefficient(const MaterialPropertyCoefficient *coeff)
{
  return (coeff && !coeff->IsExactlyZero()) ? coeff : nullptr;
}

void AddIntegrators(BilinearForm &a, const MaterialPropertyCoefficient *df,
                    const MaterialPropertyCoefficient *f,
                    const MaterialPropertyCoefficient *dfb,
                    const MaterialPropertyCoefficient *fb,
                    const MaterialPropertyCoefficient *fp,
                    const MaterialPropertyCoefficient *frb,
                    const MaterialPropertyCoefficient *fpb, bool assemble_q_data = false,
                    const std::vector<PMLIntegrator> *pml = nullptr)
{
  AddConfiguredIntegrators(a, ActiveCoefficient(df), ActiveCoefficient(f),
                           ActiveCoefficient(dfb), ActiveCoefficient(fb),
                           ActiveCoefficient(fp), ActiveCoefficient(frb),
                           ActiveCoefficient(fpb), assemble_q_data, pml);
}

void AddConfiguredAuxIntegrators(BilinearForm &a, const MaterialPropertyCoefficient *f,
                                 const MaterialPropertyCoefficient *fb,
                                 const MaterialPropertyCoefficient *frb,
                                 bool assemble_q_data = false,
                                 const std::vector<PMLIntegrator> *pml = nullptr)
{
  if (f && !f->empty())
  {
    a.AddDomainIntegrator<DiffusionIntegrator>(*f);
  }
  if (fb && !fb->empty())
  {
    a.AddBoundaryIntegrator<DiffusionIntegrator>(*fb);
  }
  if (frb && !frb->empty())
  {
    a.AddBoundaryIntegrator<RotatedVectorFEMassIntegrator>(*frb);
  }
  AddPMLIntegrators(a, pml, true);
  if (assemble_q_data)
  {
    a.AssembleQuadratureData();
  }
}

void AddAuxIntegrators(BilinearForm &a, const MaterialPropertyCoefficient *f,
                       const MaterialPropertyCoefficient *fb,
                       const MaterialPropertyCoefficient *frb, bool assemble_q_data = false,
                       const std::vector<PMLIntegrator> *pml = nullptr)
{
  AddConfiguredAuxIntegrators(a, ActiveCoefficient(f), ActiveCoefficient(fb),
                              ActiveCoefficient(frb), assemble_q_data, pml);
}

auto AssembleOperator(
    const AssemblyDomains &domains, const FiniteElementSpace &fespace,
    const MaterialPropertyCoefficient *df, const MaterialPropertyCoefficient *f,
    const MaterialPropertyCoefficient *dfb, const MaterialPropertyCoefficient *fb,
    const MaterialPropertyCoefficient *fp, const MaterialPropertyCoefficient *frb,
    const MaterialPropertyCoefficient *fpb, bool skip_zeros = false,
    bool assemble_q_data = false, const std::vector<PMLIntegrator> *pml = nullptr)
{
  BilinearForm a(fespace);
  // System matrices are only ever applied matrix-free (their sums are assembled into
  // SumOperator, which has no sparse representation), so elements where a coefficient is
  // exactly zero can be left out of the partially assembled operator.
  a.SkipZeroCoefficientElements();
  AddIntegrators(a, df, f, dfb, fb, fp, frb, fpb, assemble_q_data, pml);
  domains.Restrict(a);
  return a.Assemble(skip_zeros);
}

// The environment's term R_Γᵀ S_E(ω) R_Γ of the operator of a substructured region,
// complex symmetric (see EnvironmentOperator).
class EnvironmentTermOperator : public ComplexOperator
{
private:
  EnvironmentOperator &env;
  const double omega;

public:
  EnvironmentTermOperator(EnvironmentOperator &env, double omega, int n)
    : ComplexOperator(n), env(env), omega(omega)
  {
  }

  void Mult(const ComplexVector &x, ComplexVector &y) const override
  {
    y = 0.0;
    env.AddMult(omega, x, y);
  }

  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a) const override
  {
    if (a == 1.0)
    {
      env.AddMult(omega, x, y);
      return;
    }
    ComplexVector t(y.Size());
    t.UseDevice(true);
    t = 0.0;
    env.AddMult(omega, x, t);
    y.AXPY(a, t);
  }
};

auto AssembleOperators(
    const AssemblyDomains &domains, const FiniteElementSpaceHierarchy &fespaces,
    const MaterialPropertyCoefficient *df, const MaterialPropertyCoefficient *f,
    const MaterialPropertyCoefficient *dfb, const MaterialPropertyCoefficient *fb,
    const MaterialPropertyCoefficient *fp, const MaterialPropertyCoefficient *frb,
    const MaterialPropertyCoefficient *fpb, bool skip_zeros = false,
    bool assemble_q_data = false, bool flip_coarse_dfb = false,
    const MaterialPropertyCoefficient *dfb_flip_add = nullptr,
    const std::vector<PMLIntegrator> *pml = nullptr,
    const std::vector<PMLIntegrator> *coarse_pml = nullptr)
{
  // Optionally flip the sign of dfb on the coarse level only (see AssemblePreconditioner),
  // and add dfb_flip_add to the flipped copy. The PML integrators of the coarse level are
  // pml, or coarse_pml if given.
  std::optional<MaterialPropertyCoefficient> dfb_flip;
  const MaterialPropertyCoefficient *dfb_coarse = dfb;
  const bool flip_add = dfb_flip_add && !dfb_flip_add->empty();
  if (flip_coarse_dfb && dfb && (!dfb->empty() || flip_add))
  {
    dfb_flip.emplace(*dfb);
    *dfb_flip *= -1.0;
    if (flip_add)
    {
      dfb_flip->AddCoefficient(dfb_flip_add->GetAttributeToMaterial(),
                               dfb_flip_add->GetMaterialProperties());
    }
    dfb_coarse = &*dfb_flip;
  }

  // Keep configured coefficient support on the coarse sparse level. This is required for
  // an exact-complex coarse solve: a nonlinear eigenvalue seed λ = iω has Im(λ²) = 0, but
  // the corresponding full-domain imaginary mass term becomes nonzero once Re(λ) != 0.
  // Its off-diagonal block graph must exist before the first symbolic factorization. Use
  // the same rule for real-approximate coarse solves to keep this hierarchy assembly
  // independent of the chosen sparse-solver representation. Fine partial-assembly levels
  // have no symbolic structure to preserve, so exact-zero terms are omitted there.
  BilinearForm coarse(fespaces.GetFESpaceAtLevel(0));
  AddConfiguredIntegrators(coarse, df, f, dfb_coarse, fb, fp, frb, fpb, false,
                           coarse_pml ? coarse_pml : pml);
  domains.Restrict(coarse);
  std::vector<std::unique_ptr<Operator>> ops;
  ops.reserve(fespaces.GetNumLevels());
  ops.push_back(coarse.Assemble(skip_zeros));
  if (fespaces.GetNumLevels() > 1)
  {
    BilinearForm fine(fespaces.GetFinestFESpace());
    // Levels above the coarsest are only applied matrix-free (their smoothers use the
    // operator action and its diagonal), so elements where a coefficient is exactly zero
    // can be left out.
    fine.SkipZeroCoefficientElements();
    AddIntegrators(fine, df, f, dfb, fb, fp, frb, fpb, assemble_q_data, pml);
    domains.Restrict(fine);
    auto fine_ops = fine.Assemble(fespaces, skip_zeros, 1);
    for (auto &op : fine_ops)
    {
      ops.push_back(std::move(op));
    }
  }
  return ops;
}

auto AssembleAuxOperators(const AssemblyDomains &domains,
                          const FiniteElementSpaceHierarchy &fespaces,
                          const MaterialPropertyCoefficient *f,
                          const MaterialPropertyCoefficient *fb,
                          const MaterialPropertyCoefficient *frb, bool skip_zeros = false,
                          bool assemble_q_data = false,
                          const std::vector<PMLIntegrator> *pml = nullptr)
{
  // Match the configured coarse-support policy of the primary hierarchy above.
  BilinearForm coarse(fespaces.GetFESpaceAtLevel(0));
  AddConfiguredAuxIntegrators(coarse, f, fb, frb, false, pml);
  domains.Restrict(coarse);
  std::vector<std::unique_ptr<Operator>> ops;
  ops.reserve(fespaces.GetNumLevels());
  ops.push_back(coarse.Assemble(skip_zeros));
  if (fespaces.GetNumLevels() > 1)
  {
    BilinearForm fine(fespaces.GetFinestFESpace());
    fine.SkipZeroCoefficientElements();
    AddAuxIntegrators(fine, f, fb, frb, assemble_q_data, pml);
    domains.Restrict(fine);
    auto fine_ops = fine.Assemble(fespaces, skip_zeros, 1);
    for (auto &op : fine_ops)
    {
      ops.push_back(std::move(op));
    }
  }
  return ops;
}

// Operators with PML terms: a static stretch, at its fixed reference frequency, contributes
// to K, C, and M, while a frequency-dependent stretch, at the solve frequency, contributes
// to A2(ω).
enum class PMLFilter : char
{
  STATIC,
  FREQUENCY_DEPENDENT
};

// Append the PML integrator(s) for the complex terms c_muinv μ̃⁻¹ and/or c_eps ε̃ of the
// given kind (the prefactor not used by the kind is ignored), for each PML block whose
// stretch is selected by filter, with a frequency-dependent stretch evaluated at omega. The
// real part is appended to re and the imaginary part to im. For a real-valued
// approximation (null im), or with stretch_magnitude, the tensors use the magnitudes of the
// stretch factors (see pml::TensorPart). Nothing is appended if the terms are identically
// zero.
void AppendPML(std::vector<PMLIntegrator> &re, std::vector<PMLIntegrator> *im,
               const MaterialOperator &mat_op, PMLIntegKind kind,
               std::complex<double> c_muinv, std::complex<double> c_eps, PMLFilter filter,
               std::complex<double> omega = 0.0, bool stretch_magnitude = false)
{
  const bool floquet =
      (kind == PMLIntegKind::FLOQUET_MASS || kind == PMLIntegKind::FLOQUET_CROSS ||
       kind == PMLIntegKind::FLOQUET_DIFFUSION);
  if (!mat_op.HasPML() || (floquet && !mat_op.HasWaveVector()))
  {
    return;
  }
  if (kind == PMLIntegKind::CURL_CURL_MASS && c_muinv == 0.0)
  {
    kind = PMLIntegKind::MASS;
  }
  else if (kind == PMLIntegKind::CURL_CURL_MASS && c_eps == 0.0)
  {
    kind = PMLIntegKind::CURL_CURL;
  }
  const bool eps_only = (kind == PMLIntegKind::MASS || kind == PMLIntegKind::DIFFUSION);
  if ((eps_only ? c_eps : c_muinv) == 0.0)
  {
    return;
  }

  pml::ContextHeader header;
  header.c_muinv = eps_only ? 0.0 : c_muinv;
  header.c_eps = (eps_only || kind == PMLIntegKind::CURL_CURL_MASS) ? c_eps : 0.0;
  header.omega = omega;
  if (floquet)
  {
    const auto &kx = mat_op.GetWaveVectorCross();
    MFEM_VERIFY(kx.Height() == 3 && kx.Width() == 3,
                "Floquet PML terms require a 3D Floquet wave vector!");
    std::copy_n(kx.Data(), 9, header.wave_vector_cross.begin());
  }
  for (const auto &layer : mat_op.GetPMLLayers())
  {
    if ((filter == PMLFilter::FREQUENCY_DEPENDENT) != layer.IsFrequencyDependent())
    {
      continue;
    }
    const bool magnitude = stretch_magnitude || !im;
    for (auto [dst, part] :
         {std::pair{&re, magnitude ? pml::TensorPart::STRETCH_MAGNITUDE_REAL
                                   : pml::TensorPart::REAL},
          std::pair{im, magnitude ? pml::TensorPart::STRETCH_MAGNITUDE_IMAG
                                  : pml::TensorPart::IMAG}})
    {
      if (!dst)
      {
        continue;
      }
      header.part = part;
      auto ctx = layer.PackContext(header);
      if (!ctx.empty())
      {
        dst->push_back({std::move(ctx), kind, layer.GetAttributes()});
      }
    }
  }
}

// Shifted variant of a mass term prefactor for preconditioning, see pc_mat_shifted.
std::complex<double> ShiftedMassCoefficient(std::complex<double> a, bool shifted)
{
  return shifted ? std::complex<double>(std::abs(a.real()), a.imag()) : a;
}

// Coefficients of the PML contributions K_pml, C_pml, and M_pml to the stiffness, damping,
// and mass matrices in the combination a0 K_pml + a1 C_pml + a2 M_pml. These mirror the
// bulk terms: K_pml is the μ̃⁻¹ curl-curl term and, for a fixed Floquet wave vector, the
// Floquet terms [k ×]ᵀ μ̃⁻¹ [k ×] + i ([k ×]ᵀ μ̃⁻¹ curl ·, ·) - i (μ̃⁻¹ [k ×] ·, curl ·),
// C_pml is the Floquet cross term for a frequency-scaled wave vector k_F = ω k₀, and M_pml
// is the ε̃ mass term and, for a frequency-scaled wave vector, -[k₀ ×]ᵀ μ̃⁻¹ [k₀ ×]. The
// Floquet part of the mass term can have its own coefficient to match the bulk terms of the
// preconditioner matrices (see pc_mat_shifted).
struct PMLPencil
{
  std::complex<double> a0 = 0.0, a1 = 0.0, a2 = 0.0, a2_floquet = 0.0;
};

// Append the PML integrators for the pencil terms if the stretch is selected by filter,
// with the stretch evaluated at omega (frequency-dependent stretch) or at its reference
// frequency (static stretch). The auxiliary space terms are the gradient subspace
// projections of the mass terms. For a real-valued approximation (null im), the tensors use
// the magnitudes of the stretch factors and the Floquet cross terms are dropped. With
// stretch_magnitude, the complex-valued terms use the magnitudes of the stretch factors
// too.
//
// For a frequency-dependent stretch, A2(ω) = K_pml(ω) + iω C_pml(ω) - ω² M_pml(ω). The
// stretch factors s(ω) = κ + σ / (α + iω) and therefore the PML tensors are analytic in ω,
// so a single complex-ω evaluation covers both the real frequencies of driven simulations
// and the complex eigenfrequencies of the eigenmode nonlinear solve (ω = -iλ), where it
// yields the exact analytic continuation of A2.
void AppendPMLPencil(std::vector<PMLIntegrator> &re, std::vector<PMLIntegrator> *im,
                     std::vector<PMLIntegrator> *aux_re, std::vector<PMLIntegrator> *aux_im,
                     const MaterialOperator &mat_op, PMLFilter filter,
                     std::complex<double> omega, const PMLPencil &c,
                     bool stretch_magnitude = false)
{
  AppendPML(re, im, mat_op, PMLIntegKind::CURL_CURL_MASS, c.a0, c.a2, filter, omega,
            stretch_magnitude);
  if (aux_re)
  {
    AppendPML(*aux_re, aux_im, mat_op, PMLIntegKind::DIFFUSION, 0.0, c.a2, filter, omega);
  }
  if (mat_op.HasWaveVector())
  {
    const bool scaled = mat_op.HasFloquetFrequencyScaling();
    const auto c_mass = scaled ? -c.a2_floquet : c.a0, c_cross = scaled ? c.a1 : 1i * c.a0;
    AppendPML(re, im, mat_op, PMLIntegKind::FLOQUET_MASS, c_mass, 0.0, filter, omega,
              stretch_magnitude);
    if (im)
    {
      AppendPML(re, im, mat_op, PMLIntegKind::FLOQUET_CROSS, c_cross, 0.0, filter, omega,
                stretch_magnitude);
    }
    if (aux_re)
    {
      AppendPML(*aux_re, aux_im, mat_op, PMLIntegKind::FLOQUET_DIFFUSION, c_mass, 0.0,
                filter, omega);
    }
  }
}

// Pencil coefficients of the frequency-dependent PML terms A2(ω), with the optionally
// shifted mass term for preconditioning.
PMLPencil ExtraSystemPMLPencil(std::complex<double> omega, bool shifted = false,
                               bool shifted_floquet = false)
{
  return {1.0, 1i * omega, ShiftedMassCoefficient(-omega * omega, shifted),
          ShiftedMassCoefficient(-omega * omega, shifted_floquet)};
}

// Assemble a complex-valued operator from PML integrators only (null if empty on all
// processes).
std::unique_ptr<ComplexOperator> AssemblePMLOperator(const FiniteElementSpace &fespace,
                                                     const mfem::Array<int> &dbc_tdof_list,
                                                     std::vector<PMLIntegrator> &pml_re,
                                                     std::vector<PMLIntegrator> &pml_im,
                                                     Operator::DiagonalPolicy diag_policy)
{
  int empty[2] = {pml_re.empty(), pml_im.empty()};
  Mpi::GlobalMin(2, empty, fespace.GetComm());
  if (empty[0] && empty[1])
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  std::unique_ptr<Operator> ar, ai;
  if (!empty[0])
  {
    ar = AssembleOperator(AssemblyDomains(), fespace, nullptr, nullptr, nullptr, nullptr,
                          nullptr, nullptr, nullptr, skip_zeros, false, &pml_re);
  }
  if (!empty[1])
  {
    ai = AssembleOperator(AssemblyDomains(), fespace, nullptr, nullptr, nullptr, nullptr,
                          nullptr, nullptr, nullptr, skip_zeros, false, &pml_im);
  }
  auto A = std::make_unique<ComplexParOperator>(std::move(ar), std::move(ai), fespace);
  A->SetEssentialTrueDofs(dbc_tdof_list, diag_policy);
  return A;
}

// Add the frequency-dependent PML terms of A2(ω), with the stretch evaluated at ω, to the
// real and imaginary parts ar and ai of the other terms of A2(ω) (null if there are none).
void AddExtraSystemPMLOperators(const FiniteElementSpace &fespace,
                                const MaterialOperator &mat_op, std::complex<double> omega,
                                std::unique_ptr<Operator> &ar,
                                std::unique_ptr<Operator> &ai)
{
  std::vector<PMLIntegrator> pml_re, pml_im;
  AppendPMLPencil(pml_re, &pml_im, nullptr, nullptr, mat_op, PMLFilter::FREQUENCY_DEPENDENT,
                  omega, ExtraSystemPMLPencil(omega));
  int empty[2] = {pml_re.empty(), pml_im.empty()};
  Mpi::GlobalMin(2, empty, fespace.GetComm());
  auto AddPML = [&fespace](std::unique_ptr<Operator> &a, std::vector<PMLIntegrator> &pml)
  {
    constexpr bool skip_zeros = false;
    auto a_pml =
        AssembleOperator(AssemblyDomains(), fespace, nullptr, nullptr, nullptr, nullptr,
                         nullptr, nullptr, nullptr, skip_zeros, false, &pml);
    if (a)
    {
      a = std::make_unique<SumOperator>(std::move(a), std::move(a_pml));
    }
    else
    {
      a = std::move(a_pml);
    }
  };
  if (!empty[0])
  {
    AddPML(ar, pml_re);
  }
  if (!empty[1])
  {
    AddPML(ai, pml_im);
  }
}

// The PML terms are complex-valued (and only set up for frequency domain problems).
template <OperatorType OperType>
void CheckPMLOperatorType(const MaterialOperator &mat_op)
{
  if constexpr (!std::is_same_v<OperType, ComplexOperator>)
  {
    MFEM_VERIFY(!mat_op.HasPML(), "PML regions require complex-valued operators!");
  }
}

}  // namespace

void SpaceOperator::AssembleFrequencyDependentPermittivityA2Operators(
    std::complex<double> omega, const MaterialPropertyCoefficient &dfbr,
    const MaterialPropertyCoefficient &dfbi, const MaterialPropertyCoefficient &fbr,
    const MaterialPropertyCoefficient &fbi, const MaterialPropertyCoefficient &frbr,
    const MaterialPropertyCoefficient &frbi, const MaterialPropertyCoefficient &fpbr,
    const MaterialPropertyCoefficient &fpbi, std::unique_ptr<Operator> &ar,
    std::unique_ptr<Operator> &ai, bool include_poles)
{
  // Preserve the established boundary-only assembly path when no volume A2 model exists.
  // In particular, exact-zero poles have already been folded into ordinary conductivity.
  if (!mat_op.HasFrequencyDependentPermittivityA2())
  {
    int empty[2] = {AreExactlyZero(dfbr, fbr, frbr, fpbr),
                    AreExactlyZero(dfbi, fbi, frbi, fpbi)};
    Mpi::GlobalMin(2, empty, GetComm());
    constexpr bool skip_zeros = false;
    if (!empty[0])
    {
      ar = AssembleOperator(assembly_domains, GetNDSpace(), nullptr, nullptr, &dfbr, &fbr,
                            nullptr, &frbr, &fpbr, skip_zeros);
    }
    if (!empty[1])
    {
      ai = AssembleOperator(assembly_domains, GetNDSpace(), nullptr, nullptr, &dfbi, &fbi,
                            nullptr, &frbi, &fpbi, skip_zeros);
    }
    return;
  }

  // Always build both slots, including a zero wrapper when the contribution cancels at
  // this frequency: nonlinear interpolation and PROM fallback selection require the same
  // operator structure at every frequency.
  std::vector<double> model_A2_real, model_A2_imag;
  GetFrequencyDependentPermittivityA2Coefficients(omega, model_A2_real, model_A2_imag,
                                                  include_poles);
  constexpr bool skip_zeros = false;
  std::unique_ptr<Operator> base_r, base_i;
  if (!AreExactlyZero(dfbr, fbr, frbr, fpbr))
  {
    base_r = AssembleOperator(assembly_domains, GetNDSpace(), nullptr, nullptr, &dfbr, &fbr,
                              nullptr, &frbr, &fpbr, skip_zeros);
  }
  if (!AreExactlyZero(dfbi, fbi, frbi, fpbi))
  {
    base_i = AssembleOperator(assembly_domains, GetNDSpace(), nullptr, nullptr, &dfbi, &fbi,
                              nullptr, &frbi, &fpbi, skip_zeros);
  }
  ar = BuildFrequencyDependentPermittivityA2Operator(std::move(base_r), model_A2_real);
  ai = BuildFrequencyDependentPermittivityA2Operator(std::move(base_i), model_A2_imag);
}

bool SpaceOperator::HasSymmetricStiffnessMatrix() const
{
  if (!mat_op.HasWaveVector())
  {
    return true;
  }
  const auto &farfield_op = surface_terms.GetFarfieldOp();
  return mat_op.HasFloquetFrequencyScaling() &&
         (farfield_op.GetAttrList().Size() == 0 || farfield_op.GetOrder() < 2);
}

template <OperatorType OperType>
std::unique_ptr<OperType>
SpaceOperator::GetStiffnessMatrix(Operator::DiagonalPolicy diag_policy)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  CheckPMLOperatorType<OperType>(mat_op);
  MaterialPropertyCoefficient df(mat_op.MaxCeedAttribute()), f(mat_op.MaxCeedAttribute()),
      fb(mat_op.MaxCeedBdrAttribute()), fc(mat_op.MaxCeedAttribute()),
      fpb(mat_op.MaxCeedBdrAttribute());
  AddStiffnessCoefficients(1.0, df, f);
  AddStiffnessBdrCoefficients(1.0, fb);
  if (!mat_op.HasFloquetFrequencyScaling())
  {
    AddRealPeriodicCoefficients(1.0, f);
    AddImagPeriodicCoefficients(1.0, fc);
  }
  else
  {
    AddStiffnessPeriodicBdrCoefficients(1.0, fpb);
  }
  RestrictToAssemblyDomains({&df, &f, &fc}, {&fb, &fpb});

  // Static PML stretch (a frequency-dependent stretch contributes to GetExtraSystemMatrix).
  std::vector<PMLIntegrator> pml_re, pml_im;
  AppendPMLPencil(pml_re, &pml_im, nullptr, nullptr, mat_op, PMLFilter::STATIC, 0.0,
                  {.a0 = 1.0});

  int empty[2] = {AreExactlyZero(df, f, fb, fpb) && pml_re.empty(),
                  fc.IsExactlyZero() && pml_im.empty()};
  Mpi::GlobalMin(2, empty, GetComm());
  if (empty[0] && empty[1])
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  std::unique_ptr<Operator> kr, ki;
  if (!empty[0])
  {
    kr = AssembleOperator(assembly_domains, GetNDSpace(), &df, &f, nullptr, &fb, nullptr,
                          nullptr, &fpb, skip_zeros, false, &pml_re);
  }
  if (!empty[1])
  {
    // Imaginary branch: periodic imaginary coefficient (fc) plus PML imag μ̃⁻¹ curl-curl.
    ki = AssembleOperator(assembly_domains, GetNDSpace(), nullptr, nullptr, nullptr,
                          nullptr, &fc, nullptr, nullptr, skip_zeros, false, &pml_im);
  }
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
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

template <OperatorType OperType>
std::unique_ptr<OperType>
SpaceOperator::GetDampingMatrix(Operator::DiagonalPolicy diag_policy)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  CheckPMLOperatorType<OperType>(mat_op);
  MaterialPropertyCoefficient f(mat_op.MaxCeedAttribute()),
      fb(mat_op.MaxCeedBdrAttribute());
  MaterialPropertyCoefficient fp(mat_op.MaxCeedAttribute()),
      frb(mat_op.MaxCeedBdrAttribute());
  AddDampingCoefficients(1.0, f);
  AddDampingBdrCoefficients(1.0, fb);
  if (mat_op.HasFloquetFrequencyScaling())
  {
    AddImagPeriodicCoefficients(1.0, fp);
    AddDampingPeriodicBdrCoefficients(1.0, frb);
  }
  RestrictToAssemblyDomains({&f, &fp}, {&fb, &frb});
  // Static PML stretch.
  std::vector<PMLIntegrator> pml_re, pml_im;
  AppendPMLPencil(pml_re, &pml_im, nullptr, nullptr, mat_op, PMLFilter::STATIC, 0.0,
                  {.a1 = 1.0});
  int empty[2] = {AreExactlyZero(f, fb, fp, frb) && pml_re.empty(), pml_im.empty()};
  Mpi::GlobalMin(2, empty, GetComm());
  if (empty[0] && empty[1])
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  std::unique_ptr<Operator> cr, ci;
  if (!empty[0])
  {
    cr = AssembleOperator(assembly_domains, GetNDSpace(), nullptr, &f, nullptr, &fb, &fp,
                          &frb, nullptr, skip_zeros, false, &pml_re);
  }
  if (!empty[1])
  {
    ci = AssembleOperator(assembly_domains, GetNDSpace(), nullptr, nullptr, nullptr,
                          nullptr, nullptr, nullptr, nullptr, skip_zeros, false, &pml_im);
  }
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    auto C =
        std::make_unique<ComplexParOperator>(std::move(cr), std::move(ci), GetNDSpace());
    C->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return C;
  }
  else
  {
    MFEM_VERIFY(!ci, "Unexpected imaginary part in GetDampingMatrix<Operator>!");
    auto C = std::make_unique<ParOperator>(std::move(cr), GetNDSpace());
    C->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return C;
  }
}

template <OperatorType OperType>
std::unique_ptr<OperType> SpaceOperator::GetMassMatrix(Operator::DiagonalPolicy diag_policy)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  CheckPMLOperatorType<OperType>(mat_op);
  MaterialPropertyCoefficient fr(mat_op.MaxCeedAttribute()), fi(mat_op.MaxCeedAttribute()),
      fbr(mat_op.MaxCeedBdrAttribute()), fbi(mat_op.MaxCeedBdrAttribute());
  AddRealMassCoefficients(1.0, fr);
  AddRealMassBdrCoefficients(1.0, fbr);
  if (mat_op.HasFloquetFrequencyScaling())
  {
    AddRealPeriodicCoefficients(-1.0, fr);
  }
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    AddImagMassCoefficients(1.0, fi);
  }
  RestrictToAssemblyDomains({&fr, &fi}, {&fbr, &fbi});
  // Static PML stretch.
  std::vector<PMLIntegrator> pml_re, pml_im;
  AppendPMLPencil(pml_re, &pml_im, nullptr, nullptr, mat_op, PMLFilter::STATIC, 0.0,
                  {.a2 = 1.0, .a2_floquet = 1.0});
  int empty[2] = {AreExactlyZero(fr, fbr) && pml_re.empty(),
                  AreExactlyZero(fi, fbi) && pml_im.empty()};
  Mpi::GlobalMin(2, empty, GetComm());
  if (empty[0] && empty[1])
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  std::unique_ptr<Operator> mr, mi;
  if (!empty[0])
  {
    mr = AssembleOperator(assembly_domains, GetNDSpace(), nullptr, &fr, nullptr, &fbr,
                          nullptr, nullptr, nullptr, skip_zeros, false, &pml_re);
  }
  if (!empty[1])
  {
    mi = AssembleOperator(assembly_domains, GetNDSpace(), nullptr, &fi, nullptr, &fbi,
                          nullptr, nullptr, nullptr, skip_zeros, false, &pml_im);
  }
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
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

template <OperatorType OperType>
std::unique_ptr<OperType>
SpaceOperator::GetExtraSystemMatrix(double omega, Operator::DiagonalPolicy diag_policy)
{
  return GetExtraSystemMatrix<OperType>(omega, diag_policy, /*include_wave_ports=*/true);
}

template <OperatorType OperType>
std::unique_ptr<OperType>
SpaceOperator::GetExtraSystemMatrix(double omega, Operator::DiagonalPolicy diag_policy,
                                    bool include_wave_ports, bool include_pml)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient dfbr(mat_op.MaxCeedBdrAttribute()),
      dfbi(mat_op.MaxCeedBdrAttribute()), fbr(mat_op.MaxCeedBdrAttribute()),
      fbi(mat_op.MaxCeedBdrAttribute()), frbr(mat_op.MaxCeedBdrAttribute()),
      frbi(mat_op.MaxCeedBdrAttribute()), fpbr(mat_op.MaxCeedBdrAttribute()),
      fpbi(mat_op.MaxCeedBdrAttribute());
  AddExtraSystemBdrCoefficients(omega, dfbr, dfbi, fbr, fbi, frbr, frbi, &fpbr, &fpbi,
                                include_wave_ports);
  RestrictToAssemblyDomains({}, {&dfbr, &dfbi, &fbr, &fbi, &frbr, &frbi, &fpbr, &fpbi});
  std::unique_ptr<Operator> ar, ai;
  AssembleFrequencyDependentPermittivityA2Operators(omega, dfbr, dfbi, fbr, fbi, frbr, frbi,
                                                    fpbr, fpbi, ar, ai);

  // Frequency-dependent PML stretch, evaluated at ω.
  CheckPMLOperatorType<OperType>(mat_op);
  if (include_pml)
  {
    AddExtraSystemPMLOperators(GetNDSpace(), mat_op, omega, ar, ai);
  }
  if (!ar && !ai)
  {
    return {};
  }
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
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

std::unique_ptr<ComplexOperator>
SpaceOperator::GetExtraSystemMatrix(std::complex<double> omega,
                                    Operator::DiagonalPolicy diag_policy,
                                    bool include_permittivity_poles, bool include_pml)
{
  auto ports =
      wave_port_op.EvaluateComplex(omega, false, GetNDSpace(), nd_dbc_tdof_lists.back());
  MFEM_VERIFY(ports.available,
              "Wave port " << ports.port_idx
                           << " eigensolver did not converge at complex ω = " << omega
                           << " (" << ports.num_converged << " of " << ports.requested
                           << " modes)!");
  return BuildExtraSystemMatrix(omega, ports.kn, diag_policy, include_permittivity_poles,
                                include_pml);
}

std::unique_ptr<ComplexOperator> SpaceOperator::BuildExtraSystemMatrix(
    std::complex<double> omega,
    const std::vector<std::pair<int, std::complex<double>>> &port_kn,
    Operator::DiagonalPolicy diag_policy, bool include_permittivity_poles, bool include_pml)
{
  // Complex-ω A2(λ) for the eigenmode nonlinear solve: nonzero-pole volume terms,
  // frequency-dependent PML terms, and frequency-dependent boundary terms (2nd-order ABC,
  // surface conductivity, rational impedance, numeric wave ports) are evaluated at the
  // genuinely complex frequency (ω = -i·λ). Always returns a ComplexOperator since these
  // terms carry a real-slot contribution at complex ω.
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient dfbr(mat_op.MaxCeedBdrAttribute()),
      dfbi(mat_op.MaxCeedBdrAttribute()), fbr(mat_op.MaxCeedBdrAttribute()),
      fbi(mat_op.MaxCeedBdrAttribute()), frbr(mat_op.MaxCeedBdrAttribute()),
      frbi(mat_op.MaxCeedBdrAttribute()), fpbr(mat_op.MaxCeedBdrAttribute()),
      fpbi(mat_op.MaxCeedBdrAttribute());
  AddExtraSystemBdrCoefficients(omega, port_kn, dfbr, dfbi, fbr, fbi, frbr, frbi, &fpbr,
                                &fpbi);
  RestrictToAssemblyDomains({}, {&dfbr, &dfbi, &fbr, &fbi, &frbr, &frbi, &fpbr, &fpbi});
  std::unique_ptr<Operator> ar, ai;
  AssembleFrequencyDependentPermittivityA2Operators(omega, dfbr, dfbi, fbr, fbi, frbr, frbi,
                                                    fpbr, fpbi, ar, ai,
                                                    include_permittivity_poles);
  if (include_pml)
  {
    AddExtraSystemPMLOperators(GetNDSpace(), mat_op, omega, ar, ai);
  }
  if (!ar && !ai)
  {
    return {};
  }
  auto A = std::make_unique<ComplexParOperator>(std::move(ar), std::move(ai), GetNDSpace());
  A->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
  return A;
}

std::array<std::unique_ptr<ComplexOperator>, 3>
SpaceOperator::GetFrozenPMLMatrices(double omega, Operator::DiagonalPolicy diag_policy)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  std::array<std::unique_ptr<ComplexOperator>, 3> ops;
  const std::array<PMLPencil, 3> terms = {
      {{.a0 = 1.0}, {.a1 = 1.0}, {.a2 = 1.0, .a2_floquet = 1.0}}};
  for (std::size_t k = 0; k < terms.size(); k++)
  {
    std::vector<PMLIntegrator> pml_re, pml_im;
    AppendPMLPencil(pml_re, &pml_im, nullptr, nullptr, mat_op,
                    PMLFilter::FREQUENCY_DEPENDENT, omega, terms[k]);
    ops[k] = AssemblePMLOperator(GetNDSpace(), nd_dbc_tdof_lists.back(), pml_re, pml_im,
                                 diag_policy);
  }
  return ops;
}

std::unique_ptr<ComplexOperator>
SpaceOperator::GetExtraSystemPMLMatrix(std::complex<double> omega,
                                       Operator::DiagonalPolicy diag_policy)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  std::vector<PMLIntegrator> pml_re, pml_im;
  AppendPMLPencil(pml_re, &pml_im, nullptr, nullptr, mat_op, PMLFilter::FREQUENCY_DEPENDENT,
                  omega, ExtraSystemPMLPencil(omega));
  return AssemblePMLOperator(GetNDSpace(), nd_dbc_tdof_lists.back(), pml_re, pml_im,
                             diag_policy);
}

template <OperatorType OperType>
std::unique_ptr<OperType>
SpaceOperator::GetWavePortBoundaryMassMatrix(int port_idx,
                                             Operator::DiagonalPolicy diag_policy)
{
  // Per-port μ⁻¹ boundary mass matrix, ω-independent — see
  // WavePortOperator::AddBoundaryMassBdrCoefficients. Pure imaginary part of A2(ω) when
  // the per-ω scalar i·k_n(ω) is reattached.
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient fb(mat_op.MaxCeedBdrAttribute());
  wave_port_op.AddBoundaryMassBdrCoefficients(port_idx, fb);
  int empty = fb.IsExactlyZero();
  Mpi::GlobalMin(1, &empty, GetComm());
  if (empty)
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  auto m = AssembleOperator(assembly_domains, GetNDSpace(), nullptr, nullptr, nullptr, &fb,
                            nullptr, nullptr, nullptr, skip_zeros);
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    auto M_op = std::make_unique<ComplexParOperator>(nullptr, std::move(m), GetNDSpace());
    M_op->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M_op;
  }
  else
  {
    auto M_op = std::make_unique<ParOperator>(std::move(m), GetNDSpace());
    M_op->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M_op;
  }
}

template <OperatorType OperType>
std::unique_ptr<OperType>
SpaceOperator::GetFarfieldBoundaryCurlCurlMatrix(Operator::DiagonalPolicy diag_policy,
                                                 bool imag_slot)
{
  // ω-independent boundary curl-curl matrix M_ff for the 2nd-order farfield ABC, with unit
  // scale: the term at frequency ω is β·B = i·(0.5/ω)·M_ff, for β = 1/(2iω) and the
  // boundary curl-curl form B (see FarfieldBoundaryOperator::AddCurlCurlBdrCoefficients),
  // so M_ff = -B includes the negative sign of β and the scalings below are positive. With
  // a Floquet wave vector which does not depend on frequency, B is Hermitian and M_ff is
  // complex. By default M_ff is stored as is (the curl-curl term on the REAL slot of the
  // resulting ComplexParOperator) so that downstream BuildParSumOperator can scale it by an
  // arbitrary complex coefficient (the real-ω path uses i·(0.5/ω); the complex-λ path uses
  // -0.5/λ). With imag_slot=true i·M_ff is stored instead (the curl-curl term on the
  // IMAGINARY slot), matching the wave-port boundary-mass convention (the i in i·f(ω)·M
  // baked in) so circuit synthesis can fold it in uniformly. Returns null if the farfield
  // ABC order < 2 or it contributes no DoFs on any rank (the check is collective, so the
  // null contract is rank-uniform).
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient dfr(mat_op.MaxCeedBdrAttribute()),
      dfi(mat_op.MaxCeedBdrAttribute()), frr(mat_op.MaxCeedBdrAttribute()),
      fri(mat_op.MaxCeedBdrAttribute()), fpr(mat_op.MaxCeedBdrAttribute()),
      fpi(mat_op.MaxCeedBdrAttribute());
  surface_terms.GetFarfieldOp().AddCurlCurlBdrCoefficients(
      imag_slot ? std::complex<double>(0.0, -1.0) : -1.0, dfr, dfi, frr, fri, &fpr, &fpi);
  int empty[2] = {AreExactlyZero(dfr, frr, fpr), AreExactlyZero(dfi, fri, fpi)};
  Mpi::GlobalMin(2, empty, GetComm());
  if (empty[0] && empty[1])
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  std::unique_ptr<Operator> mr, mi;
  if (!empty[0])
  {
    mr = AssembleOperator(assembly_domains, GetNDSpace(), nullptr, nullptr, &dfr, nullptr,
                          nullptr, &frr, &fpr, skip_zeros);
  }
  if (!empty[1])
  {
    mi = AssembleOperator(assembly_domains, GetNDSpace(), nullptr, nullptr, &dfi, nullptr,
                          nullptr, &fri, &fpi, skip_zeros);
  }
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    auto M_op =
        std::make_unique<ComplexParOperator>(std::move(mr), std::move(mi), GetNDSpace());
    M_op->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M_op;
  }
  else
  {
    MFEM_VERIFY(!imag_slot,
                "imag_slot is only meaningful for the ComplexOperator instantiation of "
                "GetFarfieldBoundaryCurlCurlMatrix!");
    MFEM_VERIFY(
        !mi, "Unexpected imaginary part in GetFarfieldBoundaryCurlCurlMatrix<Operator>!");
    auto M_op = std::make_unique<ParOperator>(std::move(mr), GetNDSpace());
    M_op->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M_op;
  }
}

template <OperatorType OperType>
std::unique_ptr<OperType>
SpaceOperator::GetExtraSystemBdrMassTermMatrix(int i, Operator::DiagonalPolicy diag_policy,
                                               bool imag_slot)
{
  // ω-independent boundary mass M_i of the frequency-dependent boundary mass term i of
  // SurfaceTerms, with its unit coefficient (including crack scaling and the two surfaces
  // of interior conducting sheets). By default stored on the REAL slot, so that its
  // contribution to A2 at ω is GetExtraSystemBdrMassTermScale(i, ω)·M_i, and a caller can
  // scale it by an arbitrary complex coefficient (the NLEPS HYBRID fit-or-freeze seed
  // strategy). With imag_slot=true it is placed on the IMAGINARY slot, matching the
  // wave-port boundary-mass convention (the i in i·f(ω)·M baked in, with f(ω) =
  // scale(ω)/i), so the adaptive sweep and circuit synthesis treat it uniformly. Returns
  // null if the term contributes no DoFs on any rank (the check is collective, so the null
  // contract is rank-uniform).
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient fb(mat_op.MaxCeedBdrAttribute());
  surface_terms.AddExtraSystemBdrMassTermCoefficients(i, fb);
  int empty = fb.IsExactlyZero();
  Mpi::GlobalMin(1, &empty, GetComm());
  if (empty)
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  auto m = AssembleOperator(assembly_domains, GetNDSpace(), nullptr, nullptr, nullptr, &fb,
                            nullptr, nullptr, nullptr, skip_zeros);
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    auto M_op =
        imag_slot
            ? std::make_unique<ComplexParOperator>(nullptr, std::move(m), GetNDSpace())
            : std::make_unique<ComplexParOperator>(std::move(m), nullptr, GetNDSpace());
    M_op->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M_op;
  }
  else
  {
    MFEM_VERIFY(!imag_slot,
                "imag_slot is only meaningful for the ComplexOperator instantiation of "
                "GetExtraSystemBdrMassTermMatrix!");
    auto M_op = std::make_unique<ParOperator>(std::move(m), GetNDSpace());
    M_op->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M_op;
  }
}

std::unique_ptr<ComplexOperator>
SpaceOperator::GetFarfieldDispersiveCurlCurlMatrix(std::size_t k,
                                                   Operator::DiagonalPolicy diag_policy)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient dfr(mat_op.MaxCeedBdrAttribute()),
      dfi(mat_op.MaxCeedBdrAttribute()), frr(mat_op.MaxCeedBdrAttribute()),
      fri(mat_op.MaxCeedBdrAttribute()), fpr(mat_op.MaxCeedBdrAttribute()),
      fpi(mat_op.MaxCeedBdrAttribute());
  // i M_k with M_k = -(C + T + i(X - Xᵀ)), as for nondispersive M_ff.
  surface_terms.GetFarfieldOp().AddDispersiveCurlCurlBdrCoefficients(
      k, std::complex<double>(0.0, -1.0), dfr, dfi, frr, fri, &fpr, &fpi);
  int empty[2] = {AreExactlyZero(dfr, frr, fpr), AreExactlyZero(dfi, fri, fpi)};
  Mpi::GlobalMin(2, empty, GetComm());
  if (empty[0] && empty[1])
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  std::unique_ptr<Operator> mr, mi;
  if (!empty[0])
  {
    mr = AssembleOperator(assembly_domains, GetNDSpace(), nullptr, nullptr, &dfr, nullptr,
                          nullptr, &frr, &fpr, skip_zeros);
  }
  if (!empty[1])
  {
    mi = AssembleOperator(assembly_domains, GetNDSpace(), nullptr, nullptr, &dfi, nullptr,
                          nullptr, &fri, &fpi, skip_zeros);
  }
  auto M_op =
      std::make_unique<ComplexParOperator>(std::move(mr), std::move(mi), GetNDSpace());
  M_op->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
  return M_op;
}

std::unique_ptr<ComplexOperator> SpaceOperator::GetFarfieldDispersiveDampingFloquetMatrix(
    std::size_t k, Operator::DiagonalPolicy diag_policy)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient fr(mat_op.MaxCeedBdrAttribute());
  surface_terms.GetFarfieldOp().AddDispersiveDampingFloquetBdrCoefficients(k, 1.0, fr);
  int empty = fr.IsExactlyZero();
  Mpi::GlobalMin(1, &empty, GetComm());
  if (empty)
  {
    return {};
  }
  auto m = AssembleOperator(assembly_domains, GetNDSpace(), nullptr, nullptr, nullptr,
                            nullptr, nullptr, &fr, nullptr, false);
  auto M = std::make_unique<ComplexParOperator>(nullptr, std::move(m), GetNDSpace());
  M->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
  return M;
}

std::unique_ptr<ComplexOperator> SpaceOperator::GetFarfieldDispersiveStiffnessFloquetMatrix(
    std::size_t k, Operator::DiagonalPolicy diag_policy)
{
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MaterialPropertyCoefficient fp(mat_op.MaxCeedBdrAttribute());
  surface_terms.GetFarfieldOp().AddDispersiveStiffnessFloquetBdrCoefficients(k, 1.0, fp);
  int empty = fp.IsExactlyZero();
  Mpi::GlobalMin(1, &empty, GetComm());
  if (empty)
  {
    return {};
  }
  auto m = AssembleOperator(assembly_domains, GetNDSpace(), nullptr, nullptr, nullptr,
                            nullptr, nullptr, nullptr, &fp, false);
  auto M = std::make_unique<ComplexParOperator>(std::move(m), nullptr, GetNDSpace());
  M->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
  return M;
}

std::vector<SpaceOperator::FarfieldPoleTerm>
SpaceOperator::GetFarfieldPoleTerms(Operator::DiagonalPolicy diag_policy)
{
  // A2(ω) holds i·(0.5/ω)·M_ff (real slot) and i·(0.5/ω)/r_k(ω)·M_k, with M_k stored as
  // i·M_k (imaginary slot), so with ω = -iλ the scalars are -0.5/λ and 0.5/(ω r_k(ω)).
  std::vector<FarfieldPoleTerm> terms;
  if (auto M_ff = GetFarfieldBoundaryCurlCurlMatrix<ComplexOperator>(diag_policy))
  {
    terms.push_back(
        {std::move(M_ff), [](std::complex<double> lambda) { return -0.5 / lambda; }});
  }
  // The scalars refer to the farfield operator of this SpaceOperator, which outlives them.
  const auto &farfield_op = surface_terms.GetFarfieldOp();
  for (std::size_t k = 0; k < farfield_op.NumDispersiveMaterials(); k++)
  {
    if (auto M_k = GetFarfieldDispersiveCurlCurlMatrix(k, diag_policy))
    {
      terms.push_back({std::move(M_k), [&ff = farfield_op, k](std::complex<double> lambda)
                       {
                         const std::complex<double> omega =
                             lambda / std::complex<double>(0.0, 1.0);
                         return 0.5 / (omega * ff.GetDispersiveScale(k, omega));
                       }});
    }
    // The frequency-scaled rotated trace and mixed curl also carry 1/r_k(ω).
    // Share their matrix owners with ROM and retain the distinct ω powers in the seed.
    for (bool rotated : {true, false})
    {
      auto M = rotated ? GetFarfieldDispersiveDampingFloquetMatrix(k, diag_policy)
                       : GetFarfieldDispersiveStiffnessFloquetMatrix(k, diag_policy);
      if (M)
      {
        terms.push_back({std::move(M),
                         [&ff = farfield_op, k, rotated](std::complex<double> lambda)
                         {
                           const auto omega = lambda / std::complex<double>(0.0, 1.0);
                           return (rotated ? omega : 1.0) / ff.GetDispersiveScale(k, omega);
                         }});
      }
    }
  }
  return terms;
}

template <OperatorType OperType>
std::unique_ptr<OperType>
SpaceOperator::GetFloquetRobinBoundaryMassMatrix(int port_idx,
                                                 Operator::DiagonalPolicy diag_policy)
{
  // ω-independent µ⁻¹ boundary mass for a single Floquet port's Robin BC, on the
  // IMAGINARY slot (matching the wave-port / farfield convention: the i in i·γ₀·M is baked
  // in). The full online contribution is i·γ₀(ω)·M_floquet, with γ₀ the specular (0,0)
  // propagation constant. Returns null if the port contributes no DoFs on any rank (the
  // check is collective, so the null contract is rank-uniform).
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  const auto &port = floquet_port_op.GetPort(port_idx);
  MaterialPropertyCoefficient fb(mat_op.MaxCeedBdrAttribute());
  MaterialPropertyCoefficient muinv_func(mat_op.GetBdrAttributeToMaterial(),
                                         mat_op.GetInvPermeability());
  muinv_func.RestrictCoefficient(mat_op.GetCeedBdrAttributes(port.GetAttrList()));
  fb.AddCoefficient(muinv_func.GetAttributeToMaterial(), muinv_func.GetMaterialProperties(),
                    1.0);
  int empty = fb.IsExactlyZero();
  Mpi::GlobalMin(1, &empty, GetComm());
  if (empty)
  {
    return {};
  }
  constexpr bool skip_zeros = false;
  auto m = AssembleOperator(assembly_domains, GetNDSpace(), nullptr, nullptr, nullptr, &fb,
                            nullptr, nullptr, nullptr, skip_zeros);
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    auto M_op = std::make_unique<ComplexParOperator>(nullptr, std::move(m), GetNDSpace());
    M_op->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M_op;
  }
  else
  {
    auto M_op = std::make_unique<ParOperator>(std::move(m), GetNDSpace());
    M_op->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), diag_policy);
    return M_op;
  }
}

template <typename OperType, typename ScalarType>
std::unique_ptr<OperType>
SpaceOperator::GetSystemMatrix(ScalarType a0, ScalarType a1, ScalarType a2,
                               const OperType *K, const OperType *C, const OperType *M,
                               const OperType *A2)
{
  // When A2 is an abstract ComplexOperator (not a sparse ComplexParOperator), it cannot
  // participate in BuildParSumOperator. Build the sparse KCM sum separately and add A2 via
  // a non-owning SumComplexOperator.
  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    if (A2 && !dynamic_cast<const ComplexParOperator *>(A2))
    {
      auto A_KCM = BuildParSumOperator({a0, a1, a2}, {K, C, M});
      return std::make_unique<SumComplexOperator>(std::move(A_KCM), *A2);
    }
  }
  return BuildParSumOperator({a0, a1, a2, ScalarType{1}}, {K, C, M, A2});
}

std::unique_ptr<ComplexOperator>
SpaceOperator::GetExtraSystemOperator(double omega, Operator::DiagonalPolicy diag_policy)
{
  // Wave ports contribute to the applied operator as the sparse local boundary mass i·k_n·M
  // (scalar-admittance n×H) plus the matrix-free modal correction W = Σ (W_full −
  // W_scalar), which supplies the ∇ₜEₙ term the sparse mass drops so the port reproduces
  // the full modal n×H on its own mode. The sparse mass is retained here
  // (include_wave_ports=true) and also serves as the preconditioner's port block; the
  // correction and the Floquet DtN are summed onto the sparse A2 via non-owning/owning
  // SumComplexOperator.
  auto A2 = GetExtraSystemMatrix<ComplexOperator>(omega, diag_policy,
                                                  /*include_wave_ports=*/true);
  auto W = wave_port_op.GetModalCorrectionOperator(omega, GetNDSpace(),
                                                   nd_dbc_tdof_lists.back());
  auto F = floquet_port_op.GetExtraSystemOperator(omega);
  std::unique_ptr<ComplexOperator> op;
  auto add = [&op](std::unique_ptr<ComplexOperator> term)
  {
    if (!term)
    {
      return;
    }
    op = op ? std::make_unique<SumComplexOperator>(std::move(op), std::move(term))
            : std::move(term);
  };
  add(std::move(A2));
  add(std::move(W));
  add(std::move(F));
  if (environment)
  {
    add(std::make_unique<EnvironmentTermOperator>(*environment, omega,
                                                  GetNDSpace().GetTrueVSize()));
  }
  return op;
}

std::unique_ptr<ComplexOperator>
SpaceOperator::GetExtraSystemOperator(std::complex<double> omega,
                                      Operator::DiagonalPolicy diag_policy)
{
  auto ports =
      wave_port_op.EvaluateComplex(omega, true, GetNDSpace(), nd_dbc_tdof_lists.back());
  MFEM_VERIFY(ports.available,
              "Wave port " << ports.port_idx
                           << " eigensolver did not converge at complex ω = " << omega
                           << " (" << ports.num_converged << " of " << ports.requested
                           << " modes)!");
  return BuildExtraSystemOperator(omega, std::move(ports), diag_policy);
}

SpaceOperator::ExtraSystemOperatorResult
SpaceOperator::TryGetExtraSystemOperator(std::complex<double> omega,
                                         Operator::DiagonalPolicy diag_policy)
{
  auto ports =
      wave_port_op.EvaluateComplex(omega, true, GetNDSpace(), nd_dbc_tdof_lists.back());
  if (!ports.available)
  {
    return {false, ports.port_idx, ports.num_converged, ports.requested, nullptr};
  }
  return {true, -1, 0, 0, BuildExtraSystemOperator(omega, std::move(ports), diag_policy)};
}

std::unique_ptr<ComplexOperator>
SpaceOperator::BuildExtraSystemOperator(std::complex<double> omega,
                                        WavePortOperator::ComplexPortContribution &&ports,
                                        Operator::DiagonalPolicy diag_policy)
{
  // The sparse local wave-port mass i·k_n(ω)·M and the matrix-free modal correction W from
  // the same port modes, summed via an owning SumComplexOperator.
  auto A2 = BuildExtraSystemMatrix(omega, ports.kn, diag_policy);
  auto W = std::move(ports.W);
  if (!W)
  {
    return A2;
  }
  if (!A2)
  {
    return W;
  }
  return std::make_unique<SumComplexOperator>(std::move(A2), std::move(W));
}

std::unique_ptr<ComplexVector> SpaceOperator::GetWavePortModeVector(int port_idx,
                                                                    double omega)
{
  return wave_port_op.GetModeExcitationVector(port_idx, omega, GetNDSpace(),
                                              nd_dbc_tdof_lists.back());
}

std::vector<int> SpaceOperator::GetModalCorrectionSynthesisPorts(double omega_ref)
{
  return wave_port_op.GetModalCorrectionSynthesisPorts(omega_ref, GetNDSpace(),
                                                       nd_dbc_tdof_lists.back());
}

WavePortOperator::ModalCorrectionSample
SpaceOperator::SampleModalCorrectionVectors(int port_idx, std::complex<double> omega)
{
  return wave_port_op.SampleModalCorrectionVectors(port_idx, omega, GetNDSpace(),
                                                   nd_dbc_tdof_lists.back());
}

std::unique_ptr<Operator> SpaceOperator::GetInnerProductMatrix(double a0, double a2,
                                                               const ComplexOperator *K,
                                                               const ComplexOperator *M,
                                                               double omega)
{
  if (mat_op.HasPML())
  {
    // The real parts of the PML tensors are not positive definite (multiple stretched
    // directions in the edge and corner regions of the layer), and the frequency-dependent
    // PML terms are not in K or M. Instead, assemble the matrix with the PML tensors of the
    // magnitudes of the stretch factors, as for the real-valued preconditioner
    // approximation.
    MFEM_VERIFY(a0 >= 0.0 && a2 >= 0.0,
                "Invalid coefficients for the inner product matrix with PML!");
    MaterialPropertyCoefficient df(mat_op.MaxCeedAttribute()), f(mat_op.MaxCeedAttribute()),
        fb(mat_op.MaxCeedBdrAttribute());
    const double c0 = (K) ? a0 : 0.0, c2 = (M) ? a2 : 0.0;
    if (c0 > 0.0)
    {
      AddStiffnessCoefficients(c0, df, f);
      AddStiffnessBdrCoefficients(c0, fb);
      if (!mat_op.HasFloquetFrequencyScaling())
      {
        AddRealPeriodicCoefficients(c0, f);
      }
    }
    if (c2 > 0.0)
    {
      AddRealMassCoefficients(c2, f);
      AddRealMassBdrCoefficients(c2, fb);
      if (mat_op.HasFloquetFrequencyScaling())
      {
        AddRealPeriodicCoefficients(-c2, f);
      }
    }
    std::vector<PMLIntegrator> pml;
    const PMLPencil c = {.a0 = c0, .a2 = c2, .a2_floquet = c2};
    AppendPMLPencil(pml, nullptr, nullptr, nullptr, mat_op, PMLFilter::STATIC, 0.0, c);
    if (omega > 0.0)
    {
      AppendPMLPencil(pml, nullptr, nullptr, nullptr, mat_op,
                      PMLFilter::FREQUENCY_DEPENDENT, omega, c);
    }
    constexpr bool skip_zeros = false;
    auto b = AssembleOperator(assembly_domains, GetNDSpace(), &df, &f, nullptr, &fb,
                              nullptr, nullptr, nullptr, skip_zeros, false, &pml);
    return std::make_unique<ParOperator>(std::move(b), GetNDSpace());
  }
  const auto *PtAP_K = (K) ? dynamic_cast<const ComplexParOperator *>(K) : nullptr;
  const auto *PtAP_M = (M) ? dynamic_cast<const ComplexParOperator *>(M) : nullptr;
  return BuildParSumOperator(
      {a0, a2}, {PtAP_K ? PtAP_K->Real() : nullptr, PtAP_M ? PtAP_M->Real() : nullptr});
}

namespace
{

template <OperatorType OperType>
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

// L2-project a boundary VectorCoefficient onto the Nedelec tangent-trace space on the
// port surface: assemble the boundary mass M_bdr on port elements and CG-solve
// M_bdr * e = f. MFEM's ProjectBdrCoefficientTangent is a DOF-functional interpolant,
// not an L2 projection, so it returns a different vector (see PR #684 revert).
//
// TODO(future): Restrict the system to port-supported DOFs. The current full-T-space
// solve with DIAG_ONE on non-port DOFs is correct but pays for a trivially-identity
// block over most unknowns; a port-DOF-only system would be much smaller.
void ProjectBdrCoefficientViaMassSolve(SumVectorCoefficient &fb, const LumpedPortData &data,
                                       const MaterialOperator &mat_op,
                                       FiniteElementSpace &nd_fespace, MPI_Comm comm,
                                       Vector &result)
{
  // Assemble the boundary linear form f = ∫ φ_i · coeff dS (parallel-correct via P^T).
  mfem::LinearForm rhs(&nd_fespace.Get());
  rhs.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fb));
  rhs.UseFastAssembly(false);
  rhs.UseDevice(false);
  rhs.Assemble();
  rhs.UseDevice(true);
  nd_fespace.GetProlongationMatrix()->MultTranspose(rhs, result);

  // Assemble boundary mass matrix M_bdr = ∫ φ_i · φ_j dS on port surface.
  MaterialPropertyCoefficient fb_mass(mat_op.MaxCeedBdrAttribute());
  for (const auto &elem : data.elems)
  {
    fb_mass.AddMaterialProperty(mat_op.GetCeedBdrAttributes(elem->GetAttrList()), 1.0);
  }
  BilinearForm m_bdr(nd_fespace);
  if (!fb_mass.IsExactlyZero())
  {
    m_bdr.AddBoundaryIntegrator<VectorFEMassIntegrator>(fb_mass);
  }
  auto M_bdr = std::make_unique<ParOperator>(m_bdr.Assemble(false), nd_fespace);

  // M_bdr is zero off the port, so set non-port DOFs essential with DIAG_ONE to get a
  // full-rank SPD system on the full T-space. non_port_tdof_list must outlive M_bdr.
  mfem::Array<int> attr_list;
  for (const auto &elem : data.elems)
  {
    attr_list.Append(elem->GetAttrList());
  }
  const auto &mesh = nd_fespace.GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker;
  mesh::AttrToMarker(bdr_attr_max, attr_list, attr_marker);

  mfem::Array<int> port_tdof_list;
  nd_fespace.Get().GetEssentialTrueDofs(attr_marker, port_tdof_list);

  mfem::Array<int> non_port_tdof_list;
  {
    std::vector<bool> is_port(nd_fespace.GetTrueVSize(), false);
    for (int i = 0; i < port_tdof_list.Size(); i++)
    {
      is_port[port_tdof_list[i]] = true;
    }
    for (int i = 0; i < nd_fespace.GetTrueVSize(); i++)
    {
      if (!is_port[i])
      {
        non_port_tdof_list.Append(i);
      }
    }
  }
  M_bdr->SetEssentialTrueDofs(non_port_tdof_list, Operator::DIAG_ONE);

  // CG solve M_bdr * e = f entirely in T-vector space.
  // TODO: Make solver parameters configurable from IoData, or inherit other settings.
  auto pcg = std::make_unique<CgSolver<Operator>>(comm, 0);
  pcg->SetInitialGuess(false);
  pcg->SetRelTol(1.0e-14);
  pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
  pcg->SetMaxIter(200);
  auto jac = std::make_unique<JacobiSmoother<Operator>>(comm);
  auto ksp = std::make_unique<BaseKspSolver<Operator>>(std::move(pcg), std::move(jac));
  ksp->SetOperators(*M_bdr, *M_bdr);

  Vector sol(nd_fespace.GetTrueVSize());
  sol.UseDevice(true);
  sol = 0.0;
  ksp->Mult(result, sol);
  result = sol;
}

}  // namespace

template <typename A3Type>
void SpaceOperator::AssemblePreconditioner(
    std::complex<double> a0, std::complex<double> a1, std::complex<double> a2, A3Type a3,
    std::vector<std::unique_ptr<Operator>> &br_vec,
    std::vector<std::unique_ptr<Operator>> &br_aux_vec,
    std::vector<std::unique_ptr<Operator>> &bi_vec,
    std::vector<std::unique_ptr<Operator>> &bi_aux_vec,
    std::unique_ptr<ComplexOperator> *pml_block)
{
  constexpr bool skip_zeros = false;
  // Cache geometry/material tensors for repeated CPU smoother applications. Only the
  // fine hierarchy uses this data; the coarse sparse operator is assembled separately.
  const bool assemble_q_data = !mfem::Device::Allows(mfem::Backend::DEVICE_MASK);
  MaterialPropertyCoefficient dfr(mat_op.MaxCeedAttribute()),
      dfi(mat_op.MaxCeedAttribute()), fr(mat_op.MaxCeedAttribute()),
      fi(mat_op.MaxCeedAttribute()), dfbr(mat_op.MaxCeedBdrAttribute()),
      dfbi(mat_op.MaxCeedBdrAttribute()), fbr(mat_op.MaxCeedBdrAttribute()),
      fbi(mat_op.MaxCeedBdrAttribute()), fpi(mat_op.MaxCeedAttribute()),
      fpr(mat_op.MaxCeedAttribute()), frbr(mat_op.MaxCeedBdrAttribute()),
      frbi(mat_op.MaxCeedBdrAttribute()), fpbr(mat_op.MaxCeedBdrAttribute()),
      fpbi(mat_op.MaxCeedBdrAttribute());
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
  AddFrequencyDependentPermittivityA2Coefficients(a3, fr, fi);
  AddExtraSystemBdrCoefficients(a3, dfbr, dfbi, fbr, fbi, frbr, frbi, &fpbr, &fpbi);
  // dfbr/dfbi hold only the second-order absorbing BC term a (μ⁻¹c ∇ₛ×E, ∇ₛ×F), with
  // a = 1/(2iω) = -0.5 (Im ω + i Re ω) / |ω|²: both parts of a are nonpositive at a real ω
  // and for decaying modes (Im ω ≥ 0). The coarse-level solver is applied to the
  // real-valued Br + Bi unless ComplexCoarseSolve is set or it is Jacobi (see
  // UseRealCoarseApproximation), and there this term would cancel part of the positive
  // semidefinite curl-curl term of Br. So flip its sign on the coarse level, like
  // PCMatShifted does for the mass term: it then adds 0.5 (Re ω + Im ω) / |ω|² ≥ 0 times a
  // positive semidefinite form to Br + Bi, as the losses do. The other levels are applied
  // as the complex-valued operator and keep the correct sign, as does the system matrix.
  // Its Floquet terms keep their sign: the rotated boundary mass term is bounded by the
  // first-order absorbing BC term when the tangential Floquet wave vector is below √2 times
  // the wavenumber at the boundary, in particular for propagating waves (above, the
  // boundary condition is not passive, see
  // FarfieldBoundaryOperator::CheckFloquetPassivity), and the mixed boundary curl terms are
  // antisymmetric.
  //
  // Next to a frequency-dependent material k, the term is a_k = 0.5/(iω r_k(ω)), which has
  // no sign which the flip would fix (it is positive and real where Re ε < 0): on the
  // flipped coarse copy only, it is replaced by its magnitude, as in the real
  // preconditioner (SurfaceTerms::AddFlippedPreconditionerBdrCoefficients, which adds a_k
  // and |a_k| to the flipped -a_k). The other levels keep a_k, as the system matrix does.
  MaterialPropertyCoefficient dfbr_flip_add(mat_op.MaxCeedBdrAttribute()),
      dfbi_flip_add(mat_op.MaxCeedBdrAttribute());
  if (pc_coarse_real)
  {
    surface_terms.AddFlippedPreconditionerBdrCoefficients(a3, dfbr_flip_add, dfbi_flip_add);
  }
  if (mat_op.HasFloquetFrequencyScaling())
  {
    // k₀-based tensors: cross terms scale with a1, mass term with a2. For the second-order
    // absorbing BC, the rotated boundary mass term scales with a1 and the mixed boundary
    // curl terms with a0.
    AddImagPeriodicCoefficients(a1.imag(), fpi);
    AddImagPeriodicCoefficients(a1.real(), fpr);
    AddRealPeriodicCoefficients(-a2.real(), fr);
    AddRealPeriodicCoefficients(-a2.imag(), fi);
    AddDampingPeriodicBdrCoefficients(a1.real(), frbr);
    AddDampingPeriodicBdrCoefficients(a1.imag(), frbi);
    AddStiffnessPeriodicBdrCoefficients(a0.real(), fpbr);
    AddStiffnessPeriodicBdrCoefficients(a0.imag(), fpbi);
  }
  else
  {
    AddRealPeriodicCoefficients(a0.real(), fr);
    AddRealPeriodicCoefficients(a0.imag(), fi);
    AddImagPeriodicCoefficients(a0.real(), fpi);
    AddImagPeriodicCoefficients(-a0.imag(), fpr);
  }

  // PML terms mirror the system matrix A = a0 K + a1 C + a2 M + A2(a3): a static stretch
  // with the prefactors of the corresponding bulk terms above, or a frequency-dependent
  // stretch evaluated at ω = a3. The auxiliary space terms are the
  // gradient subspace projections of the mass terms. Note that large PML conductivities
  // (σ / ω) can slow down the multigrid convergence: the imaginary parts of the PML tensors
  // change sign between the directions normal and tangential to the layer, which degrades
  // the effectiveness of the polynomial smoothers.
  std::vector<PMLIntegrator> pml_re, pml_im, aux_pml_re, aux_pml_im;
  AppendPMLPencil(pml_re, &pml_im, &aux_pml_re, &aux_pml_im, mat_op, PMLFilter::STATIC, 0.0,
                  {a0, a1, ShiftedMassCoefficient(a2, pc_mat_shifted), a2});
  AppendPMLPencil(pml_re, &pml_im, &aux_pml_re, &aux_pml_im, mat_op,
                  PMLFilter::FREQUENCY_DEPENDENT, a3,
                  ExtraSystemPMLPencil(a3, pc_mat_shifted));

  // A real-valued coarse solve factors the sum Ar + Ai of the real and imaginary parts of
  // the coarse-level operator of a multigrid hierarchy (see ComplexCoarseSolve), which is a
  // poor approximation of the PML terms: for a stretch factor s ≈ 1 - i σ / (ωε₀), the
  // components Re{t} + Im{t} of the PML tensors of strongly absorbing layers are large and
  // of either sign, so that the coarse-level operator is strongly indefinite and the
  // multigrid preconditioner stalls. The coarse level instead uses the PML tensors with the
  // magnitudes of the stretch factors, those of a real coordinate stretch, which reduce to
  // the bulk material terms for a unit stretch. The finest level, which is the exact
  // preconditioner matrix, keeps the PML terms.
  std::optional<std::array<std::vector<PMLIntegrator>, 2>> coarse_pml;
  if (pc_coarse_real && mat_op.HasPML() && GetNDSpaces().GetNumLevels() > 1)
  {
    auto &[coarse_re, coarse_im] = coarse_pml.emplace();
    AppendPMLPencil(coarse_re, &coarse_im, nullptr, nullptr, mat_op, PMLFilter::STATIC, 0.0,
                    {a0, a1, ShiftedMassCoefficient(a2, pc_mat_shifted), a2}, true);
    AppendPMLPencil(coarse_re, &coarse_im, nullptr, nullptr, mat_op,
                    PMLFilter::FREQUENCY_DEPENDENT, a3,
                    ExtraSystemPMLPencil(a3, pc_mat_shifted), true);
  }

  int empty[2] = {
      (dfr.empty() && fr.empty() && dfbr.empty() && fbr.empty() && fpr.empty() &&
       frbr.empty() && fpbr.empty() && dfbr_flip_add.empty() && pml_re.empty()),
      (dfi.empty() && fi.empty() && dfbi.empty() && fbi.empty() && fpi.empty() &&
       frbi.empty() && fpbi.empty() && dfbi_flip_add.empty() && pml_im.empty())};
  Mpi::GlobalMin(2, empty, GetComm());
  // With an environment, the region's solver factors the complex-valued matrix of the only
  // level (see RegionSolver), so the second-order absorbing BC term keeps its sign there.
  const bool flip_coarse_dfb = pc_coarse_real && !environment;
  if (!empty[0])
  {
    br_vec = AssembleOperators(assembly_domains, GetNDSpaces(), &dfr, &fr, &dfbr, &fbr,
                               &fpr, &frbr, &fpbr, skip_zeros, assemble_q_data,
                               flip_coarse_dfb, &dfbr_flip_add, &pml_re,
                               coarse_pml ? &(*coarse_pml)[0] : nullptr);
    if (!environment)
    {
      br_aux_vec = AssembleAuxOperators(assembly_domains, GetH1Spaces(), &fr, &fbr, &frbr,
                                        skip_zeros, assemble_q_data, &aux_pml_re);
    }
  }
  if (!empty[1])
  {
    bi_vec = AssembleOperators(assembly_domains, GetNDSpaces(), &dfi, &fi, &dfbi, &fbi,
                               &fpi, &frbi, &fpbi, skip_zeros, assemble_q_data,
                               flip_coarse_dfb, &dfbi_flip_add, &pml_im,
                               coarse_pml ? &(*coarse_pml)[1] : nullptr);
    if (!environment)
    {
      bi_aux_vec = AssembleAuxOperators(assembly_domains, GetH1Spaces(), &fi, &fbi, &frbi,
                                        skip_zeros, assemble_q_data, &aux_pml_im);
    }
  }
  if (pml_block)
  {
    // The PML integrators are only assembled on the elements of the PML regions.
    const auto &fespace = GetNDSpace();
    std::unique_ptr<Operator> ar, ai;
    if (!empty[0])
    {
      ar = AssembleOperator(assembly_domains, fespace, nullptr, nullptr, &dfbr, &fbr,
                            nullptr, &frbr, &fpbr, skip_zeros, false, &pml_re);
    }
    if (!empty[1])
    {
      ai = AssembleOperator(assembly_domains, fespace, nullptr, nullptr, &dfbi, &fbi,
                            nullptr, &frbi, &fpbi, skip_zeros, false, &pml_im);
    }
    auto A = std::make_unique<ComplexParOperator>(std::move(ar), std::move(ai), fespace);
    A->SetEssentialTrueDofs(nd_dbc_tdof_lists.back(), Operator::DiagonalPolicy::DIAG_ONE);
    *pml_block = std::move(A);
  }
}

template <typename A3Type>
void SpaceOperator::AssemblePreconditioner(
    std::complex<double> a0, std::complex<double> a1, std::complex<double> a2, A3Type a3,
    std::vector<std::unique_ptr<Operator>> &br_vec,
    std::vector<std::unique_ptr<Operator>> &br_aux_vec)
{
  constexpr bool skip_zeros = false;
  const bool assemble_q_data = !mfem::Device::Allows(mfem::Backend::DEVICE_MASK);
  MaterialPropertyCoefficient dfr(mat_op.MaxCeedAttribute()), fr(mat_op.MaxCeedAttribute()),
      dfbr(mat_op.MaxCeedBdrAttribute()), fbr(mat_op.MaxCeedBdrAttribute()),
      frbr(mat_op.MaxCeedBdrAttribute());
  AddStiffnessCoefficients(a0.real(), dfr, fr);
  AddStiffnessBdrCoefficients(a0.real(), fbr);
  AddDampingCoefficients(a1.imag(), fr);
  AddDampingBdrCoefficients(a1.imag(), fbr);
  const double a2_pc = pc_mat_shifted ? std::abs(a2.real()) : a2.real();
  AddAbsMassCoefficients(a2_pc, fr);
  AddRealMassBdrCoefficients(a2_pc, fbr);
  AddFrequencyDependentPermittivityA2Coefficients(a3, fr);
  AddRealPreconditionerExtraSystemBdrCoefficients(a3, dfbr, fbr, frbr);
  if (mat_op.HasFloquetFrequencyScaling())
  {
    AddRealPeriodicCoefficients(-(pc_mat_shifted ? std::abs(a2.real()) : a2.real()), fr);
    AddDampingPeriodicBdrCoefficients(a1.imag(), frbr);
  }
  else
  {
    AddRealPeriodicCoefficients(a0.real(), fr);
  }

  // Real-valued approximation of the PML terms: Re{c} T̂ with the tensors T̂ of the
  // magnitudes of the stretch factors, which are positive definite like the bulk material
  // properties, and the prefactors of the corresponding bulk terms above.
  std::vector<PMLIntegrator> pml, aux_pml;
  AppendPMLPencil(pml, nullptr, &aux_pml, nullptr, mat_op, PMLFilter::STATIC, 0.0,
                  {a0.real(), 0.0, a2_pc, a2_pc});
  AppendPMLPencil(pml, nullptr, &aux_pml, nullptr, mat_op, PMLFilter::FREQUENCY_DEPENDENT,
                  a3, ExtraSystemPMLPencil(a3, pc_mat_shifted, pc_mat_shifted));

  int empty = (dfr.empty() && fr.empty() && dfbr.empty() && fbr.empty() && frbr.empty() &&
               pml.empty());
  Mpi::GlobalMin(1, &empty, GetComm());
  if (!empty)
  {
    br_vec = AssembleOperators(assembly_domains, GetNDSpaces(), &dfr, &fr, &dfbr, &fbr,
                               nullptr, &frbr, nullptr, skip_zeros, assemble_q_data, false,
                               nullptr, &pml);
    if (!environment)
    {
      br_aux_vec = AssembleAuxOperators(assembly_domains, GetH1Spaces(), &fr, &fbr, &frbr,
                                        skip_zeros, assemble_q_data, &aux_pml);
    }
  }
}

void SpaceOperator::AssemblePreconditioner(
    double a0, double a1, double a2, double a3,
    std::vector<std::unique_ptr<Operator>> &br_vec,
    std::vector<std::unique_ptr<Operator>> &br_aux_vec)
{
  constexpr bool skip_zeros = false;
  const bool assemble_q_data = !mfem::Device::Allows(mfem::Backend::DEVICE_MASK);
  MaterialPropertyCoefficient dfr(mat_op.MaxCeedAttribute()), fr(mat_op.MaxCeedAttribute()),
      dfbr(mat_op.MaxCeedBdrAttribute()), fbr(mat_op.MaxCeedBdrAttribute()),
      frbr(mat_op.MaxCeedBdrAttribute());
  AddStiffnessCoefficients(a0, dfr, fr);
  AddStiffnessBdrCoefficients(a0, fbr);
  AddDampingCoefficients(a1, fr);
  AddDampingBdrCoefficients(a1, fbr);
  const double a2_pc = pc_mat_shifted ? std::abs(a2) : a2;
  AddAbsMassCoefficients(a2_pc, fr);
  AddRealMassBdrCoefficients(a2_pc, fbr);
  AddFrequencyDependentPermittivityA2Coefficients(a3, fr);
  AddRealPreconditionerExtraSystemBdrCoefficients(a3, dfbr, fbr, frbr);
  if (mat_op.HasFloquetFrequencyScaling())
  {
    AddRealPeriodicCoefficients(-(pc_mat_shifted ? std::abs(a2) : a2), fr);
    AddDampingPeriodicBdrCoefficients(a1, frbr);
  }
  else
  {
    AddRealPeriodicCoefficients(a0, fr);
  }

  CheckPMLOperatorType<Operator>(mat_op);
  int empty = (dfr.empty() && fr.empty() && dfbr.empty() && fbr.empty() && frbr.empty());
  Mpi::GlobalMin(1, &empty, GetComm());
  if (!empty)
  {
    br_vec = AssembleOperators(assembly_domains, GetNDSpaces(), &dfr, &fr, &dfbr, &fbr,
                               nullptr, &frbr, nullptr, skip_zeros, assemble_q_data);
    if (!environment)
    {
      br_aux_vec = AssembleAuxOperators(assembly_domains, GetH1Spaces(), &fr, &fbr, &frbr,
                                        skip_zeros, assemble_q_data);
    }
  }
}

namespace
{

// Largest magnitude max_a |s_a| of the PML stretch factors at the vertices of each local
// element of the mesh, or zero for the elements outside of the PML regions. A
// frequency-dependent stretch is evaluated at omega. For an extruded PML, the stretch
// factors are the principal ones along and across the direction of the extrusion (see
// pml::EvaluateStretch), at the vertices of the curved elements.
std::vector<double> GetPMLStretchMagnitudes(const Mesh &mesh,
                                            const MaterialOperator &mat_op,
                                            std::complex<double> omega)
{
  const auto &pmesh = mesh.Get();
  std::vector<double> m(pmesh.GetNE(), 0.0);
  mfem::Array<int> verts;
  for (int e = 0; e < pmesh.GetNE(); e++)
  {
    const int attr = pmesh.GetAttribute(e);
    for (const auto &layer : mat_op.GetPMLLayers())
    {
      const auto &attrs = layer.GetAttributes();
      if (!std::binary_search(attrs.begin(), attrs.end(), attr))
      {
        continue;
      }
      m[e] = 1.0;
      if (layer.GetStretch().extruded)
      {
        const mfem::IntegrationRule &ir =
            *mfem::Geometries.GetVertices(pmesh.GetElementGeometry(e));
        mfem::IsoparametricTransformation T;
        pmesh.GetElementTransformation(e, &T);
        mfem::Vector x(3);
        for (int j = 0; j < ir.GetNPoints(); j++)
        {
          const mfem::IntegrationPoint &ip = ir.IntPoint(j);
          T.SetIntPoint(&ip);
          T.Transform(ip, x);
          for (const auto &s :
               pml::EvaluateStretch(layer.GetStretch(), {x(0), x(1), x(2)}, omega))
          {
            m[e] = std::max(m[e], std::abs(s));
          }
        }
        break;
      }
      pmesh.GetElementVertices(e, verts);
      for (const int v : verts)
      {
        const double *x = pmesh.GetVertex(v);
        for (const auto &s :
             pml::EvaluateStretch(layer.GetStretch(), {x[0], x[1], x[2]}, omega))
        {
          m[e] = std::max(m[e], std::abs(s));
        }
      }
      break;
    }
  }
  return m;
}

// Local true DOFs of the elements of the PML regions whose stretch magnitude (see
// GetPMLStretchMagnitudes) exceeds the threshold (all of the PML elements for a threshold
// of at most 1), which are not shared with other elements, sorted and unique. The rows and
// columns of the system matrix for these DOFs only have contributions from these elements
// and the boundary elements. The DOFs of the local elements are mapped to true DOFs on all
// processes, since a true DOF of an element may be owned by a process which has no element
// containing it (and, for nonconforming meshes, the DOFs of an element may be constrained
// by true DOFs of other elements).
std::vector<int> GetPMLTrueDofs(const FiniteElementSpace &fespace,
                                const std::vector<double> &stretch_magnitude,
                                double threshold)
{
  const auto &pfes = fespace.Get();
  const auto &mesh = fespace.GetMesh();
  MFEM_VERIFY(stretch_magnitude.size() == static_cast<std::size_t>(mesh.Get().GetNE()),
              "Invalid PML stretch magnitudes for the mesh elements!");
  auto MarkTrueDofs = [&](bool inside)
  {
    mfem::Array<int> ldof_marker(pfes.GetVSize()), marker(pfes.GetTrueVSize());
    ldof_marker = 0;
    mfem::Array<int> vdofs;
    for (int e = 0; e < mesh.Get().GetNE(); e++)
    {
      const double m = stretch_magnitude[e];
      if ((m > 0.0 && (threshold <= 1.0 || m > threshold)) != inside)
      {
        continue;
      }
      pfes.GetElementVDofs(e, vdofs);
      for (const int vdof : vdofs)
      {
        ldof_marker[(vdof >= 0) ? vdof : -1 - vdof] = 1;
      }
    }
    pfes.Dof_TrueDof_Matrix()->BooleanMultTranspose(1, ldof_marker.HostRead(), 0,
                                                    marker.HostWrite());
    return marker;
  };
  const auto in_marker = MarkTrueDofs(true), out_marker = MarkTrueDofs(false);
  std::vector<int> tdofs;
  for (int i = 0; i < in_marker.Size(); i++)
  {
    if (in_marker[i] && !out_marker[i])
    {
      tdofs.push_back(i);
    }
  }
  return tdofs;
}

}  // namespace

template <typename OperType, typename ScalarType, typename A3Type>
std::unique_ptr<OperType> SpaceOperator::GetPreconditionerMatrix(ScalarType a0,
                                                                 ScalarType a1,
                                                                 ScalarType a2, A3Type a3)
{
  // When partially assembled, the coarse operators can reuse the fine operator quadrature
  // data if the spaces correspond to the same mesh. When appropriate, we build the
  // preconditioner on all levels based on the actual complex-valued system matrix. The
  // coarse operator is always fully assembled.
  //
  // The driven solvers may apply the system operator through the finest level alone (see
  // GetSystemOperator), without assembling K, C or M, so print the header here too.
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);
  MFEM_VERIFY(GetH1Spaces().GetNumLevels() == GetNDSpaces().GetNumLevels(),
              "Multigrid hierarchy mismatch for auxiliary space preconditioning!");

  const auto n_levels = GetNDSpaces().GetNumLevels();
  std::vector<std::unique_ptr<Operator>> br_vec(n_levels), bi_vec(n_levels),
      br_aux_vec(n_levels), bi_aux_vec(n_levels);
  std::unique_ptr<ComplexOperator> pml_block;
  if (std::is_same_v<OperType, ComplexOperator> && !pc_mat_real)
  {
    AssemblePreconditioner(a0, a1, a2, a3, br_vec, br_aux_vec, bi_vec, bi_aux_vec,
                           (pml_subdomain && mat_op.HasPML()) ? &pml_block : nullptr);
  }
  else
  {
    AssemblePreconditioner(a0, a1, a2, a3, br_vec, br_aux_vec);
  }

  // Print after the assembly, which may first solve and print the wave port modes at a new
  // frequency.
  if (print_prec_hdr)
  {
    Mpi::Print("\nAssembling multigrid hierarchy:\n");
  }
  auto B = std::make_unique<BaseMultigridOperator<OperType>>(n_levels);
  for (bool aux : {false, true})
  {
    if (aux && environment)
    {
      continue;  // the RegionSolver has no auxiliary space
    }
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
      if constexpr (std::is_same_v<OperType, ComplexOperator> &&
                    std::is_same_v<A3Type, double>)
      {
        if (environment)
        {
          // The region's single level leaves S_E(ω) out: the RegionSolver adds it to the
          // region's Schur complement onto Γ at the frequency of the matrix.
          environment->SetPreconditionerFrequency(a3);
        }
      }
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

  if constexpr (std::is_same_v<OperType, ComplexOperator>)
  {
    if (pml_block)
    {
      // The PML regions are the subdomain of the subdomain correction of the multigrid
      // preconditioner (see PMLSubdomainSolver).
      SelectPMLTrueDofs(a3);
      B->SetSubdomain(std::vector<int>(pml_subdomain_tdofs), std::move(pml_block));
      if (!pml_smoother_masks.empty())
      {
        B->SetSmootherMasks(std::vector<std::vector<int>>(pml_smoother_masks),
                            std::vector<std::vector<int>>(pml_aux_smoother_masks));
      }
    }
  }

  print_prec_hdr = false;
  return B;
}

void SpaceOperator::SetEnvironment(EnvironmentOperator &env)
{
  MFEM_VERIFY(!environment && GetNDSpaces().GetNumLevels() == 1,
              "The region of a driven substructuring is set up once, on one level!");
  environment = &env;
  region_restriction.emplace(assembly_domains, env.RegionAttributes());
  auto &dbc = nd_dbc_tdof_lists.back();
  const auto &pinned = env.PinnedTrueDofs();
  std::vector<int> merged(dbc.HostRead(), dbc.HostRead() + dbc.Size());
  merged.insert(merged.end(), pinned.HostRead(), pinned.HostRead() + pinned.Size());
  std::ranges::sort(merged);
  merged.erase(std::unique(merged.begin(), merged.end()), merged.end());
  dbc.SetSize(static_cast<int>(merged.size()));
  std::copy(merged.begin(), merged.end(), dbc.HostWrite());
}

std::unique_ptr<ComplexOperator> SpaceOperator::GetSystemOperator(double omega,
                                                                  const ComplexOperator &P)
{
  if (pc_mat_real || pc_mat_shifted || GetNDSpaces().GetNumLevels() < 2)
  {
    return {};
  }
  const auto *mg = dynamic_cast<const ComplexMultigridOperator *>(&P);
  MFEM_VERIFY(mg, "GetSystemOperator requires the preconditioner from "
                  "GetPreconditionerMatrix!");

  // The finest level has the same terms as K, C, M and the sparse part of
  // GetExtraSystemOperator(ω). Add the matrix-free terms of GetExtraSystemOperator.
  auto A = std::make_unique<SumComplexOperator>(mg->GetFinestOperator());
  if (auto W = wave_port_op.GetModalCorrectionOperator(omega, GetNDSpace(),
                                                       nd_dbc_tdof_lists.back()))
  {
    A->AddOperator(std::move(W));
  }
  if (auto F = floquet_port_op.GetExtraSystemOperator(omega))
  {
    A->AddOperator(std::move(F));
  }
  return A;
}

void SpaceOperator::RestrictToAssemblyDomains(
    std::initializer_list<MaterialPropertyCoefficient *> f,
    std::initializer_list<MaterialPropertyCoefficient *> fb)
{
  const auto *attrs = assembly_domains.Get();
  if (!attrs)
  {
    return;
  }
  const auto loc_attr = GetMesh().GetCeedAttributes(*attrs);
  const auto loc_bdr_attr = GetMesh().GetCeedBdrAttributesByNeighbor(*attrs);
  for (auto *c : f)
  {
    c->RestrictCoefficient(loc_attr);
  }
  for (auto *c : fb)
  {
    c->RestrictCoefficient(loc_bdr_attr);
  }
}

void SpaceOperator::SelectPMLTrueDofs(std::complex<double> omega)
{
  // The multigrid smoothers do not converge on the PML equations of strongly stretched
  // elements: the Jacobi-scaled operator has eigenvalues near the imaginary axis and in the
  // left half-plane, which the Chebyshev polynomials amplify. The DOFs of these elements
  // are instead only updated by the coarse-level correction and the PML subdomain
  // correction, whose subdomain (all of the PML elements) contains the smoother masks.
  if (pml_dofs_selected)
  {
    return;
  }
  pml_dofs_selected = true;
  const auto &nd_fespaces = GetNDSpaces();
  const auto n_levels = nd_fespaces.GetNumLevels();
  const auto &fespace = nd_fespaces.GetFinestFESpace();
  const auto m = GetPMLStretchMagnitudes(fespace.GetMesh(), mat_op, omega);
  pml_subdomain_tdofs = GetPMLTrueDofs(fespace, m, 0.0);
  pml_smoother_masks.clear();
  pml_aux_smoother_masks.clear();
  if (pml_smoother_threshold > 0.0)
  {
    pml_smoother_masks.resize(n_levels);
    pml_aux_smoother_masks.resize(n_levels);
    for (std::size_t l = 1; l < n_levels; l++)
    {
      const auto &nd_fespace = nd_fespaces.GetFESpaceAtLevel(l);
      const auto m_l = (&nd_fespace.GetMesh() == &fespace.GetMesh())
                           ? m
                           : GetPMLStretchMagnitudes(nd_fespace.GetMesh(), mat_op, omega);
      pml_smoother_masks[l] = GetPMLTrueDofs(nd_fespace, m_l, pml_smoother_threshold);
      pml_aux_smoother_masks[l] =
          GetPMLTrueDofs(GetH1Spaces().GetFESpaceAtLevel(l), m_l, pml_smoother_threshold);
    }
  }
  HYPRE_BigInt n[3] = {
      static_cast<HYPRE_BigInt>(pml_subdomain_tdofs.size()),
      pml_smoother_masks.empty()
          ? 0
          : static_cast<HYPRE_BigInt>(pml_smoother_masks.back().size()),
      pml_aux_smoother_masks.empty()
          ? 0
          : static_cast<HYPRE_BigInt>(pml_aux_smoother_masks.back().size())};
  Mpi::GlobalSum(3, n, GetComm());
  if (pml_smoother_threshold > 0.0)
  {
    Mpi::Print(" PML elements with stretch magnitude max|s| > {:.2f}: {:d} + {:d} "
               "(auxiliary) unknowns are not updated by the smoothers\n",
               pml_smoother_threshold, n[1], n[2]);
  }
}

void SpaceOperator::AddStiffnessCoefficients(double coeff, MaterialPropertyCoefficient &df,
                                             MaterialPropertyCoefficient &f)
{
  // Contribution from material permeability. In 2D, curl is scalar so the curl-curl
  // coefficient is scalar (1x1). The bulk material coefficients exclude the PML attributes,
  // whose contributions come from the PML integrators (see AppendPML).
  df.AddCoefficient(mat_op.GetBulkAttributeToMaterial(),
                    mat_op.GetCurlCurlInvPermeability(), coeff);

  // Contribution for London superconductors.
  if (mat_op.HasLondonDepth())
  {
    f.AddCoefficient(mat_op.GetBulkAttributeToMaterial(), mat_op.GetInvLondonDepth(),
                     coeff);
  }
}

void SpaceOperator::AddStiffnessBdrCoefficients(double coeff,
                                                MaterialPropertyCoefficient &fb)
{
  // Robin BC contributions due to surface impedance and lumped ports (inductance).
  surface_terms.AddStiffnessBdrCoefficients(coeff, fb);
  lumped_port_op.AddStiffnessBdrCoefficients(coeff, fb);
}

void SpaceOperator::AddDampingCoefficients(double coeff, MaterialPropertyCoefficient &f)
{
  // Contribution for domain conductivity.
  if (mat_op.HasConductivity())
  {
    f.AddCoefficient(mat_op.GetBulkAttributeToMaterial(), mat_op.GetConductivity(), coeff);
  }
}

void SpaceOperator::AddDampingBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb)
{
  // Robin BC contributions due to surface impedance, lumped ports, and absorbing
  // boundaries (resistance).
  surface_terms.AddDampingBdrCoefficients(coeff, fb);
  lumped_port_op.AddDampingBdrCoefficients(coeff, fb);
}

void SpaceOperator::AddRealMassCoefficients(double coeff, MaterialPropertyCoefficient &f)
{
  // The bulk material coefficients exclude the PML attributes (see AppendPML).
  f.AddCoefficient(mat_op.GetBulkAttributeToMaterial(), mat_op.GetPermittivityReal(),
                   coeff);
}

void SpaceOperator::AddRealMassBdrCoefficients(double coeff,
                                               MaterialPropertyCoefficient &fb)
{
  // Robin BC contributions due to surface impedance and lumped ports (capacitance).
  surface_terms.AddMassBdrCoefficients(coeff, fb);
  lumped_port_op.AddMassBdrCoefficients(coeff, fb);
}

void SpaceOperator::AddImagMassCoefficients(double coeff, MaterialPropertyCoefficient &f)
{
  // Contribution for loss tangent: ε -> ε * (1 - i tan(δ)).
  if (mat_op.HasLossTangent())
  {
    f.AddCoefficient(mat_op.GetBulkAttributeToMaterial(), mat_op.GetPermittivityImag(),
                     coeff);
  }
}

void SpaceOperator::AddAbsMassCoefficients(double coeff, MaterialPropertyCoefficient &f)
{
  f.AddCoefficient(mat_op.GetBulkAttributeToMaterial(), mat_op.GetPermittivityAbs(), coeff);
}

void SpaceOperator::AddExtraSystemBdrCoefficients(
    double omega, MaterialPropertyCoefficient &dfbr, MaterialPropertyCoefficient &dfbi,
    MaterialPropertyCoefficient &fbr, MaterialPropertyCoefficient &fbi,
    MaterialPropertyCoefficient &frbr, MaterialPropertyCoefficient &frbi,
    MaterialPropertyCoefficient *fpbr, MaterialPropertyCoefficient *fpbi,
    bool include_wave_ports)
{
  // Contribution for second-order farfield boundaries, finite conductivity boundaries,
  // rational impedance boundaries, and absorbing boundaries next to frequency-dependent
  // materials.
  surface_terms.AddExtraSystemBdrCoefficients(omega, dfbr, dfbi, fbr, fbi, frbr, frbi, fpbr,
                                              fpbi);

  // Contribution for numeric wave ports. Skipped when the caller will apply the wave-port
  // contribution separately via per-port operators (see GetWavePortBoundaryMassMatrix).
  if (include_wave_ports)
  {
    wave_port_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);
  }

  // Contribution for Floquet ports (Robin BC).
  floquet_port_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);
}

void SpaceOperator::AddExtraSystemBdrCoefficients(
    std::complex<double> omega, MaterialPropertyCoefficient &dfbr,
    MaterialPropertyCoefficient &dfbi, MaterialPropertyCoefficient &fbr,
    MaterialPropertyCoefficient &fbi, MaterialPropertyCoefficient &frbr,
    MaterialPropertyCoefficient &frbi, MaterialPropertyCoefficient *fpbr,
    MaterialPropertyCoefficient *fpbi)
{
  surface_terms.AddExtraSystemBdrCoefficients(omega, dfbr, dfbi, fbr, fbi, frbr, frbi, fpbr,
                                              fpbi);
  wave_port_op.AddExtraSystemBdrCoefficients(omega, fbr, fbi);
}

void SpaceOperator::AddExtraSystemBdrCoefficients(
    std::complex<double> omega,
    const std::vector<std::pair<int, std::complex<double>>> &port_kn,
    MaterialPropertyCoefficient &dfbr, MaterialPropertyCoefficient &dfbi,
    MaterialPropertyCoefficient &fbr, MaterialPropertyCoefficient &fbi,
    MaterialPropertyCoefficient &frbr, MaterialPropertyCoefficient &frbi,
    MaterialPropertyCoefficient *fpbr, MaterialPropertyCoefficient *fpbi)
{
  // Complex-ω overload for the eigenmode nonlinear solve and its matching preconditioner:
  // all frequency-dependent boundary terms (2nd-order farfield ABC, surface conductivity,
  // rational surface impedance, numeric wave ports) are evaluated at the genuinely complex
  // frequency (ω = -i·λ) so the assembled A2(λ) is the exact analytic continuation. For
  // real ω the complex overloads reduce to their double overloads.
  surface_terms.AddExtraSystemBdrCoefficients(omega, dfbr, dfbi, fbr, fbi, frbr, frbi, fpbr,
                                              fpbi);
  wave_port_op.AddExtraSystemBdrCoefficients(port_kn, fbr, fbi);
}

void SpaceOperator::AddRealPreconditionerExtraSystemBdrCoefficients(
    double omega, MaterialPropertyCoefficient &dfb, MaterialPropertyCoefficient &fb,
    MaterialPropertyCoefficient &frb)
{
  surface_terms.AddRealPreconditionerExtraSystemBdrCoefficients(omega, dfb, fb, frb);
  wave_port_op.AddExtraSystemBdrCoefficients(omega, fb, fb);
  floquet_port_op.AddExtraSystemBdrCoefficients(omega, fb, fb);
}

void SpaceOperator::AddRealPreconditionerExtraSystemBdrCoefficients(
    std::complex<double> omega, MaterialPropertyCoefficient &dfb,
    MaterialPropertyCoefficient &fb, MaterialPropertyCoefficient &frb)
{
  surface_terms.AddRealPreconditionerExtraSystemBdrCoefficients(omega, dfb, fb, frb);
  wave_port_op.AddExtraSystemBdrCoefficients(omega, fb, fb);
}

bool SpaceOperator::HasExtraSystemBdrTerms() const
{
  // See AddExtraSystemBdrCoefficients for a complex ω: the second-order absorbing
  // boundaries, the frequency-dependent boundary mass terms of SurfaceTerms, and the wave
  // ports.
  const auto &farfield_op = surface_terms.GetFarfieldOp();
  return (farfield_op.GetOrder() > 1 && farfield_op.GetAttrList().Size() > 0) ||
         surface_terms.GetNumExtraSystemBdrMassTerms() > 0 || wave_port_op.Size() > 0;
}

void SpaceOperator::AddRealPeriodicCoefficients(double coeff,
                                                MaterialPropertyCoefficient &f)
{
  // Floquet periodicity contributions.
  if (mat_op.HasWaveVector())
  {
    f.AddCoefficient(mat_op.GetBulkAttributeToMaterial(), mat_op.GetFloquetMass(), coeff);
  }
}

void SpaceOperator::AddImagPeriodicCoefficients(double coeff,
                                                MaterialPropertyCoefficient &f)
{
  // Floquet periodicity contributions.
  if (mat_op.HasWaveVector())
  {
    f.AddCoefficient(mat_op.GetBulkAttributeToMaterial(), mat_op.GetFloquetCurl(), coeff);
  }
}

void SpaceOperator::AddDampingPeriodicBdrCoefficients(double coeff,
                                                      MaterialPropertyCoefficient &frb)
{
  // Floquet periodicity contributions for second-order farfield boundaries, with a
  // frequency dependent Floquet wave vector.
  surface_terms.GetFarfieldOp().AddDampingFloquetBdrCoefficients(coeff, frb);
}

void SpaceOperator::AddStiffnessPeriodicBdrCoefficients(double coeff,
                                                        MaterialPropertyCoefficient &fpb)
{
  // Floquet periodicity contributions for second-order farfield boundaries, with a
  // frequency dependent Floquet wave vector.
  surface_terms.GetFarfieldOp().AddStiffnessFloquetBdrCoefficients(coeff, fpb);
}

bool SpaceOperator::GetExcitationVector(int excitation_idx, Vector &RHS)
{
  // Time domain excitation vector.
  RHS.SetSize(GetNDSpace().GetTrueVSize());
  RHS.UseDevice(true);
  RHS = 0.0;
  bool nnz =
      AddExcitationVector1Internal(excitation_idx, RHS, nullptr, /*time_harmonic*/ false);
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
  bool nnz1 = AddExcitationVector1Internal(excitation_idx, RHS.Real(), &RHS.Imag(),
                                           /*time_harmonic*/ true);
  RHS *= 1i * omega;
  bool nnz2 = AddExcitationVector2Internal(excitation_idx, omega, RHS);
  linalg::SetSubVector(RHS, nd_dbc_tdof_lists.back(), 0.0);
  // The sources of the environment, condensed on Γ (their own entries are eliminated).
  const bool nnz_env =
      environment && environment->AddExcitation(excitation_idx, omega, RHS);
  return nnz1 || nnz2 || nnz_env;
}

void SpaceOperator::GetLumpedPortExcitationVectorPrimaryEt(int port_idx,
                                                           ComplexVector &Et_primary)
{
  const auto &data = GetLumpedPortOp().GetPort(port_idx);

  SumVectorCoefficient fb(GetMesh().SpaceDimension());
  for (const auto &elem : data.elems)
  {
    const double Rs = 1.0 * data.GetToSquare(*elem);
    const double Einc = std::sqrt(
        Rs / (elem->GetGeometryWidth() * elem->GetGeometryLength() * data.elems.size()));
    fb.AddCoefficient(elem->GetModeCoefficient(Einc));
  }

  Et_primary.SetSize(GetNDSpace().GetTrueVSize());
  Et_primary.UseDevice(true);
  Et_primary = 0.0;

  // L2-project the port mode onto the ND tangent-trace space; DOF interpolation
  // gives a different vector and breaks the W-norm = 1 invariant Einc relies on.
  ProjectBdrCoefficientViaMassSolve(fb, data, mat_op, GetNDSpace(), GetComm(),
                                    Et_primary.Real());

  linalg::SetSubVector(Et_primary.Real(), GetNDDbcTDofLists().back(), 0.0);
}

void SpaceOperator::GetLumpedPortExcitationVectorPrimaryHtcn(int port_idx,
                                                             ComplexVector &Htcn_primary)
{
  const auto &data = lumped_port_op.GetPort(port_idx);

  SumVectorCoefficient fb(GetMesh().SpaceDimension());
  for (const auto &elem : data.elems)
  {
    const double Rs = 1.0 * data.GetToSquare(*elem);
    const double Hinc = 1.0 / std::sqrt(Rs * elem->GetGeometryWidth() *
                                        elem->GetGeometryLength() * data.elems.size());
    fb.AddCoefficient(elem->GetModeCoefficient(Hinc));
  }

  Htcn_primary.SetSize(GetNDSpace().GetTrueVSize());
  Htcn_primary.UseDevice(true);
  Htcn_primary = 0.0;

  // L2-project the port mode onto the ND tangent-trace space; DOF interpolation
  // gives a different vector and breaks the W-norm = 1 invariant Einc relies on.
  ProjectBdrCoefficientViaMassSolve(fb, data, mat_op, GetNDSpace(), GetComm(),
                                    Htcn_primary.Real());

  linalg::SetSubVector(Htcn_primary.Real(), GetNDDbcTDofLists().back(), 0.0);
}

void SpaceOperator::GetWavePortFieldVectorPrimaryEt(int port_idx, double omega_ref,
                                                    ComplexVector &Et_primary)
{
  // Trigger (or refresh) the cross-section EVP at omega_ref. WavePortData::Initialize
  // caches by omega so calling this is cheap if already initialised.
  wave_port_op.GetWavePortKn(port_idx, omega_ref);
  const auto &data = wave_port_op.GetPort(port_idx);

  // Project the modal tangential E-field onto the parent ND space, restricted to the
  // port boundary attributes. Real and imaginary parts are projected separately, since
  // the full waveport mode is generally complex.
  Et_primary.SetSize(GetNDSpace().GetTrueVSize());
  Et_primary.UseDevice(true);
  Et_primary = 0.0;

  const auto &mesh = GetNDSpace().GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker;
  mesh::AttrToMarker(bdr_attr_max, data.GetAttrList(), attr_marker);

  GridFunction rhs_re(GetNDSpace());
  GridFunction rhs_im(GetNDSpace());
  rhs_re = 0.0;
  rhs_im = 0.0;
  {
    auto fb_re = data.GetModeFieldCoefficientReal(1.0);
    rhs_re.Real().ProjectBdrCoefficientTangent(*fb_re, attr_marker);
  }
  {
    auto fb_im = data.GetModeFieldCoefficientImag(1.0);
    rhs_im.Real().ProjectBdrCoefficientTangent(*fb_im, attr_marker);
  }
  GetNDSpace().GetRestrictionMatrix()->Mult(rhs_re.Real(), Et_primary.Real());
  GetNDSpace().GetRestrictionMatrix()->Mult(rhs_im.Real(), Et_primary.Imag());

  linalg::SetSubVector(Et_primary.Real(), GetNDDbcTDofLists().back(), 0.0);
  linalg::SetSubVector(Et_primary.Imag(), GetNDDbcTDofLists().back(), 0.0);

  // The phase of the port mode is fixed by its S-parameter normalization (a real
  // unconjugated reaction), which leaves a lossy mode of uniform phase, such as a TEM mode
  // of a lossy line, with a small nonzero phase φ. Remove it, with eᵀe = |e|² e^{2iφ} for
  // e = |e| e^{iφ}, so that such a mode adds a single real vector to the synthesis basis
  // rather than two parallel ones (the basis only depends on the span of the mode field).
  const double rr = linalg::Dot(GetComm(), Et_primary.Real(), Et_primary.Real());
  const double ii = linalg::Dot(GetComm(), Et_primary.Imag(), Et_primary.Imag());
  const double ri = linalg::Dot(GetComm(), Et_primary.Real(), Et_primary.Imag());
  const std::complex<double> ee(rr - ii, 2.0 * ri);
  if (std::abs(ee) > 0.0)
  {
    const std::complex<double> rot = std::polar(1.0, -0.5 * std::arg(ee));
    ComplexVector::AXPBY(rot, Et_primary.Real(), Et_primary.Imag(), 0.0, Et_primary.Real(),
                         Et_primary.Imag());
  }
}

bool SpaceOperator::GetExcitationVector1(int excitation_idx, ComplexVector &RHS1)
{
  // Assemble the frequency domain excitation term with linear frequency dependence
  // (coefficient iω, see GetExcitationVector above, is accounted for later). RHS1 is
  // complex when a lumped port has a complex excitation weight.
  RHS1.SetSize(GetNDSpace().GetTrueVSize());
  RHS1.UseDevice(true);
  RHS1 = 0.0;
  bool nnz1 = AddExcitationVector1Internal(excitation_idx, RHS1.Real(), &RHS1.Imag(),
                                           /*time_harmonic*/ true);
  linalg::SetSubVector(RHS1, nd_dbc_tdof_lists.back(), 0.0);
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
  const bool nnz_env =
      environment && environment->AddExcitation(excitation_idx, omega, RHS2);
  return nnz2 || nnz_env;
}

bool SpaceOperator::AddExcitationVector1Internal(int excitation_idx, Vector &RHS1r,
                                                 Vector *RHS1i, bool time_harmonic)
{
  // Assemble the time domain excitation -g'(t) J or frequency domain excitation -iω J.
  // The g'(t) or iω factors are not accounted for here, they are accounted for in the time
  // integration or frequency sweep later. The lumped port incident fields include the
  // complex excitation weight of each port: the real part is added to RHS1r, and the
  // imaginary part (nonzero only for a complex excitation weight) is added to RHS1i. The
  // time domain excitation is real, so RHS1i is null in that case. Surface currents and
  // current dipoles are real.
  MFEM_VERIFY(RHS1r.Size() == GetNDSpace().GetTrueVSize(),
              "Invalid T-vector size for AddExcitationVector1Internal!");
  MFEM_VERIFY(RHS1i == nullptr || RHS1i->Size() == GetNDSpace().GetTrueVSize(),
              "Invalid T-vector size for AddExcitationVector1Internal!");

  // Boundary sources
  SumVectorCoefficient fb(GetMesh().SpaceDimension()), fbi(GetMesh().SpaceDimension());
  lumped_port_op.AddExcitationBdrCoefficients(excitation_idx, fb, fbi, time_harmonic);
  surf_j_op.AddExcitationBdrCoefficients(fb);  // No excitation_idx — currently in all

  // Domain sources (current dipoles) - use integrator-based approach
  bool has_current_dipoles = !current_dipole_op.Empty();

  // The emptiness flags must agree on all ranks, since the assembly and the multiplication
  // by the prolongation matrix below are collective.
  int empty[3] = {(fb.empty()), (!has_current_dipoles), (fbi.empty())};
  Mpi::GlobalMin(3, empty, GetComm());
  if (empty[0] && empty[1] && empty[2])
  {
    return false;
  }

  // A complex excitation is not defined for a real (time domain) vector. This is checked
  // for the configuration, so it should not be reached.
  MFEM_VERIFY(empty[2] || RHS1i != nullptr,
              "A complex excitation weight (ExcitationPhase not a multiple of 180 "
              "degrees) is not supported for real (time domain) excitations!");

  if (!empty[0] || !empty[1])
  {
    mfem::ParLinearForm rhs1(&GetNDSpace().Get());

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
    GetNDSpace().GetProlongationMatrix()->AddMultTranspose(rhs1, RHS1r);
  }

  // Imaginary part of the lumped port boundary sources.
  if (!empty[2])
  {
    mfem::ParLinearForm rhs1i(&GetNDSpace().Get());
    rhs1i.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbi));
    rhs1i.UseFastAssembly(false);
    rhs1i.UseDevice(false);
    rhs1i.Assemble();
    rhs1i.UseDevice(true);
    GetNDSpace().GetProlongationMatrix()->AddMultTranspose(rhs1i, *RHS1i);
  }
  return true;
}

bool SpaceOperator::AddExcitationVector2Internal(int excitation_idx, double omega,
                                                 ComplexVector &RHS2)
{
  // Assemble the contribution of wave ports to the frequency domain excitation term at
  // the specified frequency.
  MFEM_VERIFY(RHS2.Size() == GetNDSpace().GetTrueVSize(),
              "Invalid T-vector size for AddExcitationVector2Internal!");
  // Floquet port excitation: directly adds to RHS (not via boundary linear form).
  bool nnz_floquet = floquet_port_op.AddExcitationVector(excitation_idx, omega, RHS2);

  SumVectorCoefficient fbr(GetMesh().SpaceDimension()), fbi(GetMesh().SpaceDimension());
  wave_port_op.AddExcitationBdrCoefficients(excitation_idx, omega, fbr, fbi);
  int empty = (fbr.empty() && fbi.empty());
  Mpi::GlobalMin(1, &empty, GetComm());
  if (empty)
  {
    return nnz_floquet;
  }
  // Restrict the boundary assembly to the excited wave port boundary elements (the
  // coefficients vanish elsewhere).
  mfem::Array<int> marker = wave_port_op.GetExcitationBdrMarker(excitation_idx);
  MFEM_VERIFY(marker.Size() > 0,
              "No excited wave port found for excitation index " << excitation_idx << "!");
  {
    mfem::LinearForm rhs2(&GetNDSpace().Get());
    rhs2.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbr), marker);
    rhs2.UseFastAssembly(false);
    rhs2.UseDevice(false);
    rhs2.Assemble();
    rhs2.UseDevice(true);
    GetNDSpace().GetProlongationMatrix()->AddMultTranspose(rhs2, RHS2.Real());
  }
  {
    mfem::LinearForm rhs2(&GetNDSpace().Get());
    rhs2.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbi), marker);
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

void SpaceOperator::GetRandomInitialVector(ComplexVector &v, int seed)
{
  v.SetSize(GetNDSpace().GetTrueVSize());
  v.UseDevice(true);
  linalg::SetRandom(GetNDSpace().GetComm(), v, seed);
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
SpaceOperator::GetExtraSystemMatrix(double, Operator::DiagonalPolicy, bool, bool);
template std::unique_ptr<ComplexOperator>
SpaceOperator::GetExtraSystemMatrix(double, Operator::DiagonalPolicy, bool, bool);

template std::unique_ptr<Operator>
SpaceOperator::GetWavePortBoundaryMassMatrix(int, Operator::DiagonalPolicy);
template std::unique_ptr<ComplexOperator>
SpaceOperator::GetWavePortBoundaryMassMatrix(int, Operator::DiagonalPolicy);

template std::unique_ptr<Operator>
SpaceOperator::GetFarfieldBoundaryCurlCurlMatrix(Operator::DiagonalPolicy, bool);
template std::unique_ptr<ComplexOperator>
SpaceOperator::GetFarfieldBoundaryCurlCurlMatrix(Operator::DiagonalPolicy, bool);

template std::unique_ptr<Operator>
SpaceOperator::GetExtraSystemBdrMassTermMatrix(int, Operator::DiagonalPolicy, bool);
template std::unique_ptr<ComplexOperator>
SpaceOperator::GetExtraSystemBdrMassTermMatrix(int, Operator::DiagonalPolicy, bool);

template std::unique_ptr<Operator>
SpaceOperator::GetFloquetRobinBoundaryMassMatrix(int, Operator::DiagonalPolicy);
template std::unique_ptr<ComplexOperator>
SpaceOperator::GetFloquetRobinBoundaryMassMatrix(int, Operator::DiagonalPolicy);

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
SpaceOperator::GetPreconditionerMatrix<Operator, double, double>(double, double, double,
                                                                 double);
template std::unique_ptr<ComplexOperator>
SpaceOperator::GetPreconditionerMatrix<ComplexOperator, std::complex<double>, double>(
    std::complex<double>, std::complex<double>, std::complex<double>, double);
template std::unique_ptr<ComplexOperator>
    SpaceOperator::GetPreconditionerMatrix<ComplexOperator, std::complex<double>,
                                           std::complex<double>>(std::complex<double>,
                                                                 std::complex<double>,
                                                                 std::complex<double>,
                                                                 std::complex<double>);

}  // namespace palace
