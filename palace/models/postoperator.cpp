// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "postoperator.hpp"

#include <algorithm>
#include <string>
#include "fem/coefficient.hpp"
#include "fem/errorindicator.hpp"
#include "models/curlcurloperator.hpp"
#include "models/laplaceoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfacecurrentoperator.hpp"
#include "models/waveportoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

std::string ParaviewFoldername(const ProblemType solver_t)
{
  switch (solver_t)
  {
    case ProblemType::DRIVEN:
      return "driven";
    case ProblemType::EIGENMODE:
      return "eigenmode";
    case ProblemType::ELECTROSTATIC:
      return "electrostatic";
    case ProblemType::MAGNETOSTATIC:
      return "magnetostatic";
    case ProblemType::TRANSIENT:
      return "transient";
    default:
      return "unkown";
  }
}

}  // namespace

template <ProblemType solver_t>
PostOperator<solver_t>::PostOperator(const IoData &iodata, fem_op_t<solver_t> &fem_op_)
  : fem_op(&fem_op_), units(iodata.units), post_dir(iodata.problem.output),
    post_op_csv(iodata, fem_op_),
    // dom_post_op does not have a default ctor so specialize via immediate lambda.
    dom_post_op(std::move(
        [&iodata, &fem_op_]()
        {
          if constexpr (solver_t == ProblemType::ELECTROSTATIC)
          {
            return DomainPostOperator(iodata, fem_op_.GetMaterialOp(),
                                      fem_op_.GetH1Space());
          }
          else if constexpr (solver_t == ProblemType::MAGNETOSTATIC)
          {
            return DomainPostOperator(iodata, fem_op_.GetMaterialOp(),
                                      fem_op_.GetNDSpace());
          }
          else
          {
            return DomainPostOperator(iodata, fem_op_.GetMaterialOp(), fem_op_.GetNDSpace(),
                                      fem_op_.GetRTSpace());
          }
        }())),
    surf_post_op(iodata, fem_op->GetMaterialOp(), fem_op->GetH1Space()),
    interp_op(iodata, fem_op->GetNDSpace())
{
  // Define primary grid-functions.
  if constexpr (HasVGridFunction<solver_t>())
  {
    V = std::make_unique<GridFunction>(fem_op->GetH1Space());
  }
  if constexpr (HasAGridFunction<solver_t>())
  {
    A = std::make_unique<GridFunction>(fem_op->GetNDSpace());
  }
  if constexpr (HasEGridFunction<solver_t>())
  {
    E = std::make_unique<GridFunction>(fem_op->GetNDSpace(),
                                       HasComplexGridFunction<solver_t>());
  }
  if constexpr (HasBGridFunction<solver_t>())
  {
    B = std::make_unique<GridFunction>(fem_op->GetRTSpace(),
                                       HasComplexGridFunction<solver_t>());
  }

  // Add wave port boundary mode postprocessing, if available.
  if constexpr (std::is_same_v<fem_op_t<solver_t>, SpaceOperator>)
  {
    for (const auto &[idx, data] : fem_op->GetWavePortOp())
    {
      auto ret = port_E0.emplace(idx, WavePortFieldData());
      ret.first->second.E0r = data.GetModeFieldCoefficientReal();
      ret.first->second.E0i = data.GetModeFieldCoefficientImag();
    }
  }

  // If we write paraview fields, initialize paraview files and dependant measurements. We
  // currently don't use the dependant grid functions for non-paraview measurements, so only
  // initialize if needed.
  if (solver_t == ProblemType::DRIVEN)
  {
    paraview_save_indices = iodata.solver.driven.save_indices;
  }
  else if (solver_t == ProblemType::EIGENMODE)
  {
    paraview_n_post = iodata.solver.eigenmode.n_post;
  }
  else if (solver_t == ProblemType::ELECTROSTATIC)
  {
    paraview_n_post = iodata.solver.electrostatic.n_post;
  }
  else if (solver_t == ProblemType::MAGNETOSTATIC)
  {
    paraview_n_post = iodata.solver.magnetostatic.n_post;
  }
  else if (solver_t == ProblemType::TRANSIENT)
  {
    paraview_delta_post = iodata.solver.transient.delta_post;
  }
  InitializeParaviewDataCollection();

  // Initialize CSV files for measurements.
  post_op_csv.InitializeCSVDataCollection(*this);
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperator<solver_t>::InitializeParaviewDataCollection(int ex_idx)
    -> std::enable_if_t<U == ProblemType::DRIVEN, void>
{
  fs::path sub_folder_name = "";
  auto nr_excitations = fem_op->GetPortExcitations().Size();
  if ((nr_excitations > 1) && (ex_idx > 0))
  {
    int spacing = 1 + int(std::log10(nr_excitations));
    sub_folder_name = fmt::format(FMT_STRING("excitation_{:0>{}}"), ex_idx, spacing);
  }
  InitializeParaviewDataCollection(sub_folder_name);
}

template <ProblemType solver_t>
bool PostOperator<solver_t>::write_paraview_fields(std::size_t step)
{
  return (paraview_delta_post > 0 && step % paraview_delta_post == 0) ||
         (paraview_n_post > 0 && step < paraview_n_post) ||
         std::binary_search(paraview_save_indices.cbegin(), paraview_save_indices.cend(),
                            step);
}

template <ProblemType solver_t>
void PostOperator<solver_t>::InitializeParaviewDataCollection(
    const fs::path &sub_folder_name)
{
  if (!write_paraview_fields())
  {
    return;
  }
  fs::path paraview_dir_v = post_dir / "paraview" / ParaviewFoldername(solver_t);
  fs::path paraview_dir_b =
      post_dir / "paraview" / fmt::format("{}_boundary", ParaviewFoldername(solver_t));
  if (!sub_folder_name.empty())
  {
    paraview_dir_v /= sub_folder_name;
    paraview_dir_b /= sub_folder_name;
  }
  // Set up postprocessing for output to disk.
  paraview = {paraview_dir_v, &fem_op->GetNDSpace().GetParMesh()};
  paraview_bdr = {paraview_dir_b, &fem_op->GetNDSpace().GetParMesh()};

  // Set-up grid-functions for the paraview output / measurement.
  if constexpr (HasVGridFunction<solver_t>())
  {
    V_s = std::make_unique<BdrFieldCoefficient>(V->Real());
  }

  if constexpr (HasAGridFunction<solver_t>())
  {
    A_s = std::make_unique<BdrFieldVectorCoefficient>(A->Real());
  }

  if constexpr (HasEGridFunction<solver_t>())
  {
    // Electric Energy Density.
    U_e = std::make_unique<EnergyDensityCoefficient<EnergyDensityType::ELECTRIC>>(
        *E, fem_op->GetMaterialOp());

    // Electric Boundary Field & Surface Charge.
    E_sr = std::make_unique<BdrFieldVectorCoefficient>(E->Real());
    Q_sr = std::make_unique<BdrSurfaceFluxCoefficient<SurfaceFlux::ELECTRIC>>(
        &E->Real(), nullptr, fem_op->GetMaterialOp(), true, mfem::Vector());

    if constexpr (HasComplexGridFunction<solver_t>())
    {
      E_si = std::make_unique<BdrFieldVectorCoefficient>(E->Imag());
      Q_si = std::make_unique<BdrSurfaceFluxCoefficient<SurfaceFlux::ELECTRIC>>(
          &E->Imag(), nullptr, fem_op->GetMaterialOp(), true, mfem::Vector());
    }
  }

  if constexpr (HasBGridFunction<solver_t>())
  {
    // Magnetic Energy Density.
    U_m = std::make_unique<EnergyDensityCoefficient<EnergyDensityType::MAGNETIC>>(
        *B, fem_op->GetMaterialOp());

    // Magnetic Boundary Field & Surface Current.
    B_sr = std::make_unique<BdrFieldVectorCoefficient>(B->Real());
    J_sr = std::make_unique<BdrSurfaceCurrentVectorCoefficient>(B->Real(),
                                                                fem_op->GetMaterialOp());

    if constexpr (HasComplexGridFunction<solver_t>())
    {
      B_si = std::make_unique<BdrFieldVectorCoefficient>(B->Imag());
      J_si = std::make_unique<BdrSurfaceCurrentVectorCoefficient>(B->Imag(),
                                                                  fem_op->GetMaterialOp());
    }
  }

  if constexpr (HasEGridFunction<solver_t>() && HasBGridFunction<solver_t>())
  {
    // Poynting Vector.
    S = std::make_unique<PoyntingVectorCoefficient>(*E, *B, fem_op->GetMaterialOp());
  }

  const mfem::VTKFormat format = mfem::VTKFormat::BINARY32;
#if defined(MFEM_USE_ZLIB)
  const int compress = -1;  // Default compression level
#else
  const int compress = 0;
#endif
  const bool use_ho = true;
  const int refine_ho = HasEGridFunction<solver_t>()
                            ? E->ParFESpace()->GetMaxElementOrder()
                            : B->ParFESpace()->GetMaxElementOrder();

  // Output mesh coordinate units same as input.
  paraview->SetCycle(-1);
  paraview->SetDataFormat(format);
  paraview->SetCompressionLevel(compress);
  paraview->SetHighOrderOutput(use_ho);
  paraview->SetLevelsOfDetail(refine_ho);

  paraview_bdr->SetBoundaryOutput(true);
  paraview_bdr->SetCycle(-1);
  paraview_bdr->SetDataFormat(format);
  paraview_bdr->SetCompressionLevel(compress);
  paraview_bdr->SetHighOrderOutput(use_ho);
  paraview_bdr->SetLevelsOfDetail(refine_ho);

  // Output fields @ phase = 0 and π/2 for frequency domain (rather than, for example,
  // peak phasors or magnitude = sqrt(2) * RMS). Also output fields evaluated on mesh
  // boundaries. For internal boundary surfaces, this takes the field evaluated in the
  // neighboring element with the larger dielectric permittivity or magnetic
  // permeability.
  if (E)
  {
    if (HasComplexGridFunction<solver_t>())
    {
      paraview->RegisterField("E_real", &E->Real());
      paraview->RegisterField("E_imag", &E->Imag());
      paraview_bdr->RegisterVCoeffField("E_real", E_sr.get());
      paraview_bdr->RegisterVCoeffField("E_imag", E_si.get());
    }
    else
    {
      paraview->RegisterField("E", &E->Real());
      paraview_bdr->RegisterVCoeffField("E", E_sr.get());
    }
  }
  if (B)
  {
    if (HasComplexGridFunction<solver_t>())
    {
      paraview->RegisterField("B_real", &B->Real());
      paraview->RegisterField("B_imag", &B->Imag());
      paraview_bdr->RegisterVCoeffField("B_real", B_sr.get());
      paraview_bdr->RegisterVCoeffField("B_imag", B_si.get());
    }
    else
    {
      paraview->RegisterField("B", &B->Real());
      paraview_bdr->RegisterVCoeffField("B", B_sr.get());
    }
  }
  if (V)
  {
    paraview->RegisterField("V", &V->Real());
    paraview_bdr->RegisterCoeffField("V", V_s.get());
  }
  if (A)
  {
    paraview->RegisterField("A", &A->Real());
    paraview_bdr->RegisterVCoeffField("A", A_s.get());
  }

  // Extract energy density field for electric field energy 1/2 Dᴴ E or magnetic field
  // energy 1/2 Hᴴ B. Also Poynting vector S = E x H⋆.
  if (U_e)
  {
    paraview->RegisterCoeffField("U_e", U_e.get());
    paraview_bdr->RegisterCoeffField("U_e", U_e.get());
  }
  if (U_m)
  {
    paraview->RegisterCoeffField("U_m", U_m.get());
    paraview_bdr->RegisterCoeffField("U_m", U_m.get());
  }
  if (S)
  {
    paraview->RegisterVCoeffField("S", S.get());
    paraview_bdr->RegisterVCoeffField("S", S.get());
  }

  // Extract surface charge from normally discontinuous ND E-field. Also extract surface
  // currents from tangentially discontinuous RT B-field The surface charge and surface
  // currents are single-valued at internal boundaries.
  if (Q_sr)
  {
    if (HasComplexGridFunction<solver_t>())
    {
      paraview_bdr->RegisterCoeffField("Q_s_real", Q_sr.get());
      paraview_bdr->RegisterCoeffField("Q_s_imag", Q_si.get());
    }
    else
    {
      paraview_bdr->RegisterCoeffField("Q_s", Q_sr.get());
    }
  }
  if (J_sr)
  {
    if (HasComplexGridFunction<solver_t>())
    {
      paraview_bdr->RegisterVCoeffField("J_s_real", J_sr.get());
      paraview_bdr->RegisterVCoeffField("J_s_imag", J_si.get());
    }
    else
    {
      paraview_bdr->RegisterVCoeffField("J_s", J_sr.get());
    }
  }

  // Add wave port boundary mode postprocessing when available.
  for (const auto &[idx, data] : port_E0)
  {
    paraview_bdr->RegisterVCoeffField(fmt::format("E0_{}_real", idx), data.E0r.get());
    paraview_bdr->RegisterVCoeffField(fmt::format("E0_{}_imag", idx), data.E0i.get());
  }
}

namespace
{

template <typename T>
void ScaleGridFunctions(double L, int dim, bool imag, T &E, T &B, T &V, T &A)
{
  // For fields on H(curl) and H(div) spaces, we "undo" the effect of redimensionalizing
  // the mesh which would carry into the fields during the mapping from reference to
  // physical space through the element Jacobians. No transformation for V is needed (H1
  // interpolation). Because the coefficients are always evaluating E, B in neighboring
  // elements, the Jacobian scaling is the same for the domain and boundary data
  // collections (instead of being different for B due to the dim - 1 evaluation). Wave
  // port fields also do not require rescaling since their submesh object where they are
  // evaluated remains nondimensionalized.
  if (E)
  {
    // Piola transform: J^-T
    E->Real() *= L;
    E->Real().FaceNbrData() *= L;
    if (imag)
    {
      E->Imag() *= L;
      E->Imag().FaceNbrData() *= L;
    }
  }
  if (B)
  {
    // Piola transform: J / |J|
    const auto Ld = std::pow(L, dim - 1);
    B->Real() *= Ld;
    B->Real().FaceNbrData() *= Ld;
    if (imag)
    {
      B->Imag() *= Ld;
      B->Imag().FaceNbrData() *= Ld;
    }
  }
  if (A)
  {
    // Piola transform: J^-T
    A->Real() *= L;
    A->Real().FaceNbrData() *= L;
  }
}

}  // namespace

template <ProblemType solver_t>
void PostOperator<solver_t>::WriteFields(double time, int step)
{
  if (!write_paraview_fields())
  {
    return;
  }
  BlockTimer bt(Timer::IO);

  auto mesh_Lc0 = units.GetMeshLengthRelativeScale();

  // Given the electric field and magnetic flux density, write the fields to disk for
  // visualization. Write the mesh coordinates in the same units as originally input.
  mfem::ParMesh &mesh = E ? *E->ParFESpace()->GetParMesh() : *B->ParFESpace()->GetParMesh();
  mesh::DimensionalizeMesh(mesh, mesh_Lc0);
  ScaleGridFunctions(mesh_Lc0, mesh.Dimension(), HasComplexGridFunction<solver_t>(), E, B,
                     V, A);
  paraview->SetCycle(step);
  paraview->SetTime(time);
  paraview_bdr->SetCycle(step);
  paraview_bdr->SetTime(time);
  paraview->Save();
  paraview_bdr->Save();
  mesh::NondimensionalizeMesh(mesh, mesh_Lc0);
  ScaleGridFunctions(1.0 / mesh_Lc0, mesh.Dimension(), HasComplexGridFunction<solver_t>(),
                     E, B, V, A);
  Mpi::Barrier(fem_op->GetComm());
}

template <ProblemType solver_t>
void PostOperator<solver_t>::WriteFieldsFinal(const ErrorIndicator *indicator)
{
  if (!write_paraview_fields())
  {
    return;
  }

  BlockTimer bt(Timer::IO);

  auto mesh_Lc0 = units.GetMeshLengthRelativeScale();

  // Write the mesh partitioning and (optionally) error indicators at the final step. No
  // need for these to be parallel objects, since the data is local to each process and
  // there isn't a need to ever access the element neighbors. We set the time to some
  // non-used value to make the step identifiable within the data collection.
  mfem::ParMesh &mesh = E ? *E->ParFESpace()->GetParMesh() : *B->ParFESpace()->GetParMesh();
  mesh::DimensionalizeMesh(mesh, mesh_Lc0);
  paraview->SetCycle(paraview->GetCycle() + 1);
  if (paraview->GetTime() < 1.0)
  {
    paraview->SetTime(99.0);
  }
  else
  {
    // 1 -> 99, 10 -> 999, etc.
    paraview->SetTime(
        std::pow(10.0, 2.0 + static_cast<int>(std::log10(paraview->GetTime()))) - 1.0);
  }
  auto field_map = paraview->GetFieldMap();  // Copy, so can reregister later
  for (const auto &[name, gf] : field_map)
  {
    paraview->DeregisterField(name);
  }
  auto coeff_field_map = paraview->GetCoeffFieldMap();
  for (const auto &[name, gf] : coeff_field_map)
  {
    paraview->DeregisterCoeffField(name);
  }
  auto vcoeff_field_map = paraview->GetVCoeffFieldMap();
  for (const auto &[name, gf] : vcoeff_field_map)
  {
    paraview->DeregisterVCoeffField(name);
  }
  mfem::L2_FECollection pwconst_fec(0, mesh.Dimension());
  mfem::FiniteElementSpace pwconst_fespace(&mesh, &pwconst_fec);
  std::unique_ptr<mfem::GridFunction> rank, eta;
  {
    rank = std::make_unique<mfem::GridFunction>(&pwconst_fespace);
    *rank = mesh.GetMyRank() + 1;
    paraview->RegisterField("Rank", rank.get());
  }
  if (indicator)
  {
    eta = std::make_unique<mfem::GridFunction>(&pwconst_fespace);
    MFEM_VERIFY(eta->Size() == indicator->Local().Size(),
                "Size mismatch for provided ErrorIndicator for postprocessing!");
    *eta = indicator->Local();
    paraview->RegisterField("Indicator", eta.get());
  }
  paraview->Save();
  if (rank)
  {
    paraview->DeregisterField("Rank");
  }
  if (eta)
  {
    paraview->DeregisterField("Indicator");
  }
  for (const auto &[name, gf] : field_map)
  {
    paraview->RegisterField(name, gf);
  }
  for (const auto &[name, gf] : coeff_field_map)
  {
    paraview->RegisterCoeffField(name, gf);
  }
  for (const auto &[name, gf] : vcoeff_field_map)
  {
    paraview->RegisterVCoeffField(name, gf);
  }
  mesh::NondimensionalizeMesh(mesh, mesh_Lc0);
  Mpi::Barrier(fem_op->GetComm());
}

// Measurements.

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureDomainFieldEnergy() const
{
  measurement_cache.domain_E_field_energy_i.clear();
  measurement_cache.domain_H_field_energy_i.clear();

  measurement_cache.domain_E_field_energy_i.reserve(dom_post_op.M_i.size());
  measurement_cache.domain_H_field_energy_i.reserve(dom_post_op.M_i.size());

  if constexpr (HasEGridFunction<solver_t>())
  {
    // Use V if it has it rather than E.
    auto &field = V ? *V : *E;
    auto energy = dom_post_op.GetElectricFieldEnergy(field);
    measurement_cache.domain_E_field_energy_all = energy;

    for (const auto &[idx, data] : dom_post_op.M_i)
    {
      auto energy_i = dom_post_op.GetDomainElectricFieldEnergy(idx, field);
      auto participation_ratio = std::abs(energy_i) > 0.0 ? energy_i / energy : 0.0;
      measurement_cache.domain_E_field_energy_i.emplace_back(
          Measurement::DomainData{idx, energy_i, participation_ratio});
    }
  }
  else
  {
    // Magnetic field only.
    measurement_cache.domain_E_field_energy_all = 0.0;
    for (const auto &[idx, data] : dom_post_op.M_i)
    {
      measurement_cache.domain_E_field_energy_i.emplace_back(
          Measurement::DomainData{idx, 0.0, 0.0});
    }
  }

  if (HasBGridFunction<solver_t>())
  {
    auto &field = A ? *A : *B;
    auto energy = dom_post_op.GetMagneticFieldEnergy(field);
    measurement_cache.domain_H_field_energy_all = energy;

    for (const auto &[idx, data] : dom_post_op.M_i)
    {
      auto energy_i = dom_post_op.GetDomainMagneticFieldEnergy(idx, field);
      auto participation_ratio = std::abs(energy) > 0.0 ? energy_i / energy : 0.0;
      measurement_cache.domain_H_field_energy_i.emplace_back(
          Measurement::DomainData{idx, energy_i, participation_ratio});
    }
  }
  else
  {
    // Electric field only.
    measurement_cache.domain_H_field_energy_all = 0.0;
    for (const auto &[idx, data] : dom_post_op.M_i)
    {
      measurement_cache.domain_H_field_energy_i.emplace_back(
          Measurement::DomainData{idx, 0.0, 0.0});
    }
  }

  // Log Domain Energy.
  const auto domain_E = units.Dimensionalize<Units::ValueType::ENERGY>(
      measurement_cache.domain_E_field_energy_all);
  const auto domain_H = units.Dimensionalize<Units::ValueType::ENERGY>(
      measurement_cache.domain_H_field_energy_all);
  if constexpr (HasEGridFunction<solver_t>() && !HasBGridFunction<solver_t>())
  {
    Mpi::Print(" Field energy E = {:.3e} J\n", domain_E);
  }
  else if constexpr (!HasEGridFunction<solver_t>() && HasBGridFunction<solver_t>())
  {
    Mpi::Print(" Field energy H = {:.3e} J\n", domain_H);
  }
  else if constexpr (solver_t != ProblemType::EIGENMODE)
  {
    Mpi::Print(" Field energy E ({:.3e} J) + H ({:.3e} J) = {:.3e} J\n", domain_E, domain_H,
               domain_E + domain_H);
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureLumpedPorts() const
{
  measurement_cache.lumped_port_vi.clear();
  measurement_cache.lumped_port_inductor_energy = 0.0;
  measurement_cache.lumped_port_capacitor_energy = 0.0;

  if constexpr (solver_t == ProblemType::EIGENMODE || solver_t == ProblemType::DRIVEN ||
                solver_t == ProblemType::TRANSIENT)
  {
    for (const auto &[idx, data] : fem_op->GetLumpedPortOp())
    {
      auto &vi = measurement_cache.lumped_port_vi[idx];
      vi.P = data.GetPower(*E, *B);
      vi.V = data.GetVoltage(*E);
      if constexpr (solver_t == ProblemType::EIGENMODE || solver_t == ProblemType::DRIVEN)
      {
        // Compute current from the port impedance, separate contributions for R, L, C
        // branches.
        // Get value and make real: Matches current behaviour (even for eigensolver!).
        MFEM_VERIFY(
            measurement_cache.freq.real() > 0.0,
            "Frequency domain lumped port postprocessing requires nonzero frequency!");
        vi.I_RLC[0] =
            (std::abs(data.R) > 0.0)
                ? vi.V / data.GetCharacteristicImpedance(measurement_cache.freq.real(),
                                                         LumpedPortData::Branch::R)
                : 0.0;
        vi.I_RLC[1] =
            (std::abs(data.L) > 0.0)
                ? vi.V / data.GetCharacteristicImpedance(measurement_cache.freq.real(),
                                                         LumpedPortData::Branch::L)
                : 0.0;
        vi.I_RLC[2] =
            (std::abs(data.C) > 0.0)
                ? vi.V / data.GetCharacteristicImpedance(measurement_cache.freq.real(),
                                                         LumpedPortData::Branch::C)
                : 0.0;
        vi.I = std::accumulate(vi.I_RLC.begin(), vi.I_RLC.end(),
                               std::complex<double>{0.0, 0.0});
        vi.S = data.GetSParameter(*E);

        // Add contribution due to all inductive lumped boundaries in the model:
        //                      E_ind = ∑_j 1/2 L_j I_mj².
        if (std::abs(data.L) > 0.0)
        {
          std::complex<double> I_mj = vi.I_RLC[1];
          vi.inductor_energy = 0.5 * std::abs(data.L) * std::real(I_mj * std::conj(I_mj));
          measurement_cache.lumped_port_inductor_energy += vi.inductor_energy;
        }

        // Add contribution due to all capacitive lumped boundaries in the model:
        //                      E_cap = ∑_j 1/2 C_j V_mj².
        if (std::abs(data.C) > 0.0)
        {
          std::complex<double> V_mj = vi.V;
          vi.capacitor_energy = 0.5 * std::abs(data.C) * std::real(V_mj * std::conj(V_mj));
          measurement_cache.lumped_port_capacitor_energy += vi.capacitor_energy;
        }
      }
      else
      {
        // Compute current from P = V I^* since there is no frequency & characteristic
        // impedance of the lumped element.
        vi.I = (std::abs(vi.V) > 0.0) ? std::conj(vi.P / vi.V) : 0.0;
      }
    }
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureLumpedPortsEig() const
{
  // Depends on MeasureLumpedPorts.
  if constexpr (solver_t == ProblemType::EIGENMODE)
  {
    auto freq_re = measurement_cache.freq.real();
    auto energy_electric_all = measurement_cache.domain_E_field_energy_all +
                               measurement_cache.lumped_port_capacitor_energy;
    for (const auto &[idx, data] : fem_op->GetLumpedPortOp())
    {
      // Get previously computed data: should never fail as defined by MeasureLumpedPorts.
      auto &vi = measurement_cache.lumped_port_vi.at(idx);

      // Resistive Lumped Ports:
      // Compute participation ratio of external ports (given as any port boundary with
      // nonzero resistance). Currently no reactance of the ports is supported. The κ of
      // the port follows from:
      //                          κ_mj = 1/2 R_j I_mj² / E_m
      // from which the mode coupling quality factor is computed as:
      //                              Q_mj = ω_m / κ_mj.
      if (std::abs(data.R) > 0.0)
      {
        std::complex<double> I_mj = vi.I_RLC[0];
        // Power = 1/2 R_j I_mj².
        // Note conventions: mean(I²) = (I_r² + I_i²) / 2;
        auto resistor_power = 0.5 * std::abs(data.R) * std::real(I_mj * std::conj(I_mj));
        vi.mode_port_kappa =
            std::copysign(resistor_power / energy_electric_all, I_mj.real());
        vi.quality_factor = (vi.mode_port_kappa == 0.0)
                                ? mfem::infinity()
                                : freq_re / std::abs(vi.mode_port_kappa);
      }

      // Inductive Lumped Ports:
      // Compute energy-participation ratio of junction given by index idx for the field
      // mode. We first get the port line voltage, and use lumped port circuit impedance to
      // get peak current through the inductor: I_mj = V_mj / Z_mj,  Z_mj = i ω_m L_j. E_m
      // is the total energy in mode m: E_m = E_elec + E_cap = E_mag + E_ind. The signed EPR
      // for a lumped inductive element is computed as:
      //                            p_mj = 1/2 L_j I_mj² / E_m.
      // An element with no assigned inductance will be treated as having zero admittance
      // and thus zero current.
      if (std::abs(data.L) > 0.0)
      {
        std::complex<double> I_mj = vi.I_RLC[1];
        vi.inductive_energy_participation =
            std::copysign(vi.inductor_energy / energy_electric_all, I_mj.real());
      }
    }
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureWavePorts() const
{
  measurement_cache.wave_port_vi.clear();

  if constexpr (solver_t == ProblemType::DRIVEN)
  {
    for (const auto &[idx, data] : fem_op->GetWavePortOp())
    {
      // Get value and make real: Matches current behaviour.
      auto freq_re = measurement_cache.freq.real();  // TODO: Fix
      MFEM_VERIFY(freq_re > 0.0,
                  "Frequency domain wave port postprocessing requires nonzero frequency!");
      auto &vi = measurement_cache.wave_port_vi[idx];
      vi.P = data.GetPower(*E, *B);
      vi.S = data.GetSParameter(*E);
      // vi.V = vi.I[0] = vi.I[1] = vi.I[2] = 0.0;  // Not yet implemented
      //                                            // (Z = V² / P, I = V / Z)
    }
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureSParameter() const
{
  // Depends on LumpedPorts, WavePorts.
  if constexpr (solver_t == ProblemType::DRIVEN)
  {
    using fmt::format;
    using std::complex_literals::operator""i;

    // Don't measure S-Matrix unless there is only one excitation per port. Also, we current
    // don't support mixing wave and lumped ports, because we need to fix consistent
    // conventions / de-embedding.
    if (!fem_op->GetPortExcitations().IsMultipleSimple() ||
        !((fem_op->GetLumpedPortOp().Size() > 0) xor (fem_op->GetWavePortOp().Size() > 0)))
    {
      return;
    }

    // Assumes that for single driving port the excitation index is equal to the port index.
    auto drive_port_idx = measurement_cache.ex_idx;

    // Currently S-Parameters are not calculated for mixed lumped & wave ports, so don't
    // combine output iterators.
    for (const auto &[idx, data] : fem_op->GetLumpedPortOp())
    {
      // Get previously computed data: should never fail as defined by MeasureLumpedPorts.
      auto &vi = measurement_cache.lumped_port_vi.at(idx);

      const LumpedPortData &src_data = fem_op->GetLumpedPortOp().GetPort(drive_port_idx);
      if (idx == drive_port_idx)
      {
        vi.S.real(vi.S.real() - 1.0);
      }
      // Generalized S-parameters if the ports are resistive (avoids divide-by-zero).
      if (std::abs(data.R) > 0.0)
      {
        vi.S *= std::sqrt(src_data.R / data.R);
      }

      Mpi::Print(" {0} = {1:+.3e}{2:+.3e}i, |{0}| = {3:+.3e}, arg({0}) = {4:+.3e}\n",
                 format("S[{}][{}]", idx, drive_port_idx), vi.S.real(), vi.S.imag(),
                 Measurement::Magnitude(vi.S), Measurement::Phase(vi.S));
    }
    for (const auto &[idx, data] : fem_op->GetWavePortOp())
    {
      // Get previously computed data: should never fail as defined by MeasureWavePorts.
      auto &vi = measurement_cache.wave_port_vi.at(idx);

      // Wave port modes are not normalized to a characteristic impedance so no generalized
      // S-parameters are available.
      const WavePortData &src_data = fem_op->GetWavePortOp().GetPort(drive_port_idx);
      if (idx == drive_port_idx)
      {
        vi.S.real(vi.S.real() - 1.0);
      }
      // Port de-embedding: S_demb = S exp(ikₙᵢ dᵢ) exp(ikₙⱼ dⱼ) (distance offset is default
      // 0 unless specified).
      vi.S *= std::exp(1i * src_data.kn0 * src_data.d_offset);
      vi.S *= std::exp(1i * data.kn0 * data.d_offset);

      Mpi::Print(" {0} = {1:+.3e}{2:+.3e}i, |{0}| = {3:+.3e}, arg({0}) = {4:+.3e}\n",
                 format("S[{}][{}]", idx, drive_port_idx), vi.S.real(), vi.S.imag(),
                 Measurement::Magnitude(vi.S), Measurement::Phase(vi.S));
    }
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureSurfaceFlux() const
{
  // Compute the flux through a surface as Φ_j = ∫ F ⋅ n_j dS, with F = B, F = ε D, or F =
  // E x H. The special coefficient is used to avoid issues evaluating MFEM GridFunctions
  // which are discontinuous at interior boundary elements.
  measurement_cache.surface_flux_i.clear();
  measurement_cache.surface_flux_i.reserve(surf_post_op.flux_surfs.size());
  for (const auto &[idx, data] : surf_post_op.flux_surfs)
  {
    measurement_cache.surface_flux_i.emplace_back(Measurement::FluxData{
        idx, surf_post_op.GetSurfaceFlux(idx, E.get(), B.get()), data.type});
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureInterfaceEFieldEnergy() const
{
  // Depends on Lumped Port Energy since this is used in normalization of participation
  // ratio.

  // Compute the surface dielectric participation ratio and associated quality factor for
  // the material interface given by index idx. We have:
  //                            1/Q_mj = p_mj tan(δ)_j
  // with:
  //          p_mj = 1/2 t_j Re{∫_{Γ_j} (ε_j E_m)ᴴ E_m dS} / (E_elec + E_cap).
  measurement_cache.interface_eps_i.clear();
  if constexpr (HasEGridFunction<solver_t>())
  {
    // Domain and port energies must have been measured first. E_cap returns zero if the
    // solver does not support lumped ports.
    //
    // TODO: Should this not include other types of energy too (surface impedance case)?
    auto energy_electric_all = measurement_cache.domain_E_field_energy_all +
                               measurement_cache.lumped_port_capacitor_energy;

    measurement_cache.interface_eps_i.reserve(surf_post_op.eps_surfs.size());
    for (const auto &[idx, data] : surf_post_op.eps_surfs)
    {
      auto energy = surf_post_op.GetInterfaceElectricFieldEnergy(idx, *E);

      auto energy_participation_p = energy / energy_electric_all;
      auto loss_tangent_delta = surf_post_op.GetInterfaceLossTangent(idx);
      auto quality_factor_Q = (energy_participation_p == 0.0 || loss_tangent_delta == 0.0)
                                  ? mfem::infinity()
                                  : 1.0 / (loss_tangent_delta * energy_participation_p);

      measurement_cache.interface_eps_i.emplace_back(Measurement::InterfaceData{
          idx, energy, loss_tangent_delta, energy_participation_p, quality_factor_Q});
    }
  }
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureProbes() const
{
  measurement_cache.probe_E_field.clear();
  measurement_cache.probe_B_field.clear();

#if defined(MFEM_USE_GSLIB)
  if constexpr (HasEGridFunction<solver_t>())
  {
    if (interp_op.GetProbes().size() > 0)
    {
      measurement_cache.probe_E_field = interp_op.ProbeField(*E);
    }
  }
  if constexpr (HasBGridFunction<solver_t>())
  {
    if (interp_op.GetProbes().size() > 0)
    {
      measurement_cache.probe_B_field = interp_op.ProbeField(*B);
    }
  }
#endif
}

using fmt::format;

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperator<solver_t>::MeasureAndPrintAll(int ex_idx, int step,
                                                const ComplexVector &e,
                                                const ComplexVector &b,
                                                std::complex<double> omega)
    -> std::enable_if_t<U == ProblemType::DRIVEN, double>
{
  BlockTimer bt0(Timer::POSTPRO);
  SetEGridFunction(e);
  SetBGridFunction(b);

  measurement_cache = {};
  measurement_cache.freq = omega;
  measurement_cache.ex_idx = ex_idx;
  MeasureAllImpl();

  omega = units.Dimensionalize<Units::ValueType::FREQUENCY>(omega);
  post_op_csv.PrintAllCSVData(*this, measurement_cache, omega.real(), step, ex_idx);
  if (write_paraview_fields(step))
  {
    Mpi::Print("\n");
    auto ind = 1 + std::distance(paraview_save_indices.begin(),
                                 std::lower_bound(paraview_save_indices.begin(),
                                                  paraview_save_indices.end(), step));
    WriteFields(omega.real(), ind);
    Mpi::Print(" Wrote fields to disk at step {:d}\n", step + 1);
  }
  return measurement_cache.domain_E_field_energy_all +
         measurement_cache.domain_H_field_energy_all;
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperator<solver_t>::MeasureAndPrintAll(int step, const ComplexVector &e,
                                                const ComplexVector &b,
                                                std::complex<double> omega,
                                                double error_abs, double error_bkwd,
                                                int num_conv)
    -> std::enable_if_t<U == ProblemType::EIGENMODE, double>
{
  BlockTimer bt0(Timer::POSTPRO);
  SetEGridFunction(e);
  SetBGridFunction(b);

  measurement_cache = {};
  measurement_cache.freq = omega;
  measurement_cache.eigenmode_Q =
      (omega == 0.0) ? mfem::infinity() : 0.5 * std::abs(omega) / std::abs(omega.imag());
  measurement_cache.error_abs = error_abs;
  measurement_cache.error_bkwd = error_bkwd;

  // Mini pretty-print table of eig summaries: always print with header since other
  // measurements may log their results.
  if (Mpi::Root(fem_op->GetComm()))
  {
    Table table;
    int idx_pad = 1 + static_cast<int>(std::log10(num_conv));
    table.col_options = {6, 6};
    table.insert(Column("idx", "m", idx_pad, {}, {}) << step + 1);
    table.insert(Column("f_re", "Re{f} (GHz)")
                 << units.Dimensionalize<Units::ValueType::FREQUENCY>(omega.real()));
    table.insert(Column("f_im", "Im{f} (GHz)")
                 << units.Dimensionalize<Units::ValueType::FREQUENCY>(omega.imag()));
    table.insert(Column("q", "Q") << measurement_cache.eigenmode_Q);
    table.insert(Column("err_back", "Error (Bkwd.)") << error_bkwd);
    table.insert(Column("err_abs", "Error (Abs.)") << error_abs);
    table[0].print_as_int = true;
    Mpi::Print("{}", (step == 0) ? table.format_table() : table.format_row(0));
  }
  MeasureAllImpl();

  int print_idx = step + 1;
  post_op_csv.PrintAllCSVData(*this, measurement_cache, print_idx, step);
  if (write_paraview_fields(step))
  {
    WriteFields(step, print_idx);
    Mpi::Print(" Wrote mode {:d} to disk\n", print_idx);
  }
  return measurement_cache.domain_E_field_energy_all +
         measurement_cache.domain_H_field_energy_all;
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperator<solver_t>::MeasureAndPrintAll(int step, const Vector &v, const Vector &e,
                                                int idx)
    -> std::enable_if_t<U == ProblemType::ELECTROSTATIC, double>
{
  BlockTimer bt0(Timer::POSTPRO);
  SetVGridFunction(v);
  SetEGridFunction(e);

  measurement_cache = {};
  MeasureAllImpl();

  int print_idx = step + 1;
  post_op_csv.PrintAllCSVData(*this, measurement_cache, print_idx, step);
  if (write_paraview_fields(step))
  {
    Mpi::Print("\n");
    WriteFields(step, idx);
    Mpi::Print(" Wrote fields to disk for source {:d}\n", idx);
  }
  return measurement_cache.domain_E_field_energy_all +
         measurement_cache.domain_H_field_energy_all;
}
template <ProblemType solver_t>
template <ProblemType U>
auto PostOperator<solver_t>::MeasureAndPrintAll(int step, const Vector &a, const Vector &b,
                                                int idx)
    -> std::enable_if_t<U == ProblemType::MAGNETOSTATIC, double>
{
  BlockTimer bt0(Timer::POSTPRO);
  SetAGridFunction(a);
  SetBGridFunction(b);

  measurement_cache = {};
  MeasureAllImpl();

  int print_idx = step + 1;
  post_op_csv.PrintAllCSVData(*this, measurement_cache, print_idx, step);
  if (write_paraview_fields(step))
  {
    Mpi::Print("\n");
    WriteFields(step, idx);
    Mpi::Print(" Wrote fields to disk for source {:d}\n", idx);
  }
  return measurement_cache.domain_E_field_energy_all +
         measurement_cache.domain_H_field_energy_all;
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperator<solver_t>::MeasureAndPrintAll(int step, const Vector &e, const Vector &b,
                                                double time, double J_coef)
    -> std::enable_if_t<U == ProblemType::TRANSIENT, double>
{
  BlockTimer bt0(Timer::POSTPRO);
  SetEGridFunction(e);
  SetBGridFunction(b);

  measurement_cache = {};
  measurement_cache.Jcoeff_excitation = J_coef;
  MeasureAllImpl();

  // Time must be converted before passing into csv due to the shared PrintAllCSVData
  // method.
  time = units.Dimensionalize<Units::ValueType::TIME>(time);
  post_op_csv.PrintAllCSVData(*this, measurement_cache, time, step);
  if (write_paraview_fields(step))
  {
    Mpi::Print("\n");
    WriteFields(double(step) / paraview_delta_post, time);
    Mpi::Print(" Wrote fields to disk at step {:d}\n", step + 1);
  }
  return measurement_cache.domain_E_field_energy_all +
         measurement_cache.domain_H_field_energy_all;
}

template <ProblemType solver_t>
void PostOperator<solver_t>::MeasureFinalize(const ErrorIndicator &indicator)
{
  BlockTimer bt0(Timer::POSTPRO);
  auto indicator_stats = indicator.GetSummaryStatistics(fem_op->GetComm());
  post_op_csv.PrintErrorIndicator(Mpi::Root(fem_op->GetComm()), indicator_stats);
  if (write_paraview_fields())
  {
    WriteFieldsFinal(&indicator);
  }
}

template <ProblemType solver_t>
template <ProblemType U>
auto PostOperator<solver_t>::MeasureDomainFieldEnergyOnly(const ComplexVector &e,
                                                          const ComplexVector &b)
    -> std::enable_if_t<U == ProblemType::DRIVEN, double>
{
  SetEGridFunction(e);
  SetBGridFunction(b);
  MeasureDomainFieldEnergy();
  Mpi::Barrier(fem_op->GetComm());

  // Return total domain energy for normalizing error indicator.
  return measurement_cache.domain_E_field_energy_all +
         measurement_cache.domain_H_field_energy_all;
}

// Explict template instantiation.

template class PostOperator<ProblemType::DRIVEN>;
template class PostOperator<ProblemType::EIGENMODE>;
template class PostOperator<ProblemType::ELECTROSTATIC>;
template class PostOperator<ProblemType::MAGNETOSTATIC>;
template class PostOperator<ProblemType::TRANSIENT>;

// Function explict instantiation.
// TODO(C++20): with requires, we won't need a second template.

template auto PostOperator<ProblemType::DRIVEN>::MeasureAndPrintAll<ProblemType::DRIVEN>(
    int ex_idx, int step, const ComplexVector &e, const ComplexVector &b,
    std::complex<double> omega) -> double;

template auto
PostOperator<ProblemType::EIGENMODE>::MeasureAndPrintAll<ProblemType::EIGENMODE>(
    int step, const ComplexVector &e, const ComplexVector &b, std::complex<double> omega,
    double error_abs, double error_bkwd, int num_conv) -> double;

template auto
PostOperator<ProblemType::ELECTROSTATIC>::MeasureAndPrintAll<ProblemType::ELECTROSTATIC>(
    int step, const Vector &v, const Vector &e, int idx) -> double;

template auto
PostOperator<ProblemType::MAGNETOSTATIC>::MeasureAndPrintAll<ProblemType::MAGNETOSTATIC>(
    int step, const Vector &a, const Vector &b, int idx) -> double;

template auto
PostOperator<ProblemType::TRANSIENT>::MeasureAndPrintAll<ProblemType::TRANSIENT>(
    int step, const Vector &e, const Vector &b, double t, double J_coef) -> double;

template auto
PostOperator<ProblemType::DRIVEN>::MeasureDomainFieldEnergyOnly<ProblemType::DRIVEN>(
    const ComplexVector &e, const ComplexVector &b) -> double;

template auto
PostOperator<ProblemType::DRIVEN>::InitializeParaviewDataCollection(int ex_idx) -> void;

}  // namespace palace
