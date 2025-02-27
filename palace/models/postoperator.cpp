// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "postoperator.hpp"

#include "fem/coefficient.hpp"
#include "fem/errorindicator.hpp"
#include "models/curlcurloperator.hpp"
#include "models/laplaceoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/spaceoperator.hpp"
#include "models/surfacecurrentoperator.hpp"
#include "models/waveportoperator.hpp"
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

namespace
{

std::string CreateParaviewPath(const IoData &iodata, const std::string &name)
{
  return fs::path(iodata.problem.output) / "paraview" / name;
}

}  // namespace

PostOperator::PostOperator(const IoData &iodata, SpaceOperator &space_op,
                           const std::string &name)
  : mat_op(space_op.GetMaterialOp()),
    surf_post_op(iodata, space_op.GetMaterialOp(), space_op.GetH1Space()),
    dom_post_op(iodata, space_op.GetMaterialOp(), space_op.GetNDSpace(),
                space_op.GetRTSpace()),
    interp_op(iodata, space_op.GetNDSpace()),
    E(std::make_unique<GridFunction>(space_op.GetNDSpace(),
                                     iodata.problem.type !=
                                         config::ProblemData::Type::TRANSIENT)),
    B(std::make_unique<GridFunction>(space_op.GetRTSpace(),
                                     iodata.problem.type !=
                                         config::ProblemData::Type::TRANSIENT)),
    paraview(CreateParaviewPath(iodata, name), &space_op.GetNDSpace().GetParMesh()),
    paraview_bdr(CreateParaviewPath(iodata, name) + "_boundary",
                 &space_op.GetNDSpace().GetParMesh())
{
  U_e = std::make_unique<EnergyDensityCoefficient<EnergyDensityType::ELECTRIC>>(*E, mat_op);
  U_m = std::make_unique<EnergyDensityCoefficient<EnergyDensityType::MAGNETIC>>(*B, mat_op);
  S = std::make_unique<PoyntingVectorCoefficient>(*E, *B, mat_op);

  E_sr = std::make_unique<BdrFieldVectorCoefficient>(E->Real());
  B_sr = std::make_unique<BdrFieldVectorCoefficient>(B->Real());
  J_sr = std::make_unique<BdrSurfaceCurrentVectorCoefficient>(B->Real(), mat_op);
  Q_sr = std::make_unique<BdrSurfaceFluxCoefficient<SurfaceFluxType::ELECTRIC>>(
      &E->Real(), nullptr, mat_op, true, mfem::Vector());
  if (HasImag())
  {
    E_si = std::make_unique<BdrFieldVectorCoefficient>(E->Imag());
    B_si = std::make_unique<BdrFieldVectorCoefficient>(B->Imag());
    J_si = std::make_unique<BdrSurfaceCurrentVectorCoefficient>(B->Imag(), mat_op);
    Q_si = std::make_unique<BdrSurfaceFluxCoefficient<SurfaceFluxType::ELECTRIC>>(
        &E->Imag(), nullptr, mat_op, true, mfem::Vector());
  }

  // Add wave port boundary mode postprocessing when available.
  for (const auto &[idx, data] : space_op.GetWavePortOp())
  {
    auto ret = port_E0.emplace(idx, WavePortFieldData());
    ret.first->second.E0r = data.GetModeFieldCoefficientReal();
    ret.first->second.E0i = data.GetModeFieldCoefficientImag();
  }

  // Initialize data collection objects.
  InitializeDataCollection(iodata);
}

PostOperator::PostOperator(const IoData &iodata, LaplaceOperator &laplace_op,
                           const std::string &name)
  : mat_op(laplace_op.GetMaterialOp()),
    surf_post_op(iodata, laplace_op.GetMaterialOp(), laplace_op.GetH1Space()),
    dom_post_op(iodata, laplace_op.GetMaterialOp(), laplace_op.GetH1Space()),
    E(std::make_unique<GridFunction>(laplace_op.GetNDSpace())),
    V(std::make_unique<GridFunction>(laplace_op.GetH1Space())),
    paraview(CreateParaviewPath(iodata, name), &laplace_op.GetNDSpace().GetParMesh()),
    paraview_bdr(CreateParaviewPath(iodata, name) + "_boundary",
                 &laplace_op.GetNDSpace().GetParMesh()),
    interp_op(iodata, laplace_op.GetNDSpace())
{
  // Note: When using this constructor, you should not use any of the magnetic field related
  // postprocessing functions (magnetic field energy, inductor energy, surface currents,
  // etc.), since only V and E fields are supplied.
  U_e = std::make_unique<EnergyDensityCoefficient<EnergyDensityType::ELECTRIC>>(*E, mat_op);

  E_sr = std::make_unique<BdrFieldVectorCoefficient>(E->Real());
  V_s = std::make_unique<BdrFieldCoefficient>(V->Real());
  Q_sr = std::make_unique<BdrSurfaceFluxCoefficient<SurfaceFluxType::ELECTRIC>>(
      &E->Real(), nullptr, mat_op, true, mfem::Vector());

  // Initialize data collection objects.
  InitializeDataCollection(iodata);
}

PostOperator::PostOperator(const IoData &iodata, CurlCurlOperator &curlcurl_op,
                           const std::string &name)
  : mat_op(curlcurl_op.GetMaterialOp()),
    surf_post_op(iodata, curlcurl_op.GetMaterialOp(), curlcurl_op.GetH1Space()),
    dom_post_op(iodata, curlcurl_op.GetMaterialOp(), curlcurl_op.GetNDSpace()),
    B(std::make_unique<GridFunction>(curlcurl_op.GetRTSpace())),
    A(std::make_unique<GridFunction>(curlcurl_op.GetNDSpace())),
    paraview(CreateParaviewPath(iodata, name), &curlcurl_op.GetNDSpace().GetParMesh()),
    paraview_bdr(CreateParaviewPath(iodata, name) + "_boundary",
                 &curlcurl_op.GetNDSpace().GetParMesh()),
    interp_op(iodata, curlcurl_op.GetNDSpace())
{
  // Note: When using this constructor, you should not use any of the electric field related
  // postprocessing functions (electric field energy, capacitor energy, surface charge,
  // etc.), since only the B field is supplied.
  U_m = std::make_unique<EnergyDensityCoefficient<EnergyDensityType::MAGNETIC>>(*B, mat_op);

  B_sr = std::make_unique<BdrFieldVectorCoefficient>(B->Real());
  A_s = std::make_unique<BdrFieldVectorCoefficient>(A->Real());
  J_sr = std::make_unique<BdrSurfaceCurrentVectorCoefficient>(B->Real(), mat_op);

  // Initialize data collection objects.
  InitializeDataCollection(iodata);
}

void PostOperator::InitializeDataCollection(const IoData &iodata)
{
  // Set up postprocessing for output to disk. Results are stored in a directory at
  // `iodata.problem.output/paraview`.
  const mfem::VTKFormat format = mfem::VTKFormat::BINARY32;
#if defined(MFEM_USE_ZLIB)
  const int compress = -1;  // Default compression level
#else
  const int compress = 0;
#endif
  const bool use_ho = true;
  const int refine_ho = HasE() ? E->ParFESpace()->GetMaxElementOrder()
                               : B->ParFESpace()->GetMaxElementOrder();
  mesh_Lc0 = iodata.GetMeshLengthScale();

  // Output mesh coordinate units same as input.
  paraview.SetCycle(-1);
  paraview.SetDataFormat(format);
  paraview.SetCompressionLevel(compress);
  paraview.SetHighOrderOutput(use_ho);
  paraview.SetLevelsOfDetail(refine_ho);

  paraview_bdr.SetBoundaryOutput(true);
  paraview_bdr.SetCycle(-1);
  paraview_bdr.SetDataFormat(format);
  paraview_bdr.SetCompressionLevel(compress);
  paraview_bdr.SetHighOrderOutput(use_ho);
  paraview_bdr.SetLevelsOfDetail(refine_ho);

  // Output fields @ phase = 0 and π/2 for frequency domain (rather than, for example,
  // peak phasors or magnitude = sqrt(2) * RMS). Also output fields evaluated on mesh
  // boundaries. For internal boundary surfaces, this takes the field evaluated in the
  // neighboring element with the larger dielectric permittivity or magnetic
  // permeability.
  if (E)
  {
    if (HasImag())
    {
      paraview.RegisterField("E_real", &E->Real());
      paraview.RegisterField("E_imag", &E->Imag());
      paraview_bdr.RegisterVCoeffField("E_real", E_sr.get());
      paraview_bdr.RegisterVCoeffField("E_imag", E_si.get());
    }
    else
    {
      paraview.RegisterField("E", &E->Real());
      paraview_bdr.RegisterVCoeffField("E", E_sr.get());
    }
  }
  if (B)
  {
    if (HasImag())
    {
      paraview.RegisterField("B_real", &B->Real());
      paraview.RegisterField("B_imag", &B->Imag());
      paraview_bdr.RegisterVCoeffField("B_real", B_sr.get());
      paraview_bdr.RegisterVCoeffField("B_imag", B_si.get());
    }
    else
    {
      paraview.RegisterField("B", &B->Real());
      paraview_bdr.RegisterVCoeffField("B", B_sr.get());
    }
  }
  if (V)
  {
    paraview.RegisterField("V", &V->Real());
    paraview_bdr.RegisterCoeffField("V", V_s.get());
  }
  if (A)
  {
    paraview.RegisterField("A", &A->Real());
    paraview_bdr.RegisterVCoeffField("A", A_s.get());
  }

  // Extract energy density field for electric field energy 1/2 Dᴴ E or magnetic field
  // energy 1/2 Hᴴ B. Also Poynting vector S = E x H⋆.
  if (U_e)
  {
    paraview.RegisterCoeffField("U_e", U_e.get());
    paraview_bdr.RegisterCoeffField("U_e", U_e.get());
  }
  if (U_m)
  {
    paraview.RegisterCoeffField("U_m", U_m.get());
    paraview_bdr.RegisterCoeffField("U_m", U_m.get());
  }
  if (S)
  {
    paraview.RegisterVCoeffField("S", S.get());
    paraview_bdr.RegisterVCoeffField("S", S.get());
  }

  // Extract surface charge from normally discontinuous ND E-field. Also extract surface
  // currents from tangentially discontinuous RT B-field The surface charge and surface
  // currents are single-valued at internal boundaries.
  if (Q_sr)
  {
    if (HasImag())
    {
      paraview_bdr.RegisterCoeffField("Q_s_real", Q_sr.get());
      paraview_bdr.RegisterCoeffField("Q_s_imag", Q_si.get());
    }
    else
    {
      paraview_bdr.RegisterCoeffField("Q_s", Q_sr.get());
    }
  }
  if (J_sr)
  {
    if (HasImag())
    {
      paraview_bdr.RegisterVCoeffField("J_s_real", J_sr.get());
      paraview_bdr.RegisterVCoeffField("J_s_imag", J_si.get());
    }
    else
    {
      paraview_bdr.RegisterVCoeffField("J_s", J_sr.get());
    }
  }

  // Add wave port boundary mode postprocessing when available.
  for (const auto &[idx, data] : port_E0)
  {
    paraview_bdr.RegisterVCoeffField("E0_" + std::to_string(idx) + "_real", data.E0r.get());
    paraview_bdr.RegisterVCoeffField("E0_" + std::to_string(idx) + "_imag", data.E0i.get());
  }
}

void PostOperator::SetEGridFunction(const ComplexVector &e, bool exchange_face_nbr_data)
{
  MFEM_VERIFY(HasImag(),
              "SetEGridFunction for complex-valued output called when HasImag() == false!");
  MFEM_VERIFY(E, "Incorrect usage of PostOperator::SetEGridFunction!");
  E->Real().SetFromTrueDofs(e.Real());  // Parallel distribute
  E->Imag().SetFromTrueDofs(e.Imag());
  if (exchange_face_nbr_data)
  {
    E->Real().ExchangeFaceNbrData();  // Ready for parallel comm on shared faces
    E->Imag().ExchangeFaceNbrData();
  }
  ClearAllMeasurementCache();
}

void PostOperator::SetBGridFunction(const ComplexVector &b, bool exchange_face_nbr_data)
{
  MFEM_VERIFY(HasImag(),
              "SetBGridFunction for complex-valued output called when HasImag() == false!");
  MFEM_VERIFY(B, "Incorrect usage of PostOperator::SetBGridFunction!");
  B->Real().SetFromTrueDofs(b.Real());  // Parallel distribute
  B->Imag().SetFromTrueDofs(b.Imag());
  if (exchange_face_nbr_data)
  {
    B->Real().ExchangeFaceNbrData();  // Ready for parallel comm on shared faces
    B->Imag().ExchangeFaceNbrData();
  }
  ClearAllMeasurementCache();
}

void PostOperator::SetEGridFunction(const Vector &e, bool exchange_face_nbr_data)
{
  MFEM_VERIFY(!HasImag(),
              "SetEGridFunction for real-valued output called when HasImag() == true!");
  MFEM_VERIFY(E, "Incorrect usage of PostOperator::SetEGridFunction!");
  E->Real().SetFromTrueDofs(e);
  if (exchange_face_nbr_data)
  {
    E->Real().ExchangeFaceNbrData();
  }
  ClearAllMeasurementCache();
}

void PostOperator::SetBGridFunction(const Vector &b, bool exchange_face_nbr_data)
{
  MFEM_VERIFY(!HasImag(),
              "SetBGridFunction for real-valued output called when HasImag() == true!");
  MFEM_VERIFY(B, "Incorrect usage of PostOperator::SetBGridFunction!");
  B->Real().SetFromTrueDofs(b);
  if (exchange_face_nbr_data)
  {
    B->Real().ExchangeFaceNbrData();
  }
  ClearAllMeasurementCache();
}

void PostOperator::SetVGridFunction(const Vector &v, bool exchange_face_nbr_data)
{
  MFEM_VERIFY(!HasImag(),
              "SetVGridFunction for real-valued output called when HasImag() == true!");
  MFEM_VERIFY(V, "Incorrect usage of PostOperator::SetVGridFunction!");
  V->Real().SetFromTrueDofs(v);
  if (exchange_face_nbr_data)
  {
    V->Real().ExchangeFaceNbrData();
  }
  ClearAllMeasurementCache();
}

void PostOperator::SetAGridFunction(const Vector &a, bool exchange_face_nbr_data)
{
  MFEM_VERIFY(!HasImag(),
              "SetAGridFunction for real-valued output called when HasImag() == true!");
  MFEM_VERIFY(A, "Incorrect usage of PostOperator::SetAGridFunction!");
  A->Real().SetFromTrueDofs(a);
  if (exchange_face_nbr_data)
  {
    A->Real().ExchangeFaceNbrData();
  }
  ClearAllMeasurementCache();
}

// Measurements

void PostOperator::ClearAllMeasurementCache()
{
  // Clear Cache: Save omega since this set by hand like fields E,...
  auto omega = measurement_cache.omega;
  measurement_cache = {};
  measurement_cache.omega = omega;
}

void PostOperator::MeasureAll()
{
  ClearAllMeasurementCache();

  MeasureEFieldEnergy();
  MeasureHFieldEnergy();

  MeasureSurfaceFlux();
  MeasureInterfaceEFieldEnergy();
  MeasureProbes();
}

void PostOperator::MeasureAll(const SpaceOperator &space_op)
{
  MeasureAll();
  MeasureLumpedPorts(space_op.GetLumpedPortOp());
  MeasureWavePorts(space_op.GetWavePortOp());
}

void PostOperator::SetFrequency(double omega)
{
  measurement_cache.omega = std::complex<double>(omega);
}

void PostOperator::SetFrequency(std::complex<double> omega)
{
  measurement_cache.omega = omega;
}

std::complex<double> PostOperator::GetFrequency() const
{
  return measurement_cache.omega;
}

double PostOperator::GetEFieldEnergy() const
{
  return measurement_cache.domain_E_field_energy_all;
}

double PostOperator::GetEFieldEnergy(int idx) const
{
  auto it = measurement_cache.domain_E_field_energy_i.find(idx);
  if (it == measurement_cache.domain_E_field_energy_i.end())
  {
    MFEM_ABORT(fmt::format("Could not find domain index {} for E field energy!", idx));
  }
  return it->second;
}

void PostOperator::MeasureEFieldEnergy()
{
  measurement_cache.domain_E_field_energy_i.clear();
  for (const auto &[idx, data] : dom_post_op.M_i)
  {
    measurement_cache.domain_E_field_energy_i.emplace(
        idx, (!V && !E) ? 0.0 : dom_post_op.GetDomainElectricFieldEnergy(idx, V ? *V : *E));
  }

  measurement_cache.domain_E_field_energy_all =
      (!V && !E) ? 0.0 : dom_post_op.GetElectricFieldEnergy(V ? *V : *E);
}

void PostOperator::MeasureHFieldEnergy()
{
  measurement_cache.domain_H_field_energy_i.clear();
  for (const auto &[idx, data] : dom_post_op.M_i)
  {
    measurement_cache.domain_H_field_energy_i[idx] =
        (!A && !B) ? 0.0 : dom_post_op.GetDomainMagneticFieldEnergy(idx, A ? *A : *B);
  }

  measurement_cache.domain_H_field_energy_all =
      (!A && !B) ? 0.0 : dom_post_op.GetMagneticFieldEnergy(A ? *A : *B);
}

double PostOperator::GetHFieldEnergy() const
{
  return measurement_cache.domain_H_field_energy_all;
}

double PostOperator::GetHFieldEnergy(int idx) const
{
  auto it = measurement_cache.domain_H_field_energy_i.find(idx);
  if (it == measurement_cache.domain_H_field_energy_i.end())
  {
    MFEM_ABORT(fmt::format("Could not find domain index {} for H field energy!", idx));
  }
  return it->second;
}

void PostOperator::MeasureSurfaceFlux()
{
  // Compute the flux through a surface as Φ_j = ∫ F ⋅ n_j dS, with F = B, F = ε D, or F =
  // E x H. The special coefficient is used to avoid issues evaluating MFEM GridFunctions
  // which are discontinuous at interior boundary elements.
  measurement_cache.surface_flux_i.clear();
  if (!E && !B)
  {
    return;
  }
  measurement_cache.surface_flux_i.reserve(surf_post_op.flux_surfs.size());
  for (const auto &[idx, data] : surf_post_op.flux_surfs)
  {
    measurement_cache.surface_flux_i.emplace_back(
        FluxData{idx, surf_post_op.GetSurfaceFlux(idx, E.get(), B.get()), data.type});
  }
}

PostOperator::FluxData PostOperator::GetSurfaceFlux(int idx) const
{
  auto it = std::find_if(measurement_cache.surface_flux_i.begin(),
                         measurement_cache.surface_flux_i.end(),
                         [idx](const auto &d) { return d.idx == idx; });
  if (it == measurement_cache.surface_flux_i.end())
  {
    MFEM_ABORT(fmt::format("Could not find surface index {} for flux!", idx));
  }
  return *it;
}

void PostOperator::MeasureInterfaceEFieldEnergy()
{
  // Compute the surface dielectric participation ratio and associated quality factor for
  // the material interface given by index idx. We have:
  //                            1/Q_mj = p_mj tan(δ)_j
  // with:
  //          p_mj = 1/2 t_j Re{∫_{Γ_j} (ε_j E_m)ᴴ E_m dS} /(E_elec + E_cap).
  measurement_cache.interface_eps_i.clear();
  if (!E)
  {
    return;
  }
  measurement_cache.interface_eps_i.reserve(surf_post_op.eps_surfs.size());
  for (const auto &[idx, data] : surf_post_op.eps_surfs)
  {
    measurement_cache.interface_eps_i.emplace_back(
        InterfaceData{idx, surf_post_op.GetInterfaceElectricFieldEnergy(idx, *E),
                      surf_post_op.GetInterfaceLossTangent(idx)});
  }
}

const PostOperator::InterfaceData &PostOperator::GetInterfaceEFieldEnergy(int idx) const
{
  auto it = std::find_if(measurement_cache.interface_eps_i.begin(),
                         measurement_cache.interface_eps_i.end(),
                         [idx](const auto &d) { return d.idx == idx; });
  if (it == measurement_cache.interface_eps_i.end())
  {
    MFEM_ABORT(fmt::format("Could not find surface index {} for interface energy!", idx));
  }
  return *it;
}

double PostOperator::GetInterfaceParticipation(int idx, double E_m) const
{
  // Compute the surface dielectric participation ratio and associated quality factor for
  // the material interface given by index idx. We have:
  //                            1/Q_mj = p_mj tan(δ)_j
  // with:
  //          p_mj = 1/2 t_j Re{∫_{Γ_j} (ε_j E_m)ᴴ E_m dS} /(E_elec + E_cap).
  MFEM_VERIFY(E, "Surface Q not defined, no electric field solution found!");
  auto data = GetInterfaceEFieldEnergy(idx);
  return data.energy / E_m;
}

void PostOperator::MeasureLumpedPorts(const LumpedPortOperator &lumped_port_op)
{
  measurement_cache.lumped_port_vi.clear();
  if (!E || !B)
  {
    return;
  }
  for (const auto &[idx, data] : lumped_port_op)
  {
    auto &vi = measurement_cache.lumped_port_vi[idx];
    vi.P = data.GetPower(*E, *B);
    vi.V = data.GetVoltage(*E);
    if (HasImag())
    {
      // Compute current from the port impedance, separate contributions for R, L, C
      // branches.
      // Get value and make real: Matches current behaviour (even for eigensolver!)
      auto omega = GetFrequency().real();
      MFEM_VERIFY(
          omega > 0.0,
          "Frequency domain lumped port postprocessing requires nonzero frequency!");
      vi.I[0] =
          (std::abs(data.R) > 0.0)
              ? vi.V / data.GetCharacteristicImpedance(omega, LumpedPortData::Branch::R)
              : 0.0;
      vi.I[1] =
          (std::abs(data.L) > 0.0)
              ? vi.V / data.GetCharacteristicImpedance(omega, LumpedPortData::Branch::L)
              : 0.0;
      vi.I[2] =
          (std::abs(data.C) > 0.0)
              ? vi.V / data.GetCharacteristicImpedance(omega, LumpedPortData::Branch::C)
              : 0.0;
      vi.S = data.GetSParameter(*E);
    }
    else
    {
      // Compute current from P = V I⋆ (no scattering parameter output).
      vi.I[0] = (std::abs(vi.V) > 0.0) ? std::conj(vi.P / vi.V) : 0.0;
      vi.I[1] = vi.I[2] = vi.S = 0.0;
    }
  }

  // Add contribution due to all inductive lumped boundaries in the model:
  //                      E_ind = ∑_j 1/2 L_j I_mj².
  measurement_cache.lumped_port_inductor_energy = 0.0;
  for (const auto &[idx, data] : lumped_port_op)
  {
    if (std::abs(data.L) > 0.0)
    {
      std::complex<double> I_j = GetPortCurrent(idx, LumpedPortData::Branch::L);
      measurement_cache.lumped_port_inductor_energy +=
          0.5 * std::abs(data.L) * std::real(I_j * std::conj(I_j));
    }
  }

  // Add contribution due to all capacitive lumped boundaries in the model:
  //                      E_cap = ∑_j 1/2 C_j V_mj².
  measurement_cache.lumped_port_capacitor_energy = 0.0;
  for (const auto &[idx, data] : lumped_port_op)
  {
    if (std::abs(data.C) > 0.0)
    {
      std::complex<double> V_j = GetPortVoltage(idx);
      measurement_cache.lumped_port_capacitor_energy +=
          0.5 * std::abs(data.C) * std::real(V_j * std::conj(V_j));
    }
  }
}

void PostOperator::MeasureWavePorts(const WavePortOperator &wave_port_op)
{
  measurement_cache.wave_port_vi.clear();
  // Wave ports need imaginary component. TODO: Fix this.
  if (!E || !B || !HasImag())
  {
    return;
  }
  for (const auto &[idx, data] : wave_port_op)
  {
    // Get value and make real: Matches current behaviour
    auto omega = measurement_cache.omega.real();
    MFEM_VERIFY(omega > 0.0,
                "Frequency domain wave port postprocessing requires nonzero frequency!");
    auto &vi = measurement_cache.wave_port_vi[idx];
    vi.P = data.GetPower(*E, *B);
    vi.S = data.GetSParameter(*E);
    vi.V = vi.I[0] = vi.I[1] = vi.I[2] = 0.0;  // Not yet implemented
                                               // (Z = V² / P, I = V / Z)
  }
}

double PostOperator::GetLumpedInductorEnergy() const
{
  return measurement_cache.lumped_port_inductor_energy;
}

double PostOperator::GetLumpedCapacitorEnergy() const
{
  return measurement_cache.lumped_port_capacitor_energy;
}

std::complex<double> PostOperator::GetSParameter(const LumpedPortOperator &lumped_port_op,
                                                 int idx, int source_idx) const
{
  const LumpedPortData &data = lumped_port_op.GetPort(idx);
  const LumpedPortData &src_data = lumped_port_op.GetPort(source_idx);
  const auto it = measurement_cache.lumped_port_vi.find(idx);
  MFEM_VERIFY(src_data.excitation,
              "Lumped port index " << source_idx << " is not marked for excitation!");
  MFEM_VERIFY(it != measurement_cache.lumped_port_vi.end(),
              "Could not find lumped port when calculating port S-parameters!");
  std::complex<double> S_ij = it->second.S;
  if (idx == source_idx)
  {
    S_ij.real(S_ij.real() - 1.0);
  }
  // Generalized S-parameters if the ports are resistive (avoids divide-by-zero).
  if (std::abs(data.R) > 0.0)
  {
    S_ij *= std::sqrt(src_data.R / data.R);
  }
  return S_ij;
}

std::complex<double> PostOperator::GetSParameter(const WavePortOperator &wave_port_op,
                                                 int idx, int source_idx) const
{
  // Wave port modes are not normalized to a characteristic impedance so no generalized
  // S-parameters are available.
  const WavePortData &data = wave_port_op.GetPort(idx);
  const WavePortData &src_data = wave_port_op.GetPort(source_idx);
  const auto it = measurement_cache.wave_port_vi.find(idx);
  MFEM_VERIFY(src_data.excitation,
              "Wave port index " << source_idx << " is not marked for excitation!");
  MFEM_VERIFY(it != measurement_cache.wave_port_vi.end(),
              "Could not find wave port when calculating port S-parameters!");
  std::complex<double> S_ij = it->second.S;
  if (idx == source_idx)
  {
    S_ij.real(S_ij.real() - 1.0);
  }
  // Port de-embedding: S_demb = S exp(ikₙᵢ dᵢ) exp(ikₙⱼ dⱼ) (distance offset is default 0
  // unless specified).
  S_ij *= std::exp(1i * src_data.kn0 * src_data.d_offset);
  S_ij *= std::exp(1i * data.kn0 * data.d_offset);
  return S_ij;
}

std::complex<double> PostOperator::GetPortPower(int idx) const
{
  // TODO: In multi-excitation PR we will guarantee that lumped & wave ports have unique idx
  auto it_lumped = measurement_cache.lumped_port_vi.find(idx);
  if (it_lumped != measurement_cache.lumped_port_vi.end())
  {
    return it_lumped->second.P;
  }
  auto it_wave = measurement_cache.wave_port_vi.find(idx);
  if (it_wave != measurement_cache.wave_port_vi.end())
  {
    return it_wave->second.P;
  }
  MFEM_ABORT(
      fmt::format("Port Power: Could not find a lumped or wave port with index {}!", idx));
}

std::complex<double> PostOperator::GetPortVoltage(int idx) const
{
  // TODO: In multi-excitation PR we will guarantee that lumped & wave ports have unique idx
  auto it_lumped = measurement_cache.lumped_port_vi.find(idx);
  if (it_lumped != measurement_cache.lumped_port_vi.end())
  {
    return it_lumped->second.V;
  }
  auto it_wave = measurement_cache.wave_port_vi.find(idx);
  if (it_wave != measurement_cache.wave_port_vi.end())
  {
    MFEM_ABORT("GetPortVoltage is not yet implemented for wave port boundaries!");
  }
  MFEM_ABORT(fmt::format(
      "Port Voltage: Could not find a lumped or wave port with index {}!", idx));
}

std::complex<double> PostOperator::GetPortCurrent(int idx,
                                                  LumpedPortData::Branch branch) const
{
  // TODO: In multi-excitation PR we will guarantee that lumped & wave ports have unique idx
  auto it_lumped = measurement_cache.lumped_port_vi.find(idx);
  if (it_lumped != measurement_cache.lumped_port_vi.end())
  {
    auto &I_loc = it_lumped->second.I;
    switch (branch)
    {
      case LumpedPortData::Branch::R:
        return I_loc[0];
      case LumpedPortData::Branch::L:
        return I_loc[1];
      case LumpedPortData::Branch::C:
        return I_loc[2];
      default:
        return std::accumulate(I_loc.begin(), I_loc.end(), std::complex<double>{0.0, 0.0});
    }
  }
  auto it_wave = measurement_cache.wave_port_vi.find(idx);
  if (it_wave != measurement_cache.wave_port_vi.end())
  {
    MFEM_ABORT("GetPortCurrent is not yet implemented for wave port boundaries!");
  }
  MFEM_ABORT(fmt::format(
      "Port Current: Could not find a lumped or wave port with index {}!", idx));
}

double PostOperator::GetInductorParticipation(const LumpedPortOperator &lumped_port_op,
                                              int idx, double E_m) const
{
  // Compute energy-participation ratio of junction given by index idx for the field mode.
  // We first get the port line voltage, and use lumped port circuit impedance to get peak
  // current through the inductor: I_mj = V_mj / Z_mj,  Z_mj = i ω_m L_j. E_m is the total
  // energy in mode m: E_m = E_elec + E_cap = E_mag + E_ind. The signed EPR for a lumped
  // inductive element is computed as:
  //                            p_mj = 1/2 L_j I_mj² / E_m.
  // An element with no assigned inductance will be treated as having zero admittance and
  // thus zero current.
  const LumpedPortData &data = lumped_port_op.GetPort(idx);
  std::complex<double> I_mj = GetPortCurrent(idx, LumpedPortData::Branch::L);
  return std::copysign(0.5 * std::abs(data.L) * std::real(I_mj * std::conj(I_mj)) / E_m,
                       I_mj.real());  // mean(I²) = (I_r² + I_i²) / 2
}

double PostOperator::GetExternalKappa(const LumpedPortOperator &lumped_port_op, int idx,
                                      double E_m) const
{
  // Compute participation ratio of external ports (given as any port boundary with
  // nonzero resistance). Currently no reactance of the ports is supported. The κ of the
  // port follows from:
  //                          κ_mj = 1/2 R_j I_mj² / E_m
  // from which the mode coupling quality factor is computed as:
  //                              Q_mj = ω_m / κ_mj.
  const LumpedPortData &data = lumped_port_op.GetPort(idx);
  std::complex<double> I_mj = GetPortCurrent(idx, LumpedPortData::Branch::R);
  return std::copysign(0.5 * std::abs(data.R) * std::real(I_mj * std::conj(I_mj)) / E_m,
                       I_mj.real());  // mean(I²) = (I_r² + I_i²) / 2
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

void PostOperator::WriteFields(int step, double time) const
{
  BlockTimer bt(Timer::IO);

  // Given the electric field and magnetic flux density, write the fields to disk for
  // visualization. Write the mesh coordinates in the same units as originally input.
  mfem::ParMesh &mesh =
      HasE() ? *E->ParFESpace()->GetParMesh() : *B->ParFESpace()->GetParMesh();
  mesh::DimensionalizeMesh(mesh, mesh_Lc0);
  ScaleGridFunctions(mesh_Lc0, mesh.Dimension(), HasImag(), E, B, V, A);
  paraview.SetCycle(step);
  paraview.SetTime(time);
  paraview_bdr.SetCycle(step);
  paraview_bdr.SetTime(time);
  paraview.Save();
  paraview_bdr.Save();
  mesh::NondimensionalizeMesh(mesh, mesh_Lc0);
  ScaleGridFunctions(1.0 / mesh_Lc0, mesh.Dimension(), HasImag(), E, B, V, A);

  Mpi::Barrier(GetComm());
}

void PostOperator::WriteFieldsFinal(const ErrorIndicator *indicator) const
{
  BlockTimer bt(Timer::IO);

  // Write the mesh partitioning and (optionally) error indicators at the final step. No
  // need for these to be parallel objects, since the data is local to each process and
  // there isn't a need to ever access the element neighbors. We set the time to some
  // non-used value to make the step identifiable within the data collection.
  mfem::ParMesh &mesh =
      HasE() ? *E->ParFESpace()->GetParMesh() : *B->ParFESpace()->GetParMesh();
  mesh::DimensionalizeMesh(mesh, mesh_Lc0);
  paraview.SetCycle(paraview.GetCycle() + 1);
  if (paraview.GetTime() < 1.0)
  {
    paraview.SetTime(99.0);
  }
  else
  {
    // 1 -> 99, 10 -> 999, etc.
    paraview.SetTime(
        std::pow(10.0, 2.0 + static_cast<int>(std::log10(paraview.GetTime()))) - 1.0);
  }
  mfem::DataCollection::FieldMapType field_map(paraview.GetFieldMap());  // Copy
  for (const auto &[name, gf] : field_map)
  {
    paraview.DeregisterField(name);
  }
  mfem::DataCollection::CoeffFieldMapType coeff_field_map(paraview.GetCoeffFieldMap());
  for (const auto &[name, gf] : coeff_field_map)
  {
    paraview.DeregisterCoeffField(name);
  }
  mfem::DataCollection::VCoeffFieldMapType vcoeff_field_map(paraview.GetVCoeffFieldMap());
  for (const auto &[name, gf] : vcoeff_field_map)
  {
    paraview.DeregisterVCoeffField(name);
  }
  mfem::L2_FECollection pwconst_fec(0, mesh.Dimension());
  mfem::FiniteElementSpace pwconst_fespace(&mesh, &pwconst_fec);
  std::unique_ptr<mfem::GridFunction> rank, eta;
  {
    rank = std::make_unique<mfem::GridFunction>(&pwconst_fespace);
    *rank = mesh.GetMyRank() + 1;
    paraview.RegisterField("Rank", rank.get());
  }
  if (indicator)
  {
    eta = std::make_unique<mfem::GridFunction>(&pwconst_fespace);
    MFEM_VERIFY(eta->Size() == indicator->Local().Size(),
                "Size mismatch for provided ErrorIndicator for postprocessing!");
    *eta = indicator->Local();
    paraview.RegisterField("Indicator", eta.get());
  }
  paraview.Save();
  if (rank)
  {
    paraview.DeregisterField("Rank");
  }
  if (eta)
  {
    paraview.DeregisterField("Indicator");
  }
  for (const auto &[name, gf] : field_map)
  {
    paraview.RegisterField(name, gf);
  }
  for (const auto &[name, gf] : coeff_field_map)
  {
    paraview.RegisterCoeffField(name, gf);
  }
  for (const auto &[name, gf] : vcoeff_field_map)
  {
    paraview.RegisterVCoeffField(name, gf);
  }
  mesh::NondimensionalizeMesh(mesh, mesh_Lc0);

  Mpi::Barrier(GetComm());
}

void PostOperator::MeasureProbes()
{
  if (E && interp_op.GetProbes().size() > 0)
  {
    measurement_cache.probe_E_field = interp_op.ProbeField(*E);
  }
  if (B && interp_op.GetProbes().size() > 0)
  {
    measurement_cache.probe_B_field = interp_op.ProbeField(*B);
  }
}

std::vector<std::complex<double>> PostOperator::ProbeEField() const
{
  return measurement_cache.probe_E_field;
}
std::vector<std::complex<double>> PostOperator::ProbeBField() const
{
  return measurement_cache.probe_B_field;
}

}  // namespace palace
