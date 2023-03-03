// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "postoperator.hpp"

#include "fem/curlcurloperator.hpp"
#include "fem/laplaceoperator.hpp"
#include "fem/lumpedportoperator.hpp"
#include "fem/materialoperator.hpp"
#include "fem/spaceoperator.hpp"
#include "fem/surfacecurrentoperator.hpp"
#include "fem/waveportoperator.hpp"
#include "linalg/petsc.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/mfemcoefficients.hpp"

namespace palace
{

using namespace std::complex_literals;

namespace
{

auto LocalToShared(const mfem::ParMesh &mesh)
{
  // Construct shared face mapping required for boundary coefficients.
  std::map<int, int> l2s;
  for (int i = 0; i < mesh.GetNSharedFaces(); i++)
  {
    int i_local = mesh.GetSharedFace(i);
    l2s[i_local] = i;
  }
  return l2s;
}

auto CreateParaviewPath(const IoData &iodata, const std::string &name)
{
  std::string path = iodata.problem.output;
  if (path[path.length() - 1] != '/')
  {
    path += '/';
  }
  path += "paraview/" + name;
  return path;
}

}  // namespace

PostOperator::PostOperator(const IoData &iodata, SpaceOperator &spaceop,
                           const std::string &name)
  : local_to_shared(LocalToShared(*spaceop.GetNDSpace().GetParMesh())),
    mat_op(spaceop.GetMaterialOp()),
    surf_post_op(iodata, spaceop.GetMaterialOp(), local_to_shared, spaceop.GetH1Space()),
    dom_post_op(iodata, spaceop.GetMaterialOp(), &spaceop.GetNDSpace(),
                &spaceop.GetRTSpace()),
    has_imaginary(iodata.problem.type != config::ProblemData::Type::TRANSIENT),
    E(&spaceop.GetNDSpace()), B(&spaceop.GetRTSpace()), V(std::nullopt), A(std::nullopt),
    lumped_port_init(false), wave_port_init(false),
    paraview(CreateParaviewPath(iodata, name), spaceop.GetNDSpace().GetParMesh()),
    paraview_bdr(CreateParaviewPath(iodata, name) + "_boundary",
                 spaceop.GetNDSpace().GetParMesh()),
    interp_op(iodata, *spaceop.GetNDSpace().GetParMesh())
{
  Esr = std::make_unique<BdrFieldVectorCoefficient>(E->real(), mat_op, local_to_shared);
  Bsr = std::make_unique<BdrFieldVectorCoefficient>(B->real(), mat_op, local_to_shared);
  Jsr = std::make_unique<BdrCurrentVectorCoefficient>(B->real(), mat_op, local_to_shared);
  Qsr = std::make_unique<BdrChargeCoefficient>(E->real(), mat_op, local_to_shared);
  if (has_imaginary)
  {
    Esi = std::make_unique<BdrFieldVectorCoefficient>(E->imag(), mat_op, local_to_shared);
    Bsi = std::make_unique<BdrFieldVectorCoefficient>(B->imag(), mat_op, local_to_shared);
    Jsi = std::make_unique<BdrCurrentVectorCoefficient>(B->imag(), mat_op, local_to_shared);
    Qsi = std::make_unique<BdrChargeCoefficient>(E->imag(), mat_op, local_to_shared);
    Ue = std::make_unique<EnergyDensityCoefficient<EnergyDensityType::ELECTRIC,
                                                   EnergyDensityValueType::COMPLEX>>(
        *E, mat_op, local_to_shared);
    Um = std::make_unique<EnergyDensityCoefficient<EnergyDensityType::MAGNETIC,
                                                   EnergyDensityValueType::COMPLEX>>(
        *B, mat_op, local_to_shared);
  }
  else
  {
    Ue = std::make_unique<EnergyDensityCoefficient<EnergyDensityType::ELECTRIC,
                                                   EnergyDensityValueType::REAL>>(
        E->real(), mat_op, local_to_shared);
    Um = std::make_unique<EnergyDensityCoefficient<EnergyDensityType::MAGNETIC,
                                                   EnergyDensityValueType::REAL>>(
        B->real(), mat_op, local_to_shared);
  }

  // Initialize data collection objects and register additional fields associated with wave
  // ports (only constructed in SpaceOperator).
  InitializeDataCollection(iodata);
  for (const auto &[idx, data] : spaceop.GetWavePortOp())
  {
    paraview_bdr.RegisterVCoeffField("nxH^0_" + std::to_string(idx) + "_real",
                                     data.GetModeCoefficientReal().get());
    paraview_bdr.RegisterVCoeffField("nxH^0_" + std::to_string(idx) + "_imag",
                                     data.GetModeCoefficientImag().get());
  }
}

PostOperator::PostOperator(const IoData &iodata, LaplaceOperator &laplaceop,
                           const std::string &name)
  : local_to_shared(LocalToShared(*laplaceop.GetNDSpace().GetParMesh())),
    mat_op(laplaceop.GetMaterialOp()),
    surf_post_op(iodata, laplaceop.GetMaterialOp(), local_to_shared,
                 laplaceop.GetH1Space()),
    dom_post_op(iodata, laplaceop.GetMaterialOp(), &laplaceop.GetNDSpace(), nullptr),
    has_imaginary(false), E(&laplaceop.GetNDSpace()), B(std::nullopt),
    V(&laplaceop.GetH1Space()), A(std::nullopt), lumped_port_init(false),
    wave_port_init(false),
    paraview(CreateParaviewPath(iodata, name), laplaceop.GetNDSpace().GetParMesh()),
    paraview_bdr(CreateParaviewPath(iodata, name) + "_boundary",
                 laplaceop.GetNDSpace().GetParMesh()),
    interp_op(iodata, *laplaceop.GetNDSpace().GetParMesh())
{
  // Note: When using this constructor, you should not use any of the magnetic field related
  // postprocessing functions (magnetic field energy, inductor energy, surface currents,
  // etc.), since only V and E fields are supplied.
  Esr = std::make_unique<BdrFieldVectorCoefficient>(E->real(), mat_op, local_to_shared);
  Vs = std::make_unique<BdrFieldCoefficient>(*V, mat_op, local_to_shared);
  Ue = std::make_unique<
      EnergyDensityCoefficient<EnergyDensityType::ELECTRIC, EnergyDensityValueType::REAL>>(
      E->real(), mat_op, local_to_shared);
  Qsr = std::make_unique<BdrChargeCoefficient>(E->real(), mat_op, local_to_shared);

  // Initialize data collection objects.
  InitializeDataCollection(iodata);
}

PostOperator::PostOperator(const IoData &iodata, CurlCurlOperator &curlcurlop,
                           const std::string &name)
  : local_to_shared(LocalToShared(*curlcurlop.GetNDSpace().GetParMesh())),
    mat_op(curlcurlop.GetMaterialOp()),
    surf_post_op(iodata, curlcurlop.GetMaterialOp(), local_to_shared,
                 curlcurlop.GetH1Space()),
    dom_post_op(iodata, curlcurlop.GetMaterialOp(), nullptr, &curlcurlop.GetRTSpace()),
    has_imaginary(false), E(std::nullopt), B(&curlcurlop.GetRTSpace()), V(std::nullopt),
    A(&curlcurlop.GetNDSpace()), lumped_port_init(false), wave_port_init(false),
    paraview(CreateParaviewPath(iodata, name), curlcurlop.GetNDSpace().GetParMesh()),
    paraview_bdr(CreateParaviewPath(iodata, name) + "_boundary",
                 curlcurlop.GetNDSpace().GetParMesh()),
    interp_op(iodata, *curlcurlop.GetNDSpace().GetParMesh())
{
  // Note: When using this constructor, you should not use any of the electric field related
  // postprocessing functions (electric field energy, capacitor energy, surface charge,
  // etc.), since only the B field is supplied.
  Bsr = std::make_unique<BdrFieldVectorCoefficient>(B->real(), mat_op, local_to_shared);
  As = std::make_unique<BdrFieldVectorCoefficient>(*A, mat_op, local_to_shared);
  Um = std::make_unique<
      EnergyDensityCoefficient<EnergyDensityType::MAGNETIC, EnergyDensityValueType::REAL>>(
      B->real(), mat_op, local_to_shared);
  Jsr = std::make_unique<BdrCurrentVectorCoefficient>(B->real(), mat_op, local_to_shared);

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
  const int refine_ho =
      (E) ? E->ParFESpace()->GetMaxElementOrder() : B->ParFESpace()->GetMaxElementOrder();
  const double mesh_Lc0 =
      iodata.DimensionalizeValue(IoData::ValueType::LENGTH, 1.0 / iodata.model.L0);

  // Output mesh coordinate units same as input.
  paraview.SetDataFormat(format);
  paraview.SetCompressionLevel(compress);
  paraview.SetHighOrderOutput(use_ho);
  paraview.SetLevelsOfDetail(refine_ho);
  paraview.SetLengthScale(mesh_Lc0);

  paraview_bdr.SetBoundaryOutput(true);
  paraview_bdr.SetDataFormat(format);
  paraview_bdr.SetCompressionLevel(compress);
  paraview_bdr.SetHighOrderOutput(use_ho);
  paraview_bdr.SetLevelsOfDetail(refine_ho);
  paraview_bdr.SetLengthScale(mesh_Lc0);

  // Output fields @ phase = 0 and π/2 for frequency domain (rather than, for example,
  // peak phasors or magnitude = sqrt(2) * RMS). Also output fields evaluated on mesh
  // boundaries. For internal boundary surfaces, this takes the field evaluated in the
  // neighboring element with the larger dielectric permittivity or magnetic
  // permeability.
  if (E)
  {
    if (has_imaginary)
    {
      paraview.RegisterField("E_real", &E->real());
      paraview.RegisterField("E_imag", &E->imag());
      paraview_bdr.RegisterVCoeffField("E_real", Esr.get());
      paraview_bdr.RegisterVCoeffField("E_imag", Esi.get());
    }
    else
    {
      paraview.RegisterField("E", &E->real());
      paraview_bdr.RegisterVCoeffField("E", Esr.get());
    }
  }
  if (B)
  {
    if (has_imaginary)
    {
      paraview.RegisterField("B_real", &B->real());
      paraview.RegisterField("B_imag", &B->imag());
      paraview_bdr.RegisterVCoeffField("B_real", Bsr.get());
      paraview_bdr.RegisterVCoeffField("B_imag", Bsi.get());
    }
    else
    {
      paraview.RegisterField("B", &B->real());
      paraview_bdr.RegisterVCoeffField("B", Bsr.get());
    }
  }
  if (V)
  {
    paraview.RegisterField("V", &*V);
    paraview_bdr.RegisterCoeffField("V", Vs.get());
  }
  if (A)
  {
    paraview.RegisterField("A", &*A);
    paraview_bdr.RegisterVCoeffField("A", As.get());
  }

  // Extract surface charge from normally discontinuous ND E-field. Also extract surface
  // currents from tangentially discontinuous RT B-field The surface charge and surface
  // currents are single-valued at internal boundaries.
  if (Qsr)
  {
    if (has_imaginary)
    {
      paraview_bdr.RegisterCoeffField("Qs_real", Qsr.get());
      paraview_bdr.RegisterCoeffField("Qs_imag", Qsi.get());
    }
    else
    {
      paraview_bdr.RegisterCoeffField("Qs", Qsr.get());
    }
  }
  if (Jsr)
  {
    if (has_imaginary)
    {
      paraview_bdr.RegisterVCoeffField("Js_real", Jsr.get());
      paraview_bdr.RegisterVCoeffField("Js_imag", Jsi.get());
    }
    else
    {
      paraview_bdr.RegisterVCoeffField("Js", Jsr.get());
    }
  }

  // Extract energy density field for electric field energy 1/2 Dᴴ E or magnetic field
  // energy 1/2 Bᴴ H.
  if (Ue)
  {
    paraview.RegisterCoeffField("Ue", Ue.get());
    paraview_bdr.RegisterCoeffField("Ue", Ue.get());
  }
  if (Um)
  {
    paraview.RegisterCoeffField("Um", Um.get());
    paraview_bdr.RegisterCoeffField("Um", Um.get());
  }
}

void PostOperator::GetBField(std::complex<double> omega,
                             const petsc::PetscParMatrix &NegCurl,
                             const petsc::PetscParVector &e, petsc::PetscParVector &b)
{
  // Compute B = -1/(iω) ∇ x E on the true dofs.
  MFEM_VERIFY(e.GetSize() == NegCurl.Width() && b.GetSize() == NegCurl.Height(),
              "Size mismatch error computing B-field in PostOperator!");
  NegCurl.Mult(e, b);
  b.Scale(1.0 / (1i * omega));
}

void PostOperator::GetBField(const mfem::Operator &Curl, const mfem::Vector &a,
                             mfem::Vector &b)
{
  // Compute B = ∇ x A on the true dofs.
  MFEM_VERIFY(a.Size() == Curl.Width() && b.Size() == Curl.Height(),
              "Size mismatch error computing B-field in PostOperator!");
  Curl.Mult(a, b);
}

void PostOperator::GetEField(const mfem::Operator &NegGrad, const mfem::Vector &v,
                             mfem::Vector &e)
{
  // Compute E = -∇V on the true dofs.
  MFEM_VERIFY(v.Size() == NegGrad.Width() && e.Size() == NegGrad.Height(),
              "Size mismatch error computing E-field in PostOperator!");
  NegGrad.Mult(v, e);
}

void PostOperator::SetEGridFunction(const petsc::PetscParVector &e)
{
  MFEM_VERIFY(
      has_imaginary,
      "SetEGridFunction for complex-valued output called when has_imaginary == false!");
  MFEM_VERIFY(E, "Incorrect usage of PostOperator::SetEGridFunction!");
  mfem::Vector Er(e.GetSize()), Ei(e.GetSize());
  e.GetToVectors(Er, Ei);
  E->real().SetFromTrueDofs(Er);  // Parallel distribute
  E->imag().SetFromTrueDofs(Ei);
  E->real().ExchangeFaceNbrData();  // Ready for parallel comm on shared faces
  E->imag().ExchangeFaceNbrData();
  lumped_port_init = wave_port_init = false;
}

void PostOperator::SetBGridFunction(const petsc::PetscParVector &b)
{
  MFEM_VERIFY(
      has_imaginary,
      "SetBGridFunction for complex-valued output called when has_imaginary == false!");
  MFEM_VERIFY(B, "Incorrect usage of PostOperator::SetBGridFunction!");
  mfem::Vector Br(b.GetSize()), Bi(b.GetSize());
  b.GetToVectors(Br, Bi);
  B->real().SetFromTrueDofs(Br);  // Parallel distribute
  B->imag().SetFromTrueDofs(Bi);
  B->real().ExchangeFaceNbrData();  // Ready for parallel comm on shared faces
  B->imag().ExchangeFaceNbrData();
  lumped_port_init = wave_port_init = false;
}

void PostOperator::SetEGridFunction(const mfem::Vector &e)
{
  MFEM_VERIFY(!has_imaginary,
              "SetEGridFunction for real-valued output called when has_imaginary == true!");
  MFEM_VERIFY(E, "Incorrect usage of PostOperator::SetEGridFunction!");
  E->real().SetFromTrueDofs(e);
  E->real().ExchangeFaceNbrData();
  lumped_port_init = wave_port_init = false;
}

void PostOperator::SetBGridFunction(const mfem::Vector &b)
{
  MFEM_VERIFY(!has_imaginary,
              "SetBGridFunction for real-valued output called when has_imaginary == true!");
  MFEM_VERIFY(B, "Incorrect usage of PostOperator::SetBGridFunction!");
  B->real().SetFromTrueDofs(b);
  B->real().ExchangeFaceNbrData();
  lumped_port_init = wave_port_init = false;
}

void PostOperator::SetVGridFunction(const mfem::Vector &v)
{
  MFEM_VERIFY(!has_imaginary,
              "SetVGridFunction for real-valued output called when has_imaginary == true!");
  MFEM_VERIFY(V, "Incorrect usage of PostOperator::SetVGridFunction!");
  V->SetFromTrueDofs(v);
  V->ExchangeFaceNbrData();
}

void PostOperator::SetAGridFunction(const mfem::Vector &a)
{
  MFEM_VERIFY(!has_imaginary,
              "SetAGridFunction for real-valued output called when has_imaginary == true!");
  MFEM_VERIFY(A, "Incorrect usage of PostOperator::SetAGridFunction!");
  A->SetFromTrueDofs(a);
  A->ExchangeFaceNbrData();
}

void PostOperator::UpdatePorts(const LumpedPortOperator &lumped_port_op, double omega)
{
  MFEM_VERIFY(E && B, "Incorrect usage of PostOperator::UpdatePorts!");
  if (lumped_port_init)
  {
    return;
  }
  for (const auto &[idx, data] : lumped_port_op)
  {
    auto &vi = lumped_port_vi[idx];
    if (has_imaginary)
    {
      MFEM_VERIFY(
          omega > 0.0,
          "Frequency domain lumped port postprocessing requires nonzero frequency!");
      vi.S = data.GetSParameter(*E);
      vi.P = data.GetPower(*E, *B, mat_op, local_to_shared);
      vi.V = data.GetVoltage(*E);
      vi.Z = data.GetCharacteristicImpedance(omega);
    }
    else
    {
      vi.P = data.GetPower(E->real(), B->real(), mat_op, local_to_shared);
      vi.V = data.GetVoltage(E->real());
      vi.S = vi.Z = 0.0;
    }
  }
  lumped_port_init = true;
}

void PostOperator::UpdatePorts(const WavePortOperator &wave_port_op, double omega)
{
  MFEM_VERIFY(has_imaginary && E && B, "Incorrect usage of PostOperator::UpdatePorts!");
  if (wave_port_init)
  {
    return;
  }
  for (const auto &[idx, data] : wave_port_op)
  {
    MFEM_VERIFY(omega > 0.0,
                "Frequency domain wave port postprocessing requires nonzero frequency!");
    auto &vi = wave_port_vi[idx];
    vi.S = data.GetSParameter(*E);
    vi.P = data.GetPower(*E, *B, mat_op, local_to_shared);
    vi.V = vi.Z = 0.0;  // Not yet implemented (Z = V² / P, I = V / Z)
  }
  wave_port_init = true;
}

double PostOperator::GetEFieldEnergy() const
{
  // We use a leading factor of 1/2 instead of 1/4 even though the eigenmodes are peak
  // phasors and not RMS normalized because the same peak phasors are used to compute the
  // voltages/currents which are 2x the time-averaged values. This correctly yields an EPR
  // of 1 in cases where expected.
  MFEM_VERIFY(E, "PostOperator is not configured for electric field energy calculation!");
  return has_imaginary ? dom_post_op.GetElectricFieldEnergy(*E)
                       : dom_post_op.GetElectricFieldEnergy(E->real());
}

double PostOperator::GetHFieldEnergy() const
{
  // We use a leading factor of 1/2 instead of 1/4 even though the eigenmodes are peak
  // phasors and not RMS normalized because the same peak phasors are used to compute the
  // voltages/currents which are 2x the time-averaged values. This correctly yields an EPR
  // of 1 in cases where expected.
  MFEM_VERIFY(B, "PostOperator is not configured for magnetic field energy calculation!");
  return has_imaginary ? dom_post_op.GetMagneticFieldEnergy(*B)
                       : dom_post_op.GetMagneticFieldEnergy(B->real());
}

double PostOperator::GetLumpedInductorEnergy(const LumpedPortOperator &lumped_port_op) const
{
  // Add contribution due to all capacitive lumped boundaries in the model:
  //                      E_ind = ∑_j 1/2 L_j I_mj².
  double U = 0.0;
  for (const auto &[idx, data] : lumped_port_op)
  {
    if (std::abs(data.GetL()) > 0.0)
    {
      std::complex<double> Ij = GetPortCurrent(lumped_port_op, idx);
      U += 0.5 * std::abs(data.GetL()) * std::real(Ij * std::conj(Ij));
    }
  }
  return U;
}

double
PostOperator::GetLumpedCapacitorEnergy(const LumpedPortOperator &lumped_port_op) const
{
  // Add contribution due to all capacitive lumped boundaries in the model:
  //                      E_cap = ∑_j 1/2 C_j V_mj².
  double U = 0.0;
  for (const auto &[idx, data] : lumped_port_op)
  {
    if (std::abs(data.GetC()) > 0.0)
    {
      std::complex<double> Vj = GetPortVoltage(lumped_port_op, idx);
      U += 0.5 * std::abs(data.GetC()) * std::real(Vj * std::conj(Vj));
    }
  }
  return U;
}

std::complex<double> PostOperator::GetSParameter(const LumpedPortOperator &lumped_port_op,
                                                 int idx, int source_idx) const
{
  MFEM_VERIFY(lumped_port_init,
              "Port S-parameters not defined until ports are initialized!");
  const LumpedPortData &data = lumped_port_op.GetPort(idx);
  const LumpedPortData &src_data = lumped_port_op.GetPort(source_idx);
  const auto it = lumped_port_vi.find(idx);
  MFEM_VERIFY(src_data.IsExcited(),
              "Lumped port index " << source_idx << " is not marked for excitation!");
  MFEM_VERIFY(it != lumped_port_vi.end(),
              "Could not find lumped port when calculating port S-parameters!");
  std::complex<double> Sij = it->second.S;
  if (idx == source_idx)
  {
    Sij.real(Sij.real() - 1.0);
  }
  // Generalized S-parameters if the ports are resistive (avoids divide-by-zero).
  if (std::abs(data.GetR()) > 0.0)
  {
    Sij *= std::sqrt(src_data.GetR() / data.GetR());
  }
  return Sij;
}

std::complex<double> PostOperator::GetSParameter(const WavePortOperator &wave_port_op,
                                                 int idx, int source_idx) const
{
  // Wave port modes are not normalized to a characteristic impedance so no generalized
  // S-parameters are available.
  MFEM_VERIFY(wave_port_init, "Port S-parameters not defined until ports are initialized!");
  const WavePortData &data = wave_port_op.GetPort(idx);
  const WavePortData &src_data = wave_port_op.GetPort(source_idx);
  const auto it = wave_port_vi.find(idx);
  MFEM_VERIFY(src_data.IsExcited(),
              "Wave port index " << source_idx << " is not marked for excitation!");
  MFEM_VERIFY(it != wave_port_vi.end(),
              "Could not find wave port when calculating port S-parameters!");
  std::complex<double> Sij = it->second.S;
  if (idx == source_idx)
  {
    Sij.real(Sij.real() - 1.0);
  }
  // Port de-embedding: S_demb = S exp(-ikₙᵢ dᵢ) exp(-ikₙⱼ dⱼ) (distance offset is default
  // 0 unless specified).
  Sij *= std::exp(1i * src_data.GetPropagationConstant() * src_data.GetOffsetDistance());
  Sij *= std::exp(1i * data.GetPropagationConstant() * data.GetOffsetDistance());
  return Sij;
}

std::complex<double> PostOperator::GetPortPower(const LumpedPortOperator &lumped_port_op,
                                                int idx) const
{
  MFEM_VERIFY(lumped_port_init,
              "Lumped port quantities not defined until ports are initialized!");
  const auto it = lumped_port_vi.find(idx);
  MFEM_VERIFY(it != lumped_port_vi.end(),
              "Could not find lumped port when calculating lumped port power!");
  return it->second.P;
}

std::complex<double> PostOperator::GetPortPower(const WavePortOperator &wave_port_op,
                                                int idx) const
{
  MFEM_VERIFY(wave_port_init,
              "Wave port quantities not defined until ports are initialized!");
  const auto it = wave_port_vi.find(idx);
  MFEM_VERIFY(it != wave_port_vi.end(),
              "Could not find wave port when calculating wave port power!");
  return it->second.P;
}

std::complex<double> PostOperator::GetPortVoltage(const LumpedPortOperator &lumped_port_op,
                                                  int idx) const
{
  MFEM_VERIFY(lumped_port_init,
              "Lumped port quantities not defined until ports are initialized!");
  const auto it = lumped_port_vi.find(idx);
  MFEM_VERIFY(it != lumped_port_vi.end(),
              "Could not find lumped port when calculating lumped port voltage!");
  return it->second.V;
}

std::complex<double> PostOperator::GetPortCurrent(const LumpedPortOperator &lumped_port_op,
                                                  int idx) const
{
  MFEM_VERIFY(lumped_port_init,
              "Lumped port quantities not defined until ports are initialized!");
  const auto it = lumped_port_vi.find(idx);
  MFEM_VERIFY(it != lumped_port_vi.end(),
              "Could not find lumped port when calculating lumped port current!");
  if (std::abs(it->second.Z) > 0.0)
  {
    // Compute from V = I Z when impedance is available.
    return it->second.V / it->second.Z;
  }
  else if (std::abs(it->second.V) > 0.0)
  {
    // Compute from P = V I⋆.
    return std::conj(it->second.P / it->second.V);
  }
  return 0.0;
}

double PostOperator::GetInductorParticipation(const LumpedPortOperator &lumped_port_op,
                                              int idx, double Em) const
{
  // Compute energy-participation ratio of junction given by index idx for the field mode.
  // We first get the port line voltage, and use lumped port circuit impedance to get peak
  // current through the inductor: I_mj = V_mj / Z_mj,  Z_mj = i ω_m L_j. Em is the total
  // energy in mode m: E_m = E_elec + E_cap = E_mag + E_ind. The signed EPR for a lumped
  // inductive element is computed as:
  //                            p_mj = 1/2 L_j I_mj² / E_m.
  // An element with no assigned inductance will be treated as having zero admittance and
  // thus zero current.
  const LumpedPortData &data = lumped_port_op.GetPort(idx);
  std::complex<double> Imj = GetPortCurrent(lumped_port_op, idx);
  return std::copysign(0.5 * std::abs(data.GetL()) * std::real(Imj * std::conj(Imj)) / Em,
                       Imj.real());  // mean(I²) = (I_r² + I_i²) / 2
}

double PostOperator::GetExternalKappa(const LumpedPortOperator &lumped_port_op, int idx,
                                      double Em) const
{
  // Compute participation ratio of external ports (given as any port boundary with nonzero
  // resistance). Currently no reactance of the ports is supported. The κ of the port
  // follows from:
  //                          κ_mj = 1/2 R_j I_mj² / E_m
  // from which the mode coupling quality factor is computed as:
  //                              Q_mj = ω_m / κ_mj.
  const LumpedPortData &data = lumped_port_op.GetPort(idx);
  std::complex<double> Imj = GetPortCurrent(lumped_port_op, idx);
  return std::copysign(0.5 * std::abs(data.GetR()) * std::real(Imj * std::conj(Imj)) / Em,
                       Imj.real());  // mean(I²) = (I_r² + I_i²) / 2
}

double PostOperator::GetBulkParticipation(int idx, double Em) const
{
  // Compute the bulk dielectric participation ratio material given by index idx. Here, we
  // have:
  //                     p_mj = E_elec,j / (E_elec + E_cap).
  MFEM_VERIFY(E, "Bulk Q not defined, no electric field solution found!");
  double Ebulk = has_imaginary ? dom_post_op.GetDomainElectricFieldEnergy(idx, *E)
                               : dom_post_op.GetDomainElectricFieldEnergy(idx, E->real());
  return Ebulk / Em;
}

double PostOperator::GetBulkQualityFactor(int idx, double Em) const
{
  // Compute the associated quality factor for the material given by index idx. Here, we
  // have:
  //             1/Q_mj = p_mj tan(δ)_j = tan(δ)_j E_elec,j / (E_elec + E_cap).
  MFEM_VERIFY(E, "Bulk Q not defined, no electric field solution found!");
  double Ebulki = has_imaginary
                      ? dom_post_op.GetDomainElectricFieldEnergyLoss(idx, *E)
                      : dom_post_op.GetDomainElectricFieldEnergyLoss(idx, E->real());
  return (Ebulki == 0.0) ? mfem::infinity() : Em / Ebulki;
}

double PostOperator::GetInterfaceParticipation(int idx, double Em) const
{
  // Compute the surface dielectric participation ratio and associated quality factor for
  // the material interface given by index idx. We have:
  //                            1/Q_mj = p_mj tan(δ)_j
  // with:
  //          p_mj = 1/2 t_j Re{∫_{Γ_j} (ε_j E_m)ᴴ E_m dS} /(E_elec + E_cap).
  MFEM_VERIFY(E, "Surface Q not defined, no electric field solution found!");
  double Esurf = surf_post_op.GetInterfaceElectricFieldEnergy(idx, E->real());
  if (has_imaginary)
  {
    Esurf += surf_post_op.GetInterfaceElectricFieldEnergy(idx, E->imag());
  }
  return Esurf / Em;
}

double PostOperator::GetSurfaceCharge(int idx) const
{
  // Compute the induced charge on a surface as Q_j = ∫ D ⋅ n_j dS, which correctly handles
  // two-sided internal surfaces using a special GridFunction coefficient which accounts
  // for both sides of the surface. This then yields the capacitive coupling to the
  // excitation as C_jk = Q_j / V_k where V_k is the excitation voltage.
  MFEM_VERIFY(E, "Surface capacitance not defined, no electric field solution found!");
  double Q = surf_post_op.GetSurfaceElectricCharge(idx, E->real());
  if (has_imaginary)
  {
    double Qi = surf_post_op.GetSurfaceElectricCharge(idx, E->imag());
    Q = std::copysign(std::sqrt(Q * Q + Qi * Qi), Q);
  }
  return Q;
}

double PostOperator::GetSurfaceFlux(int idx) const
{
  // Compute the magnetic flux through a surface as Φ_j = ∫ B ⋅ n_j dS. This then yields the
  // inductive coupling to the excitation as M_jk = Φ_j / I_k where I_k is the excitation
  // current. The special coefficient is used to avoid issues evaluating MFEM GridFunctions
  // which are discontinuous at interior boundary elements.
  MFEM_VERIFY(B,
              "Surface inductance not defined, no magnetic flux density solution found!");
  double Phi = surf_post_op.GetSurfaceMagneticFlux(idx, B->real());
  if (has_imaginary)
  {
    double Phii = surf_post_op.GetSurfaceMagneticFlux(idx, B->imag());
    Phi = std::copysign(std::sqrt(Phi * Phi + Phii * Phii), Phi);
  }
  return Phi;
}

void PostOperator::WriteFields(int step, double time) const
{
  // Given the electric field and magnetic flux density, write the fields to disk for
  // visualization.
  bool first_save = (paraview.GetCycle() < 0);
  paraview.SetCycle(step);
  paraview.SetTime(time);
  paraview_bdr.SetCycle(step);
  paraview_bdr.SetTime(time);
  if (first_save)
  {
    mfem::ParMesh &mesh =
        (E) ? *E->ParFESpace()->GetParMesh() : *B->ParFESpace()->GetParMesh();
    mfem::L2_FECollection pwconst_fec(0, mesh.Dimension());
    mfem::ParFiniteElementSpace pwconst_fespace(&mesh, &pwconst_fec);
    mfem::ParGridFunction rank(&pwconst_fespace);
    rank = mesh.GetMyRank() + 1;
    paraview.RegisterField("rank", &rank);
    paraview.Save();
    paraview.DeregisterField("rank");
  }
  else
  {
    paraview.Save();
  }
  paraview_bdr.Save();
}

std::vector<std::complex<double>> PostOperator::ProbeEField() const
{
  MFEM_VERIFY(E, "PostOperator is not configured for electric field probes!");
  return interp_op.ProbeField(*E, has_imaginary);
}

std::vector<std::complex<double>> PostOperator::ProbeBField() const
{
  MFEM_VERIFY(B, "PostOperator is not configured for magnetic flux density probes!");
  return interp_op.ProbeField(*B, has_imaginary);
}

const mfem::ParComplexGridFunction &PostOperator::GetE() const
{
  MFEM_VERIFY(HasE(), "E field has not been configured yet");
  return E.value();
}

const mfem::ParComplexGridFunction &PostOperator::GetB() const
{
  MFEM_VERIFY(HasE(), "B field has not been configured yet");
  return B.value();
}

const mfem::ParGridFunction &PostOperator::GetV() const
{
  MFEM_VERIFY(V.has_value(), "V field has not been configured yet");
  return V.value();
}

const mfem::ParGridFunction &PostOperator::GetA() const
{
  MFEM_VERIFY(A.has_value(), "A field has not been configured yet");
  return A.value();
}

}  // namespace palace
