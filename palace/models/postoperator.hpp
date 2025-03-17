// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_POST_OPERATOR_HPP
#define PALACE_MODELS_POST_OPERATOR_HPP

#include <complex>
#include <map>
#include <memory>
#include <optional>
#include <type_traits>
#include <vector>
#include <mfem.hpp>
#include "fem/gridfunction.hpp"
#include "fem/interpolator.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/domainpostoperator.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/postoperatorcsv.hpp"
#include "models/surfacepostoperator.hpp"
#include "utils/configfile.hpp"
#include "utils/filesystem.hpp"
#include "utils/units.hpp"

namespace palace
{

class CurlCurlOperator;
class ErrorIndicator;
class IoData;
class LaplaceOperator;
class MaterialOperator;
class SpaceOperator;
class SurfaceCurrentOperator;
class WavePortOperator;

//
// Map solver type to fem_op type
//

// Heleper struct that does not exists in C++17.
// TODO(C++20): use std::type identify or constrain T using FemOp concept.
template <class T>
struct type_identity
{
  using type = T;
};

// TODO(C++20): Replace with this lambda in unevaluated context in decltype of fem_op_t.
template <config::ProblemData::Type solver_t>
struct fem_op_map_type
{
  static constexpr auto map_type()
  {
    if constexpr (solver_t == config::ProblemData::Type::ELECTROSTATIC)
    {
      return type_identity<LaplaceOperator>{};
    }
    else if constexpr (solver_t == config::ProblemData::Type::MAGNETOSTATIC)
    {
      return type_identity<CurlCurlOperator>{};
    }
    else
    {
      return type_identity<SpaceOperator>{};
    }
  }

  using type = typename decltype(map_type())::type;
};

// Helper alias
template <config::ProblemData::Type solver_t>
using fem_op_t = typename fem_op_map_type<solver_t>::type;

// Statically check if solver type has complex fields

template <config::ProblemData::Type solver_t>
constexpr bool HasComplexGridFunction()
{
  return solver_t == config::ProblemData::Type::DRIVEN ||
         solver_t == config::ProblemData::Type::EIGENMODE;
}

// Statically check what fields solver type has
// TODO(C++20): Change these to inline consteval and use with requires.

template <config::ProblemData::Type solver_t>
constexpr bool HasVGridFunction()
{
  return solver_t == config::ProblemData::Type::ELECTROSTATIC;
}

template <config::ProblemData::Type solver_t>
constexpr bool HasAGridFunction()
{
  return solver_t == config::ProblemData::Type::MAGNETOSTATIC;
}

template <config::ProblemData::Type solver_t>
constexpr bool HasEGridFunction()
{
  return solver_t != config::ProblemData::Type::MAGNETOSTATIC;
}

template <config::ProblemData::Type solver_t>
constexpr bool HasBGridFunction()
{
  return solver_t != config::ProblemData::Type::ELECTROSTATIC;
}

//
// A class to handle solution postprocessing.
//
template <config::ProblemData::Type solver_t>
class PostOperator
{
private:
  // Pointer to operator handling discretization and FEM space appropriate to solver,
  // also contains reference to all domains, boundary contions, etc. for measurements.
  // TODO(C++20): Use std::reference_wrapper with incomplete types.
  fem_op_t<solver_t> *fem_op;

  // Copy of unit converter from IOData to scale mesh and dimensionalize measurements
  Units units;

  // Fields: Electric, Magnetic, Scalar Potential, Vector Potential.
  std::unique_ptr<GridFunction> E, B, V, A;

  // Base post-op output directory.
  fs::path post_dir;

  // ******

  // Option to write paraview fields at all
  size_t paraview_delta_post = 0;  // printing rate for paraview (DRIVEN & TRANSIENT)
  size_t paraview_n_post = 0;      // max printing for paraview (OTHER SOLVERS)
  bool write_paraview_fields() const
  {
    return (paraview_delta_post > 0) || (paraview_n_post > 0);
  }

  // Data collection for writing fields to disk for visualization.
  std::optional<mfem::ParaViewDataCollection> paraview, paraview_bdr;

  // Measurements of field solution for paraview files (full domain or surfaces).

  // Poyting Coefficient, Electric Boundary Field (re+im), Magnetic Boundary Field (re+im),
  // Vector Potential Boundary Field, Surface Current (re+im).
  std::unique_ptr<mfem::VectorCoefficient> S, E_sr, E_si, B_sr, B_si, A_s, J_sr, J_si;

  // Electric Energy Density, Magnetic Energy Density, Scalar Potential Boundary Field,
  // Surface Charge (re+im).
  std::unique_ptr<mfem::Coefficient> U_e, U_m, V_s, Q_sr, Q_si;

  // Wave port boundary mode field postprocessing.
  struct WavePortFieldData
  {
    std::unique_ptr<mfem::VectorCoefficient> E0r, E0i;
  };
  std::map<int, WavePortFieldData> port_E0;

  void InitializeParaviewDataCollection();

  // Write to disk the E- and B-fields extracted from the solution vectors. Note that
  // fields are not redimensionalized, to do so one needs to compute: B <= B * (μ₀ H₀), E
  // <= E * (Z₀ H₀), V <= V * (Z₀ H₀ L₀), etc.
  void WriteFields(double time, int step);
  void WriteFieldsFinal(const ErrorIndicator *indicator = nullptr);

  // ******

  // Small friend class for saving csv measurements.
  friend PostOperatorCSV<solver_t>;

  PostOperatorCSV<solver_t> post_op_csv;

  // Measurement from field for csv files.

  DomainPostOperator dom_post_op;           // Energy in bulk
  SurfacePostOperator surf_post_op;         // Dielectric Interface Energy and Flux
  mutable InterpolationOperator interp_op;  // E & B fields: mutates during measure

  // Mini storage for data measurements.

  struct DomainData
  {
    int idx;
    double energy;
    double participation_ratio;
  };

  struct FluxData
  {
    int idx;                   // Surface index
    std::complex<double> Phi;  // Integrated flux
    SurfaceFluxType type;      // Flux type
  };

  struct InterfaceData
  {
    int idx;                      // Interface index
    double energy;                // Surface ELectric Field Energy
    double tandelta;              // Dissipation tangent tan(δ)
    double energy_participation;  // ratio of interface energy / total_energy
    double quality_factor;        // 1 / (energy_paricipation_p * tan δ)
  };

  // For both lumped and wave port
  struct PortPostData
  {
    std::complex<double> P = 0.0;
    std::complex<double> V = 0.0;
    std::complex<double> I = 0.0;
    // Separate R, L, and C branches for current via Z
    std::array<std::complex<double>, 3> I_RLC = {0.0, 0.0, 0.0};

    // S-Parameter
    std::complex<double> S = 0.0;
    double abs_S_ij = 0.0;
    double arg_S_ij = 0.0;

    // Energies (currently only for lumped port).
    double inductor_energy = 0.0;   // E_ind = ∑_j 1/2 L_j I_mj².
    double capacitor_energy = 0.0;  // E_cap = ∑_j 1/2 C_j V_mj².

    // Resitive lumped port (only EIGENMODE)
    double mode_port_kappa = 0.0;
    double quality_factor = mfem::infinity();

    // Inductive lumped port (only EIGENMODE)
    double inductive_energy_participation = 0.0;
  };

  // Results of measurements on field. All computations should return answers in SI units.
  // i.e. Dimenzionalize<units>(value).
  struct Measurement
  {
    // "Pseudo-measurements" — input required during measurement
    std::complex<double> freq = {0.0, 0.0};  // TODO(C++20): requires driven & eigenvalue.
    // Modulation factor for input excitation:
    // I_inc(t) = J(t) I_in for transient
    // I_inc(omega) = I_in for driven, so Jcoeff_excitation = 1.0
    double Jcoeff_excitation = 1.0;  // TODO(C++20): requires transient || driven.

    double domain_E_field_energy_all = 0.0;
    double domain_H_field_energy_all = 0.0;

    std::vector<DomainData> domain_E_field_energy_i;
    std::vector<DomainData> domain_H_field_energy_i;

    double lumped_port_capacitor_energy = 0.0;
    double lumped_port_inductor_energy = 0.0;

    std::map<int, PortPostData> lumped_port_vi;
    std::map<int, PortPostData> wave_port_vi;

    std::vector<std::complex<double>> probe_E_field;
    std::vector<std::complex<double>> probe_B_field;

    std::vector<FluxData> surface_flux_i;
    std::vector<InterfaceData> interface_eps_i;

    // Eigenmode data // TODO(C++20): requires eigenmode
    double eigenmode_Q;
    double error_bkwd;
    double error_abs;
  };

  mutable Measurement measurement_cache;

  // Component measurements to fill the cache which will then be returned and printed.
  void MeasureDomainFieldEnergy() const;
  void MeasureLumpedPorts() const;
  void MeasureWavePorts() const;
  void MeasureLumpedPortsEig() const;  // Depends: DomainFieldEnergy, LumpedPorts
  void MeasureSParameter() const;      // Depends: LumpedPorts, WavePorts
  void MeasureSurfaceFlux() const;
  void MeasureInterfaceEFieldEnergy() const;  // Depends: LumpedPorts
  void MeasureProbes() const;

  void MeasureAllImpl() const
  {
    // Call order matter due to dependent measurements!
    MeasureDomainFieldEnergy();
    MeasureLumpedPorts();
    MeasureWavePorts();
    MeasureLumpedPortsEig();
    MeasureSParameter();
    MeasureSurfaceFlux();
    MeasureInterfaceEFieldEnergy();
    MeasureProbes();
  }

  // ******

  // Populate the grid function solutions for the E- and B-field using the solution vectors
  // on the true dofs. For the real-valued overload, the electric scalar potential can be
  // specified too for electrostatic simulations. The output mesh and fields are
  // nondimensionalized consistently (B ~ E (L₀ ω₀ E₀⁻¹)).
  // TODO(C++20): Switch SFINE to requires.

  template <config::ProblemData::Type U = solver_t>
  auto SetEGridFunction(const ComplexVector &e, bool exchange_face_nbr_data = true)
      -> std::enable_if_t<HasEGridFunction<U>() && HasComplexGridFunction<U>(), void>
  {
    E->Real().SetFromTrueDofs(e.Real());  // Parallel distribute
    E->Imag().SetFromTrueDofs(e.Imag());
    if (exchange_face_nbr_data)
    {
      E->Real().ExchangeFaceNbrData();  // Ready for parallel comm on shared faces
      E->Imag().ExchangeFaceNbrData();
    }
  }

  template <config::ProblemData::Type U = solver_t>
  auto SetEGridFunction(const Vector &e, bool exchange_face_nbr_data = true)
      -> std::enable_if_t<HasEGridFunction<U>() && !HasComplexGridFunction<U>(), void>
  {
    E->Real().SetFromTrueDofs(e);
    if (exchange_face_nbr_data)
    {
      E->Real().ExchangeFaceNbrData();
    }
  }

  template <config::ProblemData::Type U = solver_t>
  auto SetBGridFunction(const ComplexVector &b, bool exchange_face_nbr_data = true)
      -> std::enable_if_t<HasBGridFunction<U>() && HasComplexGridFunction<U>(), void>
  {
    B->Real().SetFromTrueDofs(b.Real());  // Parallel distribute
    B->Imag().SetFromTrueDofs(b.Imag());
    if (exchange_face_nbr_data)
    {
      B->Real().ExchangeFaceNbrData();  // Ready for parallel comm on shared faces
      B->Imag().ExchangeFaceNbrData();
    }
  }

  template <config::ProblemData::Type U = solver_t>
  auto SetBGridFunction(const Vector &b, bool exchange_face_nbr_data = true)
      -> std::enable_if_t<HasBGridFunction<U>() && !HasComplexGridFunction<U>(), void>
  {
    B->Real().SetFromTrueDofs(b);
    if (exchange_face_nbr_data)
    {
      B->Real().ExchangeFaceNbrData();
    }
  }

  template <config::ProblemData::Type U = solver_t>
  auto SetVGridFunction(const Vector &v, bool exchange_face_nbr_data = true)
      -> std::enable_if_t<HasVGridFunction<U>() && !HasComplexGridFunction<U>(), void>
  {
    V->Real().SetFromTrueDofs(v);
    if (exchange_face_nbr_data)
    {
      V->Real().ExchangeFaceNbrData();
    }
  }

  template <config::ProblemData::Type U = solver_t>
  auto SetAGridFunction(const Vector &a, bool exchange_face_nbr_data = true)
      -> std::enable_if_t<HasAGridFunction<U>() && !HasComplexGridFunction<U>(), void>
  {
    A->Real().SetFromTrueDofs(a);
    if (exchange_face_nbr_data)
    {
      A->Real().ExchangeFaceNbrData();
    }
  }

public:
  explicit PostOperator(const IoData &iodata, fem_op_t<solver_t> &fem_op,
                        int nr_expected_measurement_rows = 1);

  // Function that triggers all available measurements; specialize by solver type.
  // TODO: Upgrade SFINE to C++20 concepts to simplify static selection since we can just
  // write `MeasurePrintAll(...) requires (solver_t == Type::A)`.

  template <config::ProblemData::Type U = solver_t>
  auto MeasurePrintAll(int step, const ComplexVector &e, const ComplexVector &b,
                       std::complex<double> omega)
      -> std::enable_if_t<U == config::ProblemData::Type::DRIVEN, double>;

  template <config::ProblemData::Type U = solver_t>
  auto MeasurePrintAll(int step, const ComplexVector &e, const ComplexVector &b,
                       std::complex<double> omega, double error_abs, double error_bkwd,
                       int num_conv)
      -> std::enable_if_t<U == config::ProblemData::Type::EIGENMODE, double>;

  template <config::ProblemData::Type U = solver_t>
  auto MeasurePrintAll(int step, const Vector &v, const Vector &e, int idx)
      -> std::enable_if_t<U == config::ProblemData::Type::ELECTROSTATIC, double>;

  template <config::ProblemData::Type U = solver_t>
  auto MeasurePrintAll(int step, const Vector &a, const Vector &b, int idx)
      -> std::enable_if_t<U == config::ProblemData::Type::MAGNETOSTATIC, double>;

  template <config::ProblemData::Type U = solver_t>
  auto MeasurePrintAll(int step, const Vector &e, const Vector &b, double t, double J_coef)
      -> std::enable_if_t<U == config::ProblemData::Type::TRANSIENT, double>;

  // Write error indicator into paraview file and deregister all paraview objects
  void MeasureFinalize(const ErrorIndicator &indicator);

  // Mini-measurement needed for PROM construction where energy is needed for error
  // indicator, but no other measurement / printing should be done.
  template <config::ProblemData::Type U = solver_t>
  auto MeasureDomainFieldEnergyOnly(const ComplexVector &e, const ComplexVector &b,
                                    bool exchange_face_nbr_data = true)
      -> std::enable_if_t<U == config::ProblemData::Type::DRIVEN, double>;

  // Access grid functions for field solutions.
  // Note these are NOT const — The electrostatics and magnetostatics globallly solver
  // do measurements in the solver and use these as scratch spaces. TODO: Consider moving
  // those global measurements into MeasureFinalize?
  // TODO(C++20): Switch SFINE to requires.
  template <config::ProblemData::Type U = solver_t>
  auto GetEGridFunction() -> std::enable_if_t<HasEGridFunction<U>(), decltype(*E) &>
  {
    return *E;
  }

  template <config::ProblemData::Type U = solver_t>
  auto GetBGridFunction() -> std::enable_if_t<HasBGridFunction<U>(), decltype(*B) &>
  {
    return *B;
  }

  template <config::ProblemData::Type U = solver_t>
  auto GetVGridFunction() -> std::enable_if_t<HasVGridFunction<U>(), decltype(*V) &>
  {
    return *V;
  }

  template <config::ProblemData::Type U = solver_t>
  auto GetAGridFunction() -> std::enable_if_t<HasAGridFunction<U>(), decltype(*A) &>
  {
    return *A;
  }

  // Access to domain postprocessing objects. Use in electro & magnetostatic global
  // post-processing.
  const auto &GetDomainPostOp() const { return dom_post_op; }

  // Expose MPI communicator for Cap & Ind matrix custom processing.
  auto GetComm() const { return fem_op->GetComm(); }
};

}  // namespace palace

#endif  // PALACE_MODELS_POST_OPERATOR_HPP
