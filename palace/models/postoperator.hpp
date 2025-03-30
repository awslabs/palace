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
#include "utils/strongtype.hpp"
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

// Statically map solver (config::ProblemData::Type) to finite element operator

template <config::ProblemData::Type solver_t>
struct fem_op_map_type
{
  using type = SpaceOperator;
};
template <>
struct fem_op_map_type<config::ProblemData::Type::ELECTROSTATIC>
{
  using type = LaplaceOperator;
};
template <>
struct fem_op_map_type<config::ProblemData::Type::MAGNETOSTATIC>
{
  using type = CurlCurlOperator;
};

template <config::ProblemData::Type solver_t>
using fem_op_t = typename fem_op_map_type<solver_t>::type;

// Statically specify if solver uses real or complex fields

template <config::ProblemData::Type solver_t>
constexpr bool HasComplexGridFunction()
{
  return solver_t == config::ProblemData::Type::DRIVEN ||
         solver_t == config::ProblemData::Type::EIGENMODE;
}

// Statically specify what fields a solver uses
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
// A class to handle solution postprocessing for all solvers.
//
template <config::ProblemData::Type solver_t>
class PostOperator
{
private:
  // Pointer to operator handling discretization and FEM space appropriate to solver. It
  // also contains the reference to all domains, boundary conditions, etc. needed for
  // measurement and printing.
  // TODO(C++20): Use std::reference_wrapper with incomplete types.
  fem_op_t<solver_t> *fem_op;

  // Unit converter from IOData to scale mesh and measurements. Lightweight class so it is
  // cheap to copy, rather than keep another reference to IOData.
  Units units;

  // Base post-op output directory.
  fs::path post_dir;

  // Fields: Electric, Magnetic, Scalar Potential, Vector Potential.
  std::unique_ptr<GridFunction> E, B, V, A;

  // ParaView Measure & Print

  // Option to write ParaView fields at all and rate / number of iterations printed.
  size_t paraview_delta_post = 0;  // printing rate for ParaView (DRIVEN & TRANSIENT)
  size_t paraview_n_post = 0;      // max printing for ParaView (OTHER SOLVERS)
  bool write_paraview_fields() const
  {
    return (paraview_delta_post > 0) || (paraview_n_post > 0);
  }

  // ParaView data collection: writing fields to disk for visualization.
  // This is an optional, since ParaViewDataCollection has no default (empty) ctor,
  // and we only want initialize it if write_paraview_fields() is true.
  std::optional<mfem::ParaViewDataCollection> paraview, paraview_bdr;

  // Measurements of field solution for ParaView files (full domain or surfaces).

  // Poynting Coefficient, Electric Boundary Field (re+im), Magnetic Boundary Field (re+im),
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

public:
  // Public functions to switch paraview output to a different sub_folder and reinitialize
  // data collection. Needed in driven solver for multi-excitations and prom output.
  void InitializeParaviewDataCollection(const fs::path &sub_folder_name = "");

  // Secondary overload for the driven solver only, that takes in an excitation index and
  // sets the correct sub_folder_name path for the primary function above.
  template <config::ProblemData::Type U = solver_t>
  auto InitializeParaviewDataCollection(ExcitationIdx ex_idx)
      -> std::enable_if_t<U == config::ProblemData::Type::DRIVEN, void>;

private:
  // Write to disk the E- and B-fields extracted from the solution vectors. Note that
  // fields are not redimensionalized, to do so one needs to compute: B <= B * (μ₀ H₀), E
  // <= E * (Z₀ H₀), V <= V * (Z₀ H₀ L₀), etc.
  void WriteFields(double time, int step);
  void WriteFieldsFinal(const ErrorIndicator *indicator = nullptr);

  // CSV Measure & Print

  // PostOperatorCSV<solver_t> is a class that contains csv tables and printers of
  // measurements. Conceptually, its members could be a part of this class, like the
  // ParaView fields and functions above. It has been separated out for code readability. To
  // achieve this, it is has a pointer back to its "parent" PostOperator class and is a
  // friend class so it can access the the private measurement_cache and references of the
  // system from fem_op.
  friend PostOperatorCSV<solver_t>;

  PostOperatorCSV<solver_t> post_op_csv;

  // Helper classes that actually do some measurements that will be saved to csv files.

  DomainPostOperator dom_post_op;           // Energy in bulk
  SurfacePostOperator surf_post_op;         // Dielectric Interface Energy and Flux
  mutable InterpolationOperator interp_op;  // E & B fields: mutates during measure

  // Mini storage structs for data measurements.

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
    SurfaceFluxType type;
  };

  struct InterfaceData
  {
    int idx;                      // Interface index
    double energy;                // Surface Electric Field Energy
    double tandelta;              // Dissipation tangent tan(δ)
    double energy_participation;  // ratio of interface energy / total_energy
    double quality_factor;        // 1 / (energy_participation * tan δ)
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

    // Energies (currently only for lumped port)
    double inductor_energy = 0.0;   // E_ind = ∑_j 1/2 L_j I_mj².
    double capacitor_energy = 0.0;  // E_cap = ∑_j 1/2 C_j V_mj².

    // Resistive lumped port (only EIGENMODE)
    double mode_port_kappa = 0.0;
    double quality_factor = mfem::infinity();

    // Inductive lumped port (only EIGENMODE)
    double inductive_energy_participation = 0.0;
  };

  // Results of measurements on fields. All values here should be in SI units (i.e. we call
  // Dimenzionalize<units>(value) on the results before storage here). Not all measurements
  // are sensible to define for all solvers, but it is a hassle to statically restrict them.
  // Since this class is private we will just leave them defaulted.
  struct Measurement
  {
    // "Pseudo-measurements": input required during measurement or data which is stored here
    // in order to pass it along to the printers.

    ExcitationIdx ex_idx = ExcitationIdx(0);  // driven

    std::complex<double> freq = {0.0, 0.0};  // driven || eigenvalue.

    // Modulation factor for input excitation:
    // - I_inc(t) = J(t) I_in, for transient
    // - I_inc(omega) = I_in, for driven so that Jcoeff_excitation = 1.0
    double Jcoeff_excitation = 1.0;  // transient || driven

    // Eigenmode data including error from solver
    double eigenmode_Q = 0.0;
    double error_bkwd = 0.0;
    double error_abs = 0.0;

    // "Actual measurements"

    double domain_E_field_energy_all = 0.0;
    double domain_H_field_energy_all = 0.0;

    std::vector<DomainData> domain_E_field_energy_i;
    std::vector<DomainData> domain_H_field_energy_i;

    double lumped_port_capacitor_energy = 0.0;
    double lumped_port_inductor_energy = 0.0;

    std::map<int, PortPostData> lumped_port_vi;
    std::map<int, PortPostData> wave_port_vi;

    // Probe data is ordered as [Fx1, Fy1, Fz1, Fx2, Fy2, Fz2, ...].
    // TODO: Replace with proper matrix: mdspan (C++23) / Eigen.
    std::vector<std::complex<double>> probe_E_field;
    std::vector<std::complex<double>> probe_B_field;

    std::vector<FluxData> surface_flux_i;
    std::vector<InterfaceData> interface_eps_i;
  };

  mutable Measurement measurement_cache;

  // Individual measurements to fill the cache/workspace. Measurements functions are not
  // constrained by solver type in the signature since they are private member functions.
  // They dispatch on solver type within the function itself using `if constexpr`, and do
  // nothing if the measurement is not solver appropriate.
  void MeasureDomainFieldEnergy() const;
  void MeasureLumpedPorts() const;
  void MeasureWavePorts() const;
  void MeasureLumpedPortsEig() const;  // Depends: DomainFieldEnergy, LumpedPorts
  void MeasureSParameter() const;      // Depends: LumpedPorts, WavePorts
  void MeasureSurfaceFlux() const;
  void MeasureInterfaceEFieldEnergy() const;  // Depends: LumpedPorts
  void MeasureProbes() const;

  // Helper function called by all solvers. Has to ensure correct call order to deal with
  // dependant measurements.
  void MeasureAllImpl() const
  {
    MeasureDomainFieldEnergy();
    MeasureLumpedPorts();
    MeasureWavePorts();
    MeasureLumpedPortsEig();
    MeasureSParameter();
    MeasureSurfaceFlux();
    MeasureInterfaceEFieldEnergy();
    MeasureProbes();
  }

  // Setting grid functions.
  //
  // Populate the grid function solutions for the E- and B-field using the solution vectors
  // on the true dofs. For the real-valued overload, the electric scalar potential can be
  // specified too for electrostatic simulations. The output mesh and fields are
  // non-dimensionalized consistently (B ~ E (L₀ ω₀ E₀⁻¹)).
  //
  // These functions are private helper functions. We want to enforce that a caller passes
  // the appropriate ones as part of the MeasureAndPrintAll interface, rather than do a
  // runtime check to see that they have been set.
  //
  // TODO(C++20): Switch SFINAE to requires.

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
  // Ctor.
  // - nr_expected_measurement_rows is helper int to reserve space in csv printer tables.
  explicit PostOperator(const IoData &iodata, fem_op_t<solver_t> &fem_op,
                        int nr_expected_measurement_rows = 1);

  // MeasureAndPrintAll is the primary public interface of this class. It is specialized by
  // solver type, since each solver has different fields and extra data required. These
  // functions all:
  // 1) Set the GridFunctions which have to be passed as part of the call.
  // 2) Perform all measurements and populate measurement_cache with temporary results. This
  //    cache structure exists since measurements have dependencies; we may use some
  //    measurement results in later measurements.
  // 3) Pass the measurement cache to the csv printer which will add the appropriate
  //    rows/cols to the csv tables and print to file.
  // 4) Trigger ParaView field computation and save.
  //
  // The functions return the total domain energy which is the only thing needed in the
  // solver to normalize the error indicator. If more measurements were needed by the solver
  // loop, we could imagine passing a small struct (like Measurement above or some sub-set
  // therefore).
  //
  // The measure functions will also do logging of (some) measurements to stdout.
  //
  // TODO(C++20): Upgrade SFINAE to C++20 concepts to simplify static selection since we can
  // just write `MeasureAndPrintAll(...) requires (solver_t == Type::A)` without extra
  // template.

  template <config::ProblemData::Type U = solver_t>
  auto MeasureAndPrintAll(ExcitationIdx ex_idx, int step, const ComplexVector &e,
                          const ComplexVector &b, std::complex<double> omega)
      -> std::enable_if_t<U == config::ProblemData::Type::DRIVEN, double>;

  template <config::ProblemData::Type U = solver_t>
  auto MeasureAndPrintAll(int step, const ComplexVector &e, const ComplexVector &b,
                          std::complex<double> omega, double error_abs, double error_bkwd,
                          int num_conv)
      -> std::enable_if_t<U == config::ProblemData::Type::EIGENMODE, double>;

  template <config::ProblemData::Type U = solver_t>
  auto MeasureAndPrintAll(int step, const Vector &v, const Vector &e, int idx)
      -> std::enable_if_t<U == config::ProblemData::Type::ELECTROSTATIC, double>;

  template <config::ProblemData::Type U = solver_t>
  auto MeasureAndPrintAll(int step, const Vector &a, const Vector &b, int idx)
      -> std::enable_if_t<U == config::ProblemData::Type::MAGNETOSTATIC, double>;

  template <config::ProblemData::Type U = solver_t>
  auto MeasureAndPrintAll(int step, const Vector &e, const Vector &b, double t,
                          double J_coef)
      -> std::enable_if_t<U == config::ProblemData::Type::TRANSIENT, double>;

  // Write error indicator into ParaView file and print summary statistics to csv. Should be
  // called once at the end of the solver loop.
  void MeasureFinalize(const ErrorIndicator &indicator);

  // Measurement of the domain energy without printing. This is needed during the driven
  // simulation with PROM. There samples are taken and we need the total domain energy for
  // the error indicator, but no other measurement / printing should be done.
  //
  // TODO(C++20): SFINAE to requires.
  template <config::ProblemData::Type U = solver_t>
  auto MeasureDomainFieldEnergyOnly(
      const ComplexVector &e, const ComplexVector &b, bool exchange_face_nbr_data = true,
      std::optional<std::pair<int, double>> debug_print_paraview_opt = std::nullopt)
      -> std::enable_if_t<U == config::ProblemData::Type::DRIVEN, double>;

  // Access grid functions for field solutions. Note that these are NOT const functions. The
  // electrostatics / magnetostatics solver do measurements of the capacitance/ inductance
  // matrix globally at the end of all solves. This is done in the solver class, but uses
  // the GridFunctions in this (PostOp) class as already allocated scratch workspace.
  //
  // Future: Consider moving those cap/ind measurements into this class and MeasureFinalize?
  // Would need to store vector of V,A.
  //
  // TODO(C++20): Switch SFINAE to requires.
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

  // Access to domain postprocessing objects. Use in electrostatic & magnetostatic matrix
  // measurement (see above).
  const auto &GetDomainPostOp() const { return dom_post_op; }

  // Expose MPI communicator from fem_op for electrostatic & magnetostatic matrix processing
  // (see above).
  auto GetComm() const { return fem_op->GetComm(); }
};

}  // namespace palace

#endif  // PALACE_MODELS_POST_OPERATOR_HPP
