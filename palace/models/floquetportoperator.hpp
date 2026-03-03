// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_FLOQUET_PORT_OPERATOR_HPP
#define PALACE_MODELS_FLOQUET_PORT_OPERATOR_HPP

#include <complex>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>
#include <mfem.hpp>
#include "fem/gridfunction.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class IoData;
class MaterialOperator;
class MaterialPropertyCoefficient;

namespace config
{
struct FloquetPortData;
}

// Represents a single diffraction order (m, n) with a specific polarization (TE/TM).
struct FloquetMode
{
  int m, n;          // Lattice indices
  bool is_te;        // true = TE (s-pol), false = TM (p-pol)
  mfem::Vector B_mn;  // Transverse wavevector B_mn = m*b1 + n*b2
  mfem::Vector e_pol;  // Polarization unit vector (3D, tangential to port)
  ComplexVector v;     // Fourier projection: v_j = int_Gamma (nxnxN_j).e_p exp(-iB.r) dS
  double gamma_sq;     // gamma_mn^2 = omega^2*mu*eps - |B_mn - k_F|^2 (freq-dependent)
};

// Low-rank complex operator: F*x = sum_k g_k (v_k^H x) v_k.
// Used for the Floquet port boundary matrix.
class LowRankComplexOperator : public ComplexOperator
{
public:
  struct Term
  {
    const ComplexVector *v;   // Projection vector (not owned)
    std::complex<double> g;    // Frequency-dependent scalar weight
  };

private:
  std::vector<Term> terms;
  MPI_Comm comm;
  mutable ComplexVector tmp;

public:
  LowRankComplexOperator(MPI_Comm comm, int n) : ComplexOperator(n), comm(comm), tmp(n)
  {
    tmp.UseDevice(true);
  }

  void AddTerm(const ComplexVector *v, std::complex<double> g) { terms.push_back({v, g}); }

  void Mult(const ComplexVector &x, ComplexVector &y) const override;

  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const override;
};

// Sum of two ComplexOperators: (A + B)*x = A*x + B*x.
class SumComplexOperator : public ComplexOperator
{
private:
  std::unique_ptr<ComplexOperator> op_A, op_B;

public:
  SumComplexOperator(std::unique_ptr<ComplexOperator> &&A,
                     std::unique_ptr<ComplexOperator> &&B)
    : ComplexOperator(A->Height(), A->Width()), op_A(std::move(A)), op_B(std::move(B))
  {
  }

  void Mult(const ComplexVector &x, ComplexVector &y) const override
  {
    op_A->Mult(x, y);
    op_B->AddMult(x, y);
  }

  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a = 1.0) const override
  {
    op_A->AddMult(x, y, a);
    op_B->AddMult(x, y, a);
  }
};

//
// Data for a single Floquet port boundary.
//
class FloquetPortData
{
public:
  int excitation;
  bool active;

  [[nodiscard]] constexpr bool HasExcitation() const { return excitation != 0; }
  const auto &GetAttrList() const { return attr_list; }

  // Access the enumerated modes (for CSV column setup).
  const auto &GetModes() const { return modes; }

  // Check if mode (m, n, is_te) matches the incident mode for this port.
  [[nodiscard]] bool IsIncidentMode(int m, int n, bool is_te) const
  {
    return m == inc_m && n == inc_n && is_te == inc_te;
  }

  FloquetPortData(const config::FloquetPortData &data, const IoData &iodata,
                  const MaterialOperator &mat_op,
                  mfem::ParFiniteElementSpace &nd_fespace);

  // Update propagation constants for given omega. Cached.
  void Initialize(double omega);

  // Get the low-rank boundary operator F(omega) for this port.
  std::unique_ptr<LowRankComplexOperator> GetBoundaryOperator() const;

  // S-parameter for all propagating orders at the current frequency.
  std::map<std::tuple<int, int, bool>, std::complex<double>>
  GetAllSParameters(const GridFunction &E) const;

  // Get the number of propagating orders at current frequency.
  int NumPropagatingOrders() const;

  // Add excitation vector contribution directly to the RHS.
  // For the incident mode: RHS += -2i*omega/mu * v_inc.
  bool AddExcitationVector(double omega, ComplexVector &RHS) const;

private:
  const MaterialOperator &mat_op;
  mfem::Array<int> attr_list;

  // Lattice and reciprocal lattice vectors.
  mfem::Vector a1, a2, b1, b2;

  // Bloch wave vector (from periodic BC config).
  mfem::Vector k_F;

  // Port geometry.
  mfem::Vector port_normal;
  double port_area;

public:
  // Material properties at the port (nondimensional, from adjacent volume element).
  double mu_eps_port;  // mu_r * eps_r (for propagation constant)
  double mu_r_port;    // mu_r (for DtN coefficient and excitation)

private:

  // Diffraction order limits.
  int max_order_m, max_order_n;

  // All Floquet modes with pre-assembled projection vectors.
  std::vector<FloquetMode> modes;

  // Frequency cache.
  double omega0 = 0.0;

  // Incident mode specification.
  int inc_m, inc_n;
  bool inc_te;

  // MPI communicator.
  MPI_Comm comm;

  void ComputeReciprocalLattice();
  void EnumerateOrders();
  void AssembleFourierProjections(mfem::ParFiniteElementSpace &nd_fespace);
};

//
// Manager class for all Floquet port boundaries.
//
class FloquetPortOperator
{
private:
  const MaterialOperator &mat_op;
  std::map<int, FloquetPortData> ports;
  bool suppress_output = false;

public:
  FloquetPortOperator(const IoData &iodata, const MaterialOperator &mat_op,
                      mfem::ParFiniteElementSpace &nd_fespace);

  const FloquetPortData &GetPort(int idx) const { return ports.at(idx); }
  auto begin() const { return ports.begin(); }
  auto end() const { return ports.end(); }
  auto Size() const { return ports.size(); }
  bool Empty() const { return ports.empty(); }

  void SetSuppressOutput(bool suppress) { suppress_output = suppress; }

  // Initialize all ports for the given frequency.
  void Initialize(double omega);

  // Get the combined low-rank boundary operator F(omega) for all Floquet ports.
  // Returns nullptr if no Floquet ports are configured.
  std::unique_ptr<ComplexOperator> GetExtraSystemOperator(double omega);

  // Add Robin boundary mass coefficient (iγ₀/μ) for absorption on Floquet port faces.
  // This provides full-rank absorption (like wave ports) while the low-rank F operator
  // provides corrections for modes with γ ≠ γ₀.
  void AddExtraSystemBdrCoefficients(double omega, MaterialPropertyCoefficient &fbr,
                                     MaterialPropertyCoefficient &fbi);

  // Add excitation vector contributions for the given excitation index.
  bool AddExcitationVector(int excitation_idx, double omega, ComplexVector &RHS);

  // Collect all boundary attributes from all Floquet ports.
  mfem::Array<int> GetAttrList() const;
};

}  // namespace palace

#endif  // PALACE_MODELS_FLOQUET_PORT_OPERATOR_HPP
