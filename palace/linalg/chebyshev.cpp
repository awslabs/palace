// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "chebyshev.hpp"

#include <cmath>
#include <mfem/general/forall.hpp>

namespace palace
{

namespace
{

bool HasFiniteNonzeroDiagonal(MPI_Comm comm, const Vector &real,
                              const Vector *imag = nullptr)
{
  int valid = 1;
  if (real.Size() > 0)
  {
    const auto *R = real.HostRead();
    const auto *I = imag ? imag->HostRead() : nullptr;
    for (int i = 0; i < real.Size(); i++)
    {
      if (!(std::isfinite(R[i]) && R[i] != 0.0) ||
          (I && !(std::isfinite(I[i]) && I[i] == 0.0)))
      {
        valid = 0;
        break;
      }
    }
  }
  Mpi::GlobalMin(1, &valid, comm);
  return valid;
}

double LambdaMax(MPI_Comm comm, const Operator &A, const Vector &dinv, Vector &x0)
{
  // D⁻¹A is generally not Hermitian. With Σ = sgn(D), it is similar to ΣS for the
  // Hermitian S = |D|⁻¹ᐟ²A|D|⁻¹ᐟ², so ρ(D⁻¹A) ≤ ‖S‖₂, with equality when Σ commutes with S:
  // for a diagonal of uniform sign (A definite of either sign), and for a negative definite
  // A with essential DOFs eliminated using a unit diagonal.
  Vector dinv_sqrt(dinv);
  dinv_sqrt.Abs();
  linalg::Sqrt(dinv_sqrt);
  DiagonalOperator DinvSqrt(dinv_sqrt);
  ProductOperator ADinvSqrt(A, DinvSqrt);
  ProductOperator S(DinvSqrt, ADinvSqrt);
  return linalg::SpectralNorm(comm, S, true, &x0);
}

double GetLambdaMax(MPI_Comm comm, const Operator &A, const Vector &dinv, Vector &x0)
{
  MFEM_VERIFY(HasFiniteNonzeroDiagonal(comm, dinv),
              "Chebyshev smoother spectral estimation requires a finite, nonzero operator "
              "diagonal!");
  return LambdaMax(comm, A, dinv, x0);
}

double GetLambdaMax(MPI_Comm comm, const ComplexOperator &A, const ComplexVector &dinv,
                    ComplexVector &x0)
{
  if (A.IsReal())
  {
    MFEM_VERIFY(HasFiniteNonzeroDiagonal(comm, dinv.Real(), &dinv.Imag()),
                "Chebyshev smoother spectral estimation requires a finite, nonzero real "
                "operator diagonal!");

    // Estimate in real arithmetic, which halves the cost, starting from and returning the
    // real part of x0.
    const int n = A.Height();
    Vector x0r;
    if (x0.Size() == n)
    {
      x0r.SetSize(n);
      x0r.UseDevice(true);
      x0r = x0.Real();
    }
    const double lambda_max = LambdaMax(comm, *A.Real(), dinv.Real(), x0r);
    if (x0r.Size() == n)
    {
      x0.SetSize(n);
      x0.UseDevice(true);
      x0.Real() = x0r;
      x0.Imag() = 0.0;
    }
    return lambda_max;
  }
  ComplexDiagonalOperator Dinv(dinv);
  ComplexProductOperator DinvA(Dinv, A);
  return linalg::SpectralNorm(comm, DinvA, false, &x0);
}

// EXPERIMENT (PR 1201 review): coefficients beta_i of the optimized 4th-kind Chebyshev
// iteration of Lottes (2023), eq. (19), x_i = x_{i-1} + beta_i z_i, for orders 1 to 16,
// from the algorithm of its Appendix B. Orders above 16 use the plain 4th kind.
constexpr int max_opt_order = 16;
constexpr double opt_beta[max_opt_order][max_opt_order] = {
    {1.125000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {1.023872875703131, 1.264089053710855, 0.000000000000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {1.008425447820283, 1.088678392087309, 1.337531259096186, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {1.003913104272854, 1.040358111885930, 1.148634985462549, 1.382688692410005,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {1.002129301461647, 1.021737115492609, 1.078724331926030, 1.198100652926630,
     1.413225427916821, 0.000000000000000, 0.000000000000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {1.001285172559402, 1.013042930352334, 1.046782151241133, 1.116164894196755,
     1.238290202184445, 1.435242971067448, 0.000000000000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {1.000834643979124, 1.008439494301229, 1.030087077687138, 1.074083840920038,
     1.150361867073664, 1.271164740461399, 1.451866586493640, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {1.000572466311977, 1.005774276624156, 1.020501879229414, 1.050198034445655,
     1.101155729849411, 1.180860428068566, 1.298385853825760, 1.464860731510997,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {1.000409600728328, 1.004124395061066, 1.014602121482666, 1.035611136266718,
     1.071399725291942, 1.126882737109629, 1.207852191407293, 1.321219307167467,
     1.475296428206996, 0.000000000000000, 0.000000000000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {1.000303122296529, 1.003048406607969, 1.010770227153876, 1.026190115976401,
     1.052317249337553, 1.092557432075499, 1.150833766639722, 1.231722508708947,
     1.340608020244600, 1.483861244070113, 0.000000000000000, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {1.000230585952094, 1.002316750240285, 1.008172453963049, 1.019829865663422,
     1.039502102353247, 1.069650427005415, 1.113057542957428, 1.172908762755644,
     1.252883005767926, 1.357255799195193, 1.491016725641392, 0.000000000000000,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {1.000179472008284, 1.001801891396195, 1.006348619073076, 1.015378645663060,
     1.030569428307604, 1.053760196939435, 1.086998625920729, 1.132591830979134,
     1.193162733581726, 1.271712936751106, 1.371693379697995, 1.497084185755624,
     0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {1.000142419215596, 1.001429069326297, 1.005030289862981, 1.012169105184954,
     1.024148743427926, 1.042381588808204, 1.068420081287009, 1.103990109367600,
     1.151027482426457, 1.211718119101251, 1.288542648651281, 1.384326193809916,
     1.502294187573689, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {1.000114905382619, 1.001152463769146, 1.004053573332648, 1.009795905731531,
     1.019413004729946, 1.034014250354368, 1.054805996066293, 1.083114203018131,
     1.120408916608924, 1.168330956554461, 1.228721222882388, 1.303653057078175,
     1.395468140536785, 1.506816462095838, 0.000000000000000, 0.000000000000000},
    {1.000094047507526, 1.000942916963435, 1.003314490564448, 1.008002948338161,
     1.015842362591404, 1.027720833177050, 1.044595354228315, 1.067507612061253,
     1.097600925458897, 1.136138553661573, 1.184523614262362, 1.244320873044759,
     1.317280690833927, 1.405365438935602, 1.510778725018453, 0.000000000000000},
    {1.000077948281792, 1.000781268472533, 1.002744879744013, 1.006622910170151,
     1.013098588369713, 1.022894483293370, 1.036783214099839, 1.055598757198966,
     1.080248484055602, 1.111726071314977, 1.151125434310726, 1.199655846149737,
     1.258658417449458, 1.329624126566652, 1.414213606955761, 1.514278917303463}};

inline double OptBeta(int order, int i)
{
  return (order <= max_opt_order) ? opt_beta[order - 1][i] : 1.0;
}

template <bool Transpose = false>
inline void ApplyOp(const Operator &A, const Vector &x, Vector &y)
{
  A.Mult(x, y);
}

template <bool Transpose = false>
inline void ApplyOp(const ComplexOperator &A, const ComplexVector &x, ComplexVector &y)
{
  if constexpr (!Transpose)
  {
    A.Mult(x, y);
  }
  else
  {
    A.MultHermitianTranspose(x, y);
  }
}

template <bool Transpose = false>
inline void ApplyOp(const Operator &A, const Vector &x, Vector &y, const double a)
{
  A.AddMult(x, y, a);
}

template <bool Transpose = false>
inline void ApplyOp(const ComplexOperator &A, const ComplexVector &x, ComplexVector &y,
                    const double a)
{
  if constexpr (!Transpose)
  {
    A.AddMult(x, y, a);
  }
  else
  {
    A.AddMultHermitianTranspose(x, y, a);
  }
}

template <bool Transpose = false>
inline void ApplyOrder0(double sr, const Vector &dinv, const Vector &r, Vector &d)
{
  const bool use_dev = dinv.UseDevice() || r.UseDevice() || d.UseDevice();
  const int N = d.Size();
  const auto *DI = dinv.Read(use_dev);
  const auto *R = r.Read(use_dev);
  auto *D = d.Write(use_dev);
  mfem::forall_switch(use_dev, N,
                      [=] MFEM_HOST_DEVICE(int i) { D[i] = sr * DI[i] * R[i]; });
}

template <bool Transpose = false>
inline void ApplyOrder0(const double sr, const ComplexVector &dinv, const ComplexVector &r,
                        ComplexVector &d)
{
  const bool use_dev = dinv.UseDevice() || r.UseDevice() || d.UseDevice();
  const int N = dinv.Size();
  const auto *DIR = dinv.Real().Read(use_dev);
  const auto *DII = dinv.Imag().Read(use_dev);
  const auto *RR = r.Real().Read(use_dev);
  const auto *RI = r.Imag().Read(use_dev);
  auto *DR = d.Real().Write(use_dev);
  auto *DI = d.Imag().Write(use_dev);
  if constexpr (!Transpose)
  {
    mfem::forall_switch(use_dev, N,
                        [=] MFEM_HOST_DEVICE(int i)
                        {
                          DR[i] = sr * (DIR[i] * RR[i] - DII[i] * RI[i]);
                          DI[i] = sr * (DII[i] * RR[i] + DIR[i] * RI[i]);
                        });
  }
  else
  {
    mfem::forall_switch(use_dev, N,
                        [=] MFEM_HOST_DEVICE(int i)
                        {
                          DR[i] = sr * (DIR[i] * RR[i] + DII[i] * RI[i]);
                          DI[i] = sr * (-DII[i] * RR[i] + DIR[i] * RI[i]);
                        });
  }
}

template <bool Transpose = false>
inline void ApplyOrderK(const double sd, const double sr, const Vector &dinv,
                        const Vector &r, Vector &d)
{
  const bool use_dev = dinv.UseDevice() || r.UseDevice() || d.UseDevice();
  const int N = dinv.Size();
  const auto *DI = dinv.Read(use_dev);
  const auto *R = r.Read(use_dev);
  auto *D = d.ReadWrite(use_dev);
  mfem::forall_switch(use_dev, N, [=] MFEM_HOST_DEVICE(int i)
                      { D[i] = sd * D[i] + sr * DI[i] * R[i]; });
}

template <bool Transpose = false>
inline void ApplyOrderK(const double sd, const double sr, const ComplexVector &dinv,
                        const ComplexVector &r, ComplexVector &d)
{
  const bool use_dev = dinv.UseDevice() || r.UseDevice() || d.UseDevice();
  const int N = dinv.Size();
  const auto *DIR = dinv.Real().Read(use_dev);
  const auto *DII = dinv.Imag().Read(use_dev);
  const auto *RR = r.Real().Read(use_dev);
  const auto *RI = r.Imag().Read(use_dev);
  auto *DR = d.Real().ReadWrite(use_dev);
  auto *DI = d.Imag().ReadWrite(use_dev);
  if constexpr (!Transpose)
  {
    mfem::forall_switch(use_dev, N,
                        [=] MFEM_HOST_DEVICE(int i)
                        {
                          DR[i] = sd * DR[i] + sr * (DIR[i] * RR[i] - DII[i] * RI[i]);
                          DI[i] = sd * DI[i] + sr * (DII[i] * RR[i] + DIR[i] * RI[i]);
                        });
  }
  else
  {
    mfem::forall_switch(use_dev, N,
                        [=] MFEM_HOST_DEVICE(int i)
                        {
                          DR[i] = sd * DR[i] + sr * (DIR[i] * RR[i] + DII[i] * RI[i]);
                          DI[i] = sd * DI[i] + sr * (-DII[i] * RR[i] + DIR[i] * RI[i]);
                        });
  }
}

}  // namespace

template <OperatorType OperType>
ChebyshevSmoother<OperType>::ChebyshevSmoother(MPI_Comm comm, int smooth_it, int poly_order,
                                               double sf_max)
  : Solver<OperType>(), comm(comm), pc_it(smooth_it), order(poly_order), A(nullptr),
    lambda_max(0.0), sf_max(sf_max)
{
  MFEM_VERIFY(order > 0, "Polynomial order for Chebyshev smoothing must be positive!");
}

template <OperatorType OperType>
void ChebyshevSmoother<OperType>::SetOperator(const OperType &op)
{
  A = &op;
  d.SetSize(op.Height());
  dinv.SetSize(op.Height());
  d.UseDevice(true);
  dinv.UseDevice(true);
  op.AssembleDiagonal(dinv);
  dinv.Reciprocal();

  // Set up Chebyshev coefficients using the computed maximum eigenvalue estimate. See
  // mfem::OperatorChebyshevSmoother or Adams et al. (2003). The estimate starts from the
  // dominant vector of the previous operator, which is usually a nearby one (for example,
  // the next frequency of a driven sweep).
  lambda_max = sf_max * GetLambdaMax(comm, *A, dinv, lambda_vec);
  MFEM_VERIFY(std::isfinite(lambda_max) && lambda_max > 0.0,
              "Encountered invalid maximum eigenvalue in Chebyshev smoother!");

  this->height = op.Height();
  this->width = op.Width();
}

template <OperatorType OperType>
void ChebyshevSmoother<OperType>::Mult2(const VecType &x, VecType &y, VecType &r) const
{
  // Apply smoother: y = y + p(A) (x - A y) .
  for (int it = 0; it < pc_it; it++)
  {
    if (this->initial_guess || it > 0)
    {
      ApplyOp(*A, y, r);
      linalg::AXPBY(1.0, x, -1.0, r);
    }
    else
    {
      r = x;
      y = 0.0;
    }

    // 4th-kind Chebyshev smoother, from Phillips and Fischer or Lottes (with k -> k + 1
    // shift due to 1-based indexing).
    ApplyOrder0(4.0 / (3.0 * lambda_max), dinv, r, d);
    for (int k = 1; k < order; k++)
    {
      linalg::AXPY(OptBeta(order, k - 1), d, y);
      ApplyOp(*A, d, r, -1.0);
      const double sd = (2.0 * k - 1.0) / (2.0 * k + 3.0);
      const double sr = (8.0 * k + 4.0) / ((2.0 * k + 3.0) * lambda_max);
      ApplyOrderK(sd, sr, dinv, r, d);
    }
    linalg::AXPY(OptBeta(order, order - 1), d, y);
  }
}

template <OperatorType OperType>
ChebyshevSmoother1stKind<OperType>::ChebyshevSmoother1stKind(MPI_Comm comm, int smooth_it,
                                                             int poly_order, double sf_max,
                                                             double sf_min)
  : Solver<OperType>(), comm(comm), pc_it(smooth_it), order(poly_order), A(nullptr),
    theta(0.0), sf_max(sf_max), sf_min(sf_min)
{
  MFEM_VERIFY(order > 0, "Polynomial order for Chebyshev smoothing must be positive!");
}

template <OperatorType OperType>
void ChebyshevSmoother1stKind<OperType>::SetOperator(const OperType &op)
{
  A = &op;
  d.SetSize(op.Height());
  dinv.SetSize(op.Height());
  d.UseDevice(true);
  dinv.UseDevice(true);
  op.AssembleDiagonal(dinv);
  dinv.Reciprocal();

  // Set up Chebyshev coefficients using the computed maximum eigenvalue estimate. The
  // optimized estimate of lambda_min comes from (2.24) of Phillips and Fischer (2022). See
  // ChebyshevSmoother::SetOperator for the starting vector.
  if (sf_min <= 0.0)
  {
    sf_min = 1.69 / (std::pow(order, 1.68) + 2.11 * order + 1.98);
  }
  const double lambda_max = sf_max * GetLambdaMax(comm, *A, dinv, lambda_vec);
  MFEM_VERIFY(std::isfinite(lambda_max) && lambda_max > 0.0,
              "Encountered invalid maximum eigenvalue in Chebyshev smoother!");
  const double lambda_min = sf_min * lambda_max;
  theta = 0.5 * (lambda_max + lambda_min);
  delta = 0.5 * (lambda_max - lambda_min);

  this->height = op.Height();
  this->width = op.Width();
}

template <OperatorType OperType>
void ChebyshevSmoother1stKind<OperType>::Mult2(const VecType &x, VecType &y,
                                               VecType &r) const
{
  // Apply smoother: y = y + p(A) (x - A y) .
  for (int it = 0; it < pc_it; it++)
  {
    if (this->initial_guess || it > 0)
    {
      ApplyOp(*A, y, r);
      linalg::AXPBY(1.0, x, -1.0, r);
    }
    else
    {
      r = x;
      y = 0.0;
    }

    // 1th-kind Chebyshev smoother, from Phillips and Fischer or Adams.
    ApplyOrder0(1.0 / theta, dinv, r, d);
    double rhop = delta / theta;
    for (int k = 1; k < order; k++)
    {
      y += d;
      ApplyOp(*A, d, r, -1.0);
      const double rho = 1.0 / (2.0 * theta / delta - rhop);
      const double sd = rho * rhop;
      const double sr = 2.0 * rho / delta;
      ApplyOrderK(sd, sr, dinv, r, d);
      rhop = rho;
    }
    y += d;
  }
}

template class ChebyshevSmoother<Operator>;
template class ChebyshevSmoother<ComplexOperator>;

template class ChebyshevSmoother1stKind<Operator>;
template class ChebyshevSmoother1stKind<ComplexOperator>;

}  // namespace palace
