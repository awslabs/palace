// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "phasereference.hpp"

#include <cmath>
#include "fem/fespace.hpp"
#include "utils/geodata.hpp"

namespace palace
{

namespace
{

// Generic coefficients, the fractional parts of the multiples of the golden ratio, which
// give the reference fields no symmetry of their own.
double Coefficient(int k)
{
  return std::fmod(k * 0.6180339887498949, 1.0);
}

}  // namespace

void PhaseReference::EvalField(int k, const double *xi, int dim, double *w)
{
  if (k == 0)
  {
    // One term Πₖ∈ₜ ξₖ for every subset t of the coordinates.
    const int num_terms = 1 << dim;
    for (int a = 0; a < dim; a++)
    {
      double s = 0.0;
      for (int t = 0; t < num_terms; t++)
      {
        double m = 0.5 + Coefficient(1 + a * num_terms + t);
        for (int d = 0; d < dim; d++)
        {
          if (t & (1 << d))
          {
            m *= xi[d];
          }
        }
        s += m;
      }
      w[a] = s;
    }
  }
  else
  {
    for (int a = 0; a < dim; a++)
    {
      double s = 1.0;
      for (int d = 0; d < dim; d++)
      {
        const double kappa = 2.0 + 2.0 * Coefficient(100 + a * dim + d);
        s *= std::cos(kappa * xi[d]) + std::sin(kappa * xi[d]);
      }
      w[a] = s;
    }
  }
}

PhaseReference::PhaseReference(FiniteElementSpace &nd_fespace) : fespace(nd_fespace)
{
  MFEM_VERIFY(!fespace.IsBroken(),
              "The eigenvector phase reference needs a space which is not broken!");
  const auto &mesh = fespace.GetParMesh();
  const int dim = mesh.SpaceDimension();
  mfem::Vector bb_min, bb_max;
  mesh::GetAxisAlignedBoundingBox(mesh, bb_min, bb_max);
  center.SetSize(dim);
  half.SetSize(dim);
  for (int d = 0; d < dim; d++)
  {
    center(d) = 0.5 * (bb_min(d) + bb_max(d));
    half(d) = 0.5 * (bb_max(d) - bb_min(d));
    if (!(half(d) > 0.0))
    {
      half(d) = 1.0;
    }
  }
  GetLinearForm(0);
}

const Vector &PhaseReference::GetLinearForm(int k) const
{
  if (!assembled[k])
  {
    const int dim = center.Size();
    mfem::VectorFunctionCoefficient wk(
        dim,
        [this, k, dim](const mfem::Vector &x, mfem::Vector &w)
        {
          double xi[3];
          for (int d = 0; d < dim; d++)
          {
            xi[d] = (x(d) - center(d)) / half(d);
          }
          EvalField(k, xi, dim, w.GetData());
        });
    mfem::LinearForm lf(&fespace.Get());
    lf.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(wk));
    lf.UseFastAssembly(false);
    lf.UseDevice(false);
    lf.Assemble();
    lf.UseDevice(true);
    b[k].SetSize(fespace.GetTrueVSize());
    b[k].UseDevice(true);
    fespace.GetProlongationMatrix()->MultTranspose(lf, b[k]);
    b_norm[k] = linalg::Norml2(fespace.GetComm(), b[k]);
    assembled[k] = true;
  }
  return b[k];
}

int PhaseReference::Normalize(ComplexVector &x) const
{
  // The first reference integral, the squared norm and the sum of the degrees of freedom in
  // one reduction, and the next reference integrals only when needed.
  const MPI_Comm comm = fespace.GetComm();
  const Vector &b0 = GetLinearForm(0);
  const std::complex<double> sum = linalg::LocalSum(x);
  double sums[5] = {linalg::LocalDot(b0, x.Real()), linalg::LocalDot(b0, x.Imag()),
                    linalg::LocalDot(x.Real(), x.Real()) +
                        linalg::LocalDot(x.Imag(), x.Imag()),
                    sum.real(), sum.imag()};
  Mpi::GlobalSum(5, sums, comm);
  const double x_norm = std::sqrt(sums[2]);
  if (!(x_norm > 0.0))
  {
    return -1;
  }
  std::complex<double> z(sums[0], sums[1]);
  int ref = 0;
  while (!Accept(z, b_norm[ref], x_norm))
  {
    if (++ref == num_fields)
    {
      z = {sums[3], sums[4]};
      if (!(std::abs(z) > 0.0))
      {
        return -1;
      }
      break;
    }
    const Vector &bk = GetLinearForm(ref);
    double f[2] = {linalg::LocalDot(bk, x.Real()), linalg::LocalDot(bk, x.Imag())};
    Mpi::GlobalSum(2, f, comm);
    z = {f[0], f[1]};
  }
  x *= std::conj(z) / std::abs(z);
  return ref;
}

}  // namespace palace
