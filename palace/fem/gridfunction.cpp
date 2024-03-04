// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "gridfunction.hpp"

#include "fem/fespace.hpp"

namespace palace
{

GridFunction::GridFunction(mfem::ParFiniteElementSpace &fespace, bool complex)
  : gfr(&fespace)
{
  if (complex)
  {
    gfi.SetSpace(&fespace);
  }
}

GridFunction::GridFunction(FiniteElementSpace &fespace, bool complex)
  : GridFunction(fespace.Get(), complex)
{
}

GridFunction &GridFunction::operator=(std::complex<double> s)
{
  Real() = s.real();
  if (HasImag())
  {
    Imag() = s.imag();
  }
  else
  {
    MFEM_ASSERT(
        s.imag() == 0.0,
        "Cannot assign complex scalar to a non-complex-valued GridFunction object!");
  }
  return *this;
}

GridFunction &GridFunction::operator*=(double s)
{
  Real() *= s;
  if (HasImag())
  {
    Imag() *= s;
  }
  return *this;
}

void GridFunction::Update()
{
  Real().Update();
  if (HasImag())
  {
    Imag().Update();
  }
}

}  // namespace palace
