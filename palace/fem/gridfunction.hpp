// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_GRIDFUNCTION_HPP
#define PALACE_FEM_GRIDFUNCTION_HPP

#include <mfem.hpp>

namespace palace
{

class FiniteElementSpace;

//
// A real- or complex-valued grid function represented as two real grid functions, one for
// each component. This unifies mfem::ParGridFunction and mfem::ParComplexGridFunction, and
// replaces the latter due to some issues observed with memory aliasing on GPUs.
//
class GridFunction
{
private:
  mfem::ParGridFunction gfr, gfi;

public:
  GridFunction(mfem::ParFiniteElementSpace &fespace, bool complex = false);
  GridFunction(FiniteElementSpace &fespace, bool complex = false);

  // Get access to the real and imaginary grid function parts.
  const mfem::ParGridFunction &Real() const { return gfr; }
  mfem::ParGridFunction &Real() { return gfr; }
  const mfem::ParGridFunction &Imag() const
  {
    MFEM_ASSERT(HasImag(),
                "Invalid access of imaginary part of a real-valued GridFunction object!");
    return gfi;
  }
  mfem::ParGridFunction &Imag()
  {
    MFEM_ASSERT(HasImag(),
                "Invalid access of imaginary part of a real-valued GridFunction object!");
    return gfi;
  }

  // Check if the grid function is suited for storing complex-valued fields.
  bool HasImag() const { return (gfi.ParFESpace() != nullptr); }

  // Get access to the underlying finite element space (match MFEM interface).
  mfem::FiniteElementSpace *FESpace() { return gfr.FESpace(); }
  const mfem::FiniteElementSpace *FESpace() const { return gfr.FESpace(); }
  mfem::ParFiniteElementSpace *ParFESpace() { return gfr.ParFESpace(); }
  const mfem::ParFiniteElementSpace *ParFESpace() const { return gfr.ParFESpace(); }

  // Dimension of the vector field represented by the grid function.
  int VectorDim() const { return gfr.VectorDim(); }

  // Set all entries equal to s.
  GridFunction &operator=(std::complex<double> s);
  GridFunction &operator=(double s)
  {
    *this = std::complex<double>(s, 0.0);
    return *this;
  }

  // Scale all entries by s.
  GridFunction &operator*=(double s);

  // Transform for space update (for example on mesh change).
  void Update();

  // Get the associated MPI communicator.
  MPI_Comm GetComm() const { return ParFESpace()->GetComm(); }
};

}  // namespace palace

#endif  // PALACE_FEM_GRIDFUNCTION_HPP
