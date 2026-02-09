// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_GRIDFUNCTION_HPP
#define PALACE_FEM_GRIDFUNCTION_HPP

#include <mfem.hpp>
#include "linalg/vector.hpp"

namespace palace
{

class FiniteElementSpace;
class Mesh;

//
// A real- or complex-valued grid function represented as two real grid functions, one for
// each component. This unifies mfem::ParGridFunction and mfem::ParComplexGridFunction, and
// replaces the latter due to some issues observed with memory aliasing on GPUs.
//
class GridFunction
{
private:
  mfem::ParGridFunction gfr, gfi;

  // Back-references to palace wrappers. Set when constructed from
  // palace::FiniteElementSpace; null when constructed from raw MFEM types.
  FiniteElementSpace *fespace;
  Mesh *mesh;

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

  // Get access to the underlying palace mesh. Only available when the grid function was
  // constructed from a palace::FiniteElementSpace.
  const Mesh &GetMesh() const
  {
    MFEM_ASSERT(mesh, "GridFunction does not have a palace::Mesh reference!");
    return *mesh;
  }
  Mesh &GetMesh()
  {
    MFEM_ASSERT(mesh, "GridFunction does not have a palace::Mesh reference!");
    return *mesh;
  }
  bool HasMesh() const { return mesh != nullptr; }

  // Get access to the palace finite element space. Only available when constructed from
  // palace::FiniteElementSpace.
  const FiniteElementSpace &GetFESpace() const
  {
    MFEM_ASSERT(fespace,
                "GridFunction does not have a palace::FiniteElementSpace reference!");
    return *fespace;
  }
  FiniteElementSpace &GetFESpace()
  {
    MFEM_ASSERT(fespace,
                "GridFunction does not have a palace::FiniteElementSpace reference!");
    return *fespace;
  }
  bool HasFESpace() const { return fespace != nullptr; }

  // Get access to the underlying MFEM finite element space.
  mfem::FiniteElementSpace *FESpace() { return gfr.FESpace(); }
  const mfem::FiniteElementSpace *FESpace() const { return gfr.FESpace(); }
  mfem::ParFiniteElementSpace *ParFESpace() { return gfr.ParFESpace(); }
  const mfem::ParFiniteElementSpace *ParFESpace() const { return gfr.ParFESpace(); }

  // Dimension of the vector field represented by the grid function.
  int VectorDim() const { return gfr.VectorDim(); }

  // Convenience: max element order of the underlying FE space.
  auto GetMaxElementOrder() const { return gfr.ParFESpace()->GetMaxElementOrder(); }

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

  // Distribute true DOF values to the local grid function and optionally exchange
  // face neighbor data for parallel evaluation.
  void SetFromTrueDofs(const ComplexVector &tv, bool exchange_face_nbr = true)
  {
    Real().SetFromTrueDofs(tv.Real());
    if (HasImag())
    {
      Imag().SetFromTrueDofs(tv.Imag());
    }
    if (exchange_face_nbr)
    {
      ExchangeFaceNbrData();
    }
  }

  void SetFromTrueDofs(const Vector &tv, bool exchange_face_nbr = true)
  {
    Real().SetFromTrueDofs(tv);
    if (exchange_face_nbr)
    {
      ExchangeFaceNbrData();
    }
  }

  void ExchangeFaceNbrData()
  {
    Real().ExchangeFaceNbrData();
    if (HasImag())
    {
      Imag().ExchangeFaceNbrData();
    }
  }
};

}  // namespace palace

#endif  // PALACE_FEM_GRIDFUNCTION_HPP
