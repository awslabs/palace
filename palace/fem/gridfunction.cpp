// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "gridfunction.hpp"
#include "fem/fespace.hpp"

namespace palace
{
namespace
{
// Same volume FE and Piola map as the source; all coefficients are element-interior.
// This is an evaluation/output collection, not an MFEM restart-file format.
class DiscontinuousCollection : public mfem::FiniteElementCollection
{
  const mfem::FiniteElementCollection &source;
  const int dim;

public:
  DiscontinuousCollection(const mfem::FiniteElementCollection &source, int dim)
    : mfem::FiniteElementCollection(source.GetOrder()), source(source), dim(dim)
  {
  }
  const mfem::FiniteElement *
  FiniteElementForGeometry(mfem::Geometry::Type geom) const override
  {
    return mfem::Geometry::Dimension[geom] == dim ? source.FiniteElementForGeometry(geom)
                                                  : nullptr;
  }
  int DofForGeometry(mfem::Geometry::Type geom) const override
  {
    const auto *fe = FiniteElementForGeometry(geom);
    return fe ? fe->GetDof() : 0;
  }
  const int *DofOrderForOrientation(mfem::Geometry::Type, int) const override
  {
    return nullptr;
  }
  const char *Name() const override { return "PalaceDG"; }
  int GetContType() const override { return DISCONTINUOUS; }
};
}  // namespace

struct GridFunction::DGData
{
  FiniteElementSpace &source;
  DiscontinuousCollection collection;
  mfem::ParFiniteElementSpace space;
  Vector true_real, true_imag, local_real, local_imag;
  DGData(FiniteElementSpace &source)
    : source(source), collection(source.GetFEColl(), source.Dimension()),
      space(&source.GetParMesh(), &collection)
  {
  }
};

GridFunction::GridFunction(mfem::ParFiniteElementSpace &space, bool complex) : gfr(&space)
{
  if (complex)
  {
    gfi.SetSpace(&space);
  }
}
GridFunction::GridFunction(FiniteElementSpace &space, bool complex)
  : GridFunction(space.Get(), complex)
{
}
GridFunction::GridFunction(FiniteElementSpace &space, Storage, bool complex)
  : dg(std::make_unique<DGData>(space)), gfr(&dg->space)
{
  if (complex)
  {
    gfi.SetSpace(&dg->space);
  }
  *this = 0.0;
}
GridFunction::~GridFunction() = default;

std::unique_ptr<GridFunction> GridFunction::MakeOutput(FiniteElementSpace &space,
                                                       bool complex)
{
  return space.IsBroken() ? std::make_unique<GridFunction>(space, Storage::DG, complex)
                          : std::make_unique<GridFunction>(space, complex);
}

void GridFunction::CopyToDG(const Vector &x, mfem::ParGridFunction &snapshot,
                            Vector &local) const
{
  auto &source = dg->source;
  MFEM_VERIFY(x.Size() == source.GetTrueVSize(), "Invalid source size for DG snapshot!");
  local.SetSize(source.GetVSize());
  source.GetProlongationMatrix()->Mult(x, local);
  mfem::Array<int> dofs, output;
  mfem::DofTransformation dt;
  Vector native;
  for (int e = 0; e < source.GetMesh().GetNE(); e++)
  {
    source.GetElementDofs(e, dofs, dt);
    local.GetSubVector(dofs, native);
    if (!dt.IsIdentity())
    {
      dt.InvTransformPrimal(native);
    }
    dg->space.GetElementDofs(e, output);
    snapshot.SetSubVector(output, native);
  }
}

void GridFunction::SetFromTrueDofs(const Vector &x)
{
  if (dg)
  {
    dg->true_real = x;
    CopyToDG(x, gfr, dg->local_real);
  }
  else
  {
    gfr.SetFromTrueDofs(x);
  }
}
void GridFunction::SetFromTrueDofs(const ComplexVector &x)
{
  MFEM_VERIFY(HasImag(), "Complex field needs a complex grid function!");
  SetFromTrueDofs(x.Real());
  if (dg)
  {
    dg->true_imag = x.Imag();
    CopyToDG(x.Imag(), gfi, dg->local_imag);
  }
  else
  {
    gfi.SetFromTrueDofs(x.Imag());
  }
}
void GridFunction::GetTrueDofs(Vector &y) const
{
  if (dg)
  {
    y = dg->true_real;
  }
  else
  {
    gfr.GetTrueDofs(y);
  }
}
void GridFunction::GetTrueDofs(ComplexVector &y) const
{
  MFEM_VERIFY(HasImag(), "Complex field needs a complex grid function!");
  y.SetSize(dg ? dg->source.GetTrueVSize() : gfr.ParFESpace()->GetTrueVSize());
  GetTrueDofs(y.Real());
  if (dg)
  {
    y.Imag() = dg->true_imag;
  }
  else
  {
    gfi.GetTrueDofs(y.Imag());
  }
}
const Vector &GridFunction::GetLVector(bool imaginary) const
{
  MFEM_VERIFY(!imaginary || HasImag(), "Imaginary part of a real field!");
  if (dg)
  {
    return imaginary ? dg->local_imag : dg->local_real;
  }
  return imaginary ? Imag() : Real();
}

GridFunction &GridFunction::operator=(std::complex<double> s)
{
  MFEM_VERIFY(HasImag() || s.imag() == 0.0, "Complex scalar assigned to a real field!");
  if (dg)
  {
    Vector real(dg->source.GetTrueVSize());
    real = s.real();
    SetFromTrueDofs(real);
    if (HasImag())
    {
      dg->true_imag.SetSize(real.Size());
      dg->true_imag = s.imag();
      CopyToDG(dg->true_imag, gfi, dg->local_imag);
    }
  }
  else
  {
    Real() = s.real();
    if (HasImag())
    {
      Imag() = s.imag();
    }
  }
  return *this;
}
GridFunction &GridFunction::operator*=(double s)
{
  if (dg)
  {
    dg->true_real *= s;
    dg->local_real *= s;
    if (HasImag())
    {
      dg->true_imag *= s;
      dg->local_imag *= s;
    }
  }
  // The face-neighbor data is a cached copy of the field values on neighboring ranks,
  // so it scales together with the field to stay consistent with it.
  Real() *= s;
  Real().FaceNbrData() *= s;
  if (HasImag())
  {
    Imag() *= s;
    Imag().FaceNbrData() *= s;
  }
  return *this;
}
void GridFunction::Update()
{
  if (dg)
  {
    dg->space.Update(false);
  }
  Real().Update();
  if (HasImag())
  {
    Imag().Update();
  }
  // Source spaces have been rebuilt; there is deliberately no DG-to-source transfer.
  if (dg)
  {
    *this = 0.0;
  }
}

}  // namespace palace
