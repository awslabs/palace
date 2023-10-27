// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_FESPACE_HPP
#define PALACE_FEM_FESPACE_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"

namespace palace
{

//
// Wrapper for MFEM's ParFiniteElementSpace class, where the finite element space object
// is constructed with a unique ID associated with it. This is useful for defining equality
// operations between spaces (either different spaces on the same mesh, or the same space
// type on different meshes).
//
class FiniteElementSpace : public mfem::ParFiniteElementSpace
{
private:
  static std::size_t global_id;
  mutable std::size_t id;
  mutable long int prev_sequence;
  mutable bool init = false;

public:
  using mfem::ParFiniteElementSpace::ParFiniteElementSpace;
  FiniteElementSpace(const mfem::ParFiniteElementSpace &fespace)
    : mfem::ParFiniteElementSpace(fespace)
  {
  }

  // Get the ID associated with the instance of this class. If the underlying sequence has
  // changed (due to a mesh update, for example), regenerate the ID.
  std::size_t GetId() const;
};

//
// An AuxiliaryFiniteElement space is just a FiniteElementSpace which allows for lazy
// construction of the interpolation operator (discrete gradient or curl) from the primal
// space to this one.
//
class AuxiliaryFiniteElementSpace : public FiniteElementSpace
{
private:
  const FiniteElementSpace &primal_fespace;
  mutable std::unique_ptr<Operator> G;

  const Operator &BuildDiscreteInterpolator() const;

public:
  template <typename... T>
  AuxiliaryFiniteElementSpace(const FiniteElementSpace &primal_fespace, T &&...args)
    : FiniteElementSpace(std::forward<T>(args)...), primal_fespace(primal_fespace)
  {
  }

  // Return the discrete gradient or discrete curl matrix interpolating from the auxiliary
  // to the primal space, constructing it on the fly as necessary.
  const Operator &GetDiscreteInterpolator() const
  {
    return G ? *G : BuildDiscreteInterpolator();
  }
};

//
// A collection of FiniteElementSpace objects constructed on the same mesh with the ability
// to construct the prolongation operators between them as needed.
//
template <typename FESpace>
class Hierarchy
{
protected:
  std::vector<std::unique_ptr<FESpace>> fespaces;
  mutable std::vector<std::unique_ptr<Operator>> P;

  static_assert(std::is_convertible_v<FESpace *, FiniteElementSpace *>,
                "A hierarchy can only be constructed of FiniteElementSpaces");

  const Operator &BuildProlongationAtLevel(std::size_t l) const;

public:
  Hierarchy<FESpace>() = default;
  Hierarchy<FESpace>(std::unique_ptr<FESpace> &&fespace) { AddLevel(std::move(fespace)); }

  auto GetNumLevels() const { return fespaces.size(); }
  auto size() const { return GetNumLevels(); }
  bool empty() const { return GetNumLevels() == 0; }

  virtual void AddLevel(std::unique_ptr<FESpace> &&fespace)
  {
    fespaces.push_back(std::move(fespace));
    P.push_back(nullptr);
  }

  FESpace &GetFESpaceAtLevel(std::size_t l)
  {
    MFEM_ASSERT(l >= 0 && l < GetNumLevels(),
                "Out of bounds request for finite element space at level " << l << "!");
    return *fespaces[l];
  }
  const FESpace &GetFESpaceAtLevel(std::size_t l) const
  {
    MFEM_ASSERT(l >= 0 && l < GetNumLevels(),
                "Out of bounds request for finite element space at level " << l << "!");
    return *fespaces[l];
  }

  FESpace &GetFinestFESpace()
  {
    MFEM_ASSERT(!empty(), "Out of bounds request for finite element space at level 0!");
    return *fespaces.back();
  }
  const FESpace &GetFinestFESpace() const
  {
    MFEM_ASSERT(!empty(), "Out of bounds request for finite element space at level 0!");
    return *fespaces.back();
  }

  const Operator &GetProlongationAtLevel(std::size_t l) const
  {
    MFEM_ASSERT(l >= 0 && l < GetNumLevels() - 1,
                "Out of bounds request for finite element space prolongation at level "
                    << l << "!");
    return P[l] ? *P[l] : BuildProlongationAtLevel(l);
  }

  std::vector<const Operator *> GetProlongationOperators() const
  {
    std::vector<const Operator *> P_(GetNumLevels() - 1);
    for (std::size_t l = 0; l < P_.size(); l++)
    {
      P_[l] = &GetProlongationAtLevel(l);
    }
    return P_;
  }
};

class FiniteElementSpaceHierarchy : public Hierarchy<FiniteElementSpace>
{
public:
  using Hierarchy<FiniteElementSpace>::Hierarchy;
};

//
// A special type of FiniteElementSpaceHierarchy where all members are auxiliary finite
// element spaces.
//
class AuxiliaryFiniteElementSpaceHierarchy : public Hierarchy<AuxiliaryFiniteElementSpace>
{

public:
  using Hierarchy<AuxiliaryFiniteElementSpace>::Hierarchy;

  const Operator &GetDiscreteInterpolatorAtLevel(std::size_t l) const
  {
    return GetFESpaceAtLevel(l).GetDiscreteInterpolator();
  }

  const Operator &GetDiscreteInterpolatorAtFinestLevel() const
  {
    return GetFinestFESpace().GetDiscreteInterpolator();
  }

  std::vector<const Operator *> GetDiscreteInterpolators() const
  {
    std::vector<const Operator *> G_(GetNumLevels());
    for (std::size_t l = 0; l < G_.size(); l++)
    {
      G_[l] = &GetDiscreteInterpolatorAtLevel(l);
    }
    return G_;
  }
};

}  // namespace palace

#endif  // PALACE_FEM_FESPACE_HPP
