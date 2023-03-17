// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_QUAD_RULES_HPP
#define PALACE_QUAD_RULES_HPP

#include <mfem.hpp>

namespace palace
{

namespace utils
{
// Helper functions for creating an integration rule to exactly integrate 2p + q
// polynomials. order_increment can be used to raise or lower the order, e.g. in
// the case of derivative fes.
inline const mfem::IntegrationRule *GetDefaultRule(const mfem::FiniteElement &trial_fe,
                                                   const mfem::FiniteElement &test_fe,
                                                   mfem::ElementTransformation &Tr,
                                                   int order_increment = 0)
{
  const int ir_order =
      trial_fe.GetOrder() + test_fe.GetOrder() + Tr.OrderW() + order_increment;
  MFEM_ASSERT(ir_order >= 0, "Negative integration order not allowed");
  return &mfem::IntRules.Get(trial_fe.GetGeomType(), ir_order);
}
inline const mfem::IntegrationRule *GetDefaultRule(const mfem::FiniteElement &fe,
                                                   mfem::ElementTransformation &Tr,
                                                   int order_increment = 0)
{
  return GetDefaultRule(fe, fe, Tr, order_increment);
}

}  // namespace utils

}  // namespace palace

#endif  // PALACE_QUAD_RULES_HPP
