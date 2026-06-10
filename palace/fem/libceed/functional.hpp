// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_FUNCTIONAL_HPP
#define PALACE_LIBCEED_FUNCTIONAL_HPP

#include <string>
#include <vector>
#include "fem/libceed/ceed.hpp"
#include "fem/libceed/integrator.hpp"

namespace palace::ceed
{

// Description of a passive field input for a surface functional operator. The field
// vector is evaluated at the (mapped face) quadrature points through the provided
// element restriction and basis.
struct CeedFunctionalFieldInput
{
  // Name of the QFunction input field (QFunction inputs are added in the order the
  // inputs are provided, after the geometry data inputs).
  std::string name;

  // Field vector (L-vector layout matching the restriction).
  CeedVector vec;

  // Element restriction and basis for evaluation at quadrature points.
  CeedElemRestriction restr;
  CeedBasis basis;

  // Evaluation modes (ceed::EvalMode) for the field input.
  unsigned int ops;
};

// Construct a libCEED operator which evaluates a functional (integral of a function of
// the provided fields) over every element described by the geometry data. The QFunction
// inputs are, in order: the face geometry data, optionally the volume geometry data
// (for evaluation of volume fields on boundary elements; pass nullptr to skip), and the
// field inputs. The QFunction output "v" (size num_out_comp) is summed over quadrature
// points via an all-ones basis, yielding num_out_comp values per element in the active
// output vector through out_restr. Apply with
// CeedOperatorApplyAdd(op, CEED_VECTOR_NONE, output, ...).
void AssembleCeedSurfaceFunctional(
    const CeedQFunctionInfo &info, void *ctx, std::size_t ctx_size, Ceed ceed,
    const std::vector<CeedFunctionalFieldInput> &inputs, CeedVector face_geom_data,
    CeedElemRestriction face_geom_data_restr, CeedVector vol_geom_data,
    CeedElemRestriction vol_geom_data_restr, CeedInt num_out_comp,
    CeedElemRestriction out_restr, CeedOperator *op);

// Construct a libCEED operator which evaluates a pointwise function of the provided
// fields at arbitrary points of each element (for example the nodal points of an
// interpolatory output space), writing num_out_comp values per point through out_restr
// (CEED_EVAL_NONE, so the number of "quadrature" points of the operator must match the
// element size of the output restriction). In contrast to the surface functional
// assembly, no quadrature weighting or sum over points is performed. The QFunction
// inputs are, in order: the geometry data (for Piola transformations and attribute
// based material lookup at the evaluation points), then the field inputs. Apply with
// CeedOperatorApplyAdd(op, CEED_VECTOR_NONE, output, ...) (contributions accumulate,
// e.g. for the real and imaginary part applications of quadratic quantities).
void AssembleCeedPointEvaluator(const CeedQFunctionInfo &info, void *ctx,
                                std::size_t ctx_size, Ceed ceed,
                                const std::vector<CeedFunctionalFieldInput> &inputs,
                                CeedVector geom_data, CeedElemRestriction geom_data_restr,
                                CeedInt num_out_comp, CeedElemRestriction out_restr,
                                CeedOperator *op);

}  // namespace palace::ceed

#endif  // PALACE_LIBCEED_FUNCTIONAL_HPP
