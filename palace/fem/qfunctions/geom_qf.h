// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_GEOM_QF_H
#define PALACE_LIBCEED_GEOM_QF_H

// libCEED QFunction for building geometry factors for integration and transformations.
// At every quadrature point, compute qw * det(J) and adj(J)^T / |J| and store the result.
// in[0] is quadrature weights, shape [Q]
// in[1] is Jacobians, shape [qcomp=dim, ncomp=space_dim, Q]
// out[0] is quadrature data, stored as {attribute, Jacobian determinant, (transpose)
//        adjugate Jacobian} quadrature data, shape [ncomp=2+space_dim*dim, Q]

#include "21/geom_21_qf.h"
#include "22/geom_22_qf.h"
#include "32/geom_32_qf.h"
#include "33/geom_33_qf.h"

#endif  // PALACE_LIBCEED_GEOM_QF_H
