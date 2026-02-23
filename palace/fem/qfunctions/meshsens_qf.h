// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_MESHSENS_QF_H
#define PALACE_LIBCEED_MESHSENS_QF_H

// libCEED QFunction for mesh sensitivity of the CurlCurlMass eigenvalue integrand.
//
// Computes df/dJ at each quadrature point, where f = f_K - lambda * f_M is the combined
// stiffness (curl-curl) and mass integrand. The derivative is computed via:
//   - Enzyme reverse-mode AD (if PALACE_WITH_ENZYME), or
//   - Forward-mode dual numbers (otherwise)
//
// in[0] is element attributes, shape [Q]         (CEED_EVAL_INTERP)
// in[1] is quadrature weights, shape [Q]         (CEED_EVAL_WEIGHT)
// in[2] is Jacobian J, shape [9*Q]               (CEED_EVAL_GRAD, active)
// in[3] is eigenvector u, shape [3*Q]            (CEED_EVAL_INTERP)
// in[4] is eigenvector curl, shape [3*Q]         (CEED_EVAL_CURL)
// out[0] is df/dJ sensitivity, shape [9*Q]       (CEED_EVAL_GRAD, active)

#include "33/meshsens_33_qf.h"

#endif  // PALACE_LIBCEED_MESHSENS_QF_H
