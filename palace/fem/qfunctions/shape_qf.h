// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_SHAPE_QF_H
#define PALACE_LIBCEED_SHAPE_QF_H

// libCEED QFunctions for the derivatives of forms a(u, w) with respect to the positions of
// the mesh nodes. For a node displacement V, the derivative is ∫ S : ∇V dx, so at each
// quadrature point the output on the gradient of the mesh nodal basis is the quadrature
// weight times |J| S adj(J)^T / |J|.
// in[0] is geometry quadrature data, shape [ncomp=2+space_dim*dim, Q]
// in[1], in[2] are u and w, shape [qcomp, ncomp=1, Q]: reference values, gradients or
// curls. The QFunctions of scalar curls take the quadrature weight before them, and those
// of linear forms w only.
// out[0] is the gradient output, shape [qcomp=dim, ncomp=space_dim, Q]

#include "21/shape_21_qf.h"
#include "22/shape_22_qf.h"
#include "32/shape_32_qf.h"
#include "33/shape_33_qf.h"

#endif  // PALACE_LIBCEED_SHAPE_QF_H
