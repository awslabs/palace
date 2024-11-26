// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_L2_MASS_QF_H
#define PALACE_LIBCEED_L2_MASS_QF_H

// libCEED QFunctions for L2 + H(div) mass operators (Piola transformations u = 1 / det(J) ̂u
// and u = J / det(J) ̂u).
// Note: J / det(J) = adj(adj(J)^T / det(J))^T
// in[0] is geometry quadrature data, shape [ncomp=2+space_dim*dim, Q]
// in[1] is quadrature weights, shape [Q]
// in[2] is active vector, shape [qcomp=dim, ncomp=1, Q]
// in[3] is active vector divergence, shape [ncomp=1, Q]
// out[0] is active vector, shape [qcomp=dim, ncomp=1, Q]
// out[1] is active vector divergence, shape [ncomp=1, Q]

// Build functions assemble the quadrature point data.

#include "21/l2mass_21_qf.h"
#include "21/l2mass_build_21_qf.h"
#include "22/l2mass_22_qf.h"
#include "22/l2mass_build_22_qf.h"
#include "32/l2mass_32_qf.h"
#include "32/l2mass_build_32_qf.h"
#include "33/l2mass_33_qf.h"
#include "33/l2mass_build_33_qf.h"

#endif  // PALACE_LIBCEED_L2_MASS_QF_H
