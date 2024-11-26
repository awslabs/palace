// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HDIV_MASS_QF_H
#define PALACE_LIBCEED_HDIV_MASS_QF_H

// libCEED QFunctions for H(div) + H(curl) mass operators in 3D (Piola transformations u =
// J / det(J) ̂u and u = adj(J)^T / det(J) ̂u).
// Note: J / det(J) = adj(adj(J)^T / det(J))^T
// in[0] is geometry quadrature data, shape [ncomp=2+space_dim*dim, Q]
// in[1] is active vector, shape [qcomp=dim, ncomp=1, Q]
// in[2] is active vector curl, shape [qcomp=dim, ncomp=1, Q]
// out[0] is active vector, shape [qcomp=dim, ncomp=1, Q]
// out[1] is active vector curl, shape [qcomp=dim, ncomp=1, Q]

// In 2D, this actually uses the L2 Piola transformation on the curl (u = 1 / det(J) ̂u) and
// the curl is has qcomp=1. There is no boundary integrator support in 2D.
// in[0] is geometry quadrature data, shape [ncomp=2+space_dim*dim, Q]
// in[1] is quadrature weights, shape [Q]
// in[2] is active vector, shape [qcomp=dim, ncomp=1, Q]
// in[3] is active vector curl, shape [ncomp=1, Q]
// out[0] is active vector, shape [qcomp=dim, ncomp=1, Q]
// out[1] is active vector curl, shape [ncomp=1, Q]

// Build functions assemble the quadrature point data.

#include "22/hdivmass_22_qf.h"
#include "22/hdivmass_build_22_qf.h"
#include "32/hdivmass_32_qf.h"
#include "32/hdivmass_build_32_qf.h"
#include "33/hdivmass_33_qf.h"
#include "33/hdivmass_build_33_qf.h"

#endif  // PALACE_LIBCEED_HDIV_MASS_QF_H
