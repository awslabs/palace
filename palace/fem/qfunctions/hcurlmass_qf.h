// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_MASS_QF_H
#define PALACE_LIBCEED_HCURL_MASS_QF_H

// libCEED QFunctions for H(curl) + H1 mass operators (Piola transformation u =
// adj(J)^T / det(J) ̂u and u = ̂u).
// in[0] is geometry quadrature data, shape [ncomp=2+space_dim*dim, Q]
// in[1] is active vector, shape [ncomp=1, Q]
// in[2] is active vector gradient, shape [qcomp=dim, ncomp=1, Q]
// out[0] is active vector, shape [ncomp=1, Q]
// out[1] is active vector gradient, shape [qcomp=dim, ncomp=1, Q]

// Build functions assemble the quadrature point data.

#include "21/hcurlmass_21_qf.h"
#include "21/hcurlmass_build_21_qf.h"
#include "22/hcurlmass_22_qf.h"
#include "22/hcurlmass_build_22_qf.h"
#include "32/hcurlmass_32_qf.h"
#include "32/hcurlmass_build_32_qf.h"
#include "33/hcurlmass_33_qf.h"
#include "33/hcurlmass_build_33_qf.h"

#endif  // PALACE_LIBCEED_HCURL_MASS_QF_H
