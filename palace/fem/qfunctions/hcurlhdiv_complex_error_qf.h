// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_HDIV_COMPLEX_ERROR_QF_H
#define PALACE_LIBCEED_HCURL_HDIV_COMPLEX_ERROR_QF_H

// libCEED QFunctions for the energy-weighted error of a complex constitutive flux,
// (εE - D)ᴴ W⁻¹ (εE - D) with the complex permittivity ε = ε' + i ε'' of a frequency-
// dependent material, between an H(curl) field E and an H(div) flux D (Piola
// transformations u = adj(J)^T / det(J) ̂u and u = J / det(J) ̂u).
// in[0] is geometry quadrature data, shape [ncomp=2+space_dim*dim, Q]
// in[1], in[2] are Re{E}, Im{E}, shape [qcomp=dim, ncomp=1, Q]
// in[3], in[4] are Re{D}, Im{D}, shape [qcomp=dim, ncomp=1, Q]
// out[0] is active vector, shape [ncomp=1, Q]
// The context holds the ε', ε'', and W⁻¹ coefficient tables contiguously.

// Only for the square Jacobian case where dim = space_dim.

#include "22/hcurlhdiv_complex_error_22_qf.h"
#include "33/hcurlhdiv_complex_error_33_qf.h"

#endif  // PALACE_LIBCEED_HCURL_HDIV_COMPLEX_ERROR_QF_H
