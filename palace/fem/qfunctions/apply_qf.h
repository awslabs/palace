// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_APPLY_QF_H
#define PALACE_LIBCEED_APPLY_QF_H

// libCEED QFunctions for application of a generic operator with assembled quadrature data.
// in[0] is quadrature data, shape [ncomp=vdim*vdim, Q]
// in[1] is active vector, shape [ncomp=vdim, Q]
// out[0] is active vector, shape [ncomp=vdim, Q]

// For pairwise apply functions, the inputs and outputs come in pairs and the quadrature
// data is arranged to be applied with the first vdim*vdim components for the first
// input/output and the remainder for the second.

#include "apply/apply_12_qf.h"
#include "apply/apply_13_qf.h"
#include "apply/apply_1_qf.h"
#include "apply/apply_21_qf.h"
#include "apply/apply_22_qf.h"
#include "apply/apply_2_qf.h"
#include "apply/apply_31_qf.h"
#include "apply/apply_33_qf.h"
#include "apply/apply_3_qf.h"

#endif  // PALACE_LIBCEED_APPLY_QF_H
