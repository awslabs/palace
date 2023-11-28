// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_UTILS_QF_H
#define PALACE_LIBCEED_UTILS_QF_H

CEED_QFUNCTION_HELPER void MatUnpack22(const CeedScalar *A, const CeedInt A_stride,
                                       CeedScalar A_loc[4])
{
  A_loc[0] = A[A_stride * 0];
  A_loc[1] = A[A_stride * 1];
  A_loc[2] = A[A_stride * 2];
  A_loc[3] = A[A_stride * 3];
}

CEED_QFUNCTION_HELPER void MatUnpack33(const CeedScalar *A, const CeedInt A_stride,
                                       CeedScalar A_loc[9])
{
  A_loc[0] = A[A_stride * 0];
  A_loc[1] = A[A_stride * 1];
  A_loc[2] = A[A_stride * 2];
  A_loc[3] = A[A_stride * 3];
  A_loc[4] = A[A_stride * 4];
  A_loc[5] = A[A_stride * 5];
  A_loc[6] = A[A_stride * 6];
  A_loc[7] = A[A_stride * 7];
  A_loc[8] = A[A_stride * 8];
  A_loc[9] = A[A_stride * 9];
}

CEED_QFUNCTION_HELPER void MatUnpack21(const CeedScalar *A, const CeedInt A_stride,
                                       CeedScalar A_loc[2])
{
  A_loc[0] = A[A_stride * 0];
  A_loc[1] = A[A_stride * 1];
}

CEED_QFUNCTION_HELPER void MatUnpack32(const CeedScalar *A, const CeedInt A_stride,
                                       CeedScalar A_loc[6])
{
  A_loc[0] = A[A_stride * 0];
  A_loc[1] = A[A_stride * 1];
  A_loc[2] = A[A_stride * 2];
  A_loc[3] = A[A_stride * 3];
  A_loc[4] = A[A_stride * 4];
  A_loc[5] = A[A_stride * 5];
  A_loc[6] = A[A_stride * 6];
}

CEED_QFUNCTION_HELPER void MultAtBA22(const CeedScalar *A, const CeedInt A_stride,
                                      const CeedScalar B[3], CeedScalar C[3])
{
  // A: 0 2   B: 0 1   C: 0 1
  //    1 3      1 2      1 2
  CeedScalar A_loc[4];
  MatUnpack22(A, A_stride, A_loc);

  // First compute entries of R = B A.
  const CeedScalar R11 = B[0] * A_loc[0] + B[1] * A_loc[1];
  const CeedScalar R21 = B[1] * A_loc[0] + B[2] * A_loc[1];
  const CeedScalar R12 = B[0] * A_loc[2] + B[1] * A_loc[3];
  const CeedScalar R22 = B[1] * A_loc[2] + B[2] * A_loc[3];

  C[0] = A_loc[0] * R11 + A_loc[1] * R21;
  C[1] = A_loc[0] * R12 + A_loc[1] * R22;
  C[2] = A_loc[2] * R12 + A_loc[3] * R22;
}

CEED_QFUNCTION_HELPER void MultAtBA33(const CeedScalar *A, const CeedInt A_stride,
                                      const CeedScalar B[6], CeedScalar C[6])
{
  // A: 0 3 6   B: 0 1 2   C: 0 1 2
  //    1 4 7      1 3 4      1 3 4
  //    2 5 8      2 4 5      2 4 5
  CeedScalar A_loc[9];
  MatUnpack33(A, A_stride, A_loc);

  // First compute entries of R = B A.
  const CeedScalar R11 = B[0] * A_loc[0] + B[1] * A_loc[1] + B[2] * A_loc[2];
  const CeedScalar R21 = B[1] * A_loc[0] + B[3] * A_loc[1] + B[4] * A_loc[2];
  const CeedScalar R31 = B[2] * A_loc[0] + B[4] * A_loc[1] + B[5] * A_loc[2];
  const CeedScalar R12 = B[0] * A_loc[3] + B[1] * A_loc[4] + B[2] * A_loc[5];
  const CeedScalar R22 = B[1] * A_loc[3] + B[3] * A_loc[4] + B[4] * A_loc[5];
  const CeedScalar R32 = B[2] * A_loc[3] + B[4] * A_loc[4] + B[5] * A_loc[5];
  const CeedScalar R13 = B[0] * A_loc[6] + B[1] * A_loc[7] + B[2] * A_loc[8];
  const CeedScalar R23 = B[1] * A_loc[6] + B[3] * A_loc[7] + B[4] * A_loc[8];
  const CeedScalar R33 = B[2] * A_loc[6] + B[4] * A_loc[7] + B[5] * A_loc[8];

  C[0] = A_loc[0] * R11 + A_loc[1] * R21 + A_loc[2] * R31;
  C[1] = A_loc[0] * R12 + A_loc[1] * R22 + A_loc[2] * R32;
  C[2] = A_loc[0] * R13 + A_loc[1] * R23 + A_loc[2] * R33;
  C[3] = A_loc[3] * R12 + A_loc[4] * R22 + A_loc[5] * R32;
  C[4] = A_loc[3] * R13 + A_loc[4] * R23 + A_loc[5] * R33;
  C[5] = A_loc[6] * R13 + A_loc[7] * R23 + A_loc[8] * R33;
}

CEED_QFUNCTION_HELPER void MultAtBA21(const CeedScalar *A, const CeedInt A_stride,
                                      const CeedScalar B[3], CeedScalar C[1])
{
  // A: 0   B: 0 1   C: 0
  //    1      1 2
  CeedScalar A_loc[2];
  MatUnpack21(A, A_stride, A_loc);

  // First compute entries of R = B A.
  const CeedScalar R11 = B[0] * A_loc[0] + B[1] * A_loc[1];
  const CeedScalar R21 = B[1] * A_loc[0] + B[2] * A_loc[1];

  C[0] = A_loc[0] * R11 + A_loc[1] * R21;
}

CEED_QFUNCTION_HELPER void MultAtBA32(const CeedScalar *A, const CeedInt A_stride,
                                      const CeedScalar B[6], CeedScalar C[3])
{
  // A: 0 3   B: 0 1 2   C: 0 1
  //    1 4      1 3 4      1 2
  //    2 5      2 4 5
  CeedScalar A_loc[6];
  MatUnpack32(A, A_stride, A_loc);

  // First compute entries of R = B A.
  const CeedScalar R11 = B[0] * A_loc[0] + B[1] * A_loc[1] + B[2] * A_loc[2];
  const CeedScalar R21 = B[1] * A_loc[0] + B[3] * A_loc[1] + B[4] * A_loc[2];
  const CeedScalar R31 = B[2] * A_loc[0] + B[4] * A_loc[1] + B[5] * A_loc[2];
  const CeedScalar R12 = B[0] * A_loc[3] + B[1] * A_loc[4] + B[2] * A_loc[5];
  const CeedScalar R22 = B[1] * A_loc[3] + B[3] * A_loc[4] + B[4] * A_loc[5];
  const CeedScalar R32 = B[2] * A_loc[3] + B[4] * A_loc[4] + B[5] * A_loc[5];

  C[0] = A_loc[0] * R11 + A_loc[1] * R21 + A_loc[2] * R31;
  C[1] = A_loc[0] * R12 + A_loc[1] * R22 + A_loc[2] * R32;
  C[2] = A_loc[3] * R12 + A_loc[4] * R22 + A_loc[5] * R32;
}

CEED_QFUNCTION_HELPER void MultAtBC22(const CeedScalar *A, const CeedInt A_stride,
                                      const CeedScalar B[3], const CeedScalar *C,
                                      const CeedInt C_stride, CeedScalar D[4])
{
  // A, C: 0 2   B: 0 1   D: 0 2
  //       1 3      1 2      1 3
  CeedScalar A_loc[4];
  MatUnpack22(C, C_stride, A_loc);

  // First compute entries of R = B C.
  const CeedScalar R11 = B[0] * A_loc[0] + B[1] * A_loc[1];
  const CeedScalar R21 = B[1] * A_loc[0] + B[2] * A_loc[1];
  const CeedScalar R12 = B[0] * A_loc[2] + B[1] * A_loc[3];
  const CeedScalar R22 = B[1] * A_loc[2] + B[2] * A_loc[3];

  MatUnpack22(A, A_stride, A_loc);

  D[0] = A_loc[0] * R11 + A_loc[1] * R21;
  D[1] = A_loc[2] * R11 + A_loc[3] * R21;
  D[2] = A_loc[0] * R12 + A_loc[1] * R22;
  D[3] = A_loc[2] * R12 + A_loc[3] * R22;
}

CEED_QFUNCTION_HELPER void MultAtBC33(const CeedScalar *A, const CeedInt A_stride,
                                      const CeedScalar B[6], const CeedScalar *C,
                                      const CeedInt C_stride, CeedScalar D[9])
{
  // A, C: 0 3 6   B: 0 1 2   D: 0 3 6
  //       1 4 7      1 3 4      1 4 7
  //       2 5 8      2 4 5      2 5 8
  CeedScalar A_loc[9];
  MatUnpack33(C, C_stride, A_loc);

  // First compute entries of R = B C.
  const CeedScalar R11 = B[0] * A_loc[0] + B[1] * A_loc[1] + B[2] * A_loc[2];
  const CeedScalar R21 = B[1] * A_loc[0] + B[3] * A_loc[1] + B[4] * A_loc[2];
  const CeedScalar R31 = B[2] * A_loc[0] + B[4] * A_loc[1] + B[5] * A_loc[2];
  const CeedScalar R12 = B[0] * A_loc[3] + B[1] * A_loc[4] + B[2] * A_loc[5];
  const CeedScalar R22 = B[1] * A_loc[3] + B[3] * A_loc[4] + B[4] * A_loc[5];
  const CeedScalar R32 = B[2] * A_loc[3] + B[4] * A_loc[4] + B[5] * A_loc[5];
  const CeedScalar R13 = B[0] * A_loc[6] + B[1] * A_loc[7] + B[2] * A_loc[8];
  const CeedScalar R23 = B[1] * A_loc[6] + B[3] * A_loc[7] + B[4] * A_loc[8];
  const CeedScalar R33 = B[2] * A_loc[6] + B[4] * A_loc[7] + B[5] * A_loc[8];

  MatUnpack33(A, A_stride, A_loc);

  D[0] = A_loc[0] * R11 + A_loc[1] * R21 + A_loc[2] * R31;
  D[1] = A_loc[3] * R11 + A_loc[4] * R21 + A_loc[5] * R31;
  D[2] = A_loc[6] * R11 + A_loc[7] * R21 + A_loc[8] * R31;
  D[3] = A_loc[0] * R12 + A_loc[1] * R22 + A_loc[2] * R32;
  D[4] = A_loc[3] * R12 + A_loc[4] * R22 + A_loc[5] * R32;
  D[5] = A_loc[6] * R12 + A_loc[7] * R22 + A_loc[8] * R32;
  D[6] = A_loc[0] * R13 + A_loc[1] * R23 + A_loc[2] * R33;
  D[7] = A_loc[3] * R13 + A_loc[4] * R23 + A_loc[5] * R33;
  D[8] = A_loc[6] * R13 + A_loc[7] * R23 + A_loc[8] * R33;
}

CEED_QFUNCTION_HELPER void MultAtBC21(const CeedScalar *A, const CeedInt A_stride,
                                      const CeedScalar B[3], const CeedScalar *C,
                                      const CeedInt C_stride, CeedScalar D[1])
{
  // A, C: 0   B: 0 1   D: 0
  //       1      1 2
  CeedScalar A_loc[2];
  MatUnpack21(C, C_stride, A_loc);

  // First compute entries of R = B C.
  const CeedScalar R11 = B[0] * A_loc[0] + B[1] * A_loc[1];
  const CeedScalar R21 = B[1] * A_loc[0] + B[2] * A_loc[1];

  MatUnpack21(A, A_stride, A_loc);

  D[0] = A_loc[0] * R11 + A_loc[1] * R21;
}

CEED_QFUNCTION_HELPER void MultAtBC32(const CeedScalar *A, const CeedInt A_stride,
                                      const CeedScalar B[6], const CeedScalar *C,
                                      const CeedInt C_stride, CeedScalar D[4])
{
  // A, C: 0 3   B: 0 1 2   D: 0 2
  //       1 4      1 3 4      1 3
  //       2 5      2 4 5
  CeedScalar A_loc[6];
  MatUnpack32(C, C_stride, A_loc);

  // First compute entries of R = B C.
  const CeedScalar R11 = B[0] * A_loc[0] + B[1] * A_loc[1] + B[2] * A_loc[2];
  const CeedScalar R21 = B[1] * A_loc[0] + B[3] * A_loc[1] + B[4] * A_loc[2];
  const CeedScalar R31 = B[2] * A_loc[0] + B[4] * A_loc[1] + B[5] * A_loc[2];
  const CeedScalar R12 = B[0] * A_loc[3] + B[1] * A_loc[4] + B[2] * A_loc[5];
  const CeedScalar R22 = B[1] * A_loc[3] + B[3] * A_loc[4] + B[4] * A_loc[5];
  const CeedScalar R32 = B[2] * A_loc[3] + B[4] * A_loc[4] + B[5] * A_loc[5];

  MatUnpack32(A, A_stride, A_loc);

  D[0] = A_loc[0] * R11 + A_loc[1] * R21 + A_loc[2] * R31;
  D[1] = A_loc[3] * R11 + A_loc[4] * R21 + A_loc[5] * R31;
  D[2] = A_loc[0] * R12 + A_loc[1] * R22 + A_loc[2] * R32;
  D[3] = A_loc[3] * R12 + A_loc[4] * R22 + A_loc[5] * R32;
}

CEED_QFUNCTION_HELPER void MultBA22(const CeedScalar *A, const CeedInt A_stride,
                                    const CeedScalar B[3], CeedScalar C[4])
{
  // A: 0 2   B: 0 1   C: 0 2
  //    1 3      1 2      1 3
  CeedScalar A_loc[4];
  MatUnpack22(A, A_stride, A_loc);

  C[0] = B[0] * A_loc[0] + B[1] * A_loc[1];
  C[1] = B[1] * A_loc[0] + B[2] * A_loc[1];
  C[2] = B[0] * A_loc[2] + B[1] * A_loc[3];
  C[3] = B[1] * A_loc[2] + B[2] * A_loc[3];
}

CEED_QFUNCTION_HELPER void MultBA33(const CeedScalar *A, const CeedInt A_stride,
                                    const CeedScalar B[6], CeedScalar C[9])
{
  // A: 0 3 6   B: 0 1 2   C: 0 3 6
  //    1 4 7      1 3 4      1 4 7
  //    2 5 8      2 4 5      2 5 8
  CeedScalar A_loc[9];
  MatUnpack33(A, A_stride, A_loc);

  C[0] = B[0] * A_loc[0] + B[1] * A_loc[1] + B[2] * A_loc[2];
  C[1] = B[1] * A_loc[0] + B[3] * A_loc[1] + B[4] * A_loc[2];
  C[2] = B[2] * A_loc[0] + B[4] * A_loc[1] + B[5] * A_loc[2];
  C[3] = B[0] * A_loc[3] + B[1] * A_loc[4] + B[2] * A_loc[5];
  C[4] = B[1] * A_loc[3] + B[3] * A_loc[4] + B[4] * A_loc[5];
  C[5] = B[2] * A_loc[3] + B[4] * A_loc[4] + B[5] * A_loc[5];
  C[6] = B[0] * A_loc[6] + B[1] * A_loc[7] + B[2] * A_loc[8];
  C[7] = B[1] * A_loc[6] + B[3] * A_loc[7] + B[4] * A_loc[8];
  C[8] = B[2] * A_loc[6] + B[4] * A_loc[7] + B[5] * A_loc[8];
}

CEED_QFUNCTION_HELPER void MultBA21(const CeedScalar *A, const CeedInt A_stride,
                                    const CeedScalar B[3], CeedScalar C[2])
{
  // A: 0   B: 0 1   C: 0
  //    1      1 2      1
  CeedScalar A_loc[2];
  MatUnpack21(A, A_stride, A_loc);

  C[0] = B[0] * A_loc[0] + B[1] * A_loc[1];
  C[1] = B[1] * A_loc[0] + B[2] * A_loc[1];
}

CEED_QFUNCTION_HELPER void MultBA32(const CeedScalar *A, const CeedInt A_stride,
                                    const CeedScalar B[6], CeedScalar C[6])
{
  // A: 0 3   B: 0 1 2   C: 0 3
  //    1 4      1 3 4      1 4
  //    2 5      2 4 5      2 5
  CeedScalar A_loc[6];
  MatUnpack32(A, A_stride, A_loc);

  C[0] = B[0] * A_loc[0] + B[1] * A_loc[1] + B[2] * A_loc[2];
  C[1] = B[1] * A_loc[0] + B[3] * A_loc[1] + B[4] * A_loc[2];
  C[2] = B[2] * A_loc[0] + B[4] * A_loc[1] + B[5] * A_loc[2];
  C[3] = B[0] * A_loc[3] + B[1] * A_loc[4] + B[2] * A_loc[5];
  C[4] = B[1] * A_loc[3] + B[3] * A_loc[4] + B[4] * A_loc[5];
  C[5] = B[2] * A_loc[3] + B[4] * A_loc[4] + B[5] * A_loc[5];
}

// XX TODO WIP REMOVE

// CEED_QFUNCTION_HELPER void MultAdjJCAdjJt22(const CeedScalar *J, const CeedInt J_stride,
//                                             const CeedScalar *c, const CeedInt c_stride,
//                                             const CeedInt c_comp, const CeedScalar qw,
//                                             const CeedInt qd_stride, CeedScalar *qd)
// {
//   // Compute qw / det(J) adj(J) C adj(J)^T and store the symmetric part of the result.
//   // J: 0 2   adj(J):  J22 -J12   qd: 0 1
//   //    1 3           -J21  J11       1 2
//   const CeedScalar J11 = J[J_stride * 0];
//   const CeedScalar J21 = J[J_stride * 1];
//   const CeedScalar J12 = J[J_stride * 2];
//   const CeedScalar J22 = J[J_stride * 3];
//   const CeedScalar w = qw / (J11 * J22 - J21 * J12);
//   if (c_comp == 3)  // Matrix coefficient (symmetric)
//   {
//     // First compute entries of R = C adj(J)^T.
//     // c: 0 1
//     //    1 2
//     const CeedScalar R11 = c[c_stride * 0] * J22 - c[c_stride * 1] * J12;
//     const CeedScalar R21 = c[c_stride * 1] * J22 - c[c_stride * 2] * J12;
//     const CeedScalar R12 = -c[c_stride * 0] * J21 + c[c_stride * 1] * J11;
//     const CeedScalar R22 = -c[c_stride * 1] * J21 + c[c_stride * 2] * J11;
//     qd[qd_stride * 0] = w * (J22 * R11 - J12 * R21);
//     qd[qd_stride * 1] = w * (J11 * R21 - J21 * R11);
//     qd[qd_stride * 2] = w * (J11 * R22 - J21 * R12);
//   }
//   else if (c_comp == 2)  // Vector coefficient
//   {
//     // c: 0
//     //      1
//     qd[qd_stride * 0] = w * (c[c_stride * 1] * J12 * J12 + c[c_stride * 0] * J22 * J22);
//     qd[qd_stride * 1] = -w * (c[c_stride * 1] * J11 * J12 + c[c_stride * 0] * J21 * J22);
//     qd[qd_stride * 2] = w * (c[c_stride * 1] * J11 * J11 + c[c_stride * 0] * J21 * J21);
//   }
//   else  // Scalar coefficient
//   {
//     qd[qd_stride * 0] = w * c[c_stride * 0] * (J12 * J12 + J22 * J22);
//     qd[qd_stride * 1] = -w * c[c_stride * 0] * (J11 * J12 + J21 * J22);
//     qd[qd_stride * 2] = w * c[c_stride * 0] * (J11 * J11 + J21 * J21);
//   }
// }

// CEED_QFUNCTION_HELPER void MultAdjJCAdjJt21(const CeedScalar *J, const CeedInt J_stride,
//                                             const CeedScalar *c, const CeedInt c_stride,
//                                             const CeedInt c_comp, const CeedScalar qw,
//                                             const CeedInt qd_stride, CeedScalar *qd)
// {
//   // Compute qw / det(J) adj(J) C adj(J)^T and store the symmetric part of the result.
//   // J: 0   adj(J): 1/sqrt(J^T J) J^T   qd: 0
//   //    1
//   const CeedScalar J11 = J[J_stride * 0];
//   const CeedScalar J21 = J[J_stride * 1];
//   const CeedScalar d = J11 * J11 + J21 * J21;
//   const CeedScalar w = qw / sqrt(d);
//   if (c_comp == 3)  // Matrix coefficient (symmetric)
//   {
//     // First compute entries of R = C adj(J)^T.
//     // c: 0 1
//     //    1 2
//     const CeedScalar R11 = c[c_stride * 0] * J11 + c[c_stride * 1] * J21;
//     const CeedScalar R21 = c[c_stride * 1] * J11 + c[c_stride * 2] * J21;
//     qd[qd_stride * 0] = w * (J11 * R11 + J21 * R21) / d;
//   }
//   else if (c_comp == 2)  // Vector coefficient
//   {
//     // c: 0
//     //      1
//     qd[qd_stride * 0] = w * (c[c_stride * 0] * J11 * J11 + c[c_stride * 1] * J21 * J21) /
//     d;
//   }
//   else  // Scalar coefficient
//   {
//     qd[qd_stride * 0] = w * c[c_stride * 0];
//   }
// }

// CEED_QFUNCTION_HELPER void MultAdjJCAdjJt33(const CeedScalar *J, const CeedInt J_stride,
//                                             const CeedScalar *c, const CeedInt c_stride,
//                                             const CeedInt c_comp, const CeedScalar qw,
//                                             const CeedInt qd_stride, CeedScalar *qd)
// {
//   // Compute qw / det(J) adj(J) C adj(J)^T and store the symmetric part of the result.
//   // J: 0 3 6   qd: 0 1 2
//   //    1 4 7       1 3 4
//   //    2 5 8       2 4 5
//   const CeedScalar J11 = J[J_stride * 0];
//   const CeedScalar J21 = J[J_stride * 1];
//   const CeedScalar J31 = J[J_stride * 2];
//   const CeedScalar J12 = J[J_stride * 3];
//   const CeedScalar J22 = J[J_stride * 4];
//   const CeedScalar J32 = J[J_stride * 5];
//   const CeedScalar J13 = J[J_stride * 6];
//   const CeedScalar J23 = J[J_stride * 7];
//   const CeedScalar J33 = J[J_stride * 8];
//   const CeedScalar A11 = J22 * J33 - J23 * J32;
//   const CeedScalar A21 = J23 * J31 - J21 * J33;
//   const CeedScalar A31 = J21 * J32 - J22 * J31;
//   const CeedScalar A12 = J13 * J32 - J12 * J33;
//   const CeedScalar A22 = J11 * J33 - J13 * J31;
//   const CeedScalar A32 = J12 * J31 - J11 * J32;
//   const CeedScalar A13 = J12 * J23 - J13 * J22;
//   const CeedScalar A23 = J13 * J21 - J11 * J23;
//   const CeedScalar A33 = J11 * J22 - J12 * J21;
//   const CeedScalar w = qw / (J11 * A11 + J21 * A12 + J31 * A13);
//   if (c_comp == 6)  // Matrix coefficient (symmetric)
//   {
//     // First compute entries of R = C adj(J)^T.
//     // c: 0 1 2
//     //    1 3 4
//     //    2 4 5
//     const CeedScalar R11 =
//         c[c_stride * 0] * A11 + c[c_stride * 1] * A12 + c[c_stride * 2] * A13;
//     const CeedScalar R21 =
//         c[c_stride * 1] * A11 + c[c_stride * 3] * A12 + c[c_stride * 4] * A13;
//     const CeedScalar R31 =
//         c[c_stride * 2] * A11 + c[c_stride * 4] * A12 + c[c_stride * 5] * A13;
//     const CeedScalar R12 =
//         c[c_stride * 0] * A21 + c[c_stride * 1] * A22 + c[c_stride * 2] * A23;
//     const CeedScalar R22 =
//         c[c_stride * 1] * A21 + c[c_stride * 3] * A22 + c[c_stride * 4] * A23;
//     const CeedScalar R32 =
//         c[c_stride * 2] * A21 + c[c_stride * 4] * A22 + c[c_stride * 5] * A23;
//     const CeedScalar R13 =
//         c[c_stride * 0] * A31 + c[c_stride * 1] * A32 + c[c_stride * 2] * A33;
//     const CeedScalar R23 =
//         c[c_stride * 1] * A31 + c[c_stride * 3] * A32 + c[c_stride * 4] * A33;
//     const CeedScalar R33 =
//         c[c_stride * 2] * A31 + c[c_stride * 4] * A32 + c[c_stride * 5] * A33;
//     qd[qd_stride * 0] = w * (A11 * R11 + A12 * R21 + A13 * R31);
//     qd[qd_stride * 1] = w * (A11 * R12 + A12 * R22 + A13 * R32);
//     qd[qd_stride * 2] = w * (A11 * R13 + A12 * R23 + A13 * R33);
//     qd[qd_stride * 3] = w * (A21 * R12 + A22 * R22 + A23 * R32);
//     qd[qd_stride * 4] = w * (A21 * R13 + A22 * R23 + A23 * R33);
//     qd[qd_stride * 5] = w * (A31 * R13 + A32 * R23 + A33 * R33);
//   }
//   else if (c_comp == 3)  // Vector coefficient
//   {
//     // c: 0
//     //      1
//     //        2
//     qd[qd_stride * 0] = w * (c[c_stride * 0] * A11 * A11 + c[c_stride * 1] * A12 * A12 +
//                              c[c_stride * 2] * A13 * A13);
//     qd[qd_stride * 1] = w * (c[c_stride * 0] * A11 * A21 + c[c_stride * 1] * A12 * A22 +
//                              c[c_stride * 2] * A13 * A23);
//     qd[qd_stride * 2] = w * (c[c_stride * 0] * A11 * A31 + c[c_stride * 1] * A12 * A32 +
//                              c[c_stride * 2] * A13 * A33);
//     qd[qd_stride * 3] = w * (c[c_stride * 0] * A21 * A21 + c[c_stride * 1] * A22 * A22 +
//                              c[c_stride * 2] * A23 * A23);
//     qd[qd_stride * 4] = w * (c[c_stride * 0] * A21 * A31 + c[c_stride * 1] * A22 * A32 +
//                              c[c_stride * 2] * A23 * A33);
//     qd[qd_stride * 5] = w * (c[c_stride * 0] * A31 * A31 + c[c_stride * 1] * A32 * A32 +
//                              c[c_stride * 2] * A33 * A33);
//   }
//   else  // Scalar coefficient
//   {
//     qd[qd_stride * 0] = w * c[c_stride * 0] * (A11 * A11 + A12 * A12 + A13 * A13);
//     qd[qd_stride * 1] = w * c[c_stride * 0] * (A11 * A21 + A12 * A22 + A13 * A23);
//     qd[qd_stride * 2] = w * c[c_stride * 0] * (A11 * A31 + A12 * A32 + A13 * A33);
//     qd[qd_stride * 3] = w * c[c_stride * 0] * (A21 * A21 + A22 * A22 + A23 * A23);
//     qd[qd_stride * 4] = w * c[c_stride * 0] * (A21 * A31 + A22 * A32 + A23 * A33);
//     qd[qd_stride * 5] = w * c[c_stride * 0] * (A31 * A31 + A32 * A32 + A33 * A33);
//   }
// }

// CEED_QFUNCTION_HELPER void MultAdjJCAdjJt32(const CeedScalar *J, const CeedInt J_stride,
//                                             const CeedScalar *c, const CeedInt c_stride,
//                                             const CeedInt c_comp, const CeedScalar qw,
//                                             const CeedInt qd_stride, CeedScalar *qd)
// {
//   // Compute qw / det(J) adj(J) C adj(J)^T and store the symmetric part of the result.
//   // J: 0 3   qd: 0 1
//   //    1 4       1 2
//   //    2 5
//   const CeedScalar J11 = J[J_stride * 0];
//   const CeedScalar J21 = J[J_stride * 1];
//   const CeedScalar J31 = J[J_stride * 2];
//   const CeedScalar J12 = J[J_stride * 3];
//   const CeedScalar J22 = J[J_stride * 4];
//   const CeedScalar J32 = J[J_stride * 5];
//   const CeedScalar E = J11 * J11 + J21 * J21 + J31 * J31;
//   const CeedScalar G = J12 * J12 + J22 * J22 + J32 * J32;
//   const CeedScalar F = J11 * J12 + J21 * J22 + J31 * J32;
//   const CeedScalar d = E * G - F * F;
//   const CeedScalar w = qw / sqrt(d);
//   if (c_comp == 6)  // Matrix coefficient (symmetric)
//   {
//     // First compute entries of R = C adj(J)^T.
//     // c: 0 1 2
//     //    1 3 4
//     //    2 4 5
//     const CeedScalar R11 =
//         G * (c[c_stride * 0] * J11 + c[c_stride * 1] * J21 + c[c_stride * 2] * J31) -
//         F * (c[c_stride * 0] * J12 + c[c_stride * 1] * J22 + c[c_stride * 2] * J32);
//     const CeedScalar R21 =
//         G * (c[c_stride * 1] * J11 + c[c_stride * 3] * J21 + c[c_stride * 4] * J31) -
//         F * (c[c_stride * 1] * J12 + c[c_stride * 3] * J22 + c[c_stride * 4] * J32);
//     const CeedScalar R31 =
//         G * (c[c_stride * 2] * J11 + c[c_stride * 4] * J21 + c[c_stride * 5] * J31) -
//         F * (c[c_stride * 2] * J12 + c[c_stride * 4] * J22 + c[c_stride * 5] * J32);
//     const CeedScalar R12 =
//         E * (c[c_stride * 0] * J12 + c[c_stride * 1] * J22 + c[c_stride * 2] * J32) -
//         F * (c[c_stride * 0] * J11 + c[c_stride * 1] * J21 + c[c_stride * 2] * J31);
//     const CeedScalar R22 =
//         E * (c[c_stride * 1] * J12 + c[c_stride * 3] * J22 + c[c_stride * 4] * J32) -
//         F * (c[c_stride * 1] * J11 + c[c_stride * 3] * J21 + c[c_stride * 4] * J31);
//     const CeedScalar R32 =
//         E * (c[c_stride * 2] * J12 + c[c_stride * 4] * J22 + c[c_stride * 5] * J32) -
//         F * (c[c_stride * 2] * J11 + c[c_stride * 4] * J21 + c[c_stride * 5] * J31);
//     qd[qd_stride * 0] = w *
//                         (G * (J11 * R11 + J21 * R21 + J31 * R31) -
//                          F * (J12 * R11 + J22 * R21 + J32 * R31)) /
//                         d;
//     qd[qd_stride * 1] = w *
//                         (G * (J11 * R12 + J21 * R22 + J31 * R32) -
//                          F * (J12 * R12 + J22 * R22 + J32 * R32)) /
//                         d;
//     qd[qd_stride * 2] = w *
//                         (E * (J12 * R12 + J22 * R22 + J32 * R32) -
//                          F * (J11 * R12 + J21 * R22 + J31 * R32)) /
//                         d;
//   }
//   else if (c_comp == 3)  // Vector coefficient
//   {
//     // First compute entries of R = C adj(J)^T.
//     // c: 0
//     //      1
//     //        2
//     const CeedScalar R11 = c[c_stride * 0] * (G * J11 - F * J12);
//     const CeedScalar R21 = c[c_stride * 1] * (G * J21 - F * J22);
//     const CeedScalar R31 = c[c_stride * 2] * (G * J31 - F * J32);
//     const CeedScalar R12 = c[c_stride * 0] * (E * J12 - F * J11);
//     const CeedScalar R22 = c[c_stride * 1] * (E * J22 - F * J21);
//     const CeedScalar R32 = c[c_stride * 2] * (E * J32 - F * J31);
//     qd[qd_stride * 0] = w *
//                         (G * (J11 * R11 + J21 * R21 + J31 * R31) -
//                          F * (J12 * R11 + J22 * R21 + J32 * R31)) /
//                         d;
//     qd[qd_stride * 1] = w *
//                         (G * (J11 * R12 + J21 * R22 + J31 * R32) -
//                          F * (J12 * R12 + J22 * R22 + J32 * R32)) /
//                         d;
//     qd[qd_stride * 2] = w *
//                         (E * (J12 * R12 + J22 * R22 + J32 * R32) -
//                          F * (J11 * R12 + J21 * R22 + J31 * R32)) /
//                         d;
//   }
//   else  // Scalar coefficient
//   {
//     qd[qd_stride * 0] = w * c[c_stride * 0] * G;
//     qd[qd_stride * 1] = -w * c[c_stride * 0] * F;
//     qd[qd_stride * 2] = w * c[c_stride * 0] * E;
//   }
// }

// CEED_QFUNCTION_HELPER void MultJtCJ22(const CeedScalar *J, const CeedInt J_stride,
//                                       const CeedScalar *c, const CeedInt c_stride,
//                                       const CeedInt c_comp, const CeedScalar qw,
//                                       const CeedInt qd_stride, CeedScalar *qd)
// {
//   // Compute qw / det(J) J^T C J and store the symmetric part of the result.
//   // J: 0 2   qd: 0 1
//   //    1 3       1 2
//   const CeedScalar J11 = J[J_stride * 0];
//   const CeedScalar J21 = J[J_stride * 1];
//   const CeedScalar J12 = J[J_stride * 2];
//   const CeedScalar J22 = J[J_stride * 3];
//   const CeedScalar w = qw / (J11 * J22 - J21 * J12);
//   if (c_comp == 3)  // Matrix coefficient (symmetric)
//   {
//     // First compute entries of R = C J.
//     // c: 0 1
//     //    1 2
//     const CeedScalar R11 = c[c_stride * 0] * J11 + c[c_stride * 1] * J21;
//     const CeedScalar R21 = c[c_stride * 1] * J11 + c[c_stride * 2] * J21;
//     const CeedScalar R12 = c[c_stride * 0] * J12 + c[c_stride * 1] * J22;
//     const CeedScalar R22 = c[c_stride * 1] * J12 + c[c_stride * 2] * J22;
//     qd[qd_stride * 0] = w * (J11 * R11 + J21 * R21);
//     qd[qd_stride * 1] = w * (J11 * R12 + J21 * R22);
//     qd[qd_stride * 2] = w * (J12 * R12 + J22 * R22);
//   }
//   else if (c_comp == 2)  // Vector coefficient
//   {
//     // c: 0
//     //      1
//     qd[qd_stride * 0] = w * (c[c_stride * 0] * J11 * J11 + c[c_stride * 1] * J21 * J21);
//     qd[qd_stride * 1] = w * (c[c_stride * 0] * J11 * J12 + c[c_stride * 1] * J21 * J22);
//     qd[qd_stride * 2] = w * (c[c_stride * 0] * J12 * J12 + c[c_stride * 1] * J22 * J22);
//   }
//   else  // Scalar coefficient
//   {
//     qd[qd_stride * 0] = w * c[c_stride * 0] * (J11 * J11 + J21 * J21);
//     qd[qd_stride * 1] = w * c[c_stride * 0] * (J11 * J12 + J21 * J22);
//     qd[qd_stride * 2] = w * c[c_stride * 0] * (J12 * J12 + J22 * J22);
//   }
// }

// CEED_QFUNCTION_HELPER void MultJtCJ21(const CeedScalar *J, const CeedInt J_stride,
//                                       const CeedScalar *c, const CeedInt c_stride,
//                                       const CeedInt c_comp, const CeedScalar qw,
//                                       const CeedInt qd_stride, CeedScalar *qd)
// {
//   // Compute qw / det(J) J^T C J and store the symmetric part of the result.
//   // J: 0   qd: 0
//   //    1
//   const CeedScalar J11 = J[J_stride * 0];
//   const CeedScalar J21 = J[J_stride * 1];
//   if (c_comp == 3)  // Matrix coefficient (symmetric)
//   {
//     // First compute entries of R = C J.
//     // c: 0 1
//     //    1 2
//     const CeedScalar w = qw / sqrt(J11 * J11 + J21 * J21);
//     const CeedScalar R11 = c[c_stride * 0] * J11 + c[c_stride * 1] * J21;
//     const CeedScalar R21 = c[c_stride * 1] * J11 + c[c_stride * 2] * J21;
//     qd[qd_stride * 0] = w * (J11 * R11 + J21 * R21);
//   }
//   else if (c_comp == 2)  // Vector coefficient
//   {
//     // c: 0
//     //      1
//     const CeedScalar w = qw / sqrt(J11 * J11 + J21 * J21);
//     qd[qd_stride * 0] = w * (c[c_stride * 0] * J11 * J11 + c[c_stride * 1] * J21 * J21);
//   }
//   else  // Scalar coefficient
//   {
//     qd[qd_stride * 0] = qw * c[c_stride * 0] * sqrt(J11 * J11 + J21 * J21);
//   }
// }

// CEED_QFUNCTION_HELPER void MultJtCJ33(const CeedScalar *J, const CeedInt J_stride,
//                                       const CeedScalar *c, const CeedInt c_stride,
//                                       const CeedInt c_comp, const CeedScalar qw,
//                                       const CeedInt qd_stride, CeedScalar *qd)
// {
//   // Compute qw / det(J) J^T C J and store the symmetric part of the result.
//   // J: 0 3 6   qd: 0 1 2
//   //    1 4 7       1 3 4
//   //    2 5 8       2 4 5
//   const CeedScalar J11 = J[J_stride * 0];
//   const CeedScalar J21 = J[J_stride * 1];
//   const CeedScalar J31 = J[J_stride * 2];
//   const CeedScalar J12 = J[J_stride * 3];
//   const CeedScalar J22 = J[J_stride * 4];
//   const CeedScalar J32 = J[J_stride * 5];
//   const CeedScalar J13 = J[J_stride * 6];
//   const CeedScalar J23 = J[J_stride * 7];
//   const CeedScalar J33 = J[J_stride * 8];
//   const CeedScalar w = qw / (J11 * (J22 * J33 - J23 * J32) + J21 * (J13 * J32 - J12 *
//   J33) +
//                              J31 * (J12 * J23 - J13 * J22));
//   if (c_comp == 6)  // Matrix coefficient (symmetric)
//   {
//     // First compute entries of R = C J.
//     // c: 0 1 2
//     //    1 3 4
//     //    2 4 5
//     const CeedScalar R11 =
//         c[c_stride * 0] * J11 + c[c_stride * 1] * J21 + c[c_stride * 2] * J31;
//     const CeedScalar R21 =
//         c[c_stride * 1] * J11 + c[c_stride * 3] * J21 + c[c_stride * 4] * J31;
//     const CeedScalar R31 =
//         c[c_stride * 2] * J11 + c[c_stride * 4] * J21 + c[c_stride * 5] * J31;
//     const CeedScalar R12 =
//         c[c_stride * 0] * J12 + c[c_stride * 1] * J22 + c[c_stride * 2] * J32;
//     const CeedScalar R22 =
//         c[c_stride * 1] * J12 + c[c_stride * 3] * J22 + c[c_stride * 4] * J32;
//     const CeedScalar R32 =
//         c[c_stride * 2] * J12 + c[c_stride * 4] * J22 + c[c_stride * 5] * J32;
//     const CeedScalar R13 =
//         c[c_stride * 0] * J13 + c[c_stride * 1] * J23 + c[c_stride * 2] * J33;
//     const CeedScalar R23 =
//         c[c_stride * 1] * J13 + c[c_stride * 3] * J23 + c[c_stride * 4] * J33;
//     const CeedScalar R33 =
//         c[c_stride * 2] * J13 + c[c_stride * 4] * J23 + c[c_stride * 5] * J33;
//     qd[qd_stride * 0] = w * (J11 * R11 + J21 * R21 + J31 * R31);
//     qd[qd_stride * 1] = w * (J11 * R12 + J21 * R22 + J31 * R32);
//     qd[qd_stride * 2] = w * (J11 * R13 + J21 * R23 + J31 * R33);
//     qd[qd_stride * 3] = w * (J12 * R12 + J22 * R22 + J32 * R32);
//     qd[qd_stride * 4] = w * (J12 * R13 + J22 * R23 + J32 * R33);
//     qd[qd_stride * 5] = w * (J13 * R13 + J23 * R23 + J33 * R33);
//   }
//   else if (c_comp == 3)  // Vector coefficient
//   {
//     // c: 0
//     //      1
//     //        2
//     qd[qd_stride * 0] = w * (c[c_stride * 0] * J11 * J11 + c[c_stride * 1] * J21 * J21 +
//                              c[c_stride * 2] * J31 * J31);
//     qd[qd_stride * 1] = w * (c[c_stride * 0] * J11 * J12 + c[c_stride * 1] * J21 * J22 +
//                              c[c_stride * 2] * J31 * J32);
//     qd[qd_stride * 2] = w * (c[c_stride * 0] * J11 * J13 + c[c_stride * 1] * J21 * J23 +
//                              c[c_stride * 2] * J31 * J33);
//     qd[qd_stride * 3] = w * (c[c_stride * 0] * J12 * J12 + c[c_stride * 1] * J22 * J22 +
//                              c[c_stride * 2] * J32 * J32);
//     qd[qd_stride * 4] = w * (c[c_stride * 0] * J12 * J13 + c[c_stride * 1] * J22 * J23 +
//                              c[c_stride * 2] * J32 * J33);
//     qd[qd_stride * 5] = w * (c[c_stride * 0] * J13 * J13 + c[c_stride * 1] * J23 * J23 +
//                              c[c_stride * 2] * J33 * J33);
//   }
//   else  // Scalar coefficient
//   {
//     qd[qd_stride * 0] = w * c[c_stride * 0] * (J11 * J11 + J21 * J21 + J31 * J31);
//     qd[qd_stride * 1] = w * c[c_stride * 0] * (J11 * J12 + J21 * J22 + J31 * J32);
//     qd[qd_stride * 2] = w * c[c_stride * 0] * (J11 * J13 + J21 * J23 + J31 * J33);
//     qd[qd_stride * 3] = w * c[c_stride * 0] * (J12 * J12 + J22 * J22 + J32 * J32);
//     qd[qd_stride * 4] = w * c[c_stride * 0] * (J12 * J13 + J22 * J23 + J32 * J33);
//     qd[qd_stride * 5] = w * c[c_stride * 0] * (J13 * J13 + J23 * J23 + J33 * J33);
//   }
// }

// CEED_QFUNCTION_HELPER void MultJtCJ32(const CeedScalar *J, const CeedInt J_stride,
//                                       const CeedScalar *c, const CeedInt c_stride,
//                                       const CeedInt c_comp, const CeedScalar qw,
//                                       const CeedInt qd_stride, CeedScalar *qd)
// {
//   // Compute qw / det(J) J^T C J and store the symmetric part of the result.
//   // J: 0 3   qd: 0 1
//   //    1 4       1 2
//   //    2 5
//   const CeedScalar J11 = J[J_stride * 0];
//   const CeedScalar J21 = J[J_stride * 1];
//   const CeedScalar J31 = J[J_stride * 2];
//   const CeedScalar J12 = J[J_stride * 3];
//   const CeedScalar J22 = J[J_stride * 4];
//   const CeedScalar J32 = J[J_stride * 5];
//   const CeedScalar E = J11 * J11 + J21 * J21 + J31 * J31;
//   const CeedScalar G = J12 * J12 + J22 * J22 + J32 * J32;
//   const CeedScalar F = J11 * J12 + J21 * J22 + J31 * J32;
//   const CeedScalar w = qw / sqrt(E * G - F * F);
//   if (c_comp == 6)  // Matrix coefficient (symmetric)
//   {
//     // First compute entries of R = C J.
//     // c: 0 1 2
//     //    1 3 4
//     //    2 4 5
//     const CeedScalar R11 =
//         c[c_stride * 0] * J11 + c[c_stride * 1] * J21 + c[c_stride * 2] * J31;
//     const CeedScalar R21 =
//         c[c_stride * 1] * J11 + c[c_stride * 3] * J21 + c[c_stride * 4] * J31;
//     const CeedScalar R31 =
//         c[c_stride * 2] * J11 + c[c_stride * 4] * J21 + c[c_stride * 5] * J31;
//     const CeedScalar R12 =
//         c[c_stride * 0] * J12 + c[c_stride * 1] * J22 + c[c_stride * 2] * J32;
//     const CeedScalar R22 =
//         c[c_stride * 1] * J12 + c[c_stride * 3] * J22 + c[c_stride * 4] * J32;
//     const CeedScalar R32 =
//         c[c_stride * 2] * J12 + c[c_stride * 4] * J22 + c[c_stride * 5] * J32;
//     qd[qd_stride * 0] = w * (J11 * R11 + J21 * R21 + J31 * R31);
//     qd[qd_stride * 1] = w * (J11 * R12 + J21 * R22 + J31 * R32);
//     qd[qd_stride * 2] = w * (J12 * R12 + J22 * R22 + J32 * R32);
//   }
//   else if (c_comp == 3)  // Vector coefficient
//   {
//     // c: 0
//     //      1
//     //        2
//     qd[qd_stride * 0] = w * (c[c_stride * 0] * J11 * J11 + c[c_stride * 1] * J21 * J21 +
//                              c[c_stride * 2] * J31 * J31);
//     qd[qd_stride * 1] = w * (c[c_stride * 0] * J11 * J12 + c[c_stride * 1] * J21 * J22 +
//                              c[c_stride * 2] * J31 * J32);
//     qd[qd_stride * 2] = w * (c[c_stride * 0] * J12 * J12 + c[c_stride * 1] * J22 * J22 +
//                              c[c_stride * 2] * J32 * J32);
//   }
//   else  // Scalar coefficient
//   {
//     qd[qd_stride * 0] = w * c[c_stride * 0] * E;
//     qd[qd_stride * 1] = w * c[c_stride * 0] * F;
//     qd[qd_stride * 2] = w * c[c_stride * 0] * G;
//   }
// }

// template <bool Transpose = false>
// CEED_QFUNCTION_HELPER void MultJtCAdjJt22(const CeedScalar *J, const CeedInt J_stride,
//                                           const CeedScalar *c, const CeedInt c_stride,
//                                           const CeedInt c_comp, const CeedScalar qw,
//                                           const CeedInt qd_stride, CeedScalar *qd)
// {
//   // Compute qw / det(J) J^T C adj(J)^T and store the result.
//   // J: 0 2   adj(J):  J22 -J12   qd: 0 2
//   //    1 3           -J21  J11       1 3
//   const CeedScalar J11 = J[J_stride * 0];
//   const CeedScalar J21 = J[J_stride * 1];
//   const CeedScalar J12 = J[J_stride * 2];
//   const CeedScalar J22 = J[J_stride * 3];
//   const CeedScalar w = qw / (J11 * J22 - J21 * J12);
//   if (c_comp == 3)  // Matrix coefficient (symmetric)
//   {
//     // First compute entries of R = C adj(J)^T.
//     // c: 0 1
//     //    1 2
//     const CeedScalar R11 = c[c_stride * 0] * J22 - c[c_stride * 1] * J12;
//     const CeedScalar R21 = c[c_stride * 1] * J22 - c[c_stride * 2] * J12;
//     const CeedScalar R12 = -c[c_stride * 0] * J21 + c[c_stride * 1] * J11;
//     const CeedScalar R22 = -c[c_stride * 1] * J21 + c[c_stride * 2] * J11;
//     qd[qd_stride * 0] = w * (J11 * R11 + J21 * R21);
//     qd[qd_stride * 1] = w * (J12 * R11 + J22 * R21);
//     qd[qd_stride * 2] = w * (J11 * R12 + J21 * R22);
//     qd[qd_stride * 3] = w * (J12 * R12 + J22 * R22);
//   }
//   else if (c_comp == 2)  // Vector coefficient
//   {
//     // c: 0
//     //      1
//     qd[qd_stride * 0] = w * (c[c_stride * 0] * J11 * J22 - c[c_stride * 1] * J12 * J21);
//     qd[qd_stride * 1] = w * (c[c_stride * 0] * J12 * J22 - c[c_stride * 1] * J12 * J22);
//     qd[qd_stride * 2] = w * (-c[c_stride * 0] * J11 * J21 + c[c_stride * 1] * J11 * J21);
//     qd[qd_stride * 3] = w * (-c[c_stride * 0] * J12 * J21 + c[c_stride * 1] * J11 * J22);
//   }
//   else  // Scalar coefficient
//   {
//     qd[qd_stride * 0] = qw * c[c_stride * 0];
//     qd[qd_stride * 1] = 0.0;
//     qd[qd_stride * 2] = 0.0;
//     qd[qd_stride * 3] = qw * c[c_stride * 0];
//   }
//   if (Transpose && c_comp > 1)
//   {
//     const CeedScalar qd21 = qd[qd_stride * 1];
//     qd[qd_stride * 1] = qd[qd_stride * 2];
//     qd[qd_stride * 2] = qd21;
//   }
// }

// template <bool Transpose = false>
// CEED_QFUNCTION_HELPER void MultJtCAdjJt21(const CeedScalar *J, const CeedInt J_stride,
//                                           const CeedScalar *c, const CeedInt c_stride,
//                                           const CeedInt c_comp, const CeedScalar qw,
//                                           const CeedInt qd_stride, CeedScalar *qd)
// {
//   // Compute qw / det(J) J^T C adj(J)^T and store the result.
//   // J: 0   adj(J): 1/sqrt(J^T J) J^T   qd: 0
//   //    1
//   const CeedScalar J11 = J[J_stride * 0];
//   const CeedScalar J21 = J[J_stride * 1];
//   const CeedScalar w = qw / (J11 * J11 + J21 * J21);
//   if (c_comp == 3)  // Matrix coefficient (symmetric)
//   {
//     // First compute entries of R = C adj(J)^T.
//     // c: 0 1
//     //    1 2
//     const CeedScalar R11 = c[c_stride * 0] * J11 + c[c_stride * 1] * J21;
//     const CeedScalar R21 = c[c_stride * 1] * J11 + c[c_stride * 2] * J21;
//     qd[qd_stride * 0] = w * (J11 * R11 + J21 * R21);
//   }
//   else if (c_comp == 2)  // Vector coefficient
//   {
//     // c: 0
//     //      1
//     qd[qd_stride * 0] = w * (c[c_stride * 0] * J11 * J11 + c[c_stride * 1] * J21 * J21);
//   }
//   else  // Scalar coefficient
//   {
//     qd[qd_stride * 0] = qw * c[c_stride * 0];
//   }
// }

// template <bool Transpose = false>
// CEED_QFUNCTION_HELPER void MultJtCAdjJt33(const CeedScalar *J, const CeedInt J_stride,
//                                           const CeedScalar *c, const CeedInt c_stride,
//                                           const CeedInt c_comp, const CeedScalar qw,
//                                           const CeedInt qd_stride, CeedScalar *qd)
// {
//   // Compute qw / det(J) J^T C adj(J)^T and store the result.
//   // J: 0 3 6   qd: 0 3 6
//   //    1 4 7       1 4 7
//   //    2 5 8       2 5 8
//   const CeedScalar J11 = J[J_stride * 0];
//   const CeedScalar J21 = J[J_stride * 1];
//   const CeedScalar J31 = J[J_stride * 2];
//   const CeedScalar J12 = J[J_stride * 3];
//   const CeedScalar J22 = J[J_stride * 4];
//   const CeedScalar J32 = J[J_stride * 5];
//   const CeedScalar J13 = J[J_stride * 6];
//   const CeedScalar J23 = J[J_stride * 7];
//   const CeedScalar J33 = J[J_stride * 8];
//   const CeedScalar A11 = J22 * J33 - J23 * J32;
//   const CeedScalar A21 = J23 * J31 - J21 * J33;
//   const CeedScalar A31 = J21 * J32 - J22 * J31;
//   const CeedScalar A12 = J13 * J32 - J12 * J33;
//   const CeedScalar A22 = J11 * J33 - J13 * J31;
//   const CeedScalar A32 = J12 * J31 - J11 * J32;
//   const CeedScalar A13 = J12 * J23 - J13 * J22;
//   const CeedScalar A23 = J13 * J21 - J11 * J23;
//   const CeedScalar A33 = J11 * J22 - J12 * J21;
//   const CeedScalar w = qw / (J11 * A11 + J21 * A12 + J31 * A13);
//   if (c_comp == 6)  // Matrix coefficient (symmetric)
//   {
//     // First compute entries of R = C adj(J)^T.
//     // c: 0 1 2
//     //    1 3 4
//     //    2 4 5
//     const CeedScalar R11 =
//         c[c_stride * 0] * A11 + c[c_stride * 1] * A12 + c[c_stride * 2] * A13;
//     const CeedScalar R21 =
//         c[c_stride * 1] * A11 + c[c_stride * 3] * A12 + c[c_stride * 4] * A13;
//     const CeedScalar R31 =
//         c[c_stride * 2] * A11 + c[c_stride * 4] * A12 + c[c_stride * 5] * A13;
//     const CeedScalar R12 =
//         c[c_stride * 0] * A21 + c[c_stride * 1] * A22 + c[c_stride * 2] * A23;
//     const CeedScalar R22 =
//         c[c_stride * 1] * A21 + c[c_stride * 3] * A22 + c[c_stride * 4] * A23;
//     const CeedScalar R32 =
//         c[c_stride * 2] * A21 + c[c_stride * 4] * A22 + c[c_stride * 5] * A23;
//     const CeedScalar R13 =
//         c[c_stride * 0] * A31 + c[c_stride * 1] * A32 + c[c_stride * 2] * A33;
//     const CeedScalar R23 =
//         c[c_stride * 1] * A31 + c[c_stride * 3] * A32 + c[c_stride * 4] * A33;
//     const CeedScalar R33 =
//         c[c_stride * 2] * A31 + c[c_stride * 4] * A32 + c[c_stride * 5] * A33;
//     qd[qd_stride * 0] = w * (J11 * R11 + J21 * R21 + J31 * R31);
//     qd[qd_stride * 1] = w * (J12 * R11 + J22 * R21 + J32 * R31);
//     qd[qd_stride * 2] = w * (J13 * R11 + J23 * R21 + J33 * R31);
//     qd[qd_stride * 3] = w * (J11 * R12 + J21 * R22 + J31 * R32);
//     qd[qd_stride * 4] = w * (J12 * R12 + J22 * R22 + J32 * R32);
//     qd[qd_stride * 5] = w * (J13 * R12 + J23 * R22 + J33 * R32);
//     qd[qd_stride * 6] = w * (J11 * R13 + J21 * R23 + J31 * R33);
//     qd[qd_stride * 7] = w * (J12 * R13 + J22 * R23 + J32 * R33);
//     qd[qd_stride * 8] = w * (J13 * R13 + J23 * R23 + J33 * R33);
//   }
//   else if (c_comp == 3)  // Vector coefficient
//   {
//     // c: 0
//     //      1
//     //        2
//     qd[qd_stride * 0] = w * (c[c_stride * 0] * A11 * J11 + c[c_stride * 1] * A12 * J21 +
//                              c[c_stride * 2] * A13 * J31);
//     qd[qd_stride * 1] = w * (c[c_stride * 0] * A11 * J12 + c[c_stride * 1] * A12 * J22 +
//                              c[c_stride * 2] * A13 * J32);
//     qd[qd_stride * 2] = w * (c[c_stride * 0] * A11 * J13 + c[c_stride * 1] * A12 * J23 +
//                              c[c_stride * 2] * A13 * J33);
//     qd[qd_stride * 3] = w * (c[c_stride * 0] * A21 * J11 + c[c_stride * 1] * A22 * J21 +
//                              c[c_stride * 2] * A23 * J31);
//     qd[qd_stride * 4] = w * (c[c_stride * 0] * A21 * J12 + c[c_stride * 1] * A22 * J22 +
//                              c[c_stride * 2] * A23 * J32);
//     qd[qd_stride * 5] = w * (c[c_stride * 0] * A21 * J13 + c[c_stride * 1] * A22 * J23 +
//                              c[c_stride * 2] * A23 * J33);
//     qd[qd_stride * 6] = w * (c[c_stride * 0] * A31 * J11 + c[c_stride * 1] * A32 * J21 +
//                              c[c_stride * 2] * A33 * J31);
//     qd[qd_stride * 7] = w * (c[c_stride * 0] * A31 * J12 + c[c_stride * 1] * A32 * J22 +
//                              c[c_stride * 2] * A33 * J32);
//     qd[qd_stride * 8] = w * (c[c_stride * 0] * A31 * J13 + c[c_stride * 1] * A32 * J23 +
//                              c[c_stride * 2] * A33 * J33);
//   }
//   else  // Scalar coefficient
//   {
//     qd[qd_stride * 0] = qw * c[c_stride * 0];
//     qd[qd_stride * 1] = 0.0;
//     qd[qd_stride * 2] = 0.0;
//     qd[qd_stride * 3] = 0.0;
//     qd[qd_stride * 4] = qw * c[c_stride * 0];
//     qd[qd_stride * 5] = 0.0;
//     qd[qd_stride * 6] = 0.0;
//     qd[qd_stride * 7] = 0.0;
//     qd[qd_stride * 8] = qw * c[c_stride * 0];
//   }
//   if (Transpose && c_comp > 1)
//   {
//     {
//       const CeedScalar qd21 = qd[qd_stride * 1];
//       qd[qd_stride * 1] = qd[qd_stride * 3];
//       qd[qd_stride * 3] = qd21;
//     }
//     {
//       const CeedScalar qd31 = qd[qd_stride * 2];
//       qd[qd_stride * 2] = qd[qd_stride * 6];
//       qd[qd_stride * 6] = qd31;
//     }
//     {
//       const CeedScalar qd32 = qd[qd_stride * 5];
//       qd[qd_stride * 5] = qd[qd_stride * 7];
//       qd[qd_stride * 7] = qd32;
//     }
//   }
// }

// template <bool Transpose = false>
// CEED_QFUNCTION_HELPER void MultJtCAdjJt32(const CeedScalar *J, const CeedInt J_stride,
//                                           const CeedScalar *c, const CeedInt c_stride,
//                                           const CeedInt c_comp, const CeedScalar qw,
//                                           const CeedInt qd_stride, CeedScalar *qd)
// {
//   // Compute qw / det(J) J^T C adj(J)^T and store the result.
//   // J: 0 3   qd: 0 2
//   //    1 4       1 3
//   //    2 5
//   const CeedScalar J11 = J[J_stride * 0];
//   const CeedScalar J21 = J[J_stride * 1];
//   const CeedScalar J31 = J[J_stride * 2];
//   const CeedScalar J12 = J[J_stride * 3];
//   const CeedScalar J22 = J[J_stride * 4];
//   const CeedScalar J32 = J[J_stride * 5];
//   const CeedScalar E = J11 * J11 + J21 * J21 + J31 * J31;
//   const CeedScalar G = J12 * J12 + J22 * J22 + J32 * J32;
//   const CeedScalar F = J11 * J12 + J21 * J22 + J31 * J32;
//   const CeedScalar w = qw / (E * G - F * F);
//   if (c_comp == 6)  // Matrix coefficient (symmetric)
//   {
//     // First compute entries of R = C adj(J)^T.
//     // c: 0 1 2
//     //    1 3 4
//     //    2 4 5
//     const CeedScalar R11 =
//         G * (c[c_stride * 0] * J11 + c[c_stride * 1] * J21 + c[c_stride * 2] * J31) -
//         F * (c[c_stride * 0] * J12 + c[c_stride * 1] * J22 + c[c_stride * 2] * J32);
//     const CeedScalar R21 =
//         G * (c[c_stride * 1] * J11 + c[c_stride * 3] * J21 + c[c_stride * 4] * J31) -
//         F * (c[c_stride * 1] * J12 + c[c_stride * 3] * J22 + c[c_stride * 4] * J32);
//     const CeedScalar R31 =
//         G * (c[c_stride * 2] * J11 + c[c_stride * 4] * J21 + c[c_stride * 5] * J31) -
//         F * (c[c_stride * 2] * J12 + c[c_stride * 4] * J22 + c[c_stride * 5] * J32);
//     const CeedScalar R12 =
//         E * (c[c_stride * 0] * J12 + c[c_stride * 1] * J22 + c[c_stride * 2] * J32) -
//         F * (c[c_stride * 0] * J11 + c[c_stride * 1] * J21 + c[c_stride * 2] * J31);
//     const CeedScalar R22 =
//         E * (c[c_stride * 1] * J12 + c[c_stride * 3] * J22 + c[c_stride * 4] * J32) -
//         F * (c[c_stride * 1] * J11 + c[c_stride * 3] * J21 + c[c_stride * 4] * J31);
//     const CeedScalar R32 =
//         E * (c[c_stride * 2] * J12 + c[c_stride * 4] * J22 + c[c_stride * 5] * J32) -
//         F * (c[c_stride * 2] * J11 + c[c_stride * 4] * J21 + c[c_stride * 5] * J31);
//     qd[qd_stride * 0] = w * (J11 * R11 + J21 * R21 + J31 * R31);
//     qd[qd_stride * 1] = w * (J12 * R11 + J22 * R21 + J32 * R31);
//     qd[qd_stride * 2] = w * (J11 * R12 + J21 * R22 + J31 * R32);
//     qd[qd_stride * 3] = w * (J12 * R12 + J22 * R22 + J32 * R32);
//   }
//   else if (c_comp == 3)  // Vector coefficient
//   {
//     // First compute entries of R = C adj(J)^T.
//     // c: 0
//     //      1
//     //        2
//     const CeedScalar R11 = c[c_stride * 0] * (G * J11 - F * J12);
//     const CeedScalar R21 = c[c_stride * 1] * (G * J21 - F * J22);
//     const CeedScalar R31 = c[c_stride * 2] * (G * J31 - F * J32);
//     const CeedScalar R12 = c[c_stride * 0] * (E * J12 - F * J11);
//     const CeedScalar R22 = c[c_stride * 1] * (E * J22 - F * J21);
//     const CeedScalar R32 = c[c_stride * 2] * (E * J32 - F * J31);
//     qd[qd_stride * 0] = w * (J11 * R11 + J21 * R21 + J31 * R31);
//     qd[qd_stride * 1] = w * (J12 * R11 + J22 * R21 + J32 * R31);
//     qd[qd_stride * 2] = w * (J11 * R12 + J21 * R22 + J31 * R32);
//     qd[qd_stride * 3] = w * (J12 * R12 + J22 * R22 + J32 * R32);
//   }
//   else  // Scalar coefficient
//   {
//     qd[qd_stride * 0] = qw * c[c_stride * 0];
//     qd[qd_stride * 1] = 0.0;
//     qd[qd_stride * 2] = 0.0;
//     qd[qd_stride * 3] = qw * c[c_stride * 0];
//   }
//   if (Transpose && c_comp > 1)
//   {
//     const CeedScalar qd21 = qd[qd_stride * 1];
//     qd[qd_stride * 1] = qd[qd_stride * 2];
//     qd[qd_stride * 2] = qd21;
//   }
// }

// CEED_QFUNCTION_HELPER void MultCAdjJt22(const CeedScalar *J, const CeedInt J_stride,
//                                         const CeedScalar *c, const CeedInt c_stride,
//                                         const CeedInt c_comp, const CeedScalar qw,
//                                         const CeedInt qd_stride, CeedScalar *qd)
// {
//   // Compute qw C adj(J)^T and store the result.
//   // J: 0 2   adj(J):  J22 -J12   qd: 0 2
//   //    1 3           -J21  J11       1 3
//   const CeedScalar J11 = J[J_stride * 0];
//   const CeedScalar J21 = J[J_stride * 1];
//   const CeedScalar J12 = J[J_stride * 2];
//   const CeedScalar J22 = J[J_stride * 3];
//   if (c_comp == 3)  // Matrix coefficient (symmetric)
//   {
//     // c: 0 1
//     //    1 2
//     qd[qd_stride * 0] = qw * (c[c_stride * 0] * J22 - c[c_stride * 1] * J12);
//     qd[qd_stride * 1] = qw * (c[c_stride * 1] * J22 - c[c_stride * 2] * J12);
//     qd[qd_stride * 2] = qw * (-c[c_stride * 0] * J21 + c[c_stride * 1] * J11);
//     qd[qd_stride * 3] = qw * (-c[c_stride * 1] * J21 + c[c_stride * 2] * J11);
//   }
//   else if (c_comp == 2)  // Vector coefficient
//   {
//     // c: 0
//     //      1
//     const CeedScalar wc0 = qw * c[c_stride * 0];
//     const CeedScalar wc1 = qw * c[c_stride * 1];
//     qd[qd_stride * 0] = wc0 * J22;
//     qd[qd_stride * 1] = -wc1 * J12;
//     qd[qd_stride * 2] = -wc0 * J21;
//     qd[qd_stride * 3] = wc1 * J11;
//   }
//   else  // Scalar coefficient
//   {
//     const CeedScalar wc = qw * c[c_stride * 0];
//     qd[qd_stride * 0] = wc * J22;
//     qd[qd_stride * 1] = -wc * J12;
//     qd[qd_stride * 2] = -wc * J21;
//     qd[qd_stride * 3] = wc * J11;
//   }
// }

// CEED_QFUNCTION_HELPER void MultCAdjJt21(const CeedScalar *J, const CeedInt J_stride,
//                                         const CeedScalar *c, const CeedInt c_stride,
//                                         const CeedInt c_comp, const CeedScalar qw,
//                                         const CeedInt qd_stride, CeedScalar *qd)
// {
//   // Compute qw C adj(J)^T and store the result.
//   // J: 0   adj(J): 1/sqrt(J^T J) J^T   qd: 0
//   //    1                                   1
//   const CeedScalar J11 = J[J_stride * 0];
//   const CeedScalar J21 = J[J_stride * 1];
//   const CeedScalar w = qw / sqrt(J11 * J11 + J21 * J21);
//   if (c_comp == 3)  // Matrix coefficient (symmetric)
//   {
//     // c: 0 1
//     //    1 2
//     qd[qd_stride * 0] = w * (c[c_stride * 0] * J11 + c[c_stride * 1] * J21);
//     qd[qd_stride * 1] = w * (c[c_stride * 1] * J11 + c[c_stride * 2] * J21);
//   }
//   else if (c_comp == 2)  // Vector coefficient
//   {
//     // c: 0
//     //      1
//     qd[qd_stride * 0] = w * c[c_stride * 0] * J11;
//     qd[qd_stride * 1] = w * c[c_stride * 1] * J21;
//   }
//   else  // Scalar coefficient
//   {
//     const CeedScalar wc = w * c[c_stride * 0];
//     qd[qd_stride * 0] = wc * J11;
//     qd[qd_stride * 1] = wc * J21;
//   }
// }

// CEED_QFUNCTION_HELPER void MultCAdjJt33(const CeedScalar *J, const CeedInt J_stride,
//                                         const CeedScalar *c, const CeedInt c_stride,
//                                         const CeedInt c_comp, const CeedScalar qw,
//                                         const CeedInt qd_stride, CeedScalar *qd)
// {
//   // Compute qw C adj(J)^T and store the result.
//   // J: 0 3 6   qd: 0 3 6
//   //    1 4 7       1 4 7
//   //    2 5 8       2 5 8
//   const CeedScalar J11 = J[J_stride * 0];
//   const CeedScalar J21 = J[J_stride * 1];
//   const CeedScalar J31 = J[J_stride * 2];
//   const CeedScalar J12 = J[J_stride * 3];
//   const CeedScalar J22 = J[J_stride * 4];
//   const CeedScalar J32 = J[J_stride * 5];
//   const CeedScalar J13 = J[J_stride * 6];
//   const CeedScalar J23 = J[J_stride * 7];
//   const CeedScalar J33 = J[J_stride * 8];
//   const CeedScalar A11 = J22 * J33 - J23 * J32;
//   const CeedScalar A21 = J23 * J31 - J21 * J33;
//   const CeedScalar A31 = J21 * J32 - J22 * J31;
//   const CeedScalar A12 = J13 * J32 - J12 * J33;
//   const CeedScalar A22 = J11 * J33 - J13 * J31;
//   const CeedScalar A32 = J12 * J31 - J11 * J32;
//   const CeedScalar A13 = J12 * J23 - J13 * J22;
//   const CeedScalar A23 = J13 * J21 - J11 * J23;
//   const CeedScalar A33 = J11 * J22 - J12 * J21;
//   if (c_comp == 6)  // Matrix coefficient (symmetric)
//   {
//     // c: 0 1 2
//     //    1 3 4
//     //    2 4 5
//     qd[qd_stride * 0] =
//         qw * (c[c_stride * 0] * A11 + c[c_stride * 1] * A12 + c[c_stride * 2] * A13);
//     qd[qd_stride * 1] =
//         qw * (c[c_stride * 1] * A11 + c[c_stride * 3] * A12 + c[c_stride * 4] * A13);
//     qd[qd_stride * 2] =
//         qw * (c[c_stride * 2] * A11 + c[c_stride * 4] * A12 + c[c_stride * 5] * A13);
//     qd[qd_stride * 3] =
//         qw * (c[c_stride * 0] * A21 + c[c_stride * 1] * A22 + c[c_stride * 2] * A23);
//     qd[qd_stride * 4] =
//         qw * (c[c_stride * 1] * A21 + c[c_stride * 3] * A22 + c[c_stride * 4] * A23);
//     qd[qd_stride * 5] =
//         qw * (c[c_stride * 2] * A21 + c[c_stride * 4] * A22 + c[c_stride * 5] * A23);
//     qd[qd_stride * 6] =
//         qw * (c[c_stride * 0] * A31 + c[c_stride * 1] * A32 + c[c_stride * 2] * A33);
//     qd[qd_stride * 7] =
//         qw * (c[c_stride * 1] * A31 + c[c_stride * 3] * A32 + c[c_stride * 4] * A33);
//     qd[qd_stride * 8] =
//         qw * (c[c_stride * 2] * A31 + c[c_stride * 4] * A32 + c[c_stride * 5] * A33);
//   }
//   else if (c_comp == 3)  // Vector coefficient
//   {
//     // c: 0
//     //      1
//     //        2
//     const CeedScalar wc0 = qw * c[c_stride * 0];
//     const CeedScalar wc1 = qw * c[c_stride * 1];
//     const CeedScalar wc2 = qw * c[c_stride * 2];
//     qd[qd_stride * 0] = wc0 * A11;
//     qd[qd_stride * 1] = wc1 * A12;
//     qd[qd_stride * 2] = wc2 * A13;
//     qd[qd_stride * 3] = wc0 * A21;
//     qd[qd_stride * 4] = wc1 * A22;
//     qd[qd_stride * 5] = wc2 * A23;
//     qd[qd_stride * 6] = wc0 * A31;
//     qd[qd_stride * 7] = wc1 * A32;
//     qd[qd_stride * 8] = wc2 * A33;
//   }
//   else  // Scalar coefficient
//   {
//     const CeedScalar wc = qw * c[c_stride * 0];
//     qd[qd_stride * 0] = wc * A11;
//     qd[qd_stride * 1] = wc * A12;
//     qd[qd_stride * 2] = wc * A13;
//     qd[qd_stride * 3] = wc * A21;
//     qd[qd_stride * 4] = wc * A22;
//     qd[qd_stride * 5] = wc * A23;
//     qd[qd_stride * 6] = wc * A31;
//     qd[qd_stride * 7] = wc * A32;
//     qd[qd_stride * 8] = wc * A33;
//   }
// }

// CEED_QFUNCTION_HELPER void MultCAdjJt32(const CeedScalar *J, const CeedInt J_stride,
//                                         const CeedScalar *c, const CeedInt c_stride,
//                                         const CeedInt c_comp, const CeedScalar qw,
//                                         const CeedInt qd_stride, CeedScalar *qd)
// {
//   // Compute qw C adj(J)^T and store the result.
//   // J: 0 3   qd: 0 3
//   //    1 4       1 4
//   //    2 5       2 5
//   const CeedScalar J11 = J[J_stride * 0];
//   const CeedScalar J21 = J[J_stride * 1];
//   const CeedScalar J31 = J[J_stride * 2];
//   const CeedScalar J12 = J[J_stride * 3];
//   const CeedScalar J22 = J[J_stride * 4];
//   const CeedScalar J32 = J[J_stride * 5];
//   const CeedScalar E = J11 * J11 + J21 * J21 + J31 * J31;
//   const CeedScalar G = J12 * J12 + J22 * J22 + J32 * J32;
//   const CeedScalar F = J11 * J12 + J21 * J22 + J31 * J32;
//   const CeedScalar w = qw / sqrt(E * G - F * F);
//   if (c_comp == 6)  // Matrix coefficient (symmetric)
//   {
//     // c: 0 1 2
//     //    1 3 4
//     //    2 4 5
//     qd[qd_stride * 0] =
//         w * (G * (c[c_stride * 0] * J11 + c[c_stride * 1] * J21 + c[c_stride * 2] * J31)
//         -
//              F * (c[c_stride * 0] * J12 + c[c_stride * 1] * J22 + c[c_stride * 2] *
//              J32));
//     qd[qd_stride * 1] =
//         w * (G * (c[c_stride * 1] * J11 + c[c_stride * 3] * J21 + c[c_stride * 4] * J31)
//         -
//              F * (c[c_stride * 1] * J12 + c[c_stride * 3] * J22 + c[c_stride * 4] *
//              J32));
//     qd[qd_stride * 2] =
//         w * (G * (c[c_stride * 2] * J11 + c[c_stride * 4] * J21 + c[c_stride * 5] * J31)
//         -
//              F * (c[c_stride * 2] * J12 + c[c_stride * 4] * J22 + c[c_stride * 5] *
//              J32));
//     qd[qd_stride * 3] =
//         w * (E * (c[c_stride * 0] * J12 + c[c_stride * 1] * J22 + c[c_stride * 2] * J32)
//         -
//              F * (c[c_stride * 0] * J11 + c[c_stride * 1] * J21 + c[c_stride * 2] *
//              J31));
//     qd[qd_stride * 4] =
//         w * (E * (c[c_stride * 1] * J12 + c[c_stride * 3] * J22 + c[c_stride * 4] * J32)
//         -
//              F * (c[c_stride * 1] * J11 + c[c_stride * 3] * J21 + c[c_stride * 4] *
//              J31));
//     qd[qd_stride * 5] =
//         w * (E * (c[c_stride * 2] * J12 + c[c_stride * 4] * J22 + c[c_stride * 5] * J32)
//         -
//              F * (c[c_stride * 2] * J11 + c[c_stride * 4] * J21 + c[c_stride * 5] *
//              J31));
//   }
//   else if (c_comp == 3)  // Vector coefficient
//   {
//     // c: 0
//     //      1
//     //        2
//     const CeedScalar wc0 = w * c[c_stride * 0];
//     const CeedScalar wc1 = w * c[c_stride * 1];
//     const CeedScalar wc2 = w * c[c_stride * 2];
//     qd[qd_stride * 0] = wc0 * (G * J11 - F * J12);
//     qd[qd_stride * 1] = wc1 * (G * J21 - F * J22);
//     qd[qd_stride * 2] = wc2 * (G * J31 - F * J32);
//     qd[qd_stride * 3] = wc0 * (E * J12 - F * J11);
//     qd[qd_stride * 4] = wc1 * (E * J22 - F * J21);
//     qd[qd_stride * 5] = wc2 * (E * J32 - F * J31);
//   }
//   else  // Scalar coefficient
//   {
//     const CeedScalar wc = w * c[c_stride * 0];
//     qd[qd_stride * 0] = wc * (G * J11 - F * J12);
//     qd[qd_stride * 1] = wc * (G * J21 - F * J22);
//     qd[qd_stride * 2] = wc * (G * J31 - F * J32);
//     qd[qd_stride * 3] = wc * (E * J12 - F * J11);
//     qd[qd_stride * 4] = wc * (E * J22 - F * J21);
//     qd[qd_stride * 5] = wc * (E * J32 - F * J31);
//   }
// }

#endif  // PALACE_LIBCEED_UTILS_QF_H
