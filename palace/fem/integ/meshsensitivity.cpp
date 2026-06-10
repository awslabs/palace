// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "meshsensitivity.hpp"

#include "fem/fespace.hpp"
#include "fem/libceed/ceed.hpp"
#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"
#include "fem/mesh.hpp"
#include "utils/omp.hpp"

#include "fem/qfunctions/meshsens_qf.h"

namespace palace
{

using namespace ceed;

void MeshSensitivityIntegrator::Apply(Mesh &mesh, const FiniteElementSpace &nd_fespace,
                                      const MaterialPropertyCoefficient *Q_mass,
                                      const MaterialPropertyCoefficient *Q_stiff,
                                      double lambda, const mfem::Vector &eigenvector,
                                      mfem::Vector &sensitivity)
{
  sensitivity = 0.0;

  const mfem::ParMesh &pmesh = mesh.Get();
  const mfem::GridFunction &mesh_nodes = *pmesh.GetNodes();
  const mfem::FiniteElementSpace &mesh_fespace = *mesh_nodes.FESpace();
  const int space_dim = pmesh.SpaceDimension();
  const int dim = pmesh.Dimension();

  MFEM_VERIFY(space_dim == 3 && dim == 3,
              "MeshSensitivityIntegrator only supports 3D (space_dim=dim=3)!");

  // Build the coefficient pair context (mass first, stiffness second) and append lambda.
  auto ctx = PopulateCoefficientContext(space_dim, Q_mass, dim, Q_stiff);
  CeedIntScalar lambda_entry;
  lambda_entry.second = lambda;
  ctx.push_back(lambda_entry);

  // Convert the true-DOF eigenvector to an L-vector (local, with shared DOFs duplicated)
  // for use with CeedOperator element restrictions.
  mfem::Vector eig_lvec(nd_fespace.GetVSize());
  const auto *P = nd_fespace.GetProlongationMatrix();
  if (P)
  {
    P->Mult(eigenvector, eig_lvec);
  }
  else
  {
    eig_lvec = eigenvector;
  }

  // Iterate over Ceed contexts (one per thread in OpenMP builds).
  PalacePragmaOmp(parallel if (ceed::internal::NumCeeds() > 1))
  {
    Ceed ceed = ceed::internal::GetCeedObjects()[utils::GetThreadNum()];

    for (const auto &[geom, geom_data] : mesh.GetCeedGeomFactorData(ceed))
    {
      // Only domain elements (dim == mesh dimension), not boundary.
      if (mfem::Geometry::Dimension[geom] != dim)
      {
        continue;
      }

      const auto &indices = geom_data.indices;

      // Build mesh node basis, restriction, and CeedVector.
      CeedElemRestriction mesh_restr =
          FiniteElementSpace::BuildCeedElemRestriction(mesh_fespace, ceed, geom, indices);
      CeedBasis mesh_basis = FiniteElementSpace::BuildCeedBasis(mesh_fespace, ceed, geom);
      CeedVector mesh_nodes_vec;
      ceed::InitCeedVector(mesh_nodes, ceed, &mesh_nodes_vec);

      // Build ND (H(curl)) solution basis and restriction.
      CeedBasis soln_basis = nd_fespace.GetCeedBasis(ceed, geom);
      CeedElemRestriction soln_restr =
          nd_fespace.GetCeedElemRestriction(ceed, geom, indices);

      // Wrap the L-vector eigenvector as a CeedVector.
      CeedVector eig_vec;
      ceed::InitCeedVector(eig_lvec, ceed, &eig_vec);

      // Build attribute basis and restriction (one constant value per element).
      const std::size_t num_elem = indices.size();
      CeedInt num_qpts;
      PalaceCeedCall(ceed, CeedBasisGetNumQuadraturePoints(mesh_basis, &num_qpts));

      CeedElemRestriction attr_restr;
      CeedBasis attr_basis;
      PalaceCeedCall(ceed,
                     CeedElemRestrictionCreateStrided(ceed, num_elem, 1, 1, num_elem,
                                                      CEED_STRIDES_BACKEND, &attr_restr));
      {
        mfem::Vector Bt(num_qpts), Gt(num_qpts), qX(num_qpts), qW(num_qpts);
        Bt = 1.0;
        Gt = 0.0;
        qX = 0.0;
        qW = 0.0;
        PalaceCeedCall(ceed, CeedBasisCreateH1(ceed, CEED_TOPOLOGY_LINE, 1, 1, num_qpts,
                                               Bt.GetData(), Gt.GetData(), qX.GetData(),
                                               qW.GetData(), &attr_basis));
      }

      // Build element attribute vector (libCEED-local attributes).
      mfem::Vector elem_attr_data(num_elem);
      for (std::size_t k = 0; k < num_elem; k++)
      {
        elem_attr_data[k] = static_cast<double>(
            mesh.GetCeedAttributes().at(pmesh.GetAttribute(indices[k])));
      }
      CeedVector elem_attr_vec;
      ceed::InitCeedVector(elem_attr_data, ceed, &elem_attr_vec);

      // Build sensitivity output CeedVector (same size as mesh nodes).
      CeedVector sens_vec;
      PalaceCeedCall(ceed, CeedVectorCreate(ceed, (CeedSize)mesh_nodes.Size(), &sens_vec));
      PalaceCeedCall(ceed, CeedVectorSetValue(sens_vec, 0.0));

      // Create QFunction.
      CeedQFunction qf;
      PalaceCeedCall(ceed, CeedQFunctionCreateInterior(
                               ceed, 1, f_apply_meshsens_33,
                               PalaceQFunctionRelativePath(f_apply_meshsens_33_loc), &qf));

      // Set coefficient context.
      CeedQFunctionContext qf_ctx;
      PalaceCeedCall(ceed, CeedQFunctionContextCreate(ceed, &qf_ctx));
      PalaceCeedCall(ceed, CeedQFunctionContextSetData(
                               qf_ctx, CEED_MEM_HOST, CEED_COPY_VALUES,
                               ctx.size() * sizeof(CeedIntScalar), ctx.data()));
      PalaceCeedCall(ceed, CeedQFunctionSetContext(qf, qf_ctx));
      PalaceCeedCall(ceed, CeedQFunctionContextDestroy(&qf_ctx));

      // QFunction inputs:
      //   [0] attr        : 1 component,  CEED_EVAL_INTERP  (passive)
      //   [1] q_w         : 1 component,  CEED_EVAL_WEIGHT  (passive)
      //   [2] grad_x / J  : 9 components, CEED_EVAL_GRAD    (active input)
      //   [3] u (interp)  : 3 components, CEED_EVAL_INTERP  (passive)
      //   [4] curl_u      : 3 components, CEED_EVAL_CURL    (passive)
      PalaceCeedCall(ceed, CeedQFunctionAddInput(qf, "attr", 1, CEED_EVAL_INTERP));
      PalaceCeedCall(ceed, CeedQFunctionAddInput(qf, "q_w", 1, CEED_EVAL_WEIGHT));
      PalaceCeedCall(ceed,
                     CeedQFunctionAddInput(qf, "grad_x", space_dim * dim, CEED_EVAL_GRAD));
      PalaceCeedCall(ceed, CeedQFunctionAddInput(qf, "u", space_dim, CEED_EVAL_INTERP));
      PalaceCeedCall(ceed, CeedQFunctionAddInput(qf, "curl_u", dim, CEED_EVAL_CURL));

      // QFunction output:
      //   [0] df_dJ : 9 components, CEED_EVAL_GRAD (active output)
      PalaceCeedCall(ceed,
                     CeedQFunctionAddOutput(qf, "df_dJ", space_dim * dim, CEED_EVAL_GRAD));

      // Create operator.
      CeedOperator op;
      PalaceCeedCall(ceed, CeedOperatorCreate(ceed, qf, nullptr, nullptr, &op));
      PalaceCeedCall(ceed, CeedQFunctionDestroy(&qf));

      // Set operator fields.
      PalaceCeedCall(
          ceed, CeedOperatorSetField(op, "attr", attr_restr, attr_basis, elem_attr_vec));
      PalaceCeedCall(ceed, CeedOperatorSetField(op, "q_w", CEED_ELEMRESTRICTION_NONE,
                                                mesh_basis, CEED_VECTOR_NONE));
      PalaceCeedCall(ceed, CeedOperatorSetField(op, "grad_x", mesh_restr, mesh_basis,
                                                CEED_VECTOR_ACTIVE));
      PalaceCeedCall(ceed, CeedOperatorSetField(op, "u", soln_restr, soln_basis, eig_vec));
      PalaceCeedCall(ceed,
                     CeedOperatorSetField(op, "curl_u", soln_restr, soln_basis, eig_vec));
      PalaceCeedCall(ceed, CeedOperatorSetField(op, "df_dJ", mesh_restr, mesh_basis,
                                                CEED_VECTOR_ACTIVE));

      PalaceCeedCall(ceed, CeedOperatorCheckReady(op));

      // Apply: mesh_nodes → sensitivity.
      PalaceCeedCall(
          ceed, CeedOperatorApply(op, mesh_nodes_vec, sens_vec, CEED_REQUEST_IMMEDIATE));

      // Accumulate into the output sensitivity vector. In OpenMP builds, multiple threads
      // may contribute to shared DOFs, so synchronize the accumulation.
      {
        const CeedScalar *sens_data;
        PalaceCeedCall(ceed, CeedVectorGetArrayRead(sens_vec, CEED_MEM_HOST, &sens_data));
        PalacePragmaOmp(critical(meshsens_accum))
        for (int i = 0; i < mesh_nodes.Size(); i++)
        {
          sensitivity(i) += sens_data[i];
        }
        PalaceCeedCall(ceed, CeedVectorRestoreArrayRead(sens_vec, &sens_data));
      }

      // Cleanup.
      PalaceCeedCall(ceed, CeedOperatorDestroy(&op));
      PalaceCeedCall(ceed, CeedVectorDestroy(&sens_vec));
      PalaceCeedCall(ceed, CeedVectorDestroy(&eig_vec));
      PalaceCeedCall(ceed, CeedVectorDestroy(&mesh_nodes_vec));
      PalaceCeedCall(ceed, CeedVectorDestroy(&elem_attr_vec));
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&mesh_restr));
      PalaceCeedCall(ceed, CeedBasisDestroy(&mesh_basis));
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&attr_restr));
      PalaceCeedCall(ceed, CeedBasisDestroy(&attr_basis));
    }
  }
}

}  // namespace palace
