// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <mfem.hpp>
#include <catch2/catch_test_macros.hpp>
#include "fem/substructure.hpp"
#include "utils/communication.hpp"

namespace palace
{

TEST_CASE("Interface true DOFs are identified consistently in parallel",
          "[substructure][Serial][Parallel]")
{
  // The global interface true-DOF count must be independent of the MPI partition. Cube
  // split at x=0.5, H1: 49 interface DOFs at order 1 (7x7 plane), 169 at order 2.
  auto count = [](int order)
  {
    mfem::Mesh serial = mfem::Mesh::MakeCartesian3D(6, 6, 6, mfem::Element::HEXAHEDRON);
    for (int e = 0; e < serial.GetNE(); e++)
    {
      mfem::Vector c;
      serial.GetElementCenter(e, c);
      serial.SetAttribute(e, (c(0) < 0.5) ? 1 : 2);
    }
    serial.SetAttributes();
    mfem::ParMesh mesh(Mpi::World(), serial);
    mfem::H1_FECollection fec(order, 3);
    mfem::ParFiniteElementSpace pfes(&mesh, &fec);
    mfem::Array<int> ra(1), ea(1);
    ra[0] = 1;
    ea[0] = 2;
    mfem::Array<int> rm, em, im;
    return MarkInterfaceTrueDofs(pfes, ra, ea, rm, em, im);
  };
  CHECK(count(1) == 49);
  CHECK(count(2) == 169);
}

}  // namespace palace
