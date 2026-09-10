// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Integrate physical measures of affine or curved meshes; do not substitute a
// center Jacobian when validating a curved geometry against CAD.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <stdexcept>
#include <vector>
#include <mfem.hpp>
#include <nlohmann/json.hpp>

using Measures = std::map<std::pair<int, int>, double>;

Measures Integrate(mfem::Mesh &mesh, int order)
{
  Measures result;
  for (int dimension : {2, 3})
  {
    int count = dimension == 3 ? mesh.GetNE() : mesh.GetNBE();
    for (int i = 0; i < count; i++)
    {
      auto *T = dimension == 3 ? mesh.GetElementTransformation(i)
                               : mesh.GetBdrElementTransformation(i);
      auto geometry = dimension == 3 ? mesh.GetElementBaseGeometry(i)
                                     : mesh.GetBdrElementGeometry(i);
      int attribute = dimension == 3 ? mesh.GetAttribute(i) : mesh.GetBdrAttribute(i);
      const auto &rule = mfem::IntRules.Get(geometry, order);
      double measure = 0.0;
      for (int q = 0; q < rule.GetNPoints(); q++)
      {
        const auto &point = rule.IntPoint(q);
        T->SetIntPoint(&point);
        if (dimension == 3 && !(T->Jacobian().Det() > 0.0))
        {
          throw std::runtime_error("Nonpositive quadrature-point Jacobian");
        }
        measure += point.weight * T->Weight();
      }
      result[{dimension, attribute}] += measure;
    }
  }
  return result;
}

int main(int argc, char **argv)
{
  if (argc != 3)
  {
    return 2;
  }
  mfem::Mpi::Init(argc, argv);
  mfem::Mesh mesh(argv[1], 1, 1, true);
  if (mesh.Dimension() != 3)
  {
    return 3;
  }
  std::vector<int> boundary_count(mesh.GetNFaces(), 0);
  for (int i = 0; i < mesh.GetNBE(); i++)
  {
    int face = mesh.GetBdrElementFaceIndex(i);
    if (++boundary_count.at(face) != 1)
    {
      throw std::runtime_error("Duplicate physical assignment of a mesh face");
    }
  }
  for (int face = 0; face < mesh.GetNFaces(); face++)
  {
    int first, second;
    mesh.GetFaceElements(face, &first, &second);
    if ((second < 0 || mesh.GetAttribute(first) != mesh.GetAttribute(second)) &&
        boundary_count[face] != 1)
    {
      throw std::runtime_error("Missing exterior or material-interface boundary face");
    }
  }
  const auto coarse = Integrate(mesh, 4), fine = Integrate(mesh, 8);
  nlohmann::json result = {{"Mesh", argv[1]}, {"QuadratureOrders", {4, 8}},
                           {"BoundaryCoverageChecked", true}};
  double difference = 0.0;
  for (const auto &[key, measure] : fine)
  {
    const char *group = key.first == 3 ? "MaterialVolumes" : "BoundaryAreas";
    result[group][std::to_string(key.second)] = measure;
    difference = std::max(difference, std::abs(measure - coarse.at(key)) /
                                          std::max(std::abs(measure), 1e-300));
  }
  result["MaximumRelativeQuadratureDifference"] = difference;
  std::ofstream output(argv[2]);
  if (!output)
  {
    return 4;
  }
  output << result.dump(2) << '\n';
  return output ? 0 : 4;
}
