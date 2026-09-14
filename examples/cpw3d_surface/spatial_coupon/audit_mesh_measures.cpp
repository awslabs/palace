// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

// Integrate physical measures of affine or curved meshes; do not substitute a
// center Jacobian when validating a curved geometry against CAD.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <numeric>
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
      auto geometry =
          dimension == 3 ? mesh.GetElementBaseGeometry(i) : mesh.GetBdrElementGeometry(i);
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
  std::vector<int> component(mesh.GetNE());
  std::iota(component.begin(), component.end(), 0);
  auto Find = [&](int element)
  {
    while (component[element] != element)
    {
      component[element] = component[component[element]];
      element = component[element];
    }
    return element;
  };
  long long interior_physical_faces = 0;
  for (int face = 0; face < mesh.GetNFaces(); face++)
  {
    int first, second;
    mesh.GetFaceElements(face, &first, &second);
    if ((second < 0 || mesh.GetAttribute(first) != mesh.GetAttribute(second)) &&
        boundary_count[face] != 1)
    {
      throw std::runtime_error("Missing exterior or material-interface boundary face");
    }
    if (second >= 0 && boundary_count[face])
    {
      interior_physical_faces++;
    }
    else if (second >= 0)
    {
      component[Find(first)] = Find(second);
    }
  }
  const auto coarse = Integrate(mesh, 4), fine = Integrate(mesh, 8);
  nlohmann::json result = {
      {"Mesh", argv[1]}, {"QuadratureOrders", {4, 8}}, {"BoundaryCoverageChecked", true}};
  bool all_tetrahedra = true;
  long long component_count = 0;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    all_tetrahedra &= mesh.GetElementBaseGeometry(element) == mfem::Geometry::TETRAHEDRON;
    component_count += Find(element) == element;
  }
  if (all_tetrahedra)
  {
    // Cutting the S physical/matching faces yields C connected volume components.
    // I of those faces are two-sided. Since 4T = S + I + 2F and F >= T - C,
    // every conforming tetrahedral replacement of this same shell needs
    // T >= ceil((S + I) / 2) - C. This is a necessary bound, not achievability.
    const long long shell_incidence = mesh.GetNBE() + interior_physical_faces;
    result["FrozenSurfaceElementLowerBound"] = {
        {"BoundaryTriangles", mesh.GetNBE()},
        {"TwoSidedInterfaceTriangles", interior_physical_faces},
        {"CutVolumeComponents", component_count},
        {"MinimumTetrahedra", (shell_incidence + 1) / 2 - component_count},
        {"Scope", "Necessary count bound with every physical/matching triangle preserved"}};
  }
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
