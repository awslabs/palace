// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "singularfield.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "utils/communication.hpp"

namespace palace
{

namespace fem
{

namespace singular
{

namespace
{

void ValidateInteriorPoint(const mfem::IntegrationPoint &point)
{
  const BarycentricPoint lambda{1.0 - point.x - point.y - point.z, point.x, point.y,
                                point.z};
  if (!std::all_of(lambda.begin(), lambda.end(),
                   [](double value) { return std::isfinite(value) && value > 0.0; }))
  {
    throw std::invalid_argument(
        "Singular H1 field evaluation requires a strictly interior tetrahedron point!");
  }
}

void ValidateClosurePoint(const mfem::IntegrationPoint &point)
{
  const BarycentricPoint lambda{1.0 - point.x - point.y - point.z, point.x, point.y,
                                point.z};
  double sum = 0.0;
  for (double value : lambda)
  {
    if (!std::isfinite(value) || value < -16.0 * std::numeric_limits<double>::epsilon())
    {
      throw std::invalid_argument(
          "Singular field closure evaluation requires a point in the closed reference "
          "tetrahedron!");
    }
    sum += value;
  }
  if (std::abs(sum - 1.0) > 16.0 * std::numeric_limits<double>::epsilon())
  {
    throw std::invalid_argument(
        "Singular field closure evaluation requires tetrahedron barycentric coordinates "
        "which sum to one!");
  }
}

void ValidateInteriorTrianglePoint(const mfem::IntegrationPoint &point)
{
  const TriangleBarycentricPoint lambda{1.0 - point.x - point.y, point.x, point.y};
  if (!std::all_of(lambda.begin(), lambda.end(),
                   [](double value) { return std::isfinite(value) && value > 0.0; }))
  {
    throw std::invalid_argument(
        "Singular H1 field evaluation requires a strictly interior triangle point!");
  }
}

void ValidateInteriorTriangleBarycentricPoint(const TriangleBarycentricPoint &lambda)
{
  double sum = 0.0;
  for (double value : lambda)
  {
    if (!std::isfinite(value) || !(value > 0.0))
    {
      throw std::invalid_argument(
          "Singular field evaluation requires strictly positive triangle barycentric "
          "coordinates!");
    }
    sum += value;
  }
  if (std::abs(sum - 1.0) > 16.0 * std::numeric_limits<double>::epsilon())
  {
    throw std::invalid_argument(
        "Singular field evaluation requires triangle barycentric coordinates which sum "
        "to one!");
  }
}

void ValidateTriangleClosurePoint(const mfem::IntegrationPoint &point)
{
  const TriangleBarycentricPoint lambda{1.0 - point.x - point.y, point.x, point.y};
  double sum = 0.0;
  for (double value : lambda)
  {
    if (!std::isfinite(value) || value < -16.0 * std::numeric_limits<double>::epsilon())
    {
      throw std::invalid_argument(
          "Singular field closure evaluation requires a point in the closed reference "
          "triangle!");
    }
    sum += value;
  }
  if (std::abs(sum - 1.0) > 16.0 * std::numeric_limits<double>::epsilon())
  {
    throw std::invalid_argument(
        "Singular field closure evaluation requires triangle barycentric coordinates "
        "which sum to one!");
  }
}

double Dot(const Vector3 &left, const Vector3 &right)
{
  return left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
}

double Dot(const Vector2 &left, const Vector2 &right)
{
  return left[0] * right[0] + left[1] * right[1];
}

int DecodeInteger(double value, const char *description)
{
  if (!std::isfinite(value) || value != std::round(value) ||
      value < static_cast<double>(std::numeric_limits<int>::min()) ||
      value > static_cast<double>(std::numeric_limits<int>::max()))
  {
    throw std::runtime_error(description);
  }
  return static_cast<int>(value);
}

HigherOrderBasisFamily DecodeBasisFamily(int value)
{
  if (value < static_cast<int>(HigherOrderBasisFamily::NODE_GRADIENT) ||
      value > static_cast<int>(HigherOrderBasisFamily::EDGE_ROTATIONAL))
  {
    throw std::runtime_error("Singular trace exchange received an invalid basis family!");
  }
  return static_cast<HigherOrderBasisFamily>(value);
}

template <typename ElementMap>
struct FaceNeighborRecordTraits;

template <>
struct FaceNeighborRecordTraits<ElementDofMap>
{
  static constexpr int RecordFields = 11;

  static const std::vector<ElementDof> &GetDofs(const ElementDofMap &element, bool h1)
  {
    return h1 ? element.h1 : element.nd;
  }

  static std::vector<ElementDof> &GetDofs(ElementDofMap &element, bool h1)
  {
    return h1 ? element.h1 : element.nd;
  }

  static void Pack(const ElementDof &element_dof, double *record)
  {
    int next = 0;
    record[next++] = static_cast<int>(element_dof.basis.family);
    for (int node : element_dof.basis.nodes)
    {
      record[next++] = node;
    }
    for (int index : element_dof.basis.interpolation_indices)
    {
      record[next++] = index;
    }
    record[next++] = element_dof.basis.order;
    record[next++] = element_dof.basis.nu;
    if (next != RecordFields)
    {
      throw std::logic_error(
          "Tetrahedral singular trace topology packing is inconsistent!");
    }
  }

  static ElementDof Unpack(const double *record, std::size_t coefficient)
  {
    int next = 0;
    HigherOrderBasis basis;
    basis.family = DecodeBasisFamily(DecodeInteger(
        record[next++], "Singular trace exchange received an invalid basis family!"));
    for (int &node : basis.nodes)
    {
      node = DecodeInteger(
          record[next++],
          "Tetrahedral singular trace exchange received an invalid basis node!");
    }
    for (int &index : basis.interpolation_indices)
    {
      index = DecodeInteger(
          record[next++],
          "Tetrahedral singular trace exchange received an invalid interpolation index!");
    }
    basis.order = DecodeInteger(
        record[next++],
        "Tetrahedral singular trace exchange received an invalid basis order!");
    basis.nu = record[next++];
    if (next != RecordFields || !std::isfinite(basis.nu) || !(basis.nu > 0.0) ||
        !(basis.nu < 1.0))
    {
      throw std::runtime_error(
          "Tetrahedral singular trace exchange received invalid topology!");
    }
    return {coefficient, basis};
  }
};

template <>
struct FaceNeighborRecordTraits<TriangleElementDofMap>
{
  static constexpr int RecordFields = 6;

  static const std::vector<TriangleElementDof> &
  GetDofs(const TriangleElementDofMap &element, bool h1)
  {
    return h1 ? element.h1 : element.nd;
  }

  static std::vector<TriangleElementDof> &GetDofs(TriangleElementDofMap &element, bool h1)
  {
    return h1 ? element.h1 : element.nd;
  }

  static void Pack(const TriangleElementDof &element_dof, double *record)
  {
    int next = 0;
    record[next++] = static_cast<int>(element_dof.basis.family);
    for (int node : element_dof.basis.nodes)
    {
      record[next++] = node;
    }
    record[next++] = element_dof.basis.order;
    record[next++] = element_dof.basis.nu;
    if (next != RecordFields)
    {
      throw std::logic_error("Triangular singular trace topology packing is inconsistent!");
    }
  }

  static TriangleElementDof Unpack(const double *record, std::size_t coefficient)
  {
    int next = 0;
    TriangleBasis basis;
    basis.family = DecodeBasisFamily(DecodeInteger(
        record[next++], "Singular trace exchange received an invalid basis family!"));
    for (int &node : basis.nodes)
    {
      node = DecodeInteger(
          record[next++],
          "Triangular singular trace exchange received an invalid basis node!");
    }
    basis.order = DecodeInteger(
        record[next++],
        "Triangular singular trace exchange received an invalid basis order!");
    basis.nu = record[next++];
    if (next != RecordFields || !std::isfinite(basis.nu) || !(basis.nu > 0.0) ||
        !(basis.nu < 1.0))
    {
      throw std::runtime_error(
          "Triangular singular trace exchange received invalid topology!");
    }
    return {coefficient, basis};
  }
};

template <typename ElementMap>
void InitializeFaceNeighborEnrichment(mfem::ParFiniteElementSpace &fespace,
                                      const std::vector<ElementMap> &elements, bool h1,
                                      FaceNeighborEnrichmentData<ElementMap> &data)
{
  auto *mesh = fespace.GetParMesh();
  if (!mesh || elements.size() != static_cast<std::size_t>(mesh->GetNE()))
  {
    throw std::invalid_argument(
        "Singular trace exchange requires matching parallel element topology!");
  }
  mesh->ExchangeFaceNbrData();

  using Traits = FaceNeighborRecordTraits<ElementMap>;
  std::size_t local_maximum = 0;
  for (const auto &element : elements)
  {
    local_maximum = std::max(local_maximum, Traits::GetDofs(element, h1).size());
  }
  bool valid = local_maximum <= static_cast<std::size_t>(std::numeric_limits<int>::max());
  Mpi::GlobalAnd(1, &valid, fespace.GetComm());
  if (!valid)
  {
    throw std::overflow_error(
        "Singular trace topology exceeds supported element record dimensions!");
  }
  int maximum = static_cast<int>(local_maximum);
  Mpi::GlobalMax(1, &maximum, fespace.GetComm());
  if (maximum > (std::numeric_limits<int>::max() - 1) / Traits::RecordFields)
  {
    throw std::overflow_error(
        "Singular trace topology exceeds supported element record dimensions!");
  }
  const int neighbor_elements = mesh->GetNFaceNeighborElements();
  const int metadata_size = 1 + maximum * Traits::RecordFields;
  valid =
      maximum == 0 || (neighbor_elements <= std::numeric_limits<int>::max() / maximum &&
                       mesh->GetNE() <= std::numeric_limits<int>::max() / maximum &&
                       mesh->GetNE() <= std::numeric_limits<int>::max() / metadata_size);
  Mpi::GlobalAnd(1, &valid, fespace.GetComm());
  if (!valid)
  {
    throw std::overflow_error(
        "Singular trace topology exceeds supported finite element space dimensions!");
  }
  data.maximum_element_dofs = maximum;
  data.element_dofs.resize(neighbor_elements);
  data.coefficients.SetSize(neighbor_elements * maximum);
  data.coefficients = 0.0;
  if (maximum == 0)
  {
    return;
  }

  data.coefficient_collection =
      std::make_unique<mfem::L2_FECollection>(0, mesh->Dimension());
  data.coefficient_space = std::make_unique<mfem::ParFiniteElementSpace>(
      mesh, data.coefficient_collection.get(), maximum, mfem::Ordering::byVDIM);
  data.coefficient_field =
      std::make_unique<mfem::ParGridFunction>(data.coefficient_space.get());

  mfem::L2_FECollection metadata_collection(0, mesh->Dimension());
  mfem::ParFiniteElementSpace metadata_space(mesh, &metadata_collection, metadata_size,
                                             mfem::Ordering::byVDIM);
  mfem::ParGridFunction metadata(&metadata_space);
  metadata = 0.0;
  mfem::Array<int> vdofs;
  mfem::Vector values(metadata_size);
  for (int element = 0; element < mesh->GetNE(); element++)
  {
    values = 0.0;
    const auto &element_dofs = Traits::GetDofs(elements[element], h1);
    values[0] = static_cast<int>(element_dofs.size());
    for (std::size_t record = 0; record < element_dofs.size(); record++)
    {
      Traits::Pack(element_dofs[record],
                   values.GetData() + 1 + record * Traits::RecordFields);
    }
    metadata_space.GetElementVDofs(element, vdofs);
    if (vdofs.Size() != metadata_size)
    {
      throw std::runtime_error(
          "Singular trace metadata space has inconsistent element dimensions!");
    }
    // SetSubVector and GetElementDofValues apply MFEM's signed-VDof convention.
    metadata.SetSubVector(vdofs, values);
  }
  metadata.ExchangeFaceNbrData();

  for (int neighbor = 0; neighbor < mesh->GetNFaceNeighborElements(); neighbor++)
  {
    metadata.GetElementDofValues(mesh->GetNE() + neighbor, values);
    if (values.Size() != metadata_size)
    {
      throw std::runtime_error(
          "Singular trace exchange received inconsistent metadata dimensions!");
    }
    const int count = DecodeInteger(
        values[0], "Singular trace exchange received an invalid element record count!");
    if (count < 0 || count > maximum)
    {
      throw std::runtime_error(
          "Singular trace exchange received an invalid element record count!");
    }
    auto &element_dofs = Traits::GetDofs(data.element_dofs[neighbor], h1);
    element_dofs.reserve(count);
    for (int record = 0; record < count; record++)
    {
      const std::size_t coefficient = static_cast<std::size_t>(neighbor) * maximum + record;
      element_dofs.push_back(Traits::Unpack(
          values.GetData() + 1 + record * Traits::RecordFields, coefficient));
    }
  }
}

template <typename ElementMap>
void ExchangeFaceNeighborEnrichment(const std::vector<ElementMap> &elements,
                                    const mfem::Vector &local_coefficients, bool h1,
                                    mfem::ParFiniteElementSpace &fespace,
                                    FaceNeighborEnrichmentData<ElementMap> &data)
{
  using Traits = FaceNeighborRecordTraits<ElementMap>;
  const int maximum = data.maximum_element_dofs;
  if (maximum == 0)
  {
    return;
  }
  if (!data.coefficient_space || !data.coefficient_field ||
      data.coefficient_space->GetParMesh() != fespace.GetParMesh() ||
      elements.size() != static_cast<std::size_t>(fespace.GetNE()))
  {
    throw std::logic_error(
        "Singular trace coefficient exchange is not initialized consistently!");
  }

  auto &field = *data.coefficient_field;
  field = 0.0;
  mfem::Array<int> vdofs;
  mfem::Vector values(maximum);
  for (int element = 0; element < fespace.GetNE(); element++)
  {
    values = 0.0;
    const auto &element_dofs = Traits::GetDofs(elements[element], h1);
    if (element_dofs.size() > static_cast<std::size_t>(maximum))
    {
      throw std::logic_error(
          "Singular trace coefficient exchange exceeds its element record dimensions!");
    }
    for (std::size_t record = 0; record < element_dofs.size(); record++)
    {
      if (element_dofs[record].dof >= static_cast<std::size_t>(local_coefficients.Size()))
      {
        throw std::out_of_range(
            "Singular trace coefficient exchange references an invalid local DOF!");
      }
      values[static_cast<int>(record)] =
          local_coefficients[static_cast<int>(element_dofs[record].dof)];
    }
    data.coefficient_space->GetElementVDofs(element, vdofs);
    if (vdofs.Size() != maximum)
    {
      throw std::runtime_error(
          "Singular trace coefficient space has inconsistent element dimensions!");
    }
    // This preserves any orientation sign encoded as MFEM's -1-vdof convention.
    field.SetSubVector(vdofs, values);
  }
  field.ExchangeFaceNbrData();

  data.coefficients = 0.0;
  for (int neighbor = 0; neighbor < fespace.GetParMesh()->GetNFaceNeighborElements();
       neighbor++)
  {
    field.GetElementDofValues(fespace.GetNE() + neighbor, values);
    if (values.Size() != maximum)
    {
      throw std::runtime_error(
          "Singular trace exchange received inconsistent coefficient dimensions!");
    }
    const auto &element_dofs = Traits::GetDofs(data.element_dofs[neighbor], h1);
    for (std::size_t record = 0; record < element_dofs.size(); record++)
    {
      const std::size_t coefficient = static_cast<std::size_t>(neighbor) * maximum + record;
      data.coefficients[static_cast<int>(coefficient)] = values[static_cast<int>(record)];
    }
  }
}

template <typename ElementMap>
const ElementMap &
GetParallelElementDofs(int element, int local_elements,
                       const std::vector<ElementMap> &local_maps,
                       const FaceNeighborEnrichmentData<ElementMap> &face_neighbor)
{
  if (element < 0)
  {
    throw std::out_of_range("Singular trace element index is negative!");
  }
  if (element < local_elements)
  {
    return local_maps[element];
  }
  const int neighbor = element - local_elements;
  if (neighbor < 0 || neighbor >= static_cast<int>(face_neighbor.element_dofs.size()))
  {
    throw std::out_of_range("Singular trace face-neighbor element is out of range!");
  }
  return face_neighbor.element_dofs[neighbor];
}

mfem::ElementTransformation &
GetParallelElementTransformation(mfem::ParFiniteElementSpace &fespace, int element)
{
  const int local_elements = fespace.GetNE();
  mfem::ElementTransformation *transformation =
      element < local_elements
          ? fespace.GetElementTransformation(element)
          : fespace.GetFaceNbrElementTransformation(element - local_elements);
  if (!transformation)
  {
    throw std::runtime_error("Singular trace element has no physical transformation!");
  }
  return *transformation;
}

void GetStandardGradient(mfem::ParGridFunction &field, mfem::ParFiniteElementSpace &fespace,
                         int element, mfem::ElementTransformation &transformation,
                         const mfem::IntegrationPoint &point, mfem::Vector &gradient)
{
  if (element < fespace.GetNE())
  {
    field.GetGradient(transformation, gradient);
    return;
  }
  const int neighbor = element - fespace.GetNE();
  mfem::Array<int> dofs;
  mfem::DofTransformation dof_transformation;
  fespace.GetFaceNbrElementVDofs(neighbor, dofs, dof_transformation);
  mfem::Vector local_values;
  field.FaceNbrData().GetSubVector(dofs, local_values);
  dof_transformation.InvTransformPrimal(local_values);
  const auto *finite_element = fespace.GetFaceNbrFE(neighbor);
  if (!finite_element || finite_element->GetRangeType() != mfem::FiniteElement::SCALAR ||
      finite_element->GetMapType() != mfem::FiniteElement::VALUE ||
      finite_element->GetDof() != local_values.Size())
  {
    throw std::runtime_error(
        "Singular trace received an incompatible face-neighbor H1 element!");
  }
  transformation.SetIntPoint(&point);
  mfem::DenseMatrix derivative(finite_element->GetDof(), transformation.GetSpaceDim());
  finite_element->CalcPhysDShape(transformation, derivative);
  gradient.SetSize(transformation.GetSpaceDim());
  derivative.MultTranspose(local_values, gradient);
}

void GetStandardCurl(mfem::ParGridFunction &field, mfem::ParFiniteElementSpace &fespace,
                     int element, mfem::ElementTransformation &transformation,
                     const mfem::IntegrationPoint &point, mfem::Vector &curl)
{
  if (element < fespace.GetNE())
  {
    field.GetCurl(transformation, curl);
    return;
  }
  const int neighbor = element - fespace.GetNE();
  mfem::Array<int> dofs;
  mfem::DofTransformation dof_transformation;
  fespace.GetFaceNbrElementVDofs(neighbor, dofs, dof_transformation);
  mfem::Vector local_values;
  field.FaceNbrData().GetSubVector(dofs, local_values);
  dof_transformation.InvTransformPrimal(local_values);
  const auto *finite_element = fespace.GetFaceNbrFE(neighbor);
  if (!finite_element || finite_element->GetRangeType() != mfem::FiniteElement::VECTOR ||
      finite_element->GetDof() != local_values.Size())
  {
    throw std::runtime_error(
        "Singular trace received an incompatible face-neighbor H(curl) element!");
  }
  transformation.SetIntPoint(&point);
  mfem::DenseMatrix curl_shape(finite_element->GetDof(), finite_element->GetCurlDim());
  finite_element->CalcPhysCurlShape(transformation, curl_shape);
  curl.SetSize(curl_shape.Width());
  curl_shape.MultTranspose(local_values, curl);
}

double DistanceToLine(const mfem::Vector &point, const mfem::Vector &first,
                      const mfem::Vector &second)
{
  if (point.Size() != 3 || first.Size() != 3 || second.Size() != 3)
  {
    throw std::invalid_argument(
        "Singular edge diagnostic requires three-dimensional points!");
  }
  Vector3 direction, offset;
  for (int d = 0; d < 3; d++)
  {
    direction[d] = second[d] - first[d];
    offset[d] = point[d] - first[d];
  }
  const Vector3 cross{offset[1] * direction[2] - offset[2] * direction[1],
                      offset[2] * direction[0] - offset[0] * direction[2],
                      offset[0] * direction[1] - offset[1] * direction[0]};
  const double direction_norm = std::sqrt(Dot(direction, direction));
  if (!std::isfinite(direction_norm) || !(direction_norm > 0.0))
  {
    throw std::runtime_error("Singular edge diagnostic found a degenerate mesh edge!");
  }
  return std::sqrt(Dot(cross, cross)) / direction_norm;
}

class PotentialCoefficient final : public mfem::Coefficient
{
private:
  EnrichedH1FieldEvaluator &evaluator;

public:
  explicit PotentialCoefficient(EnrichedH1FieldEvaluator &evaluator) : evaluator(evaluator)
  {
  }

  double Eval(mfem::ElementTransformation &transformation,
              const mfem::IntegrationPoint &point) override
  {
    return evaluator.Evaluate(transformation.ElementNo, point).potential;
  }
};

class ElectricFieldCoefficient final : public mfem::VectorCoefficient
{
private:
  EnrichedH1FieldEvaluator &evaluator;

public:
  explicit ElectricFieldCoefficient(EnrichedH1FieldEvaluator &evaluator)
    : mfem::VectorCoefficient(3), evaluator(evaluator)
  {
  }

  void Eval(mfem::Vector &value, mfem::ElementTransformation &transformation,
            const mfem::IntegrationPoint &point) override
  {
    const auto field = evaluator.Evaluate(transformation.ElementNo, point);
    value.SetSize(3);
    for (int d = 0; d < 3; d++)
    {
      value[d] = -field.gradient[d];
    }
  }
};

class TrianglePotentialCoefficient final : public mfem::Coefficient
{
private:
  TriangleEnrichedH1FieldEvaluator &evaluator;

public:
  explicit TrianglePotentialCoefficient(TriangleEnrichedH1FieldEvaluator &evaluator)
    : evaluator(evaluator)
  {
  }

  double Eval(mfem::ElementTransformation &transformation,
              const mfem::IntegrationPoint &point) override
  {
    return evaluator.Evaluate(transformation.ElementNo, point).potential;
  }
};

class TriangleElectricFieldCoefficient final : public mfem::VectorCoefficient
{
private:
  TriangleEnrichedH1FieldEvaluator &evaluator;

public:
  explicit TriangleElectricFieldCoefficient(TriangleEnrichedH1FieldEvaluator &evaluator)
    : mfem::VectorCoefficient(2), evaluator(evaluator)
  {
  }

  void Eval(mfem::Vector &value, mfem::ElementTransformation &transformation,
            const mfem::IntegrationPoint &point) override
  {
    const auto field = evaluator.Evaluate(transformation.ElementNo, point);
    value.SetSize(2);
    for (int d = 0; d < 2; d++)
    {
      value[d] = -field.gradient[d];
    }
  }
};

class TriangleNDFieldCoefficient final : public mfem::VectorCoefficient
{
private:
  TriangleEnrichedNDFieldEvaluator &evaluator;

public:
  explicit TriangleNDFieldCoefficient(TriangleEnrichedNDFieldEvaluator &evaluator)
    : mfem::VectorCoefficient(2), evaluator(evaluator)
  {
  }

  void Eval(mfem::Vector &value, mfem::ElementTransformation &transformation,
            const mfem::IntegrationPoint &point) override
  {
    const auto field = evaluator.Evaluate(transformation.ElementNo, point);
    value.SetSize(2);
    for (int d = 0; d < 2; d++)
    {
      value[d] = field.value[d];
    }
  }
};

class TriangleNDCurlCoefficient final : public mfem::Coefficient
{
private:
  TriangleEnrichedNDFieldEvaluator &evaluator;

public:
  explicit TriangleNDCurlCoefficient(TriangleEnrichedNDFieldEvaluator &evaluator)
    : evaluator(evaluator)
  {
  }

  double Eval(mfem::ElementTransformation &transformation,
              const mfem::IntegrationPoint &point) override
  {
    return evaluator.Evaluate(transformation.ElementNo, point).curl;
  }
};

template <typename FieldNorm>
std::vector<H1TipSlopeDiagnostic>
FitTriangleTipSlopes(mfem::ParFiniteElementSpace &fespace,
                     const TriangleFeatureTopology &features,
                     const std::vector<GlobalVertexId> &source_vertex_ids,
                     const std::vector<GlobalVertexId> &source_element_ids,
                     const EdgeSlopeOptions &options, FieldNorm &&field_norm_at_point)
{
  if (features.elements.size() != static_cast<std::size_t>(fespace.GetNE()) ||
      source_vertex_ids.size() != static_cast<std::size_t>(fespace.GetParMesh()->GetNV()) ||
      source_element_ids.size() != static_cast<std::size_t>(fespace.GetNE()) ||
      options.sample_count < 3 || !std::isfinite(options.minimum_barycentric_radius) ||
      !std::isfinite(options.maximum_barycentric_radius) ||
      !(options.minimum_barycentric_radius > 0.0) ||
      !(options.minimum_barycentric_radius < options.maximum_barycentric_radius) ||
      !(options.maximum_barycentric_radius < 1.0))
  {
    throw std::invalid_argument(
        "Triangular singular tip slope diagnostics received inconsistent topology or "
        "options!");
  }

  const double log_minimum_radius = std::log(options.minimum_barycentric_radius);
  const double log_maximum_radius = std::log(options.maximum_barycentric_radius);
  std::vector<H1TipSlopeDiagnostic> diagnostics;
  mfem::IntegrationPoint point;
  mfem::Vector physical_point(2);
  for (int element = 0; element < fespace.GetNE(); element++)
  {
    const auto *triangle = fespace.GetParMesh()->GetElement(element);
    auto *transformation = fespace.GetElementTransformation(element);
    if (!triangle || triangle->GetGeometryType() != mfem::Geometry::TRIANGLE ||
        !transformation)
    {
      throw std::runtime_error(
          "Triangular singular tip slope diagnostics require triangular elements!");
    }
    for (const auto &node : features.elements[element].nodes)
    {
      if (node.vertex >= features.vertices.size())
      {
        throw std::invalid_argument(
            "Triangular singular tip slope diagnostics found invalid feature incidence!");
      }
      const auto &feature = features.vertices[node.vertex];
      if (feature.id != node.vertex || feature.selected_segments.empty() ||
          feature.selected_segments.size() > 2 ||
          std::any_of(feature.selected_segments.begin(), feature.selected_segments.end(),
                      [&features](std::size_t segment)
                      { return segment >= features.selected_segments.size(); }) ||
          !std::isfinite(feature.nu) || !(feature.nu > 0.0) || !(feature.nu < 1.0))
      {
        throw std::invalid_argument(
            "Triangular singular tip slope diagnostics found an invalid feature!");
      }

      std::array<bool, 3> seen{false, false, false};
      std::array<GlobalVertexId, 3> canonical_vertices;
      for (int local = 0; local < 3; local++)
      {
        const int barycentric_node = node.canonical_nodes[local];
        if (barycentric_node < 0 || barycentric_node >= 3 || seen[barycentric_node])
        {
          throw std::invalid_argument(
              "Triangular singular tip slope diagnostics found invalid canonical "
              "nodes!");
        }
        seen[barycentric_node] = true;
        const int mesh_vertex = triangle->GetVertices()[barycentric_node];
        if (mesh_vertex < 0 || mesh_vertex >= static_cast<int>(source_vertex_ids.size()))
        {
          throw std::invalid_argument(
              "Triangular singular tip slope diagnostics found invalid mesh vertices!");
        }
        canonical_vertices[local] = source_vertex_ids[mesh_vertex];
      }
      if (canonical_vertices[0] != feature.mesh_vertex)
      {
        throw std::invalid_argument(
            "Triangular singular tip incidence does not match its source vertex!");
      }

      mfem::Vector tip(2);
      point.Set2(node.canonical_nodes[0] == 1 ? 1.0 : 0.0,
                 node.canonical_nodes[0] == 2 ? 1.0 : 0.0);
      transformation->Transform(point, tip);
      std::vector<double> log_distance(options.sample_count);
      std::vector<double> log_field_norm(options.sample_count);
      std::vector<double> distance(options.sample_count);
      std::vector<double> field_norm(options.sample_count);
      bool valid = true;
      for (int sample = 0; sample < options.sample_count; sample++)
      {
        const double fraction =
            static_cast<double>(sample) / static_cast<double>(options.sample_count - 1);
        const double radius = std::exp(
            log_minimum_radius + fraction * (log_maximum_radius - log_minimum_radius));
        TriangleBarycentricPoint lambda;
        lambda[node.canonical_nodes[0]] = 1.0 - radius;
        lambda[node.canonical_nodes[1]] = 0.5 * radius;
        lambda[node.canonical_nodes[2]] = 0.5 * radius;
        point.Set2(lambda[1], lambda[2]);
        transformation->Transform(point, physical_point);
        const double dx = physical_point[0] - tip[0];
        const double dy = physical_point[1] - tip[1];
        distance[sample] = std::sqrt(dx * dx + dy * dy);
        field_norm[sample] = field_norm_at_point(element, point);
        if (!std::isfinite(distance[sample]) || !(distance[sample] > 0.0) ||
            !std::isfinite(field_norm[sample]) || !(field_norm[sample] > 0.0))
        {
          valid = false;
          break;
        }
        log_distance[sample] = std::log(distance[sample]);
        log_field_norm[sample] = std::log(field_norm[sample]);
      }

      double fitted_slope = 0.0;
      double r_squared = 0.0;
      if (valid)
      {
        long double mean_distance = 0.0L;
        long double mean_field = 0.0L;
        for (int sample = 0; sample < options.sample_count; sample++)
        {
          mean_distance += log_distance[sample];
          mean_field += log_field_norm[sample];
        }
        mean_distance /= options.sample_count;
        mean_field /= options.sample_count;
        long double distance_variance = 0.0L;
        long double covariance = 0.0L;
        for (int sample = 0; sample < options.sample_count; sample++)
        {
          const long double centered_distance = log_distance[sample] - mean_distance;
          distance_variance += centered_distance * centered_distance;
          covariance += centered_distance * (log_field_norm[sample] - mean_field);
        }
        if (!(distance_variance > 0.0L))
        {
          valid = false;
        }
        else
        {
          fitted_slope = static_cast<double>(covariance / distance_variance);
          const long double intercept = mean_field - fitted_slope * mean_distance;
          long double residual_sum = 0.0L;
          long double total_sum = 0.0L;
          for (int sample = 0; sample < options.sample_count; sample++)
          {
            const long double residual =
                log_field_norm[sample] - (intercept + fitted_slope * log_distance[sample]);
            const long double centered = log_field_norm[sample] - mean_field;
            residual_sum += residual * residual;
            total_sum += centered * centered;
          }
          r_squared = (total_sum > 0.0L)
                          ? static_cast<double>(1.0L - residual_sum / total_sum)
                          : ((residual_sum == 0.0L) ? 1.0 : 0.0);
          valid = std::isfinite(fitted_slope) && std::isfinite(r_squared);
        }
      }

      diagnostics.push_back(
          {source_element_ids[element], node.vertex, feature.selected_segments[0],
           canonical_vertices, options.sample_count, feature.nu, feature.nu - 1.0,
           fitted_slope, r_squared, valid ? distance.front() : 0.0,
           valid ? distance.back() : 0.0, valid ? field_norm.front() : 0.0,
           valid ? field_norm.back() : 0.0, valid});
    }
  }
  return diagnostics;
}

}  // namespace

H1FieldValue EvaluateElementH1Enrichment(const ElementDofMap &element_dofs,
                                         const mfem::Vector &local_coefficients,
                                         const BarycentricPoint &lambda,
                                         const BarycentricGradients &grad_lambda)
{
  H1FieldValue result;
  for (const auto &element_dof : element_dofs.h1)
  {
    if (element_dof.dof >= static_cast<std::size_t>(local_coefficients.Size()) ||
        !IsGradientFamily(element_dof.basis.family))
    {
      throw std::invalid_argument(
          "Singular H1 field evaluation received an invalid element DOF map!");
    }
    const double coefficient = local_coefficients[element_dof.dof];
    if (!std::isfinite(coefficient))
    {
      throw std::invalid_argument(
          "Singular H1 field evaluation received a nonfinite coefficient!");
    }
    const auto basis = EvaluateHigherOrderBasis(lambda, grad_lambda, element_dof.basis);
    result.potential +=
        coefficient * EvaluateHigherOrderGradientPotential(lambda, element_dof.basis);
    for (int d = 0; d < 3; d++)
    {
      result.gradient[d] += coefficient * basis.value[d];
    }
  }
  return result;
}

NDFieldValue EvaluateElementNDEnrichment(const ElementDofMap &element_dofs,
                                         const mfem::Vector &local_coefficients,
                                         const BarycentricPoint &lambda,
                                         const BarycentricGradients &grad_lambda)
{
  NDFieldValue result;
  for (const auto &element_dof : element_dofs.nd)
  {
    if (element_dof.dof >= static_cast<std::size_t>(local_coefficients.Size()))
    {
      throw std::invalid_argument(
          "Singular ND field evaluation received an invalid element DOF map!");
    }
    const double coefficient = local_coefficients[element_dof.dof];
    if (!std::isfinite(coefficient))
    {
      throw std::invalid_argument(
          "Singular ND field evaluation received a nonfinite coefficient!");
    }
    const auto basis = EvaluateHigherOrderBasis(lambda, grad_lambda, element_dof.basis);
    for (int d = 0; d < 3; d++)
    {
      result.value[d] += coefficient * basis.value[d];
      result.curl[d] += coefficient * basis.curl[d];
    }
  }
  return result;
}

Vector3 EvaluateElementNDEnrichmentValue(const ElementDofMap &element_dofs,
                                         const mfem::Vector &local_coefficients,
                                         const BarycentricPoint &lambda,
                                         const BarycentricGradients &grad_lambda)
{
  Vector3 result{};
  for (const auto &element_dof : element_dofs.nd)
  {
    if (element_dof.dof >= static_cast<std::size_t>(local_coefficients.Size()))
    {
      throw std::invalid_argument(
          "Singular ND field evaluation received an invalid element DOF map!");
    }
    const double coefficient = local_coefficients[element_dof.dof];
    if (!std::isfinite(coefficient))
    {
      throw std::invalid_argument(
          "Singular ND field evaluation received a nonfinite coefficient!");
    }
    const auto value =
        EvaluateHigherOrderBasisValue(lambda, grad_lambda, element_dof.basis);
    for (int d = 0; d < 3; d++)
    {
      result[d] += coefficient * value[d];
    }
  }
  return result;
}

TriangleH1FieldValue EvaluateElementTriangleH1Enrichment(
    const TriangleElementDofMap &element_dofs, const mfem::Vector &local_coefficients,
    const TriangleBarycentricPoint &lambda, const TriangleBarycentricGradients &grad_lambda)
{
  TriangleH1FieldValue result;
  for (const auto &element_dof : element_dofs.h1)
  {
    if (element_dof.dof >= static_cast<std::size_t>(local_coefficients.Size()) ||
        element_dof.basis.family != HigherOrderBasisFamily::NODE_GRADIENT ||
        element_dof.basis.order != 1)
    {
      throw std::invalid_argument(
          "Triangular singular H1 field evaluation received an invalid element DOF map!");
    }
    const double coefficient = local_coefficients[element_dof.dof];
    if (!std::isfinite(coefficient))
    {
      throw std::invalid_argument(
          "Triangular singular H1 field evaluation received a nonfinite coefficient!");
    }
    const auto &basis = element_dof.basis;
    const auto gradient = EvaluateTriangleNodeGradient(lambda, grad_lambda, basis.nodes[0],
                                                       basis.nodes[1], basis.nu);
    result.potential += coefficient * EvaluateTriangleNodeGradientPotential(
                                          lambda, basis.nodes[0], basis.nodes[1], basis.nu);
    for (int d = 0; d < 2; d++)
    {
      result.gradient[d] += coefficient * gradient.value[d];
    }
  }
  return result;
}

TriangleNDFieldValue EvaluateElementTriangleNDEnrichment(
    const TriangleElementDofMap &element_dofs, const mfem::Vector &local_coefficients,
    const TriangleBarycentricPoint &lambda, const TriangleBarycentricGradients &grad_lambda)
{
  TriangleNDFieldValue result;
  for (const auto &element_dof : element_dofs.nd)
  {
    if (element_dof.dof >= static_cast<std::size_t>(local_coefficients.Size()) ||
        element_dof.basis.order != 1 ||
        (element_dof.basis.family != HigherOrderBasisFamily::NODE_GRADIENT &&
         element_dof.basis.family != HigherOrderBasisFamily::NODE_ROTATIONAL))
    {
      throw std::invalid_argument(
          "Triangular singular ND field evaluation received an invalid element DOF "
          "map!");
    }
    const double coefficient = local_coefficients[element_dof.dof];
    if (!std::isfinite(coefficient))
    {
      throw std::invalid_argument(
          "Triangular singular ND field evaluation received a nonfinite coefficient!");
    }
    const auto &basis = element_dof.basis;
    const auto value =
        (basis.family == HigherOrderBasisFamily::NODE_GRADIENT)
            ? EvaluateTriangleNodeGradient(lambda, grad_lambda, basis.nodes[0],
                                           basis.nodes[1], basis.nu)
            : EvaluateTriangleNodeRotational(lambda, grad_lambda, basis.nodes[0],
                                             basis.nodes[1], basis.nodes[2], basis.nu);
    for (int d = 0; d < 2; d++)
    {
      result.value[d] += coefficient * value.value[d];
    }
    result.curl += coefficient * value.curl;
  }
  return result;
}

TriangleBoundaryModeMagneticFieldValue ReconstructTriangleBoundaryModeMagneticField(
    const Vector2 &et_real, const Vector2 &et_imag, double curl_et_real,
    double curl_et_imag, const Vector2 &grad_en_real, const Vector2 &grad_en_imag,
    std::complex<double> kn, double omega)
{
  const bool valid = std::all_of(et_real.begin(), et_real.end(),
                                 [](double value) { return std::isfinite(value); }) &&
                     std::all_of(et_imag.begin(), et_imag.end(),
                                 [](double value) { return std::isfinite(value); }) &&
                     std::all_of(grad_en_real.begin(), grad_en_real.end(),
                                 [](double value) { return std::isfinite(value); }) &&
                     std::all_of(grad_en_imag.begin(), grad_en_imag.end(),
                                 [](double value) { return std::isfinite(value); }) &&
                     std::isfinite(curl_et_real) && std::isfinite(curl_et_imag) &&
                     std::isfinite(kn.real()) && std::isfinite(kn.imag()) &&
                     std::isfinite(omega) && omega > 0.0;
  if (!valid)
  {
    throw std::invalid_argument(
        "Triangular BoundaryMode magnetic reconstruction received nonfinite fields or "
        "an invalid frequency!");
  }

  TriangleBoundaryModeMagneticFieldValue result;
  result.normal_real = curl_et_imag / omega;
  result.normal_imag = -curl_et_real / omega;
  const std::complex<double> q = -kn / omega;
  const Vector2 z_cross_et_real{-et_real[1], et_real[0]};
  const Vector2 z_cross_et_imag{-et_imag[1], et_imag[0]};
  const Vector2 grad_en_cross_z_real{grad_en_real[1], -grad_en_real[0]};
  const Vector2 grad_en_cross_z_imag{grad_en_imag[1], -grad_en_imag[0]};
  for (int d = 0; d < 2; d++)
  {
    result.transverse_real[d] = q.real() * z_cross_et_real[d] -
                                q.imag() * z_cross_et_imag[d] +
                                grad_en_cross_z_imag[d] / omega;
    result.transverse_imag[d] = q.real() * z_cross_et_imag[d] +
                                q.imag() * z_cross_et_real[d] -
                                grad_en_cross_z_real[d] / omega;
  }
  return result;
}

EnrichedH1FieldEvaluator::EnrichedH1FieldEvaluator(const DofTopology &topology,
                                                   const ParallelDofNumbering &numbering,
                                                   mfem::ParFiniteElementSpace &fespace)
  : topology(topology), numbering(numbering), fespace(fespace), standard_field(&fespace),
    enrichment_prolongation(
        BuildParallelEnrichmentProlongation(fespace.GetComm(), numbering.h1)),
    h1_exponents(topology.h1_dofs.size(), std::numeric_limits<double>::quiet_NaN()),
    initialized(false)
{
  bool valid =
      fespace.GetParMesh() &&
      topology.elements.size() == static_cast<std::size_t>(fespace.GetNE()) &&
      topology.h1_dofs.size() == numbering.h1.local_to_true.size() &&
      numbering.h1.local_size == static_cast<HYPRE_BigInt>(topology.h1_dofs.size()) &&
      numbering.h1.owned_size <=
          static_cast<HYPRE_BigInt>(std::numeric_limits<int>::max()) &&
      enrichment_prolongation &&
      enrichment_prolongation->Height() == static_cast<int>(topology.h1_dofs.size()) &&
      enrichment_prolongation->Width() == static_cast<int>(numbering.h1.owned_size);
  if (valid)
  {
    for (const auto &element : topology.elements)
    {
      for (const auto &element_dof : element.h1)
      {
        if (element_dof.dof >= h1_exponents.size() ||
            !IsGradientFamily(element_dof.basis.family) ||
            topology.h1_dofs[element_dof.dof].family != element_dof.basis.family ||
            topology.h1_dofs[element_dof.dof].order != element_dof.basis.order ||
            !std::isfinite(element_dof.basis.nu) || !(element_dof.basis.nu > 0.0) ||
            !(element_dof.basis.nu < 1.0))
        {
          valid = false;
          break;
        }
        double &exponent = h1_exponents[element_dof.dof];
        if (std::isnan(exponent))
        {
          exponent = element_dof.basis.nu;
        }
        else if (exponent != element_dof.basis.nu)
        {
          valid = false;
          break;
        }
      }
      if (!valid)
      {
        break;
      }
    }
    valid = valid && std::all_of(h1_exponents.begin(), h1_exponents.end(),
                                 [](double exponent) { return std::isfinite(exponent); });
  }
  Mpi::GlobalAnd(1, &valid, fespace.GetComm());
  if (!valid)
  {
    throw std::invalid_argument(
        "Singular H1 field evaluator received inconsistent space and DOF topology!");
  }
  local_enrichment.SetSize(static_cast<int>(topology.h1_dofs.size()));
  InitializeFaceNeighborEnrichment(fespace, topology.elements, true,
                                   face_neighbor_enrichment);
}

void EnrichedH1FieldEvaluator::SetFromTrueDofs(const mfem::Vector &combined_true_dofs)
{
  const int standard_size = fespace.GetTrueVSize();
  const int enrichment_size = static_cast<int>(numbering.h1.owned_size);
  bool valid = combined_true_dofs.Size() == standard_size + enrichment_size;
  Mpi::GlobalAnd(1, &valid, fespace.GetComm());
  if (!valid)
  {
    throw std::invalid_argument(
        "Combined singular H1 vector has inconsistent process-local dimensions!");
  }

  mfem::Vector standard_true_dofs(standard_size);
  mfem::Vector enrichment_true_dofs(enrichment_size);
  for (int i = 0; i < standard_size; i++)
  {
    standard_true_dofs[i] = combined_true_dofs[i];
  }
  for (int i = 0; i < enrichment_size; i++)
  {
    enrichment_true_dofs[i] = combined_true_dofs[standard_size + i];
  }
  standard_field.SetFromTrueDofs(standard_true_dofs);
  standard_field.ExchangeFaceNbrData();
  enrichment_prolongation->Mult(enrichment_true_dofs, local_enrichment);
  ExchangeFaceNeighborEnrichment(topology.elements, local_enrichment, true, fespace,
                                 face_neighbor_enrichment);
  initialized = true;
}

H1FieldValue EnrichedH1FieldEvaluator::Evaluate(int element,
                                                const mfem::IntegrationPoint &point)
{
  ValidateInteriorPoint(point);
  const BarycentricPoint lambda{1.0 - point.x - point.y - point.z, point.x, point.y,
                                point.z};
  return EvaluateBarycentric(element, point, lambda);
}

H1FieldValue EnrichedH1FieldEvaluator::EvaluateClosure(int element,
                                                       const mfem::IntegrationPoint &point)
{
  ValidateClosurePoint(point);
  const BarycentricPoint lambda{1.0 - point.x - point.y - point.z, point.x, point.y,
                                point.z};
  return EvaluateBarycentric(element, point, lambda);
}

std::vector<TetrahedronFaceSingularity>
EnrichedH1FieldEvaluator::GetElementFaceSingularities(
    int element, const std::array<int, 3> &face_nodes) const
{
  const auto &element_dofs = GetParallelElementDofs(
      element, fespace.GetNE(), topology.elements, face_neighbor_enrichment);
  auto sorted_face = face_nodes;
  std::sort(sorted_face.begin(), sorted_face.end());
  if (sorted_face[0] < 0 || sorted_face[2] >= 4 ||
      std::adjacent_find(sorted_face.begin(), sorted_face.end()) != sorted_face.end())
  {
    throw std::invalid_argument(
        "Singular H1 face-singularity query received invalid face nodes!");
  }
  const auto on_face = [&sorted_face](int node)
  { return std::binary_search(sorted_face.begin(), sorted_face.end(), node); };

  std::vector<TetrahedronFaceSingularity> result;
  const auto add =
      [&result](TetrahedronFaceSingularityType type, std::array<int, 2> nodes, double nu)
  {
    if (type == TetrahedronFaceSingularityType::EDGE && nodes[1] < nodes[0])
    {
      std::swap(nodes[0], nodes[1]);
    }
    const auto existing =
        std::find_if(result.begin(), result.end(), [type, &nodes](const auto &feature)
                     { return feature.type == type && feature.nodes == nodes; });
    if (existing == result.end())
    {
      result.push_back({type, nodes, nu});
    }
    else if (existing->nu != nu)
    {
      throw std::runtime_error(
          "Singular H1 face contains inconsistent exponents for one feature!");
    }
  };

  for (const auto &element_dof : element_dofs.h1)
  {
    const auto &basis = element_dof.basis;
    if (basis.family == HigherOrderBasisFamily::NODE_GRADIENT)
    {
      if (on_face(basis.nodes[0]))
      {
        add(TetrahedronFaceSingularityType::NODE, {basis.nodes[0], -1}, basis.nu);
      }
    }
    else if (basis.family == HigherOrderBasisFamily::EDGE_GRADIENT)
    {
      if (on_face(basis.nodes[0]) && on_face(basis.nodes[1]))
      {
        add(TetrahedronFaceSingularityType::EDGE, {basis.nodes[0], basis.nodes[1]},
            basis.nu);
      }
    }
  }
  std::sort(result.begin(), result.end(),
            [](const auto &left, const auto &right)
            {
              if (left.type != right.type)
              {
                return left.type < right.type;
              }
              return left.nodes < right.nodes;
            });
  return result;
}

H1FieldValue EnrichedH1FieldEvaluator::EvaluateBarycentric(
    int element, const mfem::IntegrationPoint &point, const BarycentricPoint &lambda)
{
  if (!initialized)
  {
    throw std::logic_error("Singular H1 field evaluator has no combined true-DOF vector!");
  }
  const auto &element_dofs = GetParallelElementDofs(
      element, fespace.GetNE(), topology.elements, face_neighbor_enrichment);
  auto &transformation = GetParallelElementTransformation(fespace, element);
  double jacobian_determinant;
  const auto grad_lambda =
      GetBarycentricGradients(transformation, point, jacobian_determinant);

  transformation.SetIntPoint(&point);
  H1FieldValue result;
  result.potential = standard_field.GetValue(transformation, point);
  mfem::Vector standard_gradient(3);
  GetStandardGradient(standard_field, fespace, element, transformation, point,
                      standard_gradient);
  for (int d = 0; d < 3; d++)
  {
    result.gradient[d] = standard_gradient[d];
  }

  const auto &coefficients =
      element < fespace.GetNE() ? local_enrichment : face_neighbor_enrichment.coefficients;
  const auto enrichment =
      EvaluateElementH1Enrichment(element_dofs, coefficients, lambda, grad_lambda);
  result.potential += enrichment.potential;
  for (int d = 0; d < 3; d++)
  {
    result.gradient[d] += enrichment.gradient[d];
  }
  return result;
}

AdaptiveQuadratureResult EnrichedH1FieldEvaluator::IntegrateElementGradientEnergy(
    int element, double electric_coefficient, const AdaptiveAssemblyOptions &options)
{
  if (!initialized)
  {
    throw std::logic_error("Singular H1 energy evaluator has no combined true-DOF vector!");
  }
  if (!std::isfinite(electric_coefficient) || !(electric_coefficient > 0.0))
  {
    throw std::invalid_argument(
        "Singular H1 energy integration requires a positive electric coefficient!");
  }
  const bool valid_fixed = options.fixed_subdivision && options.quadrature_order >= 1 &&
                           options.subdivisions >= 0 && options.subdivisions <= 8;
  const bool valid_adaptive =
      !options.fixed_subdivision && options.quadrature_order >= 1 &&
      std::isfinite(options.absolute_tolerance) && options.absolute_tolerance >= 0.0 &&
      std::isfinite(options.relative_tolerance) && options.relative_tolerance >= 0.0 &&
      (options.absolute_tolerance > 0.0 || options.relative_tolerance > 0.0) &&
      options.maximum_subdivisions >= 1;
  if (!valid_fixed && !valid_adaptive)
  {
    throw std::invalid_argument(
        "Singular H1 energy integration has invalid quadrature options!");
  }
  if (element < 0 || element >= fespace.GetNE())
  {
    throw std::out_of_range("Singular H1 energy integration element is out of range!");
  }
  auto *transformation = fespace.GetElementTransformation(element);
  if (!transformation)
  {
    throw std::runtime_error("Singular H1 energy element has no transformation!");
  }
  if (topology.elements[element].h1.empty())
  {
    const auto *finite_element = fespace.GetFE(element);
    if (!finite_element || finite_element->GetGeomType() != mfem::Geometry::TETRAHEDRON ||
        finite_element->GetMapType() != mfem::FiniteElement::VALUE)
    {
      throw std::runtime_error(
          "Singular H1 energy integration requires a tetrahedral value element!");
    }
    const int integration_order =
        2 * std::max(0, finite_element->GetOrder() - 1) + transformation->OrderW();
    const auto &rule = mfem::IntRules.Get(mfem::Geometry::TETRAHEDRON, integration_order);
    mfem::Vector gradient(3);
    double value = 0.0;
    for (int q = 0; q < rule.GetNPoints(); q++)
    {
      const auto &point = rule.IntPoint(q);
      transformation->SetIntPoint(&point);
      standard_field.GetGradient(*transformation, gradient);
      value += point.weight * electric_coefficient * transformation->Weight() *
               (gradient * gradient);
    }
    return {value, 0.0, 1, 0, true};
  }

  mfem::IntegrationPoint point;
  const auto integrand = [&](const BarycentricPoint &lambda)
  {
    point.Set3(lambda[1], lambda[2], lambda[3]);
    double jacobian_determinant;
    (void)GetBarycentricGradients(*transformation, point, jacobian_determinant);
    const auto value = Evaluate(element, point);
    return electric_coefficient * jacobian_determinant *
           Dot(value.gradient, value.gradient);
  };
  if (!options.fixed_subdivision)
  {
    return IntegrateReferenceTetrahedronAdaptive(
        options.quadrature_order, options.absolute_tolerance, options.relative_tolerance,
        options.maximum_subdivisions, integrand);
  }
  std::size_t leaves = 1;
  for (int level = 0; level < options.subdivisions; level++)
  {
    leaves *= 8;
  }
  return {IntegrateReferenceTetrahedron(options.quadrature_order, options.subdivisions,
                                        integrand),
          0.0, leaves, options.subdivisions, true};
}

void EnrichedH1FieldEvaluator::ProjectToDiscontinuousGridFunctions(
    mfem::ParGridFunction &potential, mfem::ParGridFunction &electric_field)
{
  const auto *potential_space = potential.ParFESpace();
  const auto *electric_space = electric_field.ParFESpace();
  if (!initialized || !potential_space || !electric_space ||
      potential_space->GetParMesh() != fespace.GetParMesh() ||
      electric_space->GetParMesh() != fespace.GetParMesh() ||
      potential_space->GetVDim() != 1 || electric_space->GetVDim() != 3 ||
      potential_space->FEColl()->GetMapType(3) != mfem::FiniteElement::VALUE ||
      electric_space->FEColl()->GetMapType(3) != mfem::FiniteElement::VALUE)
  {
    throw std::invalid_argument(
        "Singular field sampling requires initialized scalar and vector value spaces "
        "on the solve mesh!");
  }

  PotentialCoefficient potential_coefficient(*this);
  ElectricFieldCoefficient electric_coefficient(*this);
  potential.ProjectCoefficient(potential_coefficient);
  electric_field.ProjectCoefficient(electric_coefficient);
}

std::vector<H1CoefficientDiagnostic>
EnrichedH1FieldEvaluator::GetOwnedCoefficientDiagnostics() const
{
  if (!initialized)
  {
    throw std::logic_error(
        "Singular H1 coefficient diagnostics require a combined true-DOF vector!");
  }

  std::vector<H1CoefficientDiagnostic> diagnostics;
  diagnostics.reserve(static_cast<std::size_t>(numbering.h1.owned_size));
  const int rank = Mpi::Rank(fespace.GetComm());
  for (std::size_t local = 0; local < topology.h1_dofs.size(); local++)
  {
    if (numbering.h1.owner[local] != rank)
    {
      continue;
    }
    const HYPRE_BigInt true_dof = numbering.h1.local_to_true[local];
    if (true_dof < numbering.h1.owned_offset ||
        true_dof >= numbering.h1.owned_offset + numbering.h1.owned_size ||
        !std::isfinite(local_enrichment[static_cast<int>(local)]) ||
        !std::isfinite(h1_exponents[local]))
    {
      throw std::runtime_error(
          "Singular H1 coefficient diagnostics found inconsistent owned data!");
    }
    diagnostics.push_back({true_dof, topology.h1_dofs[local], h1_exponents[local],
                           local_enrichment[static_cast<int>(local)]});
  }
  std::sort(diagnostics.begin(), diagnostics.end(), [](const auto &left, const auto &right)
            { return left.true_dof < right.true_dof; });
  if (diagnostics.size() != static_cast<std::size_t>(numbering.h1.owned_size))
  {
    throw std::runtime_error(
        "Singular H1 coefficient diagnostics do not cover every owned true DOF!");
  }
  return diagnostics;
}

std::vector<H1EdgeSlopeDiagnostic> EnrichedH1FieldEvaluator::FitEdgeSlopes(
    const FeatureTopology &features, const std::vector<GlobalVertexId> &source_vertex_ids,
    const std::vector<GlobalVertexId> &source_element_ids, const EdgeSlopeOptions &options)
{
  if (!initialized)
  {
    throw std::logic_error(
        "Singular edge slope diagnostics require a combined true-DOF vector!");
  }
  if (features.elements.size() != static_cast<std::size_t>(fespace.GetNE()) ||
      source_vertex_ids.size() != static_cast<std::size_t>(fespace.GetParMesh()->GetNV()) ||
      source_element_ids.size() != static_cast<std::size_t>(fespace.GetNE()) ||
      options.sample_count < 3 || !std::isfinite(options.minimum_barycentric_radius) ||
      !std::isfinite(options.maximum_barycentric_radius) ||
      !(options.minimum_barycentric_radius > 0.0) ||
      !(options.minimum_barycentric_radius < options.maximum_barycentric_radius) ||
      !(options.maximum_barycentric_radius < 1.0))
  {
    throw std::invalid_argument(
        "Singular edge slope diagnostics received inconsistent topology or options!");
  }

  const double log_minimum_radius = std::log(options.minimum_barycentric_radius);
  const double log_maximum_radius = std::log(options.maximum_barycentric_radius);
  std::vector<H1EdgeSlopeDiagnostic> diagnostics;
  mfem::IntegrationPoint point;
  mfem::Vector physical_point(3);
  for (int element = 0; element < fespace.GetNE(); element++)
  {
    const auto *tetrahedron = fespace.GetParMesh()->GetElement(element);
    auto *transformation = fespace.GetElementTransformation(element);
    if (!tetrahedron || tetrahedron->GetGeometryType() != mfem::Geometry::TETRAHEDRON ||
        !transformation)
    {
      throw std::runtime_error(
          "Singular edge slope diagnostics require tetrahedral elements!");
    }
    for (const auto &edge : features.elements[element].edges)
    {
      if (edge.feature >= features.features.size() ||
          edge.segment >= features.segments.size() ||
          features.segments[edge.segment].feature != edge.feature)
      {
        throw std::invalid_argument(
            "Singular edge slope diagnostics found invalid feature incidence!");
      }
      std::array<bool, 4> seen{false, false, false, false};
      std::array<GlobalVertexId, 4> canonical_vertices;
      for (int local = 0; local < 4; local++)
      {
        const int node = edge.canonical_nodes[local];
        if (node < 0 || node >= 4 || seen[node])
        {
          throw std::invalid_argument(
              "Singular edge slope diagnostics found invalid canonical nodes!");
        }
        seen[node] = true;
        const int mesh_vertex = tetrahedron->GetVertices()[node];
        if (mesh_vertex < 0 || mesh_vertex >= static_cast<int>(source_vertex_ids.size()))
        {
          throw std::invalid_argument(
              "Singular edge slope diagnostics found invalid mesh vertices!");
        }
        canonical_vertices[local] = source_vertex_ids[mesh_vertex];
      }
      std::array<GlobalVertexId, 2> incidence_edge{canonical_vertices[0],
                                                   canonical_vertices[1]};
      std::array<GlobalVertexId, 2> source_edge{
          features.segments[edge.segment].mesh_vertices[0],
          features.segments[edge.segment].mesh_vertices[1]};
      std::sort(incidence_edge.begin(), incidence_edge.end());
      std::sort(source_edge.begin(), source_edge.end());
      if (incidence_edge != source_edge)
      {
        throw std::invalid_argument(
            "Singular edge slope diagnostic incidence does not match its source segment!");
      }

      mfem::Vector first(3), second(3);
      point.Set3(edge.canonical_nodes[0] == 1 ? 1.0 : 0.0,
                 edge.canonical_nodes[0] == 2 ? 1.0 : 0.0,
                 edge.canonical_nodes[0] == 3 ? 1.0 : 0.0);
      transformation->Transform(point, first);
      point.Set3(edge.canonical_nodes[1] == 1 ? 1.0 : 0.0,
                 edge.canonical_nodes[1] == 2 ? 1.0 : 0.0,
                 edge.canonical_nodes[1] == 3 ? 1.0 : 0.0);
      transformation->Transform(point, second);
      std::vector<double> log_distance(options.sample_count);
      std::vector<double> log_field_norm(options.sample_count);
      std::vector<double> distance(options.sample_count);
      std::vector<double> field_norm(options.sample_count);
      bool valid = true;
      for (int sample = 0; sample < options.sample_count; sample++)
      {
        const double fraction =
            static_cast<double>(sample) / static_cast<double>(options.sample_count - 1);
        const double radius = std::exp(
            log_minimum_radius + fraction * (log_maximum_radius - log_minimum_radius));
        BarycentricPoint lambda;
        lambda[edge.canonical_nodes[0]] = 0.5 * (1.0 - radius);
        lambda[edge.canonical_nodes[1]] = 0.5 * (1.0 - radius);
        lambda[edge.canonical_nodes[2]] = 0.5 * radius;
        lambda[edge.canonical_nodes[3]] = 0.5 * radius;
        point.Set3(lambda[1], lambda[2], lambda[3]);
        const auto field = Evaluate(element, point);
        transformation->Transform(point, physical_point);
        distance[sample] = DistanceToLine(physical_point, first, second);
        field_norm[sample] = std::sqrt(Dot(field.gradient, field.gradient));
        if (!std::isfinite(distance[sample]) || !(distance[sample] > 0.0) ||
            !std::isfinite(field_norm[sample]) || !(field_norm[sample] > 0.0))
        {
          valid = false;
          break;
        }
        log_distance[sample] = std::log(distance[sample]);
        log_field_norm[sample] = std::log(field_norm[sample]);
      }

      double fitted_slope = 0.0;
      double r_squared = 0.0;
      if (valid)
      {
        long double mean_distance = 0.0L;
        long double mean_field = 0.0L;
        for (int sample = 0; sample < options.sample_count; sample++)
        {
          mean_distance += log_distance[sample];
          mean_field += log_field_norm[sample];
        }
        mean_distance /= options.sample_count;
        mean_field /= options.sample_count;
        long double distance_variance = 0.0L;
        long double covariance = 0.0L;
        for (int sample = 0; sample < options.sample_count; sample++)
        {
          const long double centered_distance = log_distance[sample] - mean_distance;
          distance_variance += centered_distance * centered_distance;
          covariance += centered_distance * (log_field_norm[sample] - mean_field);
        }
        if (!(distance_variance > 0.0L))
        {
          valid = false;
        }
        else
        {
          fitted_slope = static_cast<double>(covariance / distance_variance);
          const long double intercept = mean_field - fitted_slope * mean_distance;
          long double residual_sum = 0.0L;
          long double total_sum = 0.0L;
          for (int sample = 0; sample < options.sample_count; sample++)
          {
            const long double residual =
                log_field_norm[sample] - (intercept + fitted_slope * log_distance[sample]);
            const long double centered = log_field_norm[sample] - mean_field;
            residual_sum += residual * residual;
            total_sum += centered * centered;
          }
          r_squared = (total_sum > 0.0L)
                          ? static_cast<double>(1.0L - residual_sum / total_sum)
                          : ((residual_sum == 0.0L) ? 1.0 : 0.0);
          valid = std::isfinite(fitted_slope) && std::isfinite(r_squared);
        }
      }

      const double exponent = features.features[edge.feature].nu;
      if (!std::isfinite(exponent) || !(exponent > 0.0) || !(exponent < 1.0))
      {
        throw std::invalid_argument(
            "Singular edge slope diagnostics found an invalid feature exponent!");
      }
      diagnostics.push_back(
          {source_element_ids[element], edge.feature, edge.segment, canonical_vertices,
           options.sample_count, exponent, exponent - 1.0, fitted_slope, r_squared,
           valid ? distance.front() : 0.0, valid ? distance.back() : 0.0,
           valid ? field_norm.front() : 0.0, valid ? field_norm.back() : 0.0, valid});
    }
  }
  return diagnostics;
}

EnrichedNDFieldEvaluator::EnrichedNDFieldEvaluator(const DofTopology &topology,
                                                   const ParallelDofNumbering &numbering,
                                                   mfem::ParFiniteElementSpace &fespace)
  : topology(topology), numbering(numbering), fespace(fespace), standard_field(&fespace),
    enrichment_prolongation(
        BuildParallelEnrichmentProlongation(fespace.GetComm(), numbering.nd)),
    nd_exponents(topology.nd_dofs.size(), std::numeric_limits<double>::quiet_NaN()),
    initialized(false)
{
  bool valid =
      fespace.GetParMesh() && fespace.GetParMesh()->Dimension() == 3 &&
      fespace.GetParMesh()->SpaceDimension() == 3 &&
      fespace.FEColl()->GetMapType(3) == mfem::FiniteElement::H_CURL &&
      topology.elements.size() == static_cast<std::size_t>(fespace.GetNE()) &&
      topology.nd_dofs.size() == numbering.nd.local_to_true.size() &&
      numbering.nd.local_size == static_cast<HYPRE_BigInt>(topology.nd_dofs.size()) &&
      numbering.nd.owned_size <=
          static_cast<HYPRE_BigInt>(std::numeric_limits<int>::max()) &&
      enrichment_prolongation &&
      enrichment_prolongation->Height() == static_cast<int>(topology.nd_dofs.size()) &&
      enrichment_prolongation->Width() == static_cast<int>(numbering.nd.owned_size);
  if (valid)
  {
    for (const auto &element : topology.elements)
    {
      for (const auto &element_dof : element.nd)
      {
        if (element_dof.dof >= nd_exponents.size() ||
            topology.nd_dofs[element_dof.dof].family != element_dof.basis.family ||
            topology.nd_dofs[element_dof.dof].order != element_dof.basis.order ||
            !std::isfinite(element_dof.basis.nu) || !(element_dof.basis.nu > 0.0) ||
            !(element_dof.basis.nu < 1.0))
        {
          valid = false;
          break;
        }
        double &exponent = nd_exponents[element_dof.dof];
        if (std::isnan(exponent))
        {
          exponent = element_dof.basis.nu;
        }
        else if (exponent != element_dof.basis.nu)
        {
          valid = false;
          break;
        }
      }
      if (!valid)
      {
        break;
      }
    }
    valid = valid && std::all_of(nd_exponents.begin(), nd_exponents.end(),
                                 [](double exponent) { return std::isfinite(exponent); });
  }
  Mpi::GlobalAnd(1, &valid, fespace.GetComm());
  if (!valid)
  {
    throw std::invalid_argument(
        "Singular ND field evaluator received inconsistent space and DOF topology!");
  }
  local_enrichment.SetSize(static_cast<int>(topology.nd_dofs.size()));
  InitializeFaceNeighborEnrichment(fespace, topology.elements, false,
                                   face_neighbor_enrichment);
  const std::size_t parallel_element_count =
      topology.elements.size() + face_neighbor_enrichment.element_dofs.size();
  barycentric_gradient_cache.resize(parallel_element_count);
  standard_element_coefficients.resize(parallel_element_count);
  standard_element_coefficients_valid.assign(parallel_element_count, false);
}

void EnrichedNDFieldEvaluator::SetFromTrueDofs(const mfem::Vector &combined_true_dofs)
{
  const int standard_size = fespace.GetTrueVSize();
  const int enrichment_size = static_cast<int>(numbering.nd.owned_size);
  bool valid = combined_true_dofs.Size() == standard_size + enrichment_size;
  Mpi::GlobalAnd(1, &valid, fespace.GetComm());
  if (!valid)
  {
    throw std::invalid_argument(
        "Combined singular ND vector has inconsistent process-local dimensions!");
  }

  mfem::Vector standard_true_dofs(standard_size);
  mfem::Vector enrichment_true_dofs(enrichment_size);
  for (int i = 0; i < standard_size; i++)
  {
    standard_true_dofs[i] = combined_true_dofs[i];
  }
  for (int i = 0; i < enrichment_size; i++)
  {
    enrichment_true_dofs[i] = combined_true_dofs[standard_size + i];
  }
  standard_field.SetFromTrueDofs(standard_true_dofs);
  standard_field.ExchangeFaceNbrData();
  enrichment_prolongation->Mult(enrichment_true_dofs, local_enrichment);
  ExchangeFaceNeighborEnrichment(topology.elements, local_enrichment, false, fespace,
                                 face_neighbor_enrichment);
  std::fill(standard_element_coefficients_valid.begin(),
            standard_element_coefficients_valid.end(), false);
  initialized = true;
}

NDFieldValue EnrichedNDFieldEvaluator::Evaluate(int element,
                                                const mfem::IntegrationPoint &point)
{
  ValidateInteriorPoint(point);
  const BarycentricPoint lambda{1.0 - point.x - point.y - point.z, point.x, point.y,
                                point.z};
  return EvaluateBarycentric(element, point, lambda);
}

NDFieldValue EnrichedNDFieldEvaluator::EvaluateClosure(int element,
                                                       const mfem::IntegrationPoint &point)
{
  ValidateClosurePoint(point);
  const BarycentricPoint lambda{1.0 - point.x - point.y - point.z, point.x, point.y,
                                point.z};
  return EvaluateBarycentric(element, point, lambda);
}

Vector3 EnrichedNDFieldEvaluator::EvaluateValueClosure(int element,
                                                       const mfem::IntegrationPoint &point)
{
  ValidateClosurePoint(point);
  const BarycentricPoint lambda{1.0 - point.x - point.y - point.z, point.x, point.y,
                                point.z};
  return EvaluateValueBarycentric(element, point, lambda);
}

NDFieldValuePair EnrichedNDFieldEvaluator::EvaluateValueClosurePair(
    EnrichedNDFieldEvaluator &second, int element, const mfem::IntegrationPoint &point)
{
  ValidateClosurePoint(point);
  const BarycentricPoint lambda{1.0 - point.x - point.y - point.z, point.x, point.y,
                                point.z};
  return EvaluateValueBarycentricPair(second, element, point, lambda);
}

void EnrichedNDFieldEvaluator::EvaluateValueClosureBatch(
    const std::vector<std::pair<EnrichedNDFieldEvaluator *, EnrichedNDFieldEvaluator *>>
        &evaluators,
    int element, const mfem::IntegrationPoint &point, std::vector<NDFieldValuePair> &values)
{
  ValidateClosurePoint(point);
  if (evaluators.empty() || values.size() != evaluators.size())
  {
    throw std::invalid_argument(
        "Batched singular ND field evaluation requires matching nonempty input and "
        "output!");
  }
  for (const auto &[real, imaginary] : evaluators)
  {
    if (!real || !imaginary || !real->initialized || !imaginary->initialized ||
        &topology != &real->topology || &topology != &imaginary->topology ||
        &numbering != &real->numbering || &numbering != &imaginary->numbering ||
        &fespace != &real->fespace || &fespace != &imaginary->fespace)
    {
      throw std::invalid_argument(
          "Batched singular ND field evaluation requires initialized matching spaces "
          "and topology!");
    }
  }

  const BarycentricPoint lambda{1.0 - point.x - point.y - point.z, point.x, point.y,
                                point.z};
  const auto &element_dofs = GetParallelElementDofs(
      element, fespace.GetNE(), topology.elements, face_neighbor_enrichment);
  auto &transformation = GetParallelElementTransformation(fespace, element);
  const auto grad_lambda = GetValueBarycentricGradients(element, transformation, point);
  const int local_elements = fespace.GetNE();
  const auto *finite_element = element < local_elements
                                   ? fespace.GetFE(element)
                                   : fespace.GetFaceNbrFE(element - local_elements);
  if (!finite_element || finite_element->GetRangeType() != mfem::FiniteElement::VECTOR ||
      finite_element->GetDof() <= 0 || finite_element->GetRangeDim() > 3)
  {
    throw std::runtime_error(
        "Batched singular ND field evaluation received an incompatible standard "
        "element!");
  }
  transformation.SetIntPoint(&point);
  standard_vector_shape.SetSize(finite_element->GetDof(), 3);
  finite_element->CalcVShape(transformation, standard_vector_shape);
  for (std::size_t field = 0; field < evaluators.size(); field++)
  {
    auto &[real, imaginary] = evaluators[field];
    double real_data[3], imaginary_data[3];
    mfem::Vector real_standard(real_data, 3), imaginary_standard(imaginary_data, 3);
    standard_vector_shape.MultTranspose(real->GetStandardElementCoefficients(element),
                                        real_standard);
    standard_vector_shape.MultTranspose(imaginary->GetStandardElementCoefficients(element),
                                        imaginary_standard);
    values[field] = {{real_standard[0], real_standard[1], real_standard[2]},
                     {imaginary_standard[0], imaginary_standard[1], imaginary_standard[2]}};
  }

  for (const auto &element_dof : element_dofs.nd)
  {
    const auto value =
        EvaluateHigherOrderBasisValue(lambda, grad_lambda, element_dof.basis);
    for (std::size_t field = 0; field < evaluators.size(); field++)
    {
      auto &[real, imaginary] = evaluators[field];
      const auto &real_coefficients = element < fespace.GetNE()
                                          ? real->local_enrichment
                                          : real->face_neighbor_enrichment.coefficients;
      const auto &imaginary_coefficients =
          element < fespace.GetNE() ? imaginary->local_enrichment
                                    : imaginary->face_neighbor_enrichment.coefficients;
      if (element_dof.dof >= static_cast<std::size_t>(real_coefficients.Size()) ||
          element_dof.dof >= static_cast<std::size_t>(imaginary_coefficients.Size()))
      {
        throw std::invalid_argument(
            "Batched singular ND field evaluation received an invalid element DOF map!");
      }
      const double real_coefficient = real_coefficients[element_dof.dof];
      const double imaginary_coefficient = imaginary_coefficients[element_dof.dof];
      if (!std::isfinite(real_coefficient) || !std::isfinite(imaginary_coefficient))
      {
        throw std::invalid_argument(
            "Batched singular ND field evaluation received a nonfinite coefficient!");
      }
      for (int d = 0; d < 3; d++)
      {
        values[field].first[d] += real_coefficient * value[d];
        values[field].second[d] += imaginary_coefficient * value[d];
      }
    }
  }
}

std::vector<TetrahedronFaceSingularity>
EnrichedNDFieldEvaluator::GetElementFaceSingularities(
    int element, const std::array<int, 3> &face_nodes) const
{
  const auto &element_dofs = GetParallelElementDofs(
      element, fespace.GetNE(), topology.elements, face_neighbor_enrichment);
  auto sorted_face = face_nodes;
  std::sort(sorted_face.begin(), sorted_face.end());
  if (sorted_face[0] < 0 || sorted_face[2] >= 4 ||
      std::adjacent_find(sorted_face.begin(), sorted_face.end()) != sorted_face.end())
  {
    throw std::invalid_argument(
        "Singular ND face-singularity query received invalid face nodes!");
  }
  const auto on_face = [&sorted_face](int node)
  { return std::binary_search(sorted_face.begin(), sorted_face.end(), node); };

  std::vector<TetrahedronFaceSingularity> result;
  const auto add =
      [&result](TetrahedronFaceSingularityType type, std::array<int, 2> nodes, double nu)
  {
    if (type == TetrahedronFaceSingularityType::EDGE && nodes[1] < nodes[0])
    {
      std::swap(nodes[0], nodes[1]);
    }
    const auto existing =
        std::find_if(result.begin(), result.end(), [type, &nodes](const auto &feature)
                     { return feature.type == type && feature.nodes == nodes; });
    if (existing == result.end())
    {
      result.push_back({type, nodes, nu});
    }
    else if (existing->nu != nu)
    {
      throw std::runtime_error(
          "Singular ND face contains inconsistent exponents for one feature!");
    }
  };

  for (const auto &element_dof : element_dofs.nd)
  {
    const auto &basis = element_dof.basis;
    if (basis.family == HigherOrderBasisFamily::NODE_GRADIENT && on_face(basis.nodes[0]))
    {
      add(TetrahedronFaceSingularityType::NODE, {basis.nodes[0], -1}, basis.nu);
    }
    else if (basis.family == HigherOrderBasisFamily::EDGE_GRADIENT &&
             on_face(basis.nodes[0]) && on_face(basis.nodes[1]))
    {
      add(TetrahedronFaceSingularityType::EDGE, {basis.nodes[0], basis.nodes[1]}, basis.nu);
    }
  }
  std::sort(result.begin(), result.end(),
            [](const auto &left, const auto &right)
            {
              if (left.type != right.type)
              {
                return left.type < right.type;
              }
              return left.nodes < right.nodes;
            });
  return result;
}

NDFieldValue EnrichedNDFieldEvaluator::EvaluateBarycentric(
    int element, const mfem::IntegrationPoint &point, const BarycentricPoint &lambda)
{
  if (!initialized)
  {
    throw std::logic_error("Singular ND field evaluator has no combined true-DOF vector!");
  }
  const auto &element_dofs = GetParallelElementDofs(
      element, fespace.GetNE(), topology.elements, face_neighbor_enrichment);
  auto &transformation = GetParallelElementTransformation(fespace, element);
  const auto grad_lambda = GetValueBarycentricGradients(element, transformation, point);

  transformation.SetIntPoint(&point);
  NDFieldValue result;
  mfem::Vector standard_value(3), standard_curl(3);
  standard_field.GetVectorValue(transformation, point, standard_value);
  transformation.SetIntPoint(&point);
  GetStandardCurl(standard_field, fespace, element, transformation, point, standard_curl);
  for (int d = 0; d < 3; d++)
  {
    result.value[d] = standard_value[d];
    result.curl[d] = standard_curl[d];
  }

  const auto &coefficients =
      element < fespace.GetNE() ? local_enrichment : face_neighbor_enrichment.coefficients;
  const auto enrichment =
      EvaluateElementNDEnrichment(element_dofs, coefficients, lambda, grad_lambda);
  for (int d = 0; d < 3; d++)
  {
    result.value[d] += enrichment.value[d];
    result.curl[d] += enrichment.curl[d];
  }
  return result;
}

const mfem::Vector &EnrichedNDFieldEvaluator::GetStandardElementCoefficients(int element)
{
  if (element < 0 || element >= static_cast<int>(standard_element_coefficients.size()))
  {
    throw std::out_of_range(
        "Singular ND value evaluation received an invalid element index!");
  }
  if (standard_element_coefficients_valid[element])
  {
    return standard_element_coefficients[element];
  }

  mfem::Array<int> dofs;
  mfem::DofTransformation dof_transformation;
  const int local_elements = fespace.GetNE();
  const mfem::FiniteElement *finite_element;
  auto &coefficients = standard_element_coefficients[element];
  if (element < local_elements)
  {
    fespace.GetElementVDofs(element, dofs, dof_transformation);
    standard_field.GetSubVector(dofs, coefficients);
    finite_element = fespace.GetFE(element);
  }
  else
  {
    const int neighbor = element - local_elements;
    fespace.GetFaceNbrElementVDofs(neighbor, dofs, dof_transformation);
    standard_field.FaceNbrData().GetSubVector(dofs, coefficients);
    finite_element = fespace.GetFaceNbrFE(neighbor);
  }
  dof_transformation.InvTransformPrimal(coefficients);
  if (!finite_element || finite_element->GetRangeType() != mfem::FiniteElement::VECTOR ||
      finite_element->GetDof() != coefficients.Size())
  {
    throw std::runtime_error(
        "Singular ND value evaluation received an incompatible standard element!");
  }
  standard_element_coefficients_valid[element] = true;
  return coefficients;
}

BarycentricGradients EnrichedNDFieldEvaluator::GetValueBarycentricGradients(
    int element, mfem::ElementTransformation &transformation,
    const mfem::IntegrationPoint &point)
{
  if (element < 0 || element >= static_cast<int>(barycentric_gradient_cache.size()))
  {
    throw std::out_of_range(
        "Singular ND value evaluation received an invalid element index!");
  }
  auto &cache = barycentric_gradient_cache[element];
  if (!cache.initialized)
  {
    cache.constant = transformation.OrderJ() == 0;
    cache.initialized = true;
    if (cache.constant)
    {
      double jacobian_determinant;
      cache.gradients =
          GetBarycentricGradients(transformation, point, jacobian_determinant);
    }
  }
  if (cache.constant)
  {
    return cache.gradients;
  }
  double jacobian_determinant;
  return GetBarycentricGradients(transformation, point, jacobian_determinant);
}

Vector3 EnrichedNDFieldEvaluator::EvaluateValueBarycentric(
    int element, const mfem::IntegrationPoint &point, const BarycentricPoint &lambda)
{
  if (!initialized)
  {
    throw std::logic_error("Singular ND field evaluator has no combined true-DOF vector!");
  }
  const auto &element_dofs = GetParallelElementDofs(
      element, fespace.GetNE(), topology.elements, face_neighbor_enrichment);
  auto &transformation = GetParallelElementTransformation(fespace, element);
  const auto grad_lambda = GetValueBarycentricGradients(element, transformation, point);

  Vector3 result{};
  transformation.SetIntPoint(&point);
  const int local_elements = fespace.GetNE();
  const auto *finite_element = element < local_elements
                                   ? fespace.GetFE(element)
                                   : fespace.GetFaceNbrFE(element - local_elements);
  if (!finite_element || finite_element->GetRangeType() != mfem::FiniteElement::VECTOR ||
      finite_element->GetDof() <= 0 || finite_element->GetRangeDim() > 3)
  {
    throw std::runtime_error(
        "Singular ND value evaluation received an incompatible standard element!");
  }
  standard_vector_shape.SetSize(finite_element->GetDof(), 3);
  finite_element->CalcVShape(transformation, standard_vector_shape);
  mfem::Vector standard_value(3);
  standard_vector_shape.MultTranspose(GetStandardElementCoefficients(element),
                                      standard_value);
  result = {standard_value[0], standard_value[1], standard_value[2]};

  const auto &coefficients =
      element < fespace.GetNE() ? local_enrichment : face_neighbor_enrichment.coefficients;
  const auto enrichment =
      EvaluateElementNDEnrichmentValue(element_dofs, coefficients, lambda, grad_lambda);
  for (int d = 0; d < 3; d++)
  {
    result[d] += enrichment[d];
  }
  return result;
}

NDFieldValuePair EnrichedNDFieldEvaluator::EvaluateValueBarycentricPair(
    EnrichedNDFieldEvaluator &second, int element, const mfem::IntegrationPoint &point,
    const BarycentricPoint &lambda)
{
  if (!initialized || !second.initialized)
  {
    throw std::logic_error(
        "Paired singular ND field evaluation requires two initialized fields!");
  }
  if (&topology != &second.topology || &numbering != &second.numbering ||
      &fespace != &second.fespace)
  {
    throw std::invalid_argument(
        "Paired singular ND field evaluation requires matching spaces and topology!");
  }
  const auto &element_dofs = GetParallelElementDofs(
      element, fespace.GetNE(), topology.elements, face_neighbor_enrichment);
  auto &transformation = GetParallelElementTransformation(fespace, element);
  const auto grad_lambda = GetValueBarycentricGradients(element, transformation, point);

  NDFieldValuePair result;
  transformation.SetIntPoint(&point);
  const int local_elements = fespace.GetNE();
  const auto *finite_element = element < local_elements
                                   ? fespace.GetFE(element)
                                   : fespace.GetFaceNbrFE(element - local_elements);
  if (!finite_element || finite_element->GetRangeType() != mfem::FiniteElement::VECTOR ||
      finite_element->GetDof() <= 0 || finite_element->GetRangeDim() > 3)
  {
    throw std::runtime_error(
        "Paired singular ND field evaluation received an incompatible standard "
        "element!");
  }
  standard_vector_shape.SetSize(finite_element->GetDof(), 3);
  finite_element->CalcVShape(transformation, standard_vector_shape);
  mfem::Vector first_standard(3), second_standard(3);
  standard_vector_shape.MultTranspose(GetStandardElementCoefficients(element),
                                      first_standard);
  standard_vector_shape.MultTranspose(second.GetStandardElementCoefficients(element),
                                      second_standard);
  result = {{first_standard[0], first_standard[1], first_standard[2]},
            {second_standard[0], second_standard[1], second_standard[2]}};

  const auto &first_coefficients =
      element < fespace.GetNE() ? local_enrichment : face_neighbor_enrichment.coefficients;
  const auto &second_coefficients = element < fespace.GetNE()
                                        ? second.local_enrichment
                                        : second.face_neighbor_enrichment.coefficients;
  for (const auto &element_dof : element_dofs.nd)
  {
    if (element_dof.dof >= static_cast<std::size_t>(first_coefficients.Size()) ||
        element_dof.dof >= static_cast<std::size_t>(second_coefficients.Size()))
    {
      throw std::invalid_argument(
          "Paired singular ND field evaluation received an invalid element DOF map!");
    }
    const double first_coefficient = first_coefficients[element_dof.dof];
    const double second_coefficient = second_coefficients[element_dof.dof];
    if (!std::isfinite(first_coefficient) || !std::isfinite(second_coefficient))
    {
      throw std::invalid_argument(
          "Paired singular ND field evaluation received a nonfinite coefficient!");
    }
    const auto value =
        EvaluateHigherOrderBasisValue(lambda, grad_lambda, element_dof.basis);
    for (int d = 0; d < 3; d++)
    {
      result.first[d] += first_coefficient * value[d];
      result.second[d] += second_coefficient * value[d];
    }
  }
  return result;
}

TriangleEnrichedH1FieldEvaluator::TriangleEnrichedH1FieldEvaluator(
    const TriangleDofTopology &topology, const ParallelDofNumbering &numbering,
    mfem::ParFiniteElementSpace &fespace)
  : topology(topology), numbering(numbering), fespace(fespace), standard_field(&fespace),
    enrichment_prolongation(
        BuildParallelEnrichmentProlongation(fespace.GetComm(), numbering.h1)),
    h1_exponents(topology.h1_dofs.size(), std::numeric_limits<double>::quiet_NaN()),
    initialized(false)
{
  bool valid =
      fespace.GetParMesh() && fespace.GetParMesh()->Dimension() == 2 &&
      fespace.GetParMesh()->SpaceDimension() == 2 &&
      topology.elements.size() == static_cast<std::size_t>(fespace.GetNE()) &&
      topology.h1_dofs.size() == numbering.h1.local_to_true.size() &&
      numbering.h1.local_size == static_cast<HYPRE_BigInt>(topology.h1_dofs.size()) &&
      numbering.h1.owned_size <=
          static_cast<HYPRE_BigInt>(std::numeric_limits<int>::max()) &&
      enrichment_prolongation &&
      enrichment_prolongation->Height() == static_cast<int>(topology.h1_dofs.size()) &&
      enrichment_prolongation->Width() == static_cast<int>(numbering.h1.owned_size);
  if (valid)
  {
    for (const auto &element : topology.elements)
    {
      for (const auto &element_dof : element.h1)
      {
        if (element_dof.dof >= h1_exponents.size() ||
            element_dof.basis.family != HigherOrderBasisFamily::NODE_GRADIENT ||
            element_dof.basis.order != 1 ||
            topology.h1_dofs[element_dof.dof].family != element_dof.basis.family ||
            topology.h1_dofs[element_dof.dof].order != element_dof.basis.order ||
            !std::isfinite(element_dof.basis.nu) || !(element_dof.basis.nu > 0.0) ||
            !(element_dof.basis.nu < 1.0))
        {
          valid = false;
          break;
        }
        double &exponent = h1_exponents[element_dof.dof];
        if (std::isnan(exponent))
        {
          exponent = element_dof.basis.nu;
        }
        else if (exponent != element_dof.basis.nu)
        {
          valid = false;
          break;
        }
      }
      if (!valid)
      {
        break;
      }
    }
    valid = valid && std::all_of(h1_exponents.begin(), h1_exponents.end(),
                                 [](double exponent) { return std::isfinite(exponent); });
  }
  Mpi::GlobalAnd(1, &valid, fespace.GetComm());
  if (!valid)
  {
    throw std::invalid_argument(
        "Triangular singular H1 field evaluator received inconsistent space and DOF "
        "topology!");
  }
  local_enrichment.SetSize(static_cast<int>(topology.h1_dofs.size()));
  InitializeFaceNeighborEnrichment(fespace, topology.elements, true,
                                   face_neighbor_enrichment);
}

void TriangleEnrichedH1FieldEvaluator::SetFromTrueDofs(
    const mfem::Vector &combined_true_dofs)
{
  const int standard_size = fespace.GetTrueVSize();
  const int enrichment_size = static_cast<int>(numbering.h1.owned_size);
  bool valid = combined_true_dofs.Size() == standard_size + enrichment_size;
  Mpi::GlobalAnd(1, &valid, fespace.GetComm());
  if (!valid)
  {
    throw std::invalid_argument(
        "Combined triangular singular H1 vector has inconsistent process-local "
        "dimensions!");
  }

  mfem::Vector standard_true_dofs(standard_size);
  mfem::Vector enrichment_true_dofs(enrichment_size);
  for (int i = 0; i < standard_size; i++)
  {
    standard_true_dofs[i] = combined_true_dofs[i];
  }
  for (int i = 0; i < enrichment_size; i++)
  {
    enrichment_true_dofs[i] = combined_true_dofs[standard_size + i];
  }
  standard_field.SetFromTrueDofs(standard_true_dofs);
  standard_field.ExchangeFaceNbrData();
  enrichment_prolongation->Mult(enrichment_true_dofs, local_enrichment);
  ExchangeFaceNeighborEnrichment(topology.elements, local_enrichment, true, fespace,
                                 face_neighbor_enrichment);
  initialized = true;
}

TriangleH1FieldValue
TriangleEnrichedH1FieldEvaluator::Evaluate(int element, const mfem::IntegrationPoint &point)
{
  ValidateInteriorTrianglePoint(point);
  const TriangleBarycentricPoint lambda{1.0 - point.x - point.y, point.x, point.y};
  return EvaluateBarycentric(element, point, lambda);
}

TriangleH1FieldValue
TriangleEnrichedH1FieldEvaluator::EvaluateClosure(int element,
                                                  const mfem::IntegrationPoint &point)
{
  ValidateTriangleClosurePoint(point);
  const TriangleBarycentricPoint lambda{1.0 - point.x - point.y, point.x, point.y};
  return EvaluateBarycentric(element, point, lambda);
}

std::vector<TriangleVectorTraceTerm>
TriangleEnrichedH1FieldEvaluator::EvaluateGradientTraceExpansion(
    int element, const mfem::IntegrationPoint &point,
    const std::array<int, 2> &endpoint_nodes)
{
  ValidateTriangleClosurePoint(point);
  if (!initialized)
  {
    throw std::logic_error(
        "Triangular singular H1 trace evaluation requires a combined true-DOF vector!");
  }
  if (endpoint_nodes[0] < 0 || endpoint_nodes[0] >= 3 || endpoint_nodes[1] < 0 ||
      endpoint_nodes[1] >= 3 || endpoint_nodes[0] == endpoint_nodes[1])
  {
    throw std::invalid_argument(
        "Triangular singular H1 trace evaluation received an invalid element edge!");
  }
  const auto &element_dofs = GetParallelElementDofs(
      element, fespace.GetNE(), topology.elements, face_neighbor_enrichment);
  const auto &coefficients =
      element < fespace.GetNE() ? local_enrichment : face_neighbor_enrichment.coefficients;
  const TriangleBarycentricPoint lambda{1.0 - point.x - point.y, point.x, point.y};
  auto &transformation = GetParallelElementTransformation(fespace, element);
  double jacobian_determinant;
  const auto grad_lambda =
      GetTriangleBarycentricGradients(transformation, point, jacobian_determinant);

  transformation.SetIntPoint(&point);
  mfem::Vector standard_gradient(2);
  GetStandardGradient(standard_field, fespace, element, transformation, point,
                      standard_gradient);
  std::vector<TriangleVectorTraceTerm> terms{
      {{0.0, 0.0}, {standard_gradient[0], standard_gradient[1]}}};
  for (const auto &element_dof : element_dofs.h1)
  {
    const auto &basis = element_dof.basis;
    const double coefficient = coefficients[static_cast<int>(element_dof.dof)];
    const int singular_node = basis.nodes[0];
    int endpoint = -1;
    for (int e = 0; e < 2; e++)
    {
      endpoint = endpoint_nodes[e] == singular_node ? e : endpoint;
    }
    if (endpoint < 0)
    {
      const auto value = EvaluateTriangleNodeGradient(lambda, grad_lambda, basis.nodes[0],
                                                      basis.nodes[1], basis.nu);
      terms.push_back(
          {{0.0, 0.0}, {coefficient * value.value[0], coefficient * value.value[1]}});
      continue;
    }

    const int j = basis.nodes[1];
    const double rho = TriangleNodeRadialCoordinate(lambda, singular_node);
    if (!(rho > 0.0))
    {
      throw std::domain_error(
          "Triangular singular H1 trace expansion cannot be evaluated at its endpoint!");
    }
    const double ratio = lambda[j] / rho;
    terms.push_back(
        {{0.0, 0.0}, {coefficient * grad_lambda[j][0], coefficient * grad_lambda[j][1]}});
    TriangleTraceExponents exponents;
    (endpoint == 0 ? exponents.left : exponents.right) = basis.nu - 1.0;
    terms.push_back(
        {exponents,
         {coefficient * (-grad_lambda[j][0] +
                         ratio * (basis.nu - 1.0) * grad_lambda[singular_node][0]),
          coefficient * (-grad_lambda[j][1] +
                         ratio * (basis.nu - 1.0) * grad_lambda[singular_node][1])}});
  }
  return terms;
}

std::vector<TriangleScalarTraceTerm>
TriangleEnrichedH1FieldEvaluator::EvaluateValueTraceExpansion(
    int element, const mfem::IntegrationPoint &point,
    const std::array<int, 2> &endpoint_nodes)
{
  ValidateTriangleClosurePoint(point);
  if (!initialized)
  {
    throw std::logic_error(
        "Triangular singular H1 trace evaluation requires a combined true-DOF vector!");
  }
  if (endpoint_nodes[0] < 0 || endpoint_nodes[0] >= 3 || endpoint_nodes[1] < 0 ||
      endpoint_nodes[1] >= 3 || endpoint_nodes[0] == endpoint_nodes[1])
  {
    throw std::invalid_argument(
        "Triangular singular H1 trace evaluation received an invalid element edge!");
  }
  const auto &element_dofs = GetParallelElementDofs(
      element, fespace.GetNE(), topology.elements, face_neighbor_enrichment);
  const auto &coefficients =
      element < fespace.GetNE() ? local_enrichment : face_neighbor_enrichment.coefficients;
  const TriangleBarycentricPoint lambda{1.0 - point.x - point.y, point.x, point.y};
  auto &transformation = GetParallelElementTransformation(fespace, element);
  transformation.SetIntPoint(&point);
  std::vector<TriangleScalarTraceTerm> terms{
      {{0.0, 0.0}, standard_field.GetValue(transformation, point)}};
  for (const auto &element_dof : element_dofs.h1)
  {
    const auto &basis = element_dof.basis;
    const double coefficient = coefficients[static_cast<int>(element_dof.dof)];
    const int singular_node = basis.nodes[0];
    int endpoint = -1;
    for (int e = 0; e < 2; e++)
    {
      endpoint = endpoint_nodes[e] == singular_node ? e : endpoint;
    }
    if (endpoint < 0)
    {
      terms.push_back(
          {{0.0, 0.0},
           coefficient * EvaluateTriangleNodeGradientPotential(lambda, basis.nodes[0],
                                                               basis.nodes[1], basis.nu)});
      continue;
    }

    const int j = basis.nodes[1];
    const double rho = TriangleNodeRadialCoordinate(lambda, singular_node);
    if (!(rho > 0.0))
    {
      throw std::domain_error(
          "Triangular singular H1 trace expansion cannot be evaluated at its endpoint!");
    }
    const double ratio = lambda[j] / rho;
    terms.push_back({{0.0, 0.0}, coefficient * lambda[j]});
    TriangleTraceExponents exponents;
    (endpoint == 0 ? exponents.left : exponents.right) = basis.nu;
    terms.push_back({exponents, -coefficient * ratio});
  }
  return terms;
}

double TriangleEnrichedH1FieldEvaluator::GetElementNodeSingularExponent(int element,
                                                                        int node) const
{
  if (node < 0 || node >= 3)
  {
    throw std::out_of_range("Triangular singular H1 node-exponent query is out of range!");
  }
  const auto &element_dofs = GetParallelElementDofs(
      element, fespace.GetNE(), topology.elements, face_neighbor_enrichment);
  double exponent = 1.0;
  for (const auto &dof : element_dofs.h1)
  {
    if (dof.basis.nodes[0] == node)
    {
      exponent = std::min(exponent, dof.basis.nu);
    }
  }
  return exponent;
}

TriangleH1FieldValue TriangleEnrichedH1FieldEvaluator::EvaluateBarycentric(
    int element, const mfem::IntegrationPoint &point,
    const TriangleBarycentricPoint &lambda)
{
  if (!initialized)
  {
    throw std::logic_error(
        "Triangular singular H1 field evaluator has no combined true-DOF vector!");
  }
  const auto &element_dofs = GetParallelElementDofs(
      element, fespace.GetNE(), topology.elements, face_neighbor_enrichment);
  auto &transformation = GetParallelElementTransformation(fespace, element);
  double jacobian_determinant;
  const auto grad_lambda =
      GetTriangleBarycentricGradients(transformation, point, jacobian_determinant);

  transformation.SetIntPoint(&point);
  TriangleH1FieldValue result;
  result.potential = standard_field.GetValue(transformation, point);
  mfem::Vector standard_gradient(2);
  GetStandardGradient(standard_field, fespace, element, transformation, point,
                      standard_gradient);
  for (int d = 0; d < 2; d++)
  {
    result.gradient[d] = standard_gradient[d];
  }

  const auto &coefficients =
      element < fespace.GetNE() ? local_enrichment : face_neighbor_enrichment.coefficients;
  const auto enrichment =
      EvaluateElementTriangleH1Enrichment(element_dofs, coefficients, lambda, grad_lambda);
  result.potential += enrichment.potential;
  for (int d = 0; d < 2; d++)
  {
    result.gradient[d] += enrichment.gradient[d];
  }
  return result;
}

AdaptiveQuadratureResult TriangleEnrichedH1FieldEvaluator::IntegrateElementGradientEnergy(
    int element, double electric_coefficient, const AdaptiveAssemblyOptions &options)
{
  if (!initialized)
  {
    throw std::logic_error(
        "Triangular singular H1 energy evaluator has no combined true-DOF vector!");
  }
  if (!std::isfinite(electric_coefficient) || !(electric_coefficient > 0.0))
  {
    throw std::invalid_argument(
        "Triangular singular H1 energy integration requires a positive electric "
        "coefficient!");
  }
  if (options.quadrature_order < 1 || !std::isfinite(options.absolute_tolerance) ||
      options.absolute_tolerance < 0.0 || !std::isfinite(options.relative_tolerance) ||
      options.relative_tolerance < 0.0 ||
      !(options.absolute_tolerance > 0.0 || options.relative_tolerance > 0.0) ||
      options.maximum_subdivisions < 1)
  {
    throw std::invalid_argument(
        "Triangular singular H1 energy integration has invalid quadrature options!");
  }
  if (element < 0 || element >= fespace.GetNE())
  {
    throw std::out_of_range(
        "Triangular singular H1 energy integration element is out of range!");
  }
  auto *transformation = fespace.GetElementTransformation(element);
  if (!transformation)
  {
    throw std::runtime_error(
        "Triangular singular H1 energy element has no transformation!");
  }
  const auto *finite_element = fespace.GetFE(element);
  if (!finite_element || finite_element->GetGeomType() != mfem::Geometry::TRIANGLE ||
      finite_element->GetMapType() != mfem::FiniteElement::VALUE)
  {
    throw std::runtime_error(
        "Triangular singular H1 energy integration requires a value triangle element!");
  }
  if (topology.elements[element].h1.empty())
  {
    const int integration_order =
        2 * std::max(0, finite_element->GetOrder() - 1) + transformation->OrderW();
    const auto &rule = mfem::IntRules.Get(mfem::Geometry::TRIANGLE, integration_order);
    mfem::Vector gradient(2);
    double value = 0.0;
    for (int q = 0; q < rule.GetNPoints(); q++)
    {
      const auto &point = rule.IntPoint(q);
      transformation->SetIntPoint(&point);
      standard_field.GetGradient(*transformation, gradient);
      value += point.weight * electric_coefficient * transformation->Weight() *
               (gradient * gradient);
    }
    return {value, 0.0, 1, 0, true};
  }

  std::vector<int> singular_nodes;
  for (const auto &dof : topology.elements[element].h1)
  {
    const int node = dof.basis.nodes[0];
    if (node < 0 || node >= 3)
    {
      throw std::runtime_error(
          "Triangular singular H1 energy found an invalid singular node!");
    }
    singular_nodes.push_back(node);
  }
  std::sort(singular_nodes.begin(), singular_nodes.end());
  singular_nodes.erase(std::unique(singular_nodes.begin(), singular_nodes.end()),
                       singular_nodes.end());

  mfem::IntegrationPoint point;
  const auto integrand = [&](const TriangleBarycentricPoint &lambda)
  {
    point.Set2(lambda[1], lambda[2]);
    double jacobian_determinant;
    (void)GetTriangleBarycentricGradients(*transformation, point, jacobian_determinant);
    const auto value = EvaluateBarycentric(element, point, lambda);
    return electric_coefficient * jacobian_determinant *
           Dot(value.gradient, value.gradient);
  };
  const auto partition_weight =
      [&singular_nodes](const TriangleBarycentricPoint &lambda, int aligned_node)
  {
    if (singular_nodes.size() == 1)
    {
      return 1.0;
    }
    double minimum_rho = 1.0;
    for (int node : singular_nodes)
    {
      const double rho = TriangleNodeRadialCoordinate(lambda, node);
      if (!(rho > 0.0))
      {
        throw std::domain_error(
            "Partitioned triangle energy quadrature requires interior points!");
      }
      minimum_rho = std::min(minimum_rho, rho);
    }
    long double denominator = 0.0L;
    long double numerator = 0.0L;
    for (int node : singular_nodes)
    {
      const double rho = TriangleNodeRadialCoordinate(lambda, node);
      const long double weight =
          std::pow(minimum_rho / rho, MultiFeatureDuffyPartitionPower);
      denominator += weight;
      if (node == aligned_node)
      {
        numerator = weight;
      }
    }
    if (!(denominator > 0.0L) || !(numerator > 0.0L))
    {
      throw std::runtime_error(
          "Partitioned triangle energy quadrature produced an invalid weight!");
    }
    return static_cast<double>(numerator / denominator);
  };
  const auto integrate = [&](int order)
  {
    long double result = 0.0L;
    for (int node : singular_nodes)
    {
      result += IntegrateReferenceTriangleNodeDuffy(
          order, node, TriangleDuffyRadialPower, [&](const TriangleBarycentricPoint &lambda)
          { return partition_weight(lambda, node) * integrand(lambda); });
    }
    return static_cast<double>(result);
  };

  const int high_order = std::max(H1DuffyReferenceOrder, 2 * options.quadrature_order + 15);
  const int comparison_order = high_order - 8;
  const double value = integrate(high_order);
  const double comparison = integrate(comparison_order);
  const double scale = std::max({1.0, std::abs(value), std::abs(comparison)});
  const double error = H1DuffyErrorSafetyFactor * std::abs(value - comparison) +
                       64.0 * std::numeric_limits<double>::epsilon() * scale;
  const double tolerance =
      options.absolute_tolerance + options.relative_tolerance * std::abs(value);
  const bool converged = std::isfinite(value) && std::isfinite(comparison) &&
                         std::isfinite(error) && error >= 0.0 && error <= tolerance;
  const std::size_t charts = singular_nodes.size();
  const std::size_t point_count =
      charts * (ReferenceTriangleNodeDuffyQuadraturePointCount(high_order) +
                ReferenceTriangleNodeDuffyQuadraturePointCount(comparison_order));
  return {value, error, point_count, 0, converged};
}

void TriangleEnrichedH1FieldEvaluator::ProjectToDiscontinuousGridFunctions(
    mfem::ParGridFunction &potential, mfem::ParGridFunction &electric_field)
{
  const auto *potential_space = potential.ParFESpace();
  const auto *electric_space = electric_field.ParFESpace();
  if (!initialized || !potential_space || !electric_space ||
      potential_space->GetParMesh() != fespace.GetParMesh() ||
      electric_space->GetParMesh() != fespace.GetParMesh() ||
      potential_space->GetVDim() != 1 || electric_space->GetVDim() != 2 ||
      potential_space->FEColl()->GetMapType(2) != mfem::FiniteElement::VALUE ||
      electric_space->FEColl()->GetMapType(2) != mfem::FiniteElement::VALUE)
  {
    throw std::invalid_argument(
        "Triangular singular field sampling requires initialized scalar and vector value "
        "spaces on the solve mesh!");
  }

  TrianglePotentialCoefficient potential_coefficient(*this);
  TriangleElectricFieldCoefficient electric_coefficient(*this);
  potential.ProjectCoefficient(potential_coefficient);
  electric_field.ProjectCoefficient(electric_coefficient);
}

std::vector<H1CoefficientDiagnostic>
TriangleEnrichedH1FieldEvaluator::GetOwnedCoefficientDiagnostics() const
{
  if (!initialized)
  {
    throw std::logic_error(
        "Triangular singular H1 coefficient diagnostics require a combined true-DOF "
        "vector!");
  }

  std::vector<H1CoefficientDiagnostic> diagnostics;
  diagnostics.reserve(static_cast<std::size_t>(numbering.h1.owned_size));
  const int rank = Mpi::Rank(fespace.GetComm());
  for (std::size_t local = 0; local < topology.h1_dofs.size(); local++)
  {
    if (numbering.h1.owner[local] != rank)
    {
      continue;
    }
    const HYPRE_BigInt true_dof = numbering.h1.local_to_true[local];
    if (true_dof < numbering.h1.owned_offset ||
        true_dof >= numbering.h1.owned_offset + numbering.h1.owned_size ||
        !std::isfinite(local_enrichment[static_cast<int>(local)]) ||
        !std::isfinite(h1_exponents[local]))
    {
      throw std::runtime_error(
          "Triangular singular H1 coefficient diagnostics found inconsistent owned "
          "data!");
    }
    diagnostics.push_back({true_dof, topology.h1_dofs[local], h1_exponents[local],
                           local_enrichment[static_cast<int>(local)]});
  }
  std::sort(diagnostics.begin(), diagnostics.end(), [](const auto &left, const auto &right)
            { return left.true_dof < right.true_dof; });
  if (diagnostics.size() != static_cast<std::size_t>(numbering.h1.owned_size))
  {
    throw std::runtime_error(
        "Triangular singular H1 coefficient diagnostics do not cover every owned true "
        "DOF!");
  }
  return diagnostics;
}

std::vector<H1TipSlopeDiagnostic> TriangleEnrichedH1FieldEvaluator::FitTipSlopes(
    const TriangleFeatureTopology &features,
    const std::vector<GlobalVertexId> &source_vertex_ids,
    const std::vector<GlobalVertexId> &source_element_ids, const EdgeSlopeOptions &options)
{
  if (!initialized)
  {
    throw std::logic_error(
        "Triangular singular tip slope diagnostics require a combined true-DOF vector!");
  }
  return FitTriangleTipSlopes(fespace, features, source_vertex_ids, source_element_ids,
                              options,
                              [this](int element, const mfem::IntegrationPoint &point)
                              {
                                const auto field = Evaluate(element, point);
                                return std::sqrt(Dot(field.gradient, field.gradient));
                              });
}

TriangleEnrichedNDFieldEvaluator::TriangleEnrichedNDFieldEvaluator(
    const TriangleDofTopology &topology, const ParallelDofNumbering &numbering,
    mfem::ParFiniteElementSpace &fespace)
  : topology(topology), numbering(numbering), fespace(fespace), standard_field(&fespace),
    enrichment_prolongation(
        BuildParallelEnrichmentProlongation(fespace.GetComm(), numbering.nd)),
    nd_exponents(topology.nd_dofs.size(), std::numeric_limits<double>::quiet_NaN()),
    initialized(false)
{
  bool valid =
      fespace.GetParMesh() && fespace.GetParMesh()->Dimension() == 2 &&
      fespace.GetParMesh()->SpaceDimension() == 2 &&
      fespace.FEColl()->GetMapType(2) == mfem::FiniteElement::H_CURL &&
      topology.elements.size() == static_cast<std::size_t>(fespace.GetNE()) &&
      topology.nd_dofs.size() == numbering.nd.local_to_true.size() &&
      numbering.nd.local_size == static_cast<HYPRE_BigInt>(topology.nd_dofs.size()) &&
      numbering.nd.owned_size <=
          static_cast<HYPRE_BigInt>(std::numeric_limits<int>::max()) &&
      enrichment_prolongation &&
      enrichment_prolongation->Height() == static_cast<int>(topology.nd_dofs.size()) &&
      enrichment_prolongation->Width() == static_cast<int>(numbering.nd.owned_size);
  if (valid)
  {
    for (const auto &element : topology.elements)
    {
      for (const auto &element_dof : element.nd)
      {
        const auto family = element_dof.basis.family;
        if (element_dof.dof >= nd_exponents.size() || element_dof.basis.order != 1 ||
            (family != HigherOrderBasisFamily::NODE_GRADIENT &&
             family != HigherOrderBasisFamily::NODE_ROTATIONAL) ||
            topology.nd_dofs[element_dof.dof].family != family ||
            topology.nd_dofs[element_dof.dof].order != element_dof.basis.order ||
            !std::isfinite(element_dof.basis.nu) || !(element_dof.basis.nu > 0.0) ||
            !(element_dof.basis.nu < 1.0))
        {
          valid = false;
          break;
        }
        double &exponent = nd_exponents[element_dof.dof];
        if (std::isnan(exponent))
        {
          exponent = element_dof.basis.nu;
        }
        else if (exponent != element_dof.basis.nu)
        {
          valid = false;
          break;
        }
      }
      if (!valid)
      {
        break;
      }
    }
    valid = valid && std::all_of(nd_exponents.begin(), nd_exponents.end(),
                                 [](double exponent) { return std::isfinite(exponent); });
  }
  Mpi::GlobalAnd(1, &valid, fespace.GetComm());
  if (!valid)
  {
    throw std::invalid_argument(
        "Triangular singular ND field evaluator received inconsistent space and DOF "
        "topology!");
  }
  local_enrichment.SetSize(static_cast<int>(topology.nd_dofs.size()));
  InitializeFaceNeighborEnrichment(fespace, topology.elements, false,
                                   face_neighbor_enrichment);
}

void TriangleEnrichedNDFieldEvaluator::SetFromTrueDofs(
    const mfem::Vector &combined_true_dofs)
{
  const int standard_size = fespace.GetTrueVSize();
  const int enrichment_size = static_cast<int>(numbering.nd.owned_size);
  bool valid = combined_true_dofs.Size() == standard_size + enrichment_size;
  Mpi::GlobalAnd(1, &valid, fespace.GetComm());
  if (!valid)
  {
    throw std::invalid_argument(
        "Combined triangular singular ND vector has inconsistent process-local "
        "dimensions!");
  }

  mfem::Vector standard_true_dofs(standard_size);
  mfem::Vector enrichment_true_dofs(enrichment_size);
  for (int i = 0; i < standard_size; i++)
  {
    standard_true_dofs[i] = combined_true_dofs[i];
  }
  for (int i = 0; i < enrichment_size; i++)
  {
    enrichment_true_dofs[i] = combined_true_dofs[standard_size + i];
  }
  standard_field.SetFromTrueDofs(standard_true_dofs);
  standard_field.ExchangeFaceNbrData();
  enrichment_prolongation->Mult(enrichment_true_dofs, local_enrichment);
  ExchangeFaceNeighborEnrichment(topology.elements, local_enrichment, false, fespace,
                                 face_neighbor_enrichment);
  initialized = true;
}

TriangleNDFieldValue
TriangleEnrichedNDFieldEvaluator::Evaluate(int element, const mfem::IntegrationPoint &point)
{
  ValidateInteriorTrianglePoint(point);
  const TriangleBarycentricPoint lambda{1.0 - point.x - point.y, point.x, point.y};
  return EvaluateBarycentric(element, point, lambda);
}

TriangleNDFieldValue
TriangleEnrichedNDFieldEvaluator::EvaluateClosure(int element,
                                                  const mfem::IntegrationPoint &point)
{
  ValidateTriangleClosurePoint(point);
  const TriangleBarycentricPoint lambda{1.0 - point.x - point.y, point.x, point.y};
  return EvaluateBarycentric(element, point, lambda);
}

std::vector<TriangleVectorTraceTerm>
TriangleEnrichedNDFieldEvaluator::EvaluateTraceExpansion(
    int element, const mfem::IntegrationPoint &point,
    const std::array<int, 2> &endpoint_nodes)
{
  ValidateTriangleClosurePoint(point);
  if (!initialized)
  {
    throw std::logic_error(
        "Triangular singular ND trace evaluation requires a combined true-DOF vector!");
  }
  if (endpoint_nodes[0] < 0 || endpoint_nodes[0] >= 3 || endpoint_nodes[1] < 0 ||
      endpoint_nodes[1] >= 3 || endpoint_nodes[0] == endpoint_nodes[1])
  {
    throw std::invalid_argument(
        "Triangular singular ND trace evaluation received an invalid element edge!");
  }
  const auto &element_dofs = GetParallelElementDofs(
      element, fespace.GetNE(), topology.elements, face_neighbor_enrichment);
  const auto &coefficients =
      element < fespace.GetNE() ? local_enrichment : face_neighbor_enrichment.coefficients;
  const TriangleBarycentricPoint lambda{1.0 - point.x - point.y, point.x, point.y};
  auto &transformation = GetParallelElementTransformation(fespace, element);
  double jacobian_determinant;
  const auto grad_lambda =
      GetTriangleBarycentricGradients(transformation, point, jacobian_determinant);

  transformation.SetIntPoint(&point);
  mfem::Vector standard_value(2);
  standard_field.GetVectorValue(transformation, point, standard_value);
  std::vector<TriangleVectorTraceTerm> terms{
      {{0.0, 0.0}, {standard_value[0], standard_value[1]}}};
  for (const auto &element_dof : element_dofs.nd)
  {
    const auto &basis = element_dof.basis;
    const double coefficient = coefficients[static_cast<int>(element_dof.dof)];
    const int singular_node = basis.nodes[0];
    int endpoint = -1;
    for (int e = 0; e < 2; e++)
    {
      endpoint = endpoint_nodes[e] == singular_node ? e : endpoint;
    }
    if (endpoint < 0)
    {
      const auto value =
          basis.family == HigherOrderBasisFamily::NODE_GRADIENT
              ? EvaluateTriangleNodeGradient(lambda, grad_lambda, basis.nodes[0],
                                             basis.nodes[1], basis.nu)
              : EvaluateTriangleNodeRotational(lambda, grad_lambda, basis.nodes[0],
                                               basis.nodes[1], basis.nodes[2], basis.nu);
      terms.push_back(
          {{0.0, 0.0}, {coefficient * value.value[0], coefficient * value.value[1]}});
      continue;
    }

    if (basis.family == HigherOrderBasisFamily::NODE_GRADIENT)
    {
      const int j = basis.nodes[1];
      const double rho = TriangleNodeRadialCoordinate(lambda, singular_node);
      if (!(rho > 0.0))
      {
        throw std::domain_error(
            "Triangular singular ND trace expansion cannot be evaluated at its endpoint!");
      }
      const double ratio = lambda[j] / rho;
      terms.push_back(
          {{0.0, 0.0}, {coefficient * grad_lambda[j][0], coefficient * grad_lambda[j][1]}});
      TriangleTraceExponents exponents;
      (endpoint == 0 ? exponents.left : exponents.right) = basis.nu - 1.0;
      terms.push_back(
          {exponents,
           {coefficient * (-grad_lambda[j][0] +
                           ratio * (basis.nu - 1.0) * grad_lambda[singular_node][0]),
            coefficient * (-grad_lambda[j][1] +
                           ratio * (basis.nu - 1.0) * grad_lambda[singular_node][1])}});
    }
    else
    {
      const auto edge =
          EvaluateTriangleStandardEdge(lambda, grad_lambda, basis.nodes[1], basis.nodes[2]);
      terms.push_back(
          {{0.0, 0.0}, {-coefficient * edge.value[0], -coefficient * edge.value[1]}});
      TriangleTraceExponents exponents;
      (endpoint == 0 ? exponents.left : exponents.right) = basis.nu;
      terms.push_back(
          {exponents, {coefficient * edge.value[0], coefficient * edge.value[1]}});
    }
  }
  return terms;
}

double TriangleEnrichedNDFieldEvaluator::GetElementNodeSingularExponent(int element,
                                                                        int node) const
{
  if (node < 0 || node >= 3)
  {
    throw std::out_of_range("Triangular singular ND node-exponent query is out of range!");
  }
  const auto &element_dofs = GetParallelElementDofs(
      element, fespace.GetNE(), topology.elements, face_neighbor_enrichment);
  double exponent = 1.0;
  for (const auto &dof : element_dofs.nd)
  {
    if (dof.basis.nodes[0] == node)
    {
      exponent = std::min(exponent, dof.basis.nu);
    }
  }
  return exponent;
}

TriangleNDFieldValue TriangleEnrichedNDFieldEvaluator::EvaluateBarycentric(
    int element, const mfem::IntegrationPoint &point,
    const TriangleBarycentricPoint &lambda)
{
  if (!initialized)
  {
    throw std::logic_error(
        "Triangular singular ND field evaluator has no combined true-DOF vector!");
  }
  const auto &element_dofs = GetParallelElementDofs(
      element, fespace.GetNE(), topology.elements, face_neighbor_enrichment);
  auto &transformation = GetParallelElementTransformation(fespace, element);
  double jacobian_determinant;
  const auto grad_lambda =
      GetTriangleBarycentricGradients(transformation, point, jacobian_determinant);

  transformation.SetIntPoint(&point);
  TriangleNDFieldValue result;
  mfem::Vector standard_value(2), standard_curl(1);
  standard_field.GetVectorValue(transformation, point, standard_value);
  transformation.SetIntPoint(&point);
  GetStandardCurl(standard_field, fespace, element, transformation, point, standard_curl);
  for (int d = 0; d < 2; d++)
  {
    result.value[d] = standard_value[d];
  }
  result.curl = standard_curl[0];

  const auto &coefficients =
      element < fespace.GetNE() ? local_enrichment : face_neighbor_enrichment.coefficients;
  const auto enrichment =
      EvaluateElementTriangleNDEnrichment(element_dofs, coefficients, lambda, grad_lambda);
  for (int d = 0; d < 2; d++)
  {
    result.value[d] += enrichment.value[d];
  }
  result.curl += enrichment.curl;
  return result;
}

AdaptiveQuadratureResult TriangleEnrichedNDFieldEvaluator::IntegrateElementFieldEnergy(
    int element, double coefficient, const AdaptiveAssemblyOptions &options)
{
  if (!initialized)
  {
    throw std::logic_error(
        "Triangular singular ND energy evaluator has no combined true-DOF vector!");
  }
  if (!std::isfinite(coefficient) || !(coefficient > 0.0))
  {
    throw std::invalid_argument(
        "Triangular singular ND energy integration requires a positive coefficient!");
  }
  if (options.quadrature_order < 1 || !std::isfinite(options.absolute_tolerance) ||
      options.absolute_tolerance < 0.0 || !std::isfinite(options.relative_tolerance) ||
      options.relative_tolerance < 0.0 ||
      !(options.absolute_tolerance > 0.0 || options.relative_tolerance > 0.0) ||
      options.maximum_subdivisions < 1)
  {
    throw std::invalid_argument(
        "Triangular singular ND energy integration has invalid quadrature options!");
  }
  if (element < 0 || element >= fespace.GetNE())
  {
    throw std::out_of_range(
        "Triangular singular ND energy integration element is out of range!");
  }
  auto *transformation = fespace.GetElementTransformation(element);
  if (!transformation)
  {
    throw std::runtime_error(
        "Triangular singular ND energy element has no transformation!");
  }
  const auto *finite_element = fespace.GetFE(element);
  if (!finite_element || finite_element->GetGeomType() != mfem::Geometry::TRIANGLE ||
      finite_element->GetMapType() != mfem::FiniteElement::H_CURL)
  {
    throw std::runtime_error(
        "Triangular singular ND energy integration requires an H(curl) triangle "
        "element!");
  }
  if (topology.elements[element].nd.empty())
  {
    const int integration_order =
        2 * std::max(1, finite_element->GetOrder()) + transformation->OrderW();
    const auto &rule = mfem::IntRules.Get(mfem::Geometry::TRIANGLE, integration_order);
    double value = 0.0;
    for (int q = 0; q < rule.GetNPoints(); q++)
    {
      const auto &point = rule.IntPoint(q);
      const auto field = Evaluate(element, point);
      transformation->SetIntPoint(&point);
      value += point.weight * coefficient * transformation->Weight() *
               Dot(field.value, field.value);
    }
    return {value, 0.0, 1, 0, true};
  }

  std::vector<int> singular_nodes;
  for (const auto &dof : topology.elements[element].nd)
  {
    const int node = dof.basis.nodes[0];
    if (node < 0 || node >= 3)
    {
      throw std::runtime_error(
          "Triangular singular ND energy found an invalid singular node!");
    }
    singular_nodes.push_back(node);
  }
  std::sort(singular_nodes.begin(), singular_nodes.end());
  singular_nodes.erase(std::unique(singular_nodes.begin(), singular_nodes.end()),
                       singular_nodes.end());

  mfem::IntegrationPoint point;
  const auto integrand = [&](const TriangleBarycentricPoint &lambda)
  {
    point.Set2(lambda[1], lambda[2]);
    double jacobian_determinant;
    (void)GetTriangleBarycentricGradients(*transformation, point, jacobian_determinant);
    const auto value = EvaluateBarycentric(element, point, lambda);
    return coefficient * jacobian_determinant * Dot(value.value, value.value);
  };
  const auto partition_weight =
      [&singular_nodes](const TriangleBarycentricPoint &lambda, int aligned_node)
  {
    if (singular_nodes.size() == 1)
    {
      return 1.0;
    }
    double minimum_rho = 1.0;
    for (int node : singular_nodes)
    {
      const double rho = TriangleNodeRadialCoordinate(lambda, node);
      if (!(rho > 0.0))
      {
        throw std::domain_error(
            "Partitioned triangle ND energy quadrature requires interior points!");
      }
      minimum_rho = std::min(minimum_rho, rho);
    }
    long double denominator = 0.0L;
    long double numerator = 0.0L;
    for (int node : singular_nodes)
    {
      const double rho = TriangleNodeRadialCoordinate(lambda, node);
      const long double weight =
          std::pow(minimum_rho / rho, MultiFeatureDuffyPartitionPower);
      denominator += weight;
      if (node == aligned_node)
      {
        numerator = weight;
      }
    }
    if (!(denominator > 0.0L) || !(numerator > 0.0L))
    {
      throw std::runtime_error(
          "Partitioned triangle ND energy quadrature produced an invalid weight!");
    }
    return static_cast<double>(numerator / denominator);
  };
  const auto integrate = [&](int order)
  {
    long double result = 0.0L;
    for (int node : singular_nodes)
    {
      result += IntegrateReferenceTriangleNodeDuffy(
          order, node, TriangleDuffyRadialPower, [&](const TriangleBarycentricPoint &lambda)
          { return partition_weight(lambda, node) * integrand(lambda); });
    }
    return static_cast<double>(result);
  };

  const int high_order = std::max(H1DuffyReferenceOrder, 2 * options.quadrature_order + 15);
  const int comparison_order = high_order - 8;
  const double value = integrate(high_order);
  const double comparison = integrate(comparison_order);
  const double scale = std::max({1.0, std::abs(value), std::abs(comparison)});
  const double error = H1DuffyErrorSafetyFactor * std::abs(value - comparison) +
                       64.0 * std::numeric_limits<double>::epsilon() * scale;
  const double tolerance =
      options.absolute_tolerance + options.relative_tolerance * std::abs(value);
  const bool converged = std::isfinite(value) && std::isfinite(comparison) &&
                         std::isfinite(error) && error >= 0.0 && error <= tolerance;
  const std::size_t charts = singular_nodes.size();
  const std::size_t point_count =
      charts * (ReferenceTriangleNodeDuffyQuadraturePointCount(high_order) +
                ReferenceTriangleNodeDuffyQuadraturePointCount(comparison_order));
  return {value, error, point_count, 0, converged};
}

void TriangleEnrichedNDFieldEvaluator::ProjectToDiscontinuousGridFunctions(
    mfem::ParGridFunction &field, mfem::ParGridFunction &curl)
{
  const auto *field_space = field.ParFESpace();
  const auto *curl_space = curl.ParFESpace();
  if (!initialized || !field_space || !curl_space ||
      field_space->GetParMesh() != fespace.GetParMesh() ||
      curl_space->GetParMesh() != fespace.GetParMesh() || field_space->GetVDim() != 2 ||
      curl_space->GetVDim() != 1 ||
      field_space->FEColl()->GetMapType(2) != mfem::FiniteElement::VALUE ||
      curl_space->FEColl()->GetMapType(2) != mfem::FiniteElement::VALUE)
  {
    throw std::invalid_argument(
        "Triangular singular ND sampling requires initialized vector and scalar value "
        "spaces on the solve mesh!");
  }

  TriangleNDFieldCoefficient field_coefficient(*this);
  TriangleNDCurlCoefficient curl_coefficient(*this);
  field.ProjectCoefficient(field_coefficient);
  curl.ProjectCoefficient(curl_coefficient);
}

std::vector<H1CoefficientDiagnostic>
TriangleEnrichedNDFieldEvaluator::GetOwnedCoefficientDiagnostics() const
{
  if (!initialized)
  {
    throw std::logic_error(
        "Triangular singular ND coefficient diagnostics require a combined true-DOF "
        "vector!");
  }

  std::vector<H1CoefficientDiagnostic> diagnostics;
  diagnostics.reserve(static_cast<std::size_t>(numbering.nd.owned_size));
  const int rank = Mpi::Rank(fespace.GetComm());
  for (std::size_t local = 0; local < topology.nd_dofs.size(); local++)
  {
    if (numbering.nd.owner[local] != rank)
    {
      continue;
    }
    const HYPRE_BigInt true_dof = numbering.nd.local_to_true[local];
    if (true_dof < numbering.nd.owned_offset ||
        true_dof >= numbering.nd.owned_offset + numbering.nd.owned_size ||
        !std::isfinite(local_enrichment[static_cast<int>(local)]) ||
        !std::isfinite(nd_exponents[local]))
    {
      throw std::runtime_error(
          "Triangular singular ND coefficient diagnostics found inconsistent owned "
          "data!");
    }
    diagnostics.push_back({true_dof, topology.nd_dofs[local], nd_exponents[local],
                           local_enrichment[static_cast<int>(local)]});
  }
  std::sort(diagnostics.begin(), diagnostics.end(), [](const auto &left, const auto &right)
            { return left.true_dof < right.true_dof; });
  if (diagnostics.size() != static_cast<std::size_t>(numbering.nd.owned_size))
  {
    throw std::runtime_error(
        "Triangular singular ND coefficient diagnostics do not cover every owned true "
        "DOF!");
  }
  return diagnostics;
}

std::vector<H1TipSlopeDiagnostic> TriangleEnrichedNDFieldEvaluator::FitComplexTipSlopes(
    TriangleEnrichedNDFieldEvaluator &imaginary, const TriangleFeatureTopology &features,
    const std::vector<GlobalVertexId> &source_vertex_ids,
    const std::vector<GlobalVertexId> &source_element_ids, const EdgeSlopeOptions &options)
{
  if (!initialized || !imaginary.initialized)
  {
    throw std::logic_error(
        "Complex triangular singular tip slope diagnostics require initialized real "
        "and imaginary fields!");
  }
  if (&topology != &imaginary.topology || &numbering != &imaginary.numbering ||
      &fespace != &imaginary.fespace)
  {
    throw std::invalid_argument(
        "Complex triangular singular tip slope diagnostics require matching field "
        "evaluators!");
  }
  return FitTriangleTipSlopes(
      fespace, features, source_vertex_ids, source_element_ids, options,
      [this, &imaginary](int element, const mfem::IntegrationPoint &point)
      {
        const auto real = Evaluate(element, point);
        const auto imag = imaginary.Evaluate(element, point);
        return std::sqrt(Dot(real.value, real.value) + Dot(imag.value, imag.value));
      });
}

}  // namespace singular

}  // namespace fem

}  // namespace palace
