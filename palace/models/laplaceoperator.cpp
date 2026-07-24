// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "laplaceoperator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fem/multigrid.hpp"
#include "linalg/hypre.hpp"
#include "linalg/rap.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

namespace
{

class TracePotentialCoefficient : public mfem::Coefficient
{
private:
  struct Sample
  {
    std::array<double, 3> point;
    double value;
  };

  std::vector<Sample> curve;
  std::vector<std::array<Sample, 3>> triangles;
  int dimension;

  static double Dot(const std::array<double, 3> &a, const std::array<double, 3> &b)
  {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  }

  static std::array<double, 3> Subtract(const std::array<double, 3> &a,
                                        const std::array<double, 3> &b)
  {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
  }

  static std::pair<double, double> EvaluateTriangle(const std::array<Sample, 3> &triangle,
                                                    const std::array<double, 3> &point)
  {
    // Closest-point regions and barycentric coordinates from
    // "Real-Time Collision Detection", Christer Ericson, section 5.1.5.
    const auto ab = Subtract(triangle[1].point, triangle[0].point);
    const auto ac = Subtract(triangle[2].point, triangle[0].point);
    const auto ap = Subtract(point, triangle[0].point);
    const double d1 = Dot(ab, ap);
    const double d2 = Dot(ac, ap);
    double u, v, w;
    if (d1 <= 0.0 && d2 <= 0.0)
    {
      u = 1.0;
      v = w = 0.0;
    }
    else
    {
      const auto bp = Subtract(point, triangle[1].point);
      const double d3 = Dot(ab, bp);
      const double d4 = Dot(ac, bp);
      if (d3 >= 0.0 && d4 <= d3)
      {
        v = 1.0;
        u = w = 0.0;
      }
      else
      {
        const double vc = d1 * d4 - d3 * d2;
        if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
        {
          v = d1 / (d1 - d3);
          u = 1.0 - v;
          w = 0.0;
        }
        else
        {
          const auto cp = Subtract(point, triangle[2].point);
          const double d5 = Dot(ab, cp);
          const double d6 = Dot(ac, cp);
          if (d6 >= 0.0 && d5 <= d6)
          {
            w = 1.0;
            u = v = 0.0;
          }
          else
          {
            const double vb = d5 * d2 - d1 * d6;
            if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
            {
              w = d2 / (d2 - d6);
              u = 1.0 - w;
              v = 0.0;
            }
            else
            {
              const double va = d3 * d6 - d5 * d4;
              if (va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0)
              {
                w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
                u = 0.0;
                v = 1.0 - w;
              }
              else
              {
                const double denominator = 1.0 / (va + vb + vc);
                v = vb * denominator;
                w = vc * denominator;
                u = 1.0 - v - w;
              }
            }
          }
        }
      }
    }

    std::array<double, 3> closest;
    for (int d = 0; d < 3; d++)
    {
      closest[d] =
          u * triangle[0].point[d] + v * triangle[1].point[d] + w * triangle[2].point[d];
    }
    const auto delta = Subtract(point, closest);
    const double value =
        u * triangle[0].value + v * triangle[1].value + w * triangle[2].value;
    return {Dot(delta, delta), value};
  }

public:
  TracePotentialCoefficient(const std::string &path, int dimension_,
                            double mesh_coordinate_scale, double voltage_scale)
    : dimension(dimension_)
  {
    std::ifstream input(path);
    MFEM_VERIFY(input,
                "Unable to open prescribed potential trace file \"" << path << "\"!");

    std::string line;
    bool have_data = false;
    int columns = 0;
    int line_number = 0;
    std::map<int, std::vector<Sample>> triangle_samples;
    while (std::getline(input, line))
    {
      line_number++;
      auto first = line.find_first_not_of(" \t\r");
      if (first == std::string::npos || line[first] == '#')
      {
        continue;
      }
      std::replace(line.begin(), line.end(), ',', ' ');
      std::istringstream row(line);
      std::vector<double> data;
      double entry;
      while (row >> entry)
      {
        data.push_back(entry);
      }
      if (data.size() != 4 && data.size() != 5)
      {
        MFEM_VERIFY(!have_data, "Could not parse prescribed potential trace file \""
                                    << path << "\" at line " << line_number << "!");
        continue;  // Optional header before the first data row.
      }
      if (columns == 0)
      {
        columns = static_cast<int>(data.size());
      }
      MFEM_VERIFY(columns == static_cast<int>(data.size()),
                  "Mixed curve and surface rows in prescribed potential trace file \""
                      << path << "\" at line " << line_number << "!");
      Sample sample{{data[0], data[1], data[2]}, data[3]};
      MFEM_VERIFY(std::isfinite(sample.point[0]) && std::isfinite(sample.point[1]) &&
                      std::isfinite(sample.point[2]) && std::isfinite(sample.value),
                  "Non-finite value in prescribed potential trace file \""
                      << path << "\" at line " << line_number << "!");
      for (auto &x : sample.point)
      {
        x /= mesh_coordinate_scale;
      }
      sample.value /= voltage_scale;
      if (columns == 4)
      {
        curve.push_back(sample);
      }
      else
      {
        const int triangle = static_cast<int>(std::llround(data[4]));
        MFEM_VERIFY(triangle > 0 && static_cast<double>(triangle) == data[4],
                    "Invalid triangle index in prescribed potential surface file \""
                        << path << "\" at line " << line_number << "!");
        triangle_samples[triangle].push_back(sample);
      }
      have_data = true;
    }
    MFEM_VERIFY(have_data,
                "Prescribed potential trace file \"" << path << "\" contains no data!");
    if (columns == 4)
    {
      MFEM_VERIFY(curve.size() >= 3, "Prescribed potential curve file \""
                                         << path << "\" must contain at least 3 points!");
    }
    else
    {
      MFEM_VERIFY(dimension == 3, "Prescribed potential surface file \""
                                      << path << "\" requires a three-dimensional mesh!");
      triangles.reserve(triangle_samples.size());
      for (const auto &[index, samples] : triangle_samples)
      {
        MFEM_VERIFY(samples.size() == 3,
                    "Triangle " << index << " in prescribed potential surface \"" << path
                                << "\" must contain exactly 3 rows!");
        std::array<Sample, 3> triangle{samples[0], samples[1], samples[2]};
        const auto ab = Subtract(triangle[1].point, triangle[0].point);
        const auto ac = Subtract(triangle[2].point, triangle[0].point);
        const std::array<double, 3> cross = {ab[1] * ac[2] - ab[2] * ac[1],
                                             ab[2] * ac[0] - ab[0] * ac[2],
                                             ab[0] * ac[1] - ab[1] * ac[0]};
        MFEM_VERIFY(Dot(cross, cross) > 0.0,
                    "Degenerate triangle " << index << " in prescribed potential surface \""
                                           << path << "\"!");
        triangles.push_back(std::move(triangle));
      }
    }
  }

  double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override
  {
    double x_data[3] = {0.0, 0.0, 0.0};
    mfem::Vector x(x_data, dimension);
    T.Transform(ip, x);
    const std::array<double, 3> point{x_data[0], x_data[1], x_data[2]};

    double closest_distance_sq = std::numeric_limits<double>::infinity();
    double closest_value = 0.0;
    if (!triangles.empty())
    {
      for (const auto &triangle : triangles)
      {
        const auto [distance_sq, value] = EvaluateTriangle(triangle, point);
        if (distance_sq < closest_distance_sq)
        {
          closest_distance_sq = distance_sq;
          closest_value = value;
        }
      }
    }
    else
    {
      for (std::size_t i = 0; i < curve.size(); i++)
      {
        const std::size_t j = (i + 1) % curve.size();
        double length_sq = 0.0;
        double projection = 0.0;
        for (int d = 0; d < dimension; d++)
        {
          const double delta = curve[j].point[d] - curve[i].point[d];
          length_sq += delta * delta;
          projection += (x[d] - curve[i].point[d]) * delta;
        }
        if (length_sq <= 0.0)
        {
          continue;
        }
        const double t = std::clamp(projection / length_sq, 0.0, 1.0);
        double distance_sq = 0.0;
        for (int d = 0; d < dimension; d++)
        {
          const double delta =
              x[d] - (curve[i].point[d] + t * (curve[j].point[d] - curve[i].point[d]));
          distance_sq += delta * delta;
        }
        if (distance_sq < closest_distance_sq)
        {
          closest_distance_sq = distance_sq;
          closest_value = curve[i].value + t * (curve[j].value - curve[i].value);
        }
      }
    }
    MFEM_VERIFY(std::isfinite(closest_distance_sq),
                "Prescribed potential trace contains no valid interpolation elements!");
    return closest_value;
  }
};

}  // namespace

LaplaceOperator::LaplaceOperator(const config::BoundaryData &boundaries,
                                 const config::SolverData &solver,
                                 const std::vector<config::MaterialData> &materials,
                                 ProblemType problem_type,
                                 const std::vector<std::unique_ptr<Mesh>> &mesh)
  : print_hdr(true),
    dbc_attr(SetUpBoundaryProperties(boundaries.pec, boundaries.terminal,
                                     boundaries.prescribed_potential, *mesh.back())),
    h1_fecs(fem::ConstructFECollections<mfem::H1_FECollection>(
        solver.order, mesh.back()->Dimension(), solver.linear.mg_max_levels,
        solver.linear.mg_coarsening, false)),
    nd_fec(std::make_unique<mfem::ND_FECollection>(solver.order, mesh.back()->Dimension())),
    rt_fecs(fem::ConstructFECollections<mfem::RT_FECollection>(
        solver.order - 1, mesh.back()->Dimension(),
        solver.linear.estimator_mg ? solver.linear.mg_max_levels : 1,
        solver.linear.mg_coarsening, false)),
    h1_fespaces(fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
        solver.linear.mg_max_levels, mesh, h1_fecs, &dbc_attr, &dbc_tdof_lists)),
    nd_fespace(*mesh.back(), nd_fec.get()),
    rt_fespaces(fem::ConstructFiniteElementSpaceHierarchy<mfem::RT_FECollection>(
        solver.linear.estimator_mg ? solver.linear.mg_max_levels : 1, mesh, rt_fecs)),
    mat_op(materials, boundaries.periodic, problem_type, *mesh.back()),
    source_attr_lists(
        ConstructSources(boundaries.terminal, boundaries.prescribed_potential)),
    source_data_files(ConstructSourceDataFiles(boundaries.prescribed_potential)),
    source_terminal_attr_lists(
        ConstructSourceTerminalAttributes(boundaries.prescribed_potential)),
    mesh_coordinate_scale(1.0), voltage_scale(1.0)
{
  // Print essential BC information.
  if (dbc_attr.Size())
  {
    Mpi::Print("\nConfiguring Dirichlet BC at attributes:\n");
    utils::PrettyPrint(dbc_attr);
  }
}

LaplaceOperator::LaplaceOperator(const IoData &iodata,
                                 const std::vector<std::unique_ptr<Mesh>> &mesh)
  : LaplaceOperator(iodata.boundaries, iodata.solver, iodata.domains.materials,
                    iodata.problem.type, mesh)
{
  mesh_coordinate_scale = iodata.units.GetMeshLengthRelativeScale();
  voltage_scale = iodata.units.GetScaleFactor<Units::ValueType::VOLTAGE>();
}

mfem::Array<int> LaplaceOperator::SetUpBoundaryProperties(
    const config::PecBoundaryData &pec, const std::map<int, config::TerminalData> &terminal,
    const std::map<int, config::PrescribedPotentialData> &potential,
    const mfem::ParMesh &mesh)
{
  // Check that boundary attributes have been specified correctly.
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> bdr_attr_marker;
  if (!pec.empty() || !terminal.empty() || !potential.empty())
  {
    bdr_attr_marker.SetSize(bdr_attr_max);
    bdr_attr_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    std::set<int> bdr_warn_list;
    for (auto attr : pec.attributes)
    {
      if (attr <= 0 || attr > bdr_attr_max || !bdr_attr_marker[attr - 1])
      {
        bdr_warn_list.insert(attr);
      }
    }
    if (!bdr_warn_list.empty())
    {
      Mpi::Print("\n");
      Mpi::Warning("Unknown ground boundary attributes!\nSolver will just ignore them!");
      utils::PrettyPrint(bdr_warn_list, "Boundary attribute list:");
      Mpi::Print("\n");
    }
    for (const auto &[idx, data] : terminal)
    {
      for (auto attr : data.attributes)
      {
        MFEM_VERIFY(
            attr > 0 && attr <= bdr_attr_max,
            "Terminal boundary attribute tags must be non-negative and correspond to "
            "attributes in the mesh!");
        MFEM_VERIFY(bdr_attr_marker[attr - 1] > 0,
                    "Unknown terminal boundary attribute " << attr << "!");
      }
    }
    for (const auto &[idx, data] : potential)
    {
      for (auto attr : data.attributes)
      {
        MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
                    "Prescribed potential boundary attribute tags must be non-negative and "
                    "correspond to attributes in the mesh!");
        MFEM_VERIFY(bdr_attr_marker[attr - 1] > 0,
                    "Unknown prescribed potential boundary attribute " << attr << "!");
      }
      for (auto attr : data.terminal_attributes)
      {
        MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
                    "Prescribed potential terminal boundary attribute tags must be "
                    "non-negative and correspond to attributes in the mesh!");
        MFEM_VERIFY(bdr_attr_marker[attr - 1] > 0,
                    "Unknown prescribed potential terminal boundary attribute " << attr
                                                                                << "!");
      }
    }
  }

  // Mark selected boundary attributes from the mesh as essential (Dirichlet).
  mfem::Array<int> dbc_bcs;
  dbc_bcs.Reserve(static_cast<int>(pec.attributes.size()) +
                  static_cast<int>(terminal.size()));
  for (auto attr : pec.attributes)
  {
    if (attr <= 0 || attr > bdr_attr_max || !bdr_attr_marker[attr - 1])
    {
      continue;
    }
    dbc_bcs.Append(attr);
  }
  for (const auto &[idx, data] : terminal)
  {
    for (auto attr : data.attributes)
    {
      dbc_bcs.Append(attr);
    }
  }
  for (const auto &[idx, data] : potential)
  {
    for (auto attr : data.attributes)
    {
      dbc_bcs.Append(attr);
    }
    for (auto attr : data.terminal_attributes)
    {
      dbc_bcs.Append(attr);
    }
  }
  dbc_bcs.Sort();
  dbc_bcs.Unique();
  MFEM_VERIFY(dbc_bcs.Size() > 0,
              "Electrostatic problem is ill-posed without any Dirichlet boundaries!");
  return dbc_bcs;
}

std::map<int, mfem::Array<int>> LaplaceOperator::ConstructSources(
    const std::map<int, config::TerminalData> &terminal,
    const std::map<int, config::PrescribedPotentialData> &potential)
{
  // Construct mapping from terminal index to list of associated attributes.
  std::map<int, mfem::Array<int>> attr_lists;
  for (const auto &[idx, data] : terminal)
  {
    mfem::Array<int> &attr_list = attr_lists[idx];
    attr_list.Reserve(static_cast<int>(data.attributes.size()));
    for (auto attr : data.attributes)
    {
      attr_list.Append(attr);
    }
  }
  for (const auto &[idx, data] : potential)
  {
    mfem::Array<int> &attr_list = attr_lists[idx];
    attr_list.Reserve(
        static_cast<int>(data.attributes.size() + data.terminal_attributes.size()));
    for (auto attr : data.attributes)
    {
      attr_list.Append(attr);
    }
    for (auto attr : data.terminal_attributes)
    {
      attr_list.Append(attr);
    }
  }
  return attr_lists;
}

std::map<int, std::string> LaplaceOperator::ConstructSourceDataFiles(
    const std::map<int, config::PrescribedPotentialData> &potential)
{
  std::map<int, std::string> data_files;
  for (const auto &[idx, data] : potential)
  {
    data_files.emplace(idx, data.data_file);
  }
  return data_files;
}

std::map<int, mfem::Array<int>> LaplaceOperator::ConstructSourceTerminalAttributes(
    const std::map<int, config::PrescribedPotentialData> &potential)
{
  std::map<int, mfem::Array<int>> attributes;
  for (const auto &[idx, data] : potential)
  {
    if (!data.terminal_attributes.empty())
    {
      auto &list = attributes[idx];
      list.Reserve(static_cast<int>(data.terminal_attributes.size()));
      for (const int attribute : data.terminal_attributes)
      {
        list.Append(attribute);
      }
    }
  }
  return attributes;
}

namespace
{
void PrintHeader(const mfem::ParFiniteElementSpace &h1_fespace,
                 const mfem::ParFiniteElementSpace &nd_fespace,
                 const mfem::ParFiniteElementSpace &rt_fespace, bool &print_hdr)
{
  if (print_hdr)
  {
    Mpi::Print("\nAssembling system matrices, number of global unknowns:\n"
               " H1 (p = {:d}): {:d}, ND (p = {:d}): {:d}, RT (p = {:d}): {:d}\n Operator "
               "assembly level: {}\n",
               h1_fespace.GetMaxElementOrder(), h1_fespace.GlobalTrueVSize(),
               nd_fespace.GetMaxElementOrder(), nd_fespace.GlobalTrueVSize(),
               rt_fespace.GetMaxElementOrder(), rt_fespace.GlobalTrueVSize(),
               (h1_fespace.GetMaxElementOrder() >= BilinearForm::pa_order_threshold)
                   ? "Partial"
                   : "Full");

    const auto &mesh = *h1_fespace.GetParMesh();
    const auto geom_types = mesh::CheckElements(mesh).GetGeomTypes(mesh.Dimension());
    Mpi::Print(" Mesh geometries:\n");
    for (auto geom : geom_types)
    {
      const auto *fe = nd_fespace.FEColl()->FiniteElementForGeometry(geom);
      MFEM_VERIFY(fe, "MFEM does not support ND spaces on geometry = "
                          << mfem::Geometry::Name[geom] << "!");
      const int q_order = fem::DefaultIntegrationOrder::Get(mesh, geom);
      Mpi::Print("  {}: P = {:d}, Q = {:d} (quadrature order = {:d}){}\n",
                 mfem::Geometry::Name[geom], fe->GetDof(),
                 mfem::IntRules.Get(geom, q_order).GetNPoints(), q_order,
                 (geom == geom_types.back()) ? "" : ",");
    }

    Mpi::Print("\nAssembling multigrid hierarchy:\n");
  }
}

}  // namespace

std::unique_ptr<Operator> LaplaceOperator::GetStiffnessMatrix()
{
  // When partially assembled, the coarse operators can reuse the fine operator quadrature
  // data if the spaces correspond to the same mesh.
  PrintHeader(GetH1Space(), GetNDSpace(), GetRTSpace(), print_hdr);

  constexpr bool skip_zeros = false;
  MaterialPropertyCoefficient epsilon_func(mat_op.GetAttributeToMaterial(),
                                           mat_op.GetPermittivityReal());
  BilinearForm k(GetH1Space());
  k.AddDomainIntegrator<DiffusionIntegrator>(epsilon_func);
  // k.AssembleQuadratureData();
  auto k_vec = k.Assemble(GetH1Spaces(), skip_zeros);
  auto K = std::make_unique<MultigridOperator>(GetH1Spaces().GetNumLevels());
  for (std::size_t l = 0; l < GetH1Spaces().GetNumLevels(); l++)
  {
    const auto &h1_fespace_l = GetH1Spaces().GetFESpaceAtLevel(l);
    if (print_hdr)
    {
      Mpi::Print(" Level {:d} (p = {:d}): {:d} unknowns", l,
                 h1_fespace_l.GetMaxElementOrder(), h1_fespace_l.GlobalTrueVSize());
      if (const auto *k_spm = dynamic_cast<const hypre::HypreCSRMatrix *>(k_vec[l].get()))
      {
        HYPRE_BigInt nnz = k_spm->NNZ();
        Mpi::GlobalSum(1, &nnz, h1_fespace_l.GetComm());
        Mpi::Print(", {:d} NNZ\n", nnz);
      }
      else
      {
        Mpi::Print("\n");
      }
    }
    auto K_l = std::make_unique<ParOperator>(std::move(k_vec[l]), h1_fespace_l);
    K_l->SetEssentialTrueDofs(dbc_tdof_lists[l], Operator::DiagonalPolicy::DIAG_ONE);
    K->AddOperator(std::move(K_l));
  }

  print_hdr = false;
  return K;
}

void LaplaceOperator::GetExcitationVector(int idx, const Operator &K, Vector &X,
                                          Vector &RHS)
{
  // Apply the Dirichlet BCs to the solution vector: V = 1 on terminal boundaries with the
  // given index, V = 0 on all ground and other terminal boundaries.
  mfem::ParGridFunction x(&GetH1Space().Get());
  x = 0.0;

  // Get a marker of all boundary attributes with the given source surface index.
  const mfem::ParMesh &mesh = GetMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> source_marker = mesh::AttrToMarker(bdr_attr_max, source_attr_lists[idx]);
  if (auto it = source_data_files.find(idx); it != source_data_files.end())
  {
    TracePotentialCoefficient trace(it->second, mesh.SpaceDimension(),
                                    mesh_coordinate_scale, voltage_scale);
    x.ProjectBdrCoefficient(trace, source_marker);  // Values only correct on master
    if (auto terminal = source_terminal_attr_lists.find(idx);
        terminal != source_terminal_attr_lists.end())
    {
      auto terminal_marker = mesh::AttrToMarker(bdr_attr_max, terminal->second);
      mfem::ConstantCoefficient one_volt(1.0 / voltage_scale);
      x.ProjectBdrCoefficient(one_volt,
                              terminal_marker);  // Values only correct on master
    }
  }
  else
  {
    mfem::ConstantCoefficient one(1.0);
    x.ProjectBdrCoefficient(one, source_marker);  // Values only correct on master
  }

  // Eliminate the essential BC to get the RHS vector.
  X.SetSize(GetH1Space().GetTrueVSize());
  RHS.SetSize(GetH1Space().GetTrueVSize());
  X.UseDevice(true);
  RHS.UseDevice(true);
  X = 0.0;
  RHS = 0.0;
  x.ParallelProject(X);  // Restrict to the true dofs
  const auto *mg_K = dynamic_cast<const MultigridOperator *>(&K);
  const auto *PtAP_K = mg_K ? dynamic_cast<const ParOperator *>(&mg_K->GetFinestOperator())
                            : dynamic_cast<const ParOperator *>(&K);
  MFEM_VERIFY(PtAP_K, "LaplaceOperator requires ParOperator for RHS elimination!");
  PtAP_K->EliminateRHS(X, RHS);
}

}  // namespace palace
