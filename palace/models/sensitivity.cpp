// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "sensitivity.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <fmt/format.h>
#include <nlohmann/json.hpp>
#include "fem/mesh.hpp"
#include "linalg/vector.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/outputdir.hpp"
#include "utils/tablecsv.hpp"

namespace palace
{

namespace
{

using Point = std::array<double, 3>;

// The velocities of the design parameters at nodes of the input mesh, in mesh length units
// per unit of the parameter.
struct VelocityTable
{
  std::vector<std::string> names;
  std::vector<Point> points;

  // The velocity of parameter p at point k is velocities[p][k].
  std::vector<std::vector<Point>> velocities;
};

VelocityTable ReadVelocityTable(const std::string &path)
{
  mfem::named_ifgzstream fi(path);
  MFEM_VERIFY(fi.good(), "Unable to open the design velocity table " << path << "!");
  std::stringstream buffer;
  buffer << fi.rdbuf();
  Table table(buffer.str());
  const auto n_cols = table.n_cols();
  MFEM_VERIFY(n_cols >= 6 && n_cols % 3 == 0,
              "The design velocity table "
                  << path
                  << " must have the columns x, y, z, then <name>.x, <name>.y, <name>.z "
                     "for each design parameter!");
  const std::size_t n_rows = table[0].data.size();
  for (std::size_t j = 0; j < n_cols; j++)
  {
    MFEM_VERIFY(table[j].data.size() == n_rows,
                "Column \"" << table[j].header_text << "\" of the design velocity table "
                            << path << " is missing values!");
  }
  constexpr std::array<const char *, 3> axes = {"x", "y", "z"};
  for (int c = 0; c < 3; c++)
  {
    MFEM_VERIFY(table[c].header_text == axes[c],
                "Column " << c + 1 << " of the design velocity table " << path
                          << " must be \"" << axes[c] << "\", not \""
                          << table[c].header_text << "\"!");
  }
  VelocityTable data;
  for (std::size_t j = 3; j < n_cols; j += 3)
  {
    const auto &header = table[j].header_text;
    const auto name = header.substr(0, header.size() - 2);
    for (int c = 0; c < 3; c++)
    {
      const auto expected = fmt::format("{}.{}", name, axes[c]);
      MFEM_VERIFY(!name.empty() && table[j + c].header_text == expected,
                  "Column " << j + c + 1 << " of the design velocity table " << path
                            << " must be the " << axes[c]
                            << " component of the velocity of a design parameter, as \""
                            << (name.empty() ? "<name>" : name) << "." << axes[c]
                            << "\", not \"" << table[j + c].header_text << "\"!");
    }
    MFEM_VERIFY(std::ranges::find(data.names, name) == data.names.end(),
                "Design parameter \"" << name << "\" appears twice in the design velocity "
                                      << "table " << path << "!");
    data.names.push_back(name);
  }
  data.points.resize(n_rows);
  data.velocities.assign(data.names.size(), std::vector<Point>(n_rows));
  for (std::size_t k = 0; k < n_rows; k++)
  {
    for (int c = 0; c < 3; c++)
    {
      data.points[k][c] = table[c].data[k];
      for (std::size_t p = 0; p < data.names.size(); p++)
      {
        data.velocities[p][k][c] = table[3 * (p + 1) + c].data[k];
      }
    }
  }
  return data;
}

// The key of a cell of a grid. Cells with equal keys share their points.
std::uint64_t CellKey(const std::array<std::int64_t, 3> &cell)
{
  std::uint64_t key = 1469598103934665603u;
  for (const auto c : cell)
  {
    key = (key ^ static_cast<std::uint64_t>(c)) * 1099511628211u;
  }
  return key;
}

std::array<std::int64_t, 3> Cell(const Point &x, double h)
{
  return {static_cast<std::int64_t>(std::floor(x[0] / h)),
          static_cast<std::int64_t>(std::floor(x[1] / h)),
          static_cast<std::int64_t>(std::floor(x[2] / h))};
}

}  // namespace

InputVelocities::InputVelocities(const IoData &iodata, MPI_Comm comm)
{
#if !defined(MFEM_USE_GSLIB)
  MFEM_ABORT("\"Sensitivity\" \"Velocities\" requires building with GSLIB!");
#endif
  if (Mpi::Root(comm))
  {
    const auto &path = iodata.sensitivity.velocities;
    const auto table = ReadVelocityTable(path);
    names = table.names;
    input_mesh = mesh::ReadInputMesh(iodata);
    input_mesh->EnsureNodes();
    SetNodes();

    // The velocity fields, from the table.
    const auto &fespace = *input_mesh->GetNodes()->FESpace();
    const int sdim = input_mesh->SpaceDimension();
    fields.assign(names.size(),
                  mfem::GridFunction(const_cast<mfem::FiniteElementSpace *>(&fespace)));
    std::vector<double *> V;
    for (auto &field : fields)
    {
      V.push_back(field.HostWrite());
      std::fill_n(V.back(), field.Size(), 0.0);
    }
    std::size_t unmatched = 0;
    std::string first_unmatched;
    for (std::size_t k = 0; k < table.points.size(); k++)
    {
      const auto found = FindNodes(table.points[k]);
      if (found.empty() && unmatched++ == 0)
      {
        first_unmatched = fmt::format("({:.9e}, {:.9e}, {:.9e})", table.points[k][0],
                                      table.points[k][1], table.points[k][2]);
      }
      for (const int i : found)
      {
        for (std::size_t p = 0; p < V.size(); p++)
        {
          for (int c = 0; c < sdim; c++)
          {
            V[p][fespace.DofToVDof(i, c)] = table.velocities[p][k][c];
          }
        }
      }
    }
    MFEM_VERIFY(unmatched == 0,
                unmatched << " of the " << table.points.size() << " rows of " << path
                          << " are not nodes of the input mesh, the first at "
                          << first_unmatched
                          << ": the table must be at the nodes of the mesh file, in its "
                             "length units, with full precision!");

    // Without its curvature, the mesh moves with its vertices: interpolate their
    // velocities.
    if (iodata.model.remove_curvature && fespace.GetMaxElementOrder() > 1)
    {
      std::vector<std::vector<Point>> vertex(fields.size(),
                                             std::vector<Point>(input_mesh->GetNV()));
      mfem::Array<int> dofs;
      for (int v = 0; v < input_mesh->GetNV(); v++)
      {
        fespace.GetVertexDofs(v, dofs);
        for (std::size_t p = 0; p < fields.size(); p++)
        {
          for (int c = 0; c < sdim; c++)
          {
            vertex[p][v][c] = V[p][fespace.DofToVDof(dofs[0], c)];
          }
        }
      }
      input_mesh->SetCurvature(1);
      SetNodes();
      const auto &linear = *input_mesh->GetNodes()->FESpace();
      for (std::size_t p = 0; p < fields.size(); p++)
      {
        fields[p].SetSpace(const_cast<mfem::FiniteElementSpace *>(&linear));
        auto *h_V = fields[p].HostWrite();
        for (int v = 0; v < input_mesh->GetNV(); v++)
        {
          linear.GetVertexDofs(v, dofs);
          for (int c = 0; c < sdim; c++)
          {
            h_V[linear.DofToVDof(dofs[0], c)] = vertex[p][v][c];
          }
        }
      }
    }
#if defined(MFEM_USE_GSLIB)
    finder = std::make_unique<mfem::FindPointsGSLIB>(MPI_COMM_SELF);
    finder->Setup(*input_mesh, 0.01, 1.0e-12);
    finder->SetDefaultInterpolationValue(0.0);
    finder->SetL2AvgType(mfem::FindPointsGSLIB::NONE);
#endif
  }
  int n = static_cast<int>(names.size());
  Mpi::Broadcast(1, &n, 0, comm);
  names.resize(n);
  for (auto &name : names)
  {
    int length = static_cast<int>(name.size());
    Mpi::Broadcast(1, &length, 0, comm);
    name.resize(length);
    Mpi::Broadcast(length, name.data(), 0, comm);
  }
}

InputVelocities::~InputVelocities()
{
#if defined(MFEM_USE_GSLIB)
  if (finder)
  {
    finder->FreeData();
  }
#endif
}

void InputVelocities::SetNodes()
{
  // The positions of the nodes, by their cells in a grid of a spacing relative to the size
  // of the mesh.
  const auto &fespace = *input_mesh->GetNodes()->FESpace();
  const auto *X = input_mesh->GetNodes()->HostRead();
  const int sdim = input_mesh->SpaceDimension();
  nodes.assign(fespace.GetNDofs(), Point{0.0, 0.0, 0.0});
  Point lo, hi;
  lo.fill(std::numeric_limits<double>::max());
  hi.fill(std::numeric_limits<double>::lowest());
  for (int i = 0; i < fespace.GetNDofs(); i++)
  {
    for (int c = 0; c < sdim; c++)
    {
      nodes[i][c] = X[fespace.DofToVDof(i, c)];
      lo[c] = std::min(lo[c], nodes[i][c]);
      hi[c] = std::max(hi[c], nodes[i][c]);
    }
  }
  double extent = 0.0;
  for (int c = 0; c < sdim; c++)
  {
    extent = std::max(extent, hi[c] - lo[c]);
  }
  tol = std::max(1.0e-9 * extent, std::numeric_limits<double>::min());
  cells.clear();
  for (std::size_t i = 0; i < nodes.size(); i++)
  {
    cells[CellKey(Cell(nodes[i], tol))].push_back(static_cast<int>(i));
  }
}

std::vector<int> InputVelocities::FindNodes(const Point &x) const
{
  // The nodes within tol of x in each coordinate are in its cell or a neighbor.
  std::vector<int> found;
  const auto cell = Cell(x, tol);
  for (std::int64_t i = -1; i <= 1; i++)
  {
    for (std::int64_t j = -1; j <= 1; j++)
    {
      for (std::int64_t l = -1; l <= 1; l++)
      {
        const auto it = cells.find(CellKey({cell[0] + i, cell[1] + j, cell[2] + l}));
        if (it == cells.end())
        {
          continue;
        }
        for (const int k : it->second)
        {
          if (std::abs(nodes[k][0] - x[0]) <= tol && std::abs(nodes[k][1] - x[1]) <= tol &&
              std::abs(nodes[k][2] - x[2]) <= tol)
          {
            found.push_back(k);
          }
        }
      }
    }
  }
  std::ranges::sort(found);
  found.erase(std::ranges::unique(found).begin(), found.end());
  return found;
}

std::vector<std::vector<double>>
InputVelocities::Evaluate(const std::vector<double> &x) const
{
  // At a node of the input mesh, the velocity is its value if the nodal basis interpolates.
  // Elsewhere, it is the interpolant in the element containing the point.
  const auto &fespace = *fields[0].FESpace();
  const int sdim = input_mesh->SpaceDimension();
  const bool interpolatory =
      dynamic_cast<const mfem::NodalFiniteElement *>(fespace.GetFE(0)) != nullptr;
  std::vector<const double *> V;
  for (const auto &field : fields)
  {
    V.push_back(field.HostRead());
  }
  std::vector<std::vector<double>> v(fields.size(), std::vector<double>(x.size(), 0.0));
  std::vector<std::size_t> search;
  for (std::size_t k = 0; k < x.size() / 3; k++)
  {
    const auto found = interpolatory ? FindNodes({x[3 * k], x[3 * k + 1], x[3 * k + 2]})
                                     : std::vector<int>{};
    if (found.empty())
    {
      search.push_back(k);
      continue;
    }
    for (std::size_t p = 0; p < V.size(); p++)
    {
      for (int c = 0; c < sdim; c++)
      {
        v[p][3 * k + c] = V[p][fespace.DofToVDof(found[0], c)];
      }
    }
  }
  if (search.empty())
  {
    return v;
  }
#if defined(MFEM_USE_GSLIB)
  const int m = static_cast<int>(search.size());
  mfem::Vector xyz(m * sdim);
  for (int i = 0; i < m; i++)
  {
    for (int c = 0; c < sdim; c++)
    {
      xyz(c * m + i) = x[3 * search[i] + c];
    }
  }
  finder->FindPoints(xyz, mfem::Ordering::byNODES);
  const auto &code = finder->GetCode();
  const auto *dist = finder->GetDist().HostRead();
  for (int i = 0; i < m; i++)
  {
    MFEM_VERIFY(code[i] < 2 && dist[i] <= tol,
                fmt::format("The mesh node at ({:.9e}, {:.9e}, {:.9e}) is not in the input "
                            "mesh, so the design velocities are not defined there!",
                            x[3 * search[i]], x[3 * search[i] + 1], x[3 * search[i] + 2]));
  }
  mfem::Vector vals;
  for (std::size_t p = 0; p < fields.size(); p++)
  {
    finder->Interpolate(fields[p], vals, mfem::Ordering::byVDIM);
    const auto *h_vals = vals.HostRead();
    for (int i = 0; i < m; i++)
    {
      for (int c = 0; c < sdim; c++)
      {
        v[p][3 * search[i] + c] = h_vals[i * sdim + c];
      }
    }
  }
#else
  MFEM_ABORT("Design velocities require MFEM built with GSLIB!");
#endif
  return v;
}

DesignVelocities::DesignVelocities(const InputVelocities *input, const IoData &iodata,
                                   const Mesh &mesh)
  : mesh(mesh.Get())
{
  if (!input)
  {
    return;
  }
  names = input->GetNames();
  const auto comm = mesh.GetComm();
  const auto &nodes = *this->mesh.GetNodes();
  const auto *nodal_fespace =
      dynamic_cast<const mfem::ParFiniteElementSpace *>(nodes.FESpace());
  MFEM_VERIFY(nodal_fespace, "Design velocities require a parallel mesh with nodes!");
  const int sdim = this->mesh.SpaceDimension();

  // The owned true DOFs of the mesh nodes and their positions, in mesh length units.
  // Constrained (hanging) nodes are not true DOFs.
  const double scale = iodata.units.GetMeshLengthRelativeScale();
  std::vector<std::array<int, 3>> tdofs;
  std::vector<double> x;
  const auto *h_nodes = nodes.HostRead();
  for (int i = 0; i < nodal_fespace->GetNDofs(); i++)
  {
    std::array<int, 3> tdof = {-1, -1, -1};
    std::array<double, 3> xi = {0.0, 0.0, 0.0};
    for (int c = 0; c < sdim; c++)
    {
      const int ldof = nodal_fespace->DofToVDof(i, c);
      tdof[c] = nodal_fespace->GetLocalTDofNumber(ldof);
      xi[c] = h_nodes[ldof] * scale;
    }
    if (tdof[0] < 0)
    {
      continue;
    }
    tdofs.push_back(tdof);
    x.insert(x.end(), xi.begin(), xi.end());
  }

  // Evaluate the velocities at the nodes on the root, and scatter them one parameter at a
  // time.
  const auto x_root = Mpi::GatherOnRoot(x, comm);
  const auto v_root =
      Mpi::Root(comm) ? input->Evaluate(x_root) : std::vector<std::vector<double>>{};
  for (std::size_t p = 0; p < names.size(); p++)
  {
    const auto v = Mpi::ScatterFromRoot(Mpi::Root(comm) ? v_root[p] : std::vector<double>{},
                                        static_cast<int>(x.size()), comm);
    auto &V = velocities.emplace_back(nodal_fespace->GetTrueVSize());
    V.UseDevice(true);
    V = 0.0;
    auto *h_V = V.HostReadWrite();
    for (std::size_t k = 0; k < tdofs.size(); k++)
    {
      for (int c = 0; c < sdim; c++)
      {
        h_V[tdofs[k][c]] = v[3 * k + c] / scale;
      }
    }
  }
  long long int n_nodes = tdofs.size();
  Mpi::GlobalSum(1, &n_nodes, comm);
  Mpi::Print(comm, "\nDesign velocities of {:d} parameter{} at {:d} mesh nodes\n",
             names.size(), (names.size() > 1) ? "s" : "", n_nodes);
}

DesignVelocities::~DesignVelocities() = default;

std::vector<double> DesignVelocities::Contract(const Vector &g) const
{
  std::vector<double> d(Size());
  for (std::size_t p = 0; p < Size(); p++)
  {
    d[p] = linalg::LocalDot(g, velocities[p]);
  }
  Mpi::GlobalSum(static_cast<int>(d.size()), d.data(), mesh.GetComm());
  return d;
}

std::vector<std::vector<double>> DesignVelocities::Differentiate(
    const std::function<std::vector<double>(const mfem::ParMesh &)> &f) const
{
  std::vector<std::vector<double>> d(Size());
  if (Size() == 0)
  {
    return d;
  }
  if (!moved_mesh)
  {
    moved_mesh = std::make_unique<mfem::ParMesh>(mesh);
  }
  const auto comm = mesh.GetComm();
  const auto &X = *mesh.GetNodes();
  auto &nodes = *moved_mesh->GetNodes();
  const auto &nodal_fespace = *X.FESpace();
  const int sdim = mesh.SpaceDimension();

  // Steps which move the nodes by a small fraction of the size of the mesh.
  std::array<double, 3> lo, hi;
  lo.fill(std::numeric_limits<double>::max());
  hi.fill(std::numeric_limits<double>::lowest());
  const auto *h_X = X.HostRead();
  for (int i = 0; i < nodal_fespace.GetNDofs(); i++)
  {
    for (int c = 0; c < sdim; c++)
    {
      const double x = h_X[nodal_fespace.DofToVDof(i, c)];
      lo[c] = std::min(lo[c], x);
      hi[c] = std::max(hi[c], x);
    }
  }
  Mpi::GlobalMin(sdim, lo.data(), comm);
  Mpi::GlobalMax(sdim, hi.data(), comm);
  double extent = 0.0;
  for (int c = 0; c < sdim; c++)
  {
    extent = std::max(extent, hi[c] - lo[c]);
  }

  mfem::Vector V_l(nodal_fespace.GetVSize());
  for (std::size_t p = 0; p < Size(); p++)
  {
    double v_max = velocities[p].Normlinf();
    Mpi::GlobalMax(1, &v_max, comm);
    const double h = (v_max > 0.0) ? 1.0e-6 * extent / v_max : 1.0;
    nodal_fespace.GetProlongationMatrix()->Mult(velocities[p], V_l);
    auto Evaluate = [&](double t)
    {
      nodes = X;
      nodes.Add(t, V_l);
      moved_mesh->NodesUpdated();
      return f(*moved_mesh);
    };
    const auto fp = Evaluate(h), fm = Evaluate(-h);
    MFEM_ASSERT(fp.size() == fm.size(), "Geometric quantities changed in number!");
    d[p].resize(fp.size());
    for (std::size_t k = 0; k < fp.size(); k++)
    {
      d[p][k] = (fp[k] - fm[k]) / (2.0 * h);
    }
  }
  nodes = X;
  moved_mesh->NodesUpdated();
  return d;
}

mfem::DenseMatrix GramMatrix(const mfem::Operator &op, const std::vector<Vector> &u,
                             MPI_Comm comm)
{
  const int n = static_cast<int>(u.size());
  mfem::DenseMatrix G(n);
  G = 0.0;
  Vector y(op.Height());
  y.UseDevice(true);
  for (int i = 0; i < n; i++)
  {
    op.Mult(u[i], y);
    for (int j = i; j < n; j++)
    {
      G(i, j) = linalg::LocalDot(u[j], y);
    }
  }
  Mpi::GlobalSum(n * n, G.Data(), comm);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < i; j++)
    {
      G(i, j) = G(j, i);
    }
  }
  return G;
}

nlohmann::json ToJson(const mfem::DenseMatrix &A, double scale)
{
  auto rows = nlohmann::json::array();
  for (int i = 0; i < A.Height(); i++)
  {
    auto &row = rows.emplace_back(nlohmann::json::array());
    for (int j = 0; j < A.Width(); j++)
    {
      row.push_back(A(i, j) * scale);
    }
  }
  return rows;
}

std::string MaterialDerivativeName(int m, const std::string &property)
{
  return fmt::format("Materials[{:d}].{}", m + 1, property);
}

std::string LumpedPortDerivativeName(int idx, const std::string &element)
{
  return fmt::format("LumpedPorts[{:d}].{}", idx, element);
}

void WriteSensitivities(const nlohmann::json &sensitivities,
                        const std::filesystem::path &post_dir, MPI_Comm comm)
{
  WriteRootOutputFile(post_dir / "sensitivity.json", comm,
                      [&](std::ostream &os) { os << sensitivities.dump(2) << '\n'; });
}

}  // namespace palace
