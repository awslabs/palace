// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "sensitivity.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
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
  std::ifstream fi(path);
  MFEM_VERIFY(fi.is_open(), "Unable to open the design velocity table " << path << "!");
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

// Points, found by their coordinates to within a tolerance.
class PointHash
{
private:
  struct KeyHash
  {
    std::size_t operator()(const std::array<std::int64_t, 3> &key) const
    {
      std::size_t h = 0;
      for (const auto k : key)
      {
        h = h * 1000003u ^ std::hash<std::int64_t>{}(k);
      }
      return h;
    }
  };

  const std::vector<Point> &points;
  double tol;
  std::unordered_map<std::array<std::int64_t, 3>, std::vector<int>, KeyHash> cells;

  std::array<std::int64_t, 3> Key(const Point &x) const
  {
    return {static_cast<std::int64_t>(std::floor(x[0] / tol)),
            static_cast<std::int64_t>(std::floor(x[1] / tol)),
            static_cast<std::int64_t>(std::floor(x[2] / tol))};
  }

public:
  PointHash(const std::vector<Point> &points, double tol) : points(points), tol(tol)
  {
    for (std::size_t k = 0; k < points.size(); k++)
    {
      cells[Key(points[k])].push_back(static_cast<int>(k));
    }
  }

  // The indices of the points within the tolerance of x in each coordinate.
  std::vector<int> Find(const Point &x) const
  {
    std::vector<int> found;
    const auto key = Key(x);
    for (std::int64_t i = -1; i <= 1; i++)
    {
      for (std::int64_t j = -1; j <= 1; j++)
      {
        for (std::int64_t l = -1; l <= 1; l++)
        {
          const auto it = cells.find({key[0] + i, key[1] + j, key[2] + l});
          if (it == cells.end())
          {
            continue;
          }
          for (const int k : it->second)
          {
            if (std::abs(points[k][0] - x[0]) <= tol &&
                std::abs(points[k][1] - x[1]) <= tol &&
                std::abs(points[k][2] - x[2]) <= tol)
            {
              found.push_back(k);
            }
          }
        }
      }
    }
    return found;
  }
};

// The elements of a serial mesh containing points, from a grid of their bounding boxes.
class ElementLocator
{
private:
  mfem::Mesh &mesh;
  int dim;
  Point lo, h;
  std::array<int, 3> n;
  std::vector<std::vector<int>> cells;

  std::array<int, 3> Cell(const double *x) const
  {
    std::array<int, 3> cell = {0, 0, 0};
    for (int c = 0; c < dim; c++)
    {
      cell[c] =
          std::clamp(static_cast<int>(std::floor((x[c] - lo[c]) / h[c])), 0, n[c] - 1);
    }
    return cell;
  }

  int Index(const std::array<int, 3> &cell) const
  {
    return (cell[2] * n[1] + cell[1]) * n[0] + cell[0];
  }

public:
  ElementLocator(mfem::Mesh &mesh) : mesh(mesh), dim(mesh.SpaceDimension())
  {
    // The bounding box of each element, from its nodes, with a margin for curved elements.
    const auto &nodes = *mesh.GetNodes();
    const auto &fespace = *nodes.FESpace();
    std::vector<std::array<Point, 2>> boxes(mesh.GetNE());
    mfem::Array<int> dofs;
    lo.fill(std::numeric_limits<double>::max());
    Point hi;
    hi.fill(std::numeric_limits<double>::lowest());
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      fespace.GetElementDofs(e, dofs);
      auto &[bmin, bmax] = boxes[e];
      bmin.fill(std::numeric_limits<double>::max());
      bmax.fill(std::numeric_limits<double>::lowest());
      for (const int d : dofs)
      {
        for (int c = 0; c < dim; c++)
        {
          const double x = nodes(fespace.DofToVDof(d, c));
          bmin[c] = std::min(bmin[c], x);
          bmax[c] = std::max(bmax[c], x);
        }
      }
      double size = 0.0;
      for (int c = 0; c < dim; c++)
      {
        size = std::max(size, bmax[c] - bmin[c]);
      }
      for (int c = 0; c < dim; c++)
      {
        bmin[c] -= 0.1 * size;
        bmax[c] += 0.1 * size;
        lo[c] = std::min(lo[c], bmin[c]);
        hi[c] = std::max(hi[c], bmax[c]);
      }
    }

    // About one cell per element.
    const int n_cells = std::max(1, static_cast<int>(std::pow(mesh.GetNE(), 1.0 / dim)));
    n.fill(1);
    h.fill(1.0);
    for (int c = 0; c < dim; c++)
    {
      n[c] = n_cells;
      h[c] = std::max(hi[c] - lo[c], std::numeric_limits<double>::min()) / n[c];
    }
    cells.resize(n[0] * n[1] * n[2]);
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      const auto c0 = Cell(boxes[e][0].data()), c1 = Cell(boxes[e][1].data());
      for (int k = c0[2]; k <= c1[2]; k++)
      {
        for (int j = c0[1]; j <= c1[1]; j++)
        {
          for (int i = c0[0]; i <= c1[0]; i++)
          {
            cells[Index({i, j, k})].push_back(e);
          }
        }
      }
    }
  }

  // Find an element containing x, or the closest one within tol of it, and the reference
  // coordinates of x in it.
  bool Locate(const Point &x, double tol, int &elem, mfem::IntegrationPoint &ip) const
  {
    mfem::Vector xv(const_cast<double *>(x.data()), dim), y(dim);
    double closest = std::numeric_limits<double>::max();
    for (const int e : cells[Index(Cell(x.data()))])
    {
      auto &T = *mesh.GetElementTransformation(e);
      mfem::InverseElementTransformation inv(&T);
      inv.SetInitialGuessType(mfem::InverseElementTransformation::ClosestPhysNode);
      mfem::IntegrationPoint ipe;
      const int result = inv.Transform(xv, ipe);
      if (result == mfem::InverseElementTransformation::Inside)
      {
        elem = e;
        ip = ipe;
        return true;
      }
      if (result == mfem::InverseElementTransformation::Outside)
      {
        // The reference point is projected onto the element.
        T.Transform(ipe, y);
        y -= xv;
        if (y.Normlinf() < closest)
        {
          closest = y.Normlinf();
          elem = e;
          ip = ipe;
        }
      }
    }
    return closest <= tol;
  }
};

// The velocities [p][3 k + c] of the design parameters at the points x_k, in mesh length
// units: the input mesh's interpolant of the table, zero at the nodes it does not list.
std::vector<std::vector<double>> EvaluateVelocities(const IoData &iodata,
                                                    const VelocityTable &table,
                                                    const std::vector<double> &x)
{
  auto imesh = mesh::ReadInputMesh(iodata);
  imesh->EnsureNodes();
  const int sdim = imesh->SpaceDimension();
  const auto &nodes = *imesh->GetNodes();
  const auto &fespace = *nodes.FESpace();

  // The positions of the nodes of the input mesh, and a tolerance relative to its size.
  std::vector<Point> node_points(fespace.GetNDofs(), Point{0.0, 0.0, 0.0});
  Point lo, hi;
  lo.fill(std::numeric_limits<double>::max());
  hi.fill(std::numeric_limits<double>::lowest());
  for (int i = 0; i < fespace.GetNDofs(); i++)
  {
    for (int c = 0; c < sdim; c++)
    {
      node_points[i][c] = nodes(fespace.DofToVDof(i, c));
      lo[c] = std::min(lo[c], node_points[i][c]);
      hi[c] = std::max(hi[c], node_points[i][c]);
    }
  }
  double extent = 0.0;
  for (int c = 0; c < sdim; c++)
  {
    extent = std::max(extent, hi[c] - lo[c]);
  }
  const PointHash node_hash(node_points,
                            std::max(1.0e-9 * extent, std::numeric_limits<double>::min()));

  // The velocity fields on the input mesh, from the table.
  std::vector<mfem::GridFunction> V(
      table.names.size(),
      mfem::GridFunction(const_cast<mfem::FiniteElementSpace *>(&fespace)));
  for (auto &v : V)
  {
    v = 0.0;
  }
  std::size_t unmatched = 0;
  std::string first_unmatched;
  for (std::size_t k = 0; k < table.points.size(); k++)
  {
    const auto found = node_hash.Find(table.points[k]);
    if (found.empty())
    {
      if (unmatched++ == 0)
      {
        first_unmatched = fmt::format("({:.9e}, {:.9e}, {:.9e})", table.points[k][0],
                                      table.points[k][1], table.points[k][2]);
      }
      continue;
    }
    for (const int i : found)
    {
      for (std::size_t p = 0; p < V.size(); p++)
      {
        for (int c = 0; c < sdim; c++)
        {
          V[p](fespace.DofToVDof(i, c)) = table.velocities[p][k][c];
        }
      }
    }
  }
  MFEM_VERIFY(
      unmatched == 0,
      unmatched << " of the " << table.points.size()
                << " points of the design velocity table are not nodes of the input "
                   "mesh, the first at "
                << first_unmatched
                << ": the table must be at the nodes of the mesh file, in its "
                   "length units, with full precision!");

  // Evaluate them at each point: at a node of the input mesh, its value if the nodal basis
  // interpolates, and otherwise the interpolant in the element containing the point.
  const bool interpolatory =
      dynamic_cast<const mfem::NodalFiniteElement *>(fespace.GetFE(0)) != nullptr;
  const std::size_t n_points = x.size() / 3;
  std::vector<std::vector<double>> v(V.size(), std::vector<double>(x.size(), 0.0));
  const ElementLocator locator(*imesh);
  mfem::Vector value(sdim);
  for (std::size_t k = 0; k < n_points; k++)
  {
    const Point xk = {x[3 * k], x[3 * k + 1], x[3 * k + 2]};
    if (const auto found = interpolatory ? node_hash.Find(xk) : std::vector<int>{};
        !found.empty())
    {
      for (std::size_t p = 0; p < V.size(); p++)
      {
        for (int c = 0; c < sdim; c++)
        {
          v[p][3 * k + c] = V[p](fespace.DofToVDof(found[0], c));
        }
      }
      continue;
    }
    int e;
    mfem::IntegrationPoint ip;
    MFEM_VERIFY(locator.Locate(xk, 1.0e-9 * extent, e, ip),
                fmt::format("The mesh node at ({:.9e}, {:.9e}, {:.9e}) is not in the input "
                            "mesh, so the design velocities are not defined there!",
                            xk[0], xk[1], xk[2]));
    for (std::size_t p = 0; p < V.size(); p++)
    {
      V[p].GetVectorValue(e, ip, value);
      for (int c = 0; c < sdim; c++)
      {
        v[p][3 * k + c] = value(c);
      }
    }
  }
  return v;
}

}  // namespace

DesignVelocities::DesignVelocities(const IoData &iodata, const Mesh &mesh)
  : mesh(mesh.Get())
{
  const auto &path = iodata.sensitivity.velocities;
  if (path.empty())
  {
    return;
  }
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

  // Read the table and evaluate the velocities at the nodes on the root, then scatter
  // them one parameter at a time.
  const auto x_root = Mpi::GatherOnRoot(x, comm);
  VelocityTable table;
  std::vector<std::vector<double>> v_root;
  if (Mpi::Root(comm))
  {
    table = ReadVelocityTable(path);
    v_root = EvaluateVelocities(iodata, table, x_root);
  }
  int n_params = static_cast<int>(table.names.size());
  Mpi::Broadcast(1, &n_params, 0, comm);
  for (int p = 0; p < n_params; p++)
  {
    int length = Mpi::Root(comm) ? static_cast<int>(table.names[p].size()) : 0;
    Mpi::Broadcast(1, &length, 0, comm);
    std::string name = Mpi::Root(comm) ? table.names[p] : std::string(length, ' ');
    Mpi::Broadcast(length, name.data(), 0, comm);
    names.push_back(name);

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
  Mpi::Print(
      comm,
      "\nDesign velocities of {:d} parameter{} from {}, evaluated at {:d} mesh nodes\n",
      names.size(), (names.size() > 1) ? "s" : "", path, n_nodes);
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
  if (!Mpi::Root(comm))
  {
    return;
  }
  const auto path = post_dir / "sensitivity.json";
  std::ofstream fo(path);
  MFEM_VERIFY(fo.is_open(), "Unable to open sensitivity file " << path << "!");
  fo << sensitivities.dump(2) << '\n';
  MFEM_VERIFY(fo.good(), "Failed to write sensitivity file " << path << "!");
}

}  // namespace palace
