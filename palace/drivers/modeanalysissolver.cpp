// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "modeanalysissolver.hpp"

#include <algorithm>
#include <array>
#include <complex>
#include <map>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "fem/errorindicator.hpp"
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/vector.hpp"
#include "models/boundarymodeoperator.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/modeanalysisoperator.hpp"
#include "models/postoperator.hpp"
#include "models/surfaceconductivityoperator.hpp"
#include "models/surfaceimpedanceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

std::pair<ErrorIndicator, long long int>
ModeAnalysisSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  const auto &ma_data = iodata.solver.mode_analysis;
  const bool use_submesh = !ma_data.attributes.empty();

  if (use_submesh)
  {
    MFEM_VERIFY(mesh.back()->Dimension() == 3,
                "ModeAnalysis with \"Attributes\" requires a 3D mesh!");
  }
  else
  {
    MFEM_VERIFY(mesh.back()->Dimension() == 2,
                "ModeAnalysis solver requires a 2D mesh (waveguide cross-section), "
                "or a 3D mesh with \"Attributes\" specifying the cross-section boundary!");
  }

  const double freq_GHz = ma_data.freq;
  const int num_modes = ma_data.n;
  const double tol = ma_data.tol;
  const double omega =
      2.0 * M_PI * iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_GHz);

  Mpi::Print("\nConfiguring 2D waveguide mode analysis at f = {:.3e} GHz "
             "(omega = {:.6e})\n",
             freq_GHz, omega);

  BlockTimer bt0(Timer::CONSTRUCT);

  std::unique_ptr<Mesh> submesh_holder;
  std::unique_ptr<MaterialOperator> owned_mat_op;
  mfem::Vector submesh_centroid, submesh_e1, submesh_e2;
  std::vector<std::array<int, 3>> pec_internal_edges;  // {v1, v2, attr}
  Mesh *solve_mesh;

  if (use_submesh)
  {
    Mpi::Print(" Extracting 2D submesh from 3D boundary attributes...\n");
    const auto &parent_mesh = mesh.back()->Get();
    MPI_Comm comm = parent_mesh.GetComm();
    mfem::Array<int> attr_list;
    attr_list.Append(ma_data.attributes.data(), ma_data.attributes.size());

    // Step 1: Extract ParSubMesh from 3D boundary.
    auto par_submesh = mfem::ParSubMesh::CreateFromBoundary(parent_mesh, attr_list);

    // Step 2: Compute surface normal (MPI collective — all ranks participate).
    mfem::Vector surface_normal = mesh::GetSurfaceNormal(par_submesh);
    Mpi::Print(" Surface normal = ({:+.3e}, {:+.3e}, {:+.3e})\n", surface_normal(0),
               surface_normal(1), surface_normal(2));

    // Step 3: Remap domain attributes from neighboring 3D elements.
    mesh::RemapSubMeshAttributes(par_submesh);

    // Step 4: Remap boundary attributes from parent boundary faces.
    mesh::RemapSubMeshBdrAttributes(par_submesh, attr_list);

    // Step 5a: Collect PEC internal edges from the parent mesh that exist on the submesh.
    // PEC boundary faces in the 3D parent are distributed across all MPI ranks. Each rank
    // collects its local PEC parent edge indices, gathers them to rank 0 (which has the
    // submesh), and rank 0 maps them to submesh vertex pairs for boundary element creation.
    {
      // Collect PEC-related boundary attributes from config.
      std::vector<int> pec_attrs;
      for (auto a : iodata.boundaries.pec.attributes)
      {
        pec_attrs.push_back(a);
      }
      for (auto a : iodata.boundaries.auxpec.attributes)
      {
        pec_attrs.push_back(a);
      }
      for (const auto &d : iodata.boundaries.conductivity)
      {
        for (auto a : d.attributes)
        {
          pec_attrs.push_back(a);
        }
      }

      // Each rank collects PEC parent edge indices and attributes from its local
      // parent boundary elements. Use GLOBAL edge numbering for cross-rank matching.
      // ParMesh doesn't expose global edge numbers directly, but we can use the
      // edge-to-vertex map with global vertex numbers for matching instead.
      // Approach: each rank sends {global_v0, global_v1, attr} for PEC edges.
      // Get global vertex numbering from the parent ParMesh (collective MPI operation).
      mfem::Array<HYPRE_BigInt> pvert_gi;
      parent_mesh.GetGlobalVertexIndices(pvert_gi);

      struct GVEdge
      {
        long long gv0, gv1;
        int attr;
      };
      std::vector<GVEdge> local_pec_gv_edges;
      {
        mfem::Array<int> edges, orientations, ev;
        for (int be = 0; be < parent_mesh.GetNBE(); be++)
        {
          int attr = parent_mesh.GetBdrAttribute(be);
          if (std::find(pec_attrs.begin(), pec_attrs.end(), attr) == pec_attrs.end())
          {
            continue;
          }
          parent_mesh.GetBdrElementEdges(be, edges, orientations);
          for (int j = 0; j < edges.Size(); j++)
          {
            parent_mesh.GetEdgeVertices(edges[j], ev);
            long long gv0 = pvert_gi[ev[0]];
            long long gv1 = pvert_gi[ev[1]];
            local_pec_gv_edges.push_back(
                {std::min(gv0, gv1), std::max(gv0, gv1), attr});
          }
        }
      }

      // Gather all PEC global-vertex edges on rank 0.
      int local_count = static_cast<int>(local_pec_gv_edges.size());
      int nprocs_pec = Mpi::Size(comm);
      std::vector<int> all_counts(nprocs_pec);
      MPI_Gather(&local_count, 1, MPI_INT, all_counts.data(), 1, MPI_INT, 0, comm);

      int total_count = 0;
      std::vector<int> displs(nprocs_pec, 0);
      if (Mpi::Root(comm))
      {
        for (int i = 0; i < nprocs_pec; i++)
        {
          displs[i] = total_count;
          total_count += all_counts[i];
        }
      }

      // Pack as flat array: {gv0, gv1, attr} per edge (3 longs per entry, but use int
      // for simplicity — global vertex indices fit in int for typical meshes).
      std::vector<int> local_flat(local_count * 3);
      for (int i = 0; i < local_count; i++)
      {
        local_flat[3 * i + 0] = static_cast<int>(local_pec_gv_edges[i].gv0);
        local_flat[3 * i + 1] = static_cast<int>(local_pec_gv_edges[i].gv1);
        local_flat[3 * i + 2] = local_pec_gv_edges[i].attr;
      }
      std::vector<int> all_flat;
      std::vector<int> flat_counts(nprocs_pec), flat_displs(nprocs_pec);
      if (Mpi::Root(comm))
      {
        for (int i = 0; i < nprocs_pec; i++)
        {
          flat_counts[i] = all_counts[i] * 3;
          flat_displs[i] = displs[i] * 3;
        }
        all_flat.resize(total_count * 3);
      }
      MPI_Gatherv(local_flat.data(), local_count * 3, MPI_INT, all_flat.data(),
                  flat_counts.data(), flat_displs.data(), MPI_INT, 0, comm);

      // Build the submesh edge → global vertex pair mapping on ALL ranks (each rank
      // only has submesh edges for its local portion). Gather to rank 0 for matching.
      // Each submesh edge maps parent_edge → (global_v0, global_v1) → (sub_v0, sub_v1).
      std::vector<int> local_sub_edges;  // {gv0, gv1, sv0, sv1} per edge
      {
        const mfem::Array<int> &parent_edge_map = par_submesh.GetParentEdgeIDMap();
        for (int i = 0; i < parent_edge_map.Size(); i++)
        {
          int parent_edge = parent_edge_map[i];
          mfem::Array<int> pev;
          parent_mesh.GetEdgeVertices(parent_edge, pev);
          int gv0 = static_cast<int>(pvert_gi[pev[0]]);
          int gv1 = static_cast<int>(pvert_gi[pev[1]]);
          mfem::Array<int> sev;
          par_submesh.GetEdgeVertices(i, sev);
          local_sub_edges.insert(local_sub_edges.end(),
                                 {std::min(gv0, gv1), std::max(gv0, gv1),
                                  sev[0], sev[1]});
        }
      }

      // Gather submesh edge mappings on rank 0.
      int local_sub_count = static_cast<int>(local_sub_edges.size()) / 4;
      std::vector<int> sub_counts(nprocs_pec);
      MPI_Gather(&local_sub_count, 1, MPI_INT, sub_counts.data(), 1, MPI_INT, 0, comm);

      std::vector<int> all_sub_flat;
      {
        int sub_total = 0;
        std::vector<int> sub_flat_counts(nprocs_pec), sub_flat_displs(nprocs_pec);
        for (int i = 0; i < nprocs_pec; i++)
        {
          sub_flat_displs[i] = sub_total * 4;
          sub_flat_counts[i] = sub_counts[i] * 4;
          sub_total += sub_counts[i];
        }
        if (Mpi::Root(comm))
        {
          all_sub_flat.resize(sub_total * 4);
        }
        MPI_Gatherv(local_sub_edges.data(), local_sub_count * 4, MPI_INT,
                    all_sub_flat.data(), sub_flat_counts.data(), sub_flat_displs.data(),
                    MPI_INT, 0, comm);
      }

      // On rank 0: match PEC global-vertex edges against submesh edges.
      if (Mpi::Root(comm))
      {
        // Build global vertex pair → submesh vertex pair map from gathered data.
        std::map<std::pair<int, int>, std::pair<int, int>> gvpair_to_submesh_verts;
        int sub_total = static_cast<int>(all_sub_flat.size()) / 4;
        for (int i = 0; i < sub_total; i++)
        {
          auto key = std::make_pair(all_sub_flat[4 * i + 0], all_sub_flat[4 * i + 1]);
          gvpair_to_submesh_verts[key] = {all_sub_flat[4 * i + 2],
                                          all_sub_flat[4 * i + 3]};
        }

        // Match PEC edges against submesh edges.
        std::set<std::pair<int, int>> seen;
        for (int i = 0; i < total_count; i++)
        {
          int gv0 = all_flat[3 * i + 0];
          int gv1 = all_flat[3 * i + 1];
          int attr = all_flat[3 * i + 2];
          auto key = std::make_pair(std::min(gv0, gv1), std::max(gv0, gv1));
          auto it = gvpair_to_submesh_verts.find(key);
          if (it != gvpair_to_submesh_verts.end())
          {
            auto sv = it->second;
            int sv0 = std::min(sv.first, sv.second);
            int sv1 = std::max(sv.first, sv.second);
            if (seen.insert({sv0, sv1}).second)
            {
              pec_internal_edges.push_back({sv.first, sv.second, attr});
            }
          }
        }
      }
      Mpi::Print(" Found {:d} PEC internal edges on the cross-section\n",
                 static_cast<int>(pec_internal_edges.size()));
    }

    // Step 5b: Gather the parallel submesh into a serial mesh, project to 2D, redistribute
    // across all MPI processes. PrintAsOne replaces domain element attributes with
    // (rank + 1), so we gather the real attributes separately and restore them.
    {
      // Gather element and boundary element attributes via MPI.
      int ne_local = par_submesh.GetNE();
      int nbe_local = par_submesh.GetNBE();
      std::vector<int> local_elem_attrs(ne_local), local_bdr_attrs(nbe_local);
      for (int i = 0; i < ne_local; i++)
      {
        local_elem_attrs[i] = par_submesh.GetAttribute(i);
      }
      for (int i = 0; i < nbe_local; i++)
      {
        local_bdr_attrs[i] = par_submesh.GetBdrAttribute(i);
      }

      int nprocs = Mpi::Size(comm);
      std::vector<int> ne_counts(nprocs), nbe_counts(nprocs);
      MPI_Gather(&ne_local, 1, MPI_INT, ne_counts.data(), 1, MPI_INT, 0, comm);
      MPI_Gather(&nbe_local, 1, MPI_INT, nbe_counts.data(), 1, MPI_INT, 0, comm);

      std::vector<int> ne_displs(nprocs, 0), nbe_displs(nprocs, 0);
      int ne_total = 0, nbe_total = 0;
      for (int i = 0; i < nprocs; i++)
      {
        ne_displs[i] = ne_total;
        ne_total += ne_counts[i];
        nbe_displs[i] = nbe_total;
        nbe_total += nbe_counts[i];
      }

      std::vector<int> global_elem_attrs(ne_total), global_bdr_attrs(nbe_total);
      MPI_Gatherv(local_elem_attrs.data(), ne_local, MPI_INT,
                  global_elem_attrs.data(), ne_counts.data(), ne_displs.data(),
                  MPI_INT, 0, comm);
      MPI_Gatherv(local_bdr_attrs.data(), nbe_local, MPI_INT,
                  global_bdr_attrs.data(), nbe_counts.data(), nbe_displs.data(),
                  MPI_INT, 0, comm);

      // Broadcast gathered attributes to all ranks.
      MPI_Bcast(&ne_total, 1, MPI_INT, 0, comm);
      MPI_Bcast(&nbe_total, 1, MPI_INT, 0, comm);
      if (Mpi::Rank(comm) != 0)
      {
        global_elem_attrs.resize(ne_total);
        global_bdr_attrs.resize(nbe_total);
      }
      MPI_Bcast(global_elem_attrs.data(), ne_total, MPI_INT, 0, comm);
      MPI_Bcast(global_bdr_attrs.data(), nbe_total, MPI_INT, 0, comm);

      // PrintAsOne gathers mesh topology on rank 0 (collective).
      std::ostringstream oss;
      par_submesh.PrintAsOne(oss);
      std::string serial_str = oss.str();
      int str_size = static_cast<int>(serial_str.size());
      MPI_Bcast(&str_size, 1, MPI_INT, 0, comm);
      if (Mpi::Rank(comm) != 0)
      {
        serial_str.resize(str_size);
      }
      MPI_Bcast(serial_str.data(), str_size, MPI_CHAR, 0, comm);

      // Read the serial mesh and restore correct attributes.
      std::istringstream iss(serial_str);
      mfem::Mesh serial_mesh(iss);
      MFEM_VERIFY(serial_mesh.GetNE() == ne_total,
                  "PrintAsOne element count mismatch!");
      for (int i = 0; i < ne_total; i++)
      {
        serial_mesh.SetAttribute(i, global_elem_attrs[i]);
      }
      for (int i = 0; i < std::min(nbe_total, serial_mesh.GetNBE()); i++)
      {
        serial_mesh.SetBdrAttribute(i, global_bdr_attrs[i]);
      }
      serial_mesh.SetAttributes();

      // Broadcast PEC internal edges to all ranks (collected on rank 0 only).
      {
        int n_pec = static_cast<int>(pec_internal_edges.size());
        MPI_Bcast(&n_pec, 1, MPI_INT, 0, comm);
        if (Mpi::Rank(comm) != 0)
        {
          pec_internal_edges.resize(n_pec);
        }
        MPI_Bcast(pec_internal_edges.data(), n_pec * 3, MPI_INT, 0, comm);
      }

      // Add PEC internal edges as boundary elements in the serial mesh. These are
      // material interface edges that correspond to metal traces. Adding them before
      // repartitioning ensures the ParMesh has proper boundary elements for PEC.
      if (!pec_internal_edges.empty())
      {
        int dim_s = serial_mesh.Dimension();
        int sdim_s = serial_mesh.SpaceDimension();
        int nv_s = serial_mesh.GetNV();
        int ne_s = serial_mesh.GetNE();
        int nbe_s = serial_mesh.GetNBE();
        int n_pec = static_cast<int>(pec_internal_edges.size());

        mfem::Mesh new_serial(dim_s, nv_s, ne_s, nbe_s + n_pec, sdim_s);
        for (int v = 0; v < nv_s; v++)
        {
          new_serial.AddVertex(serial_mesh.GetVertex(v));
        }
        for (int e = 0; e < ne_s; e++)
        {
          new_serial.AddElement(serial_mesh.GetElement(e)->Duplicate(&new_serial));
        }
        for (int be = 0; be < nbe_s; be++)
        {
          new_serial.AddBdrElement(
              serial_mesh.GetBdrElement(be)->Duplicate(&new_serial));
        }
        for (const auto &edge : pec_internal_edges)
        {
          int vs[2] = {edge[0], edge[1]};
          new_serial.AddBdrSegment(vs, edge[2]);
        }
        new_serial.FinalizeTopology();
        // For linear meshes (order 1), vertex coordinates are the mesh nodes.
        // FinalizeTopology preserved them from AddVertex. Call EnsureNodes to
        // create the Nodes GridFunction from vertex data.
        new_serial.EnsureNodes();
        serial_mesh = std::move(new_serial);
      }

      // Debug: check element attributes.
      if (Mpi::Root(comm))
      {
        std::map<int, int> attr_counts;
        for (int i = 0; i < serial_mesh.GetNE(); i++)
        {
          attr_counts[serial_mesh.GetAttribute(i)]++;
        }
        for (const auto &[a, c] : attr_counts)
        {
          Mpi::Print(" Domain attr {:d}: {:d} elements\n", a, c);
        }
        std::map<int, int> bdr_counts;
        for (int i = 0; i < serial_mesh.GetNBE(); i++)
        {
          bdr_counts[serial_mesh.GetBdrAttribute(i)]++;
        }
        for (const auto &[a, c] : bdr_counts)
        {
          Mpi::Print(" Boundary attr {:d}: {:d} elements\n", a, c);
        }
      }

      // Project the serial mesh from 3D to 2D coordinates.
      {
        // Build orthonormal tangent frame from the surface normal.
        submesh_e1.SetSize(3);
        submesh_e2.SetSize(3);
        {
          int min_idx = 0;
          double min_val = std::abs(surface_normal(0));
          for (int d = 1; d < 3; d++)
          {
            if (std::abs(surface_normal(d)) < min_val)
            {
              min_val = std::abs(surface_normal(d));
              min_idx = d;
            }
          }
          mfem::Vector axis(3);
          axis = 0.0;
          axis(min_idx) = 1.0;
          submesh_e1(0) = axis(1) * surface_normal(2) - axis(2) * surface_normal(1);
          submesh_e1(1) = axis(2) * surface_normal(0) - axis(0) * surface_normal(2);
          submesh_e1(2) = axis(0) * surface_normal(1) - axis(1) * surface_normal(0);
          submesh_e1 /= submesh_e1.Norml2();
          submesh_e2(0) =
              surface_normal(1) * submesh_e1(2) - surface_normal(2) * submesh_e1(1);
          submesh_e2(1) =
              surface_normal(2) * submesh_e1(0) - surface_normal(0) * submesh_e1(2);
          submesh_e2(2) =
              surface_normal(0) * submesh_e1(1) - surface_normal(1) * submesh_e1(0);
          submesh_e2 /= submesh_e2.Norml2();
        }

        // Compute centroid.
        submesh_centroid.SetSize(3);
        submesh_centroid = 0.0;
        int nv = serial_mesh.GetNV();
        for (int i = 0; i < nv; i++)
        {
          const double *v = serial_mesh.GetVertex(i);
          for (int d = 0; d < 3; d++)
          {
            submesh_centroid(d) += v[d];
          }
        }
        if (nv > 0)
        {
          submesh_centroid /= nv;
        }

        // Project node coordinates to 2D.
        int mesh_order = 1;
        if (serial_mesh.GetNodes())
        {
          mesh_order = serial_mesh.GetNodes()->FESpace()->GetMaxElementOrder();
        }
        mfem::GridFunction *nodes3d = serial_mesh.GetNodes();
        const int sdim = serial_mesh.SpaceDimension();
        if (nodes3d)
        {
          const int npts = nodes3d->Size() / sdim;
          std::vector<std::array<double, 2>> projected(npts);
          for (int i = 0; i < npts; i++)
          {
            double coords[3] = {0.0, 0.0, 0.0};
            for (int d = 0; d < sdim; d++)
            {
              int idx = (nodes3d->FESpace()->GetOrdering() == mfem::Ordering::byNODES)
                            ? d * npts + i
                            : i * sdim + d;
              coords[d] = (*nodes3d)(idx);
            }
            double dx = coords[0] - submesh_centroid(0);
            double dy = coords[1] - submesh_centroid(1);
            double dz = coords[2] - submesh_centroid(2);
            projected[i][0] = dx * submesh_e1(0) + dy * submesh_e1(1) + dz * submesh_e1(2);
            projected[i][1] = dx * submesh_e2(0) + dy * submesh_e2(1) + dz * submesh_e2(2);
          }
          serial_mesh.SetCurvature(mesh_order, false, 2, mfem::Ordering::byNODES);
          mfem::GridFunction *nodes2d = serial_mesh.GetNodes();
          for (int i = 0; i < npts; i++)
          {
            (*nodes2d)(0 * npts + i) = projected[i][0];
            (*nodes2d)(1 * npts + i) = projected[i][1];
          }
        }
      }

      // Repartition across all MPI ranks and construct a distributed ParMesh.
      auto *partitioning = serial_mesh.GeneratePartitioning(nprocs);
      auto mesh_2d = std::make_unique<mfem::ParMesh>(comm, serial_mesh, partitioning);
      delete[] partitioning;

      submesh_holder = std::make_unique<Mesh>(std::move(mesh_2d));
      solve_mesh = submesh_holder.get();
    }

    // Construct MaterialOperator from the standalone 2D mesh.
    owned_mat_op = std::make_unique<MaterialOperator>(iodata, *solve_mesh);

    // Rotate material tensors from global frame to the local tangent frame.
    owned_mat_op->RotateMaterialTensors(iodata, submesh_e1, submesh_e2, surface_normal);


  }
  else
  {
    solve_mesh = mesh.back().get();
    owned_mat_op = std::make_unique<MaterialOperator>(iodata, *mesh.back());
  }
  MaterialOperator &mat_op = *owned_mat_op;

  // Construct FE spaces: ND for tangential E, H1 for normal (out-of-plane) E.
  auto nd_fec = std::make_unique<mfem::ND_FECollection>(iodata.solver.order,
                                                        solve_mesh->Dimension());
  auto h1_fec = std::make_unique<mfem::H1_FECollection>(iodata.solver.order,
                                                        solve_mesh->Dimension());
  FiniteElementSpace nd_fespace(*solve_mesh, nd_fec.get());
  FiniteElementSpace h1_fespace(*solve_mesh, h1_fec.get());

  // Essential (Dirichlet) BCs: PEC + AuxPEC + conductivity boundaries. For submesh mode
  // analysis, other wave port boundaries on the cross-section are also treated as PEC.
  mfem::Array<int> nd_dbc_tdof_list, h1_dbc_tdof_list;
  {
    const auto &pmesh = solve_mesh->Get();
    int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
    mfem::Array<int> dbc_bcs;
    for (auto attr : iodata.boundaries.pec.attributes)
    {
      if (attr > 0 && attr <= bdr_attr_max)
      {
        dbc_bcs.Append(attr);
      }
    }
    for (auto attr : iodata.boundaries.auxpec.attributes)
    {
      if (attr > 0 && attr <= bdr_attr_max)
      {
        dbc_bcs.Append(attr);
      }
    }
    for (const auto &data : iodata.boundaries.conductivity)
    {
      for (auto attr : data.attributes)
      {
        if (attr > 0 && attr <= bdr_attr_max)
        {
          dbc_bcs.Append(attr);
        }
      }
    }
    if (use_submesh)
    {
      // Other wave port boundaries become PEC on this cross-section.
      for (const auto &[idx, data] : iodata.boundaries.waveport)
      {
        for (auto attr : data.attributes)
        {
          if (std::find(ma_data.attributes.begin(), ma_data.attributes.end(), attr) !=
              ma_data.attributes.end())
          {
            continue;
          }
          if (attr > 0 && attr <= bdr_attr_max)
          {
            dbc_bcs.Append(attr);
          }
        }
      }
    }
    dbc_bcs.Sort();
    dbc_bcs.Unique();

    auto dbc_marker = mesh::AttrToMarker(bdr_attr_max, dbc_bcs);
    nd_fespace.Get().GetEssentialTrueDofs(dbc_marker, nd_dbc_tdof_list);
    h1_fespace.Get().GetEssentialTrueDofs(dbc_marker, h1_dbc_tdof_list);

  }

  const int nd_size = nd_fespace.GetTrueVSize();
  const int h1_size = h1_fespace.GetTrueVSize();
  Mpi::Print(" ND space: {:d} DOFs, H1 space: {:d} DOFs, total: {:d}\n",
             nd_fespace.GlobalTrueVSize(), h1_fespace.GlobalTrueVSize(),
             nd_fespace.GlobalTrueVSize() + h1_fespace.GlobalTrueVSize());

  // Material operator and PostOperator.
  ModeAnalysisOperator mode_op(mat_op, nd_fespace, h1_fespace, *solve_mesh,
                               iodata.solver.order);
  PostOperator<ProblemType::MODEANALYSIS> post_op(iodata, mode_op);

  // Project impedance path coordinates from 3D to the 2D local frame.
  if (use_submesh)
  {
    post_op.ProjectImpedancePaths(submesh_centroid, submesh_e1, submesh_e2);
  }

  // Error estimator setup (skip for submesh for now — isolating boundary dim issue).
  std::unique_ptr<mfem::RT_FECollection> rt_fec;
  std::unique_ptr<FiniteElementSpaceHierarchy> nd_fespaces_h, rt_fespaces_h,
      h1_fespaces_est;
  std::unique_ptr<TimeDependentFluxErrorEstimator<ComplexVector>> estimator;
  if (!use_submesh)
  {
    rt_fec = std::make_unique<mfem::RT_FECollection>(iodata.solver.order - 1,
                                                     solve_mesh->Dimension());
    nd_fespaces_h = std::make_unique<FiniteElementSpaceHierarchy>(
        std::make_unique<FiniteElementSpace>(*solve_mesh, nd_fec.get()));
    rt_fespaces_h = std::make_unique<FiniteElementSpaceHierarchy>(
        std::make_unique<FiniteElementSpace>(*solve_mesh, rt_fec.get()));
    h1_fespaces_est = std::make_unique<FiniteElementSpaceHierarchy>(
        std::make_unique<FiniteElementSpace>(*solve_mesh, h1_fec.get()));
    estimator = std::make_unique<TimeDependentFluxErrorEstimator<ComplexVector>>(
        mat_op, *nd_fespaces_h, *rt_fespaces_h, iodata.solver.linear.estimator_tol,
        iodata.solver.linear.estimator_max_it, 0, iodata.solver.linear.estimator_mg,
        &mode_op.GetCurlSpace(), h1_fespaces_est.get());
  }
  ErrorIndicator indicator;

  // Determine kn_target from user-specified n_eff target or material properties.
  double kn_target;
  if (ma_data.target > 0.0)
  {
    kn_target = ma_data.target * omega;
    Mpi::Print(" Target n_eff = {:.6e}, kn_target = {:.6e}\n", ma_data.target, kn_target);
  }
  else
  {
    double c_min = mat_op.GetLightSpeedMax().Min();
    Mpi::GlobalMin(1, &c_min, nd_fespace.GetComm());
    MFEM_VERIFY(c_min > 0.0 && c_min < mfem::infinity(),
                "Invalid material speed of light!");
    kn_target = omega / c_min * std::sqrt(1.1);
    Mpi::Print(" Auto kn_target = {:.6e} (from c_min = {:.6e})\n", kn_target, c_min);
  }

  // Construct boundary operators for impedance, absorbing, and conductivity BCs.
  // For submesh mode analysis, these operators can't be constructed on 1D boundary edges
  // with 2D material tensors. The BoundaryModeOperator handles impedance/absorbing BCs
  // internally through material coefficients instead.
  std::unique_ptr<SurfaceImpedanceOperator> surf_z_op;
  std::unique_ptr<FarfieldBoundaryOperator> farfield_op;
  std::unique_ptr<SurfaceConductivityOperator> surf_sigma_op;
  if (!use_submesh)
  {
    surf_z_op = std::make_unique<SurfaceImpedanceOperator>(iodata, mat_op, *solve_mesh);
    farfield_op = std::make_unique<FarfieldBoundaryOperator>(iodata, mat_op, *solve_mesh);
    surf_sigma_op =
        std::make_unique<SurfaceConductivityOperator>(iodata, mat_op, *solve_mesh);
  }

  // Configure solver.
  BoundaryModeOperatorConfig config;
  config.attr_to_material = &mat_op.GetAttributeToMaterial();
  config.inv_permeability = &mat_op.GetInvPermeability();
  config.curlcurl_inv_permeability = &mat_op.GetCurlCurlInvPermeability();
  config.permittivity_real = &mat_op.GetPermittivityReal();
  config.permittivity_scalar = &mat_op.GetPermittivityScalar();
  config.normal = nullptr;  // 2D domain mesh (projected for submesh case)
  config.permittivity_imag =
      mat_op.HasLossTangent() ? &mat_op.GetPermittivityImag() : nullptr;
  config.permittivity_imag_scalar =
      mat_op.HasLossTangent() ? &mat_op.GetPermittivityImagScalar() : nullptr;
  config.has_loss_tangent = mat_op.HasLossTangent();
  config.conductivity = mat_op.HasConductivity() ? &mat_op.GetConductivity() : nullptr;
  config.has_conductivity = mat_op.HasConductivity();
  config.inv_london_depth =
      mat_op.HasLondonDepth() ? &mat_op.GetInvLondonDepth() : nullptr;
  config.has_london_depth = mat_op.HasLondonDepth();
  config.mat_op = &mat_op;
  config.surf_z_op = surf_z_op.get();
  config.farfield_op = farfield_op.get();
  config.surf_sigma_op = surf_sigma_op.get();
  config.num_modes = num_modes;
  config.num_vec = -1;
  config.eig_tol = tol;
  config.which_eig = (ma_data.target > 0.0)
                         ? EigenvalueSolver::WhichType::LARGEST_MAGNITUDE
                         : EigenvalueSolver::WhichType::LARGEST_REAL;
  config.linear = &iodata.solver.linear;
  config.eigen_backend = ma_data.type;
  config.verbose = iodata.problem.verbose;

  // Build combined dbc_tdof_list for the block system (PEC only on both ND and H1).
  mfem::Array<int> dbc_tdof_list;
  dbc_tdof_list.Append(nd_dbc_tdof_list);
  for (int i = 0; i < h1_dbc_tdof_list.Size(); i++)
  {
    dbc_tdof_list.Append(nd_size + h1_dbc_tdof_list[i]);
  }

  // Construct the boundary mode solver.
  BoundaryModeOperator mode_solver(config, nd_fespace, h1_fespace, dbc_tdof_list,
                                   solve_mesh->GetComm());

  // Store Btt for impedance postprocessing.
  mode_op.SetBttMatrix(std::make_unique<mfem::HypreParMatrix>(*mode_solver.GetBtt()));

  // Solve the GEP. The eigenvalue is lambda = -kn^2 (shift-and-invert with
  // sigma = -kn_target^2).
  BlockTimer bt1(Timer::EPS);
  Mpi::Print("\nSolving GEP for {:d} propagation mode(s)...\n", num_modes);

  double sigma = -kn_target * kn_target;
  auto result = mode_solver.Solve(omega, sigma);
  int num_conv = result.num_converged;
  Mpi::Print(" Found {:d} converged eigenvalue{} (sigma = {:.6e})\n", num_conv,
             (num_conv != 1) ? "s" : "", sigma);

  // Print ALL converged eigenvalues. The shift-and-invert eigenvalue lambda is related
  // to the original eigenvalue -kn^2 by: kn^2 = -sigma - 1/lambda.
  for (int i = 0; i < num_conv; i++)
  {
    std::complex<double> lambda = mode_solver.GetEigenvalue(i);
    std::complex<double> kn = std::sqrt(-sigma - 1.0 / lambda);
    Mpi::Print(" eig {:d}: kn = {:.6e}{:+.6e}i, n_eff = {:.6e}{:+.6e}i\n", i,
               kn.real(), kn.imag(), kn.real() / omega, kn.imag() / omega);
  }

  // Postprocessing.
  BlockTimer bt2(Timer::POSTPRO);
  if (const auto *ksp = mode_solver.GetLinearSolver())
  {
    SaveMetadata(*ksp);
  }
  Mpi::Print("\nComputing solution error estimates and performing postprocessing\n\n");

  auto &l2_curl_fespace = mode_op.GetCurlSpace();
  const auto &CurlOp = l2_curl_fespace.GetDiscreteInterpolator(nd_fespace);
  const int l2_size = l2_curl_fespace.GetTrueVSize();

  const int n_print = std::min(num_conv, num_modes);
  for (int i = 0; i < n_print; i++)
  {
    // Recover kn from eigenvalue (shift-and-invert recovery).
    std::complex<double> lambda = mode_solver.GetEigenvalue(i);
    std::complex<double> kn = std::sqrt(-sigma - 1.0 / lambda);
    double error_bkwd = mode_solver.GetError(i, EigenvalueSolver::ErrorType::BACKWARD);
    double error_abs = mode_solver.GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE);

    // Extract et (ND) and en_tilde (H1) from eigenvector.
    ComplexVector et(nd_size), en(h1_size);
    {
      ComplexVector e0(nd_size + h1_size);
      mode_solver.GetEigenvector(i, e0);
      std::copy_n(e0.Real().begin(), nd_size, et.Real().begin());
      std::copy_n(e0.Imag().begin(), nd_size, et.Imag().begin());
      std::copy_n(e0.Real().begin() + nd_size, h1_size, en.Real().begin());
      std::copy_n(e0.Imag().begin() + nd_size, h1_size, en.Imag().begin());
    }

    // Power-normalize eigenvector.
    {
      const auto *Btt = mode_op.GetBtt();
      if (Btt)
      {
        Vector Btt_etr(nd_size), Btt_eti(nd_size);
        Btt->Mult(et.Real(), Btt_etr);
        Btt->Mult(et.Imag(), Btt_eti);
        double p_rr = mfem::InnerProduct(nd_fespace.GetComm(), et.Real(), Btt_etr);
        double p_ii = mfem::InnerProduct(nd_fespace.GetComm(), et.Imag(), Btt_eti);
        double p_ri = mfem::InnerProduct(nd_fespace.GetComm(), et.Real(), Btt_eti);
        double p_ir = mfem::InnerProduct(nd_fespace.GetComm(), et.Imag(), Btt_etr);
        std::complex<double> etH_Btt_et(p_rr + p_ii, p_ri - p_ir);
        std::complex<double> P = 0.5 * std::conj(kn) / omega * etH_Btt_et;
        double P_abs = std::abs(P);
        if (P_abs > 0.0)
        {
          double scale = 1.0 / std::sqrt(P_abs);
          et.Real() *= scale;
          et.Imag() *= scale;
          en.Real() *= scale;
          en.Imag() *= scale;
        }
        else
        {
          double norm2 = p_rr + p_ii;
          if (norm2 > 0.0)
          {
            double scale = 1.0 / std::sqrt(norm2);
            et.Real() *= scale;
            et.Imag() *= scale;
            en.Real() *= scale;
            en.Imag() *= scale;
          }
        }
      }
    }

    auto total_domain_energy =
        post_op.MeasureAndPrintAll(i, et, en, kn, omega, error_abs, error_bkwd, n_print);

    const bool is_propagating =
        std::abs(kn.imag()) < 0.1 * std::abs(kn.real()) && std::abs(kn.real()) > 0.0;
    if (i < num_modes && is_propagating && estimator)
    {
      ComplexVector bz(l2_size);
      {
        Vector curl_etr(l2_size), curl_eti(l2_size);
        CurlOp.Mult(et.Real(), curl_etr);
        CurlOp.Mult(et.Imag(), curl_eti);
        bz.Real() = curl_eti;
        bz.Real() *= 1.0 / omega;
        bz.Imag() = curl_etr;
        bz.Imag() *= -1.0 / omega;
      }
      estimator->AddErrorIndicator(et, bz, total_domain_energy, indicator);
    }
  }
  Mpi::Print("\n");

  post_op.MeasureFinalize(indicator);
  return {indicator, nd_fespace.GlobalTrueVSize() + h1_fespace.GlobalTrueVSize()};
}

}  // namespace palace
