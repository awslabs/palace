// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "facenbrexchange.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <memory>
#include <mutex>
#include "fem/fespace.hpp"
#include "fem/libceed/basis.hpp"
#include "fem/libceed/ceed.hpp"
#include "fem/libceed/functional.hpp"
#include "fem/libceed/restriction.hpp"
#include "fem/mesh.hpp"
#include "utils/communication.hpp"
#include "utils/diagnostic.hpp"

PalacePragmaDiagnosticPush
PalacePragmaDiagnosticDisableUnused

#include "fem/qfunctions/22/eval_22_qf.h"
#include "fem/qfunctions/33/eval_33_qf.h"

PalacePragmaDiagnosticPop

namespace palace
{

namespace
{

// Key identifying one export point-evaluator group. The request supplies an
// integer/topological point_key (reference-face topology/orientation/subface identity),
// so grouping never depends on rounded physical or reference point coordinates. Ad-hoc
// requests with an empty point_key receive a unique request-order key.
using PointConfigKey = std::vector<long long>;

// Message tags for the setup (payload size, payload) and evaluation exchanges. A
// single evaluation tag is sufficient: all processes perform the evaluation calls of
// every exchange object in the same order, so per-pair messages match in order (MPI
// non-overtaking guarantee).
constexpr int TAG_SETUP_SIZE = 1741, TAG_SETUP_PAYLOAD = 1742, TAG_EVAL = 1743;

// Registry of evaluation point integration rules with application lifetime (as in
// output_functionals.cpp): mfem::FiniteElement::GetDofToQuad caches tabulations keyed by
// the IntegrationRule pointer inside the (global, shared) FiniteElement objects, so
// destroying an IntegrationRule which was used for tabulation would leave a dangling
// cache entry.
const mfem::IntegrationRule *GetRegisteredIr(const PointConfigKey &key,
                                             const std::vector<mfem::IntegrationPoint> &pts)
{
  static std::map<PointConfigKey, std::unique_ptr<mfem::IntegrationRule>> registry;
  static std::mutex registry_mutex;
  std::lock_guard<std::mutex> lock(registry_mutex);
  auto it = registry.find(key);
  if (it == registry.end())
  {
    auto ir = std::make_unique<mfem::IntegrationRule>(static_cast<int>(pts.size()));
    for (std::size_t q = 0; q < pts.size(); q++)
    {
      ir->IntPoint(static_cast<int>(q)) = pts[q];
    }
    it = registry.emplace(key, std::move(ir)).first;
  }
  return it->second.get();
}

bool CeedSupportsNonTensorAtPoints(Ceed ceed)
{
  const char *resource;
  PalaceCeedCall(ceed, CeedGetResource(ceed, &resource));
  return !std::getenv("PALACE_SURFACE_DISABLE_ATPOINTS") &&
         (std::strstr(resource, "/gpu/cuda/ref") ||
          std::strstr(resource, "/gpu/cuda/magma"));
}

int TetNumModes(int degree)
{
  return (degree + 1) * (degree + 2) * (degree + 3) / 6;
}

mfem::IntegrationRule MakeTetLatticeRule(int degree)
{
  mfem::IntegrationRule ir(TetNumModes(degree));
  int q = 0;
  if (degree == 0)
  {
    ir.IntPoint(q).Set3(0.25, 0.25, 0.25);
    ir.IntPoint(q++).weight = 1.0;
  }
  else
  {
    for (int total = 0; total <= degree; total++)
    {
      for (int i = 0; i <= total; i++)
      {
        for (int j = 0; j <= total - i; j++)
        {
          const int k = total - i - j;
          ir.IntPoint(q).Set3(static_cast<double>(i) / degree,
                              static_cast<double>(j) / degree,
                              static_cast<double>(k) / degree);
          ir.IntPoint(q++).weight = 1.0;
        }
      }
    }
  }
  return ir;
}

void InitTetBasisForAtPoints(const mfem::FiniteElement &fe, bool grad_only,
                             CeedInt num_comp, Ceed ceed, CeedBasis *basis)
{
  MFEM_VERIFY(fe.GetGeomType() == mfem::Geometry::TETRAHEDRON,
              "FaceNbrFieldExchange AtPoints export currently supports tetrahedral "
              "volume elements only!");
  // MAGMA's hardened non-tensor AtPoints basis construction (libCEED
  // cuda-nontensor-atpoints branch) requires tabulation points that overdetermine the
  // complete polynomial space for non-H1 spaces; square tabulations are rejected. Bump the
  // lattice by one degree, matching the InitTetBasisForAtPoints copy in
  // output_functionals.cpp, since this exchange builds bases for the same ND/RT field
  // spaces.
  const int degree = std::max(0, fe.GetOrder() - (grad_only ? 1 : 0) + 1);
  const mfem::IntegrationRule ir = MakeTetLatticeRule(degree);
  ceed::InitBasisAtPoints(fe, ir, num_comp, ceed, basis);
}

void CreateSequentialPointRestriction(Ceed ceed, std::size_t num_elem, int nq, int num_comp,
                                      CeedSize l_size, CeedElemRestriction *restr)
{
  const CeedInt total_pts = static_cast<CeedInt>(num_elem * nq);
  std::vector<CeedInt> offsets(num_elem + 1 + total_pts);
  for (std::size_t e = 0; e <= num_elem; e++)
  {
    offsets[e] = static_cast<CeedInt>(num_elem + 1 + e * nq);
  }
  for (CeedInt i = 0; i < total_pts; i++)
  {
    offsets[num_elem + 1 + i] = i;
  }
  PalaceCeedCall(ceed, CeedElemRestrictionCreateAtPoints(
                           ceed, static_cast<CeedInt>(num_elem), total_pts, num_comp,
                           l_size, CEED_MEM_HOST, CEED_COPY_VALUES, offsets.data(), restr));
}

}  // namespace

FaceNbrFieldExchange::FaceNbrFieldExchange(
    const Mesh &mesh,
    const std::array<const mfem::ParFiniteElementSpace *, MaxSources> &fespaces,
    const std::vector<Request> &requests)
  : comm(mesh.GetComm())
{
  const mfem::ParMesh &pmesh = mesh.Get();
  const int num_nbr = pmesh.GetNFaceNeighbors();
  const int value_dim = pmesh.SpaceDimension();
  MFEM_VERIFY(value_dim == 2 || value_dim == 3,
              "FaceNbrFieldExchange requires 2D or 3D physical-space fields!");
  MFEM_VERIFY(requests.empty() || num_nbr > 0,
              "FaceNbrFieldExchange requires face neighbor data "
              "(ParMesh::ExchangeFaceNbrData)!");

  // Route each request to its face neighbor process: ghost element index fn falls in
  // the contiguous block [offset[i], offset[i+1]) of face neighbor i, and the ghost at
  // position p = fn - offset[i] of the block is the element at position p of the
  // neighbor's send_face_nbr_elements row for this process. This is the same
  // correspondence ParMesh::ExchangeFaceNbrData (and the face neighbor dof data of
  // ParFiniteElementSpace) relies on, for both conformal and nonconforming meshes.
  const auto &elem_offsets = pmesh.face_nbr_elements_offset;
  std::vector<std::vector<int>> nbr_reqs(num_nbr);
  for (std::size_t r = 0; r < requests.size(); r++)
  {
    const auto &req = requests[r];
    MFEM_VERIFY(req.face_nbr_elem >= 0 &&
                    req.face_nbr_elem < pmesh.GetNFaceNeighborElements() &&
                    req.source_mask != 0 && !req.pts.empty(),
                "Invalid face neighbor field exchange request!");
    const int *it =
        std::upper_bound(elem_offsets.begin(), elem_offsets.end(), req.face_nbr_elem);
    nbr_reqs[static_cast<int>(it - elem_offsets.begin()) - 1].push_back(
        static_cast<int>(r));
  }

  // Serialize the requests per neighbor (as doubles: position in the neighbor's send
  // element list, number of points, source mask, integer point-key length and entries,
  // then the point coordinates), and assign the import offsets (per neighbor, in
  // request construction order, source slots in ascending order; the serving process
  // lays out the reply values in exactly
  // this order).
  import_offsets.resize(requests.size());
  for (auto &offsets : import_offsets)
  {
    offsets.fill(-1);
  }
  std::vector<std::vector<double>> send_payload(num_nbr);
  int import_size = 0;
  for (int i = 0; i < num_nbr; i++)
  {
    if (nbr_reqs[i].empty())
    {
      continue;
    }
    auto &payload = send_payload[i];
    const int block_start = import_size;
    for (const int r : nbr_reqs[i])
    {
      const auto &req = requests[r];
      const int nq = static_cast<int>(req.pts.size());
      payload.push_back(static_cast<double>(req.face_nbr_elem - elem_offsets[i]));
      payload.push_back(static_cast<double>(nq));
      payload.push_back(static_cast<double>(req.source_mask));
      payload.push_back(static_cast<double>(req.point_key.size()));
      for (auto v : req.point_key)
      {
        payload.push_back(static_cast<double>(v));
      }
      for (const auto &ip : req.pts)
      {
        payload.push_back(ip.x);
        payload.push_back(ip.y);
        payload.push_back(ip.z);
      }
      for (int s = 0; s < MaxSources; s++)
      {
        if (req.source_mask & (1u << s))
        {
          MFEM_VERIFY(fespaces[s],
                      "Missing finite element space for requested source slot!");
          import_offsets[r][s] = import_size;
          import_size += value_dim * nq;
        }
      }
    }
    recv_msgs.push_back({pmesh.GetFaceNbrRank(i), block_start, import_size - block_start});
  }

  // Symmetric setup exchange with all face neighbors (every pair exchanges, possibly
  // empty): payload sizes, then the payloads.
  std::vector<int> send_len(num_nbr), recv_len(num_nbr);
  {
    std::vector<MPI_Request> mpi_reqs(2 * num_nbr);
    for (int i = 0; i < num_nbr; i++)
    {
      send_len[i] = static_cast<int>(send_payload[i].size());
      MPI_Irecv(&recv_len[i], 1, MPI_INT, pmesh.GetFaceNbrRank(i), TAG_SETUP_SIZE, comm,
                &mpi_reqs[2 * i]);
      MPI_Isend(&send_len[i], 1, MPI_INT, pmesh.GetFaceNbrRank(i), TAG_SETUP_SIZE, comm,
                &mpi_reqs[2 * i + 1]);
    }
    MPI_Waitall(static_cast<int>(mpi_reqs.size()), mpi_reqs.data(), MPI_STATUSES_IGNORE);
  }
  std::vector<std::vector<double>> recv_payload(num_nbr);
  {
    std::vector<MPI_Request> mpi_reqs;
    mpi_reqs.reserve(2 * num_nbr);
    for (int i = 0; i < num_nbr; i++)
    {
      if (recv_len[i] > 0)
      {
        recv_payload[i].resize(recv_len[i]);
        MPI_Irecv(recv_payload[i].data(), recv_len[i], MPI_DOUBLE, pmesh.GetFaceNbrRank(i),
                  TAG_SETUP_PAYLOAD, comm, &mpi_reqs.emplace_back());
      }
      if (send_len[i] > 0)
      {
        MPI_Isend(send_payload[i].data(), send_len[i], MPI_DOUBLE, pmesh.GetFaceNbrRank(i),
                  TAG_SETUP_PAYLOAD, comm, &mpi_reqs.emplace_back());
      }
    }
    MPI_Waitall(static_cast<int>(mpi_reqs.size()), mpi_reqs.data(), MPI_STATUSES_IGNORE);
  }

  Ceed ceed = ceed::internal::GetCeedObjects()[0];
  const bool use_at_points = CeedSupportsNonTensorAtPoints(ceed);

  // Parse the received requests, assigning export offsets with the same layout rules
  // as the import offsets above. On CUDA backends, tetrahedral export evaluators use
  // libCEED AtPoints so requests with the same source/geometry/point count can share one
  // operator even when their finite NC trace maps differ. Other backends keep the
  // mapped-integration-rule grouping by integer/topological point-key data; empty keys
  // intentionally fall back to unique request-order groups.
  struct ExportGroup
  {
    bool at_points = false;
    std::vector<mfem::IntegrationPoint> pts;  // Representative IR or concatenated AtPoints
    std::vector<int> elems;
    std::vector<int> bases;  // Export vector base offset per element entry
  };
  std::map<PointConfigKey, ExportGroup> export_map;
  int export_group_id = 0;
  int export_size = 0;
  for (int i = 0; i < num_nbr; i++)
  {
    const auto &payload = recv_payload[i];
    const int block_start = export_size;
    std::size_t k = 0;
    while (k < payload.size())
    {
      const int p = static_cast<int>(std::llround(payload[k++]));
      const int nq = static_cast<int>(std::llround(payload[k++]));
      const auto mask = static_cast<unsigned int>(std::llround(payload[k++]));
      const int point_key_size = static_cast<int>(std::llround(payload[k++]));
      MFEM_VERIFY(point_key_size >= 0,
                  "Invalid face neighbor point-key size in received field exchange "
                  "request!");
      std::vector<long long> point_key(point_key_size);
      for (int q = 0; q < point_key_size; q++)
      {
        point_key[q] = std::llround(payload[k++]);
      }
      MFEM_VERIFY(p >= 0 && p < pmesh.send_face_nbr_elements.RowSize(i),
                  "Invalid face neighbor element position in received field exchange "
                  "request!");
      const int elem = pmesh.send_face_nbr_elements.GetRow(i)[p];
      std::vector<mfem::IntegrationPoint> pts(nq);
      for (int q = 0; q < nq; q++)
      {
        pts[q].Set3(payload[k], payload[k + 1], payload[k + 2]);
        pts[q].weight = 1.0;
        k += 3;
      }
      const auto geom = pmesh.GetElementGeometry(elem);
      for (int s = 0; s < MaxSources; s++)
      {
        if (mask & (1u << s))
        {
          MFEM_VERIFY(fespaces[s],
                      "Missing finite element space for received source slot!");
          const bool at_points_group = use_at_points && geom == mfem::Geometry::TETRAHEDRON;
          PointConfigKey key;
          key.reserve(5 + point_key.size());
          key.push_back(s);
          key.push_back(static_cast<long long>(geom));
          key.push_back(nq);
          key.push_back(static_cast<long long>(at_points_group));
          if (!at_points_group)
          {
            if (point_key.empty())
            {
              key.push_back(export_group_id++);
            }
            else
            {
              key.insert(key.end(), point_key.begin(), point_key.end());
            }
          }
          auto &group = export_map[key];
          group.at_points = at_points_group;
          if (at_points_group)
          {
            group.pts.insert(group.pts.end(), pts.begin(), pts.end());
          }
          else if (group.pts.empty())
          {
            group.pts = pts;
          }
          group.elems.push_back(elem);
          group.bases.push_back(export_size);
          export_size += value_dim * nq;
        }
      }
    }
    if (export_size > block_start)
    {
      send_msgs.push_back(
          {pmesh.GetFaceNbrRank(i), block_start, export_size - block_start});
    }
  }

  imported.SetSize(import_size);
  imported.UseDevice(true);
  imported = 0.0;
  exported.SetSize(export_size);
  exported.UseDevice(true);
  if (export_size == 0)
  {
    return;
  }

  // Assemble a libCEED point evaluator for each export group, writing the
  // physical-space field values (space-dimension components per point, point-major)
  // into the exported vector at the assigned offsets.
  int max_vsize = 0;
  for (const auto *fespace : fespaces)
  {
    max_vsize = std::max(max_vsize, fespace ? fespace->GetVSize() : 0);
  }
  field_staging.SetSize(max_vsize);
  field_staging.UseDevice(true);
  field_staging = 0.0;
  const mfem::FiniteElementSpace &mesh_fespace = *pmesh.GetNodes()->FESpace();
  for (const auto &[key, group] : export_map)
  {
    const int s = static_cast<int>(key[0]);
    const auto geom = static_cast<mfem::Geometry::Type>(key[1]);
    const int nq = static_cast<int>(key[2]);
    const std::size_t num_elem = group.elems.size();
    const auto &fespace = *fespaces[s];
    MFEM_VERIFY(!group.at_points ||
                    group.pts.size() == num_elem * static_cast<std::size_t>(nq),
                "Invalid AtPoints export point layout for face neighbor exchange!");
    const mfem::IntegrationRule *ir =
        group.at_points ? nullptr : GetRegisteredIr(key, group.pts);

    // Inputs: mesh node gradients (on-the-fly geometry for the Piola transformations)
    // and the field, both evaluated at the points.
    std::vector<ceed::CeedFunctionalFieldInput> inputs;
    std::vector<std::pair<std::string, int>> field_sources;
    CeedElemRestriction points_restr = nullptr;
    CeedVector points_vec = nullptr;
    if (group.at_points)
    {
      auto &points = export_attrs.emplace_back(3 * num_elem * nq);
      for (std::size_t e = 0; e < num_elem; e++)
      {
        for (int q = 0; q < nq; q++)
        {
          const auto &ip = group.pts[e * nq + q];
          const std::size_t off = 3 * (e * nq + q);
          points[off + 0] = ip.x;
          points[off + 1] = ip.y;
          points[off + 2] = ip.z;
        }
      }
      CreateSequentialPointRestriction(ceed, num_elem, nq, 3, points.Size(), &points_restr);
      ceed::InitCeedVector(points, ceed, &points_vec);
    }
    CeedElemRestriction mesh_restr = FiniteElementSpace::BuildCeedElemRestriction(
        mesh_fespace, ceed, geom, group.elems, /*is_interp*/ true);
    const mfem::FiniteElement *mesh_fe =
        mesh_fespace.FEColl()->FiniteElementForGeometry(geom);
    CeedBasis mesh_basis;
    if (group.at_points)
    {
      InitTetBasisForAtPoints(*mesh_fe, /*grad_only*/ true, mesh_fespace.GetVDim(), ceed,
                              &mesh_basis);
    }
    else
    {
      ceed::InitBasisAtPoints(*mesh_fe, *ir, mesh_fespace.GetVDim(), ceed, &mesh_basis);
    }
    CeedVector mesh_nodes_vec;
    ceed::InitCeedVector(*mesh_fespace.GetMesh()->GetNodes(), ceed, &mesh_nodes_vec);
    inputs.push_back({"x", mesh_nodes_vec, mesh_restr, mesh_basis, ceed::EvalMode::Grad});
    CeedElemRestriction field_restr;
    ceed::InitRestriction(fespace, group.elems, false, /*is_interp*/ true, false, ceed,
                          &field_restr);
    const mfem::FiniteElement *fe = fespace.FEColl()->FiniteElementForGeometry(geom);
    MFEM_VERIFY(fe, "Unable to get field finite element for face neighbor exchange!");
    CeedBasis field_basis;
    if (group.at_points)
    {
      InitTetBasisForAtPoints(*fe, /*grad_only*/ false, fespace.GetVDim(), ceed,
                              &field_basis);
    }
    else
    {
      ceed::InitBasisAtPoints(*fe, *ir, fespace.GetVDim(), ceed, &field_basis);
    }
    CeedVector field_vec;
    ceed::InitCeedVector(field_staging, ceed, &field_vec);
    inputs.push_back({"u_1", field_vec, field_restr, field_basis, ceed::EvalMode::Interp});
    field_sources.emplace_back("u_1", s);

    // Output restriction: one physical-space vector per point at the assigned export
    // offsets.
    std::vector<CeedInt> offsets(num_elem * nq);
    for (std::size_t e = 0; e < num_elem; e++)
    {
      for (int j = 0; j < nq; j++)
      {
        offsets[e * nq + j] = group.bases[e] + value_dim * j;
      }
    }
    CeedElemRestriction out_restr;
    PalaceCeedCall(ceed, CeedElemRestrictionCreate(ceed, static_cast<CeedInt>(num_elem), nq,
                                                   value_dim, 1, (CeedSize)export_size,
                                                   CEED_MEM_HOST, CEED_COPY_VALUES,
                                                   offsets.data(), &out_restr));

    // The reply contains physical-space field values: the Piola transformation
    // (H(curl) or H(div) depending on the source space) is applied here so the
    // requester needs no neighbor element geometry.
    const auto map_type = fespace.FEColl()->GetMapType(pmesh.Dimension());
    MFEM_VERIFY(map_type == mfem::FiniteElement::H_CURL ||
                    map_type == mfem::FiniteElement::H_DIV,
                "FaceNbrFieldExchange requires H(curl) or H(div) source spaces!");
    ceed::CeedQFunctionInfo info;
    if (map_type == mfem::FiniteElement::H_CURL)
    {
      if (value_dim == 2)
      {
        info.apply_qf = f_eval_probe_hcurl_22;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_eval_probe_hcurl_22_loc);
      }
      else
      {
        info.apply_qf = f_eval_probe_hcurl_33;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_eval_probe_hcurl_33_loc);
      }
    }
    else
    {
      if (value_dim == 2)
      {
        info.apply_qf = f_eval_probe_hdiv_22;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_eval_probe_hdiv_22_loc);
      }
      else
      {
        info.apply_qf = f_eval_probe_hdiv_33;
        info.apply_qf_path = PalaceQFunctionRelativePath(f_eval_probe_hdiv_33_loc);
      }
    }
    CeedOperator op;
    if (group.at_points)
    {
      ceed::AssembleCeedPointEvaluatorAtPoints(info, nullptr, 0, ceed, inputs, points_restr,
                                               points_vec, value_dim, out_restr, &op);
    }
    else
    {
      ceed::AssembleCeedPointEvaluator(info, nullptr, 0, ceed, inputs, value_dim, out_restr,
                                       &op);
    }
    export_groups.push_back({ceed, op, std::move(field_sources)});

    // Cleanup (the assembled operator holds its own references).
    if (points_vec)
    {
      PalaceCeedCall(ceed, CeedVectorDestroy(&points_vec));
    }
    if (points_restr)
    {
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&points_restr));
    }
    PalaceCeedCall(ceed, CeedVectorDestroy(&mesh_nodes_vec));
    PalaceCeedCall(ceed, CeedVectorDestroy(&field_vec));
    PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&mesh_restr));
    PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&field_restr));
    PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&out_restr));
    PalaceCeedCall(ceed, CeedBasisDestroy(&mesh_basis));
    PalaceCeedCall(ceed, CeedBasisDestroy(&field_basis));
  }
}

FaceNbrFieldExchange::~FaceNbrFieldExchange()
{
  fem::DestroyGroupOperators(export_groups);
}

void FaceNbrFieldExchange::Exchange(
    const std::array<const Vector *, MaxSources> &srcs) const
{
  // Evaluate the values requested by the neighbors (on the device).
  if (exported.Size() > 0)
  {
    exported = 0.0;
    fem::ApplyAddGroupOperators(export_groups, srcs, exported);
  }
  if (recv_msgs.empty() && send_msgs.empty())
  {
    return;
  }

  // Exchange: receive directly into the imported vector (contiguous per neighbor by
  // construction), send directly from the exported vector.
  std::vector<MPI_Request> mpi_reqs;
  mpi_reqs.reserve(recv_msgs.size() + send_msgs.size());
  double *dst = (imported.Size() > 0) ? imported.HostWrite() : nullptr;
  for (const auto &msg : recv_msgs)
  {
    MPI_Irecv(dst + msg.start, msg.size, MPI_DOUBLE, msg.rank, TAG_EVAL, comm,
              &mpi_reqs.emplace_back());
  }
  const double *src = (exported.Size() > 0) ? exported.HostRead() : nullptr;
  for (const auto &msg : send_msgs)
  {
    MPI_Isend(src + msg.start, msg.size, MPI_DOUBLE, msg.rank, TAG_EVAL, comm,
              &mpi_reqs.emplace_back());
  }
  MPI_Waitall(static_cast<int>(mpi_reqs.size()), mpi_reqs.data(), MPI_STATUSES_IGNORE);
}

}  // namespace palace
