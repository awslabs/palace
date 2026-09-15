// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_FACE_NBR_EXCHANGE_HPP
#define PALACE_FEM_FACE_NBR_EXCHANGE_HPP

#include <array>
#include <vector>
#include <mfem.hpp>
#include "fem/ceed_group_operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class Mesh;

//
// Exchanges field values evaluated at prescribed reference points of face neighbor
// (ghost) elements between processes. This enables libCEED surface functionals to
// evaluate two-sided integrands on interior boundaries crossing parallel interfaces
// without exchanging field dofs (the legacy approach of
// ParGridFunction::ExchangeFaceNbrData): each process posts requests against ghost
// elements (by face neighbor index, with evaluation points in the neighbor element's
// reference space, e.g. mapped through FaceElementTransformations::Loc2). At setup, the
// requests are routed to the owning processes through the pairwise face neighbor
// communication structure of the ParMesh (the same correspondence
// ParMesh::ExchangeFaceNbrData relies on, valid for conformal and nonconforming
// meshes). At evaluation, the owning process evaluates the requested fields at the
// points with libCEED point evaluators (on the device) and replies with physical-space
// field values (Piola transformations applied by the evaluating process for H(curl) and
// H(div) sources, scalar VALUE/INTEGRAL sources copied directly), so the requester needs
// no neighbor element geometry; material properties are applied by the requester using
// the ghost element attributes).
//
class FaceNbrFieldExchange
{
public:
  // Number of source field slots (matching the source vectors of the surface
  // functional evaluation calls).
  static constexpr int MaxSources = 4;

  struct Request
  {
    // Face neighbor (ghost) element index, in [0, ParMesh::GetNFaceNeighborElements()),
    // e.g. FaceElementTransformations::Elem2No - ParMesh::GetNE() for a shared
    // interior face.
    int face_nbr_elem;

    // Bitmask over the source field slots to evaluate (bit s set: source slot s).
    unsigned int source_mask;

    // Canonical identity of the fixed reference point rule. SurfaceFunctional fills this
    // with its reference-face/subface transformation key, including exact point-matrix
    // coefficient bits. Requests with the same key, source slot, element geometry, and
    // point count share one ordinary mapped-rule libCEED evaluator. Empty keys are
    // allowed for ad-hoc requests and use exact coordinate bits for identity.
    std::vector<long long> point_key;

    // Evaluation point rule in the neighbor element's reference space (the ghost
    // element reference space coincides with the owning process' local element reference
    // space, as the ghost elements preserve the neighbor's vertex ordering). Exact
    // coordinate and weight bits are included in the fixed-rule registry identity.
    std::vector<mfem::IntegrationPoint> pts;
  };

private:
  MPI_Comm comm;

  // Pairwise evaluation-time messages with face neighbor ranks. Receives fill
  // [start, start + size) of the imported vector; sends read [start, start + size) of
  // the exported vector (both contiguous per neighbor by construction).
  struct Message
  {
    int rank;
    int start, size;
  };
  std::vector<Message> recv_msgs, send_msgs;

  // Number of imported components for each source slot. H(curl)/H(div) fields use the
  // mesh space dimension; scalar VALUE/INTEGRAL fields use one component.
  std::array<int, MaxSources> source_num_comp;

  // For request r (construction order), base offset of the values for source slot s in
  // the imported vector (-1 when not requested). Values are point-major with
  // source_num_comp[s] components per point.
  std::vector<std::array<int, MaxSources>> import_offsets;

  // Assembled libCEED point evaluators serving the requests of neighboring processes,
  // writing into the exported vector.
  std::vector<fem::CeedGroupOperator> export_groups;

  // Exported and imported field-value vectors.
  mutable Vector exported, imported;

public:
  // Construct the exchange for the given requests against ghost elements of the mesh.
  // fespaces[s] provides the evaluation space for source slot s (may be nullptr for
  // unused slots; all processes must pass the same spaces). Collective on the mesh
  // communicator: the requests are routed to the owning processes following the face
  // neighbor communication structure of ParMesh::ExchangeFaceNbrData (which must have
  // been called for the mesh already, see fem/mesh.cpp).
  FaceNbrFieldExchange(
      const Mesh &mesh,
      const std::array<const mfem::ParFiniteElementSpace *, MaxSources> &fespaces,
      const std::vector<Request> &requests);
  ~FaceNbrFieldExchange();

  FaceNbrFieldExchange(const FaceNbrFieldExchange &) = delete;
  FaceNbrFieldExchange &operator=(const FaceNbrFieldExchange &) = delete;

  // Size of the imported values vector.
  int ImportSize() const { return imported.Size(); }

  // Number of point-major components imported for source slot s.
  int ImportNumComp(int s) const { return source_num_comp[s]; }

  // Base offset in the imported vector of the values of source slot s for request r
  // (construction order), or -1 when source s was not requested. Layout per request
  // and source: ImportNumComp(s) components per point, point-major.
  int ImportOffset(int r, int s) const { return import_offsets[r][s]; }

  // The imported values vector (filled by Exchange()).
  const Vector &Imported() const { return imported; }

  // Detach borrowed source arrays from the export operators before caller-owned warm-up
  // storage is released. The next Exchange() reattaches its source arrays directly.
  void DetachFieldVectors() const;

  // Evaluate the exported field values for the neighbors' requests from the given
  // source vectors (L-vectors of the corresponding source spaces), exchange with the
  // face neighbor processes, and fill the imported vector with this process' requested
  // values. All processes must call in the same order with the corresponding sources
  // (point-to-point communication with the face neighbor ranks; processes with no
  // requests and no exports perform no communication).
  void Exchange(const std::array<const Vector *, MaxSources> &srcs) const;
};

}  // namespace palace

#endif  // PALACE_FEM_FACE_NBR_EXCHANGE_HPP
