// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SURFACE_CURL_SOLVER_HPP
#define PALACE_MODELS_SURFACE_CURL_SOLVER_HPP

#include <map>
#include <memory>
#include <vector>
#include "linalg/vector.hpp"
#include "utils/labels.hpp"

namespace palace
{

class Mesh;
class FiniteElementSpace;
class MaterialOperator;

// Forward declarations
class SurfaceFluxData;
class CurlCurlOperator;
template <ProblemType T>
class PostOperator;

// London cohomology generators kept across adaptive refinements, keyed by flux loop index.
// One built before hanging edges appear is exact and stays curl-free on the film when
// prolonged. Kept as a discontinuous linear vector field (exact for ND_1), since MFEM's ND
// rebalance transfer misindexes the DOFs of migrated elements.
class CohomologyGeneratorCache
{
private:
  struct Entry
  {
    const mfem::ParMesh *mesh;
    long sequence;
    std::unique_ptr<mfem::L2_FECollection> fec;
    std::unique_ptr<mfem::ParFiniteElementSpace> fespace;
    std::unique_ptr<mfem::ParGridFunction> a;
  };
  std::map<int, Entry> entries;

public:
  CohomologyGeneratorCache();
  ~CohomologyGeneratorCache();

  // Stored generator for this loop if it is up to date on the given mesh, else nullptr.
  const mfem::ParGridFunction *Find(int key, const mfem::ParMesh &mesh) const;

  // Keep the given ND_1 generator for this loop.
  void Store(int key, const mfem::ParGridFunction &a);

  // Transfer the stored generators after a single modification (refinement or rebalance)
  // of the given mesh. Entries that cannot be transferred exactly are dropped.
  void Update(mfem::ParMesh &mesh);
};

// Build the curl-free cut cohomology generator (±Φ across a cut through the hole centroid)
// used as the London shifted-penalty drive a_h, on the 3D ND space (see
// BuildCutCohomologyGenerator). With a cache, the generator is kept across adaptive
// refinements under the given key.
Vector SolveSurfaceCurlProblem(const SurfaceFluxData &flux_data, const Mesh &mesh,
                               const FiniteElementSpace &nd_fespace,
                               PostOperator<ProblemType::MAGNETOSTATIC> &post_op,
                               CohomologyGeneratorCache *cache = nullptr, int key = 0);

void SolveSurfaceCurlProblem(const SurfaceFluxData &flux_data, const Mesh &mesh,
                             const FiniteElementSpace &nd_fespace,
                             PostOperator<ProblemType::MAGNETOSTATIC> &post_op,
                             Vector &result, CohomologyGeneratorCache *cache = nullptr,
                             int key = 0);

// Integrate the magnetic flux density B ⋅ n over the given boundary attributes, with the
// surface normal oriented according to flux_direction. Returns the global (summed) flux.
double ComputeFluxThroughSurface(const mfem::ParGridFunction &B_gf,
                                 const std::vector<int> &attributes, const Mesh &mesh,
                                 const MaterialOperator &mat_op,
                                 const mfem::Vector &flux_direction, MPI_Comm comm);

}  // namespace palace

#endif  // PALACE_MODELS_SURFACE_CURL_SOLVER_HPP
