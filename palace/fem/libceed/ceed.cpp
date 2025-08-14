// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "ceed.hpp"

#include <string_view>
#include "utils/omp.hpp"

namespace palace::ceed
{

namespace internal
{

static std::vector<Ceed> ceeds;

const std::vector<Ceed> &GetCeedObjects()
{
  return ceeds;
}

std::size_t NumCeeds()
{
  return GetCeedObjects().size();
}

}  // namespace internal

void Initialize(const char *resource, const char *jit_source_dir)
{
  PalacePragmaOmp(parallel)
  {
    PalacePragmaOmp(master)
    {
      // Only parallelize libCEED operators over threads when not using the GPU.
      const int nt = !std::string_view(resource).compare(0, 4, "/cpu")
                         ? utils::GetNumActiveThreads()
                         : 1;
      internal::ceeds.resize(nt, nullptr);
    }
  }

  // Master thread initializes all Ceed objects (ineherently sequential anyway due to shared
  // resources).
  for (std::size_t i = 0; i < internal::ceeds.size(); i++)
  {
    int ierr = CeedInit(resource, &internal::ceeds[i]);
    MFEM_VERIFY(!ierr, "Failed to initialize libCEED with resource " << resource << "!");
    Ceed ceed = internal::ceeds[i];

    // Configure error handling (allow errors to be handled by PalaceCeedCallBackend or
    // PalaceCeedCall).
    PalaceCeedCall(ceed, CeedSetErrorHandler(ceed, CeedErrorStore));

    // Configure QFunction search path.
    if (jit_source_dir)
    {
      PalaceCeedCall(ceed, CeedAddJitSourceRoot(ceed, jit_source_dir));
    }
  }
}

void Finalize()
{
  // Destroy Ceed context(s).
  for (std::size_t i = 0; i < internal::ceeds.size(); i++)
  {
    int ierr = CeedDestroy(&internal::ceeds[i]);
    MFEM_VERIFY(!ierr, "Failed to finalize libCEED!");
  }
  internal::ceeds.clear();
}

std::string Print()
{
  MFEM_VERIFY(internal::GetCeedObjects().size() > 0,
              "libCEED must be initialized before querying the active backend!");
  Ceed ceed = internal::GetCeedObjects()[0];
  const char *ceed_resource;
  PalaceCeedCall(ceed, CeedGetResource(ceed, &ceed_resource));
  return std::string(ceed_resource);
}

void InitCeedVector(const mfem::Vector &v, Ceed ceed, CeedVector *cv, bool init)
{
  CeedMemType mem;
  PalaceCeedCall(ceed, CeedGetPreferredMemType(ceed, &mem));
  if (!mfem::Device::Allows(mfem::Backend::DEVICE_MASK) && mem == CEED_MEM_DEVICE)
  {
    mem = CEED_MEM_HOST;
  }
  const auto *data = v.Read(mem == CEED_MEM_DEVICE);
  if (init)
  {
    PalaceCeedCall(ceed, CeedVectorCreate(ceed, v.Size(), cv));
  }
  else
  {
    PalaceCeedCall(ceed, CeedVectorTakeArray(*cv, mem, nullptr));
  }
  PalaceCeedCall(
      ceed, CeedVectorSetArray(*cv, mem, CEED_USE_POINTER, const_cast<CeedScalar *>(data)));
}

CeedElemTopology GetCeedTopology(mfem::Geometry::Type geom)
{
  switch (geom)
  {
    case mfem::Geometry::SEGMENT:
      return CEED_TOPOLOGY_LINE;
    case mfem::Geometry::TRIANGLE:
      return CEED_TOPOLOGY_TRIANGLE;
    case mfem::Geometry::SQUARE:
      return CEED_TOPOLOGY_QUAD;
    case mfem::Geometry::TETRAHEDRON:
      return CEED_TOPOLOGY_TET;
    case mfem::Geometry::CUBE:
      return CEED_TOPOLOGY_HEX;
    case mfem::Geometry::PRISM:
      return CEED_TOPOLOGY_PRISM;
    case mfem::Geometry::PYRAMID:
      return CEED_TOPOLOGY_PYRAMID;
    default:
      MFEM_ABORT("This type of element is not supported!");
      return CEED_TOPOLOGY_LINE;  // Silence compiler warning
  }
}

mfem::Geometry::Type GetMfemTopology(CeedElemTopology geom)
{
  switch (geom)
  {
    case CEED_TOPOLOGY_LINE:
      return mfem::Geometry::SEGMENT;
    case CEED_TOPOLOGY_TRIANGLE:
      return mfem::Geometry::TRIANGLE;
    case CEED_TOPOLOGY_QUAD:
      return mfem::Geometry::SQUARE;
    case CEED_TOPOLOGY_TET:
      return mfem::Geometry::TETRAHEDRON;
    case CEED_TOPOLOGY_HEX:
      return mfem::Geometry::CUBE;
    case CEED_TOPOLOGY_PRISM:
      return mfem::Geometry::PRISM;
    case CEED_TOPOLOGY_PYRAMID:
      return mfem::Geometry::PYRAMID;
    default:
      MFEM_ABORT("This type of element is not supported!");
      return mfem::Geometry::SEGMENT;  // Silence compiler warning
  }
}

}  // namespace palace::ceed
