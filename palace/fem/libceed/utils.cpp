// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "utils.hpp"

#include "basis.hpp"
#include "restriction.hpp"

#if defined(MFEM_USE_OPENMP)
#include <omp.h>
#endif

namespace palace::ceed
{

namespace internal
{

std::vector<Ceed> ceed;

}  // namespace internal

void Initialize(const char *resource, const char *jit_source_dir)
{
  PalacePragmaOmp(parallel)
  {
    PalacePragmaOmp(master)
    {
#if defined(MFEM_USE_OPENMP)
      const int nt = omp_get_num_threads();
#else
      const int nt = 1;
#endif
      internal::ceed.resize(nt, nullptr);
    }
  }

  // Master thread initializes all Ceed objects (ineherently sequential anyway due to shared
  // resources).
  for (std::size_t i = 0; i < internal::ceed.size(); i++)
  {
    int ierr = CeedInit(resource, &internal::ceed[i]);
    MFEM_VERIFY(!ierr, "Failed to initialize libCEED with resource " << resource << "!");
    Ceed ceed = internal::ceed[i];

    // Configure QFunction search path.
    if (jit_source_dir)
    {
      PalaceCeedCall(ceed, CeedAddJitSourceRoot(ceed, jit_source_dir));
    }

    // XX TODO: Do this always even not in debug?
    // Configure error handling.
#ifdef MFEM_DEBUG
    PalaceCeedCall(ceed, CeedSetErrorHandler(ceed, CeedErrorAbort));
#endif
  }
}

void Finalize()
{
  // Destroy global basis and element restriction caches.
  ClearBasisRestrictionCache();

  // Destroy Ceed context(s).
  for (std::size_t i = 0; i < internal::ceed.size(); i++)
  {
    int ierr = CeedDestroy(&internal::ceed[i]);
    MFEM_VERIFY(!ierr, "Failed to finalize libCEED!");
  }
  internal::ceed.clear();
}

std::string Print()
{
  MFEM_VERIFY(internal::ceed.size() > 0,
              "libCEED must be initialized before querying the active backend!");
  Ceed ceed = internal::ceed[0];
  const char *ceed_resource;
  PalaceCeedCall(ceed, CeedGetResource(ceed, &ceed_resource));
  return std::string(ceed_resource);
}

void ClearBasisRestrictionCache()
{
  for (auto [k, v] : internal::basis_map)
  {
    Ceed ceed;
    PalaceCeedCallBackend(CeedBasisGetCeed(v, &ceed));
    PalaceCeedCall(ceed, CeedBasisDestroy(&v));
  }
  for (auto [k, v] : internal::restr_map)
  {
    Ceed ceed;
    PalaceCeedCallBackend(CeedElemRestrictionGetCeed(v, &ceed));
    PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&v));
  }
  internal::basis_map.clear();
  internal::restr_map.clear();
}

void InitCeedVector(const mfem::Vector &v, Ceed ceed, CeedVector *cv)
{
  CeedMemType mem;
  const CeedScalar *data;
  PalaceCeedCall(ceed, CeedVectorCreate(ceed, v.Size(), cv));
  PalaceCeedCall(ceed, CeedGetPreferredMemType(ceed, &mem));
  if (mfem::Device::Allows(mfem::Backend::DEVICE_MASK) && mem == CEED_MEM_DEVICE)
  {
    data = v.Read();
  }
  else
  {
    data = v.HostRead();
    mem = CEED_MEM_HOST;
  }
  PalaceCeedCall(
      ceed, CeedVectorSetArray(*cv, mem, CEED_USE_POINTER, const_cast<CeedScalar *>(data)));
}

}  // namespace palace::ceed
