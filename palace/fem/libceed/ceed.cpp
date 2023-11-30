// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "ceed.hpp"

#include "utils/omp.hpp"

#if defined(MFEM_USE_OPENMP)
#include <omp.h>
#endif

namespace palace::ceed
{

namespace internal
{

static std::vector<Ceed> ceeds;

const std::vector<Ceed> &GetCeedObjects()
{
  return ceeds;
}

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
