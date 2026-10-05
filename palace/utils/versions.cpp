// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "versions.hpp"

#include <cstdint>

#include <Eigen/Core>
#include <ceed.h>
#include <metis.h>
#include <mpi.h>
#include <parmetis.h>
#include <fmt/format.h>
#include <scn/scan.h>
#include <mfem.hpp>
#include <nlohmann/json.hpp>

#if defined(MFEM_USE_SUPERLU)
#include <superlu_defs.h>
#endif
#if defined(MFEM_USE_STRUMPACK)
#include <StrumpackConfig.hpp>
#endif
#if defined(MFEM_USE_SUNDIALS)
#include <sundials/sundials_version.h>
#endif
#if defined(PALACE_WITH_SLEPC)
#include <petscversion.h>
#include <slepcversion.h>
#endif
#if defined(MFEM_USE_CUDA)
#include <cuda_runtime.h>
#elif defined(MFEM_USE_HIP)
#include <hip/hip_runtime.h>
#endif
#if defined(PALACE_WITH_UMPIRE)
#include <camp/config.hpp>
#include <umpire/Umpire.hpp>
#endif
#if defined(PALACE_WITH_CUDSS)
#include <cudss.h>
#endif
#if defined(PALACE_WITH_ZLIB)
#include <zlib.h>
#endif
#if defined(PALACE_WITH_ZFP)
#include <zfp.h>
#endif

// SLATE, BLAS++, and LAPACK++ version queries, declared here to avoid pulling in their full
// headers (which require e.g. OpenMP).
#if defined(PALACE_WITH_SLATE)
namespace slate
{
int version();
}  // namespace slate
#endif
#if defined(PALACE_WITH_BLASPP)
namespace blas
{
int blaspp_version();
}  // namespace blas
#endif
#if defined(PALACE_WITH_LAPACKPP)
namespace lapack
{
int lapackpp_version();
}  // namespace lapack
#endif

// LAPACK and vendor BLAS version queries, declared here since their headers are
// vendor-specific. Integer arguments are 64-bit and zero-initialized so they read correctly
// from both LP64 and ILP64 libraries (little-endian).
extern "C"
{
  void ilaver_(std::int64_t *major, std::int64_t *minor, std::int64_t *patch);
#if defined(PALACE_BLAS_OPENBLAS)
  char *openblas_get_config();
#elif defined(PALACE_BLAS_MKL)
  void MKL_Get_Version_String(char *buffer, int len);
#elif defined(PALACE_BLAS_ARMPL)
  void armplversion(std::int64_t *major, std::int64_t *minor, std::int64_t *patch,
                    const char **tag);
#elif defined(PALACE_BLAS_AOCL)
  const char *bli_info_get_version_str();
#endif
}

namespace palace
{

namespace
{

// Decode an integer version encoded as major * major_scale + minor * minor_scale + patch.
std::string DecodeVersion(long long version, long long major_scale, long long minor_scale)
{
  return fmt::format("{}.{}.{}", version / major_scale,
                     (version % major_scale) / minor_scale, version % minor_scale);
}

// Format a date-based version encoded as yyyymmrr (SLATE, BLAS++, LAPACK++).
[[maybe_unused]] std::string DateVersion(int version)
{
  return fmt::format("{:04}.{:02}.{:02}", version / 10000, (version / 100) % 100,
                     version % 100);
}

}  // namespace

std::vector<std::pair<std::string, std::string>> GetDependencyVersions()
{
  std::vector<std::pair<std::string, std::string>> versions;

  char mpi_version[MPI_MAX_LIBRARY_VERSION_STRING];
  int mpi_len;
  MPI_Get_library_version(mpi_version, &mpi_len);
  std::string mpi(mpi_version, mpi_len);
  versions.emplace_back("MPI", mpi.substr(0, mpi.find_first_of("\r\n")));

  versions.emplace_back("MFEM", mfem::GetVersionStr());

  int major, minor, patch;
  bool release;
  CeedGetVersion(&major, &minor, &patch, &release);
  versions.emplace_back(
      "libCEED", fmt::format("{}.{}.{}{}", major, minor, patch, release ? "" : "-dev"));

  versions.emplace_back("HYPRE", HYPRE_RELEASE_VERSION);

  std::int64_t lapack_major = 0, lapack_minor = 0, lapack_patch = 0;
  ilaver_(&lapack_major, &lapack_minor, &lapack_patch);
  versions.emplace_back("LAPACK API",
                        fmt::format("{}.{}.{}", lapack_major, lapack_minor, lapack_patch));
#if defined(PALACE_BLAS_OPENBLAS)
  versions.emplace_back("OpenBLAS", openblas_get_config());
#elif defined(PALACE_BLAS_MKL)
  char mkl_version[256];
  MKL_Get_Version_String(mkl_version, sizeof(mkl_version));
  versions.emplace_back("Intel MKL", mkl_version);
#elif defined(PALACE_BLAS_ARMPL)
  std::int64_t armpl_major = 0, armpl_minor = 0, armpl_patch = 0;
  const char *armpl_tag = nullptr;
  armplversion(&armpl_major, &armpl_minor, &armpl_patch, &armpl_tag);
  versions.emplace_back("Arm PL",
                        fmt::format("{}.{}.{}", armpl_major, armpl_minor, armpl_patch));
#elif defined(PALACE_BLAS_AOCL)
  versions.emplace_back("AOCL-BLAS", bli_info_get_version_str());
#endif

#if defined(PALACE_WITH_SLEPC)
  versions.emplace_back("SLEPc", fmt::format("{}.{}.{}{}", SLEPC_VERSION_MAJOR,
                                             SLEPC_VERSION_MINOR, SLEPC_VERSION_SUBMINOR,
                                             SLEPC_VERSION_RELEASE ? "" : "-dev"));
  versions.emplace_back("PETSc", fmt::format("{}.{}.{}{}", PETSC_VERSION_MAJOR,
                                             PETSC_VERSION_MINOR, PETSC_VERSION_SUBMINOR,
                                             PETSC_VERSION_RELEASE ? "" : "-dev"));
#endif
#if defined(PALACE_WITH_ARPACK)
#if defined(PALACE_ARPACK_VERSION)
  versions.emplace_back("ARPACK-NG", PALACE_ARPACK_VERSION);
#else
  versions.emplace_back("ARPACK-NG", "unknown");
#endif
#endif
#if defined(MFEM_USE_SUPERLU)
  superlu_dist_GetVersionNumber(&major, &minor, &patch);
  versions.emplace_back("SuperLU_DIST", fmt::format("{}.{}.{}", major, minor, patch));
#endif
#if defined(MFEM_USE_STRUMPACK)
  versions.emplace_back("STRUMPACK",
                        fmt::format("{}.{}.{}", STRUMPACK_VERSION_MAJOR,
                                    STRUMPACK_VERSION_MINOR, STRUMPACK_VERSION_PATCH));
#endif
#if defined(PALACE_BUTTERFLYPACK_VERSION)
  versions.emplace_back("ButterflyPACK", PALACE_BUTTERFLYPACK_VERSION);
#endif
#if defined(PALACE_WITH_ZFP)
  versions.emplace_back("ZFP", zfp_version_string);
#endif
#if defined(PALACE_WITH_SLATE)
  versions.emplace_back("SLATE", DateVersion(slate::version()));
#endif
#if defined(PALACE_WITH_BLASPP)
  versions.emplace_back("BLAS++", DateVersion(blas::blaspp_version()));
#endif
#if defined(PALACE_WITH_LAPACKPP)
  versions.emplace_back("LAPACK++", DateVersion(lapack::lapackpp_version()));
#endif
#if defined(MFEM_USE_MUMPS)
  versions.emplace_back("MUMPS", MUMPS_VERSION);
#endif
#if defined(MFEM_USE_SUNDIALS)
  char sundials_version[64];
  SUNDIALSGetVersion(sundials_version, sizeof(sundials_version));
  versions.emplace_back("SUNDIALS", sundials_version);
#endif
#if defined(PALACE_WITH_CUDSS)
  major = minor = patch = 0;
  cudssGetProperty(MAJOR_VERSION, &major);
  cudssGetProperty(MINOR_VERSION, &minor);
  cudssGetProperty(PATCH_LEVEL, &patch);
  versions.emplace_back("cuDSS", fmt::format("{}.{}.{}", major, minor, patch));
#endif
#if defined(PALACE_GSLIB_VERSION)
  versions.emplace_back("GSLIB", PALACE_GSLIB_VERSION);
#endif
  versions.emplace_back("METIS", fmt::format("{}.{}.{}", METIS_VER_MAJOR, METIS_VER_MINOR,
                                             METIS_VER_SUBMINOR));
  versions.emplace_back("ParMETIS",
                        fmt::format("{}.{}.{}", PARMETIS_MAJOR_VERSION,
                                    PARMETIS_MINOR_VERSION, PARMETIS_SUBMINOR_VERSION));

#if defined(MFEM_USE_CUDA)
  int cuda_version = 0;
  if (cudaRuntimeGetVersion(&cuda_version) == cudaSuccess)
  {
    versions.emplace_back("CUDA runtime", fmt::format("{}.{}", cuda_version / 1000,
                                                      (cuda_version % 1000) / 10));
  }
#elif defined(MFEM_USE_HIP)
  int hip_version = 0;
  if (hipRuntimeGetVersion(&hip_version) == hipSuccess)
  {
    versions.emplace_back("HIP runtime", DecodeVersion(hip_version, 10000000, 100000));
  }
#endif
#if defined(PALACE_WITH_UMPIRE)
  versions.emplace_back("Umpire", fmt::format("{}.{}.{}", umpire::get_major_version(),
                                              umpire::get_minor_version(),
                                              umpire::get_patch_version()));
  versions.emplace_back("camp", fmt::format("{}.{}.{}", CAMP_VERSION_MAJOR,
                                            CAMP_VERSION_MINOR, CAMP_VERSION_PATCH));
#endif
#if defined(PALACE_LIBXSMM_VERSION)
  versions.emplace_back("libxsmm", PALACE_LIBXSMM_VERSION);
#endif
#if defined(PALACE_MAGMA_VERSION)
  versions.emplace_back("MAGMA", PALACE_MAGMA_VERSION);
#endif
#if defined(PALACE_WITH_ZLIB)
  versions.emplace_back("zlib", zlibVersion());
#endif

#if defined(EIGEN_VERSION_STRING)
  versions.emplace_back("Eigen", EIGEN_VERSION_STRING);
#else
  versions.emplace_back("Eigen", fmt::format("{}.{}.{}", EIGEN_WORLD_VERSION,
                                             EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION));
#endif
  versions.emplace_back("fmt", DecodeVersion(FMT_VERSION, 10000, 100));
  versions.emplace_back("scn", DecodeVersion(SCN_VERSION, 10000000, 10000));
#if defined(PALACE_JSON_SCHEMA_VALIDATOR_VERSION)
  versions.emplace_back("nlohmann/json-schema-validator",
                        PALACE_JSON_SCHEMA_VALIDATOR_VERSION);
#endif
  versions.emplace_back("nlohmann/json",
                        fmt::format("{}.{}.{}", NLOHMANN_JSON_VERSION_MAJOR,
                                    NLOHMANN_JSON_VERSION_MINOR,
                                    NLOHMANN_JSON_VERSION_PATCH));

  return versions;
}

}  // namespace palace
