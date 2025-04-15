# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import (
    CMakePackage,
    CudaPackage,
    ROCmPackage,
    maintainers,
    version,
    variant,
    depends_on,
    when,
    conflicts,
)


class Palace(CMakePackage, CudaPackage, ROCmPackage):
    """3D finite element solver for computational electromagnetics"""

    tags = ["cem", "fem", "finite-elements", "hpc", "solver"]

    homepage = "https://github.com/awslabs/palace"
    git = "https://github.com/awslabs/palace.git"

    maintainers("sebastiangrimberg")

    version("develop", branch="main")

    # Note: 'cuda' and 'cuda_arch' variants are added by the CudaPackage
    # Note: 'rocm' and 'amdgpu_target' variants are added by the ROCmPackage
    variant("shared", default=True, description="Build shared libraries")
    variant("int64", default=False, description="Use 64 bit integers")
    variant(
        "openmp", default=False, description="Use OpenMP for shared-memory parallelism"
    )
    variant(
        "superlu-dist",
        default=True,
        description="Build with SuperLU_DIST sparse direct solver",
    )
    variant(
        "strumpack",
        default=False,
        description="Build with STRUMPACK sparse direct solver",
    )
    variant(
        "sundials",
        default=True,
        description="Build with SUNDIALS differential/algebraic equations solver",
    )
    variant("mumps", default=False, description="Build with MUMPS sparse direct solver")
    variant("slepc", default=True, description="Build with SLEPc eigenvalue solver")
    variant("arpack", default=False, description="Build with ARPACK eigenvalue solver")
    variant("magma", default=True, description="Build with MAGMA backend for libCEED")
    variant(
        "gslib",
        default=True,
        description="Build with GSLIB library for high-order field interpolation",
    )

    # TODO: Apply patches for all packages...
    # TODO: We should actually use these as externals...

    # These are our hard Dependencies
    # TODO: Need to specify @git.v4.8-rc0=develop (maybe only in spack.yaml)
    # depends_on("mfem@git.v4.8-rc0=develop")
    depends_on("mfem@develop+metis+zlib~fms~libceed")
    depends_on("metis@5:")
    depends_on("hypre~complex")
    depends_on("gslib+mpi")

    # superlu-dist isn't a hard dep, but some variants are
    depends_on("superlu-dist+parmetis~openmp~cuda~rocm", when="+superlu-dist")

    # LibCEED is a core dep
    # TODO: We need to specify @git.v0.13.0-rc.1=develop (maybe only in spack.yaml)
    depends_on("libceed@develop")
    depends_on("libceed+magma", when="+magma")
    # Spack says that libxsmm isn't available on Darwin...
    # Are there other operating systems that we can add support to (windows)?
    # depends_on("libceed~libxsmm", when="platform=darwin")
    # depends_on("libceed+libxsmm", when="platform=linux")
    # NOTE: @=main != @main since libxsmm has a version main-2023-22
    # depends_on("libxsmm@=main~shared blas=0", when="platform=linux")
    # depends_on("libxsmm@=main~shared blas=0")

    depends_on("cmake@3.21:", type="build")
    depends_on("pkgconfig", type="build")
    depends_on("mpi")
    depends_on("zlib-api")
    depends_on("nlohmann-json")
    depends_on("fmt")
    depends_on("eigen")

    # Conditional base dependencies
    depends_on("slepc", when="+slepc")
    depends_on("strumpack+butterflypack+zfp+parmetis", when="+strumpack")
    depends_on("mumps+metis+parmetis", when="+mumps")
    depends_on("petsc+mpi+double+complex", when="+slepc")
    depends_on("arpack-ng+mpi+icb@develop", when="+arpack")

    # Further propagate variants.
    for pkg in ["mumps", "strumpack", "superlu-dist", "gslib", "sundials"]:
        depends_on(f"mfem+{pkg}", when=f"+{pkg}")

    with when("build_type=Debug"):
        depends_on("mfem+libunwind")
        depends_on("libxsmm+debug")

    # Magma is our GPU backend, so we need it when gpus are enabled
    conflicts("~magma", when="+cuda")
    conflicts("~magma", when="+rocm")
    conflicts(
        "+cuda+rocm", msg="PALACE_WITH_CUDA is not compatible with PALACE_WITH_HIP"
    )

    # Basic constraints of the package
    conflicts("~arpack~slepc", msg="At least one eigenvalue solver is required")
    conflicts(
        "~superlu-dist~strumpack~sundials~mumps",
        msg="Need at least one sparse direct solver",
    )

    # More dependency variant conflicts
    conflicts(
        "^hypre+int64", msg="Palace uses HYPRE's mixedint option for 64 bit integers"
    )
    conflicts("^mumps+int64", msg="Palace requires MUMPS without 64 bit integers")
    conflicts("^slepc+arpack", msg="Palace requires SLEPc without ARPACK")

    # Propogate important variants
    # First element is what we depend on
    # Second is when we depend on it. If no val, always depend on it / no variant controls it
    for pkg in [
        ("metis", ""),
        ("hypre", ""),
        ("strumpack", "+strumpack"),
        ("superlu-dist", "+superlu-dist"),
        ("sundials", "+sundials"),
        ("mumps", "+mumps"),
        ("petsc", "+slepc"),  # Need PETSc when we use slepc
        ("arpack-ng", "+arpack"),
        ("magma", "+magma"),
        ("mfem", ""),
        ("gslib", ""),
    ]:
        depends_on(f"{pkg[0]}+shared", when=f"{pkg[1]}+shared")
        depends_on(f"{pkg[0]}~shared", when=f"{pkg[1]}~shared")

        # For complex / int64
        if pkg[0] in ["metis", "superlu-dist", "petsc"]:
            depends_on(f"{pkg[0]}+int64", when=f"{pkg[1]}+int64")
            depends_on(f"{pkg[0]}~int64", when=f"{pkg[1]}~int64")
        # Hypre is special
        elif pkg[0] == "hypre~complex":
            depends_on(f"{pkg[0]}+mixedint", when=f"{pkg[1]}+int64")
            depends_on(f"{pkg[0]}~mixedint", when=f"{pkg[1]}~int64")

        # OpenMP
        if pkg[0] in [
            "hypre",
            "strumpack",
            "sundials",
            "mumps",
            "petsc",
            "mfem",
        ]:
            depends_on(f"{pkg[0]}+openmp", when=f"{pkg[1]}+openmp")
            depends_on(f"{pkg[0]}~openmp", when=f"{pkg[1]}~openmp")

    # Now for GPU targets
    conflicts(
        "cuda_arch=none",
        when="+cuda",
        msg="palace: Please specify a CUDA arch value / values",
    )
    conflicts(
        "amdgpu_target=none",
        when="+rocm",
        msg="palace: Please specify an AMD GPU target / targets",
    )

    # Magma is at the core of our GPU backend, so that's our ~/+gpu variant...
    with when("+magma"):
        for gpu_pkg in [
            ("hypre", ""),
            ("strumpack", "+strumpack"),
            ("sundials", "+sundials"),
            ("slepc", "+slepc"),
            ("petsc", "+slepc"),  # Need PETSc when we use slepc
            ("magma", "+magma"),
            ("mfem", ""),
            ("libceed", ""),
        ]:
            with when("+cuda"):
                for arch in CudaPackage.cuda_arch_values:
                    cuda_variant = f"+cuda cuda_arch={arch}"
                    depends_on(
                        f"{gpu_pkg[0]}{cuda_variant}",
                        when=f"{gpu_pkg[1]}{cuda_variant}",
                    )
            with when("+rocm"):
                for arch in ROCmPackage.amdgpu_targets:
                    rocm_variant = f"+rocm amdgpu_target={arch}"
                    depends_on(
                        f"{gpu_pkg[0]}{rocm_variant}",
                        when=f"{gpu_pkg[1]}{rocm_variant}",
                    )

    def cmake_args(self):
        args = [
            self.define_from_variant("BUILD_SHARED_LIBS", "shared"),
            self.define_from_variant("PALACE_WITH_64BIT_INT", "int64"),
            self.define_from_variant("PALACE_WITH_OPENMP", "openmp"),
            self.define_from_variant("PALACE_WITH_SUPERLU", "superlu-dist"),
            self.define_from_variant("PALACE_WITH_STRUMPACK", "strumpack"),
            self.define_from_variant("PALACE_WITH_SUNDIALS", "sundials"),
            self.define_from_variant("PALACE_WITH_MUMPS", "mumps"),
            self.define_from_variant("PALACE_WITH_SLEPC", "slepc"),
            self.define_from_variant("PALACE_WITH_ARPACK", "arpack"),
            self.define("PALACE_WITH_LIBXSMM", True),
            self.define_from_variant("PALACE_WITH_MAGMA", "magma"),
            self.define_from_variant("PALACE_WITH_GSLIB", "gslib"),
            self.define("libCEED_DIR", self.spec["libceed"].prefix),
            self.define("PALACE_BUILD_EXTERNAL_DEPS", False),
            self.define_from_variant("PALACE_WITH_CUDA", "cuda"),
            self.define_from_variant("PALACE_WITH_HIP", "rocm"),
            # This is an experimental flag while we transition our meta-build
            self.define("PALACE_WITH_SPACK", True),
        ]

        # We guarantee that there are arch specs with conflicts above
        if "+cuda" in self.spec:
            args.append(
                self.define(
                    "CMAKE_CUDA_ARCHITECTURES", self.spec.variants["cuda_arch"].value
                )
            )
        if "+rocm" in self.spec:
            args.append(
                self.define(
                    "CMAKE_HIP_ARCHITECTURES", self.spec.variants["amdgpu_target"].value
                )
            )

        # HYPRE is always built with external BLAS/LAPACK
        args.extend(
            [
                self.define("HYPRE_REQUIRED_PACKAGES", "LAPACK;BLAS"),
                self.define("BLAS_LIBRARIES", self.spec["blas"].libs),  # type: ignore
                self.define("LAPACK_LIBRARIES", self.spec["lapack"].libs),  # type: ignore
            ]
        )

        # MPI compiler wrappers are not required, but MFEM test builds need to know to link
        # against MPI libraries.
        # Eventually these will all be external spack based dependencies.
        if "+superlu-dist" in self.spec:
            args.append(self.define("SuperLUDist_REQUIRED_PACKAGES", "LAPACK;BLAS;MPI"))
        if "+sundials" in self.spec:
            args.append(self.define("SUNDIALS_REQUIRED_PACKAGES", "LAPACK;BLAS;MPI"))
        if "+strumpack" in self.spec:
            args.append(
                self.define(
                    "STRUMPACK_REQUIRED_PACKAGES", "LAPACK;BLAS;MPI;MPI_Fortran"
                )
            )
        if "+mumps" in self.spec:
            args.append(
                self.define("MUMPS_REQUIRED_PACKAGES", "LAPACK;BLAS;MPI;MPI_Fortran")
            )

        # Configure libCEED build
        if "+libxsmm" in self.spec:
            args.append(self.define("LIBXSMM_DIR", self.spec["libxsmm"].prefix))  # type: ignore
        if "+magma" in self.spec:
            args.append(self.define("MAGMA_DIR", self.spec["magma"].prefix))  # type: ignore

        return args

    def install(self, spec, prefix):
        # No install phase for Palace (always performed during build)
        pass
