# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage
from spack_repo.builtin.build_systems.cuda import CudaPackage
from spack_repo.builtin.build_systems.rocm import ROCmPackage

from spack.package import *


class Palace(CMakePackage, CudaPackage, ROCmPackage):
    """3D finite element solver for computational electromagnetics"""

    tags = ["cem", "fem", "finite-elements", "hpc", "solver"]

    homepage = "https://github.com/awslabs/palace"
    git = "https://github.com/awslabs/palace.git"
    license("Apache-2.0")

    maintainers("hughcars", "simlap", "cameronrutherford", "sbozzolo")

    version("develop", branch="main")
    version("0.15.0", tag="v0.15.0", commit="b6762777d85a06072fdf4cc96e8a365da73df170")
    version("0.14.0", tag="v0.14.0", commit="a428a3a32dbbd6a2a6013b3b577016c3e9425abc")
    version("0.13.0", tag="v0.13.0", commit="a61c8cbe0cacf496cde3c62e93085fae0d6299ac")
    version("0.12.0", tag="v0.12.0", commit="8c192071206466638d5818048ee712e1fada386f")
    version("0.11.2", tag="v0.11.2", commit="6c3aa5f84a934a6ddd58022b2945a1bdb5fa329d")

    # Note: 'cuda' and 'cuda_arch' variants are added by the CudaPackage
    # Note: 'rocm' and 'amdgpu_target' variants are added by the ROCmPackage
    variant("shared", default=True, description="Build shared libraries")
    variant("int64", default=False, description="Use 64 bit integers")
    variant("openmp", default=False, description="Use OpenMP for shared-memory parallelism")
    variant(
        "superlu-dist", default=True, description="Build with SuperLU_DIST sparse direct solver"
    )
    variant("strumpack", default=False, description="Build with STRUMPACK sparse direct solver")
    variant(
        "sundials",
        default=True,
        description="Build with SUNDIALS differential/algebraic equations solver",
        when="@0.14:",
    )
    variant("mumps", default=False, description="Build with MUMPS sparse direct solver")
    variant("slepc", default=True, description="Build with SLEPc eigenvalue solver")
    variant("arpack", default=False, description="Build with ARPACK eigenvalue solver")
    variant("libxsmm", default=True, description="Build with libxsmm backend for libCEED")
    variant(
        "gslib",
        default=True,
        description="Build with GSLIB library for high-order field interpolation",
    )

    # Fix API mismatch between libxsmm@main and internal libceed build
    patch("palace-0.12.0.patch", when="@0.12")

    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("cmake@3.21:", type="build", when="@0.14:0.15")
    depends_on("cmake@3.24:", type="build", when="@0.16:")
    depends_on("pkgconfig", type="build")
    depends_on("mpi")
    depends_on("blas")
    depends_on("lapack")
    depends_on("zlib-api")
    depends_on("nlohmann-json")
    depends_on("fmt+shared", when="+shared")
    depends_on("fmt~shared", when="~shared")
    depends_on("scnlib+shared", when="+shared@0.14:")
    depends_on("scnlib~shared", when="~shared@0.14:")
    depends_on("eigen")

    conflicts("~superlu-dist~strumpack~mumps", msg="Need at least one sparse direct solver")

    conflicts("^mumps+int64", msg="Palace requires MUMPS without 64 bit integers")
    with when("+mumps"):
        depends_on("fortran", type="build")
        depends_on("mumps+metis+parmetis")
        depends_on("mumps+shared", when="+shared")
        depends_on("mumps~shared", when="~shared")
        depends_on("mumps+openmp", when="+openmp")
        depends_on("mumps~openmp", when="~openmp")

    with when("+superlu-dist"):
        depends_on("superlu-dist+parmetis")
        depends_on("superlu-dist+shared", when="+shared")
        depends_on("superlu-dist~shared", when="~shared")
        depends_on("superlu-dist+int64", when="+int64")
        depends_on("superlu-dist~int64", when="~int64")
        depends_on("superlu-dist+openmp", when="+openmp")
        depends_on("superlu-dist~openmp", when="~openmp")

    with when("+strumpack"):
        depends_on("fortran", type="build")
        depends_on("strumpack+butterflypack+zfp+parmetis")
        depends_on("strumpack+shared", when="+shared")
        depends_on("strumpack~shared", when="~shared")
        depends_on("strumpack+openmp", when="+openmp")
        depends_on("strumpack~openmp", when="~openmp")

    conflicts("~arpack~slepc", msg="At least one eigenvalue solver is required")

    with when("+slepc"):
        depends_on("slepc~arpack")
        depends_on("petsc+mpi+double+complex")
        depends_on("petsc+shared", when="+shared")
        depends_on("petsc~shared", when="~shared")
        depends_on("petsc+int64", when="+int64")
        depends_on("petsc~int64", when="~int64")
        depends_on("petsc+openmp", when="+openmp")
        depends_on("petsc~openmp", when="~openmp")

    with when("+arpack"):
        depends_on("fortran", type="build")
        depends_on("arpack-ng+mpi+icb")
        depends_on("arpack-ng+shared", when="+shared")
        depends_on("arpack-ng~shared", when="~shared")

    with when("+gslib @0.14:"):
        depends_on("gslib+mpi")
        depends_on("gslib+shared", when="+shared")
        depends_on("gslib~shared", when="~shared")

    depends_on("metis@5:")
    depends_on("metis+shared", when="+shared")
    depends_on("metis~shared", when="~shared")
    depends_on("metis+int64", when="+int64")
    depends_on("metis~int64", when="~int64")

    conflicts("^hypre+int64", msg="Palace uses HYPRE's mixedint option for 64 bit integers")
    depends_on("hypre@:2", when="@:0.15.0")
    depends_on("hypre@3:")
    depends_on("hypre~complex")
    depends_on("hypre~unified-memory")
    depends_on("hypre+shared", when="+shared")
    depends_on("hypre~shared", when="~shared")
    depends_on("hypre+mixedint", when="+int64")
    depends_on("hypre~mixedint", when="~int64")
    depends_on("hypre+openmp", when="+openmp")
    depends_on("hypre~openmp", when="~openmp")
    depends_on("hypre+gpu-aware-mpi", when="^mpi+cuda")
    # Use external blas/lapack with hypre
    depends_on("hypre+lapack")

    # NOTE: hypre+gpu-profiling is also useful: it adds NVTX annotations, which
    # are great for GPU profiling with Nsight.

    with when("@0.16:"):
        # +lapack means: use external lapack
        depends_on(
            "mfem+mpi+metis+lapack cxxstd=17 commit=bf1b25d82e232909e1d3aa1f416dbdae9931a550",
            patches=[
                "patch_mesh_vis_dev.diff",
                "patch_par_tet_mesh_fix_dev.diff",
                "patch_gmsh_parser_performance.diff",
                "patch_race_condition_fix.diff",
            ],
        )
        depends_on("mfem+shared", when="+shared")
        depends_on("mfem~shared", when="~shared")
        depends_on("mfem+openmp", when="+openmp")
        depends_on("mfem+threadsafe", when="+openmp")
        depends_on("mfem~openmp", when="~openmp")
        depends_on("mfem+superlu-dist", when="+superlu-dist")
        depends_on("mfem~superlu-dist", when="~superlu-dist")
        depends_on("mfem+strumpack", when="+strumpack")
        depends_on("mfem~strumpack", when="~strumpack")
        depends_on("mfem+mumps", when="+mumps")
        depends_on("mfem~mumps", when="~mumps")
        depends_on("mfem+sundials", when="+sundials")
        depends_on("mfem~sundials", when="~sundials")
        depends_on("mfem+gslib", when="+gslib")
        depends_on("mfem~gslib", when="~gslib")
        depends_on("mfem+exceptions", type="test")

        depends_on("mfem+libunwind", when="build_type=Debug")

    with when("+libxsmm"):
        # NOTE: @=main != @main since libxsmm has a version main-2023-22
        depends_on("libxsmm@=main blas=0")
        depends_on("libxsmm+debug", when="build_type=Debug")
        depends_on("libceed+libxsmm", when="@0.14:")
        # NOTE: libxsmm builds on MacOS have linker issues
        # https://github.com/libxsmm/libxsmm/issues/883
        depends_on("libxsmm+shared")

    with when("@0.14:"):
        depends_on("libceed@0.13:")
        depends_on("libceed+openmp", when="+openmp")
        depends_on("libceed~openmp", when="~openmp")
        depends_on("libceed+shared", when="+shared")
        depends_on("libceed~shared", when="~shared")

    with when("+sundials @0.14:"):
        depends_on("sundials+mpi+lapack~examples~examples-install")
        depends_on("sundials+shared", when="+shared")
        depends_on("sundials~shared", when="~shared")
        depends_on("sundials+openmp", when="+openmp")
        depends_on("sundials~openmp", when="~openmp")

    conflicts("+cuda", when="@:0.13", msg="CUDA is only supported for Palace versions after 0.13")
    conflicts("+rocm", when="@:0.13", msg="ROCm is only supported for Palace versions after 0.13")
    conflicts("+cuda+rocm", msg="PALACE_WITH_CUDA is not compatible with PALACE_WITH_HIP")
    conflicts(
        "cuda_arch=none", when="+cuda", msg="palace: Please specify a CUDA arch value / values"
    )
    conflicts(
        "amdgpu_target=none",
        when="+rocm",
        msg="palace: Please specify an AMD GPU target / targets",
    )

    with when("+cuda"):
        depends_on("magma+shared", when="+shared")
        depends_on("magma~shared", when="~shared")
        depends_on("libceed+magma", when="@0.14:")

    with when("+rocm"):
        depends_on("magma+shared", when="+shared")
        depends_on("magma~shared", when="~shared")
        depends_on("libceed+magma", when="@0.14:")

    with when("+cuda"):
        for arch in CudaPackage.cuda_arch_values:
            cuda_variant = f"+cuda cuda_arch={arch}"

            # We need https://github.com/llnl/blt/pull/735, which is not available
            # in <= 0.7.1
            depends_on("blt@0.7.2:")

            depends_on(f"umpire{cuda_variant}", when=f"{cuda_variant}")
            depends_on(f"hypre+umpire{cuda_variant}", when=f"{cuda_variant}")
            depends_on(f"mfem+umpire{cuda_variant}", when=f"{cuda_variant}")
            depends_on(f"magma{cuda_variant}", when=f"{cuda_variant}")
            depends_on(f"libceed{cuda_variant}", when=f"{cuda_variant} @0.14:")
            depends_on(f"sundials{cuda_variant}", when=f"+sundials{cuda_variant} @0.14:")
            depends_on(f"slepc{cuda_variant}", when=f"+slepc{cuda_variant}")
            depends_on(f"petsc{cuda_variant}", when=f"+slepc{cuda_variant}")
            depends_on(f"superlu-dist{cuda_variant}", when=f"+superlu-dist{cuda_variant}")
            depends_on(f"strumpack{cuda_variant}", when=f"+strumpack{cuda_variant}")

    with when("+rocm"):
        for arch in ROCmPackage.amdgpu_targets:
            rocm_variant = f"+rocm amdgpu_target={arch}"
            depends_on(f"umpire{rocm_variant}", when=f"{rocm_variant}")
            depends_on(f"hypre+umpire{rocm_variant}", when=f"{rocm_variant}")
            depends_on(f"mfem+umpire{rocm_variant}", when=f"{rocm_variant}")
            depends_on(f"magma{rocm_variant}", when=f"{rocm_variant}")
            depends_on(f"libceed{rocm_variant}", when=f"{rocm_variant} @0.14:")
            depends_on(f"sundials{rocm_variant}", when=f"+sundials{rocm_variant} @0.14:")
            depends_on(f"slepc{rocm_variant}", when=f"+slepc{rocm_variant}")
            depends_on(f"petsc{rocm_variant}", when=f"+slepc{rocm_variant}")
            depends_on(f"superlu-dist{rocm_variant}", when=f"+superlu-dist{rocm_variant}")
            depends_on(f"strumpack{rocm_variant}", when=f"+strumpack{rocm_variant}")

    depends_on("catch2@3:", type="test")

    def cmake_args(self):
        args = [
            self.define_from_variant("BUILD_SHARED_LIBS", "shared"),
            self.define_from_variant("PALACE_WITH_64BIT_INT", "int64"),
            self.define_from_variant("PALACE_WITH_ARPACK", "arpack"),
            self.define_from_variant("PALACE_WITH_CUDA", "cuda"),
            self.define_from_variant("PALACE_WITH_GSLIB", "gslib"),
            self.define_from_variant("PALACE_WITH_HIP", "rocm"),
            self.define_from_variant("PALACE_WITH_LIBXSMM", "libxsmm"),
            self.define_from_variant("PALACE_WITH_MUMPS", "mumps"),
            self.define_from_variant("PALACE_WITH_OPENMP", "openmp"),
            self.define_from_variant("PALACE_WITH_SLEPC", "slepc"),
            self.define_from_variant("PALACE_WITH_STRUMPACK", "strumpack"),
            self.define_from_variant("PALACE_WITH_SUNDIALS", "sundials"),
            self.define_from_variant("PALACE_WITH_SUPERLU", "superlu-dist"),
            self.define("PALACE_BUILD_EXTERNAL_DEPS", False),
            self.define("PALACE_MFEM_USE_EXCEPTIONS", self.run_tests),
            self.define("PALACE_WITH_GPU_AWARE_MPI", self.spec.satisfies("^mpi+cuda")),
        ]

        if self.spec.satisfies("@0.16:"):
            args.append(self.define("MFEM_DIR", self.spec["mfem"].prefix))
            if self.spec.satisfies("+mumps"):
                args.append(self.define("MUMPS_DIR", self.spec["mumps"].prefix))
            if self.spec.satisfies("+strumpack"):
                args.append(self.define("STRUMPACK_DIR", self.spec["strumpack"].prefix))
            if self.spec.satisfies("+superlu-dist"):
                args.append(self.define("SUPERLU_DIST_DIR", self.spec["superlu-dist"].prefix))
            args.append(self.define("METIS_DIR", self.spec["metis"].prefix))
            args.append(self.define("PARMETIS_DIR", self.spec["parmetis"].prefix))
            args.append(self.define("HYPRE_DIR", self.spec["hypre"].prefix))

        # We guarantee that there are arch specs with conflicts above
        if self.spec.satisfies("+cuda"):
            args.append(
                self.define(
                    "CMAKE_CUDA_ARCHITECTURES", ";".join(self.spec.variants["cuda_arch"].value)
                )
            )
        if self.spec.satisfies("+rocm"):
            args.append(
                self.define(
                    "CMAKE_HIP_ARCHITECTURES", ";".join(self.spec.variants["amdgpu_target"].value)
                )
            )

        # Pass down external BLAS/LAPACK
        args.extend(
            [
                self.define("BLAS_LIBRARIES", self.spec["blas"].libs.joined(";")),
                self.define("LAPACK_LIBRARIES", self.spec["lapack"].libs.joined(";")),
            ]
        )

        if self.spec.satisfies("@:0.13"):
            # In v0.13 and prior libCEED and gslib were internally built and required the libxsmm
            # and magma build information be passed in.
            if self.spec.satisfies("+libxsmm"):
                args.append(self.define("LIBXSMM_DIR", self.spec["libxsmm"].prefix))
            if self.spec.satisfies("+cuda") or self.spec.satisfies("+rocm"):
                args.append(self.define("MAGMA_DIR", self.spec["magma"].prefix))
        else:
            # After v 0.13 gslib and libceed are built externally and
            # so the directories for these are passed explicitly.
            args.append(self.define("LIBCEED_DIR", self.spec["libceed"].prefix))
            if self.spec.satisfies("+gslib"):
                args.append(self.define("GSLIB_DIR", self.spec["gslib"].prefix))

        return args

    def build(self, spec, prefix):
        with working_dir(self.build_directory):
            if self.run_tests:
                make("palace-tests")
            else:
                make()

    def install(self, spec, prefix):
        # No install phase for Palace (always performed during build)
        pass
