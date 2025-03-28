#!/bin/bash

if [[ ! -f ${PWD}/spack/scripts/$(basename $0) ]]; then
  echo "Please run this script from the palace root dir."
  exit 1
fi

if [[ -z ${SPACK_ROOT} ]]; then
  echo "You have not configured a system spack!"
  echo "Please set SPACK_ROOT and re-run to continue"
  echo "The script will clone spack to SPACK ROOT if no spack present."
  exit 1
fi

# Clone spack if it doesn't exist
if [ ! -f ${SPACK_ROOT}/share/spack/setup-env.sh ]; then
  mkdir -p ${SPACK_ROOT}
  git clone -c feature.manyFiles=true --depth=2 https://github.com/spack/spack.git $SPACK_ROOT
fi

source ${SPACK_ROOT}/share/spack/setup-env.sh || exit 1

SPACK_ENV=${PWD}/$(spack arch --generic)

# Can easily change compiler / MPI / blas spec
# NOTES:
#   - intel-oneapi-mkl only works on Linux / x86_64
# PALACE_SPEC="local.palace@develop ^intel-oneapi-mkl ^openmpi"
PALACE_SPEC="local.palace@develop ^openblas ^openmpi"

FRESH_INSTALL=false

while [[ $# -gt 0 ]]; do
  case "$1" in
  -f)
    echo "-f set. Doing a full env rebuild"
    FRESH_INSTALL=true
    shift
    ;;
  *)
    echo "Error: invalid option: $1"
    return 1
    ;;
  esac
done

# Prevents loading ~/.spack
export SPACK_DISABLE_LOCAL_CONFIG=1

# Configured temporary directories for building
# Will use `/dev/shm` for now, but might need more space...
export TMP=/tmp/spack-tmp
export TMPDIR=${TMP}
export tempdir=${TMP}
mkdir -p $TMP

# Assumes that Python is installed / present on the system
SPACK_PYTHON=$(which python3)

# Should do a fresh install if no spack.yaml in env
if [ ! -f ${SPACK_ENV}/spack.yaml ]; then
  FRESH_INSTALL=true
fi

if [ ${FRESH_INSTALL} = "true" ]; then
  if [ -d ${SPACK_ENV} ]; then
    rm -rfd ${SPACK_ENV}
  fi
  spack env create -d ${SPACK_ENV}

  # We don't need to clean every time, but might as well to avoid issues...
  # spack -e ${SPACK_ENV} clean -m

  # Add our repo-local definition of palace spack package
  spack -e ${SPACK_ENV} repo add spack/local

  # Add our package
  spack -e ${SPACK_ENV} add ${PALACE_SPEC}
  spack -e ${SPACK_ENV} develop --path=${PWD} ${PALACE_SPEC}

  # Configure externals / concretizer
  spack -e ${SPACK_ENV} compiler find
  spack -e ${SPACK_ENV} external find

  __IS_MAC=${__IS_MAC:-$(test $(uname -s) == "Darwin" && echo 'true')}

  if [ -n "${__IS_MAC}" ]; then
    # Assumes that you have an openblas / openmpi installation you want to use
    # Install with brew if you would like to use this
    PACKAGES=(
      "openblas"
      "openmpi"
      "python"
      "cmake"
      "curl"
    )
    for PKG in "${PACKAGES[@]}"; do
      PKG_PATH=$(brew --prefix ${PKG})
      spack -e ${SPACK_ENV} external find --path ${PKG_PATH} ${PKG}
      # We could also add specification of virtual providers here...
      spack -e ${SPACK_ENV} config add packages:${PKG}:buildable:false
    done
  fi

  # Views are only for environments being used with many core specs
  spack -e ${SPACK_ENV} config add view:false
  # We want to re-use externals/buildcache
  # However, the concretization kept giving new-deps, so disabling...
  # Maybe reuse:dependencies would be better
  spack -e ${SPACK_ENV} config add concretizer:reuse:false
  # Unify doesn't really apply until you have more than one core spec
  spack -e ${SPACK_ENV} config add concretizer:unify:true
  # Allow splicing of MPI / Compilers from build cache for faster builds
  # While splicing does speed up builds, it often creates wild DAGs
  spack -e ${SPACK_ENV} config add concretizer:splice:automatic:false
  # Try to prevent duplcate packages in concretization
  spack -e ${SPACK_ENV} config add concretizer:duplicates:strategy:none
  # Target the most generic microarch to encourage re-use
  spack -e ${SPACK_ENV} config add concretizer:targets:granularity:generic
  # Lets make sure that spack uses our tmp directories
  spack -e ${SPACK_ENV} config add config:source_cache:${TMP}
  spack -e ${SPACK_ENV} config add config:misc_cache:${TMP}
  spack -e ${SPACK_ENV} config add config:build_stage:${TMP}
  # We have some requirements for other package variants
  spack -e ${SPACK_ENV} config add packages:petsc:require:~hdf5

  # Add public mirror to help with build times
  spack -e ${SPACK_ENV} mirror add develop https://binaries.spack.io/develop
  spack -e ${SPACK_ENV} buildcache keys --install --trust

  # Verify concretization is correct before installing
  spack -e ${SPACK_ENV} concretize -f
  echo ""
  echo "NOTE: Even though it doesn't say palace.local, it's using that version."
  echo "If you are happy with this concretization, press Enter to continue..."
  read -r
  echo "Installing..."
  echo ""
fi

spack -e ${SPACK_ENV} install -j $(nproc 2>/dev/null || sysctl -n hw.ncpu) --fail-fast --only-concrete

DEV_PATH=$(spack -e ${SPACK_ENV} location --build-dir ${PALACE_SPEC})

echo
echo "Installation done / cancelled / failed. Feel free to re-run build with:"
echo
echo "    cd ${DEV_PATH}/spack-build"
echo "    make -j"
echo
echo "OR re-run the script without -f as so:"
echo "    ${PWD}/spack/scripts/$(basename $0)"
echo
echo "NOTE: The path in the cd command is also in ./build-$(spack arch)-<hash>"
echo "      which is symlinked to the one output above."
echo
