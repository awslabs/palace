#!/bin/bash

if [[ ! -z $SPACK_COMMAND ]]; then
  echo "..."
  echo
else
  SPACK_COMMAND=spack
fi

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

# Add optional bootstrap of Python if using ubuntu (in a container...)
# This is made irrelevant by running with ubuntu/python instead...
python3 --version >/dev/null 2>&1
if [ ! $? -eq 0 ]; then
  echo "Python not found!"
  OS_VER=$(awk -F= '/^NAME/{print $2}' /etc/os-release)
  echo "OS_VER == " ${OS_VER}
  # The quotes are messy here and I don't like it
  if [ "${OS_VER}" = '"Ubuntu"' ]; then
    apt-get update
    apt-get install -y python3
    apt-get clean
    rm -rfd /var/lib/apt/lists/*
  else
    exit 1
  fi
fi

# Add this optional bootstrap of gcc if using Ubuntu (in a container...)
gcc --version >/dev/null 2>&1
if [ ! $? -eq 0 ]; then
  echo "GCC not found!"
  OS_VER=$(awk -F= '/^NAME/{print $2}' /etc/os-release)
  echo "OS_VER == " ${OS_VER}
  # The quotes are messy here and I don't like it
  if [ "${OS_VER}" = '"Debian GNU/Linux"' ]; then
    apt-get update
    apt-get install -y gcc g++ gfortran
    apt-get clean
    rm -rfd /var/lib/apt/lists/*
  else
    exit 1
  fi

fi

if [[ ! -z $SPACK_COMMAND ]]; then
  if [[ $SPACK_COMMAND == "spack" ]]; then
    SPACK_PYTHON=$(which python3)
    source ${SPACK_ROOT}/share/spack/setup-env.sh
  fi
fi

# A quick check to make sure spack is working...
echo Using spack@$(${SPACK_COMMAND} --version) || exit 1

# Maybe this could be a one-liner...
GARCH_STR="${SPACK_COMMAND} arch --generic"
GARCH="$(${GARCH_STR})"
SPACK_ENV=${PWD}/${GARCH}

# Can easily change compiler / MPI / blas spec
# NOTES:
#   - intel-oneapi-mkl only works on Linux / x86_64
# PALACE_SPEC="local.palace@develop ^intel-oneapi-mkl ^openmpi"
PALACE_SPEC="local.palace@develop ^openblas ^openmpi"

# Somtimes we want to default to true
FRESH_INSTALL=false
# FRESH_INSTALL=true

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

# Should do a fresh install if no spack.yaml in env
if [ ! -f ${SPACK_ENV}/spack.yaml ]; then
  FRESH_INSTALL=true
fi

if [ ${FRESH_INSTALL} = "true" ]; then
  if [ -d ${SPACK_ENV} ]; then
    rm -rfd ${SPACK_ENV}
  fi
  ${SPACK_COMMAND} env create -d ${SPACK_ENV} || exit 1

  # We don't need to clean every time, but might as well to avoid issues...
  ${SPACK_COMMAND} -e ${SPACK_ENV} clean -m

  # Add our repo-local definition of palace spack package
  ${SPACK_COMMAND} -e ${SPACK_ENV} repo add spack/local

  # Add our package
  ${SPACK_COMMAND} -e ${SPACK_ENV} add ${PALACE_SPEC}
  ${SPACK_COMMAND} -e ${SPACK_ENV} develop --path=${PWD} ${PALACE_SPEC}

  # Configure externals / concretizer
  ${SPACK_COMMAND} -e ${SPACK_ENV} compiler find
  ${SPACK_COMMAND} -e ${SPACK_ENV} external find --all

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
      ${SPACK_COMMAND} -e ${SPACK_ENV} external find --path ${PKG_PATH} ${PKG}
      # We could also add specification of virtual providers here...
      ${SPACK_COMMAND} -e ${SPACK_ENV} config add packages:${PKG}:buildable:false
    done
  fi

  # Views are only for environments being used with many core specs
  ${SPACK_COMMAND} -e ${SPACK_ENV} config add view:false
  # We want to re-use externals/buildcache
  # However, the concretization kept giving new-deps, so disabling...
  # Maybe reuse:dependencies would be better
  ${SPACK_COMMAND} -e ${SPACK_ENV} config add concretizer:reuse:false
  # Unify doesn't really apply until you have more than one core spec
  ${SPACK_COMMAND} -e ${SPACK_ENV} config add concretizer:unify:true
  # Allow splicing of MPI / Compilers from build cache for faster builds
  # While splicing does speed up builds, it often creates wild DAGs
  ${SPACK_COMMAND} -e ${SPACK_ENV} config add concretizer:splice:automatic:false
  # Try to prevent duplcate packages in concretization
  ${SPACK_COMMAND} -e ${SPACK_ENV} config add concretizer:duplicates:strategy:none
  # Target the most generic microarch to encourage re-use
  ${SPACK_COMMAND} -e ${SPACK_ENV} config add concretizer:targets:granularity:generic
  # Lets make sure that ${SPACK_COMMAND} uses our tmp directories
  ${SPACK_COMMAND} -e ${SPACK_ENV} config add config:source_cache:${TMP}
  ${SPACK_COMMAND} -e ${SPACK_ENV} config add config:misc_cache:${TMP}
  ${SPACK_COMMAND} -e ${SPACK_ENV} config add config:build_stage:${TMP}
  # We have some requirements for other package variants
  ${SPACK_COMMAND} -e ${SPACK_ENV} config add packages:petsc:require:~hdf5

  # Add public mirror to help with build times
  ${SPACK_COMMAND} -e ${SPACK_ENV} mirror add develop https://binaries.spack.io/develop
  ${SPACK_COMMAND} -e ${SPACK_ENV} buildcache keys --install --trust

  # Verify concretization is correct before installing
  ${SPACK_COMMAND} -e ${SPACK_ENV} concretize -f || exit 1
  echo ""
  echo "NOTE: Even though it doesn't say palace.local, it's using that version."
  echo "If you are happy with this concretization, press Enter to continue..."
  read -r
  echo "Installing..."
  echo ""
fi

${SPACK_COMMAND} -e ${SPACK_ENV} install -j $(nproc 2>/dev/null || sysctl -n hw.ncpu) --fail-fast --only-concrete

DEV_PATH=$(${SPACK_COMMAND} -e ${SPACK_ENV} location --build-dir ${PALACE_SPEC})

echo
echo "Installation done / cancelled / failed. Feel free to re-run build with:"
echo
echo "    cd ${DEV_PATH}/spack-build"
echo "    make -j"
echo
echo "OR re-run the script without -f as so:"
echo "    ${PWD}/spack/scripts/$(basename $0)"
echo
echo "NOTE: The path in the cd command is also in ./build-$(${SPACK_COMMAND} arch)-<hash>"
echo "      which is symlinked to the one output above."
echo
echo "NOTE: Just re-run the script if you are using the container script"
