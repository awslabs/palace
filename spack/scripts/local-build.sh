#!/bin/bash

set -e

if [[ ! -z $SPACK_COMMAND ]]; then
  echo "..."
  echo
else
  # SPACK_COMMAND="spack --debug"
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

if [[ $SPACK_COMMAND == "spack" ]]; then
  SPACK_PYTHON=$(which python3)
  source ${SPACK_ROOT}/share/spack/setup-env.sh
fi

# A quick check to make sure spack is working...
# echo Using spack@$(${SPACK_COMMAND} --version) || exit 1

# Maybe this could be a one-liner...
GARCH_STR="${SPACK_COMMAND} arch --generic"
GARCH="$(${GARCH_STR})"
if [[ $SPACK_COMMAND == "spack" ]]; then
  export SPACK_ENV=${PWD}/${GARCH}
else
  export SPACK_ENV=${CONTAINER_ROOT}/${GARCH}
fi

echo $SPACK_ENV

FRESH_INSTALL=false

while [[ $# -gt 0 ]]; do
  case "$1" in
  -f)
    echo "-f set. Doing a full env rebuild"
    FRESH_INSTALL=true
    shift
    ;;
  -ff)
    echo "-ff set. Running additional spack clean"
    FRESH_INSTALL=true
    FORCE_FRESH_INSTALL=true
    shift
    ;;
  -cuda)
    echo "-cuda set. Running with cuda build"
    CUDA=true
    shift
    ;;
  *)
    echo "Error: invalid option: $1"
    return 1
    ;;
  esac
done

# Can easily change compiler / MPI / blas spec
# NOTES:
#   - intel-oneapi-mkl only works on Linux / x86_64
BLAS_SPEC="openblas"

if [ ${CUDA} = "true" ]; then
  GPU="+cuda cuda_arch=89"
else
  GPU=""
fi
PALACE_SPEC="local.palace@develop${GPU} ^${BLAS_SPEC} ^openmpi"

# Prevents loading ~/.spack
export SPACK_DISABLE_LOCAL_CONFIG=1

# Should do a fresh install if no spack.yaml in env
if [ ! -f ${PWD}/${GARCH}/spack.yaml ]; then
  FRESH_INSTALL=true
fi

if [ ${FRESH_INSTALL} = "true" ]; then
  if [[ "${SPACK_COMMAND}" == "spack" ]]; then
    if [ -d ${SPACK_ENV} ]; then
      rm -rfd ${SPACK_ENV}
    fi
    mkdir -p ${SPACK_ENV}
    TMP_SPACK_ENV=${SPACK_ENV}
  else
    TMP_SPACK_ENV=${PWD}/${GARCH}
    if [ -d ${TMP_SPACK_ENV} ]; then
      rm -rfd ${TMP_SPACK_ENV}
    fi
    mkdir -p ${TMP_SPACK_ENV}
  fi

  # To enable / disable using a build cache aggressively, toggle:
  #   - reuse: true
  #   - splce: automatic: true
  # TODO: What versions of Palace support corresponding MFEM versions?
  cat <<EOF >${TMP_SPACK_ENV}/spack.yaml
  spack:
    specs: 
      - ${PALACE_SPEC}
      - local.gslib+shared
      - local.libceed
    repos:
    - ${SPACK_ENV}/../spack/local
    develop:
      palace:
        spec: ${PALACE_SPEC}
        path: ${SPACK_ENV}/..
    view: false
    concretizer:
      reuse: false
      unify: true
      splice:
        automatic: false
      duplicates:
        strategy: none
      targets:
        granularity: generic
    packages:
      petsc:
        require: ~hdf5
      rocblas:
        require: ~tensile
      openmpi:
        require: ~cuda
    mirrors:
        develop: https://binaries.spack.io/develop
EOF

  if [[ ${FORCE_FRESH_INSTALL} = "true" ]]; then
    ${SPACK_COMMAND} clean -abm
    rm -rfd ~/.spack*
  fi

  # We don't need to clean every time, but might as well to avoid issues...
  # ${SPACK_COMMAND} -e ${SPACK_ENV} clean -abm
  # ${SPACK_COMMAND} -e ${SPACK_ENV} clean -m
  # *** The next one removes everything spack has installed ***
  # ${SPACK_COMMAND} -e ${SPACK_ENV} gc -by

  # Configure externals / compiler
  ${SPACK_COMMAND} -e ${SPACK_ENV} external find --all \
    --exclude openssl \
    --exclude openssh \
    --exclude python \
    --exclude ncurses \
    --exclude bzip2 \
    --exclude xz \
    --exclude curl

  # MAC OS(?)
  # ${SPACK_COMMAND} -e ${SPACK_ENV} external find --all --exclude curl
  if [[ "${SPACK_COMMAND}" == "spack" ]]; then
    # Assumes that you have an openblas / openmpi installation you want to use
    # Install with brew if you would like to use this
    if command -v brew 2>&1 >/dev/null; then
      PACKAGES=(
        "openblas"
        "curl"
      )
      for PKG in "${PACKAGES[@]}"; do
        PKG_PATH=$(brew --prefix ${PKG})
        ${SPACK_COMMAND} -e ${SPACK_ENV} external find --path ${PKG_PATH} ${PKG}
        # We could also add specification of virtual providers here...
        ${SPACK_COMMAND} -e ${SPACK_ENV} config add packages:${PKG}:buildable:false
      done
    fi
  fi

  # Add public mirror to help with build times
  if [[ ${FORCE_FRESH_INSTALL} = "true" ]]; then
    ${SPACK_COMMAND} -e ${SPACK_ENV} buildcache keys --install --trust
  fi

  # TODO: FIX
  # if [[ ! "${SPACK_COMMAND}" == "spack" ]]; then
  #   # ${SPACK_COMMAND} -e ${SPACK_ENV} mirror add \
  #   #   --oci-password-variable GITHUB_PAT \
  #   #   --oci-username cameronrutherford \
  #   #   github oci://ghcr.io/awslabs/palace
  # fi

  # TODO: Add config to add GitHub cache
  # TODO: Configure secrets in a secure way...
  # TODO: We will have many caches...

  # Verify concretization is correct before installing
  ${SPACK_COMMAND} -e ${SPACK_ENV} concretize -f || exit 1
  echo ""
  echo "NOTE: Even though it doesn't say palace.local, it's using that version."
  echo "If you are happy with this concretization, press Enter to continue..."
  read -r
fi

echo
echo "Installing Palace..."
echo

${SPACK_COMMAND} -e ${SPACK_ENV} install --fail-fast --only-concrete --keep-stage --only dependencies --show-log-on-error
${SPACK_COMMAND} -e ${SPACK_ENV} install --fail-fast --only-concrete --keep-stage --verbose --show-log-on-error
DEV_PATH=$(${SPACK_COMMAND} -e ${SPACK_ENV} -j $(nproc 2>/dev/null || sysctl -n hw.ncpu) location --build-dir ${PALACE_SPEC})

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

${SPACK_COMMAND} -e ${SPACK_ENV} gc -byE
