#!/bin/bash

# Assumes that you have finch installed locally

if [[ ! -f ${PWD}/spack/scripts/$(basename $0) ]]; then
  echo "Please run this script from the palace root dir."
  exit 1
fi

if [[ -z ${SPACK_ROOT} ]]; then
  echo "You have not configured a system spack!"
  echo "Please set SPACK_ROOT and re-run to continue"
  echo "This is so the container can share your system spack to cache / install."
  exit 1
fi

# This is technically Debian and not ubuntu, but I didn't want to wait
# for Python to install each time while debugging.
CONTAINER=spack/ubuntu-noble
CONTAINER_ROOT=/scratch
CONTAINER_SPACK_ROOT=/opt/local-spack
# Changing tmpdir requires changing local-build.sh as well.
TMPDIR=/tmp/spack-tmp

finch vm start >/dev/null 2>&1 || true

# This is our magic command that:
#   - Runs spack in the ubuntu container
#   - Using our local clone of spack
#   - Storing files in the local spack install
#
# Since this literally only has spack bootstrapped, you can't run a container
# directly. I wanted to iterate quickly, so I went with this.
#
# In the long term, devcontainers are nicer (once they work).
#
# To iterate quickly, and to cache well, these scripts are nice.
finch run --rm \
  -v ${PWD}:${CONTAINER_ROOT} \
  -v ${SPACK_ROOT}:${CONTAINER_SPACK_ROOT} \
  -v ${TMPDIR}:${TMPDIR} \
  --env SPACK_ROOT=${CONTAINER_SPACK_ROOT} \
  -w ${CONTAINER_ROOT} \
  ${CONTAINER} \
  "--version"

SPACK_COMMAND="\
finch run --rm \
  -v ${PWD}:${CONTAINER_ROOT} \
  -v ${SPACK_ROOT}:${CONTAINER_SPACK_ROOT} \
  -v ${TMPDIR}:${TMPDIR} \
  --env SPACK_ROOT=${CONTAINER_SPACK_ROOT} \
  -w ${CONTAINER_ROOT} \
  ${CONTAINER}"

echo Using:
echo
echo $SPACK_COMMAND
echo
echo As the pseudo-spack command for the local build script.
echo

SPACK_COMMAND=$SPACK_COMMAND ./spack/scripts/local-build.sh
