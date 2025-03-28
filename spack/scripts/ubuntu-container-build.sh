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

UBUNTU_TAG=oracular
FINCH_TEST="finch run --rm ubuntu:${UBUNTU_TAG} echo Finch runs Ubuntu..."
CONTAINER_ROOT=/scratch
CONTAINER_SPACK_ROOT=${CONTAINER_ROOT}/local-clone-of-spack
# Changing tmpdir requires changing local-build.sh as well.
TMPDIR=/tmp/spack-tmp

${FINCH_TEST} || (finch vm start && ${FINCH_TEST})

finch run --rm \
  -v ${PWD}:${CONTAINER_ROOT} \
  -v ${SPACK_ROOT}:${CONTAINER_SPACK_ROOT} \
  -v ${TMPDIR}:${TMPDIR} \
  --env SPACK_ROOT=${CONTAINER_SPACK_ROOT} \
  -w ${CONTAINER_ROOT} \
  ubuntu:${UBUNTU_TAG} \
  ${CONTAINER_ROOT}/spack/scripts/local-build.sh
