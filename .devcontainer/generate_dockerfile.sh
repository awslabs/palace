#!/bin/bash

set -e

# Configures a Dockerfile for use with our Devcontainer setup
# Assumes that you have a local clone of spack around to generate this

if [ -z "${SPACK_ROOT}" ]; then
  echo "Environment var SPACK_ROOT is unset!"
  echo "Please make a local clone of spack and set SPACK_ROOT before re-running."
  exit 1
fi

if [ ! -f ${PWD}/.devcontainer/generate_dockerfile.sh ]; then
  echo "Please re-run from the Palace root directory!"
  exit 1
fi

export PATH=${PATH}:${SPACK_ROOT}/bin

# This is just for reference in logs
# Note that the Python version might not be the same as you expect...
spack debug report

# Going with a heredoc in this case, just to make the config a little more readable.
# Other options include:
#   - Setting everything through Spack CLI commands
#   - Adding this file explicitly to the repo
cat >${PWD}/.devcontainer/spack.yaml <<EOF
spack:
  specs:
    - local.palace@develop
  concretizer:
    unify: true
    reuse: true
  mirrors:
    develop: https://binaries.spack.io/develop
    palace: oci://ghcr.io/awslabs/palace
  container:
    format: docker
    images:
      os: ubuntu:24.04
      spack:
        ref: v1.0.0-alpha.4
    os_packages:
      build:
        - autoconf
        - gcc-14
        - g++-14
      final:
        - gfortran
  packages:
    petsc:
      require: ~hdf5
    all:
      providers:
        compiler: [gcc]
        mpi: [openmpi]
        blas: [intel-oneapi-mkl]
        lapack: [intel-oneapi-mkl]
EOF

cd .devcontainer
spack containerize >${PWD}/Dockerfile
cd -

# While this auto-generated Dockerfile is great, it's not quite what we want.
# Ideally one builds this container, and then inherits from it to make the devcontainer.
# Since we aren't quite there yet (build issues...), this is the best we can do.

# Configure environment for VSCode User in devcontainer
echo "# Make sure devcontainer user gets spack packages" >>./.devcontainer/Dockerfile
echo "RUN echo \"source /entrypoint.sh\" >> /home/vscode/.bashrc" >>./.devcontainer/Dockerfile

# Make sure vscode user is correct
echo "" >>./.devcontainer/Dockerfile
echo "# Configure user for container" >>./.devcontainer/Dockerfile
echo "USER vscode" >>./.devcontainer/Dockerfile

# This adds some code right before installing to trust / configure source caches
# NOTE: sed commands are notorious for working differently on different OSes
# MAC SED COMMAND vvv
# sed -i "" "s|# Install the software|# Do this separate of install to cache keys...\nRUN --mount=type=secret,id=.secrets,target=.secrets export \$(cat .secrets \| xargs) \&\& cd /opt/spack-environment \&\& spack -e . buildcache keys --install --trust \&\& spack -e . mirror set --oci-username __token__ --oci-password-variable GITHUB_PAT palace \&\& spack -e . concretize -f \&\& spack -e . mirror create -a\n\n# Install the software|" .devcontainer/Dockerfile
sed -i 's|# Install the software|# Do this separately to cache keys\nRUN cd /opt/spack-environment \&\& \\\n\t spack -e . buildcache keys --install --trust\n\n# Install the software|' .devcontainer/Dockerfile
sed -i 's|# Install the software|# Do this to configure private mirror\nRUN --mount=type=bind,source=.,destination=/home/app,readwrite --mount=type=secret,id=.secrets,target=.secrets export \$(cat .secrets \| xargs) \&\& \\\n\t cd /opt/spack-environment \&\& \\\n\t spack -e . mirror set --oci-username-variable GITHUB_USERNAME --oci-password-variable GITHUB_PAT palace \&\& \\\n\t spack -e . repo add /home/app/spack/local \&\& \\\n\t spack -e . develop --path=/home/app local.palace@=develop \&\& \\\n\t spack -e . concretize -f \&\& \\\n\t spack -e . mirror create --all --dependencies --exclude-specs local.palace \n\n# Install the software|' .devcontainer/Dockerfile

# Lets make sure that we configure externals / compilers properly...
# NOTE: This relies on the comment string '# Do this' from the previous sed command
sed -i 's|# Do this separately to cache keys|# Lets make sure our externals / compilers are correct\nRUN cd /opt/spack-environment \&\& \\\n\t spack -e . external find --all --exclude python \&\& \\\n\t spack -e . compiler list\n\n# Do this separately to cache keys|' .devcontainer/Dockerfile

# Minor spack optimization...
# Since we concretize above, install command definitely shouldn't
# MAC SED COMMAND vvv
# sed -i "" "s|--fail-fast|--fail-fast --only-concrete|" .devcontainer/Dockerfile
sed -i 's|--fail-fast|--fail-fast --only-concrete|' .devcontainer/Dockerfile

# Bind-mount the PWD at build time to use the source code here
# We could also bind-mount the spack cache dirs / installation dirs for speed...
# MAC SED COMMAND vvv (this didn't work, and I didn't debug...)
# sed -i "" "s|cd /opt/spack-environment && spack env activate .|--mount=type=bind,source=.,destination=/home/app cd /opt/spack-environment && spack env activate .|" .devcontainer/Dockerfile
# This command looks different as I only want to run it once.
sed -i 's|cd /opt/spack-environment \&\& spack env activate .|--mount=type=bind,source=.,destination=/home/app,readwrite cd /opt/spack-environment \&\& \\\n\t spack env activate .|' .devcontainer/Dockerfile

# Cleanup temp files
rm .devcontainer/spack.yaml
rm -rfd .devconatiner/.spack-env

