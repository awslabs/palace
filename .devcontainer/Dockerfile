FROM ubuntu:24.04 AS bootstrap

ENV SPACK_ROOT=/opt/spack \
    CURRENTLY_BUILDING_DOCKER_IMAGE=1 \
    container=docker

ENV DEBIAN_FRONTEND=noninteractive   \
    LANGUAGE=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8

RUN apt-get -yqq update \
 && apt-get -yqq upgrade \
 && apt-get -yqq install --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        file \
        g++ \
        gcc \
        gfortran \
        git \
        gnupg2 \
        iproute2 \
        locales \
        make \
        mercurial \
        subversion \
        python3 \
        python3-boto3 \
        unzip \
        zstd \
 && locale-gen en_US.UTF-8 \
 && rm -rf /var/lib/apt/lists/*

RUN mkdir $SPACK_ROOT && cd $SPACK_ROOT && \
    git init --quiet && git remote add origin https://github.com/spack/spack.git && git fetch --depth=1 origin v1.0.0-alpha.4 && git checkout --detach FETCH_HEAD && \
    mkdir -p $SPACK_ROOT/opt/spack

RUN ln -s $SPACK_ROOT/share/spack/docker/entrypoint.bash \
          /usr/local/bin/docker-shell \
 && ln -s $SPACK_ROOT/share/spack/docker/entrypoint.bash \
          /usr/local/bin/interactive-shell \
 && ln -s $SPACK_ROOT/share/spack/docker/entrypoint.bash \
          /usr/local/bin/spack-env

RUN mkdir -p /root/.spack \
 && cp $SPACK_ROOT/share/spack/docker/modules.yaml \
        /root/.spack/modules.yaml \
 && rm -rf /root/*.* /run/nologin

# [WORKAROUND]
# https://superuser.com/questions/1241548/
#     xubuntu-16-04-ttyname-failed-inappropriate-ioctl-for-device#1253889
RUN [ -f ~/.profile ]                                               \
 && sed -i 's/mesg n/( tty -s \\&\\& mesg n || true )/g' ~/.profile \
 || true


WORKDIR /root
SHELL ["docker-shell"]

# Creates the package cache
RUN spack bootstrap now \
    && spack bootstrap status --optional \
    && spack spec hdf5+mpi

ENTRYPOINT ["/bin/bash", "/opt/spack/share/spack/docker/entrypoint.bash"]
CMD ["interactive-shell"]

# Build stage with Spack pre-installed and ready to be used
FROM bootstrap AS builder

# Install OS packages needed to build the software
RUN apt-get -yqq update && apt-get -yqq upgrade \
 && apt-get -yqq install autoconf gcc-14 g++-14 \
 && rm -rf /var/lib/apt/lists/*

# What we want to install and how we want to install it
# is specified in a manifest file (spack.yaml)
RUN mkdir -p /opt/spack-environment && \
set -o noclobber \
&&  (echo spack: \
&&   echo '  specs:' \
&&   echo '  - local.palace@develop' \
&&   echo '  concretizer:' \
&&   echo '    unify: true' \
&&   echo '    reuse: true' \
&&   echo '  mirrors:' \
&&   echo '    develop: https://binaries.spack.io/develop' \
&&   echo '    palace: oci://ghcr.io/awslabs/palace' \
&&   echo '  packages:' \
&&   echo '    petsc:' \
&&   echo '      require: ~hdf5' \
&&   echo '    all:' \
&&   echo '      providers:' \
&&   echo '        compiler: [gcc]' \
&&   echo '        mpi: [openmpi]' \
&&   echo '        blas: [intel-oneapi-mkl]' \
&&   echo '        lapack: [intel-oneapi-mkl]' \
&&   echo '  config:' \
&&   echo '    install_tree: /opt/software' \
&&   echo '  view: /opt/views/view') > /opt/spack-environment/spack.yaml

# Lets make sure our externals / compilers are correct
RUN cd /opt/spack-environment && \
	 spack -e . external find --all --exclude python && \
	 spack -e . compiler list

# Do this separately to cache keys
RUN cd /opt/spack-environment && \
	 spack -e . buildcache keys --install --trust

# Do this to configure private mirror
RUN --mount=type=bind,source=.,destination=/home/app,readwrite --mount=type=secret,id=.secrets,target=.secrets export $(cat .secrets | xargs) && \
	 cd /opt/spack-environment && \
	 spack -e . mirror set --oci-username-variable GITHUB_USERNAME --oci-password-variable GITHUB_PAT palace && \
	 spack -e . repo add /home/app/spack/local && \
	 spack -e . develop --path=/home/app local.palace@=develop && \
	 spack -e . concretize -f && \
	 spack -e . mirror create --all --dependencies --exclude-specs local.palace 

# Install the software, remove unnecessary deps
RUN --mount=type=bind,source=.,destination=/home/app,readwrite cd /opt/spack-environment && \
	 spack env activate . && spack install --fail-fast --only-concrete && spack gc -y

# Strip all the binaries
RUN find -L /opt/views/view/* -type f -exec readlink -f '{}' \; | \
    xargs file -i | \
    grep 'charset=binary' | \
    grep 'x-executable\|x-archive\|x-sharedlib' | \
    awk -F: '{print $1}' | xargs strip

# Modifications to the environment that are necessary to run
RUN cd /opt/spack-environment && \
    spack env activate --sh -d . > activate.sh


# Bare OS image to run the installed executables
FROM ubuntu:24.04

COPY --from=builder /opt/spack-environment /opt/spack-environment
COPY --from=builder /opt/software /opt/software

# paths.view is a symlink, so copy the parent to avoid dereferencing and duplicating it
COPY --from=builder /opt/views /opt/views

RUN { \
      echo '#!/bin/sh' \
      && echo '.' /opt/spack-environment/activate.sh \
      && echo 'exec "$@"'; \
    } > /entrypoint.sh \
&& chmod a+x /entrypoint.sh \
&& ln -s /opt/views/view /opt/view


RUN apt-get -yqq update && apt-get -yqq upgrade \
 && apt-get -yqq install gfortran \
 && rm -rf /var/lib/apt/lists/*
ENTRYPOINT [ "/entrypoint.sh" ]
CMD [ "/bin/bash" ]

# Make sure devcontainer user gets spack packages
RUN echo "source /entrypoint.sh" >> /home/vscode/.bashrc

# Configure user for container
USER vscode
