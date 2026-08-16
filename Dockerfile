# syntax=docker/dockerfile:1

FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive
ARG BUILD_JOBS=4

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# -----------------------------------------------------------------------------
# System build dependencies
# -----------------------------------------------------------------------------

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        cmake \
        curl \
        git \
        libnuma-dev \
        pkg-config \
        && \
    rm -rf /var/lib/apt/lists/*

# -----------------------------------------------------------------------------
# Gurobi
#
# The defaults use Gurobi 12.0.1 to match the current README and paper
# environment. To use 12.0.3 instead, override both build arguments:
#
#   --build-arg GUROBI_VERSION=12.0.3
#   --build-arg GUROBI_DIR_VERSION=1203
# -----------------------------------------------------------------------------

ARG GUROBI_VERSION=12.0.1
ARG GUROBI_DIR_VERSION=1201
ARG GUROBI_DOWNLOAD_SERIES=12.0

ENV GUROBI_HOME="/opt/gurobi${GUROBI_DIR_VERSION}/linux64"
ENV PATH="${GUROBI_HOME}/bin:${PATH}"
ENV LD_LIBRARY_PATH="${GUROBI_HOME}/lib"

RUN gurobi_archive="gurobi${GUROBI_VERSION}_linux64.tar.gz" && \
    curl -fsSL --retry 3 \
        "https://packages.gurobi.com/${GUROBI_DOWNLOAD_SERIES}/${gurobi_archive}" \
        -o "/tmp/${gurobi_archive}" && \
    tar -xzf "/tmp/${gurobi_archive}" -C /opt && \
    rm -f "/tmp/${gurobi_archive}" && \
    test -f "${GUROBI_HOME}/include/gurobi_c++.h" && \
    test -f "${GUROBI_HOME}/lib/libgurobi120.so" && \
    echo "${GUROBI_HOME}/lib" > /etc/ld.so.conf.d/gurobi.conf && \
    mkdir -p /opt/gurobi && \
    ldconfig

# Do not COPY gurobi.lic into the image. Mount it read-only at runtime.

# -----------------------------------------------------------------------------
# LEMON 1.3.1
#
# The artifact CMakeLists.txt directly uses both the source directory and the
# build directory because lemon/config.h is generated during configuration.
# -----------------------------------------------------------------------------

ARG LEMON_VERSION=1.3.1

ENV LEMON_SOURCE_DIR="/opt/lemon-${LEMON_VERSION}"
ENV LEMON_BUILD_DIR="/opt/lemon-${LEMON_VERSION}/build"

RUN curl -fsSL --retry 3 \
        "https://lemon.cs.elte.hu/pub/sources/lemon-${LEMON_VERSION}.tar.gz" \
        -o "/tmp/lemon-${LEMON_VERSION}.tar.gz" && \
    tar -xzf "/tmp/lemon-${LEMON_VERSION}.tar.gz" -C /opt && \
    rm -f "/tmp/lemon-${LEMON_VERSION}.tar.gz" && \
    cmake \
        -S "${LEMON_SOURCE_DIR}" \
        -B "${LEMON_BUILD_DIR}" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DCMAKE_INSTALL_LIBDIR=lib && \
    cmake --build "${LEMON_BUILD_DIR}" --parallel "${BUILD_JOBS}" && \
    cmake --install "${LEMON_BUILD_DIR}" && \
    test -f "${LEMON_SOURCE_DIR}/lemon/list_graph.h" && \
    test -f "${LEMON_BUILD_DIR}/lemon/config.h" && \
    echo "/usr/local/lib" > /etc/ld.so.conf.d/lemon.conf && \
    ldconfig

# -----------------------------------------------------------------------------
# Python environment
#
# This is needed because the README runs:
#   conda run -n rtss26-figures ...
# -----------------------------------------------------------------------------

ARG MINIFORGE_INSTALLER=Miniforge3-Linux-x86_64.sh

ENV CONDA_DIR=/opt/conda

RUN curl -fsSL --retry 3 \
        "https://github.com/conda-forge/miniforge/releases/latest/download/${MINIFORGE_INSTALLER}" \
        -o /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p "${CONDA_DIR}" && \
    rm -f /tmp/miniforge.sh

ENV PATH="${CONDA_DIR}/bin:${PATH}"

COPY environment.yml /tmp/environment.yml

RUN conda env create --file /tmp/environment.yml && \
    conda clean --all --yes && \
    rm -f /tmp/environment.yml && \
    conda run -n rtss26-figures python --version && \
    conda run -n rtss26-figures jupyter lab --version

ENV CONDA_DEFAULT_ENV=rtss26-figures
ENV PATH="${CONDA_DIR}/envs/rtss26-figures/bin:${PATH}"

# Preserve the artifact environment when commands use "bash -lc".
RUN printf '%s\n' \
        "export GUROBI_HOME=\"${GUROBI_HOME}\"" \
        "export CONDA_DIR=\"${CONDA_DIR}\"" \
        'export PATH="${CONDA_DIR}/envs/rtss26-figures/bin:${CONDA_DIR}/bin:${GUROBI_HOME}/bin:${PATH}"' \
        'export LD_LIBRARY_PATH="${GUROBI_HOME}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"' \
        > /etc/profile.d/rtss26.sh

# -----------------------------------------------------------------------------
# Project
#
# The Dockerfile and environment.yml are expected to be in the repository root.
# -----------------------------------------------------------------------------

ENV PROJECT_DIR=/opt/RTSS26-Multi-Telescope-Transient-Follow-Up-Search

WORKDIR ${PROJECT_DIR}

COPY . .

RUN cmake \
        -S . \
        -B build \
        -DCMAKE_BUILD_TYPE=Release \
        -DLEMON_SOURCE_DIR="${LEMON_SOURCE_DIR}" \
        -DLEMON_BUILD_DIR="${LEMON_BUILD_DIR}" && \
    cmake --build build --parallel "${BUILD_JOBS}" && \
    test -x build/ts_maxp && \
    test -x build/ts_mink

WORKDIR ${PROJECT_DIR}

CMD ["/bin/bash"]