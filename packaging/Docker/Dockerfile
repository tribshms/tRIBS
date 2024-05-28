FROM alpine:3.19

# Define version numbers as arguments
ARG BUILD_BASE_VERSION=0.5-r3
ARG CMAKE_VERSION=3.27.8-r0
ARG BASH_VERSION=5.2.21-r0
ARG ZSH_VERSION=5.9-r2
ARG PERL_VERSION=5.38.2-r0
ARG PYTHON_VERSION=3.11.9-r0
ARG OPENMPI_VERSION=5.0.2
ARG SHADOW_VERSION=4.14.2-r0

# Update and install dependencies
RUN apk update && \
    apk add --no-cache \
        build-base=${BUILD_BASE_VERSION} \
        cmake=${CMAKE_VERSION} \
        bash=${BASH_VERSION} \
        zsh=${ZSH_VERSION} \
        perl=${PERL_VERSION} \
        python3=${PYTHON_VERSION} \
        shadow=${SHADOW_VERSION}\
        openssh\
        zlib-dev\
    && rm -rf /var/cache/apk/*

# Download and install OpenMPI
ADD https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-${OPENMPI_VERSION}.tar.bz2 .
RUN tar xf openmpi-${OPENMPI_VERSION}.tar.bz2 \
    && cd openmpi-${OPENMPI_VERSION} \
    && ./configure \
    && make -j4 all \
    && make install \
    && cd .. && rm -rf openmpi-${OPENMPI_VERSION} openmpi-${OPENMPI_VERSION}.tar.bz2 /tmp/*

# Set environment variables for OpenMPI
ENV PATH="/usr/local/bin:$PATH"
ENV LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"
ENV MPICXX=/usr/local/bin/mpicxx

# Set the working directory
WORKDIR /tribs

# Copy source code
COPY src/ ./src/

# Copy CMakeLists.txt for parallel build
COPY CMakeLists.txt ./

# Build parallel
WORKDIR /tribs
RUN cmake -B cmake-build-parallel -S . -Dparallel=ON -DCMAKE_BUILD_TYPE=Release && \
    cmake --build cmake-build-parallel --target all

# Build serial
RUN cmake -B cmake-build-serial -S . -Dparallel=OFF -DCMAKE_BUILD_TYPE=Release && \
    cmake --build cmake-build-serial --target all

# Create a data directory
WORKDIR /tribs/
RUN mkdir "/tribs/shared"&& mkdir "/tribs/bin" &&\
    cp /tribs/cmake-build-serial/tRIBS /tribs/bin/ \
    && cp /tribs/cmake-build-parallel/tRIBSpar /tribs/bin/ \
    && rm -rf /tribs/cmake-build-parallel/ && rm -rf /tribs/cmake-build-serial/ && rm -rf /tribs/CMakeLists.txt && rm -rf /tribs/src/

# Switch to the non-root user
# add user for non-root access of openmpi
RUN useradd -ms /bin/bash tRIBSuser
USER tRIBSuser