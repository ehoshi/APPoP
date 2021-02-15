#!/bin/sh

module purge
#module load armadillo/9.800.3-gcc/8.3.1-cuda10_2
module load openmpi/4.0.3-gcc/8.3.1-ucx
#module load openmpi/3.1.5-gcc/7.1.0-ucx
#module load gcc/8.2.0 openmpi/4.0.1-gcc82-ucx

mkdir build
cd build
cmake -D "CMAKE_CXX_FLAGS=-Og" \
      -D "CMAKE_BUILD_TYPE=Debug" \
      -D "USE_BUNDLED_ARMADILLO=On" \
      ..
make "$@"
