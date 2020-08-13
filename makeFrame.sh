#!/bin/sh

module load cuda/10.2.89-gcc/8.3.1
module load cmake/3.17.3-gcc/8.3.1
module load openmpi/3.1.6-gcc/8.3.1-cuda10_2-ucx
module load armadillo/7.950.1-gcc/8.3.1-cuda10_2
#OR module load gcc/8.2.0 openmpi/4.0.1-gcc82-ucx

mkdir build
cd build
cmake -D "CMAKE_PREFIX_PATH=$ARMADILLO_DIR" \
      -D "CMAKE_CXX_FLAGS=-Og" \
      -D "CMAKE_BUILD_TYPE=Debug" \
      ..
make "$@"
