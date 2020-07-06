#!/bin/sh

module load cmake/3.13.1
module load armadillo
module add gcc/6.1.0 openmpi/1.6.4-eth
#OR module load gcc/8.2.0 openmpi/4.0.1-gcc82-ucx

mkdir build
cd build
cmake -D "CMAKE_PREFIX_PATH=$ARMADILLO_DIR" \
      -D "CMAKE_CXX_FLAGS=-Og" \
      -D "CMAKE_BUILD_TYPE=Debug" \
      ..
make "$@"
