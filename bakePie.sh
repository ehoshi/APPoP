#!/bin/sh

module purge
module load openmpi/4.0.3-gcc/8.3.1-ucx
rm -r build

mkdir build
cd build
cmake -D "CMAKE_CXX_FLAGS=-Og" \
      -D "CMAKE_BUILD_TYPE=Debug" \
      -D "USE_BUNDLED_ARMADILLO=On" \
      ..
make "$@"
