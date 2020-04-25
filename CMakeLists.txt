cmake_minimum_required(VERSION 3.10)

project(appop LANGUAGES CXX)

find_package(Armadillo REQUIRED)
find_package(MPI REQUIRED)

add_executable(frame
  frame.cpp
  frm_FXN.cpp
  frm_MAS.cpp
  frm_SIM.cpp
  frm_WRK.cpp)
target_compile_options(frame PRIVATE
  $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:
  -Werror -Wall -pedantic-errors -Wextra -Wconversion -Wsign-conversion>)
target_include_directories(frame SYSTEM PRIVATE ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(frame
  PRIVATE ${ARMADILLO_LIBRARIES}
  PRIVATE MPI::MPI_CXX)
