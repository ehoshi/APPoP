cmake_minimum_required(VERSION 3.10)

project(appop LANGUAGES CXX)
if(NOT USE_BUNDLED_ARMADILLO)
    find_package(Armadillo REQUIRED)
endif()
find_package(MPI REQUIRED)

add_executable(frame
  frame.cpp
  frm_FXN.cpp
  frm_MAS.cpp
  frm_SIM.cpp
  frm_WRK.cpp)
if(USE_BUNDLED_ARMADILLO)
    set(ARMADILLO_INCLUDE_DIRS armadillo-9.800.3/include)
endif()
target_compile_options(frame PRIVATE
  $<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:
  -Werror -Wall -pedantic-errors -Wextra -Wconversion -Wsign-conversion>)
target_include_directories(frame SYSTEM PRIVATE ${ARMADILLO_INCLUDE_DIRS})
set_target_properties(frame PROPERTIES
  CXX_STANDARD 11
  CXX_STANDARD_REQUIRED ON)
target_link_libraries(frame
  PRIVATE ${ARMADILLO_LIBRARIES}
  PRIVATE MPI::MPI_CXX)
