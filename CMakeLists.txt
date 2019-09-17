cmake_minimum_required(VERSION 2.8)

project(rbf LANGUAGES Fortran)

if (DEFINED ENV{GSL} AND (NOT DEFINED ENV{GSL_ROOT}))
  set(ENV{GSL_ROOT} $ENV{GSL})
endif ()
if (DEFINED ENV{GSL_ROOT})
  set(USE_GSL ON)
  add_definitions(-DUSE_GSL)
  include_directories("$ENV{GSL_ROOT}/include")
  link_directories("$ENV{GSL_ROOT}/lib")
endif ()

set(sources
  src/const_mod.F90
  src/node_placing_mod.F90
  src/kernel_mod.F90
  src/rbf_mod.F90
)

add_library(rbf ${sources})
if (USE_GSL)
  target_link_libraries(rbf gsl)
endif ()
