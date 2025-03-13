##
#  CMake script for the step-3 tutorial program:
##

# Usually, you will not need to modify anything beyond this point...

cmake_minimum_required(VERSION 3.5.0)


SET(TARGET "sample")
SET(TARGET_SRC "sample.cc")

list(APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/cmake/modules
  )

find_package(deal.II 9.5.0
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
if(NOT ${deal.II_FOUND})
  message(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
endif()

set(WITH_NETCDF ON CACHE BOOL "Check if the user wants to compile with the NetCDF libraries.")

if(WITH_NETCDF)
  find_package(NETCDF)
  if(${NETCDF_FOUND})
    message(STATUS "Using WITH_NETCDF = '${ASPECT_WITH_NETCDF}'")
    message(STATUS "  NETCDF_INCLUDE_DIR: ${NETCDF_INCLUDE_DIR}")
    message(STATUS "  NETCDF_LIBRARY: ${NETCDF_LIBRARY}")
    message(STATUS "  NETCDF_VERSION: ${NETCDF_VERSION}")
  else()
    message(STATUS "NetCDF not found. Disabling WITH_NETCDF. You can specify a hint to your installation directory with NETCDF_DIR.")
    set(WITH_NETCDF OFF CACHE BOOL "" FORCE)
  endif()
else()
  message(STATUS "Using WITH_NETCDF = 'OFF'")
endif()



deal_ii_initialize_cached_variables()

project(${TARGET})
deal_ii_invoke_autopilot()
