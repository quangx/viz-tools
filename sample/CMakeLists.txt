##
#  CMake script for the step-3 tutorial program:
##

project(me)
# Usually, you will not need to modify anything beyond this point...

cmake_minimum_required(VERSION 3.3.0)

find_package(deal.II 9.6.0
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )

if(NOT ${deal.II_FOUND})
  message(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
endif()


FIND_PACKAGE(Aspect 2.4.0 QUIET HINTS ${Aspect_DIR} ../ ../../ $ENV{ASPECT_DIR})

IF (NOT Aspect_FOUND)
  MESSAGE(FATAL_ERROR "\n"
    "Could not find a valid ASPECT build/installation directory. "
    "Please specify the directory where you are building ASPECT by passing\n"
    "   -D Aspect_DIR=<path to ASPECT>\n"
    "to cmake or by setting the environment variable ASPECT_DIR in your shell "
    "before calling cmake. See the section 'How to write a plugin' in the "
    "manual for more information.")
ENDIF ()

deal_ii_initialize_cached_variables()

SET(TARGET "sample")
PROJECT(${TARGET})

# ADD_LIBRARY(${TARGET} SHARED sample.cc)
add_executable(sample sample.cc)
deal_ii_setup_target(sample)
ASPECT_SETUP_PLUGIN(sample)