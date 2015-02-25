# Copyright (c) 2012-2014 CNRS
# Author: Olivier Roussel
#
# This file is part of the qserl package.
# qserl is free software: you can redistribute it
# and/or modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation, either version
# 3 of the License, or (at your option) any later version.
#
# qserl is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Lesser Public License for more details.  You should have
# received a copy of the GNU Lesser General Public License along with
# qserl.  If not, see
# <http://www.gnu.org/licenses/>.

include_directories("${PROJECT_SOURCE_DIR}/src")

# -------------------------------------
# Dependencies
# -------------------------------------

# Use Boost
if(MSVC)
	add_definitions(-DBOOST_ALL_DYN_LINK)
endif(MSVC)
if(IS_ICPC)
    set(Boost_USE_STATIC_LIBS ON CACHE STRING "Use statically linked Boost libraries")
else(IS_ICPC)
    # Ensure dynamic linking with boost unit_test_framework
    add_definitions(-DBOOST_TEST_DYN_LINK)
endif(IS_ICPC)
find_package(Boost 1.55 COMPONENTS chrono date_time filesystem program_options regex serialization system thread unit_test_framework REQUIRED)
if (Boost_FOUND)
  include_directories("${Boost_INCLUDE_DIR}")
  set(LIBS ${LIBS} ${Boost_LIBRARIES}) 
else()
  message("Could not find Boost")
endif()

# Use Eigen3 (with LGSM module)
if (NOT $ENV{EIGEN3_INCLUDE_DIR} STREQUAL "")
  set(EIGEN3_INCLUDE_DIR $ENV{EIGEN3_INCLUDE_DIR})
endif()
find_package(Eigen3 REQUIRED)
if (EIGEN3_FOUND)
  include_directories(${EIGEN3_INCLUDE_DIR})
else()
  message("Could not find Eigen3")
endif()
# Set EIGEN_DONT_ALIGN definition for 32bits platforms
# XXX crashes if not set also under MSVC10 Win64 _ 
#if(NOT CMAKE_CL_64)
#	add_definitions(-DEIGEN_DONT_ALIGN)
#endif()

# -------------------------------------
# Setting up target
# -------------------------------------
 
set(EXECUTABLE_NAME qserl-tests)

# Getting sources
SOURCE_FROM_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR})
set(${EXECUTABLE_NAME}_SOURCES ${SOURCE_FILES})
set(${EXECUTABLE_NAME}_HEADERS ${HEADER_FILES})

get_property(${PROJECT_NAME}_LIBRARY_RELEASE TARGET qserl PROPERTY LOCATION_RELEASE)
get_property(${PROJECT_NAME}_LIBRARY_DEBUG TARGET qserl PROPERTY LOCATION_DEBUG)
set(${PROJECT_NAME}_LIBRARIES ${${PROJECT_NAME}_LIBRARIES}
	"debug" ${${PROJECT_NAME}_LIBRARY_DEBUG})
set(${PROJECT_NAME}_LIBRARIES ${${PROJECT_NAME}_LIBRARIES}
	"optimized" ${${PROJECT_NAME}_LIBRARY_RELEASE})

add_executable(${EXECUTABLE_NAME} 
  ${${EXECUTABLE_NAME}_SOURCES} 
  ${${EXECUTABLE_NAME}_HEADERS} )

target_link_libraries(${EXECUTABLE_NAME} ${LIBS} ${${PROJECT_NAME}_LIBRARIES})

include(CMakeUnitTests.txt)

add_custom_command(TARGET ${EXECUTABLE_NAME} POST_BUILD 
	COMMAND ${EXECUTABLE_NAME} 
	--log_level=test_suite
	--result_code=no 
	--report_level=confirm 
	--detect_memory_leak=1) 

