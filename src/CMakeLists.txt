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

# include self source
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

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
	add_definitions(-DEIGEN_DONT_ALIGN)
#endif()

# -------------------------------------
# Setting up target
# -------------------------------------

set(LIB_NAME ${PROJECT_NAME})

file(GLOB ${LIB_NAME}_HEADERS 
	${PROJECT_SOURCE_DIR}/include/${${PROJECT_NAME}_HEADER_DIR}/*.h)
file(GLOB ${LIB_NAME}_PRIV_HEADERS 
	${CMAKE_CURRENT_SOURCE_DIR}/*.h)
file(GLOB ${LIB_NAME}_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/*.cc)

if(QSERL_BUILD_SHARED)
	add_definitions(-DQSERL_DLL_EXPORTS)
	add_definitions(-DQSERL_USE_SHARED)
	add_library(${LIB_NAME} SHARED
		${${LIB_NAME}_SOURCES}
		${${LIB_NAME}_HEADERS}
		${${LIB_NAME}_PRIV_HEADERS})
else(QSERL_BUILD_SHARED)
	add_library(${LIB_NAME} STATIC
		${${LIB_NAME}_SOURCES}
		${${LIB_NAME}_HEADERS}
		${${LIB_NAME}_PRIV_HEADERS})
endif(QSERL_BUILD_SHARED)

set_target_properties(${LIB_NAME} PROPERTIES 
	VERSION ${PROJECT_VERSION}
	SOVERSION ${PROJECT_SOVERSION})

target_link_libraries(${LIB_NAME} LINK_INTERFACE_LIBRARIES ${LIBS})

if(MSVC)
	# debug suffix already added at this point
	set_target_properties(${PROJECT_NAME} PROPERTIES 
		DEBUG_OUTPUT_NAME "${PROJECT_NAME}d-${PLATFORM_COMPILER}" 
		RELEASE_OUTPUT_NAME "${PROJECT_NAME}-${PLATFORM_COMPILER}")

else(MSVC)
	set_target_properties(${PROJECT_NAME} PROPERTIES 
		VERSION "${PROJECT_VERSION}" 
		SOVERSION "${PROJECT_SOVERSION}")
endif(MSVC)

# -------------------------------------
# Installation rules
# -------------------------------------

# Set /bin and /lib directories
set(OUTPUT_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin")
make_directory(${OUTPUT_BIN_DIR})
set(OUTPUT_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib")
make_directory(${OUTPUT_LIB_DIR})
set(OUTPUT_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include/${PROJECT_NAME}")
make_directory(${OUTPUT_INCLUDE_DIR})

if(WIN32)
  install(TARGETS ${LIB_NAME} 
    RUNTIME DESTINATION ${OUTPUT_BIN_DIR}
    LIBRARY DESTINATION ${OUTPUT_BIN_DIR}
    ARCHIVE DESTINATION ${OUTPUT_LIB_DIR})
else(WIN32)
  install(TARGETS ${LIB_NAME} 
    RUNTIME DESTINATION ${OUTPUT_BIN_DIR}
    LIBRARY DESTINATION ${OUTPUT_LIB_DIR}
    ARCHIVE DESTINATION ${OUTPUT_LIB_DIR})
endif(WIN32)

install(FILES ${${LIB_NAME}_HEADERS} DESTINATION ${OUTPUT_INCLUDE_DIR})

