# Copyright (c) 2012-2018 CNRS
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

cmake_minimum_required(VERSION 3.8.2)

project(libqserl VERSION 1.1.0 LANGUAGES CXX)

# Use cmake/ dir for user modules
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)

# Option for building shared lib
# XXX Building shared DLL library is disabled since the library exposes directly interfaces
# using the STL library
# FIXME
option(QSERL_BUILD_SHARED ON "Build shared libraries")

# Option for building tests
option(QSERL_BUILD_TEST ON "Build tests")

#------------------------------------------------------------------------------
# Dependencies
#------------------------------------------------------------------------------

find_package(Boost 1.55 REQUIRED MODULE)
find_package(EigenLgsm REQUIRED CONFIG)
find_package(cereal REQUIRED CONFIG)

#------------------------------------------------------------------------------
# Setting up target
#------------------------------------------------------------------------------

set(qserl_SOURCES
  src/rod2d/analytic_dqda.cc
  src/rod2d/analytic_energy.cc
  src/rod2d/analytic_mu.cc
  src/rod2d/analytic_q.cc
  src/rod2d/costate_system.cc
  src/rod2d/inverse_geometry.cc
  src/rod2d/jacobian_system.cc
  src/rod2d/rod.cc
  src/rod2d/state_system.cc
  src/rod2d/workspace_integrated_state.cc
  src/rod2d/workspace_state.cc
  src/rod3d/costate_system.cc
  src/rod3d/jacobian_system.cc
  src/rod3d/parameters.cc
  src/rod3d/rod.cc
  src/rod3d/state_system.cc
  src/rod3d/workspace_integrated_state.cc
  src/rod3d/workspace_state.cc
  src/util/lie_algebra_utils.cc
  src/util/timer.cc
  src/util/utils.cc
  )

if(QSERL_BUILD_SHARED)
	add_definitions(-DQSERL_DLL_EXPORTS)
	add_definitions(-DQSERL_USE_SHARED)
	add_library(qserl SHARED
    ${qserl_SOURCES})
else()
	add_library(qserl STATIC
    ${qserl_SOURCES})
endif()

target_include_directories(qserl
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
  )

target_compile_features(qserl PRIVATE cxx_std_11)
target_compile_options(qserl PRIVATE -Wall -Wextra)

target_link_libraries(qserl
  PUBLIC
    Boost::boost
    cereal
    EigenLgsm
  )

#------------------------------------------------------------------------------
# Installation rules
#------------------------------------------------------------------------------

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/Qserl)

install(TARGETS qserl
  EXPORT qserl-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

set_target_properties(PROPERTIES EXPORT_NAME Qserl)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Export the targets to a script
install(EXPORT qserl-targets
  FILE
  QserlTargets.cmake
  DESTINATION
  ${INSTALL_CONFIGDIR}
  )

# Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/QserlConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
  )

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/QserlConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/QserlConfig.cmake
  INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
  )

# Install the config, configversion and custom find modules
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/QserlConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/QserlConfigVersion.cmake
  DESTINATION ${INSTALL_CONFIGDIR}
  )

export(EXPORT qserl-targets FILE ${CMAKE_CURRENT_BINARY_DIR}/QserlTargets.cmake)

# Register package in user's package registry
export(PACKAGE Qserl)

# Uninstall target
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY
  )

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
  )
endif()

#------------------------------------------------------------------------------
# Add subdirectories
#------------------------------------------------------------------------------

add_subdirectory(examples)

if(QSERL_BUILD_TEST)
  add_subdirectory(test)
endif()