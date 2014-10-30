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

# Check for minimum version
cmake_minimum_required(VERSION 2.8)

include(cmake/detect_msvc_versions.cmake)
include(cmake/common.cmake)
include(cmake/install.cmake)

set(PROJECT_NAME "qserl")
set(PROJECT_VERSION 1.0)
set(PROJECT_SOVERSION 0)
set(PROJECT_DESCRIPTION
  "Quasi-static Elastic Rods model library.")
set(PROJECT_URL "")

project(${PROJECT_NAME})

if(MSVC)
	DETECT_MSVC_VERSION()
endif(MSVC)

# Define & include header directory
set(${PROJECT_NAME}_HEADER_DIR qserl)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

# Use cmake/ dir for user modules
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Auto-install
set(CMAKE_INSTALL_PREFIX "$ENV{DEVEL_DIR}/install/unstable/" CACHE PATH "Prefix prepended to install directories." FORCE)

# Option for building shared lib
# XXX Building shared DLL library is disabled since the library exposes directly interfaces
# using the STL library
#option(QSERL_BUILD_SHARED OFF "Build shared libraries")

# Option for building tests
option(QSERL_BUILD_TEST ON "Build tests")

# -------------------------------------
# Visit subdirectories
# -------------------------------------

if(NOT SUBMODULE_MODE)
	# Standalone mode
	set(SUBMODULE_MODE ON)
	
	set(PROJECT_GROUP "Modules")
	
	# add submodules if any

	# Source directory
	set(PROJECT_GROUP "")
	add_subdirectory(src)
	
	if (QSERL_BUILD_TEST)
		enable_testing()
		add_subdirectory(test)
	endif(QSERL_BUILD_TEST)

else()
	# Source directory
	add_subdirectory(src)
endif()