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

# Locate Qserl library.
#
# This script defines:
#   QSERL_FOUND, set to 1 if found
#   QSERL_LIBRARIES
#   QSERL_INCLUDE_DIR
#
# This script will look in standard locations for installed QSERL. However, 
# if you have installed it into a non-standard location, you can use the 
# QSERL_ROOT (environment or CMake) variable to specify the location.
#
# Inputs:
#   QSERL_DEBUG  Set to true if you want debug messages
#
# Warning: If using MSVC, the variable PLATFORM_COMPILER must have been set
# 	to the appropriate used Msvc compiler (e.g. vc100 for Msvc10).
#	Use the DETECT_MSVC_VERSION() macro from the detect_msvc_versions.cmake file.

if(MSVC)
	if(PLATFORM_COMPILER)
		set(QSERL_LIBRARY_SUFFIX "-${PLATFORM_COMPILER}")
	else(PLATFORM_COMPILER)
		message("The variable PLATFORM_COMPILER must be set under MSVC.")
	endif(PLATFORM_COMPILER)
else(MSVC)
	set(QSERL_LIBRARY_SUFFIX "")
endif(MSVC)

unset( QSERL_INCLUDE_DIR CACHE )
unset( QSERL_LIBRARY CACHE )
unset( QSERL_LIBRARY_DEBUG CACHE )
unset( QSERL_LIBRARIES CACHE )
mark_as_advanced( QSERL_INCLUDE_DIR )
mark_as_advanced( QSERL_LIBRARY )
mark_as_advanced( QSERL_LIBRARY_DEBUG )
find_path( QSERL_INCLUDE_DIR "qserl/exports.h"
  PATHS
      ${QSERL_ROOT}
      $ENV{QSERL_ROOT}
      /usr
      /usr/include
      /usr/local
      /usr/local/include
  PATH_SUFFIXES
      /include
  )

if(QSERL_DEBUG)
  message(STATUS "QSERL_INCLUDE_DIR=${QSERL_INCLUDE_DIR}")
  message(STATUS "Looking for library file: qserl${QSERL_LIBRARY_SUFFIX}.lib")
endif(QSERL_DEBUG)
find_library( QSERL_LIBRARY
  NAMES
      qserl${QSERL_LIBRARY_SUFFIX}
  PATHS
      ${QSERL_ROOT}
      $ENV{QSERL_ROOT}
      /lib
      /lib64
  PATH_SUFFIXES
      /lib
  )
    
if(QSERL_DEBUG)
  message(STATUS "QSERL_INCLUDE_DIR=${QSERL_INCLUDE_DIR}")
  message(STATUS "Looking for library file: qserld${QSERL_LIBRARY_SUFFIX}.lib")
endif(QSERL_DEBUG)

find_library( QSERL_LIBRARY_DEBUG
  NAMES
      qserld${QSERL_LIBRARY_SUFFIX}
  PATHS
      ${QSERL_ROOT}
      $ENV{QSERL_ROOT}
      /lib
      /lib64
  PATH_SUFFIXES
      /lib
  )
  
if(QSERL_DEBUG)
  message(STATUS "QSERL_LIBRARY=${QSERL_LIBRARY}")
  message(STATUS "QSERL_LIBRARY_DEBUG=${QSERL_LIBRARY_DEBUG}")
endif(QSERL_DEBUG)

if( QSERL_LIBRARY )
  set( QSERL_LIBRARIES ${QSERL_LIBRARIES}
  "optimized" ${QSERL_LIBRARY}
)
endif( QSERL_LIBRARY )
if( QSERL_LIBRARY_DEBUG )
  set( QSERL_LIBRARIES ${QSERL_LIBRARIES}
      "debug" ${QSERL_LIBRARY_DEBUG}
  )
endif( QSERL_LIBRARY_DEBUG )

if(QSERL_DEBUG)
  message(STATUS "QSERL_LIBRARIES=${QSERL_LIBRARIES}")
endif(QSERL_DEBUG)

set( QSERL_FOUND 0 )
if( QSERL_INCLUDE_DIR AND QSERL_LIBRARIES )
  set( QSERL_FOUND 1 )
endif( QSERL_INCLUDE_DIR AND QSERL_LIBRARIES )

if(QSERL_DEBUG)
  message(STATUS "QSERL_FOUND=${QSERL_FOUND}")
endif(QSERL_DEBUG)