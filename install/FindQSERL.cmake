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

UNSET( QSERL_INCLUDE_DIR CACHE )
UNSET( QSERL_LIBRARY CACHE )
UNSET( QSERL_LIBRARY_DEBUG CACHE )
UNSET( QSERL_LIBRARIES CACHE )
MARK_AS_ADVANCED( QSERL_INCLUDE_DIR )
MARK_AS_ADVANCED( QSERL_LIBRARY )
MARK_AS_ADVANCED( QSERL_LIBRARY_DEBUG )
FIND_PATH( QSERL_INCLUDE_DIR "qserl/exports.h"
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

IF(QSERL_DEBUG)
  MESSAGE(STATUS "QSERL_INCLUDE_DIR=${QSERL_INCLUDE_DIR}")
ENDIF(QSERL_DEBUG)
FIND_LIBRARY( QSERL_LIBRARY
  NAMES
      qsserl
  PATHS
      ${QSERL_ROOT}
      $ENV{QSERL_ROOT}
      /lib
      /lib64
  PATH_SUFFIXES
      /lib
  )
    
FIND_LIBRARY( QSERL_LIBRARY_DEBUG
  NAMES
      qsserld
  PATHS
      ${QSERL_ROOT}
      $ENV{QSERL_ROOT}
      /lib
      /lib64
  PATH_SUFFIXES
      /lib
  )
  
IF(QSERL_DEBUG)
  MESSAGE(STATUS "QSERL_LIBRARY=${QSERL_LIBRARY}")
  MESSAGE(STATUS "QSERL_LIBRARY_DEBUG=${QSERL_LIBRARY_DEBUG}")
ENDIF(QSERL_DEBUG)

IF( QSERL_LIBRARY )
  SET( QSERL_LIBRARIES ${QSERL_LIBRARIES}
  "optimized" ${QSERL_LIBRARY}
)
ENDIF( QSERL_LIBRARY )
IF( QSERL_LIBRARY_DEBUG )
  SET( QSERL_LIBRARIES ${QSERL_LIBRARIES}
      "debug" ${QSERL_LIBRARY_DEBUG}
  )
ENDIF( QSERL_LIBRARY_DEBUG )

IF(QSERL_DEBUG)
  MESSAGE(STATUS "QSERL_LIBRARIES=${QSERL_LIBRARIES}")
ENDIF(QSERL_DEBUG)

SET( QSERL_FOUND 0 )
IF( QSERL_INCLUDE_DIR AND QSERL_LIBRARIES )
  SET( QSERL_FOUND 1 )
ENDIF( QSERL_INCLUDE_DIR AND QSERL_LIBRARIES )

IF(QSERL_DEBUG)
  MESSAGE(STATUS "QSERL_FOUND=${QSERL_FOUND}")
ENDIF(QSERL_DEBUG)