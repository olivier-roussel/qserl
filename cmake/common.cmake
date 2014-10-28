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

macro(SUBDIRLIST result curdir)
  file(GLOB children RELATIVE ${curdir} ${curdir}/*)
  set(dirlist "")
  foreach(child ${children})
    if(IS_DIRECTORY ${curdir}/${child})
        set(dirlist ${dirlist} ${child})
    endif()
  endforeach()
  set(${result} ${dirlist})
endmacro()

# Recursive subroutine, for internal use.
# Use SOURCE_FROM_DIRECTORY instead.
macro(SOURCE_FROM_DIRECTORY_REC ROOTDIR SUBPATH)
if ("${SUBPATH}" STREQUAL "")
	set(FULL_CUR_DIR ${ROOTDIR})
else()
	set(FULL_CUR_DIR ${ROOTDIR}/${SUBPATH})
endif()

file(GLOB SUBDIR_HEADERS
	${FULL_CUR_DIR}/*.h 
	${FULL_CUR_DIR}/*.hh 
	${FULL_CUR_DIR}/*.hpp)

set(HEADER_FILES ${HEADER_FILES} ${SUBDIR_HEADERS})

file(GLOB SUBDIR_SOURCES
	${FULL_CUR_DIR}/*.cc 
	${FULL_CUR_DIR}/*.cpp)

set(SOURCE_FILES ${SOURCE_FILES} ${SUBDIR_SOURCES})

# add source group if not in root
if (NOT "${SUBPATH}" STREQUAL "")
	STRING(REGEX REPLACE "/" "\\\\" FULL_SUB_PATH_REPL ${SUBPATH})
	source_group("Header Files\\${FULL_SUB_PATH_REPL}" FILES ${SUBDIR_HEADERS}) 
	source_group("Source Files\\${FULL_SUB_PATH_REPL}" FILES ${SUBDIR_SOURCES}) 
endif()

SUBDIRLIST(SOURCES_SUBDIRS ${FULL_CUR_DIR})
FOREACH(SUBDIR ${SOURCES_SUBDIRS})
	#message("SUBDIR=${SUBDIR}")
	if ("${SUBPATH}" STREQUAL "")
		set(FULL_SUB_DIR ${SUBDIR})
	else()
		set(FULL_SUB_DIR ${SUBPATH}/${SUBDIR})
	endif()
	SOURCE_FROM_DIRECTORY_REC(${ROOTDIR} ${FULL_SUB_DIR})
ENDFOREACH()
endmacro()

macro(SOURCE_FROM_DIRECTORY ROOTDIR)
	set(HEADER_FILES "")
	set(SOURCE_FILES "")
	SOURCE_FROM_DIRECTORY_REC(${ROOTDIR} "")
endmacro()