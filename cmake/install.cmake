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

# BASE_DIR is the directory parent to the headers files given in HEADER_LIST as relative paths
MACRO(INSTALL_HEADERS_WITH_DIRECTORY BASE_DIR HEADER_LIST)
foreach(HEADER ${HEADER_LIST})
	string(REGEX MATCH "(.*)[/\\]" DIR ${HEADER})
	install(FILES ${BASE_DIR}/${HEADER} DESTINATION include/${DIR})
endforeach(HEADER)

ENDMACRO(INSTALL_HEADERS_WITH_DIRECTORY)