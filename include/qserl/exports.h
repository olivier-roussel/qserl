/**
* Copyright (c) 2012-2018 CNRS
* Author: Olivier Roussel
*
* This file is part of the qserl package.
* qserl is free software: you can redistribute it
* and/or modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation, either version
* 3 of the License, or (at your option) any later version.
*
* qserl is distributed in the hope that it will be
* useful, but WITHOUT ANY WARRANTY; without even the implied warranty
* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Lesser Public License for more details.  You should have
* received a copy of the GNU Lesser General Public License along with
* qserl.  If not, see
* <http://www.gnu.org/licenses/>.
**/

/** 
* \file exports.h
* \deprecated Not used as we do not want yet to export DLL symbols 
* since we expose STL interfaces
* \brief DLL symbol import / export  definitions for qserl
*/

#ifndef QSERL_EXPORTS_H_
#define QSERL_EXPORTS_H_

#if defined(WIN32) && defined(QSERL_USE_SHARED)
	#if defined(QSERL_DLL_EXPORTS)
		#undef QSERL_DLL_EXPORTS                      
		#define QSERL_EXPORT __declspec(dllexport)
		#define QSER_EXPIMP_TEMPLATE
	#else
		#define QSERL_EXPORT __declspec(dllimport)
		#define QSER_EXPIMP_TEMPLATE extern
	#endif
#else
	#undef QSERL_DLL_EXPORTS
	#define QSERL_EXPORT
#endif

#endif // QSERL_EXPORTS_H_

