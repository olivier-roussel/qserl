/**
* Copyright (c) 2012-2014 CNRS
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

#ifndef QSERL_FORWARD_CLASS_H_
#define QSERL_FORWARD_CLASS_H_

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

#define DECLARE_CLASS(C)																							 \
    class C;                                                           \
    typedef boost::shared_ptr<C> C##ShPtr;														 \
    typedef boost::weak_ptr<C> C##WkPtr;															 \
    typedef boost::shared_ptr<const C> C##ConstShPtr;									 \
    typedef boost::weak_ptr<const C> C##ConstWkPtr				 

#endif // QSERL_FORWARD_CLASS_H_
