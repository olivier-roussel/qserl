/**
* Copyright (c) 2019 CNRS
* Author: Joseph Mirabel
* Inspired by file unittest/explog.cc of Pinocchio library.
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


#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

#include <qserl/util/explog.h>

using Eigen::Vector3d;
using Eigen::Matrix3d;

typedef Eigen::Matrix<double,6,1> Vector6d;
typedef Eigen::Matrix<double,4,4> Matrix4d;

BOOST_AUTO_TEST_SUITE ( explog )

BOOST_AUTO_TEST_CASE(exp)
{
  Matrix4d M(Matrix4d::Random());

  Vector6d v(Vector6d::Random());

  Matrix3d R = qserl::exp3(v.head<3>());
  BOOST_CHECK(R.isApprox(Eigen::AngleAxis<double>(v.head<3>().norm(), v.head<3>().normalized()).matrix()));
  BOOST_CHECK(v.head<3>().isApprox (qserl::log3(R)));

  Matrix3d R0 = qserl::exp3(Vector3d::Zero());
  BOOST_CHECK(R0.isIdentity());
  
  M = qserl::exp6(v);
  
  BOOST_CHECK(R.isApprox(M.topLeftCorner<3,3>()));
  BOOST_CHECK(v.isApprox (qserl::log6(M)));
}

BOOST_AUTO_TEST_SUITE_END()
