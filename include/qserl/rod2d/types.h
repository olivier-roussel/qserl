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

#ifndef QSERL_2D_TYPES_H_
#define QSERL_2D_TYPES_H_

#include "qserl/exports.h"

#include "qserl/rod3d/types.h"

#include <array>

namespace qserl {
namespace rod2d {

/**< For convenience, a 2-dimensional wrench (with reference to wrenches in screw theory
* for 3-dimensional bodies) is stored in an array where index storage represent:
*	0: Force along X axis
*	1: Force along Y axis
*	2: Torque
* \warning This convention is different from the 3D case and the one used in Eigen/Lgsm, where 
*   rotational components are before translationnal ones. To be changed in next major release.
*/
typedef Eigen::Matrix<double, 3, 1> Wrench2D;

/**< For convenience, a 2-dimensional displacement (with reference to displacement in screw theory
* for 3-dimensional bodies is stored in an array where index storage represent:
*	0: Position along X axis
*	1: Position along Y axis
*	2: Rotation
* \warning This convention is different from the 3D case and the one used in Eigen/Lgsm, where 
*   rotational components are before translationnal ones. To be changed in next major release.
*/
typedef Eigen::Matrix<double, 3, 1> Displacement2D;

/**< Helper transforming a 2-dimensional pseudo-displacement into a 3D displacement
* in the XY plane at z=0.
*/
inline Eigen::Matrix4d
toHomogeneousMatrix(const Displacement2D& i_disp)
{
//  return se3::SE3(Eigen::Vector3d(i_disp[0], i_disp[1], 0.),
//                              Eigen::AngleAxisd(i_disp[2], Eigen::Vector3d::UnitZ()));
  const double cos_theta = cos(i_disp[2]);
  const double sin_theta = sin(i_disp[2]);
  Eigen::Matrix4d out;
  out << cos_theta, -sin_theta, 0., i_disp[0],
         sin_theta, cos_theta, 0., i_disp[1],
         0., 0., 1., 0.,
         0., 0., 0., 1.;
  return out;
}

/**< Helper transforming a 2-dimensional pseudo-displacement into a homogeneous 3x3 matrix.
*/
inline void
toHomogeneousMatrix(const Displacement2D& i_disp,
                    Eigen::Matrix3d& o_mat)
{
  const double cos_theta = cos(i_disp[2]);
  const double sin_theta = sin(i_disp[2]);
  o_mat << cos_theta, -sin_theta, i_disp[0],
      sin_theta, cos_theta, i_disp[1],
      0., 0., 1.;
}

/**< Helper transforming a 2-dimensional pseudo-wrench into a 3D wrench.
*/
inline rod3d::Wrench
toWrench3D(const Wrench2D& i_wrench)
{
  rod3d::Wrench w;
  w << 0., 0., i_wrench[2], i_wrench[0], i_wrench[1], 0.;
  return w;
}

/**< \brief Comparator for 2D wrenches.
* Returns true if w1 is strictly less than w2, false otherwise.
*/
template<typename T, typename U>
inline
bool
isLess(const T& i_w1,
       const U& i_w2)
{
  bool res = true;
  for(int k = 0; k < 3 && res; ++k)
  {
    res = i_w1[k] < i_w2[k];
  }
  return res;
}

}  // namespace rod2d
}  // namespace qserl

#endif // QSERL_2D_TYPES_H_
