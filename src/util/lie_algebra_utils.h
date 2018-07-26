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


#ifndef QSERL_UTIL_LIE_ALGEBRA_UTILS_H_
#define QSERL_UTIL_LIE_ALGEBRA_UTILS_H_

//#pragma warning( push, 0 )
#include <Eigen/Lgsm>
//#pragma warning( pop )

namespace qserl {
namespace util {

// TODO: pass matrix as output parameter instead of returned value
inline Eigen::Matrix3d
hat(const Eigen::Vector3d& i_v)
{
  Eigen::Matrix3d out;
  out << 0., -i_v[2], i_v[1],
      i_v[2], 0., -i_v[0],
      -i_v[1], i_v[0], 0.;
  return out;
}

inline void
hat_SE2(const Eigen::Vector3d& i_v,
        Eigen::Matrix3d& o_mat)
{
  o_mat << 0., -i_v[2], i_v[0],
      i_v[2], 0., i_v[1],
      0., 0., 0.;
}

Eigen::Matrix4d
GetHomogenousMatrix(const Eigen::Displacementd& disp);

Eigen::Matrix4d
GetTransformationMatrix(const Eigen::Matrix4d& H,
                        const Eigen::Vector3d& scale);

Eigen::Matrix4d
GetTransformationMatrix(const Eigen::Displacementd& disp,
                        const Eigen::Vector3d& scale);

Eigen::Vector3d
GetScaleFromMatrix(const Eigen::Matrix4d& TM);

Eigen::Vector3d
GetTranslationFromMatrix(const Eigen::Matrix4d& TM);

Eigen::Rotation3d
GetRotationFromMatrix(const Eigen::Matrix4d& TM);

Eigen::Matrix3d
rotationMatrixFromBryantAngles(double rx,
                               double ry,
                               double rz);

Eigen::Matrix4d
transformationMatrixFromBryantAngles(double x,
                                     double y,
                                     double z,
                                     double rx,
                                     double ry,
                                     double rz);

void
transformationMatrixToBryantAngles(const Eigen::Matrix4d& trans,
                                   double& x,
                                   double& y,
                                   double& z,
                                   double& rx,
                                   double& ry,
                                   double& rz);

void
rotationMatrixToBryantAngles(const Eigen::Matrix3d& rot,
                             double& rx,
                             double& ry,
                             double& rz);

} // namespace util
} // namespace qserl

#endif // QSERL_UTIL_LIE_ALGEBRA_UTILS_H_
