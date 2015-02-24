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

#include "lie_algebra_utils.h"

namespace qserl {
namespace util {

Eigen::Matrix4d GetHomogenousMatrix(const Eigen::Displacementd & disp)
{
	Eigen::Matrix4d res = Eigen::Matrix4d::Identity();
	res.block<3,1>(0,3) = disp.getTranslation();
	res.block<3,3>(0,0) = disp.getRotation().adjoint();
	return res;
}

Eigen::Matrix4d GetTransformationMatrix(const Eigen::Matrix4d & H, const Eigen::Vector3d & scale)
{
	Eigen::Matrix4d res = H;
	res.block<3,1>(0,0) *= scale.x();
	res.block<3,1>(0,1) *= scale.y();
	res.block<3,1>(0,2) *= scale.z();
	return res;
}

Eigen::Matrix4d GetTransformationMatrix(const Eigen::Displacementd & disp, const Eigen::Vector3d & scale)
{
	return GetTransformationMatrix(GetHomogenousMatrix(disp), scale);
}

Eigen::Vector3d GetScaleFromMatrix(const Eigen::Matrix4d &TM)
{
	double sign = 1.;
	if( TM.block<3,3>(0,0).determinant() < 0) sign = -1.; //inverting scaling
	return sign * Eigen::Vector3d(Eigen::Vector3d(TM(0,0),TM(1,0),TM(2,0)).norm(),
		Eigen::Vector3d(TM(0,1),TM(1,1),TM(2,1)).norm(),
		Eigen::Vector3d(TM(0,2),TM(1,2),TM(2,2)).norm());
}

Eigen::Vector3d GetTranslationFromMatrix(const Eigen::Matrix4d &TM)
{
	return Eigen::Vector3d(TM.block<3,1>(0,3));
}

Eigen::Rotation3d GetRotationFromMatrix(const Eigen::Matrix4d & TM)
{
	Eigen::Vector3d scale = GetScaleFromMatrix(TM);
	Eigen::Matrix3d rot;
	rot << TM(0,0)/scale[0], TM(0,1)/scale[1], TM(0,2)/scale[2],
		TM(1,0)/scale[0], TM(1,1)/scale[1], TM(1,2)/scale[2],
		TM(2,0)/scale[0], TM(2,1)/scale[1], TM(2,2)/scale[2];

	return Eigen::Rotation3d(Eigen::Quaterniond(rot));
}

Eigen::Matrix3d rotationMatrixFromBryantAngles(double rx, double ry, double rz)
{
	const double c_rx = cos(rx);
	const double s_rx = sin(rx);
	const double c_ry = cos(ry);
	const double s_ry = sin(ry);
	const double c_rz = cos(rz);
	const double s_rz = sin(rz);
	Eigen::Matrix3d rot;
	rot << c_ry*c_rz, -c_ry*s_rz, s_ry, 
		s_rx*s_ry*c_rz + c_rx*s_rz, -s_rx*s_ry*s_rz + c_rx*c_rz, -s_rx*c_ry,
		-c_rx*s_ry*c_rz + s_rx*s_rz, c_rx*s_ry*s_rz + s_rx*c_rz, c_rx*c_ry;
	return rot;
}

Eigen::Matrix4d transformationMatrixFromBryantAngles(double x, double y, double z, double rx, double ry, double rz)
{
	Eigen::Matrix4d rot = Eigen::Matrix4d::Zero();
	rot.block<3,3>(0,0) = rotationMatrixFromBryantAngles(rx, ry, rz);
	rot(0,3) = x;
	rot(1,3) = y;
	rot(2,3) = z;
	rot(3,3) = 1.;
	return rot;
}

void transformationMatrixToBryantAngles(const Eigen::Matrix4d& trans, double& x, double& y, double& z, double& rx, double& ry, double& rz)
{
	x = trans(0,3);
	y = trans(1,3);
	z = trans(2,3);
	rotationMatrixToBryantAngles(trans.block<3,3>(0,0), rx, ry, rz);
}

void rotationMatrixToBryantAngles(const Eigen::Matrix3d& rot, double& rx, double& ry, double& rz)
{
	static const double kCosTolSing = 1.e-6;
	const double s_ry = rot(0,2);
	ry = asin(s_ry);
	const double c_ry = cos(ry);
	if (abs(c_ry) > kCosTolSing)
	{
		const double inv_c_ry = 1. / c_ry;
		const double s_rx = -rot(1,2) * inv_c_ry;
		const double c_rx = rot(2,2) * inv_c_ry;
		const double s_rz = -rot(0,1) * inv_c_ry;
		const double c_rz = rot(0,0) * inv_c_ry;
		rx = atan2(s_rx, c_rx);
		rz = atan2(s_rz, c_rz);
	}else{
		rx = 0.;
		ry = 0.;
		rz = 0.;
	}
}

} // namespace util
} // namespace qserl
