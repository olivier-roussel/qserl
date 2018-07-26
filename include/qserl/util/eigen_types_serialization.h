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

#ifndef QSERL_UTIL_EIGEN_TYPES_SERIALIZATION_H_
#define QSERL_UTIL_EIGEN_TYPES_SERIALIZATION_H_

#include <boost/serialization/serialization.hpp>
#pragma warning( push, 0 )	
#include <Eigen/Lgsm>
#pragma warning( pop )	

// Seriliazation for Eigen types

namespace boost {
namespace serialization {

template<class Archive, typename Scalar>
void serialize(Archive & ar, Eigen::Matrix<Scalar, 3, 1>& v, const unsigned int version)
{
	ar & boost::serialization::make_nvp("x", v.x()) & 
		boost::serialization::make_nvp("y", v.y()) & 
		boost::serialization::make_nvp("z", v.z());
}

template<class Archive, typename Scalar>
void serialize(Archive & ar, Eigen::Twist<Scalar>& w, const unsigned int version)
{
	ar & boost::serialization::make_nvp("rx", w.rx()) & 
			boost::serialization::make_nvp("ry", w.ry()) & 
			boost::serialization::make_nvp("rz", w.rz()) & 
			boost::serialization::make_nvp("vx", w.vx()) & 
			boost::serialization::make_nvp("vy", w.vy()) & 
			boost::serialization::make_nvp("vz", w.vz());
}

template<class Archive, typename Scalar>
void serialize(Archive & ar, Eigen::Matrix<Scalar, 6, 1>& v, const unsigned int version)
{
  for (int idxCol = 0 ; idxCol < 6 ; ++idxCol)
  {
    const std::string name = std::string("c") + std::to_string(idxCol);
    ar & boost::serialization::make_nvp(name.c_str(), v[idxCol]);
  }
}

template<class Archive, typename Scalar>
void serialize(Archive & ar, Eigen::Wrench<Scalar>& w, const unsigned int version)
{
	ar & boost::serialization::make_nvp("tx", w.tx()) & 
			boost::serialization::make_nvp("ty", w.ty()) & 
			boost::serialization::make_nvp("tz", w.tz()) & 
			boost::serialization::make_nvp("fx", w.fx()) & 
			boost::serialization::make_nvp("fy", w.fy()) & 
			boost::serialization::make_nvp("fz", w.fz());
}

template<class Archive, typename Scalar>
void serialize(Archive & ar, Eigen::Displacement<Scalar>& d, const unsigned int version)
{
	ar & boost::serialization::make_nvp("x", d.x()) & 
			boost::serialization::make_nvp("y", d.y()) & 
			boost::serialization::make_nvp("z", d.z()) & 
			boost::serialization::make_nvp("qx", d.qx()) & 
			boost::serialization::make_nvp("qy", d.qy()) & 
			boost::serialization::make_nvp("qz", d.qz()) & 
			boost::serialization::make_nvp("qw", d.qw());
}

} // namespace serialization
} // namespace boost

#endif // QSERL_UTIL_EIGEN_TYPES_SERIALIZATION_H_
