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

#ifndef QSERL_2D_PARAMETERS_H_
#define QSERL_2D_PARAMETERS_H_

#include "qserl/exports.h"

#include <Eigen/Core>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/version.hpp>


namespace qserl {
namespace rod2d {

/**
* \brief Static parameters of a Kirchhoff 2D (planar) rod (homogeneous, istropic, linear elastic material).
*/
struct QSERL_EXPORT Parameters
{
	Parameters() :
		radius(0.01), 
		length(1.), 
		integrationTime(1.),
		delta_t(0.01),
		rodModel(RM_INEXTENSIBLE)
	{}

	enum RodModelT
	{
		RM_INEXTENSIBLE = 0,
		//RM_EXTENSIBLE,	/**< Not implemented yet */
		RM_NUMBER_OF_ROD_MODELS
	};

	static const std::string getRodModelName(RodModelT i_model)
	{
		const static char * const model_names_array[] = { "INEXTENSIBLE" };
		const static std::vector<std::string> v_model_names_array(model_names_array, model_names_array + RM_NUMBER_OF_ROD_MODELS);
		return v_model_names_array[static_cast<int>(i_model)];
	}

	double														radius;							/**< Rod radius _ Default is 1e-2. */
	double														length;							/**< Rod length _ Must be set to 1 _ Default is 1. */

	double														integrationTime;		/**< Internal use _ Default is 1. */
	double														delta_t;						/**< Integration step time, defines the resolution of the rod discretization. 
																												* Default is 1e-2. */

	RodModelT													rodModel;

	//int																numNodes;				/** Deprecated. Use numberOfNodes() instead. */
	int	numberOfNodes() const
	{ return static_cast<int>(std::floor(integrationTime/delta_t)) + 1; }

	/** Serialization */
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & boost::serialization::make_nvp("radius", radius) & 
			boost::serialization::make_nvp("length", length) & 
			boost::serialization::make_nvp("integrationTime", integrationTime) &
			boost::serialization::make_nvp("delta_t", delta_t) &
			boost::serialization::make_nvp("rodModel", rodModel);
	}
};

}	// namespace rod2d
}	// namespace qserl

#endif // QSERL_2D_PARAMETERS_H_
