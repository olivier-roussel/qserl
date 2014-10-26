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

#ifndef QSERL_QUASI_STATIC_ROD_PARAMETERS_H_
#define QSERL_QUASI_STATIC_ROD_PARAMETERS_H_

#include "exports.h"

#include <Eigen/Core>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/version.hpp>

#include "util/eigen_types_serialization.h"

namespace qserl {

/**
* \brief Static parameters of a Kirchhoff 3D circular rod (homogeneous, istropic, linear elastic material).
*/
struct QSERL_EXPORT Parameters
{
	Parameters() :
		radius(0.01), 
		length(1.), 
		youngModulus(15.4e6),  /** Default Young modulus of rubber: 15.4 MPa */
		shearModulus(5.13e6),  /** Default Shear modulus of rubber: 5.13 MPa */
		density(1.1 * 10e3),   /** 1.10 kg/dm3 -> kg/m3, */ 
		integrationTime(1.),
		rodModel(RM_EXTENSIBLE_SHEARABLE),
		numNodes(100)
	{}

	enum RodModelT
	{
		RM_INEXTENSIBLE = 0,
		RM_EXTENSIBLE_SHEARABLE,
		RM_NUMBER_OF_ROD_MODELS
	};

	static const std::string getRodModelName(RodModelT i_model)
	{
		const static char * const model_names_array[] = { "INEXTENSIBLE", "EXTENSIBLE_SHEARABLE" };
		const static std::vector<std::string> v_model_names_array(model_names_array, model_names_array + RM_NUMBER_OF_ROD_MODELS);
		return v_model_names_array[static_cast<int>(i_model)];
	}

	double														radius; 
	double														length;
	double														youngModulus;			/** So called E Young modulus. */
	double														shearModulus;			/** So called G shear modulus. */
	double														density;		

	/**< These two parameters are only used when coupling inextensible/non-shearable static rod model
	* with extensible/shearable one.
	* They are required to express the ratio between stiffness coefficients for the supposed infinite values
	* of extension & shear elasticity values. 
	* \deprecated
	*/
	//double														extensionRatio;		/**< Ratio for the young modulus (thus for extension). */
	//double														shearingRatio;		/**< Ratio for the shear modulus (thus for shearing). */			

	/**< Internal use only. */
	double														integrationTime;	/**< Should be kept to 1 (default value). */

	RodModelT													rodModel;

	int																numNodes;

	/** Serialization */
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & boost::serialization::make_nvp("radius", radius) & 
			boost::serialization::make_nvp("length", length) & 
			boost::serialization::make_nvp("youngModulus", youngModulus) & 
			boost::serialization::make_nvp("shearModulus", shearModulus) & 
			boost::serialization::make_nvp("density", density) & 
			//boost::serialization::make_nvp("extensionRatio", extensionRatio) & 
			//boost::serialization::make_nvp("shearingRatio", shearingRatio) & 
			boost::serialization::make_nvp("integrationTime", integrationTime) &
			boost::serialization::make_nvp("rodModel", rodModel) &
			boost::serialization::make_nvp("numNodes", numNodes);
	}
};

}	// namespace qserl

#endif // QSERL_QUASI_STATIC_ROD_PARAMETERS_H_
