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

#ifndef QSERL_3D_PARAMETERS_H_
#define QSERL_3D_PARAMETERS_H_

#include "qserl/exports.h"

#include <Eigen/Core>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/version.hpp>

#include "qserl/util/eigen_types_serialization.h"

namespace qserl {
namespace rod3d {

/**
* \brief Static parameters of a Kirchhoff 3D circular rod (homogeneous, istropic, linear elastic material).
*/
struct QSERL_EXPORT Parameters
{
	Parameters() :
		radius(0.01), 
		length(1.), 
		//youngModulus(15.4e6),  /** Default Young modulus of rubber: 15.4 MPa */
		//shearModulus(5.13e6),  /** Default Shear modulus of rubber: 5.13 MPa */
    stiffnessCoefficients(Eigen::Matrix<double, 6, 1>::Ones()),
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
  Eigen::Matrix<double, 6, 1>       stiffnessCoefficients;
	double														density;		        /** Unused. */

	/**< These two parameters are only used when coupling inextensible/non-shearable static rod model
	* with extensible/shearable one.
	* They are required to express the ratio between stiffness coefficients for the supposed infinite values
	* of extension & shear elasticity values. 
	* \deprecated
	*/
	//double														extensionRatio;		/**< Ratio for the young modulus (thus for extension). */
	//double														shearingRatio;		/**< Ratio for the shear modulus (thus for shearing). */			

	//double														youngModulus;			/** So called E Young modulus. */
	//double														shearModulus;			/** So called G shear modulus. */

	/**< Internal use only. */
	double														integrationTime;	/**< Should be kept to 1 (default value). */

	RodModelT													rodModel;

	int																numNodes;

  /**
  * \brief Returns equivalent Young modulus from stiffness coefficients.
  * \pre Stiffness coefficients must represent transversal isotropic elasticity, i.e.
  *   stiffnessCoefficients[1] == stiffnessCoefficients[2]    and
  *   stiffnessCoefficients[4] == stiffnessCoefficients[5]
  */
  double getIsotropicYoungModulus() const;

  /**
  * \brief Returns equivalent Shear modulus from stiffness coefficients.
  * \pre Stiffness coefficients must represent transversal isotropic elasticity, i.e.
  *   stiffnessCoefficients[1] == stiffnessCoefficients[2]    and
  *   stiffnessCoefficients[4] == stiffnessCoefficients[5]
  */
  double getIsotropicShearModulus() const;

  /**
  * \brief Set transversal isotropic stiffness coefficients from elasticity parameters described by
  * Young modulus and shear modulus. 
  */
  void setIsotropicStiffnessCoefficientsFromElasticityParameters(double i_youngModulus, double i_shearModulus);

	/** Serialization */
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & boost::serialization::make_nvp("radius", radius) & 
			boost::serialization::make_nvp("length", length) & 
			boost::serialization::make_nvp("stiffnessCoefficients", stiffnessCoefficients) & 
			//boost::serialization::make_nvp("youngModulus", youngModulus) & 
			//boost::serialization::make_nvp("shearModulus", shearModulus) & 
			boost::serialization::make_nvp("density", density) & 
			//boost::serialization::make_nvp("extensionRatio", extensionRatio) & 
			//boost::serialization::make_nvp("shearingRatio", shearingRatio) & 
			boost::serialization::make_nvp("integrationTime", integrationTime) &
			boost::serialization::make_nvp("rodModel", rodModel) &
			boost::serialization::make_nvp("numNodes", numNodes);
	}
};

}	// namespace rod3d
}	// namespace qserl

#endif // QSERL_3D_PARAMETERS_H_
