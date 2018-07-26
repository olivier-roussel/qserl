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

#ifndef QSERL_3D_PARAMETERS_H_
#define QSERL_3D_PARAMETERS_H_

#include "qserl/exports.h"

#include <Eigen/Core>
#include <cereal/access.hpp>
#include <cereal/archives/xml.hpp>

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
      stiffnessCoefficients(Eigen::Matrix<double, 6, 1>::Ones()),
      //density(1.1 * 10e3),   /** 1.10 kg/dm3 -> kg/m3, */
      rodModel(RM_EXTENSIBLE_SHEARABLE),
      numNodes(100),
      integrationTime(1.)
  {
  }

  enum RodModelT
  {
    RM_INEXTENSIBLE = 0,
    RM_EXTENSIBLE_SHEARABLE,
    RM_NUMBER_OF_ROD_MODELS
  };

  static const std::string
  getRodModelName(RodModelT i_model)
  {
    const static char* const model_names_array[] = {"INEXTENSIBLE", "EXTENSIBLE_SHEARABLE"};
    const static std::vector<std::string> v_model_names_array(model_names_array,
                                                              model_names_array + RM_NUMBER_OF_ROD_MODELS);
    return v_model_names_array[static_cast<int>(i_model)];
  }

  /**
  * \brief Returns equivalent Young modulus from stiffness coefficients.
  * \pre Stiffness coefficients must represent transversal isotropic elasticity, i.e.
  *   stiffnessCoefficients[1] == stiffnessCoefficients[2]    and
  *   stiffnessCoefficients[4] == stiffnessCoefficients[5]
  */
  double
  getIsotropicYoungModulus() const;

  /**
  * \brief Returns equivalent Shear modulus from stiffness coefficients.
  * \pre Stiffness coefficients must represent transversal isotropic elasticity, i.e.
  *   stiffnessCoefficients[1] == stiffnessCoefficients[2]    and
  *   stiffnessCoefficients[4] == stiffnessCoefficients[5]
  */
  double
  getIsotropicShearModulus() const;

  /**
  * \brief Set transversal isotropic stiffness coefficients from elasticity parameters described by
  * Young modulus and shear modulus. 
  */
  void
  setIsotropicStiffnessCoefficientsFromElasticityParameters(double i_youngModulus,
                                                            double i_shearModulus);
  /**
   * Attributes
   */
  double                        radius;
  double                        length;
  Eigen::Matrix<double, 6, 1>   stiffnessCoefficients;
  RodModelT                     rodModel;
  int                           numNodes;       /** Number of discretization nodes. Related to the delta_t field
                                                  * used for 2D rod by delta_t = 1 / (numNodes - 1) */
  /**< Internal use only. */
  double                        integrationTime;  /**< Should be kept to 1 (default value). */

  /**
   * Serialization
   */
  template<class Archive>
  void
  serialize(Archive& ar,
            const unsigned int version)
  {
    ar(cereal::make_nvp("radius", radius),
       cereal::make_nvp("length", length),
       cereal::make_nvp("stiffnessCoefficients", stiffnessCoefficients),
       cereal::make_nvp("rodModel", rodModel),
       cereal::make_nvp("numNodes", numNodes),
        //cereal::make_nvp("density", density) &
       cereal::make_nvp("integrationTime", integrationTime)
    );
  }
};

}  // namespace rod3d
}  // namespace qserl

#endif // QSERL_3D_PARAMETERS_H_
