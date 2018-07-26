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

#include "qserl/rod3d/parameters.h"

#include "qserl/util/constants.h"
#include "util/utils.h"

namespace qserl {
namespace rod3d {

/************************************************************************/
/*		 setIsotropicStiffnessCoefficientsFromElasticityParameters				*/
/************************************************************************/
void
Parameters::setIsotropicStiffnessCoefficientsFromElasticityParameters(double i_youngModulus,
                                                                      double i_shearModulus)
{
  const double sectionArea = util::sqr(radius) * constants::pi;
  const double I = 0.25 * util::sqr(util::sqr(radius)) * constants::pi;
  const double J = 0.5 * util::sqr(util::sqr(radius)) * constants::pi;

  stiffnessCoefficients[0] = i_shearModulus * J;
  stiffnessCoefficients[1] = i_youngModulus * I;
  stiffnessCoefficients[2] = stiffnessCoefficients[1];
  stiffnessCoefficients[3] = i_youngModulus * sectionArea;
  stiffnessCoefficients[4] = i_shearModulus * sectionArea;
  stiffnessCoefficients[5] = stiffnessCoefficients[4];
}

/************************************************************************/
/*												getIsotropicYoungModulus											*/
/************************************************************************/
double
Parameters::getIsotropicYoungModulus() const
{
  assert(stiffnessCoefficients[1] == stiffnessCoefficients[2] && stiffnessCoefficients[4] == stiffnessCoefficients[5] &&
         "rod stifness constants must represent (transversal) isotropy elasticity");

  const double sectionArea = util::sqr(radius) * constants::pi;
  return stiffnessCoefficients[3] / sectionArea;
}

/************************************************************************/
/*												getIsotropicShearModulus											*/
/************************************************************************/
double
Parameters::getIsotropicShearModulus() const
{
  assert(stiffnessCoefficients[1] == stiffnessCoefficients[2] && stiffnessCoefficients[4] == stiffnessCoefficients[5] &&
         "rod stifness constants must represent (transversal) isotropy elasticity");

  const double sectionArea = util::sqr(radius) * constants::pi;
  return stiffnessCoefficients[4] / sectionArea;
}

}  // namespace rod3d
}  // namespace qserl
