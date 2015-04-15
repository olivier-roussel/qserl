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

/** Helper for analytic computation of mu values. 
* Just for experimentation, should not be used as long term.
*/


#ifndef QSERL_ROD2_INVERSE_GEOMETRY_Q_H_
#define QSERL_ROD2_INVERSE_GEOMETRY_Q_H_

#include "qserl/exports.h"

#pragma warning( push, 0 )	
#include <unsupported/Eigen/Lgsm>
#pragma warning( pop )	

namespace qserl {
namespace rod2d {

/**
* \param[int] i_q_des Rod geometry q(1) (i.e. at position t=1) where:
*             i_q_des[0] = theta
*             i_q_des[1] = x
*             i_q_des[2] = y
* \param[int] i_a0 Initial guess of rod parameterization in A-space where:
*             i_a0[0] is rod base torque
*             i_a0[1] is rod base force along x
*             i_a0[2] is rod base force along y
* \param[in]  i_alpha scaling factor for the descent step.
* \param[out] o_a Resulting rod parameterization in A-space where:
*             a[0] is rod base torque
*             a[1] is rod base force along x
*             a[2] is rod base force along y
* \return true if successfull, false if unstability or singularity found.
*/
QSERL_EXPORT bool inverseGeometry_Newton(const Eigen::Vector3d& i_q_des, int i_maxIter, double i_maxNormError, 
  const Eigen::Vector3d& i_a0, double i_alpha, Eigen::Vector3d& o_a);

}	// namespace rod2d
}	// namespace qserl

#endif // QSERL_ROD2_INVERSE_GEOMETRY_Q_H_