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

/** Helper for analytic computation of da / da  values. 
* Just for experimentation, should not be used as long term.
*/


#ifndef QSERL_2D_ANALYTIC_DQDA_H_
#define QSERL_2D_ANALYTIC_DQDA_H_

#include "qserl/exports.h"

#include "qserl/rod2d/analytic_q.h"

#pragma warning( push, 0 )	
#include <Eigen/Lgsm>
#pragma warning( pop )	

namespace qserl {
namespace rod2d {

struct MotionConstantsDqDa
{
  // general constants
  MotionConstantsQ qc;

  // jacobian specific
  Eigen::Vector3d dlambda2_da;
  Eigen::Vector3d dlambda4_da;
  Eigen::Vector3d ddelta_da;
  Eigen::Matrix3d dalpha_da;
  Eigen::Vector3d dm_da;
  Eigen::Vector3d dn_da;
  Eigen::Vector3d dr_da;
  Eigen::Vector3d deta_da;
  Eigen::Vector3d dtau_da;
  Eigen::Vector3d dgamma_0_da;
  Eigen::Vector3d dsn_gamma_0_da;
  Eigen::Vector3d dcn_gamma_0_da;
  Eigen::Vector3d ddn_gamma_0_da;
  Eigen::Vector3d dE_am_gamma_0_da;
  Eigen::Vector3d dbeta1_0_da;
  Eigen::Vector3d dbeta2_0_da;

  // precomputed constants 
  //double sqrt_alpha3; // = sqrt(alpha[2])
  //double  inv_r;
};

/**
* \brief Compute constants of motion for the rod from its parameterization in A-space.
* \param[in]  i_a Rod parameterization in A-space where:
*             a[0] is rod base torque
*             a[1] is rod base force along x
*             a[2] is rod base force along y
* \param[out] o_motionConstants Constants of motion for the rod.
* \return false if a(i) (i=1..3) values corresponds to a unhandled special case (then motion constants 
*   will not be valid).
*/
QSERL_EXPORT bool computeMotionConstantsDqDa(const Eigen::Vector3d& i_a, MotionConstantsDqDa& o_mc);

/**
* \pre a(i) (i=3..5) values corresponding to given motion constants must respect the unhandled following cases:
*   - if Case I (includes Case III) i.e. lambda4 >= 0, then a3 != 0 and a5 != 0
*   - if Case II i.e. lambda4 < 0, then a5 != 0 (a3 cannot be null in this case)
*/
QSERL_EXPORT bool computeDqDaAtPositionT(double i_t, const MotionConstantsDqDa& i_mc, 
  Eigen::Matrix3d& o_dqda);

}	// namespace rod2d
}	// namespace qserl

#endif // QSERL_2D_ANALYTIC_DQDA_H_