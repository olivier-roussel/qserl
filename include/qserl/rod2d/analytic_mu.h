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


#ifndef ANALYTIC_MU_H_
#define ANALYTIC_MU_H_

#include "qserl/exports.h"

#pragma warning( push, 0 )	
#include <Eigen/Lgsm>
#pragma warning( pop )	

namespace qserl {
namespace rod2d {

struct MotionConstantsMu
{
  double lambda[4];
  double delta;
  double alpha[3];
  double k; // defined as the modulus m = k^2, as Boost elliptic functions impl. uses k instead of m
  double m;
  double n;
  double r;
  double eta;
  double tau;
  signed char epsilon_tau;
  signed char epsilon_k;

  // precomputed constants 
  //double sqrt_alpha3; // = sqrt(alpha[2])
};

/**
* \brief Compute constants of motion for the rod from its parameterization in A-space.
* \param[in]  i_a Rod parameterization in A-space where:
*             a[0] is rod base torque
*             a[1] is rod base force along x
*             a[2] is rod base force along y
* \param[out] o_motionConstants Constants of motion for the rod
* \return false if a(i) (i=1..3) values corresponds to a unhandled special case (then motion constants 
*   will not be valid).
*/
QSERL_EXPORT bool computeMotionConstantsMu(const Eigen::Vector3d& i_a, MotionConstantsMu& o_motionConstants);

/**
* \brief Compute internal rod wrenches (in body frame) at the position t along the rod.
* \param[in]  i_t Normalized position t along the rod between [0,1] for computing the mu(t) values
* \param[in]  i_motionConstants Constants of motion for the rod
* \param[out] o_mu Internal rod wrenches (body frame) where:
*             o_mu[0] = k(t)                             (torque)
*             o_mu[1] = 0.5 * ( k(t)^2 + lambda2)        (force along x)
*             o_mu[2] = -k_dot(t)                        (force along y)
* \pre a(i) (i=3..5) values corresponding to given motion constants must respect the unhandled following cases:
*   - if Case I (includes Case III) i.e. lambda4 >= 0, then a3 != 0 and a5 != 0
*/
QSERL_EXPORT bool computeMuAtPositionT(double i_t, const MotionConstantsMu& i_motionConstants, Eigen::Vector3d& o_mu);

}	// namespace rod2d
}	// namespace qserl

#endif // ANALYTIC_MU_H_