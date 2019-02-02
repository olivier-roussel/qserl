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

/** Helper for analytic computation of mu values. 
* Just for experimentation, should not be used as long term.
*/


#ifndef QSERL_2D_ANALYTIC_Q_H_
#define QSERL_2D_ANALYTIC_Q_H_

#include "qserl/exports.h"

////#pragma warning( push, 0 )
//#include <Eigen/Lgsm>
////#pragma warning( pop )
#include <Eigen/Core>

namespace qserl {
namespace rod2d {

struct MotionConstantsQ
{
  //MotionConstantsMu muc;	// TODO factorize motion constants
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
  double gamma_0;
  double sn_gamma_0;
  double cn_gamma_0;
  double dn_gamma_0;
  double am_gamma_0;
  double E_am_gamma_0;
  double beta1_0;
  double beta2_0;

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
QSERL_EXPORT bool
computeMotionConstantsQ(const Eigen::Vector3d& i_a,
                        MotionConstantsQ& o_motionConstants);

/**
* \brief Compute rod geometry q (in body frame) at the position t along the rod.
* \param[in]  i_t Normalized position t along the rod between [0,1] for computing the q(t) values
* \param[in]  i_motionConstants Constants of motion for the rod
* \param[out] o_qdot Derivative of rod geometry dq(t) / dt where:
*             o_qdot[0] = dtheta / dt
*             o_qdot[1] = dx / dt
*             o_qdot[2] = dy / dt
* \param[out] o_q Rod geometry q(t) (at position t) where:
*             o_q[0] = theta
*             o_q[1] = x
*             o_q[2] = y
* \pre a(i) (i=3..5) values corresponding to given motion constants must respect the unhandled following cases:
*   - if Case I (includes Case III) i.e. lambda4 >= 0, then a3 != 0 and a5 != 0
*/
QSERL_EXPORT bool
computeQAtPositionT(double i_t,
                    const Eigen::Vector3d& i_a,
                    const MotionConstantsQ& i_mc,
                    Eigen::Vector3d& o_qdot,
                    Eigen::Vector3d& o_q);

}  // namespace rod2d
}  // namespace qserl

#endif // QSERL_2D_ANALYTIC_Q_H_