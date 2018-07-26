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

/** Helper for computation of mu values. 
* Just for experimentation, should not be used as long term.
*/

#include "qserl/rod2d/analytic_q.h"

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/acosh.hpp>
#include "util/jacobi_elliptic.h"

#include "util/utils.h"

namespace qserl {
namespace rod2d {


bool
computeMotionConstantsQ(const Eigen::Vector3d& i_a,
                        MotionConstantsQ& o_mc)
{
  static const double kEpsilonNullTorque = 1.e-10;
  static const double kEpsilonNullForce = 1.e-10;

  const double sqrd_a3 = util::sqr(i_a[0]);
  o_mc.lambda[0] = 0.; // lambda1
  o_mc.lambda[1] = sqrd_a3 + 2 * i_a[1]; // lambda2
  o_mc.lambda[2] = 0.; // lambda3
  o_mc.lambda[3] = util::sqr(i_a[2]) - sqrd_a3 * (0.25 * sqrd_a3 + i_a[1]); // lambda4

  o_mc.epsilon_tau = util::sigpos(i_a[0]) * util::sigpos(i_a[2]);
  o_mc.epsilon_k = util::sigpos(i_a[0]);

  o_mc.delta = util::sqr(o_mc.lambda[1]) + 4 * o_mc.lambda[3];
  const double sqrt_delta = sqrt(o_mc.delta);

  if(o_mc.lambda[3] >= 0.)
  {
    // unhandled special case corresponding to a3 = a5 = 0 
    // which is equivalent to lambda4 = 0 and lambda2 < 0 when a4 < 0 (Case III.2)
    if(abs(i_a[0]) < kEpsilonNullTorque && abs(i_a[2]) < kEpsilonNullForce)
    {
      return false;
    }

    // case I : lambda4 > 0 (also embeds case III where lambda4 == 0)
    if(!(o_mc.lambda[3] == 0. && o_mc.lambda[1] < 0.))
    {
      o_mc.alpha[0] = -(o_mc.lambda[1] - sqrt_delta);
      o_mc.alpha[1] = 0.;
      o_mc.alpha[2] = o_mc.lambda[1] + sqrt_delta;

      o_mc.m = 0.5 + (o_mc.lambda[1] * 0.5 / sqrt_delta);
      o_mc.k = sqrt(o_mc.m);
      o_mc.n = 1.;
      o_mc.r = sqrt(sqrt_delta * 0.5);

      const double eta_sqrd = util::clamp(1. - sqrd_a3 / o_mc.alpha[2], 0., 1.);
      o_mc.eta = sqrt(eta_sqrd);
      o_mc.tau = o_mc.epsilon_tau * boost::math::ellint_1(o_mc.k, asin(o_mc.eta)) / o_mc.r;

      o_mc.gamma_0 = o_mc.r * o_mc.tau;
      o_mc.sn_gamma_0 = boost::math::jacobi_elliptic(o_mc.k, o_mc.gamma_0, &o_mc.cn_gamma_0,
                                                     &o_mc.dn_gamma_0, &o_mc.am_gamma_0);
      o_mc.E_am_gamma_0 = boost::math::ellint_2(o_mc.k, o_mc.am_gamma_0);

      o_mc.beta1_0 = 2. * util::sqr(o_mc.dn_gamma_0) - 1.;
      o_mc.beta2_0 = o_mc.sn_gamma_0 * o_mc.dn_gamma_0;

    }
    else if(o_mc.lambda[3] == 0. && o_mc.lambda[1] == 0)
    {
      // lambda4 == 0 and lambda2 == 0 => all a_i are nulls
      assert (abs(i_a[0]) < kEpsilonNullTorque);
      assert (abs(i_a[1]) < kEpsilonNullForce);
      assert (abs(i_a[2]) < kEpsilonNullForce);

      o_mc.alpha[0] = 0.;
      o_mc.alpha[1] = 0.;
      o_mc.alpha[2] = 0.;

      o_mc.k = 0.; // undefined ?
      o_mc.n = 0.; // undefined ?
      o_mc.r = 0.;
      o_mc.eta = 1.; // to be verified
      o_mc.epsilon_tau = 0;
      o_mc.epsilon_k = 0;
      o_mc.tau = 0.; // undefined ?
      o_mc.gamma_0 = 0.;
      o_mc.sn_gamma_0 = 0.;
      o_mc.cn_gamma_0 = 1.;
      o_mc.dn_gamma_0 = 1.;
      o_mc.am_gamma_0 = 0.;
      o_mc.E_am_gamma_0 = 0.;
      o_mc.beta1_0 = 1.;
      o_mc.beta2_0 = 0.;
    }
  }
  else if(o_mc.lambda[3] < 0.)
  {
    // case II : lambda4 < 0
    if(abs(i_a[1]) > kEpsilonNullForce || abs(i_a[2]) > kEpsilonNullForce)
    {
      // case II.1 : a4 != 0  and a5 != 0 (i.e. m != 0)
      o_mc.alpha[0] = 0.;
      o_mc.alpha[1] = o_mc.lambda[1] - sqrt_delta;
      o_mc.alpha[2] = o_mc.lambda[1] + sqrt_delta;

      const double sqrt_alpha3 = sqrt(o_mc.alpha[2]);

      o_mc.m = 2 * sqrt_delta / (o_mc.lambda[1] + sqrt_delta);
      o_mc.k = sqrt(o_mc.m);
      o_mc.n = o_mc.m;
      o_mc.r = 0.5 * sqrt_alpha3;

      const double eta_sqrd = util::clamp((1. - sqrd_a3 / o_mc.alpha[2]) / o_mc.n, 0., 1.);
      o_mc.eta = sqrt(eta_sqrd);
      o_mc.tau = o_mc.epsilon_tau * boost::math::ellint_1(o_mc.k, asin(o_mc.eta)) / o_mc.r;

      o_mc.gamma_0 = o_mc.r * o_mc.tau;
      o_mc.sn_gamma_0 = boost::math::jacobi_elliptic(o_mc.k, o_mc.gamma_0, &o_mc.cn_gamma_0,
                                                     &o_mc.dn_gamma_0, &o_mc.am_gamma_0);
      o_mc.E_am_gamma_0 = boost::math::ellint_2(o_mc.k, o_mc.am_gamma_0);

      o_mc.beta1_0 = 1. - 2 * util::sqr(o_mc.sn_gamma_0);
      o_mc.beta2_0 = o_mc.cn_gamma_0 * o_mc.sn_gamma_0;;
    }
    else
    {
      // case II.2 : a4 == 0  and a5 == 0 (i.e. m == 0) and a3 != 0
      assert (abs(i_a[0]) > kEpsilonNullTorque);
      assert (abs(i_a[1]) < kEpsilonNullForce);
      assert (abs(i_a[2]) < kEpsilonNullForce);

      o_mc.alpha[0] = 0.;
      o_mc.alpha[1] = util::sqr(i_a[0]);
      o_mc.alpha[2] = o_mc.alpha[1];

      o_mc.k = 0;
      o_mc.n = 0;
      o_mc.r = 0.5 * i_a[0];
      o_mc.eta = 0.;
      o_mc.tau = 0.;
      o_mc.gamma_0 = 0.;
      o_mc.sn_gamma_0 = 0.;
      o_mc.cn_gamma_0 = 1.;
      o_mc.dn_gamma_0 = 1.;
      o_mc.am_gamma_0 = 0.;
      o_mc.E_am_gamma_0 = 0.;
      o_mc.beta1_0 = 1.;
      o_mc.beta2_0 = 0.;
    }
  }
  return true;
}


bool
computeQAtPositionT(double i_t,
                    const Eigen::Vector3d& i_a,
                    const MotionConstantsQ& i_mc,
                    Eigen::Vector3d& o_qdot,
                    Eigen::Vector3d& o_q)
{
  static const double kEpsilonNullTorque = 1.e-10;
  static const double kEpsilonNullForce = 1.e-10;

  if(i_mc.lambda[3] >= 0.)
  {
    // unhandled special case corresponding to a3 = a5 = 0 
    // which is equivalent to lambda4 = 0 and lambda2 < 0 when a4 < 0 (Case III.2)
    //if (abs(i_a[0]) < kEpsilonNullTorque && abs(i_a[2]) < kEpsilonNullForce)
    //  return false;

    // case I : lambda4 > 0 (also embeds case III where lambda4 == 0)
    if(!(i_mc.lambda[3] == 0. && i_mc.lambda[1] < 0.))
    {
      // pre-computed constants
      const double sqrt_alpha3 = sqrt(i_mc.alpha[2]);
      const double inv_r = 1. / i_mc.r;

      const double gamma_t = i_mc.r * (i_t + i_mc.tau);
      double cn_gamma_t, dn_gamma_t, am_gamma_t;
      const double sn_gamma_t = boost::math::jacobi_elliptic(i_mc.k, gamma_t, &cn_gamma_t,
                                                             &dn_gamma_t, &am_gamma_t);
      const double E_am_gamma_t = boost::math::ellint_2(i_mc.k, am_gamma_t);

      const double beta1_t = 2 * util::sqr(dn_gamma_t) - 1.;
      const double beta2_t = sn_gamma_t * dn_gamma_t;
      const double int_beta1 = (2 * inv_r) * (E_am_gamma_t - i_mc.E_am_gamma_0) - i_t;
      const double int_beta2 = -inv_r * (cn_gamma_t - i_mc.cn_gamma_0);

      o_qdot[0] = i_mc.epsilon_k * sqrt_alpha3 * cn_gamma_t;
      o_qdot[1] = i_mc.beta1_0 * beta1_t + 4 * i_mc.m * i_mc.beta2_0 * beta2_t;
      o_qdot[2] = 2 * i_mc.epsilon_k * i_mc.k * (i_mc.beta1_0 * beta2_t - beta1_t * i_mc.beta2_0);

      //o_q[0] = i_mc.epsilon_k * 2 * (acos(dn_gamma_t) - acos(i_mc.dn_gamma_0)); WRONG !
      o_q[0] = atan2(o_qdot[2], o_qdot[1]);
      o_q[1] = i_mc.beta1_0 * int_beta1 + 4 * i_mc.m * i_mc.beta2_0 * int_beta2;
      o_q[2] = 2 * i_mc.epsilon_k * i_mc.k * (i_mc.beta1_0 * int_beta2 - i_mc.beta2_0 * int_beta1);

    }
    else if(i_mc.lambda[3] == 0. && i_mc.lambda[1] == 0)
    {
      // lambda4 == 0 and lambda2 == 0 => all a_i are nulls
      (void) kEpsilonNullTorque;
      assert (abs(i_a[0]) < kEpsilonNullTorque);
      assert (abs(i_a[1]) < kEpsilonNullForce);
      assert (abs(i_a[2]) < kEpsilonNullForce);

      // straight line configuration 
      o_qdot[0] = 0.;
      o_qdot[1] = 1.;
      o_qdot[2] = 0.;

      o_q[0] = 0.;
      o_q[1] = i_t;
      o_q[2] = 0.;
    }
  }
  else if(i_mc.lambda[3] < 0.)
  {
    // case II : lambda4 < 0
    if(abs(i_a[1]) > kEpsilonNullForce || abs(i_a[2]) > kEpsilonNullForce)
    {
      // case II.1 : a4 != 0  and a5 != 0 (i.e. m != 0)
      assert (i_mc.k > 0. && i_mc.m > 0.);

      // pre-computed constants
      const double sqrt_alpha3 = sqrt(i_mc.alpha[2]);
      const double inv_r = 1. / i_mc.r;
      const double inv_m = 1. / i_mc.m;

      const double gamma_t = i_mc.r * (i_t + i_mc.tau);
      double cn_gamma_t, dn_gamma_t, am_gamma_t;
      const double sn_gamma_t = boost::math::jacobi_elliptic(i_mc.k, gamma_t, &cn_gamma_t,
                                                             &dn_gamma_t, &am_gamma_t);
      const double E_am_gamma_t = boost::math::ellint_2(i_mc.k, am_gamma_t);

      const double beta1_t = 1. - 2 * util::sqr(sn_gamma_t);
      const double beta2_t = cn_gamma_t * sn_gamma_t;
      const double int_beta1 = inv_m * (i_t * (i_mc.m - 2.) + 2 * inv_r *
                                                              (E_am_gamma_t - i_mc.E_am_gamma_0));
      const double int_beta2 = -inv_m * inv_r * (dn_gamma_t - i_mc.dn_gamma_0);

      o_qdot[0] = i_mc.epsilon_k * sqrt_alpha3 * dn_gamma_t;
      o_qdot[1] = i_mc.beta1_0 * beta1_t + 4 * beta2_t * i_mc.beta2_0;
      o_qdot[2] = 2 * i_mc.epsilon_k * (beta2_t * i_mc.beta1_0 - i_mc.beta2_0 * beta1_t);

      // o_q[0] = i_mc.epsilon_k * 2 * (asin(sn_gamma_t) - asin(i_mc.sn_gamma_0));
      o_q[0] = i_mc.epsilon_k * 2 * (am_gamma_t - i_mc.am_gamma_0);
      o_q[1] = i_mc.beta1_0 * int_beta1 + 4 * i_mc.beta2_0 * int_beta2;
      o_q[2] = 2 * i_mc.epsilon_k * (i_mc.beta1_0 * int_beta2 - i_mc.beta2_0 * int_beta1);
    }
    else
    {
      // case II.2 : a4 == 0  and a5 == 0 (i.e. m == 0) and a3 != 0
      assert (abs(i_a[0]) > kEpsilonNullTorque);
      assert (abs(i_a[1]) < kEpsilonNullForce);
      assert (abs(i_a[2]) < kEpsilonNullForce);

      const double inv_a3 = 1. / i_a[2];

      o_qdot[0] = i_a[2];
      o_qdot[1] = cos(i_a[2] * i_t);
      o_qdot[2] = sin(i_a[2] * i_t);

      o_q[0] = i_a[2] * i_t;
      o_q[1] = inv_a3 * o_qdot[2];
      o_q[2] = -inv_a3 * (o_qdot[1] - 1.);
    }
  }

  return true;
}

}  // namespace rod2d
}  // namespace qserl
