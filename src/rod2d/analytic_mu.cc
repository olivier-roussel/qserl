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

/** Helper for computation of mu values. 
* Just for experimentation, should not be used as long term.
*/

#include "qserl/rod2d/analytic_mu.h"

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/acosh.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>

#include "util/utils.h"

namespace qserl {
namespace rod2d {

  
bool computeMotionConstantsMu(const Eigen::Vector3d& i_a, MotionConstantsMu& o_mc)
{
  static const double kEpsilonNullTorque = 1.e-10;
  static const double kEpsilonNullForce = 1.e-10;
  static const double kEpsilonLambda = 1.e-12;

  const double sqrd_a3 = util::sqr(i_a[0]);
  o_mc.lambda[0] = 0.; // lambda1
  o_mc.lambda[1] = sqrd_a3 + 2*i_a[1]; // lambda2
  o_mc.lambda[2] = 0.; // lambda3
  o_mc.lambda[3] = util::sqr(i_a[2]) - sqrd_a3*(0.25*sqrd_a3 + i_a[1]); // lambda4

  o_mc.epsilon_tau = util::sigpos(i_a[0]) * util::sigpos(i_a[2]);
  o_mc.epsilon_k = util::sigpos(i_a[0]);

  o_mc.delta = util::sqr(o_mc.lambda[1]) + 4*o_mc.lambda[3];
  const double sqrt_delta = sqrt(o_mc.delta);

  if (o_mc.lambda[3] >= 0.)
  {
    // unhandled special case corresponding to a3 = a5 = 0 
    // which is equivalent to lambda4 = 0 and lambda2 < 0 when a4 < 0 (Case III.2)
    if (abs(i_a[0]) < kEpsilonNullTorque && abs(i_a[2]) < kEpsilonNullForce)
      return false;

    // case I : lambda4 > 0 (also embeds case III where lambda4 == 0)
    if (!(o_mc.lambda[3] == 0. && o_mc.lambda[1] < 0.))
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
    }else if (o_mc.lambda[3] == 0. && o_mc.lambda[1] == 0){
      // lambda4 == 0 and lambda2 == 0 => all a_i are nulls
      o_mc.alpha[0] = 0.;
      o_mc.alpha[1] = 0.;
      o_mc.alpha[2] = 0.;

      o_mc.m = 0.; // undefined ?
      o_mc.k = 0.; // undefined ?
      o_mc.n = 0.; // undefined ?
      o_mc.r = 0.;
      o_mc.eta = 1.; // to be verified

      o_mc.epsilon_tau = 0;
      o_mc.epsilon_k = 0;

      o_mc.tau = 0.; // undefined ?
    }
  }else if (o_mc.lambda[3] < 0.){
    // case II : lambda4 < 0
    if (abs(i_a[1]) > kEpsilonNullForce || abs(i_a[2]) > kEpsilonNullForce)
    {
      // case II.1 : a4 != 0  and a5 != 0 (i.e. m != 0)
      o_mc.alpha[0] = 0.;
      o_mc.alpha[1] = o_mc.lambda[1] - sqrt_delta;
      o_mc.alpha[2] = o_mc.lambda[1] + sqrt_delta;

      const double sqrt_alpha3 = sqrt(o_mc.alpha[2]);

      o_mc.m = 2*sqrt_delta / (o_mc.lambda[1] + sqrt_delta);
      o_mc.k = sqrt(o_mc.m);
      o_mc.n = o_mc.m;
      o_mc.r = 0.5 * sqrt_alpha3;

      const double eta_sqrd = util::clamp((1. - sqrd_a3 / o_mc.alpha[2]) / o_mc.n, 0., 1.);
      o_mc.eta = sqrt(eta_sqrd);
      o_mc.tau = o_mc.epsilon_tau * boost::math::ellint_1(o_mc.k, asin(o_mc.eta)) / o_mc.r;
    }else{
      // case II.2 : a4 == 0  and a5 == 0 (i.e. m == 0)
      o_mc.alpha[0] = 0.;
      o_mc.alpha[1] = util::sqr(i_a[0]);
      o_mc.alpha[2] = o_mc.alpha[1];

      o_mc.m = 0.;
      o_mc.k = 0.;
      o_mc.n = 0.;
      o_mc.r = 0.5 * i_a[0];

      o_mc.eta = 0.;
      o_mc.tau = 0.;
    }
  }
  return true;
}


bool computeMuAtPositionT(double i_t, const MotionConstantsMu& i_mc, Eigen::Vector3d& o_mu)
{
  static const double kEpsilonLambda = 1.e-12;

  double k_t = 0.;      // k(t) is rod curvature at position t
  double k_dot_t = 0.;
  if (i_mc.lambda[3] >= 0.)
  {
    // unhandled special case corresponding to a3 = a5 = 0 
    // which is equivalent to lambda4 = 0 and lambda2 < 0 when a4 < 0 (Case III.2)
    //if (abs(i_a[0]) < kEpsilonNullTorque && abs(i_a[2]) < kEpsilonNullForce)
    //  return false;

    // case I : lambda4 > 0 (also embeds case III where lambda4 == 0)
    if (!(i_mc.lambda[3] == 0. && i_mc.lambda[1] < 0.))
    {
      const double sqrt_alpha3 = sqrt(i_mc.alpha[2]);
      const double gamma_t = i_mc.r*(i_t + i_mc.tau);
      // as we need all three elliptic jacobi functions sn, cn and dn, it is faster to compute the three together
      double cn_gamma_t, dn_gamma_t;
      const double sn_gamma_t = boost::math::jacobi_elliptic(i_mc.k, gamma_t, &cn_gamma_t, &dn_gamma_t); 
      k_t = i_mc.epsilon_k * sqrt_alpha3 * cn_gamma_t;
      k_dot_t = -i_mc.epsilon_k * i_mc.r * sqrt_alpha3 * sn_gamma_t * dn_gamma_t;
    }else if (i_mc.lambda[3] == 0. && i_mc.lambda[1] == 0){
      // lambda4 == 0 and lambda2 == 0 => all a_i are nulls
      // k_t = 0 and k_dot_t = 0 => no op
    }
  }else if (i_mc.lambda[3] < 0.){
    // case II : lambda4 < 0
    if (i_mc.k != 0.)
    {
      // case II.1 : a4 != 0  and a5 != 0 (i.e. m != 0)
      const double sqrt_alpha3 = sqrt(i_mc.alpha[2]);
      const double gamma_t = i_mc.r*(i_t + i_mc.tau);
      // as we need all three elliptic jacobi functions sn, cn and dn, it is faster to compute the three together
      double cn_gamma_t, dn_gamma_t;
      const double sn_gamma_t = boost::math::jacobi_elliptic(i_mc.k, gamma_t, &cn_gamma_t, &dn_gamma_t); 
      k_t = i_mc.epsilon_k * sqrt_alpha3 * dn_gamma_t;
      k_dot_t = -i_mc.epsilon_k * util::sqr(i_mc.k) * i_mc.alpha[2] * 0.5 * sn_gamma_t * cn_gamma_t;
    }else{
      // case II.2 : a4 == 0  and a5 == 0 (i.e. m == 0)
      k_t = 2. * i_mc.r; // == a3
      k_dot_t = 0.;
    }
  }

  o_mu[0] = k_t;
  o_mu[1] = 0.5 * (-util::sqr(k_t) + i_mc.lambda[1]);
  o_mu[2] = -k_dot_t;

  return true;
}

}	// namespace rod2d
}	// namespace qserl
