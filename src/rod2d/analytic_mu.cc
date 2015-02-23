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

  
bool computeMotionConstants(const Eigen::Vector3d& i_a, MotionConstants& o_motionConstants)
{
  static const double kEpsilonNullTorque = 1.e-10;
  static const double kEpsilonNullForce = 1.e-10;
  static const double kEpsilonLambda = 1.e-12;

  const double sqrd_a3 = util::sqr(i_a[0]);
  o_motionConstants.lambda[0] = 0.; // lambda1
  o_motionConstants.lambda[1] = sqrd_a3 + 2*i_a[1]; // lambda2
  o_motionConstants.lambda[2] = 0.; // lambda3
  o_motionConstants.lambda[3] = util::sqr(i_a[2]) - sqrd_a3*(0.25*sqrd_a3 + i_a[1]); // lambda4

  o_motionConstants.epsilon_tau = util::sigpos(i_a[0]) * util::sigpos(i_a[2]);
  o_motionConstants.epsilon_k = util::sigpos(i_a[0]);

  o_motionConstants.delta = util::sqr(o_motionConstants.lambda[1]) + 4*o_motionConstants.lambda[3];
  const double sqrt_delta = sqrt(o_motionConstants.delta);

  if (o_motionConstants.lambda[3] >= 0.)
  {
    // case I : lambda4 > 0 (also embeds case III where lambda4 == 0)
    if (!(o_motionConstants.lambda[3] == 0. && o_motionConstants.lambda[1] < 0.))
    {
      o_motionConstants.alpha[0] = -(o_motionConstants.lambda[1] - sqrt_delta);
      o_motionConstants.alpha[1] = 0.;
      o_motionConstants.alpha[2] = o_motionConstants.lambda[1] + sqrt_delta;

      o_motionConstants.k = sqrt(0.5 + (o_motionConstants.lambda[1] * 0.5 / sqrt_delta));
      //o_motionConstants.m = 0.5 + (o_motionConstants.lambda[1] * 0.5 / sqrt_delta);
      o_motionConstants.n = 1.;
      o_motionConstants.r = sqrt(sqrt_delta * 0.5);

      const double eta_sqrd = util::clamp(1. - sqrd_a3 / o_motionConstants.alpha[2], 0., 1.);
      o_motionConstants.eta = sqrt(eta_sqrd);
      o_motionConstants.tau = o_motionConstants.epsilon_tau * 
        boost::math::ellint_1(o_motionConstants.k, asin(o_motionConstants.eta)) / o_motionConstants.r;
    }else if (o_motionConstants.lambda[3] == 0. && o_motionConstants.lambda[1] == 0){
      // lambda4 == 0 and lambda2 == 0 => all a_i are nulls
      o_motionConstants.alpha[0] = 0.;
      o_motionConstants.alpha[1] = 0.;
      o_motionConstants.alpha[2] = 0.;

      o_motionConstants.k = 0.; // undefined ?
      o_motionConstants.n = 0.; // undefined ?
      o_motionConstants.r = 0.;
      o_motionConstants.eta = 1.; // to be verified

      o_motionConstants.epsilon_tau = 0;
      o_motionConstants.epsilon_k = 0;

      o_motionConstants.tau = 0.; // undefined ?
    }else{
      // unhandled case corresponding to a3 == a5 == 0 (abnormal case)
      return false;
    }
  }else if (o_motionConstants.lambda[3] < 0.){
    // case II : lambda4 < 0
    if (abs(i_a[1]) > kEpsilonNullForce || abs(i_a[2]) > kEpsilonNullForce)
    {
      // case II.1 : a4 != 0  and a5 != 0 (i.e. m != 0)
      o_motionConstants.alpha[0] = 0.;
      o_motionConstants.alpha[1] = o_motionConstants.lambda[1] - sqrt_delta;
      o_motionConstants.alpha[2] = o_motionConstants.lambda[1] + sqrt_delta;

      const double sqrt_alpha3 = sqrt(o_motionConstants.alpha[2]);

      o_motionConstants.k = sqrt(2*sqrt_delta / (o_motionConstants.lambda[1] + sqrt_delta));
      o_motionConstants.n = util::sqr(o_motionConstants.k);
      o_motionConstants.r = 0.5 * sqrt_alpha3;

      const double eta_sqrd = util::clamp((1. - sqrd_a3 / o_motionConstants.alpha[2]) / o_motionConstants.n, 0., 1.);
      o_motionConstants.eta = sqrt(eta_sqrd);
      o_motionConstants.tau = o_motionConstants.epsilon_tau * 
        boost::math::ellint_1(o_motionConstants.k, asin(o_motionConstants.eta)) / o_motionConstants.r;
    }else{
      // case II.2 : a4 == 0  and a5 == 0 (i.e. m == 0)
      o_motionConstants.alpha[0] = 0.;
      o_motionConstants.alpha[1] = util::sqr(i_a[0]);
      o_motionConstants.alpha[2] = o_motionConstants.alpha[1];

      o_motionConstants.k = 0;
      o_motionConstants.n = 0;
      o_motionConstants.r = 0.5 * i_a[0];

      o_motionConstants.eta = 0.;
      o_motionConstants.tau = 0.;
    }
  }
  return true;
}


bool computeMuAtPositionT(double i_t, const MotionConstants& i_motionConstants, Eigen::Vector3d& o_mu)
{
  static const double kEpsilonLambda = 1.e-12;

  double k_t = 0.;      // k(t) is rod curvature at position t
  double k_dot_t = 0.;
  if (i_motionConstants.lambda[3] >= 0.)
  {
    // case I : lambda4 > 0 (also embeds case III where lambda4 == 0)
    if (!(i_motionConstants.lambda[3] == 0. && i_motionConstants.lambda[1] < 0.))
    {
      const double sqrt_alpha3 = sqrt(i_motionConstants.alpha[2]);
      const double r_t_tau = i_motionConstants.r*(i_t + i_motionConstants.tau);
      // as we need all three elliptic jacobi functions sn, cn and dn, it is faster to compute the three together
      double cn_rttau, dn_rttau;
      const double sn_rttau = boost::math::jacobi_elliptic(i_motionConstants.k, r_t_tau, &cn_rttau, &dn_rttau); 
      k_t = i_motionConstants.epsilon_k * sqrt_alpha3 * cn_rttau;
      k_dot_t = -i_motionConstants.epsilon_k * i_motionConstants.r * sqrt_alpha3 * 
        sn_rttau * dn_rttau;
    }else if (i_motionConstants.lambda[3] == 0. && i_motionConstants.lambda[1] == 0){
      // lambda4 == 0 and lambda2 == 0 => all a_i are nulls
      // k_t = 0 and k_dot_t = 0 => no op
    }else{
      // unhandled case corresponding to a3 == a5 == 0 (abnormal case)
      return false;
    }
  }else if (i_motionConstants.lambda[3] < 0.){
    // case II : lambda4 < 0
    if (i_motionConstants.k != 0.)
    {
      // case II.1 : a4 != 0  and a5 != 0 (i.e. m != 0)
      const double sqrt_alpha3 = sqrt(i_motionConstants.alpha[2]);
      const double r_t_tau = i_motionConstants.r*(i_t + i_motionConstants.tau);
      // as we need all three elliptic jacobi functions sn, cn and dn, it is faster to compute the three together
      double cn_rttau, dn_rttau;
      const double sn_rttau = boost::math::jacobi_elliptic(i_motionConstants.k, r_t_tau, &cn_rttau, &dn_rttau); 
      k_t = i_motionConstants.epsilon_k * sqrt_alpha3 * dn_rttau;
      k_dot_t = -i_motionConstants.epsilon_k * util::sqr(i_motionConstants.k) * i_motionConstants.alpha[2] * 0.5 *
        sn_rttau * cn_rttau;
    }else{
      // case II.2 : a4 == 0  and a5 == 0 (i.e. m == 0)
      k_t = 2. * i_motionConstants.r; // == a3
      k_dot_t = 0.;
    }
  }

  o_mu[0] = k_t;
  o_mu[1] = 0.5 * (-util::sqr(k_t) + i_motionConstants.lambda[1]);
  o_mu[2] = -k_dot_t;

  return true;
}

void computeMotionConstants_old(const Eigen::Vector3d& i_a, MotionConstants& o_motionConstants)
{
  static const double kEpsilonNullTorque = 1.e-10;
  static const double kEpsilonNullForce = 1.e-10;
  static const double kEpsilonLambda = 1.e-12;

  const double sqrd_a3 = util::sqr(i_a[0]);
  o_motionConstants.lambda[0] = 0.; // lambda1
  o_motionConstants.lambda[1] = sqrd_a3 + 2*i_a[1]; // lambda2
  o_motionConstants.lambda[2] = 0.; // lambda3
  o_motionConstants.lambda[3] = util::sqr(i_a[2]) - sqrd_a3*(0.25*sqrd_a3 + i_a[1]); // lambda4

  o_motionConstants.epsilon_tau = util::sigpos(i_a[0]) * util::sigpos(i_a[2]);
  o_motionConstants.epsilon_k = util::sigpos(i_a[0]);

  o_motionConstants.delta = util::sqr(o_motionConstants.lambda[1]) + 4*o_motionConstants.lambda[3];
  const double sqrt_delta = sqrt(o_motionConstants.delta);

  if (o_motionConstants.lambda[3] > kEpsilonLambda)
  {
    // case I : lambda4 > 0
    o_motionConstants.alpha[0] = -(o_motionConstants.lambda[1] - sqrt_delta);
    o_motionConstants.alpha[1] = 0.;
    o_motionConstants.alpha[2] = o_motionConstants.lambda[1] + sqrt_delta;

    o_motionConstants.k = sqrt(0.5 + (o_motionConstants.lambda[1] * 0.5 / sqrt_delta));
    //o_motionConstants.m = 0.5 + (o_motionConstants.lambda[1] * 0.5 / sqrt_delta);
    o_motionConstants.n = 1.;
    o_motionConstants.r = sqrt(sqrt_delta * 0.5);

    o_motionConstants.tau = o_motionConstants.epsilon_tau * 
      boost::math::ellint_1(o_motionConstants.k, asin(sqrt(1. - sqrd_a3 / o_motionConstants.alpha[2])))
      / o_motionConstants.r;

  }else if (o_motionConstants.lambda[3] < -kEpsilonLambda){
    // case II : lambda4 < 0
    if (abs(i_a[1]) > kEpsilonNullForce || abs(i_a[2]) > kEpsilonNullForce)
    {
      // case II.1 : a4 != 0  and a5 != 0 (i.e. m != 0)
      o_motionConstants.alpha[0] = 0.;
      o_motionConstants.alpha[1] = o_motionConstants.lambda[1] - sqrt_delta;
      o_motionConstants.alpha[2] = o_motionConstants.lambda[1] + sqrt_delta;

      const double sqrt_alpha3 = sqrt(o_motionConstants.alpha[2]);

      o_motionConstants.k = sqrt(2*sqrt_delta / (o_motionConstants.lambda[1] + sqrt_delta));
      o_motionConstants.n = util::sqr(o_motionConstants.k);
      o_motionConstants.r = 0.5 * sqrt_alpha3;

      o_motionConstants.epsilon_tau = util::sigpos(i_a[0] * i_a[2]);
      o_motionConstants.epsilon_k = util::sigpos(i_a[0]);

      o_motionConstants.tau = o_motionConstants.epsilon_tau * 
        boost::math::ellint_1(o_motionConstants.k, asin(sqrt((1. - sqrd_a3 / 
        o_motionConstants.alpha[2])/o_motionConstants.n))) / o_motionConstants.r;
    }else{
      // case II.2 : a4 == 0  and a5 == 0 (i.e. m == 0)
      o_motionConstants.alpha[0] = 0.;
      o_motionConstants.alpha[1] = util::sqr(i_a[0]);
      o_motionConstants.alpha[2] = o_motionConstants.alpha[1];

      o_motionConstants.k = 0;
      o_motionConstants.n = 0;
      o_motionConstants.r = 0.5 * i_a[0];

      o_motionConstants.tau = 0.;
    }
  }else{
    // case III : lambda4 = 0
    if (o_motionConstants.lambda[1] > kEpsilonLambda)
    {
      // III.1 : lambda2 > 0
      o_motionConstants.alpha[0] = 0.;
      o_motionConstants.alpha[1] = 0.;
      o_motionConstants.alpha[2] = 2*o_motionConstants.lambda[1];

      const double sqrt_alpha3 = sqrt(o_motionConstants.alpha[2]);

      o_motionConstants.k = 1.;
      o_motionConstants.n = 1.;
      o_motionConstants.r = 0.5 * sqrt_alpha3;

      o_motionConstants.epsilon_tau = util::sigpos(i_a[0] * i_a[2]);
      o_motionConstants.epsilon_k = util::sigpos(i_a[0]);

      // a3 != 0
      if (abs(i_a[0]) > kEpsilonNullTorque)
      {
        o_motionConstants.tau = /*o_motionConstants.epsilon_tau * */
          boost::math::acosh(o_motionConstants.epsilon_k * sqrt_alpha3 / i_a[0]) / o_motionConstants.r;
      }else{
        o_motionConstants.tau = std::numeric_limits<double>::infinity();
      }
    }else if(o_motionConstants.lambda[1] < kEpsilonLambda){
      // III.2 : lambda2 < 0
      o_motionConstants.alpha[0] = 0.;
      o_motionConstants.alpha[1] = 0.;
      o_motionConstants.alpha[2] = -2*o_motionConstants.lambda[1];

      const double sqrt_alpha3 = sqrt(o_motionConstants.alpha[2]);

      o_motionConstants.k = 1.;
      o_motionConstants.n = 1.;
      o_motionConstants.r = 0.5 * sqrt_alpha3;

      o_motionConstants.epsilon_tau = util::sigpos(i_a[0] * i_a[2]);
      o_motionConstants.epsilon_k = util::sigpos(i_a[0]);

        // a3 != 0
      if (abs(i_a[0]) > kEpsilonNullTorque)
      {
        o_motionConstants.tau = /*o_motionConstants.epsilon_tau * */
          boost::math::acosh(o_motionConstants.epsilon_k * sqrt_alpha3 / i_a[0]) / o_motionConstants.r;
      }else{
        o_motionConstants.tau = std::numeric_limits<double>::infinity();

      }
    }else{
      // III.3 : lambda2 = 0 ( a(i) = 0 and P(y) has triple root alpha = 0 )
      o_motionConstants.alpha[0] = 0.;
      o_motionConstants.alpha[1] = 0.;
      o_motionConstants.alpha[2] = 0.;

      o_motionConstants.k = 0.; // undefined ?
      o_motionConstants.n = 0.; // undefined ?
      o_motionConstants.r = 0.;

      o_motionConstants.epsilon_tau = 0;
      o_motionConstants.epsilon_k = 0;

      o_motionConstants.tau = 0.; // undefined ?
    }
  }
}

void computeMuAtPositionT_old(double i_t, const MotionConstants& i_motionConstants, Eigen::Vector3d& o_mu)
{
  static const double kEpsilonLambda = 1.e-12;

  double k_t = 0.;      // k(t) is rod curvature at position t
  double k_dot_t = 0.;
  if (i_motionConstants.lambda[3] > kEpsilonLambda)
  {
    // case I : lambda4 > 0
    const double sqrt_alpha3 = sqrt(i_motionConstants.alpha[2]);
    const double r_t_tau = i_motionConstants.r*(i_t + i_motionConstants.tau);
    // as we need all three elliptic jacobi functions sn, cn and dn, it is faster to compute the three together
    double cn_rttau, dn_rttau;
    const double sn_rttau = boost::math::jacobi_elliptic(i_motionConstants.k, r_t_tau, &cn_rttau, &dn_rttau); 
    k_t = i_motionConstants.epsilon_k * sqrt_alpha3 * cn_rttau;
    k_dot_t = -i_motionConstants.epsilon_k * i_motionConstants.r * sqrt_alpha3 * 
      sn_rttau * dn_rttau;
  }else if (i_motionConstants.lambda[3] < kEpsilonLambda){
    // case II : lambda4 < 0
    if (i_motionConstants.k != 0.)
    {
      // case II.1 : a4 != 0  and a5 != 0 (i.e. m != 0)
      const double sqrt_alpha3 = sqrt(i_motionConstants.alpha[2]);
      const double r_t_tau = i_motionConstants.r*(i_t + i_motionConstants.tau);
      // as we need all three elliptic jacobi functions sn, cn and dn, it is faster to compute the three together
      double cn_rttau, dn_rttau;
      const double sn_rttau = boost::math::jacobi_elliptic(i_motionConstants.k, r_t_tau, &cn_rttau, &dn_rttau); 
      k_t = i_motionConstants.epsilon_k * sqrt_alpha3 * dn_rttau;
      k_dot_t = -i_motionConstants.epsilon_k * util::sqr(i_motionConstants.k) * i_motionConstants.alpha[2] * 0.5 *
        sn_rttau * cn_rttau;
    }else{
      // case II.2 : a4 == 0  and a5 == 0 (i.e. m == 0)
      k_t = 2. * i_motionConstants.r; // == a3
      k_dot_t = 0.;
    }
  }else{
    // case III : lambda4 = 0
    if (abs(i_motionConstants.lambda[1]) > kEpsilonLambda)
    {
      // III.1 : lambda2 > 0 and // III.2 : lambda2 < 0
      // a3 != 0
      if (i_motionConstants.tau == std::numeric_limits<double>::infinity())
      //if (abs(i_a[2]) > kEpsilonNullTorque)
      {
        const double sqrt_alpha3 = sqrt(i_motionConstants.alpha[2]);
        const double r_t_tau = i_motionConstants.r*(i_t + i_motionConstants.tau);
        const double cosh_rttau = cosh(r_t_tau);
        k_t = i_motionConstants.epsilon_k * sqrt_alpha3 / cosh_rttau;
        k_dot_t = -i_motionConstants.epsilon_k * i_motionConstants.lambda[1] / (cosh_rttau * tanh(r_t_tau));
      }
      /*else{
        // a3 == 0
        // no op _ k_t = k_dot_t = 0
        k_t = 0.;
        k_dot_t = 0.;
      }*/
    }else{
      // III.3 : lambda2 = 0 ( a(i) = 0 and P(y) has triple root alpha = 0 )
      // no op _ k_t = k_dot_t = 0 (no curvature _ straight line configuration)
    }
  }

  o_mu[0] = k_t;
  o_mu[1] = 0.5 * (-util::sqr(k_t) + i_motionConstants.lambda[1]);
  o_mu[2] = -k_dot_t;
}

}	// namespace rod2d
}	// namespace qserl
