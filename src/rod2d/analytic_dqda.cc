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

/** Helper for computation of dq / da values. 
* Just for experimentation, should not be used as long term.
*/

#include "qserl/rod2d/analytic_dqda.h"

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/acosh.hpp>
#include "util/jacobi_elliptic.h"
#include "util/utils.h"

namespace qserl {
namespace rod2d {

bool computeMotionConstantsDqDa(const Eigen::Vector3d& i_a, MotionConstantsDqDa& o_mc)
{
  static const double kEpsilonNullTorque = 1.e-10;
  static const double kEpsilonNullForce = 1.e-10;
  static const double kEpsilonLambda = 1.e-12;

  static const double inv_sqrt_2 = 1. / sqrt(2.);

  const double sqrd_a3 = util::sqr(i_a[0]);
  o_mc.qc.lambda[0] = 0.; // lambda1
  o_mc.qc.lambda[1] = sqrd_a3 + 2*i_a[1]; // lambda2
  o_mc.qc.lambda[2] = 0.; // lambda3
  o_mc.qc.lambda[3] = util::sqr(i_a[2]) - sqrd_a3*(0.25*sqrd_a3 + i_a[1]); // lambda4
  o_mc.dlambda2_da = Eigen::Vector3d(2*i_a[0], 2., 0.);
  o_mc.dlambda4_da = Eigen::Vector3d(-i_a[0] * o_mc.qc.lambda[1], -sqrd_a3, 
    2*i_a[2]);

  o_mc.qc.epsilon_tau = util::sigpos(i_a[0]) * util::sigpos(i_a[2]);
  o_mc.qc.epsilon_k = util::sigpos(i_a[0]);

  o_mc.qc.delta = util::sqr(o_mc.qc.lambda[1]) + 4*o_mc.qc.lambda[3];
  const double sqrt_delta = sqrt(o_mc.qc.delta);
  const double sqrt_4_delta = sqrt(sqrt_delta);
  const double inv_delta = 1. / o_mc.qc.delta;
  const double inv_sqrt_delta = 1. / sqrt_delta;
  const double inv_sqrt_4_delta = 1. / sqrt_4_delta;
  o_mc.ddelta_da = Eigen::Vector3d(0., 8.*i_a[1], 8.*i_a[2]);

  if (o_mc.qc.lambda[3] >= 0.)
  {
    if (abs(i_a[0]) < kEpsilonNullTorque)
      return false;

    // a5 = 0 => unhandled singularity for the deta_da here (eta = 0 => eta^-1 goes inf)
    if (abs(i_a[2]) < kEpsilonNullForce)
      return false;

    // case I : lambda4 > 0 (also embeds case III where lambda4 == 0)
    if (!(o_mc.qc.lambda[3] == 0. && o_mc.qc.lambda[1] < 0.))
    {
      // compute alphas (solutions of the cubic) and its derivatives
      o_mc.qc.alpha[0] = -(o_mc.qc.lambda[1] - sqrt_delta);
      o_mc.qc.alpha[1] = 0.;
      o_mc.qc.alpha[2] = o_mc.qc.lambda[1] + sqrt_delta;
      o_mc.dalpha_da.col(0) = -o_mc.dlambda2_da + 0.5 * inv_sqrt_delta * o_mc.ddelta_da;
      o_mc.dalpha_da.col(1) = Eigen::Vector3d::Zero();
      o_mc.dalpha_da.col(2) = o_mc.dlambda2_da + 0.5 * inv_sqrt_delta * o_mc.ddelta_da;

      // compute m, n and r (elliptic parameters) and its derivatives
      o_mc.qc.m = 0.5 + (o_mc.qc.lambda[1] * 0.5 * inv_sqrt_delta);
      o_mc.qc.k = sqrt(o_mc.qc.m);
      o_mc.dm_da = 0.5 * inv_sqrt_delta * (o_mc.dlambda2_da - o_mc.qc.lambda[1] * 0.5 * inv_delta * 
        o_mc.ddelta_da);
      o_mc.qc.n = 1.;
      o_mc.dn_da = Eigen::Vector3d::Zero();
      o_mc.qc.r = sqrt_4_delta * inv_sqrt_2;
      o_mc.dr_da = 0.25 * inv_sqrt_2 * inv_sqrt_4_delta * inv_sqrt_delta * o_mc.ddelta_da;

      // pre-compute some usefull constants
      const double m1 = 1. - o_mc.qc.m;
      const double inv_r = 1. / o_mc.qc.r;
      const double inv_alpha3 = 1. / o_mc.qc.alpha[2];

      // compute the phase tau
      const double eta_sqrd = util::clamp(1. - sqrd_a3 * inv_alpha3, 0., 1.);
      o_mc.qc.eta = sqrt(eta_sqrd);
      const double inv_eta = 1. / o_mc.qc.eta;
      // note that arcsn_eta == F_arcsin_eta
      const double arcsin_eta = asin(o_mc.qc.eta);
      const double arcsn_eta = boost::math::ellint_1(o_mc.qc.k, arcsin_eta);
      const double E_arcsin_eta = boost::math::ellint_2(o_mc.qc.k, arcsin_eta);
      o_mc.qc.tau = o_mc.qc.epsilon_tau * arcsn_eta * inv_r;

      o_mc.deta_da = 0.5 * inv_eta * (sqrd_a3 * util::sqr(inv_alpha3)) * (o_mc.dalpha_da.col(2) - 
        Eigen::Vector3d(2*o_mc.qc.alpha[2] / i_a[0], 0., 0.));
      const double darcsn_eta_deta = 1. / (sqrt(1. - eta_sqrd) * sqrt(1. - o_mc.qc.m*eta_sqrd));
      double cn_F_arcsin_eta, dn_F_arcsin_eta, dummy;
      const double sn_F_arcsin_eta = boost::math::jacobi_elliptic(o_mc.qc.k, arcsn_eta, 
        &cn_F_arcsin_eta, &dn_F_arcsin_eta, &dummy); 
      const double cd_F_arcsin_eta = cn_F_arcsin_eta / dn_F_arcsin_eta;
      const double darcsn_eta_dm = (E_arcsin_eta - m1*arcsn_eta - o_mc.qc.m*o_mc.qc.eta*cd_F_arcsin_eta) /
        (2.*m1*o_mc.qc.m);
      const Eigen::Vector3d darcsn_eta_da = darcsn_eta_deta * o_mc.deta_da + darcsn_eta_dm * o_mc.dm_da;
      o_mc.dtau_da = (o_mc.qc.epsilon_tau * inv_r) * (-inv_r * o_mc.dr_da * arcsn_eta + darcsn_eta_da);

      // compute last intermediate constants of motions
      o_mc.qc.gamma_0 = o_mc.qc.r * o_mc.qc.tau;
      o_mc.dgamma_0_da = o_mc.dr_da * o_mc.qc.tau + o_mc.dtau_da * o_mc.qc.r;
      //double am_gamma_0;
      o_mc.qc.sn_gamma_0 = boost::math::jacobi_elliptic(o_mc.qc.k, o_mc.qc.gamma_0, &o_mc.qc.cn_gamma_0, 
        &o_mc.qc.dn_gamma_0, &o_mc.qc.am_gamma_0); 
      const double cd_gamma_0 = o_mc.qc.cn_gamma_0 / o_mc.qc.dn_gamma_0;
      const double sc_gamma_0 = o_mc.qc.sn_gamma_0 / o_mc.qc.cn_gamma_0;
      //const double am_gamma_0 = asin(o_mc.sn_gamma_0);   // should have been kept from calculation of elliptic functions...

      //const double F_am_gamma_0 = boost::math::ellint_1(o_mc.k, am_gamma_0);
      const double F_am_gamma_0 = o_mc.qc.gamma_0;
      // TODO elliptics integrals of the 1st and second kind code can be factorized for the same parameters
      o_mc.qc.E_am_gamma_0 = boost::math::ellint_2(o_mc.qc.k, o_mc.qc.am_gamma_0);

      // ddn_gamma_0_da
      const double ddn_gamma_0_dgamma_0 = -o_mc.qc.m * o_mc.qc.sn_gamma_0 * o_mc.qc.cn_gamma_0;
      const double ddn_gamma_0_dm = (1./(2.*m1)) * (o_mc.qc.sn_gamma_0 * o_mc.qc.cn_gamma_0 * ((o_mc.qc.m - 1)*o_mc.qc.gamma_0 + 
        o_mc.qc.E_am_gamma_0 - o_mc.qc.dn_gamma_0 * sc_gamma_0));
      o_mc.ddn_gamma_0_da = ddn_gamma_0_dgamma_0 * o_mc.dgamma_0_da + ddn_gamma_0_dm * o_mc.dm_da;

      // dsn_gamma_0_da
      const double dsn_gamma_0_dgamma_0 = o_mc.qc.cn_gamma_0 * o_mc.qc.dn_gamma_0;
      const double dsn_gamma_0_dm = (1. / (2. * o_mc.qc.m * m1)) * (o_mc.qc.dn_gamma_0 * o_mc.qc.cn_gamma_0 * 
        (m1 * o_mc.qc.gamma_0 - o_mc.qc.E_am_gamma_0 + o_mc.qc.m * cd_gamma_0 * o_mc.qc.sn_gamma_0));
      o_mc.dsn_gamma_0_da = dsn_gamma_0_dgamma_0 * o_mc.dgamma_0_da + dsn_gamma_0_dm * o_mc.dm_da;

      // dcn_gamma_0_da (will be used for computing dint_beta2_da)
      const double dcn_gamma_0_dgamma_0 = -o_mc.qc.sn_gamma_0 * o_mc.qc.dn_gamma_0;
      const double dcn_gamma_0_dm = (1. / (2. * o_mc.qc.m * m1)) * (o_mc.qc.sn_gamma_0 * o_mc.qc.dn_gamma_0 * 
        ((o_mc.qc.m - 1) * o_mc.qc.gamma_0 + o_mc.qc.E_am_gamma_0 - o_mc.qc.m * o_mc.qc.sn_gamma_0 * cd_gamma_0));
      o_mc.dcn_gamma_0_da = dcn_gamma_0_dgamma_0 * o_mc.dgamma_0_da + dcn_gamma_0_dm * o_mc.dm_da;

      // dE_am_gamma_0_da (will be used for computing dint_beta1_da)
      const double dE_am_gamma_0_dam_gamma_0 = o_mc.qc.dn_gamma_0;
      const double dam_gamma_0_dgamma_0 = o_mc.qc.dn_gamma_0;
      const double dam_gamma_0_dm = (1. / (2. * o_mc.qc.m * (o_mc.qc.m - 1))) * (((o_mc.qc.m - 1) * 
        o_mc.qc.gamma_0 + o_mc.qc.E_am_gamma_0) * o_mc.qc.dn_gamma_0 - o_mc.qc.m * o_mc.qc.cn_gamma_0 * 
        o_mc.qc.sn_gamma_0);
      const Eigen::Vector3d dam_gamma_0_da = dam_gamma_0_dgamma_0 * o_mc.dgamma_0_da + dam_gamma_0_dm * o_mc.dm_da;
      const double dE_am_gamma_0_dm = (o_mc.qc.E_am_gamma_0 - F_am_gamma_0) / (2. * o_mc.qc.m);
      o_mc.dE_am_gamma_0_da = dE_am_gamma_0_dam_gamma_0 * dam_gamma_0_da + dE_am_gamma_0_dm * o_mc.dm_da;

      o_mc.qc.beta1_0 = 2.*util::sqr(o_mc.qc.dn_gamma_0) - 1.;
      o_mc.qc.beta2_0 = o_mc.qc.sn_gamma_0 * o_mc.qc.dn_gamma_0;
      o_mc.dbeta1_0_da = 4. * o_mc.qc.dn_gamma_0 * o_mc.ddn_gamma_0_da;
      o_mc.dbeta2_0_da = o_mc.dsn_gamma_0_da * o_mc.qc.dn_gamma_0 + o_mc.ddn_gamma_0_da * o_mc.qc.sn_gamma_0;

    }
    //else if (o_mc.lambda[3] == 0. && o_mc.lambda[1] == 0){
    //  // lambda4 == 0 and lambda2 == 0 => all a_i are nulls
    //  // TODO
    //  return false;
    //}
    //else{
    //  // unhandled case corresponding to a3 == a5 == 0 (abnormal case)
    //  // TODO
    //  return false;
    //}
  }else if (o_mc.qc.lambda[3] < 0.){
    // a5 = 0 => unhandled singularity for the deta_da here (eta = 0 => eta^-1 goes inf)
    if (abs(i_a[2]) < kEpsilonNullForce)
      return false;

    // case II : lambda4 < 0
    if (abs(i_a[1]) > kEpsilonNullForce || abs(i_a[2]) > kEpsilonNullForce)
    {
      // case II.1 : a4 != 0  and a5 != 0 (i.e. m != 0)
      //assert (abs(i_a[0]) > kEpsilonNullTorque);
      //assert (abs(i_a[1]) > kEpsilonNullForce || abs(i_a[2]) > kEpsilonNullForce);

      // compute alphas (solutions of the cubic) and its derivatives
      o_mc.qc.alpha[0] = 0.;
      o_mc.qc.alpha[1] = o_mc.qc.lambda[1] - sqrt_delta;
      o_mc.qc.alpha[2] = o_mc.qc.lambda[1] + sqrt_delta;
      o_mc.dalpha_da.col(0) = Eigen::Vector3d::Zero();
      o_mc.dalpha_da.col(1) = o_mc.dlambda2_da - 0.5 * inv_sqrt_delta * o_mc.ddelta_da;
      o_mc.dalpha_da.col(2) = o_mc.dlambda2_da + 0.5 * inv_sqrt_delta * o_mc.ddelta_da;

      // pre-compute some usefull constants
      const double sqrt_alpha3 = sqrt(o_mc.qc.alpha[2]);
      const double inv_alpha3 = 1. / o_mc.qc.alpha[2];
      const double inv_sqrd_alpha3 = util::sqr(inv_alpha3);

      // compute m, n and r (elliptic parameters) and its derivatives
      o_mc.qc.m = (o_mc.qc.alpha[2] - o_mc.qc.alpha[1]) * inv_alpha3;
      o_mc.qc.k = sqrt(o_mc.qc.m);
      o_mc.dm_da = inv_sqrd_alpha3 * ((o_mc.qc.alpha[2] * inv_sqrt_delta - 1.) * o_mc.ddelta_da
        - 2. * sqrt_delta * o_mc.dlambda2_da);
      o_mc.qc.n = o_mc.qc.m;
      o_mc.dn_da = o_mc.dm_da;
      o_mc.qc.r = 0.5 * sqrt_alpha3;
      o_mc.dr_da = 0.25 * (1. / sqrt_alpha3) * o_mc.dalpha_da.col(2);

      // pre-compute some usefull constants
      const double inv_n = 1. / o_mc.qc.n;
      const double inv_r = 1. / o_mc.qc.r;
      const double m1 = 1. - o_mc.qc.m;

       // compute the phase tau
      const double eta_sqrd = util::clamp(inv_n * (1. - sqrd_a3 * inv_alpha3), 0., 1.);
      o_mc.qc.eta = sqrt(eta_sqrd);
      const double inv_eta = 1. / o_mc.qc.eta;
      // note that arcsn_eta == F_arcsin_eta
      const double arcsin_eta = asin(o_mc.qc.eta);
      const double arcsn_eta = boost::math::ellint_1(o_mc.qc.k, arcsin_eta);
      o_mc.qc.tau = o_mc.qc.epsilon_tau * arcsn_eta * inv_r;

      o_mc.deta_da = 0.5 * inv_eta * inv_n * ( sqrd_a3 * inv_sqrd_alpha3 * 
            (o_mc.dalpha_da.col(2) - Eigen::Vector3d(2. * o_mc.qc.alpha[2] / i_a[0], 0., 0.)) - 
            inv_n * (1. - sqrd_a3 * inv_alpha3) * o_mc.dn_da);
      const double darcsn_eta_deta = 1. / (sqrt(1. - eta_sqrd) * sqrt(1. - o_mc.qc.m * eta_sqrd));
      // TODO elliptics integrals of the 1st and second kind code can be factorized for the same parameters
      const double E_arcsin_eta = boost::math::ellint_2(o_mc.qc.k, arcsin_eta);
      double cn_F_arcsin_eta, dn_F_arcsin_eta, dummy;
      const double sn_F_arcsin_eta = boost::math::jacobi_elliptic(o_mc.qc.k, arcsn_eta, 
        &cn_F_arcsin_eta, &dn_F_arcsin_eta, &dummy); 
      const double cd_F_arcsin_eta = cn_F_arcsin_eta / dn_F_arcsin_eta;
      const double darcsn_eta_dm = (E_arcsin_eta - m1*arcsn_eta - o_mc.qc.m * o_mc.qc.eta * cd_F_arcsin_eta) 
        / (2. * m1 * o_mc.qc.m);
      const Eigen::Vector3d darcsn_eta_da = darcsn_eta_deta * o_mc.deta_da + darcsn_eta_dm * o_mc.dm_da;
      o_mc.dtau_da = o_mc.qc.epsilon_tau * inv_r * (-inv_r * o_mc.dr_da * arcsn_eta + darcsn_eta_da);
    
      // compute last intermediate constants of motions
      o_mc.qc.gamma_0 = o_mc.qc.r * o_mc.qc.tau;
      o_mc.dgamma_0_da = o_mc.dr_da * o_mc.qc.tau + o_mc.dtau_da * o_mc.qc.r;
      //double am_gamma_0;
      o_mc.qc.sn_gamma_0 = boost::math::jacobi_elliptic(o_mc.qc.k, o_mc.qc.gamma_0, &o_mc.qc.cn_gamma_0, 
        &o_mc.qc.dn_gamma_0, &o_mc.qc.am_gamma_0); 
      const double cd_gamma_0 = o_mc.qc.cn_gamma_0 / o_mc.qc.dn_gamma_0;
      const double sc_gamma_0 = o_mc.qc.sn_gamma_0 / o_mc.qc.cn_gamma_0;
      //const double am_gamma_0 = asin(o_mc.sn_gamma_0);   // should have been kept from calculation of elliptic functions...

      //const double F_am_gamma_0 = boost::math::ellint_1(o_mc.k, am_gamma_0);
      const double F_am_gamma_0 = o_mc.qc.gamma_0;
      // TODO elliptics integrals of the 1st and second kind code can be factorized for the same parameters
      o_mc.qc.E_am_gamma_0 = boost::math::ellint_2(o_mc.qc.k, o_mc.qc.am_gamma_0);

      // ddn_gamma_0_da
      const double ddn_gamma_0_dgamma_0 = -o_mc.qc.m * o_mc.qc.sn_gamma_0 * o_mc.qc.cn_gamma_0;
      const double ddn_gamma_0_dm = (1. / (2. * m1)) * (o_mc.qc.sn_gamma_0 * o_mc.qc.cn_gamma_0 * 
        ((o_mc.qc.m - 1.) * o_mc.qc.gamma_0 + o_mc.qc.E_am_gamma_0 - o_mc.qc.dn_gamma_0 * sc_gamma_0));
      o_mc.ddn_gamma_0_da = ddn_gamma_0_dgamma_0 * o_mc.dgamma_0_da + ddn_gamma_0_dm * o_mc.dm_da;

      // dcn_gamma_0_da
      const double dcn_gamma_0_dgamma_0 = -o_mc.qc.sn_gamma_0 * o_mc.qc.dn_gamma_0;
      const double dcn_gamma_0_dm = (1. / (2. * o_mc.qc.m * m1)) * (o_mc.qc.sn_gamma_0 * o_mc.qc.dn_gamma_0 * 
        ((o_mc.qc.m - 1.) * o_mc.qc.gamma_0 + o_mc.qc.E_am_gamma_0 - o_mc.qc.m * o_mc.qc.sn_gamma_0 * cd_gamma_0));
      o_mc.dcn_gamma_0_da = dcn_gamma_0_dgamma_0 * o_mc.dgamma_0_da + dcn_gamma_0_dm * o_mc.dm_da;

      // dsn_gamma_0_da
      const double dsn_gamma_0_dgamma_0 = o_mc.qc.cn_gamma_0 * o_mc.qc.dn_gamma_0;
      const double dsn_gamma_0_dm = (1. / (2. * o_mc.qc.m * m1)) * (o_mc.qc.dn_gamma_0 * o_mc.qc.cn_gamma_0 * 
        (m1 * o_mc.qc.gamma_0 - o_mc.qc.E_am_gamma_0 + o_mc.qc.m * cd_gamma_0 * o_mc.qc.sn_gamma_0));
      o_mc.dsn_gamma_0_da = dsn_gamma_0_dgamma_0 * o_mc.dgamma_0_da + dsn_gamma_0_dm * o_mc.dm_da;

      // dE_am_gamma_0_da (will be used for computing dint_beta1_da)
      const double dE_am_gamma_0_dam_gamma_0 = o_mc.qc.dn_gamma_0;
      const double dam_gamma_0_dgamma_0 = o_mc.qc.dn_gamma_0;
      const double dam_gamma_0_dm = (1. / (2. * o_mc.qc.m * (o_mc.qc.m - 1.))) * (((o_mc.qc.m - 1.) * 
        o_mc.qc.gamma_0 + o_mc.qc.E_am_gamma_0) * o_mc.qc.dn_gamma_0 - o_mc.qc.m * o_mc.qc.cn_gamma_0 * 
        o_mc.qc.sn_gamma_0);
      const Eigen::Vector3d dam_gamma_0_da = dam_gamma_0_dgamma_0 * o_mc.dgamma_0_da + dam_gamma_0_dm * o_mc.dm_da;
      const double dE_am_gamma_0_dm = (o_mc.qc.E_am_gamma_0 - F_am_gamma_0) / (2. * o_mc.qc.m);
      o_mc.dE_am_gamma_0_da = dE_am_gamma_0_dam_gamma_0 * dam_gamma_0_da + dE_am_gamma_0_dm * o_mc.dm_da;

      o_mc.qc.beta1_0 = 1. - 2.*util::sqr(o_mc.qc.sn_gamma_0);
      o_mc.qc.beta2_0 = o_mc.qc.sn_gamma_0 * o_mc.qc.cn_gamma_0;
      o_mc.dbeta1_0_da = -4. * o_mc.qc.sn_gamma_0 * o_mc.dsn_gamma_0_da;
      o_mc.dbeta2_0_da = o_mc.dsn_gamma_0_da * o_mc.qc.cn_gamma_0 + o_mc.dcn_gamma_0_da * o_mc.qc.sn_gamma_0;
    }
    //else{
    //  // case II.2 : a4 == 0  and a5 == 0 (i.e. m == 0)
    //  // TODO
    //  return false;
    //}
  }
  return true;
}

bool computeDqDaAtPositionT(double i_t, const MotionConstantsDqDa& i_mc, Eigen::Matrix3d& o_dqda)
{
  static const double kEpsilonLambda = 1.e-12;

  double k_t = 0.;      // k(t) is rod curvature at position t
  double k_dot_t = 0.;
  if (i_mc.qc.lambda[3] >= 0.)
  {
    //if (abs(i_a[0]) < kEpsilonNullTorque)
    //  return false;

    //// a5 = 0 => unhandled singularity for the deta_da here (eta = 0 => eta^-1 goes inf)
    //if (abs(i_a[2]) < kEpsilonNullForce)
    //  return false;

    // case I : lambda4 > 0 (also embeds case III where lambda4 == 0)
    if (!(i_mc.qc.lambda[3] == 0. && i_mc.qc.lambda[1] < 0.))
    {
      // pre-computed vars
      const double m1 = 1. - i_mc.qc.m;
      const double inv_r = 1. / i_mc.qc.r;

      // gamma(t) and elliptic functions valued
      const double gamma_t = i_mc.qc.r * (i_t + i_mc.qc.tau);
      const Eigen::Vector3d dgamma_t_da = i_mc.dr_da * (i_t + i_mc.qc.tau) + i_mc.qc.r * i_mc.dtau_da;
      double cn_gamma_t, dn_gamma_t, am_gamma_t;
      const double sn_gamma_t = boost::math::jacobi_elliptic(i_mc.qc.k, gamma_t, &cn_gamma_t, 
        &dn_gamma_t, &am_gamma_t); 
      const double cd_gamma_t = cn_gamma_t / dn_gamma_t;
      const double sc_gamma_t = sn_gamma_t / cn_gamma_t;
      //const double am_gamma_t = asin(sn_gamma_t);   // should have been kept from calculation of elliptic functions...

      //const double F_am_gamma_t = boost::math::ellint_1(i_mc.k, am_gamma_t);
      const double F_am_gamma_t = gamma_t;
      const double E_am_gamma_t = boost::math::ellint_2(i_mc.qc.k, am_gamma_t);

      // derivatives
      // ddn_gamma_t_da
      const double ddn_gamma_t_dgamma_t = -i_mc.qc.m * sn_gamma_t * cn_gamma_t;
      const double ddn_gamma_t_dm = (1. / (2. * m1)) * (sn_gamma_t * cn_gamma_t * ((i_mc.qc.m - 1.) * 
        gamma_t + E_am_gamma_t - dn_gamma_t * sc_gamma_t));
      const Eigen::Vector3d ddn_gamma_t_da = ddn_gamma_t_dgamma_t * dgamma_t_da + ddn_gamma_t_dm * i_mc.dm_da;

      // dcn_gamma_t_da
      const double dcn_gamma_t_dgamma_t = -sn_gamma_t * dn_gamma_t;
      const double dcn_gamma_t_dm = (1. / (2. * i_mc.qc.m * m1)) * (sn_gamma_t * dn_gamma_t * 
        ((i_mc.qc.m - 1.) * gamma_t + E_am_gamma_t - i_mc.qc.m * sn_gamma_t * cd_gamma_t));
      const Eigen::Vector3d dcn_gamma_t_da = dcn_gamma_t_dgamma_t * dgamma_t_da + dcn_gamma_t_dm * i_mc.dm_da;

      // dE_am_gamma_t_da
      const double dE_am_gamma_t_dam_gamma_t = dn_gamma_t;
      const double dam_gamma_t_dgamma_t = dn_gamma_t;
      const double dam_gamma_t_dm = (1. / (2. * i_mc.qc.m * (i_mc.qc.m - 1.))) * (((i_mc.qc.m - 1.) *
        gamma_t + E_am_gamma_t) * dn_gamma_t - i_mc.qc.m * cn_gamma_t * sn_gamma_t);
      const Eigen::Vector3d dam_gamma_t_da = dam_gamma_t_dgamma_t * dgamma_t_da + dam_gamma_t_dm * i_mc.dm_da;
      const double dE_am_gamma_t_dm = (E_am_gamma_t - F_am_gamma_t) / (2. * i_mc.qc.m);
      const Eigen::Vector3d dE_am_gamma_t_da = dE_am_gamma_t_dam_gamma_t * dam_gamma_t_da + dE_am_gamma_t_dm * i_mc.dm_da;

      const double int_beta1 = (2. * inv_r) * (E_am_gamma_t - i_mc.qc.E_am_gamma_0) - i_t;
      const double int_beta2 = -inv_r * (cn_gamma_t - i_mc.qc.cn_gamma_0);
      const Eigen::Vector3d dint_beta1_da = inv_r * (-(1. / (2. * i_mc.qc.delta)) * i_mc.ddelta_da * 
        (E_am_gamma_t - i_mc.qc.E_am_gamma_0) + 2. * (dE_am_gamma_t_da - i_mc.dE_am_gamma_0_da));
      const Eigen::Vector3d dint_beta2_da = inv_r * (inv_r * i_mc.dr_da * (cn_gamma_t - i_mc.qc.cn_gamma_0) - 
            dcn_gamma_t_da + i_mc.dcn_gamma_0_da);

      // dq1 / da (q1 is angle theta)
      o_dqda.row(0) = (2. * i_mc.qc.epsilon_k * (( 1. / (i_mc.qc.k * i_mc.qc.sn_gamma_0) ) * i_mc.ddn_gamma_0_da - 
        ( 1. / (i_mc.qc.k * sn_gamma_t) ) * ddn_gamma_t_da)).transpose();

      // dq2 / da (q2 is x position)
      o_dqda.row(1) = (i_mc.dbeta1_0_da * int_beta1 + dint_beta1_da * i_mc.qc.beta1_0 + 
        4 * (int_beta2 * (i_mc.dm_da * i_mc.qc.beta2_0 + i_mc.dbeta2_0_da * i_mc.qc.m ) + dint_beta2_da * 
        i_mc.qc.m * i_mc.qc.beta2_0)).transpose();

      // dq3 / da (q3 is y position)
      o_dqda.row(2) = (i_mc.qc.epsilon_k * ((1. / i_mc.qc.k) * i_mc.dm_da * (i_mc.qc.beta1_0 * int_beta2 - 
        i_mc.qc.beta2_0 * int_beta1) + 2 * i_mc.qc.k * (i_mc.dbeta1_0_da * int_beta2 + dint_beta2_da * 
        i_mc.qc.beta1_0 - i_mc.dbeta2_0_da * int_beta1 - dint_beta1_da * i_mc.qc.beta2_0))).transpose();
    }
    //else if (i_mc.lambda[3] == 0. && i_mc.lambda[1] == 0){
    //  // lambda4 == 0 and lambda2 == 0 => all a_i are nulls
    //  // TODO
    //  return false;
    //}
  }else if (i_mc.qc.lambda[3] < 0.){
    // a5 = 0 => unhandled singularity for the deta_da here (eta = 0 => eta^-1 goes inf)
    //if (abs(i_a[2]) < kEpsilonNullForce)
    //  return false;

    // case II : lambda4 < 0
    if (i_mc.qc.k != 0.)
    {
      // pre-computed vars
      const double m1 = 1. - i_mc.qc.m;
      const double inv_m = 1. / i_mc.qc.m;
      const double inv_r = 1. / i_mc.qc.r;

      // gamma(t) and elliptic functions valued
      const double gamma_t = i_mc.qc.r * (i_t + i_mc.qc.tau);
      const Eigen::Vector3d dgamma_t_da = i_mc.dr_da * (i_t + i_mc.qc.tau) + i_mc.qc.r * i_mc.dtau_da;
      double cn_gamma_t, dn_gamma_t, am_gamma_t;
      const double sn_gamma_t = boost::math::jacobi_elliptic(i_mc.qc.k, gamma_t, &cn_gamma_t, &dn_gamma_t, &am_gamma_t); 
      const double cd_gamma_t = cn_gamma_t / dn_gamma_t;
      const double sc_gamma_t = sn_gamma_t / cn_gamma_t;
      //const double am_gamma_t = asin(sn_gamma_t);   // should have been kept from calculation of elliptic functions...

      //const double F_am_gamma_t = boost::math::ellint_1(i_mc.k, am_gamma_t);
      const double F_am_gamma_t = gamma_t;
      const double E_am_gamma_t = boost::math::ellint_2(i_mc.qc.k, am_gamma_t);

      const double int_beta1_p = i_t * (i_mc.qc.m - 2.) + 2. * inv_r * (E_am_gamma_t - i_mc.qc.E_am_gamma_0);
      const double int_beta1 = inv_m * int_beta1_p;
      const double int_beta2 = -inv_m * inv_r * (dn_gamma_t - i_mc.qc.dn_gamma_0);

      // derivatives
      // ddn_gamma_t_da
      const double ddn_gamma_t_dgamma_t = -i_mc.qc.m * sn_gamma_t * cn_gamma_t;
      const double ddn_gamma_t_dm = (1. / (2. * m1)) * (sn_gamma_t * cn_gamma_t * ((i_mc.qc.m - 1.) * 
        gamma_t + E_am_gamma_t - dn_gamma_t * sc_gamma_t));
      const Eigen::Vector3d ddn_gamma_t_da = ddn_gamma_t_dgamma_t * dgamma_t_da + ddn_gamma_t_dm * i_mc.dm_da;

      // dsn_gamma_t_da
      const double dsn_gamma_t_dgamma_t = cn_gamma_t * dn_gamma_t;
      const double dsn_gamma_t_dm = (1. / (2. * i_mc.qc.m * m1)) * (dn_gamma_t * cn_gamma_t * (m1 *gamma_t - 
        E_am_gamma_t + i_mc.qc.m * cd_gamma_t * sn_gamma_t));
      const Eigen::Vector3d dsn_gamma_t_da = dsn_gamma_t_dgamma_t * dgamma_t_da + dsn_gamma_t_dm * i_mc.dm_da;
        
      // dE_am_gamma_t_da
      const double dE_am_gamma_t_dam_gamma_t = dn_gamma_t;
      const double dam_gamma_t_dgamma_t = dn_gamma_t;
      const double dam_gamma_t_dm = (1. / (2. * i_mc.qc.m * (i_mc.qc.m - 1))) * (((i_mc.qc.m - 1.) * 
        gamma_t + E_am_gamma_t) * dn_gamma_t - i_mc.qc.m * cn_gamma_t * sn_gamma_t);
      const Eigen::Vector3d dam_gamma_t_da = dam_gamma_t_dgamma_t * dgamma_t_da + dam_gamma_t_dm * i_mc.dm_da;
      const double dE_am_gamma_t_dm = (E_am_gamma_t - F_am_gamma_t) / (2. * i_mc.qc.m);
      const Eigen::Vector3d dE_am_gamma_t_da = dE_am_gamma_t_dam_gamma_t * dam_gamma_t_da + 
        dE_am_gamma_t_dm * i_mc.dm_da;

      // dint_beta1_da
      const Eigen::Vector3d dint_beta1_p_da = i_t * i_mc.dm_da + 2. * inv_r * (-inv_r * i_mc.dr_da * 
        (E_am_gamma_t - i_mc.qc.E_am_gamma_0) + dE_am_gamma_t_da - i_mc.dE_am_gamma_0_da);
      const Eigen::Vector3d dint_beta1_da = inv_m * (-inv_m * i_mc.dm_da * int_beta1_p + dint_beta1_p_da);

      // dint_beta2_da
      const Eigen::Vector3d dint_beta2_da = inv_m * inv_r * ((dn_gamma_t - i_mc.qc.dn_gamma_0) *
        (inv_m * i_mc.dm_da + inv_r * i_mc.dr_da) - (ddn_gamma_t_da - i_mc.ddn_gamma_0_da));

      // dq1 / da (q1 is angle theta)
      o_dqda.row(0) = (2 * i_mc.qc.epsilon_k * ( (1. / cn_gamma_t) * dsn_gamma_t_da - 
        (1. / i_mc.qc.cn_gamma_0) * i_mc.dsn_gamma_0_da)).transpose();

      // dq2 / da (q2 is x position)
      o_dqda.row(1) = (i_mc.dbeta1_0_da * int_beta1 + i_mc.qc.beta1_0 * dint_beta1_da + 
        4. * (i_mc.dbeta2_0_da * int_beta2 + i_mc.qc.beta2_0 * dint_beta2_da)).transpose();

      // dq3 / da (q3 is y position)
      o_dqda.row(2) = (2 * i_mc.qc.epsilon_k * (i_mc.dbeta1_0_da * int_beta2 + i_mc.qc.beta1_0 * 
        dint_beta2_da - i_mc.dbeta2_0_da * int_beta1 - i_mc.qc.beta2_0 * dint_beta1_da)).transpose();
    }
    //else{
    //  // case II.2 : a4 == 0  and a5 == 0 (i.e. m == 0)

    //  // TODO
    //}
  }
  return true;
}

}	// namespace rod2d
}	// namespace qserl
