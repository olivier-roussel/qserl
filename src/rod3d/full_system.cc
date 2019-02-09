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

#include "full_system.h"
#include <Eigen/Geometry>
#include "util/lie_algebra_utils.h"

namespace qserl {
namespace rod3d {

FullSystem::state_type
FullSystem::defaultState()
{
  state_type defaultStateArray;
  defaultStateArray.fill(0.);
  return defaultStateArray;
}

FullSystem::FullSystem(const Parameters& i_params,
                             double i_dt) :
  m_inv_c(i_params.stiffnessCoefficients.cwiseInverse()),
  m_rodParameters(i_params),
  m_dt(i_dt),
  m_stability_threshold(1.e-5),
  m_stability_tolerance(1.e-12)
{
  using namespace std::placeholders;

  m_b[0] = m_inv_c[2] - m_inv_c[1];
  m_b[1] = m_inv_c[0] - m_inv_c[2];
  m_b[2] = m_inv_c[1] - m_inv_c[0];
  m_b[3] = m_inv_c[5] - m_inv_c[4];
  m_b[4] = m_inv_c[3] - m_inv_c[5];
  m_b[5] = m_inv_c[4] - m_inv_c[3];

  if(m_rodParameters.rodModel == Parameters::RM_INEXTENSIBLE)
  {
    m_evaluationCallback = std::bind(&FullSystem::evaluateInextensible, this, _1, _2, _3);
  }
  else if(m_rodParameters.rodModel == Parameters::RM_EXTENSIBLE_SHEARABLE)
  {
    m_evaluationCallback = std::bind(&FullSystem::evaluateExtensibleShearable, this, _1, _2, _3);
  }
  else if(m_rodParameters.rodModel == Parameters::RM_INEXTENSIBLE_WITH_GRAVITY)
  {
    // XXX Note that w will be pointing to the opposite direction of gravity
    const Eigen::Vector3d w = -m_rodParameters.gravity * m_rodParameters.unitaryMass;
    m_w_x_0 = Eigen::Vector4d{w[0], w[1], w[2], 0.};
    m_evaluationCallback = std::bind(&FullSystem::evaluateInextensibleWithGravity, this, _1, _2, _3);
  }
  else
    assert(false && "invalid rod model");
}

FullSystem::~FullSystem()
{
}

void
FullSystem::operator()(const state_type& i_x,
                          state_type& o_dxdt,
                          double i_t)
{
  return m_evaluationCallback(i_x, o_dxdt, i_t);
}

void
FullSystem::evaluateInextensible(const state_type& i_x,
                                    state_type& o_dxdt,
                                    double /*i_t*/)
{
  // ----------------------
  // costate
  static const Eigen::Matrix<double, 3, 1> ke_1 = Eigen::Matrix3d::Identity().col(0);
  const Eigen::Map<const Eigen::Matrix<double, 3, 1> > m_e(i_x.data() + mu_index());
  const Eigen::Map<const Eigen::Matrix<double, 3, 1> > f_e(i_x.data() + mu_index() + 3);
  const Eigen::Matrix<double, 3, 1> u = m_e.cwiseProduct(m_inv_c.block<3, 1>(0, 0));

  Eigen::Map<Eigen::Matrix<double, 3, 1> > dmdt_e(o_dxdt.data() + mu_index());
  Eigen::Map<Eigen::Matrix<double, 3, 1> > dfdt_e(o_dxdt.data() + mu_index() + 3);

  dmdt_e = -u.cross(m_e) - (ke_1).cross(f_e);
  dfdt_e = -u.cross(f_e);

  // ----------------------
  // state
  Eigen::Matrix4d u_hat_h;
  u_hat_h << 0, -u[2], u[1], 1.,
    u[2], 0., -u[0], 0.,
    -u[1], u[0], 0., 0.,
    0., 0., 0., 0.;

  const Eigen::Map<const Eigen::Matrix4d> q_e(i_x.data() + q_index());
  Eigen::Map<Eigen::Matrix4d> dqdt_e(o_dxdt.data() + q_index());

  dqdt_e = q_e * u_hat_h;

  // ----------------------
  // Jacobians
  const Eigen::Map<const Eigen::Matrix<double, 6, 1> > mu_k(i_x.data() + mu_index());

  // F matrix
  Eigen::Matrix<double, 6, 6> F;
  F << 0., mu_k[2] * m_b[0], mu_k[1] * m_b[0], 0., 0., 0.,
    mu_k[2] * m_b[1], 0, mu_k[0] * m_b[1], 0., 0., 1.,
    mu_k[1] * m_b[2], mu_k[0] * m_b[2], 0., 0., -1., 0.,
    0., -mu_k[5] * m_inv_c[1], mu_k[4] * m_inv_c[2], 0., u[2], -u[1],
    mu_k[5] * m_inv_c[0], 0, -mu_k[3] * m_inv_c[2], -u[2], 0., u[0],
    -mu_k[4] * m_inv_c[0], mu_k[3] * m_inv_c[1], 0, u[1], -u[0], 0.;

  // G matrix
  Eigen::Matrix<double, 6, 6> G;
  G.setZero();
  G.diagonal() << m_inv_c[0], m_inv_c[1], m_inv_c[2], 0, 0, 0;

  // H matrix
  Eigen::Matrix<double, 6, 6> H;
  H << 0, u[2], -u[1], 0, 0, 0,
    -u[2], 0, u[0], 0, 0, 0,
    u[1], -u[0], 0, 0, 0, 0,
    0, 0, 0, 0, u[2], -u[1],
    0, 0, 1, -u[2], 0, u[0],
    0, -1, 0, u[1], -u[0], 0;

  // create mapping between mj array and M & J eigen matrices
  const Eigen::Map<const Eigen::Matrix<double, 6, 6> > M_e(i_x.data() + MJ_index());
  const Eigen::Map<const Eigen::Matrix<double, 6, 6> > J_e(i_x.data() + MJ_index() + 36);

  // create mapping between dmjdt array and dMdt & dJdt eigen matrices
  Eigen::Map<Eigen::Matrix<double, 6, 6> > dMdt_e(o_dxdt.data() + MJ_index());
  Eigen::Map<Eigen::Matrix<double, 6, 6> > dJdt_e(o_dxdt.data() + MJ_index() + 36);

  dMdt_e = F * M_e;
  dJdt_e = G * M_e + H * J_e;
}

void
FullSystem::evaluateInextensibleWithGravity(const state_type& i_x,
                                               state_type& o_dxdt,
                                               double /*i_t*/)
{
  // ----------------------
  // costate
  static const Eigen::Matrix<double, 3, 1> ke_1 = Eigen::Matrix3d::Identity().col(0);
  const Eigen::Map<const Eigen::Matrix<double, 3, 1> > m_e(i_x.data() + mu_index());
  const Eigen::Map<const Eigen::Matrix<double, 3, 1> > f_e(i_x.data() + mu_index() + 3);
  const Eigen::Matrix<double, 3, 1> u = m_e.cwiseProduct(m_inv_c.block<3, 1>(0, 0));

  Eigen::Map<Eigen::Matrix<double, 3, 1> > dmdt_e(o_dxdt.data() + mu_index());
  Eigen::Map<Eigen::Matrix<double, 3, 1> > dfdt_e(o_dxdt.data() + mu_index() + 3);

  const Eigen::Map<const Eigen::Matrix4d> q_e(i_x.data() + q_index());
  const Eigen::Vector3d w_x = (q_e * m_w_x_0).block<3, 1>(0, 0);

  dmdt_e = -u.cross(m_e) - (ke_1).cross(f_e);
  dfdt_e = -u.cross(f_e) + w_x;

  // ----------------------
  // state
  Eigen::Matrix4d u_hat_h;
  u_hat_h << 0, -u[2], u[1], 1.,
    u[2], 0., -u[0], 0.,
    -u[1], u[0], 0., 0.,
    0., 0., 0., 0.;

  Eigen::Map<Eigen::Matrix4d> dqdt_e(o_dxdt.data() + q_index());

  dqdt_e = q_e * u_hat_h;

  // ----------------------
  // Jacobians
  const Eigen::Map<const Eigen::Matrix<double, 6, 1> > mu_k(i_x.data() + mu_index());

  // F matrix
  Eigen::Matrix<double, 6, 6> F;
  F << 0., mu_k[2] * m_b[0], mu_k[1] * m_b[0], 0., 0., 0.,
    mu_k[2] * m_b[1], 0, mu_k[0] * m_b[1], 0., 0., 1.,
    mu_k[1] * m_b[2], mu_k[0] * m_b[2], 0., 0., -1., 0.,
    0., -mu_k[5] * m_inv_c[1], mu_k[4] * m_inv_c[2], 0., u[2], -u[1],
    mu_k[5] * m_inv_c[0], 0, -mu_k[3] * m_inv_c[2], -u[2], 0., u[0],
    -mu_k[4] * m_inv_c[0], mu_k[3] * m_inv_c[1], 0, u[1], -u[0], 0.;

  // G matrix
  Eigen::Matrix<double, 6, 6> G;
  G.setZero();
  G.diagonal() << m_inv_c[0], m_inv_c[1], m_inv_c[2], 0, 0, 0;

  // H matrix
  Eigen::Matrix<double, 6, 6> H;
  H << 0, u[2], -u[1], 0, 0, 0,
    -u[2], 0, u[0], 0, 0, 0,
    u[1], -u[0], 0, 0, 0, 0,
    0, 0, 0, 0, u[2], -u[1],
    0, 0, 1, -u[2], 0, u[0],
    0, -1, 0, u[1], -u[0], 0;

  // K matrix
  Eigen::Matrix<double, 6, 6> K;
  K.setZero();
  K.block<3, 3>(3, 0) = util::hat(w_x);

  // create mapping between mj array and M & J eigen matrices
  const Eigen::Map<const Eigen::Matrix<double, 6, 6> > M_e(i_x.data() + MJ_index());
  const Eigen::Map<const Eigen::Matrix<double, 6, 6> > J_e(i_x.data() + MJ_index() + 36);

  // create mapping between dmjdt array and dMdt & dJdt eigen matrices
  Eigen::Map<Eigen::Matrix<double, 6, 6> > dMdt_e(o_dxdt.data() + MJ_index());
  Eigen::Map<Eigen::Matrix<double, 6, 6> > dJdt_e(o_dxdt.data() + MJ_index() + 36);

  dMdt_e = F * M_e - K * J_e;
  dJdt_e = G * M_e + H * J_e;
}

void
FullSystem::evaluateExtensibleShearable(const state_type& i_x,
                                           state_type& o_dxdt,
                                           double /*i_t*/)
{
  // ----------------------
  // costate
  static const Eigen::Matrix<double, 3, 1> ke_1 = Eigen::Matrix3d::Identity().col(0);
  const Eigen::Map<const Eigen::Matrix<double, 3, 1> > m_e(i_x.data() + mu_index());
  const Eigen::Map<const Eigen::Matrix<double, 3, 1> > f_e(i_x.data() + mu_index() + 3);
  const Eigen::Matrix<double, 3, 1> u_m = m_e.cwiseProduct(m_inv_c.block<3, 1>(0, 0));
  const Eigen::Matrix<double, 3, 1> u_f = f_e.cwiseProduct(m_inv_c.block<3, 1>(3, 0));

  Eigen::Map<Eigen::Matrix<double, 3, 1> > dmdt_e(o_dxdt.data() + mu_index());
  Eigen::Map<Eigen::Matrix<double, 3, 1> > dfdt_e(o_dxdt.data() + mu_index() + 3);

  dmdt_e = -u_m.cross(m_e) - (u_f + ke_1).cross(f_e);
  dfdt_e = -u_m.cross(f_e);

  // ----------------------
  // state
  Eigen::Matrix4d u_hat_h;
  u_hat_h << 0, -u_m[2], u_m[1], (1 + u_f[0]),
    u_m[2], 0., -u_m[0], u_f[1],
    -u_m[1], u_m[0], 0., u_f[2],
    0., 0., 0., 0.;

  const Eigen::Map<const Eigen::Matrix4d> q_e(i_x.data() + q_index());
  Eigen::Map<Eigen::Matrix4d> dqdt_e(o_dxdt.data() + q_index());

  dqdt_e = q_e * u_hat_h;

  // ----------------------
  // Jacobians
  const Eigen::Map<const Eigen::Matrix<double, 6, 1> > mu_k(i_x.data() + mu_index());

  // F matrix
  Eigen::Matrix<double, 6, 6> F;
  F << 0., mu_k[2] * m_b[0], mu_k[1] * m_b[0], 0., mu_k[5] * m_b[3], mu_k[4] * m_b[3],
    mu_k[2] * m_b[1], 0, mu_k[0] * m_b[1], mu_k[5] * m_b[4], 0., 1 + mu_k[3] * m_b[4],
    mu_k[1] * m_b[2], mu_k[0] * m_b[2], 0., mu_k[4] * m_b[5], -1 + mu_k[3] * m_b[5], 0.,
    0., -mu_k[5] * m_inv_c[1], mu_k[4] * m_inv_c[2], 0., u_m[2], -u_m[1],
    mu_k[5] * m_inv_c[0], 0, -mu_k[3] * m_inv_c[2], -u_m[2], 0., u_m[0],
    -mu_k[4] * m_inv_c[0], mu_k[3] * m_inv_c[1], 0, u_m[1], -u_m[0], 0.;

  // G matrix
  Eigen::Matrix<double, 6, 6> G;
  G.setZero();
  G.diagonal() << m_inv_c[0], m_inv_c[1], m_inv_c[2], m_inv_c[3], m_inv_c[4], m_inv_c[5];

  // H matrix
  Eigen::Matrix<double, 6, 6> H;
  H << 0, u_m[2], -u_m[1], 0, 0, 0,
    -u_m[2], 0, u_m[0], 0, 0, 0,
    u_m[1], -u_m[0], 0, 0, 0, 0,
    0, u_f[2], -u_f[1], 0, u_m[2], -u_m[1],
    -u_f[2], 0, 1 + u_f[0], -u_m[2], 0, u_m[0],
    u_f[1], -1 - u_f[0], 0, u_m[1], -u_m[0], 0;

  // create mapping between mj array and M & J eigen matrices
  const Eigen::Map<const Eigen::Matrix<double, 6, 6> > M_e(i_x.data() + MJ_index());
  const Eigen::Map<const Eigen::Matrix<double, 6, 6> > J_e(i_x.data() + MJ_index() + 36);

  // create mapping between dmjdt array and dMdt & dJdt eigen matrices
  Eigen::Map<Eigen::Matrix<double, 6, 6> > dMdt_e(o_dxdt.data() + MJ_index());
  Eigen::Map<Eigen::Matrix<double, 6, 6> > dJdt_e(o_dxdt.data() + MJ_index() + 36);

  dMdt_e = F * M_e;
  dJdt_e = G * M_e + H * J_e;
}

double
FullSystem::jacobianStabilityThreshold() const
{
  return m_stability_threshold;
}

void
FullSystem::jacobianStabilityThreshold(double stability_threshold)
{
  m_stability_threshold = stability_threshold;
}

double
FullSystem::jacobianStabilityTolerance() const
{
  return m_stability_tolerance;
}

void
FullSystem::jacobianStabilityTolerance(double stability_tolerance)
{
  m_stability_tolerance = stability_tolerance;
}

}  // namespace rod3d
}  // namespace qserl
