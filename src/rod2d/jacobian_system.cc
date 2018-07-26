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

#include "jacobian_system.h"

#include <boost/bind.hpp>

namespace qserl {
namespace rod2d {

const double JacobianSystem::kStabilityThreshold = 1.e-8;        /** Threshold for Jacobian determinant. */
const double JacobianSystem::kStabilityTolerance = 1.e-12;      /** Tolerance for which Jacobian determinant vanishes. */

JacobianSystem::state_type
JacobianSystem::defaultState()
{
  state_type defaultStateArray;
  defaultStateArray.assign(0.);
  return defaultStateArray;
}

JacobianSystem::JacobianSystem(double i_inv_stiffness,
                               double i_dt,
                               const std::vector<WorkspaceIntegratedState::costate_type>& i_mu,
                               Parameters::RodModelT i_rodModel) :
    m_inv_c(i_inv_stiffness),
    m_dt(i_dt),
    m_mu(i_mu),
    m_rodModel(i_rodModel)
{
  assert (m_dt > 0. && "integration step time must be positive.");

  if(m_rodModel == Parameters::RM_INEXTENSIBLE)
  {
    m_evaluationCallback = boost::bind(&JacobianSystem::evaluateInextensible, this, _1, _2, _3);
  }
  else
    assert(false && "invalid rod model");
}

JacobianSystem::~JacobianSystem()
{
}


void
JacobianSystem::operator()(const state_type& i_MJ,
                           state_type& o_dMJdt,
                           double i_t)
{
  return m_evaluationCallback(i_MJ, o_dMJdt, i_t);
}

void
JacobianSystem::evaluateInextensible(const state_type& i_MJ,
                                     state_type& o_dMJdt,
                                     double i_t)
{
  const size_t k = static_cast<size_t>(i_t / m_dt);
  assert (k >= 0 && k < m_mu.size() && "Given mu values array not consistent with current integration parameters.");

  const WorkspaceIntegratedState::costate_type& mu_k = m_mu[k];

  // F matrix
  Eigen::Matrix<double, 3, 3> F;
  F << 0., mu_k[2], mu_k[1],
      -mu_k[2], 0, -mu_k[0],
      0., -1., 0.;

  // G matrix
  Eigen::Matrix<double, 3, 3> G;
  G.setZero();
  G.diagonal() << 0., 0., m_inv_c;

  // H matrix
  Eigen::Matrix<double, 3, 3> H;
  H << 0., mu_k[2], 0.,
      -mu_k[2], 0., 1.,
      0., 0., 0.;

  // create mapping between mj array and M & J eigen matrices
  const Eigen::Map<const Eigen::Matrix<double, 3, 3> > M_e(i_MJ.data());
  const Eigen::Map<const Eigen::Matrix<double, 3, 3> > J_e(i_MJ.data() + 9);

  // create mapping between dmjdt array and dMdt & dJdt eigen matrices
  Eigen::Map<Eigen::Matrix<double, 3, 3> > dMdt_e(o_dMJdt.data());
  Eigen::Map<Eigen::Matrix<double, 3, 3> > dJdt_e(o_dMJdt.data() + 9);

  dMdt_e = F * M_e;
  dJdt_e = G * M_e + H * J_e;
}

}  // namespace rod2d
}  // namespace qserl
