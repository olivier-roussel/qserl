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

#include "costate_system.h"
#include <Eigen/Geometry>

namespace qserl {
namespace rod3d {

CostateSystem::state_type
CostateSystem::defaultState()
{
  state_type defaultStateArray;
  defaultStateArray.fill(0.);
  return defaultStateArray;
}

CostateSystem::CostateSystem(const Eigen::Matrix<double, 6, 1>& i_inv_stiffness,
                             double i_length,
                             Parameters::RodModelT i_rodModel) :
    m_inv_c(i_inv_stiffness),
    m_length(i_length),
    m_rodModel(i_rodModel)
{
  using namespace std::placeholders;

  if(m_rodModel == Parameters::RM_INEXTENSIBLE)
  {
    m_evaluationCallback = std::bind(&CostateSystem::evaluateInextensible, this, _1, _2, _3);
  }
  else if(m_rodModel == Parameters::RM_EXTENSIBLE_SHEARABLE)
  {
    m_evaluationCallback = std::bind(&CostateSystem::evaluateExtensibleShearable, this, _1, _2, _3);
  }
  else
    assert(false && "invalid rod model");
}

CostateSystem::~CostateSystem()
{
}

void
CostateSystem::operator()(const state_type& i_mu,
                          state_type& o_dmudt,
                          double i_t)
{
  return m_evaluationCallback(i_mu, o_dmudt, i_t);
}

void
CostateSystem::evaluateInextensible(const state_type& i_mu,
                                    state_type& o_dmudt,
                                    double /*i_t*/)
{
  static const Eigen::Matrix<double, 3, 1> ke_1 = Eigen::Matrix3d::Identity().col(0);
  const Eigen::Map<const Eigen::Matrix<double, 3, 1> > m_e(i_mu.data());
  const Eigen::Map<const Eigen::Matrix<double, 3, 1> > f_e(i_mu.data() + 3);
  const Eigen::Matrix<double, 3, 1> u = m_e.cwiseProduct(m_inv_c.block<3, 1>(0, 0));

  Eigen::Map<Eigen::Matrix<double, 3, 1> > dmdt_e(o_dmudt.data());
  Eigen::Map<Eigen::Matrix<double, 3, 1> > dfdt_e(o_dmudt.data() + 3);

  dmdt_e = -u.cross(m_e) - (ke_1).cross(f_e);
  dfdt_e = -u.cross(f_e);
}

void
CostateSystem::evaluateExtensibleShearable(const state_type& i_mu,
                                           state_type& o_dmudt,
                                           double /*i_t*/)
{
  static const Eigen::Matrix<double, 3, 1> ke_1 = Eigen::Matrix3d::Identity().col(0);
  const Eigen::Map<const Eigen::Matrix<double, 3, 1> > m_e(i_mu.data());
  const Eigen::Map<const Eigen::Matrix<double, 3, 1> > f_e(i_mu.data() + 3);
  const Eigen::Matrix<double, 3, 1> u_m = m_e.cwiseProduct(m_inv_c.block<3, 1>(0, 0));
  const Eigen::Matrix<double, 3, 1> u_f = f_e.cwiseProduct(m_inv_c.block<3, 1>(3, 0));

  Eigen::Map<Eigen::Matrix<double, 3, 1> > dmdt_e(o_dmudt.data());
  Eigen::Map<Eigen::Matrix<double, 3, 1> > dfdt_e(o_dmudt.data() + 3);

  dmdt_e = -u_m.cross(m_e) - (u_f + ke_1).cross(f_e);
  dfdt_e = -u_m.cross(f_e);
}

}  // namespace rod3d
}  // namespace qserl
