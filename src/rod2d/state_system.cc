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

#include "state_system.h"

#include <boost/bind.hpp>

namespace qserl {
namespace rod2d {

const StateSystem::state_type StateSystem::kDefaultState = { {  0. } };

StateSystem::StateSystem(double i_inv_stiffness, double i_length, 
	double i_dt, const std::vector<WorkspaceIntegratedState::costate_type>& i_mu, 
	 Parameters::RodModelT i_rodModel) :
	m_inv_c(i_inv_stiffness),
	m_dt(i_dt),
	m_mu(i_mu),
	m_length(i_length),
	m_rodModel(i_rodModel)
{
	assert (m_dt > 0. && "integration step time must be positive.");

	if (m_rodModel == Parameters::RM_INEXTENSIBLE)
		m_evaluationCallback = boost::bind(&StateSystem::evaluateInextensible, this, _1, _2, _3);
	else
		assert(false && "invalid rod model");
}

StateSystem::~StateSystem()
{
}

void StateSystem::operator() (const state_type& i_q, state_type& o_dqdt, double i_t)
{
	return m_evaluationCallback(i_q, o_dqdt, i_t);
}

void StateSystem::evaluateInextensible(const state_type& i_q, state_type& o_dqdt, double i_t)
{
	const size_t k = static_cast<size_t>(i_t / m_dt);
	assert (k >= 0 && k < m_mu.size() && "Given mu costate array is not consistent with current integration parameters.");

	const WorkspaceIntegratedState::costate_type& mu_k = m_mu[k];
	
	const double u = mu_k[2] * m_inv_c;

	o_dqdt[0] = cos(i_q[2]);
	o_dqdt[1] = sin(i_q[2]);
	o_dqdt[2] = u;
}

}	// namespace rod2d
}	// namespace qserl
