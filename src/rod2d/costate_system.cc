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
	
#include "costate_system.h"

#include <boost/bind.hpp>

namespace qserl {
namespace rod2d {

const CostateSystem::state_type CostateSystem::kDefaultState = { {  0., 0., 0. } };

CostateSystem::CostateSystem(double i_inv_stiffness, double i_length,
	Parameters::RodModelT i_rodModel) :
	m_inv_c(i_inv_stiffness),
	m_length(i_length),
	m_rodModel(i_rodModel)
	{
		if (m_rodModel == Parameters::RM_INEXTENSIBLE)
			m_evaluationCallback = boost::bind(&CostateSystem::evaluateInextensible, this, _1, _2, _3);
		else
			assert(false && "invalid rod model");
	}

CostateSystem::~CostateSystem()
{
}

void CostateSystem::operator() (const state_type& i_mu, state_type& o_dmudt, double i_t)
{
	return m_evaluationCallback(i_mu, o_dmudt, i_t);
}

void CostateSystem::evaluateInextensible(const state_type& i_mu, state_type& o_dmudt, double i_t)
{
	const double u =  i_mu[2] * m_inv_c;

	o_dmudt[0] = i_mu[1] * u;
	o_dmudt[1] = -i_mu[0] * u;
	o_dmudt[2] = -i_mu[1];
}

}	// namespace rod2d
}	// namespace qserl
