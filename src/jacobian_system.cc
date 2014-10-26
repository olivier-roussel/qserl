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
	
#include "jacobian_system.h"

#include <boost/bind.hpp>

namespace qserl {

const JacobianSystem::state_type JacobianSystem::kDefaultState = { {  0. } };

const double JacobianSystem::kStabilityThreshold = 1.e-5;				/** Threshold for Jacobian determinant. */
//const double JacobianSystem::kStabilityThreshold = 1.e-8;				/** Threshold for Jacobian determinant. */
const double JacobianSystem::kStabilityTolerance = 1.e-12;			/** Tolerance for which Jacobian determinant vanishes. */
//const double JacobianSystem::kStabilityTolerance = 1.e-12;			/** Tolerance for which Jacobian determinant vanishes. */

JacobianSystem::JacobianSystem(const Eigen::Matrix<double, 6, 1>& i_inv_stiffness, double i_dt, 
	const std::vector<WorkspaceIntegratedState::costate_type>& i_mu, Parameters::RodModelT i_rodModel) :
	m_inv_c(i_inv_stiffness),
	m_dt(i_dt),
	m_mu(i_mu),
	m_rodModel(i_rodModel)
{
	assert (m_dt > 0. && "integration step time must be positive.");

	m_b[0] = m_inv_c[2] - m_inv_c[1];
	m_b[1] = m_inv_c[0] - m_inv_c[2];
	m_b[2] = m_inv_c[1] - m_inv_c[0];
	m_b[3] = m_inv_c[5] - m_inv_c[4];
	m_b[4] = m_inv_c[3] - m_inv_c[5];
	m_b[5] = m_inv_c[4] - m_inv_c[3];

	if (m_rodModel == Parameters::RM_INEXTENSIBLE)
		m_evaluationCallback = boost::bind(&JacobianSystem::evaluateInextensible, this, _1, _2, _3);
	else if (m_rodModel == Parameters::RM_EXTENSIBLE_SHEARABLE)
		m_evaluationCallback = boost::bind(&JacobianSystem::evaluateExtensibleShearable, this, _1, _2, _3);
	else
		assert(false && "invalid rod model");
}

JacobianSystem::~JacobianSystem()
{
}


void JacobianSystem::operator() (const state_type& i_MJ, state_type& o_dMJdt, double i_t)
{
	return m_evaluationCallback(i_MJ, o_dMJdt, i_t);
}

void JacobianSystem::evaluateInextensible(const state_type& i_MJ, state_type& o_dMJdt, double i_t)
{
	const size_t k = static_cast<size_t>(i_t / m_dt);
	assert (k >= 0 && k < m_mu.size() && "Given mu values array not consistent with current integration parameters.");

	const WorkspaceIntegratedState::costate_type& mu_k = m_mu[k];

	// pre-compute u
	const Eigen::Map<const Eigen::Matrix<double, 3, 1> > mu_k_e(mu_k.data());
	const Eigen::Matrix<double, 3, 1> u = mu_k_e.cwiseProduct(m_inv_c.block<3,1>(0,0));

	// F matrix
	Eigen::Matrix<double, 6, 6> F;
	F << 0., mu_k[2]*m_b[0], mu_k[1]*m_b[0], 0., 0., 0.,
			 mu_k[2]*m_b[1], 0, mu_k[0]*m_b[1], 0., 0., 1.,
			 mu_k[1]*m_b[2], mu_k[0]*m_b[2], 0., 0., -1., 0.,
			0., -mu_k[5]*m_inv_c[1], mu_k[4]*m_inv_c[2], 0., u[2], -u[1],
			mu_k[5]*m_inv_c[0], 0, -mu_k[3]*m_inv_c[2], -u[2], 0., u[0],
			-mu_k[4]*m_inv_c[0], mu_k[3]*m_inv_c[1], 0, u[1], -u[0], 0.;

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
	const Eigen::Map<const Eigen::Matrix<double, 6, 6> > M_e(i_MJ.data());
	const Eigen::Map<const Eigen::Matrix<double, 6, 6> > J_e(i_MJ.data()+36);
	
	// create mapping between dmjdt array and dMdt & dJdt eigen matrices
	Eigen::Map<Eigen::Matrix<double, 6, 6> > dMdt_e(o_dMJdt.data());
	Eigen::Map<Eigen::Matrix<double, 6, 6> > dJdt_e(o_dMJdt.data()+36);

	dMdt_e = F*M_e;
	dJdt_e = G*M_e + H*J_e;
}

void JacobianSystem::evaluateExtensibleShearable(const state_type& i_MJ, state_type& o_dMJdt, double i_t)
{
	const size_t k = static_cast<size_t>(i_t / m_dt);
	assert (k >= 0 && k < m_mu.size() && "Given mu values array not consistent with current integration parameters.");

	const WorkspaceIntegratedState::costate_type& mu_k = m_mu[k];

	// pre-compute u
	const Eigen::Map<const Eigen::Matrix<double, 6, 1> > mu_k_e(mu_k.data());
	const Eigen::Matrix<double, 6, 1> u = mu_k_e.cwiseProduct(m_inv_c);

	// F matrix
	Eigen::Matrix<double, 6, 6> F;
	F << 0.,								mu_k[2]*m_b[0],			mu_k[1]*m_b[0],				0.,								mu_k[5]*m_b[3],			mu_k[4]*m_b[3],
			 mu_k[2]*m_b[1],		0,									mu_k[0]*m_b[1],				mu_k[5]*m_b[4],		0.,									1+mu_k[3]*m_b[4],
			 mu_k[1]*m_b[2],		mu_k[0]*m_b[2],			0.,										mu_k[4]*m_b[5],		-1+mu_k[3]*m_b[5],	0.,
			0.,									-mu_k[5]*m_inv_c[1], mu_k[4]*m_inv_c[2],	0.,								u[2],								-u[1],
			mu_k[5]*m_inv_c[0], 0,									-mu_k[3]*m_inv_c[2],	-u[2],						0.,									u[0],
			-mu_k[4]*m_inv_c[0], mu_k[3]*m_inv_c[1], 0,										u[1],							-u[0],							0.;

	// G matrix
	Eigen::Matrix<double, 6, 6> G;
	G.setZero();
	G.diagonal() << m_inv_c[0], m_inv_c[1], m_inv_c[2], m_inv_c[3],  m_inv_c[4],  m_inv_c[5];

	// H matrix
	Eigen::Matrix<double, 6, 6> H;
	H << 0,			u[2],			-u[1],		0,			0,			0,
			-u[2],	0,				u[0],			0,			0,			0,
			u[1],		-u[0],		0,				0,			0,			0, 
			0,			u[5],			-u[4],		0,			u[2],		-u[1],
			-u[5],	0,				1+u[3],		-u[2],	0,			u[0],
			u[4],		-1-u[3],	0,				u[1],		-u[0],	0;

	// create mapping between mj array and M & J eigen matrices
	const Eigen::Map<const Eigen::Matrix<double, 6, 6> > M_e(i_MJ.data());
	const Eigen::Map<const Eigen::Matrix<double, 6, 6> > J_e(i_MJ.data()+36);
	
	// create mapping between dmjdt array and dMdt & dJdt eigen matrices
	Eigen::Map<Eigen::Matrix<double, 6, 6> > dMdt_e(o_dMJdt.data());
	Eigen::Map<Eigen::Matrix<double, 6, 6> > dJdt_e(o_dMJdt.data()+36);

	dMdt_e = F*M_e;
	dJdt_e = G*M_e + H*J_e;
}

}	// namespace qserl
