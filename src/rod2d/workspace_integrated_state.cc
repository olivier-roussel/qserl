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

#define NOMINMAX

#include "qserl/rod2d/workspace_integrated_state.h"

#pragma warning( push, 0 )	
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#pragma warning( pop )	
#include <boost/numeric/odeint.hpp>

#include "qserl/rod2d/rod.h"
#include "state_system.h"
#include "costate_system.h"
#include "jacobian_system.h"

namespace qserl {
namespace rod2d {

/************************************************************************/
/*													Constructor																	*/
/************************************************************************/
WorkspaceIntegratedState::WorkspaceIntegratedState(unsigned int i_nnodes, 
	const Displacement2D& i_basePosition,	const Parameters& i_rodParams):
WorkspaceState(std::vector<Displacement2D>(), i_basePosition, i_rodParams),
	m_integrationOptions() // initialize to default values
{
  assert ( i_nnodes > 1 && "rod number of nodes must be greater or equal to 2" );
	m_numNodes = i_nnodes;
}

/************************************************************************/
/*												 Destructor																		*/
/************************************************************************/
WorkspaceIntegratedState::~WorkspaceIntegratedState()
{
}

/************************************************************************/
/*														create																		*/
/************************************************************************/
WorkspaceIntegratedStateShPtr WorkspaceIntegratedState::create(const Wrench2D& i_baseWrench, unsigned int i_nnodes, 
	const Displacement2D& i_basePosition, const Parameters& i_rodParams)
{
	WorkspaceIntegratedStateShPtr shPtr(new WorkspaceIntegratedState(i_nnodes, i_basePosition, i_rodParams));

	if (!shPtr->init(i_baseWrench))
		shPtr.reset();

	return shPtr;
}

/************************************************************************/
/*														createCopy  															*/
/************************************************************************/
WorkspaceIntegratedStateShPtr WorkspaceIntegratedState::createCopy(const WorkspaceIntegratedStateConstShPtr& i_other)
{
	WorkspaceIntegratedStateShPtr shPtr(new WorkspaceIntegratedState(*i_other));

	return shPtr;
}

/************************************************************************/
/*														init																			*/
/************************************************************************/
bool WorkspaceIntegratedState::init(const Wrench2D& i_wrench)
{
	bool success = true;

	m_isStable = false;

	m_mu.resize(1);
	m_mu[0] = i_wrench;

	return success;
}

/************************************************************************/
/*														 clone																		*/
/************************************************************************/
WorkspaceStateShPtr WorkspaceIntegratedState::clone() const
{
	return WorkspaceStateShPtr(new WorkspaceIntegratedState(*this));;
}

/************************************************************************/
/*															integrate																*/
/************************************************************************/
bool WorkspaceIntegratedState::integrate()
{
	return integrateFromBaseWrench(m_mu[0]);
}

/************************************************************************/
/*												integrateFromBaseWrench												*/
/************************************************************************/
bool WorkspaceIntegratedState::integrateFromBaseWrench(const Wrench2D& i_wrench)
{
	static const double ktstart = 0.;													// Start integration time
	const double ktend = m_rodParameters.integrationTime;			// End integration time
	const double dt = (ktend - ktstart) / static_cast<double>(m_numNodes-1);	// Integration time step

	if (Rod::isConfigurationSingular(i_wrench))
		return false;

	// 1. solve the costate system to find mu
	const double stiffnessCoefficient = Rod::getStiffnessCoefficients(m_rodParameters);
	const double invStiffness = 1. / stiffnessCoefficient;
	CostateSystem costate_system(invStiffness, m_rodParameters.length, m_rodParameters.rodModel);

	// init mu(0) = a					(base DLO wrench)
	costate_type mu_t = i_wrench;

	// until C++0x, there is no convenient way to release memory of a vector 
	// so we use temporary allocated vectors if we do not want our instance
	// memory print to explode.
	std::vector<costate_type>* mu_buffer;
	if (m_integrationOptions.keepMuValues)
	{
		mu_buffer = &m_mu;
		mu_buffer->assign(m_numNodes, CostateSystem::kDefaultState);
		m_mu[0] = mu_t;				// store mu_0
	}else{
		mu_buffer = new std::vector<costate_type>(m_numNodes, CostateSystem::kDefaultState);
		m_mu.assign(1, mu_t); // store mu_0
	}

	// integrator for costates mu
	boost::numeric::odeint::runge_kutta4< costate_type > css_stepper;

	size_t step_idx = 1;
	for (double t = ktstart; step_idx < m_numNodes ; ++step_idx, t+=dt)
	{
		css_stepper.do_step(costate_system, mu_t, t, dt);
		(*mu_buffer)[step_idx] = mu_t;
	}

	// 2. solve the state system to find q 
	StateSystem state_system(invStiffness, m_rodParameters.length, dt, *mu_buffer, m_rodParameters.rodModel);
	boost::numeric::odeint::runge_kutta4< state_type > sss_stepper;
	m_nodes.resize(m_numNodes);

	// init q_0 to identity
	state_type q_t = StateSystem::kDefaultState;
	m_nodes[0] = q_t;

	step_idx = 1;
	for (double t = ktstart; step_idx < m_numNodes ; ++step_idx, t+=dt)
	{
		sss_stepper.do_step(state_system, q_t, t, dt);
		m_nodes[step_idx] = q_t;
	}

	// 3. Solve the jacobian system (and check non-degenerescence of matrix J)
	JacobianSystem jacobianSystem(invStiffness, dt, *mu_buffer, m_rodParameters.rodModel);
	boost::numeric::odeint::runge_kutta4< jacobian_state_type > jacobianStepper;
	std::vector<Eigen::Matrix<double, 3, 3> >* M_buffer;
	if (m_integrationOptions.keepMMatrices)
	{
		M_buffer = &m_M;
		M_buffer->assign(m_numNodes, Eigen::Matrix<double, 3, 3>::Zero());
	}else{
		M_buffer = new std::vector<Eigen::Matrix<double, 3, 3> >(m_numNodes, Eigen::Matrix<double, 3, 3>::Zero());
	}

	std::vector<Eigen::Matrix<double, 3, 3> >* J_buffer;
	if (m_integrationOptions.keepJMatrices)
	{
		J_buffer = &m_J;
		J_buffer->assign(m_numNodes, Eigen::Matrix<double, 3, 3>::Zero());
	}else{
		J_buffer = new std::vector<Eigen::Matrix<double, 3, 3> >(m_numNodes, Eigen::Matrix<double, 3, 3>::Zero());
	}

	// init M_0 to identity and J_0 to zero
	jacobian_state_type jacobian_t;
	Eigen::Map<Eigen::Matrix<double, 3, 3> > M_t_e(jacobian_t.data());
	Eigen::Map<Eigen::Matrix<double, 3, 3> > J_t_e(jacobian_t.data()+9);
	M_t_e.setIdentity();
	J_t_e.setZero();
	(*M_buffer)[0] = M_t_e;
	(*J_buffer)[0] = J_t_e;

	step_idx = 1;
	m_isStable = true;
	bool isThresholdOn = false;

	std::vector<double>* J_det_buffer;
	if (m_integrationOptions.keepJdet)
	{
		J_det_buffer = &m_J_det;
		J_det_buffer->assign(m_numNodes, 0.);
	}else
		J_det_buffer = new std::vector<double>(m_numNodes, 0.);

	for (double t = ktstart; step_idx < m_numNodes && (!m_integrationOptions.stop_if_unstable || m_isStable) ; 
		++step_idx, t+=dt)
	{
		jacobianStepper.do_step(jacobianSystem, jacobian_t, t, dt);
		(*M_buffer)[step_idx] = Eigen::Map<Eigen::Matrix<double, 3, 3> >(jacobian_t.data());
		(*J_buffer)[step_idx] = Eigen::Map<Eigen::Matrix<double, 3, 3> >(jacobian_t.data()+9);
		// check if stable
		double& J_det = (*J_det_buffer)[step_idx];
		J_det = (*J_buffer)[step_idx].determinant();
		if (abs(J_det) > JacobianSystem::kStabilityThreshold)
			isThresholdOn = true;
		if (isThresholdOn && ( abs(J_det) < JacobianSystem::kStabilityTolerance ||
			J_det * (*J_det_buffer)[step_idx-1] < 0.) )	// zero crossing
			m_isStable = false;
	}

	if (!m_integrationOptions.keepMuValues)
		delete mu_buffer;

	if (!m_integrationOptions.keepJdet)
		delete J_det_buffer;

	if (!m_integrationOptions.keepMMatrices)
		delete M_buffer;

	if (!m_integrationOptions.keepJMatrices)
		delete J_buffer;

	return true; 
}

/************************************************************************/
/*																isStable															*/
/************************************************************************/
bool WorkspaceIntegratedState::isStable() const
{
	return m_isStable;
}

/************************************************************************/
/*																wrench																*/
/************************************************************************/
Wrench2D WorkspaceIntegratedState::wrench(size_t i_idxNode) const
{
	return m_mu[i_idxNode];
}

/************************************************************************/
/*																mu																		*/
/************************************************************************/
const std::vector<WorkspaceIntegratedState::costate_type>& WorkspaceIntegratedState::mu() const
{
	return m_mu;
}

/************************************************************************/
/*																getMMatrix()													*/
/************************************************************************/
const Eigen::Matrix<double, 3, 3>& WorkspaceIntegratedState::getMMatrix(size_t i_nodeIdx) const
{
	assert (i_nodeIdx >= 0 && i_nodeIdx < m_numNodes && "invalid node index");
	return m_M[i_nodeIdx];
}

/************************************************************************/
/*																getJMatrix()													*/
/************************************************************************/
const Eigen::Matrix<double, 3, 3>& WorkspaceIntegratedState::getJMatrix(size_t i_nodeIdx) const
{
	assert (i_nodeIdx >= 0 && i_nodeIdx < m_numNodes && "invalid node index");
	return m_J[i_nodeIdx];
}

/************************************************************************/
/*																J_det																	*/
/************************************************************************/
const std::vector<double>& WorkspaceIntegratedState::J_det() const
{
	return m_J_det;
}

/************************************************************************/
/*																	memUsage														*/
/************************************************************************/
size_t WorkspaceIntegratedState::memUsage() const
{
	return WorkspaceState::memUsage() +
		sizeof(m_isStable) + 
		m_mu.capacity() * sizeof(costate_type) + 
		m_M.capacity() * sizeof(Eigen::Matrix<double, 3, 3>) + 
		m_J.capacity() * sizeof(Eigen::Matrix<double, 3, 3>) + 
		m_J_det.capacity() * sizeof(double) +
		sizeof(m_integrationOptions);
}

/************************************************************************/
/*												integrationOptions														*/
/************************************************************************/
void WorkspaceIntegratedState::integrationOptions(const WorkspaceIntegratedState::IntegrationOptions& i_integrationOptions)
{
	m_integrationOptions = i_integrationOptions;
}

/************************************************************************/
/*												integrationOptions														*/
/************************************************************************/
const WorkspaceIntegratedState::IntegrationOptions& WorkspaceIntegratedState::integrationOptions() const
{
	return m_integrationOptions;
}

/************************************************************************/
/*										integrateUntilConjugatePoint											*/
/************************************************************************/
bool WorkspaceIntegratedState::integrateUntilConjugatePoint(int i_maxIter, double& o_tconj)
{
	static const double ktstart = 0.;													// Start integration time
	const double ktend = m_rodParameters.integrationTime;			// End integration time
	const double dt = (ktend - ktstart) / static_cast<double>(m_numNodes-1);	// Integration time step

	if (Rod::isConfigurationSingular(m_mu[0]))
		return false;

	const double stiffnessCoefficient = Rod::getStiffnessCoefficients(m_rodParameters);
	const double invStiffness = 1. / stiffnessCoefficient;
	CostateSystem costate_system(invStiffness, m_rodParameters.length, m_rodParameters.rodModel);
	boost::numeric::odeint::runge_kutta4< costate_type > css_stepper;

	// init mu(0) = a					(base DLO wrench)
	costate_type mu_t = m_mu[0];

		// until C++0x, there is no convenient way to release memory of a vector 
	// so we use temporary allocated vectors if we do not want our instance
	// memory print to explode.
	std::vector<costate_type>* mu_buffer;
	if (m_integrationOptions.keepMuValues)
	{
		mu_buffer = &m_mu;
	}else{
		mu_buffer = new std::vector<costate_type>();
		(*mu_buffer).push_back(m_mu[0]);
	}

	// init state integrator and q(0)
	StateSystem state_system(invStiffness, m_rodParameters.length, dt, *mu_buffer, m_rodParameters.rodModel);
	boost::numeric::odeint::runge_kutta4< state_type > sss_stepper;
	state_type q_t = StateSystem::kDefaultState;
	m_nodes.clear();
	m_nodes.push_back(q_t);
	//m_nodes.reserve(m_numNodes);

	// init jacobian integarator and M_0 to identity and J_0 to zero
	JacobianSystem jacobianSystem(invStiffness, dt, *mu_buffer, m_rodParameters.rodModel);
	boost::numeric::odeint::runge_kutta4< jacobian_state_type > jacobianStepper;
	jacobian_state_type jacobian_t;
	Eigen::Map<Eigen::Matrix<double, 3, 3> > M_t_e(jacobian_t.data());
	Eigen::Map<Eigen::Matrix<double, 3, 3> > J_t_e(jacobian_t.data()+9);
	M_t_e.setIdentity();
	J_t_e.setZero();
	if (m_integrationOptions.keepMMatrices)
		m_M.push_back(M_t_e);
	if (m_integrationOptions.keepJMatrices)
		m_J.push_back(J_t_e);
	double Jdet_cur = 0., Jdet_prev = 0.;
	m_J_det.clear();
	m_J_det.push_back(Jdet_cur);

	// integrate until unstability or max iteration reached
	m_isStable = true;
	bool isThresholdOn = false;
	int iter = 1;
	double t = 0.;
	while (m_isStable && iter < i_maxIter)
	{
		// integrate co-state
		css_stepper.do_step(costate_system, mu_t, t, dt);
		(*mu_buffer).push_back(mu_t);
		// integrate jacobian
		jacobianStepper.do_step(jacobianSystem, jacobian_t, t, dt);
		Eigen::Map<Eigen::Matrix<double, 3, 3> > J_cur = Eigen::Map<Eigen::Matrix<double, 3, 3> >(jacobian_t.data()+9);
		// compute jacobian and check stability
		Jdet_prev = Jdet_cur;
		Jdet_cur = J_cur.determinant();
		if (abs(Jdet_cur) > JacobianSystem::kStabilityThreshold)
				isThresholdOn = true;
		if (isThresholdOn && ( abs(Jdet_cur) < JacobianSystem::kStabilityTolerance ||
			Jdet_cur * Jdet_prev < 0.) )	// zero crossing
				m_isStable = false;
		if (m_isStable)
		{
			// integrate state
			sss_stepper.do_step(state_system, q_t, t, dt);
			if (m_integrationOptions.keepMMatrices)
				m_M.push_back(Eigen::Map<Eigen::Matrix<double, 3, 3> >(jacobian_t.data()));
			if (m_integrationOptions.keepJMatrices)
				m_J.push_back(J_cur);
			if (m_integrationOptions.keepJdet)
				m_J_det.push_back(Jdet_cur);
			m_nodes.push_back(q_t);
			++iter;
			t+=dt;
		}
	}

	if (!m_integrationOptions.keepMuValues)
		delete mu_buffer;

	if (!m_isStable)
	{
		// conjugate point found
		o_tconj = t;
	}else{
		// conjugate point not found
		o_tconj = -1.;
	}
	return true;
}

/************************************************************************/
/*									IntegrationOptions::Constructor											*/
/************************************************************************/
WorkspaceIntegratedState::IntegrationOptions::IntegrationOptions() :
	stop_if_unstable(true),
	keepMuValues(false),
	keepJdet(false),
	keepMMatrices(false),
	keepJMatrices(false)
{}

}	// namespace rod2d
}	// namespace qserl
