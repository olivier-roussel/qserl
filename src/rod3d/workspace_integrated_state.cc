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

#include "qserl/rod3d/workspace_integrated_state.h"

#pragma warning( push, 0 )	
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#pragma warning( pop )	
#include <boost/numeric/odeint.hpp>

#include "qserl/rod3d/rod.h"
#include "state_system.h"
#include "costate_system.h"
#include "jacobian_system.h"

namespace qserl {
namespace rod3d {

/************************************************************************/
/*													Constructor																	*/
/************************************************************************/
WorkspaceIntegratedState::WorkspaceIntegratedState(unsigned int i_nnodes, const Eigen::Displacementd& i_basePosition,
	const Parameters& i_rodParams):
WorkspaceState(std::vector<Eigen::Displacementd>(), i_basePosition, i_rodParams),
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
WorkspaceIntegratedStateShPtr WorkspaceIntegratedState::create(const Eigen::Wrenchd& i_baseWrench, unsigned int i_nnodes, 
	const Eigen::Displacementd& i_basePosition, const Parameters& i_rodParams)
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
bool WorkspaceIntegratedState::init(const Eigen::Wrenchd& i_wrench)
{
	bool success = true;

	m_isStable = false;

	m_mu.resize(1);
	for (int i = 0 ; i < 6 ; ++i)
		m_mu[0][i] = i_wrench[i];		// XXX for an Eigen::Wrench, torque are the 3 first components and force the 3 last

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
	Eigen::Wrenchd mu_0(Eigen::Matrix<double, 6, 1>(m_mu[0].data()));
	return integrateFromBaseWrench(mu_0);
}

/************************************************************************/
/*												integrateFromBaseWrench												*/
/************************************************************************/
bool WorkspaceIntegratedState::integrateFromBaseWrench(const Eigen::Wrenchd& i_wrench)
{

	static const double ktstart = 0.;													// Start integration time
	const double ktend = m_rodParameters.integrationTime;			// End integration time
	//const double ktend = 1.;
	const double dt = (ktend - ktstart) / static_cast<double>(m_numNodes-1);	// Integration time step

	if (Rod::isConfigurationSingular(i_wrench))
		return false;

	// 1. solve the costate system to find mu
	const Eigen::Matrix<double, 6, 1> stiffnessCoefficients = Rod::getStiffnessCoefficients(m_rodParameters);
	const Eigen::Matrix<double, 6, 1> invStiffness = stiffnessCoefficients.cwiseInverse();
	CostateSystem costate_system(invStiffness, m_rodParameters.length, m_rodParameters.rodModel);

	costate_type mu_t;
	// init mu(0) = a					(base DLO wrench)
	for (int i = 0 ; i < 6 ; ++i)
		mu_t[i] = i_wrench[i];		// XXX for an Eigen::Wrench, torque are the 3 first components and force the 3 last, so we can just copy simply

	// until C++0x, there is no convenient way to release memory of a vector 
	// so we use temporary allocated vectors if we do not want our instance
	// memory print to explode.
	std::vector<costate_type>* mu_buffer;
	if (m_integrationOptions.keepMuValues)
	{
		mu_buffer = &m_mu;
		mu_buffer->assign(m_numNodes, CostateSystem::defaultState());
		m_mu[0] = mu_t;				// store mu_0
	}else{
		mu_buffer = new std::vector<costate_type>(m_numNodes, CostateSystem::defaultState());
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
	std::vector<state_type> q_array(m_numNodes, StateSystem::defaultState());

	// init q_0 to identity
	state_type q_t;
	Eigen::Matrix4d::Map(q_t.data()).setIdentity();
	q_array[0] = q_t;			// store q_0

	step_idx = 1;
	for (double t = ktstart; step_idx < m_numNodes ; ++step_idx, t+=dt)
	{
		sss_stepper.do_step(state_system, q_t, t, dt);
		q_array[step_idx] = q_t;
	}

	// update state
	m_nodes.assign(q_array.size(), Eigen::Displacementd::Identity());
	for (size_t i = 0 ; i < q_array.size() ; ++i)
	{
		const Eigen::Map<const Eigen::Matrix4d> q_e(q_array[i].data());
		m_nodes[i] = Eigen::Displacementd(q_e);
	}

	// 3. Solve the jacobian system (and check non-degenerescence of matrix J)
	JacobianSystem jacobianSystem(invStiffness, dt, *mu_buffer, m_rodParameters.rodModel);
	boost::numeric::odeint::runge_kutta4< jacobian_state_type > jacobianStepper;
	std::vector<Eigen::Matrix<double, 6, 6> >* M_buffer;
	if (m_integrationOptions.keepMMatrices)
	{
		M_buffer = &m_M;
		M_buffer->assign(m_numNodes, Eigen::Matrix<double, 6, 6>::Zero());
	}else{
		M_buffer = new std::vector<Eigen::Matrix<double, 6, 6> >(m_numNodes, Eigen::Matrix<double, 6, 6>::Zero());
	}

	std::vector<Eigen::Matrix<double, 6, 6> >* J_buffer;
	if (m_integrationOptions.keepJMatrices)
	{
		J_buffer = &m_J;
		J_buffer->assign(m_numNodes, Eigen::Matrix<double, 6, 6>::Zero());
	}else{
		J_buffer = new std::vector<Eigen::Matrix<double, 6, 6> >(m_numNodes, Eigen::Matrix<double, 6, 6>::Zero());
	}

	// init M_0 to identity and J_0 to zero
	jacobian_state_type jacobian_t;
	Eigen::Map<Eigen::Matrix<double, 6, 6> > M_t_e(jacobian_t.data());
	Eigen::Map<Eigen::Matrix<double, 6, 6> > J_t_e(jacobian_t.data()+36);
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
		(*M_buffer)[step_idx] = Eigen::Map<Eigen::Matrix<double, 6, 6> >(jacobian_t.data());
		(*J_buffer)[step_idx] = Eigen::Map<Eigen::Matrix<double, 6, 6> >(jacobian_t.data()+36);
		// check if stable
		double& J_det = (*J_det_buffer)[step_idx];
		J_det = (*J_buffer)[step_idx].determinant();
		if (abs(J_det) > JacobianSystem::kStabilityThreshold)
			isThresholdOn = true;
		if (isThresholdOn && ( abs(J_det) < JacobianSystem::kStabilityTolerance ||
			J_det * (*J_det_buffer)[step_idx-1] < 0.) )	// zero crossing
			m_isStable = false;
	}

	// compute J nu part singular values
	if ((!m_integrationOptions.stop_if_unstable || m_isStable) && m_integrationOptions.computeJ_nu_sv)
	{
		m_J_nu_sv.assign(m_numNodes, Eigen::Vector3d::Zero());
		for (size_t idxNode = 1 ; idxNode < m_numNodes ; ++idxNode)
		{
			Eigen::JacobiSVD<Eigen::Matrix<double, 3, 6> > svd_J_nu((*J_buffer)[idxNode].block<3,6>(3,0));
			m_J_nu_sv[idxNode] = svd_J_nu.singularValues();
		}
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
/*									computeJacobianNuSingularValues											*/
/************************************************************************/
//void WorkspaceIntegratedState::computeJacobianNuSingularValues(bool i_compute)
//{
//	m_computeJ_nu_sv = i_compute;
//}

/************************************************************************/
/*									computeJacobianNuSingularValues											*/
/************************************************************************/
//bool WorkspaceIntegratedState::computeJacobianNuSingularValues() const
//{
//	return m_computeJ_nu_sv;
//}

/************************************************************************/
/*																isStable															*/
/************************************************************************/
bool WorkspaceIntegratedState::isStable() const
{
	return m_isStable;
}

/************************************************************************/
/*																baseWrench														*/
/************************************************************************/
Eigen::Wrenchd WorkspaceIntegratedState::baseWrench() const
{
	const Eigen::Map<const Eigen::Matrix<double, 6, 1> > mu_0(m_mu[0].data());
	return Eigen::Wrenchd(mu_0);
}

/************************************************************************/
/*																	tipWrench														*/
/************************************************************************/
Eigen::Wrenchd WorkspaceIntegratedState::tipWrench() const
{
	const Eigen::Map<const Eigen::Matrix<double, 6, 1> > mu_n(m_mu.back().data());
	return Eigen::Wrenchd(mu_n);
}

/************************************************************************/
/*																wrench																*/
/************************************************************************/
Eigen::Wrenchd WorkspaceIntegratedState::wrench(size_t i_idxNode) const
{
	const Eigen::Map<const Eigen::Matrix<double, 6, 1> > mu(m_mu[i_idxNode].data());
	return Eigen::Wrenchd(mu);
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
const Eigen::Matrix<double, 6, 6>& WorkspaceIntegratedState::getMMatrix(size_t i_nodeIdx) const
{
	assert (i_nodeIdx >= 0 && i_nodeIdx < m_numNodes && "invalid node index");
	return m_M[i_nodeIdx];
}

/************************************************************************/
/*																getJMatrix()													*/
/************************************************************************/
const Eigen::Matrix<double, 6, 6>& WorkspaceIntegratedState::getJMatrix(size_t i_nodeIdx) const
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
/*																J_nu_sv																*/
/************************************************************************/
const Eigen::Vector3d& WorkspaceIntegratedState::J_nu_sv(size_t i_nodeIdx) const
{
	assert (i_nodeIdx >= 0 && i_nodeIdx < m_numNodes && "invalid node index");
	return m_J_nu_sv[i_nodeIdx];
}

/************************************************************************/
/*											computeLinearizedNodePositions									*/
/************************************************************************/
WorkspaceStateShPtr WorkspaceIntegratedState::approximateLinearlyNeighbourState(const Eigen::Wrenchd& i_da, const Eigen::Displacementd& i_neighbBase) const
{
	//assert (m_state && "current state must be initialized to compute linearized positions. ");

	WorkspaceStateShPtr neighbApproxState = WorkspaceState::create(std::vector<Eigen::Displacementd>(m_numNodes), i_neighbBase, m_rodParameters);

	for (int nodeIdx = 0 ; nodeIdx < static_cast<int>(m_numNodes) ; ++nodeIdx)
	{
		const Eigen::Matrix<double, 6, 6>& J_mat = getJMatrix(nodeIdx);
		const Eigen::Twistd xi = J_mat * i_da;
		Eigen::Displacementd localDisp = xi.exp();
		// apply DLO scale
		localDisp.getTranslation() = localDisp.getTranslation() * m_rodParameters.length;
		neighbApproxState->m_nodes[nodeIdx] = m_nodes[nodeIdx] * localDisp;
	}
	return neighbApproxState;
}


/************************************************************************/
/*																	memUsage														*/
/************************************************************************/
size_t WorkspaceIntegratedState::memUsage() const
{
	return WorkspaceState::memUsage() +
		sizeof(m_isStable) + 
		m_mu.capacity() * sizeof(costate_type) + 
		//m_MJ.capacity() * sizeof(JacobianSystem::state_type) + 
		m_M.capacity() * sizeof(Eigen::Matrix<double, 6, 6>) + 
		m_J.capacity() * sizeof(Eigen::Matrix<double, 6, 6>) + 
		m_J_det.capacity() * sizeof(double) +
		m_J_nu_sv.capacity() * sizeof(Eigen::Vector3d) +
		sizeof(m_integrationOptions);/* +
		sizeof(m_weakPtr);*/
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
/*									IntegrationOptions::Constructor											*/
/************************************************************************/
WorkspaceIntegratedState::IntegrationOptions::IntegrationOptions() :
	computeJ_nu_sv(false),
	stop_if_unstable(true),
	keepMuValues(false),
	keepJdet(false),
	keepMMatrices(false),
	keepJMatrices(true)
{}

}	// namespace rod3d
}	// namespace qserl
