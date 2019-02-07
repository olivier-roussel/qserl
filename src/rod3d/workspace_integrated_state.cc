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

#define NOMINMAX

#include "qserl/rod3d/workspace_integrated_state.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <boost/numeric/odeint.hpp>

#include "qserl/rod3d/rod.h"
#include "full_system.h"

namespace qserl {
namespace rod3d {

/************************************************************************/
/*													Constructor																	*/
/************************************************************************/
WorkspaceIntegratedState::WorkspaceIntegratedState(unsigned int i_nnodes,
                                                   const Displacement& i_basePosition,
                                                   const Parameters& i_rodParams) :
    WorkspaceState(Displacements(), i_basePosition, i_rodParams),
    m_isInitialized{false},
    m_isStable{false},
    m_mu{},
    m_M{},
    m_J{},
    m_J_det{},
    m_J_nu_sv{},
    m_integrationOptions{} // initialize to default values
{
  assert (i_nnodes > 1 && "rod number of nodes must be greater or equal to 2");
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
WorkspaceIntegratedStateShPtr
WorkspaceIntegratedState::create(const Wrench& i_baseWrench,
                                 unsigned int i_nnodes,
                                 const Displacement& i_basePosition,
                                 const Parameters& i_rodParams)
{
  WorkspaceIntegratedStateShPtr shPtr(new WorkspaceIntegratedState(i_nnodes, i_basePosition, i_rodParams));

  if(!shPtr->init(i_baseWrench))
  {
    shPtr.reset();
  }

  return shPtr;
}

/************************************************************************/
/*														createCopy  															*/
/************************************************************************/
WorkspaceIntegratedStateShPtr
WorkspaceIntegratedState::createCopy(const WorkspaceIntegratedStateConstShPtr& i_other)
{
  WorkspaceIntegratedStateShPtr shPtr(new WorkspaceIntegratedState(*i_other));

  return shPtr;
}

/************************************************************************/
/*														init																			*/
/************************************************************************/
bool
WorkspaceIntegratedState::init(const Wrench& i_wrench)
{
  bool success = true;

  m_isStable = false;
  m_isInitialized = false;

  m_mu.resize(1);
  m_mu[0] = i_wrench;

  return success;
}

/************************************************************************/
/*														 clone																		*/
/************************************************************************/
WorkspaceStateShPtr
WorkspaceIntegratedState::clone() const
{
  return WorkspaceStateShPtr(new WorkspaceIntegratedState(*this));;
}

/************************************************************************/
/*															integrate																*/
/************************************************************************/
WorkspaceIntegratedState::IntegrationResultT
WorkspaceIntegratedState::integrate()
{
  const Wrench mu_0(Eigen::Matrix<double, 6, 1>(m_mu[0].data()));
  return integrateFromBaseWrenchRK4(mu_0);
}

/************************************************************************/
/*								integrateFromBaseWrenchRK4												*/
/************************************************************************/
WorkspaceIntegratedState::IntegrationResultT
WorkspaceIntegratedState::integrateFromBaseWrenchRK4(const Wrench& i_wrench)
{

  static const double ktstart = 0.;                          // Start integration time
  const double ktend = m_rodParameters.integrationTime;      // End integration time
  const double dt = (ktend - ktstart) / static_cast<double>(m_numNodes - 1);  // Integration time step

  m_isInitialized = true;

  if(Rod::isConfigurationSingular(i_wrench))
  {
    return IR_SINGULAR;
  }

  // 1. solve the costate system to find mu
  FullSystem full_system(m_rodParameters, dt);
  boost::numeric::odeint::runge_kutta4<FullSystem::state_type> fss_stepper;

  FullSystem::state_type x_t = FullSystem::defaultState();

  // Set initial state
  // init mu(0) = a	(base DLO wrench)
  for(int i = 0; i < 6; ++i)
  {
    (x_t.data() + FullSystem::mu_index())[i] = i_wrench[i];   // order in wrench is angular then linear
  }
  // init q_0 to identity
  Eigen::Map<Eigen::Matrix<double, 4, 4> > q_t_e(x_t.data() + FullSystem::q_index());
  q_t_e.setIdentity();
  // init M_0 to identity and J_0 to zero
  Eigen::Map<Eigen::Matrix<double, 6, 6> > M_t_e(x_t.data() + FullSystem::MJ_index());
  Eigen::Map<Eigen::Matrix<double, 6, 6> > J_t_e(x_t.data() + FullSystem::MJ_index() + 36);
  M_t_e.setIdentity();
  J_t_e.setZero();

  // setup internal memory
  m_nodes.resize(m_numNodes);
  m_nodes[0] = q_t_e;      // store q_0
  if(m_integrationOptions.keepMuValues)
  {
    m_mu.resize(m_numNodes);
    // store mu_0
    m_mu[0] = Eigen::Map<Wrench>(x_t.data() + FullSystem::mu_index());
  }
  else
  {
    m_mu.clear();
  }
  if(m_integrationOptions.keepMMatrices)
  {
    m_M.resize(m_numNodes);
    m_M[0] = M_t_e;
  }
  else
  {
    m_M.clear();
  }
  if(m_integrationOptions.keepJMatrices)
  {
    m_J.resize(m_numNodes);
    m_J[0] = J_t_e;
  }
  else
  {
    m_J.clear();
  }
  if(m_integrationOptions.keepJdet)
  {
    m_J_det.resize(m_numNodes);
    m_J_det[0] = 0.;
  }
  else
  {
    m_J_det.clear();
  }

  m_isStable = true;
  bool isThresholdOn = false;

  size_t step_idx = 1;
  double prev_det_J = 0.;
  double det_J = 0.;
  for(double t = ktstart; step_idx < m_numNodes; ++step_idx, t += dt)
  {
    fss_stepper.do_step(full_system, x_t, t, dt);
    // save state
    if(m_integrationOptions.keepMuValues)
    {
      m_mu[step_idx] = Eigen::Map<Wrench>(x_t.data() + FullSystem::mu_index());
    }
    m_nodes[step_idx] = Eigen::Map<const Eigen::Matrix4d>(x_t.data() + FullSystem::q_index());
    if(m_integrationOptions.keepMMatrices)
    {
      m_M[step_idx] = Eigen::Map<Eigen::Matrix<double, 6, 6> >(x_t.data() + FullSystem::MJ_index());
    }
    auto J_mat = Eigen::Map<Eigen::Matrix<double, 6, 6> >(x_t.data() + FullSystem::MJ_index() + 36);
    // check stability
    prev_det_J = det_J;
    det_J =  Eigen::Map<Eigen::Matrix<double, 6, 6> >(x_t.data() + FullSystem::MJ_index() + 36).determinant();
    if(m_integrationOptions.keepJMatrices)
    {
      m_J[step_idx] = J_mat;
    }
    if(m_integrationOptions.keepJdet)
    {
      m_J_det[step_idx] = det_J;
    }
    if(std::abs(det_J) > full_system.jacobianStabilityThreshold())
    {
      isThresholdOn = true;
    }
    if(isThresholdOn and (std::abs(det_J) < full_system.jacobianStabilityTolerance() or
      det_J * prev_det_J < 0.))
    {  // zero crossing
      m_isStable = false;
    }
  }

  // compute J nu part singular values
  if((!m_integrationOptions.stop_if_unstable || m_isStable) && m_integrationOptions.computeJ_nu_sv)
  {
    m_J_nu_sv.assign(m_numNodes, Eigen::Vector3d::Zero());
    for(size_t idxNode = 1; idxNode < m_numNodes; ++idxNode)
    {
      Eigen::JacobiSVD<Eigen::Matrix<double, 3, 6> > svd_J_nu(m_J[idxNode].block<3, 6>(3, 0));
      m_J_nu_sv[idxNode] = svd_J_nu.singularValues();
    }
  }

  if(not m_isStable)
  {
    return IR_UNSTABLE;
  }

  return IR_VALID;
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
bool
WorkspaceIntegratedState::isStable() const
{
  assert(m_isInitialized && "the state must be integrated first");
  return m_isStable;
}

/************************************************************************/
/*																baseWrench														*/
/************************************************************************/
Wrench
WorkspaceIntegratedState::baseWrench() const
{
  assert(m_isInitialized && "the state must be integrated first");
  return m_mu.front();
}

/************************************************************************/
/*																	tipWrench														*/
/************************************************************************/
Wrench
WorkspaceIntegratedState::tipWrench() const
{
  assert(m_isInitialized && "the state must be integrated first");
  return m_mu.back();
}

/************************************************************************/
/*																wrench																*/
/************************************************************************/
Wrench
WorkspaceIntegratedState::wrench(size_t i_idxNode) const
{
  assert(m_isInitialized && "the state must be integrated first");
  assert(i_idxNode < m_mu.size() && "invalid node index");
  return m_mu[i_idxNode];
}

/************************************************************************/
/*																mu																		*/
/************************************************************************/
const Wrenches&
WorkspaceIntegratedState::mu() const
{
  assert(m_isInitialized && "the state must be integrated first");
  return m_mu;
}

/************************************************************************/
/*																getMMatrix()													*/
/************************************************************************/
const Matrix6d&
WorkspaceIntegratedState::getMMatrix(size_t i_nodeIdx) const
{
  assert(m_isInitialized && "the state must be integrated first");
  assert(i_nodeIdx < m_numNodes && "invalid node index");
  return m_M[i_nodeIdx];
}

/************************************************************************/
/*																getJMatrix()													*/
/************************************************************************/
const Matrix6d&
WorkspaceIntegratedState::getJMatrix(size_t i_nodeIdx) const
{
  assert(m_isInitialized && "the state must be integrated first");
  assert(i_nodeIdx < m_numNodes && "invalid node index");
  return m_J[i_nodeIdx];
}

/************************************************************************/
/*																J_det																	*/
/************************************************************************/
const std::vector<double>&
WorkspaceIntegratedState::J_det() const
{
  assert(m_isInitialized && "the state must be integrated first");
  return m_J_det;
}

/************************************************************************/
/*																J_nu_sv																*/
/************************************************************************/
const Eigen::Vector3d&
WorkspaceIntegratedState::J_nu_sv(size_t i_nodeIdx) const
{
  assert(i_nodeIdx < m_numNodes && "invalid node index");
  return m_J_nu_sv[i_nodeIdx];
}

/************************************************************************/
/*											computeLinearizedNodePositions									*/
/************************************************************************/
//WorkspaceStateShPtr
//WorkspaceIntegratedState::approximateLinearlyNeighbourState(const Wrench& i_da,
//                                                            const Displacement& i_neighbBase) const
//{
//  //assert (m_state && "current state must be initialized to compute linearized positions. ");
//
//  WorkspaceStateShPtr neighbApproxState = WorkspaceState::create(Displacements(m_numNodes),
//                                                                 i_neighbBase,
//                                                                 m_rodParameters);
//
//  for(int nodeIdx = 0; nodeIdx < static_cast<int>(m_numNodes); ++nodeIdx)
//  {
//    const Eigen::Matrix<double, 6, 6>& J_mat = getJMatrix(nodeIdx);
//    const Eigen::Twistd xi = J_mat * i_da;
//    Displacement localDisp = xi.exp();
//    // apply DLO scale
//    localDisp.getTranslation() = localDisp.getTranslation() * m_rodParameters.length;
//    neighbApproxState->m_nodes[nodeIdx] = m_nodes[nodeIdx] * localDisp;
//  }
//  return neighbApproxState;
//}


/************************************************************************/
/*																	memUsage														*/
/************************************************************************/
size_t
WorkspaceIntegratedState::memUsage() const
{
  return WorkspaceState::memUsage() +
         sizeof(m_isInitialized) +
         sizeof(m_isStable) +
         m_mu.capacity() * sizeof(Wrench) +
         m_M.capacity() * sizeof(Matrix6d) +
         m_J.capacity() * sizeof(Matrix6d) +
         m_J_det.capacity() * sizeof(double) +
         m_J_nu_sv.capacity() * sizeof(Eigen::Vector3d) +
         sizeof(m_integrationOptions);
}

/************************************************************************/
/*												integrationOptions														*/
/************************************************************************/
void
WorkspaceIntegratedState::integrationOptions(const WorkspaceIntegratedState::IntegrationOptions& i_integrationOptions)
{
  m_integrationOptions = i_integrationOptions;
}

/************************************************************************/
/*												integrationOptions														*/
/************************************************************************/
const WorkspaceIntegratedState::IntegrationOptions&
WorkspaceIntegratedState::integrationOptions() const
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
{
}

}  // namespace rod3d
}  // namespace qserl
