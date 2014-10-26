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

#include "qserl/workspace_state.h"

#include <Eigen/Core>

namespace qserl {

/************************************************************************/
/*													Constructor																	*/
/************************************************************************/
WorkspaceState::WorkspaceState(const std::vector<Eigen::Displacementd>& i_nodes, const Eigen::Displacementd& i_basePosition, const Parameters& i_rodParams):
	m_numNodes(i_nodes.size()),
	m_nodes(i_nodes),
	m_base(i_basePosition),
	m_rodParameters(i_rodParams)
{
	//assert(m_numNodes > 0 && "number of rod nodes must be strictly positive");
}

/************************************************************************/
/*												 Destructor																		*/
/************************************************************************/
WorkspaceState::~WorkspaceState()
{
}

/************************************************************************/
/*														create																		*/
/************************************************************************/
WorkspaceStateShPtr WorkspaceState::create(const std::vector<Eigen::Displacementd>& i_nodes, const Eigen::Displacementd& i_basePosition, const Parameters& i_rodParams)
{
	WorkspaceStateShPtr shPtr(new WorkspaceState(i_nodes, i_basePosition, i_rodParams));

	//if (!shPtr->init(shPtr))
	//	shPtr.reset();

	return shPtr;
}

/************************************************************************/
/*														 clone																		*/
/************************************************************************/
WorkspaceStateShPtr WorkspaceState::clone() const
{
	return WorkspaceStateShPtr(new WorkspaceState(*this));;
}

/************************************************************************/
/*														init																			*/
/************************************************************************/
//bool WorkspaceState::init(const WorkspaceStateWkPtr& i_weakPtr)
//{
//	bool success = true;
//
//	if (success)
//		m_weakPtr = i_weakPtr;
//
//	return success;
//}

/************************************************************************/
/*																numNodes															*/
/************************************************************************/
size_t WorkspaceState::numNodes() const
{
	return m_numNodes;
}

/************************************************************************/
/*																	nodes																*/
/************************************************************************/
const std::vector<Eigen::Displacementd>& WorkspaceState::nodes() const
{
	return m_nodes;
}

/************************************************************************/
/*																base																	*/
/************************************************************************/
const Eigen::Displacementd& WorkspaceState::base() const
{
	return m_base;
}

/************************************************************************/
/*																base																	*/
/************************************************************************/
void WorkspaceState::base(const Eigen::Displacementd& i_base)
{
	m_base = i_base;
}

/************************************************************************/
/*														nodesAbsolutePositions										*/
/************************************************************************/
std::vector<Eigen::Vector3d> WorkspaceState::nodesAbsolutePositions() const
{
	std::vector<Eigen::Vector3d> pos;
	pos.reserve(m_numNodes);
	for (size_t i = 0 ; i < m_numNodes ; ++i)
		pos.push_back(m_base * m_nodes[i].getTranslation());
	return pos;
}

/************************************************************************/
/*														nodesAbsolute6DPositions										*/
/************************************************************************/
std::vector<Eigen::Displacementd> WorkspaceState::nodesAbsolute6DPositions() const
{
  std::vector<Eigen::Displacementd> pos;
	pos.reserve(m_numNodes);
	for (size_t i = 0 ; i < m_numNodes ; ++i)
		pos.push_back(m_base * m_nodes[i]);
	return pos;
}

/************************************************************************/
/*													staticParameters														*/
/************************************************************************/
const Parameters& WorkspaceState::staticParameters() const
{
	return m_rodParameters;
}

/************************************************************************/
/*																	memUsage														*/
/************************************************************************/
size_t WorkspaceState::memUsage() const
{
	return sizeof(m_numNodes) + 
		m_nodes.capacity() * sizeof(Eigen::Displacementd) + 
		sizeof(m_base) + 
		sizeof(m_rodParameters);/* + 
		sizeof(m_weakPtr);*/
}

/************************************************************************/
/*													torsionalRotation														*/
/************************************************************************/
double WorkspaceState::torsionalRotation() const
{
	double torsionalRotation = 0.;
	for (size_t idxNode = 0 ; idxNode < m_numNodes - 1 ; ++idxNode)
	{
		const Eigen::AngularVelocityd w = (m_nodes[idxNode].getRotation().inverse() * m_nodes[idxNode+1].getRotation()).log();
		torsionalRotation += w[0];
	}
	return torsionalRotation;
}

}	// namespace qserl
