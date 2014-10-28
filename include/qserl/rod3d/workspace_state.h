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

#ifndef QSERL_3D_WORKSPACE_STATE_H_
#define QSERL_3D_WORKSPACE_STATE_H_

#include "qserl/exports.h"

#pragma warning( push, 0 )	
#include <Eigen/Lgsm>
#pragma warning( pop )	

#include "qserl/rod3d/parameters.h"
#include "util/forward_class.h"

namespace qserl {
namespace rod3d {

DECLARE_CLASS( WorkspaceState );

DECLARE_CLASS( WorkspaceIntegratedState );

class QSERL_EXPORT WorkspaceState
{
public:

	/**
	\brief Destructor.
	*/
	virtual ~WorkspaceState();

	/**
	\brief Constructor.
	* Compute rod directly from its geometry and its static parameters.
	*/
	static WorkspaceStateShPtr create(const std::vector<Eigen::Displacementd>& i_nodes, const Eigen::Displacementd& i_basePosition, 
		const Parameters& i_rodParams);

	/**
	* \brief Returns a copy of itself.
	*/
	virtual WorkspaceStateShPtr clone() const;

	/**
	* \brief Returns the number of nodes of the rod.
	*/
	size_t numNodes() const;

	/**
	* \brief Returns a vector of rod nodes positions, if initialized, in <b>base</b> frame.
	* Else returns an empty vector.
	*/
	const std::vector<Eigen::Displacementd>& nodes() const;

	/**
	* \brief Accessor to rod base position (in world frame).
	*/
	const Eigen::Displacementd& base() const;

	/**
	* \brief Seter for the rod base position (in world frame).
	* ? deprecated ?
	*/
	void base(const Eigen::Displacementd& i_base);

	/**
	* \brief Returns a vector of rod nodes positions, if initialized, in <b>world</b> frame.
	* Else returns an empty vector.
	*/
	std::vector<Eigen::Vector3d> nodesAbsolutePositions() const;

  /**
	* \brief Returns a vector of rod nodes 6D positions, if initialized, in <b>world</b> frame.
	* Else returns an empty vector.
	*/
	std::vector<Eigen::Displacementd> nodesAbsolute6DPositions() const;

	/** \brief Const accessor to rod static paramaters. */
	const Parameters& staticParameters() const;

	size_t memUsage() const;

	/**
	* \brief Returns the total torsional rotation along the rod (from base to tip),
	* i.e. rotation arround X axis.
	* \pre resolution must be sufficient such two successive rotations does not exceed PI.
	*/
	double torsionalRotation() const;

	/** \note WorkspaceIntegratedState needs a way to access protected members of other 
	* instances of this class. 
	*/
	friend class WorkspaceIntegratedState;

protected:
	
	/**
	\brief Constructor
	*/
	WorkspaceState(const std::vector<Eigen::Displacementd>& i_nodes, const Eigen::Displacementd& i_basePosition, const Parameters& i_rodParams);

	size_t																									m_numNodes;		/**< Number of nodes N. */
	std::vector<Eigen::Displacementd>												m_nodes;			/**< Position of each node (size N), in base frame. */
	Eigen::Displacementd																		m_base;				/**< DLO base position, in world frame (absolute). */

	Parameters  																						m_rodParameters;
};

}	// namespace rod3d
}	// namespace qserl

#endif // QSERL_3D_WORKSPACE_STATE_H_
