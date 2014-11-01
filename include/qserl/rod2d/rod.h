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

#ifndef QSERL_2D_ROD_H_
#define QSERL_2D_ROD_H_

#include "qserl/exports.h"

#include <boost/thread.hpp>
#include <unordered_map>

#pragma warning( push, 0 )	
#include <Eigen/Lgsm>
#pragma warning( pop )	

#include "qserl/rod2d/parameters.h"
#include "qserl/rod2d/types.h"
#include "qserl/rod2d/workspace_integrated_state.h"
#include "util/forward_class.h"


namespace qserl {
namespace rod2d {

DECLARE_CLASS( WorkspaceState );

DECLARE_CLASS( Rod );

class QSERL_EXPORT Rod
{
public:

	/**
	* \brief Destructor.
	*/
	virtual ~Rod();

	/**
	* \brief Constructor.
	*/
	static RodShPtr create(const Parameters& parameters);

	/**
	* \brief Returns the radius of the DLO.
	* \deprecated Kept for compatibilty
	*/
	double radius() const;

	/**
	* \brief Returns the length of the DLO.
	* \deprecated Kept for compatibilty
	*/
	double length() const;

	/**
	* \brief Returns the set of parameters of the DLO.
	*/
	const Parameters& parameters() const;

	/************************************************************************/
	/*										State accessors & modifiers                       */
	/************************************************************************/

	/**
	* \brief Accessor to rod state.
	*/
	WorkspaceStateShPtr state() const;

	/**
	* \brief Accessor to rod integrated state if exists, or null pointer otherwise
	*/
	WorkspaceIntegratedStateShPtr integratedState() const;


	/**
	* \brief Compute rod state from its base wrench.
	* Rod base is independant from this as node positions are computed in local base frame.
	* Returns false if given configuration was singular.
	*/
	bool integrateStateFromBaseWrench(const Wrench2D& i_wrench, /*unsigned int i_nnodes, */
		const Displacement2D& i_basePos, const WorkspaceIntegratedState::IntegrationOptions& i_integrationOptions);

	/************************************************************************/
	/*														Static members														*/
	/************************************************************************/

	/**
	* \brief Returns true iif rod configuration is singular, 
	* (i.e. (a0[1], a0[2]) = (0, 0)
	*/
	static bool isConfigurationSingular(const Wrench2D& i_cfgWrench);
	
	/**
	* \brief Returns corresponding stifness coefficients a an extensible DLO from elasticity parameters. 
	*/
	static double getStiffnessCoefficients(const Parameters& i_parameters);

protected:
	
	/**
	\brief Constructor
	*/
	Rod(const Parameters& parameters);

	/**
	\brief Init function
	\param i_weakPtr The weak pointer to the Rod.
	*/
	bool init(const RodWkPtr& i_weakPtr);

private:
	RodWkPtr								m_weakPtr;

	Parameters							m_staticParameters;

	WorkspaceStateShPtr			m_state;		

};

}	// namespace rod2d
}	// namespace qserl

#endif // QSERL_2D_ROD_H_
