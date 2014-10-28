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

#include "qserl/rod2d/rod.h"

#include <boost/math/constants/constants.hpp>

#include "qserl/rod2d/workspace_integrated_state.h"
#include "util/utils.h"

namespace qserl {
namespace rod2d {

/************************************************************************/
/*													Constructor																	*/
/************************************************************************/
Rod::Rod(const Parameters& i_parameters) :
	m_state(),
	m_staticParameters(i_parameters)
{
}

/************************************************************************/
/*												 Destructor																		*/
/************************************************************************/
Rod::~Rod()
{
}

/************************************************************************/
/*														create																		*/
/************************************************************************/
RodShPtr Rod::create(const Parameters& i_parameters)
{
	RodShPtr shPtr(new Rod(i_parameters));

	if (!shPtr->init(shPtr))
		shPtr.reset();

	return shPtr;
}

/************************************************************************/
/*														init																			*/
/************************************************************************/
bool Rod::init(const RodWkPtr& i_weakPtr)
{
	bool success = true;

	if (success)
		m_weakPtr = i_weakPtr;

	return success;
}

/************************************************************************/
/*												integrateStateFromBaseWrench									*/
/************************************************************************/
bool Rod::integrateStateFromBaseWrench(const Wrench2D& i_wrench, unsigned int i_nnodes, 
	const Displacement2D& i_basePos, const WorkspaceIntegratedState::IntegrationOptions& i_integrationOptions)
{
	WorkspaceIntegratedStateShPtr intState = WorkspaceIntegratedState::create(i_wrench, i_nnodes, i_basePos, m_staticParameters);
	intState->integrationOptions(i_integrationOptions);
	bool success = intState->integrate();
	if (success)
		m_state = intState;
	return success;
}

/************************************************************************/
/*												radius																				*/
/************************************************************************/
double Rod::radius() const
{
	return m_staticParameters.radius;
}

/************************************************************************/
/*												length																				*/
/************************************************************************/
double Rod::length() const
{
	return m_staticParameters.length;
}

/************************************************************************/
/*													parameters																	*/
/************************************************************************/
const Parameters& Rod::parameters() const
{
	return m_staticParameters;
}

/************************************************************************/
/*													state																				*/
/************************************************************************/
WorkspaceStateShPtr Rod::state() const
{
	return m_state;
}

/************************************************************************/
/*													integratedState															*/
/************************************************************************/
WorkspaceIntegratedStateShPtr Rod::integratedState() const
{
	WorkspaceIntegratedStateShPtr integratedState = boost::dynamic_pointer_cast<WorkspaceIntegratedState>(m_state);
	return integratedState;
}

/************************************************************************/
/*											isConfigurationSingular													*/
/************************************************************************/
bool Rod::isConfigurationSingular(const Wrench2D& i_cfgWrench)
{
	static const double kWrenchTolerance = 1.e-6;

		if (fabs(i_cfgWrench[1]) > kWrenchTolerance ||
			fabs(i_cfgWrench[2]) > kWrenchTolerance)
			return false;
		else
			return true;
}

/************************************************************************/
/*												getStiffnessCoefficients											*/
/************************************************************************/
double Rod::getStiffnessCoefficients(const Parameters& i_parameters)
{
	return 1.;
}

}	// namespace rod2d
}	// namespace qserl
