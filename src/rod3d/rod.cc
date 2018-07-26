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

#include "qserl/rod3d/rod.h"

#include <boost/math/constants/constants.hpp>

#include "qserl/rod3d/workspace_integrated_state.h"
#include "util/utils.h"

namespace qserl {
namespace rod3d {

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
WorkspaceIntegratedState::IntegrationResultT Rod::integrateStateFromBaseWrench(const Eigen::Wrenchd& i_wrench, 
	const Eigen::Displacementd& i_basePos, const WorkspaceIntegratedState::IntegrationOptions& i_integrationOptions)
{
	WorkspaceIntegratedStateShPtr intState = WorkspaceIntegratedState::create(i_wrench, m_staticParameters.numNodes, 
    i_basePos, m_staticParameters);
	intState->integrationOptions(i_integrationOptions);
	WorkspaceIntegratedState::IntegrationResultT success = intState->integrate();
	if (success == WorkspaceIntegratedState::IR_VALID)
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
	WorkspaceIntegratedStateShPtr integratedState = std::dynamic_pointer_cast<WorkspaceIntegratedState>(m_state);
	return integratedState;
}

/************************************************************************/
/*											isConfigurationSingular													*/
/************************************************************************/
bool Rod::isConfigurationSingular(const Eigen::Wrenchd& i_cfgWrench)
{
	static const double kWrenchTolerance = 1.e-6;

		if (fabs(i_cfgWrench[1]) > kWrenchTolerance ||
			fabs(i_cfgWrench[2]) > kWrenchTolerance || 
			fabs(i_cfgWrench[4]) > kWrenchTolerance ||
			fabs(i_cfgWrench[5]) > kWrenchTolerance)
			return false;
		else
			return true;
}

/************************************************************************/
/*												getStiffnessCoefficients											*/
/************************************************************************/
Eigen::Matrix<double, 6, 1> Rod::getStiffnessCoefficients(const Parameters& i_parameters)
{
  return i_parameters.stiffnessCoefficients;
}

}	// namespace rod3d
}	// namespace qserl
