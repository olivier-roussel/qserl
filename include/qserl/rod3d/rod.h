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

#ifndef QSERL_3D_ROD_H_
#define QSERL_3D_ROD_H_

#include "qserl/exports.h"

#include "qserl/rod3d/parameters.h"
#include "qserl/rod3d/types.h"
#include "qserl/rod3d/workspace_integrated_state.h"
#include "qserl/util/forward_class.h"


namespace qserl {
namespace rod3d {

DECLARE_CLASS(WorkspaceState);

DECLARE_CLASS(Rod);

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
  static RodShPtr
  create(const Parameters& parameters);

  /**
  * \brief Returns the radius of the DLO.
  * \deprecated Kept for compatibilty
  */
  double
  radius() const;

  /**
  * \brief Returns the length of the DLO.
  * \deprecated Kept for compatibilty
  */
  double
  length() const;

  /**
  * \brief Returns the set of parameters of the DLO.
  */
  const Parameters&
  parameters() const;

  /************************************************************************/
  /*										State accessors & modifiers                       */
  /************************************************************************/

  /**
  * \brief Accessor to rod state.
  */
  WorkspaceStateShPtr
  state() const;

  /**
  * \brief Accessor to rod integrated state if exists, or null pointer otherwise
  */
  WorkspaceIntegratedStateShPtr
  integratedState() const;

  /**
  * \brief Setter for the rod state.
  */
  void
  state(const WorkspaceStateShPtr& i_state);

  /**
  * \brief Compute rod state from its base wrench.
  * Rod base is independant from this as node positions are computed in local base frame.
  * The corresponding rod state will be updated only if the result of integration leads to
  * WorkspaceIntegratedState::IR_VALID (see enum WorkspaceIntegratedState::IntegrationResultT).
  * \return The corresponding integration result status (see enum WorkspaceIntegratedState::IntegrationResultT).
  *	Note that IR_OUT_OF_WRENCH_BOUNDS cannot be returned, as out of bounds detection for internal
  * rod wrenches is not implemented yet.
  */
  WorkspaceIntegratedState::IntegrationResultT
  integrateStateFromBaseWrench(const Wrench& i_wrench,
                               const Displacement& i_basePos,
                               const WorkspaceIntegratedState::IntegrationOptions& i_integrationOptions);

  /************************************************************************/
  /*														Static members														*/
  /************************************************************************/

  /**
  * \brief Returns true iif rod configuration is singular,
  * (i.e. (a0.y, a0.z, a0.ry, a0.rz) = (0, 0, 0, 0)
  */
  static bool
  isConfigurationSingular(const Wrench& i_cfgWrench);

  /**
  * \brief Returns corresponding stifness coefficients a an extensible DLO from elasticity parameters.
  * \warning Deprecated. Use Parameters::getStiffnessCoefficients() instead.
  */
  static Eigen::Matrix<double, 6, 1>
  getStiffnessCoefficients(const Parameters& i_parameters);

protected:

  /**
  \brief Constructor
  */
  Rod(const Parameters& parameters);

  /**
  \brief Init function
  \param i_weakPtr The weak pointer to the Rod.
  */
  bool
  init(const RodWkPtr& i_weakPtr);

private:
  RodWkPtr m_weakPtr;

  Parameters m_staticParameters;

  WorkspaceStateShPtr m_state;

};

}  // namespace rod3d
}  // namespace qserl

#endif // QSERL_3D_ROD_H_
