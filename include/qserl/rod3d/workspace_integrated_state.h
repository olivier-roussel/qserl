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
#ifndef QSERL_3D_WORKSPACE_INTEGRATED_STATE_H_
#define QSERL_3D_WORKSPACE_INTEGRATED_STATE_H_

#include "qserl/exports.h"

#include <boost/array.hpp>
#include <boost/function.hpp>

#include "qserl/rod3d/workspace_state.h"
#include "qserl/rod3d/parameters.h"
#include "qserl/util/forward_class.h"

namespace qserl {
namespace rod3d {

DECLARE_CLASS( WorkspaceIntegratedState );

class QSERL_EXPORT WorkspaceIntegratedState : public WorkspaceState
{
public:
	
	/** Rod state types. */
	typedef boost::array<double, 16>	state_type;						/**< Type of 3D elastic rod states q(t) at position t, where q(t) is an element of the Lie Group SE(3). */
	typedef boost::array<double, 6>		costate_type;					/**< Type of 3D elastic rod co-states mu(t) at position t, where mu(t) is an element of the dual Lie algebra se(3)*. */
	typedef boost::array<double, 72>	jacobian_state_type;	/**< Type of 3D elastic rod co-state and state derviates M(t) ( resp. J(t) ) at position t, where:
																													     - M(t) is the 6x6 Jacobian matrix of the co-state mu(t) w.r.t. initial conditions (i.e. mu(0)) (first 36 elements),
																															 - J(t) is the 6x6 Jacobian matrix of the state q(t) w.r.t. initial conditions (i.e. mu(0)) (last 36 elements). */

  /**< \brief Descriptors of possible status result returned by the rod integration process.*/
	enum IntegrationResultT{
		IR_VALID = 0,													/**< The rod configuration is valid w.r.t to given criterias (stability, max wrench, etc...) */
		IR_SINGULAR,													/**< The rod configuration is singular, i.e. a[1] = a[2] = 0. */
		IR_UNSTABLE,													/**< The rod configuration is unstable. */
		IR_OUT_OF_WRENCH_BOUNDS,							/**< The rod configuration is out of maximum allowed wrench. */
		IR_NUMBER_OF_INTEGRATION_RESULTS
	};

	/**
	* \brief Destructor.
	*/
	virtual ~WorkspaceIntegratedState();

	/**
	* \brief Constructor.
	* Rod base is independant from this as node positions are computed in local base frame.
	*/
	static WorkspaceIntegratedStateShPtr create(const Eigen::Wrenchd& i_baseWrench, unsigned int i_nnodes, 
		const Eigen::Displacementd& i_basePosition, const Parameters& i_rodParams);

	/**
	* \brief Copy constructor.
	*/
	static WorkspaceIntegratedStateShPtr createCopy(const WorkspaceIntegratedStateConstShPtr& i_other);

	/**
	* \brief Returns a copy of itself.
	*/
	virtual WorkspaceStateShPtr clone() const;

	/**
	* \brief Compute rod state from its base wrench by integration.
	* \return The corresponding integration result status (see enum IntegrationResultT).
	*	Note that IR_OUT_OF_WRENCH_BOUNDS cannot be returned, as out of bounds detection for internal
	* rod wrenches is not implemented yet.
	*/
	IntegrationResultT integrate();

	/**
	*\brief Approximate nodes positions by linearization for a neighboring state of this.
	* The local neihbour state is given by this plus given da, the base wrench perturbation.
	* Output positions are in local base frame.
	* \pre this state is initialized.
	*/
	WorkspaceStateShPtr approximateLinearlyNeighbourState(const Eigen::Wrenchd& i_da, const Eigen::Displacementd& i_neighbBase) const;

	/** \brief Set whether the singular values of the linear speed nu part of the Jacobian
	* (the 3x6 block starting at (3,0) ) should be computed during integration or not.
	* Default is set to true.
	* \warning Singular values decomposition if slow (takes approximately as long as the whole
	* integration process itself. 
	*/
	//void computeJacobianNuSingularValues(bool i_compute);

	/** \brief Returns whether the singular values of the linear speed nu part of the Jacobian
	* (the 3x6 block starting at (3,0) ) should be computed during integration or not.
	* Default is set to true.
	*/
	//bool computeJacobianNuSingularValues() const;

	/**
	* \brief Returns true if rod configuration is in a stable quasi-static confiugration.
	* \pre Rod must be initialized.
	*/
	bool isStable() const;

	/**
	* \brief Returns the wrench at the rod base.
	* \note This is equivalent to access through mu()[0]
	* \deprecated Use wrench(size_t i_idxNode) instead.
	*/
	Eigen::Wrenchd baseWrench() const;

	/**
	* \brief Returns the wrench at the rod tip.
	* \note This is equivalent to access through mu()[N-1]
	* \deprecated Use wrench(size_t i_idxNode) instead.
	*/
	Eigen::Wrenchd tipWrench() const;

	/**
	* \brief Returns the wrench at the rod given node.
	*/
	Eigen::Wrenchd wrench(size_t i_idxNode) const;

	/**
	* \brief Const accessor to the wrenches (costate) of the rod for each node. 
	*   \warning Only accessible if the keepMuValues() has been set to true,
	*		with the exception of mu(0) (i.e. the base wrench) which is kept in any case.
	*/
	const std::vector<costate_type>& mu() const;

	/** \brief Returns the M matrix (i.e. dmu(t) / dmu(0) ).
	*   \warning Only accessible if the keepMMatrices integration option has been set to true.
	*/
	const Eigen::Matrix<double, 6, 6>& getMMatrix(size_t i_nodeIdx) const;

	/** \brief Returns the J matrix (i.e. dq(t) / dmu(0) ).
	*   \warning Only accessible if the keepJMatrices integration option has been set to true.
	*/
	const Eigen::Matrix<double, 6, 6>& getJMatrix(size_t i_nodeIdx) const;

	/**
	* \brief Returns the values for each node of the jacobian determinant.
	* \warning If the DLO is detected as unstable, the determinant will be 0 from
	* the instability point.
	* \warning Only accessible if the keepJdet() has been set to true.
	*/
	const std::vector<double>& J_det() const;

	/**
	* \brief Const accessor to Jacobian linear speed part nu singular values.
	* \warning Only accessible if the computeJacobianNuSingularValues() has been set to true.
	*/
	const Eigen::Vector3d& J_nu_sv(size_t i_nodeIdx) const;

	/**
	* \brief Returns the memory usage of this instance.
	*/
	size_t memUsage() const;

	/**
	* \brief Integration computation options.
	*/
	struct IntegrationOptions
	{
		/** 
		* Constructor. 
		* Initialize to default values.
		*/
		IntegrationOptions();

		bool									computeJ_nu_sv;			/**< True if linear speed nu part of Jacobian matrix singular values should be computed. */
		bool									stop_if_unstable;		/**< True if integration process should be stop if configuration is detected as not stable. */
		bool									keepMuValues;
		bool									keepJdet;
		bool									keepMMatrices;
		bool									keepJMatrices;
	};

	/**
	* \brief Set integration options.
	* Must be done beofre integrate() call
	*/
	void integrationOptions(const IntegrationOptions& i_integrationOptions);

	/**
	* \brief Accessor to integration options.
	*/
	const IntegrationOptions& integrationOptions() const;

protected:
	
	/**
	\brief Constructor
	*/
	WorkspaceIntegratedState(unsigned int i_nnodes, const Eigen::Displacementd& i_basePosition, 
		const Parameters& i_rodParams);

	/**
	\brief Init function
	*/
	bool init(const Eigen::Wrenchd& i_wrench);

	/** \brief Integrates rod state from given base wrench.. 
      Numerical integration is done through a 4-th order Runge-Kutta with constant step. */
	IntegrationResultT integrateFromBaseWrenchRK4(const Eigen::Wrenchd& i_wrench);


	bool																										m_isInitialized;/**< True if the state has been integrated.*/
	bool																										m_isStable;		/**< True if DLO state is stable. */
	std::vector<costate_type>																m_mu;					/**< Wrenches at each nodes (size N). */
	std::vector<Eigen::Matrix<double, 6, 6> >								m_M;
	std::vector<Eigen::Matrix<double, 6, 6> >								m_J;

	std::vector<double>																			m_J_det;
	std::vector<Eigen::Vector3d>														m_J_nu_sv;			/**< Singular values of the linear speed nu part of the Jacobian matrix. */

	IntegrationOptions																			m_integrationOptions;
};

}	// namespace rod3d
}	// namespace qserl

#endif // QSERL_3D_WORKSPACE_INTEGRATED_STATE_H_
