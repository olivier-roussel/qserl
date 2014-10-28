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

#ifndef QSERL_3D_JACOBIAN_SYSTEM_H_
#define QSERL_3D_JACOBIAN_SYSTEM_H_

#include "qserl/exports.h"

#include <boost/function.hpp>

#pragma warning( push, 0 )	
#include <Eigen/Lgsm>
#pragma warning( pop )	

#include "qserl/rod3d/workspace_integrated_state.h"

namespace qserl {
namespace rod3d {

class QSERL_EXPORT JacobianSystem
{
public:
	typedef WorkspaceIntegratedState::jacobian_state_type state_type;

	static const state_type kDefaultState;				/**< Deprecated */

	static const double kStabilityThreshold;			/** Threshold for Jacobian determinant. */
	static const double kStabilityTolerance;			/** Tolerance for which Jacobian determinant vanishes. */

	/**
	* Constructors, destructors
	*/
	JacobianSystem(const Eigen::Matrix<double, 6, 1>& i_inv_stiffness, double i_dt, 
		const std::vector<WorkspaceIntegratedState::costate_type>& i_mu, Parameters::RodModelT i_rodModel);
  virtual ~JacobianSystem();

	/**
	* Derivative evaluation at time t.
	*/
	void operator() (const state_type& MJ, state_type& dMJdt, double t);

private:

	Eigen::Matrix<double, 6, 1>																			m_inv_c;	/**< inverse stiffness coefficients */
	Eigen::Matrix<double, 6, 1>																			m_b;			/**<	Precomputed values from inverse stiffness coefficients, where:
																																									b(1) = inv_c(3) - inv_c(2) 
																																									b(2) = inv_c(1) - inv_c(3) 
																																									b(3) = inv_c(2) - inv_c(1) 
																																									b(4) = inv_c(6) - inv_c(5) 
																																									b(5) = inv_c(4) - inv_c(6) 
																																									b(6) = inv_c(5) - inv_c(4) 
																																									*/
	double																													m_dt;			
	const std::vector<WorkspaceIntegratedState::costate_type>&			m_mu;
	Parameters::RodModelT																						m_rodModel;

	boost::function<void(const state_type&, state_type&, double)>		m_evaluationCallback;

	/** 
	* Derivative evaluation at time t for the inextensible (RM_INEXTENSIBLE) rod model.
	*/
	void evaluateInextensible(const state_type& MJ, state_type& dMJdt, double t);
	
	/** 
	* Derivative evaluation at time t for the inextensible (RM_EXTENSIBLE_SHEARABLE) rod model.
	*/
	void evaluateExtensibleShearable(const state_type& MJ, state_type& dMJdt, double t);
};

}	// namespace rod3d

}	// namespace qserl

#endif // QSERL_3D_JACOBIAN_SYSTEM_H_
