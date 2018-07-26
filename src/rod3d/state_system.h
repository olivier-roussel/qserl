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

#ifndef QSERL_3D_STATE_SYSTEM_H_
#define QSERL_3D_STATE_SYSTEM_H_

#include "qserl/exports.h"

#include <boost/function.hpp>

#pragma warning( push, 0 )	
#include <Eigen/Lgsm>
#pragma warning( pop )	

#include "qserl/rod3d/workspace_integrated_state.h"

namespace qserl {
namespace rod3d {

class QSERL_EXPORT StateSystem
{
public:
	typedef WorkspaceIntegratedState::state_type state_type;

	//static const state_type kDefaultState;

	/**
	* Constructors, destructors
	*/
	StateSystem(const Eigen::Matrix<double,6, 1>& i_inv_stiffness, double i_length,  
		double i_dt, const std::vector<WorkspaceIntegratedState::costate_type>& i_mu,
		Parameters::RodModelT i_rodModel);
  virtual ~StateSystem();

	void operator() (const state_type& i_q, state_type& o_dqdt, double i_t);

	/** Returns default state value. */
	static state_type defaultState();

private:

	Eigen::Matrix<double,6, 1>																				m_inv_c;			/**< Inverse stiffness coefficients */
	double																														m_dt;			
	const std::vector<WorkspaceIntegratedState::costate_type>&				m_mu;
	double																														m_length;
	Parameters::RodModelT																							m_rodModel;

	boost::function<void(const state_type&, state_type&, double)>			m_evaluationCallback;

	/** 
	* Derivative evaluation at time t for the inextensible (RM_INEXTENSIBLE) rod model.
	*/
	void evaluateInextensible(const state_type& i_q, state_type& o_dqdt, double i_t);
	
	/** 
	* Derivative evaluation at time t for the inextensible (RM_EXTENSIBLE_SHEARABLE) rod model.
	*/
	void evaluateExtensibleShearable(const state_type& i_q, state_type& o_dqdt, double i_t);
};

}	// namespace rod3d
}	// namespace qserl

#endif // QSERL_3D_STATE_SYSTEM_H_

