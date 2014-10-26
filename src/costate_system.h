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

#ifndef QSERL_QUASI_STATIC_ROD_COSTATE_SYSTEM_H_
#define QSERL_QUASI_STATIC_ROD_COSTATE_SYSTEM_H_

#include "qserl/exports.h"

#include <boost/function.hpp>

#pragma warning( push, 0 )	
#include <Eigen/Lgsm>
#pragma warning( pop )	

#include "qserl/workspace_integrated_state.h"

namespace qserl {

class QSERL_EXPORT CostateSystem
{
public:
	typedef WorkspaceIntegratedState::costate_type state_type;

	static const state_type kDefaultState;

	/**
	* Constructors, destructors
	*/
	CostateSystem(const Eigen::Matrix<double,6,1>& i_inv_stiffness, double i_length, 
		Parameters::RodModelT i_rodModel);
  virtual ~CostateSystem();

	void operator() (const state_type& i_mu, state_type& o_dmudt, double i_t);

private:

	Eigen::Matrix<double,6,1>																				m_inv_c;		/**< Inverse stiffness coefficients */
	double																													m_length;
	Parameters::RodModelT																						m_rodModel;

	boost::function<void(const state_type&, state_type&, double)>		m_evaluationCallback;
	
	/** 
	* Derivative evaluation at time t for the inextensible (RM_INEXTENSIBLE) rod model.
	*/
	void evaluateInextensible(const state_type& i_mu, state_type& o_dmudt, double i_t);
	
	/** 
	* Derivative evaluation at time t for the inextensible (RM_EXTENSIBLE_SHEARABLE) rod model.
	*/
	void evaluateExtensibleShearable(const state_type& i_mu, state_type& o_dmudt, double i_t);
};

}	// namespace qserl

#endif // QSERL_QUASI_STATIC_ROD_COSTATE_SYSTEM_H_
