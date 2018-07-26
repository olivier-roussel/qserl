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

/** Helper for computation of mu values. 
* Just for experimentation, should not be used as long term.
*/

#include "qserl/rod2d/analytic_energy.h"

#include <boost/math/special_functions/ellint_2.hpp>
#include "util/jacobi_elliptic.h"

#include "util/utils.h"

namespace qserl {
namespace rod2d {

bool computeTotalElasticEnergy(const MotionConstantsQ& i_mc, double& o_energy)
{
	const double gamma_1 = i_mc.r * (1. + i_mc.tau);

	double am_gamma_1;
	double cn_gamma_1_dummy, dn_gamma_1_dummy;
	boost::math::jacobi_elliptic(i_mc.k, gamma_1, &cn_gamma_1_dummy, &dn_gamma_1_dummy, &am_gamma_1);
	const double E_am_gamma_1 = boost::math::ellint_2(i_mc.k, am_gamma_1);

	if (i_mc.lambda[3] >= 0.)
	{
		// Case I: lambda_4 > 0 (includes case III : lambda_4 = 0 )
		o_energy = i_mc.alpha[2] * 0.5 / (i_mc.m * i_mc.r) * ( E_am_gamma_1 - i_mc.E_am_gamma_0 - i_mc.r * (1 - i_mc.m) );
	}else{
		// Case II: lambda_4 < 0 
		o_energy = i_mc.alpha[2] * 0.5 / (i_mc.r) * ( E_am_gamma_1 - i_mc.E_am_gamma_0 );
	}
	return true;
}
 
}	// namespace rod2d
}	// namespace qserl
