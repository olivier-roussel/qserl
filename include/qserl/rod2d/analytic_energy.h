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


#ifndef QSERL_2D_ANALYTIC_ENERGY_H_
#define QSERL_2D_ANALYTIC_ENERGY_H_

#include "qserl/exports.h"

#include "qserl/rod2d/analytic_q.h"

////#pragma warning( push, 0 )
#include <Eigen/Lgsm>
////#pragma warning( pop )

namespace qserl {
namespace rod2d {


/**
* \brief Compute the total elastic energy of the rod.
* \param[in]  i_motionConstants Constants of motion for the rod
* \pre a(i) (i=3..5) values corresponding to given motion constants must respect the unhandled following cases:
*   - if Case I (includes Case III) i.e. lambda4 >= 0, then a3 != 0 and a5 != 0
*/
QSERL_EXPORT bool
computeTotalElasticEnergy(const MotionConstantsQ& i_motionConstants,
                          double& o_energy);

}  // namespace rod2d
}  // namespace qserl

#endif // QSERL_2D_ANALYTIC_ENERGY_H_