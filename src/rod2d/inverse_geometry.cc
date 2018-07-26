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

#include "qserl/rod2d/inverse_geometry.h"

#include "qserl/rod2d/analytic_dqda.h"
#include "qserl/util/timer.h"
#include "util/utils.h"

namespace qserl {
namespace rod2d {

bool inverseGeometry_Newton(const Eigen::Vector3d& i_q_des, int i_maxIter, double i_maxNormError, 
  const Eigen::Vector3d& i_a0, double i_alpha, Eigen::Vector3d& o_a)
{
  static const double t = 1.;
  static const double kJacobianDetNullTolerance = 1.e-9;

	util::TimePoint startSolveTime = util::getTimePoint();

  const double sqrdMaxNormError = util::sqr(i_maxNormError);
  bool isQ1Close = false;
  bool isUnstable = false;
  bool singularityFound = false;
  MotionConstantsDqDa mc_a_k;
  Eigen::Matrix3d dq_da_k;  // jacobian dq(t) / da for a_k
  Eigen::Matrix3d inv_dq_da_k;  // inverse of jacobian dq(t) / da for a_k
  Eigen::Vector3d q_a_k, q_dot_a_k;  // q(t) for a_k and q_dot(t) for a_k (unused)
  Eigen::Vector3d a_k = i_a0;
  
  int iter = 0;
  while (!isQ1Close && iter < i_maxIter && !singularityFound && !isUnstable)
  {
    // compute common motion constants to dq / da and q for a_k 
    singularityFound = !computeMotionConstantsDqDa(a_k, mc_a_k) ||
      !computeQAtPositionT(t, a_k, mc_a_k.qc, q_dot_a_k, q_a_k) ||
      !computeDqDaAtPositionT(t, mc_a_k, dq_da_k);
    if (!singularityFound)
    {
      const Eigen::Vector3d r_a_k = q_a_k - i_q_des;
      std::cout << "[DEBUG] [iter=" << iter << "] a_k = " << a_k << std::endl;
      std::cout << "[DEBUG] [iter=" << iter << "] q(a_k) = " << q_a_k << std::endl;
      std::cout << "[DEBUG] [iter=" << iter << "] ||r(a_k)|| = " << r_a_k.norm() << std::endl;
      if (r_a_k.squaredNorm() > sqrdMaxNormError)
      {
        // inverse the jacobian and check rank
        bool isJacobianInvertible;
        double det_J_k;
        dq_da_k.computeInverseAndDetWithCheck(inv_dq_da_k, det_J_k, isJacobianInvertible, kJacobianDetNullTolerance);
        if (isJacobianInvertible)
        {
          const Eigen::Vector3d p_k = -inv_dq_da_k * r_a_k * i_alpha;
          a_k = a_k + p_k;
        }else{
          std::cout << "[WARNING] inverseGeometry_Newton(): unstability point found at iter = " << iter << std::endl;
          isUnstable = true;
        }
      }else{
        // solution found
        isQ1Close = true;
      }
    }else{
      std::cout << "[WARNING] inverseGeometry_Newton(): singularity point found at iter = " << iter << std::endl;
    }
    ++iter;
  }
	double solveTimeUs = static_cast<double>(util::getElapsedTimeUsec(startSolveTime).count());

  const bool success = !singularityFound && !isUnstable && iter < i_maxIter;
  if (success)
  {
    o_a = a_k;
    std::cout << "[PROGRESS] inverseGeometry_Newton(): succefully solved after " << iter << " iterations" <<
    " _ took " << solveTimeUs << "us" << std::endl;
  }
  return success;
}


}	// namespace rod2d
}	// namespace qserl
