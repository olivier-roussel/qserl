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

#include <boost/test/unit_test.hpp>

//#include "qserl/rod2d/analytic_mu.h"
#include "qserl/rod2d/analytic_dqda.h"
//#include "qserl/rod2d/workspace_integrated_state.h"
//#include "qserl/rod2d/rod.h"
#include "util/timer.h"

/* ------------------------------------------------------------------------- */
/* AnalyticVsNumericTests  																						 */
/* ------------------------------------------------------------------------- */
BOOST_AUTO_TEST_SUITE(AnalyticVsNumericTests)

BOOST_AUTO_TEST_CASE(BenchmarkAnalyticDqDa)
{
  static const size_t numRuns = 1000000;

  // set base A-space bounds
	Eigen::Matrix<double, 3, 1> aSpaceUpperBounds, aSpaceLowerBounds;
	static const double maxTorque = 6.28;
	static const double maxForce = 100;

  aSpaceUpperBounds[0] = maxTorque;
  aSpaceLowerBounds[0] = -maxTorque;
	for (int k = 1 ; k < 3 ; ++k)
	{
		aSpaceUpperBounds[k] = maxForce;
		aSpaceLowerBounds[k] = -maxForce;
	}

  Eigen::Vector3d wrench;
  Eigen::Matrix3d dqda;
  qserl::rod2d::MotionConstantsDqDa motionConstants;
  int successfullMotionConstants = 0;
  int successfullDqDa = 0;

	util::TimePoint startBenchTime = util::getTimePoint();
  for (size_t idxRun = 0; idxRun < numRuns; ++idxRun)
  {
    for (int k = 0; k < 3; ++k)
    {
      wrench[k] = ((rand() % 10000) * ( aSpaceUpperBounds[k] - aSpaceLowerBounds[k]) ) / 1.e4 + 
        aSpaceUpperBounds[k];
    }

    if (qserl::rod2d::computeMotionConstantsDqDa(wrench, motionConstants))
    {
      ++successfullMotionConstants;
      static const double t = 1.;
      //if (qserl::rod2d::computeDqDaAtPositionT(t, motionConstants, dqda))
      //{
      //  ++successfullDqDa;
      //}
    }
  }
	double benchTimeMs = static_cast<double>(util::getElapsedTimeMsec(startBenchTime).count());
	BOOST_TEST_MESSAGE( "Num runs: " << numRuns << " / success Motion Constants:" << successfullMotionConstants 
    << " / success Jacobian:" << successfullDqDa );
	BOOST_TEST_MESSAGE( "Benchmarking total time: " << benchTimeMs << "ms for "
    << numRuns << " analytic jacobians dq / da computations" );
	double benchTimePerJacobianUs = benchTimeMs * 1.e3 / static_cast<double>(numRuns);
	BOOST_TEST_MESSAGE( "  Avg. computation time per jacobian = " << benchTimePerJacobianUs << "us" );
	
}

BOOST_AUTO_TEST_SUITE_END();
//#endif