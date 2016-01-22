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

#include "qserl/rod2d/analytic_mu.h"
#include "qserl/rod2d/workspace_integrated_state.h"
#include "qserl/rod2d/rod.h"
#include "qserl/util/timer.h"
#include "util/lie_algebra_utils.h"

/* ------------------------------------------------------------------------- */
/* Analytic vs. Numeric Tests for mu values              									   */
/* ------------------------------------------------------------------------- */
void compareAnalyticAndNumericMu(const Eigen::Vector3d& i_wrench, double i_errorTolerance,
  const qserl::rod2d::Parameters& i_rodParameters)
{
  // Integrate numerically mu(t)
  // ----------------------------------

	static const qserl::rod2d::Displacement2D identityDisp = { { 0., 0., 0. } };
  // set integration options
  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
  integrationOptions.stop_if_unstable = false;
	integrationOptions.keepMuValues = true;
	integrationOptions.keepJdet = true;
	integrationOptions.keepMMatrices = false;
	integrationOptions.keepJMatrices = true;

  qserl::rod2d::Wrench2D wrench2D;
  wrench2D[0] = i_wrench[1];
  wrench2D[1] = i_wrench[2];
  wrench2D[2] = i_wrench[0];

  qserl::rod2d::WorkspaceIntegratedStateShPtr rodState = qserl::rod2d::WorkspaceIntegratedState::create(wrench2D,
    identityDisp, i_rodParameters);
		rodState->integrationOptions(integrationOptions);
  qserl::rod2d::WorkspaceIntegratedState::IntegrationResultT integrationStatus = rodState->integrate();
  BOOST_CHECK( integrationStatus == qserl::rod2d::WorkspaceIntegratedState::IR_VALID );	

  // Compute analytically mu(t)
  // ----------------------------------
  qserl::rod2d::MotionConstantsMu motionConstants;
  bool motionConstantsSuccess = qserl::rod2d::computeMotionConstantsMu(i_wrench, motionConstants);

  BOOST_CHECK( motionConstantsSuccess );	

  if (motionConstantsSuccess)
  {
    const size_t numNodes = i_rodParameters.numberOfNodes();
    double maxError = 0.;
    for (size_t idxNode = 0; idxNode < numNodes; ++idxNode)
    {
      const double t =  static_cast<double>(idxNode) / static_cast<double>(numNodes - 1);
      Eigen::Vector3d mu_al;
      bool muSuccess = qserl::rod2d::computeMuAtPositionT(t, motionConstants, mu_al);

      BOOST_CHECK( muSuccess );	

      // Compare analytic vs. numerically integrated mu(t)
      // ----------------------------------

      const qserl::rod2d::Wrench2D wrench2D_num = rodState->wrench(idxNode);

      for (size_t j = 0; j < 3; ++j)
      {
        size_t jNum = j == 0 ? 2 : j-1; 
        const double err_mu_j = abs(mu_al[j] - wrench2D_num[jNum]);
        maxError = std::max(maxError, err_mu_j);
        BOOST_CHECK_SMALL( err_mu_j, i_errorTolerance );
      }
    }	

    BOOST_TEST_MESSAGE( "  Max error = " << maxError << " (must be less than " << i_errorTolerance << ")" );
  }
}

BOOST_AUTO_TEST_SUITE(AnalyticVsNumericTests_Mu)

BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Mu_CaseI_set1)
{
  Eigen::Vector3d wrench(2.3777, -49.6303, -9.8917);
  static const double errorTolerance = 1.e-8; // tolerance on the error of dq(i) / da(j) between
                                              // analytical and numerically integrated expressions
  qserl::rod2d::Parameters rodParameters;
	rodParameters.radius = 1.;
	rodParameters.length = 1.;
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
	rodParameters.delta_t = 1.e-3;            // integration resolution will impact on divergence with 
                                            // analytical forms

  compareAnalyticAndNumericMu(wrench, errorTolerance, rodParameters);
}

BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Mu_CaseI_set2)
{
  Eigen::Vector3d wrench(-1.2339, -21.8067, -12.0168);
  static const double errorTolerance = 1.e-8; // tolerance on the error of dq(i) / da(j) between
                                              // analytical and numerically integrated expressions
  qserl::rod2d::Parameters rodParameters;
	rodParameters.radius = 1.;
	rodParameters.length = 1.;
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
	rodParameters.delta_t = 1.e-3;              // integration resolution will impact on divergence with 
                                              // analytical forms 

  compareAnalyticAndNumericMu(wrench, errorTolerance, rodParameters);
}

// singularity with a3 = 0 _ not handled TODO
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Mu_CaseI_singular_a3_1)
//{
//  Eigen::Vector3d wrench(0., 4., 8.);
//  static const double errorTolerance = 1.e-8; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-3;
//  compareAnalyticAndNumericMu(wrench, errorTolerance, rodParameters);
//}

// singularity with a4 = 0
BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Mu_CaseI_singular_a4_1)
{
  Eigen::Vector3d wrench(2., 0., 6.);
  // XXX high deviation on dq1 / da here from numerical solutions 
  static const double errorTolerance = 1.e-8; // tolerance on the error of dq(i) / da(j) between
                                              // analytical and numerically integrated expressions
  qserl::rod2d::Parameters rodParameters;
	rodParameters.radius = 1.;
	rodParameters.length = 1.;
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
	rodParameters.delta_t = 1.e-3;              // integration resolution will impact on divergence with 
                                              // analytical forms

  compareAnalyticAndNumericMu(wrench, errorTolerance, rodParameters);
}

// singularity with a5 = 0 
BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Mu_CaseI_singular_a5_1)
{
  Eigen::Vector3d wrench(-1., -8., 0.);
  static const double errorTolerance = 1.e-8; // tolerance on the error of dq(i) / da(j) between
                                              // analytical and numerically integrated expressions
  qserl::rod2d::Parameters rodParameters;
	rodParameters.radius = 1.;
	rodParameters.length = 1.;
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
	rodParameters.delta_t = 1.e-3;
  compareAnalyticAndNumericMu(wrench, errorTolerance, rodParameters);
}

// singularity with a3 = a4 = 0 _ not handled TODO
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Mu_CaseI_singular_a3_a4_1)
//{
//  Eigen::Vector3d wrench(0., 0., -8.);
//  static const double errorTolerance = 1.e-8; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-3;
//  compareAnalyticAndNumericMu(wrench, errorTolerance, rodParameters);
//}

// singularity with a3 = a5 = 0 _ SINGULAR case (abnormal) not handled TODO
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Mu_CaseI_singular_a3_a5_1)
//{
//  Eigen::Vector3d wrench(0., 5., 0.);
//  static const double errorTolerance = 1.e-8; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-3;
//  compareAnalyticAndNumericMu(wrench, errorTolerance, rodParameters);
//}

BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Mu_CaseI_close_singular_a3_1)
{
  Eigen::Vector3d wrench(-1.e-4, 1.5, 17.);
  static const double errorTolerance = 1.e-8; // tolerance on the error of dq(i) / da(j) between
                                              // analytical and numerically integrated expressions
  qserl::rod2d::Parameters rodParameters;
	rodParameters.radius = 1.;
	rodParameters.length = 1.;
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
	rodParameters.delta_t = 1.e-3;

  compareAnalyticAndNumericMu(wrench, errorTolerance, rodParameters);
}

BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Mu_CaseII_set1)
{
  Eigen::Vector3d wrench(-4.1337, 87.7116, 18.0966);
  static const double errorTolerance = 1.e-6; 

  qserl::rod2d::Parameters rodParameters;
	rodParameters.radius = 1.;
	rodParameters.length = 1.;
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
	rodParameters.delta_t = 1.e-3;

  compareAnalyticAndNumericMu(wrench, errorTolerance, rodParameters);
}

BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Mu_CaseII_set2)
{
  Eigen::Vector3d wrench(4.1748, 23.4780, 4.0259);
  static const double errorTolerance = 1.e-8; 

  qserl::rod2d::Parameters rodParameters;
  rodParameters.radius = 1.;
  rodParameters.length = 1.;
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
  rodParameters.delta_t = 1.e-3;

  compareAnalyticAndNumericMu(wrench, errorTolerance, rodParameters);
}

BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Mu_CaseII_set3)
{
  Eigen::Vector3d wrench(1., 4., 0.5);
  static const double errorTolerance = 1.e-8; 

  qserl::rod2d::Parameters rodParameters;
	rodParameters.radius = 1.;
	rodParameters.length = 1.;
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
	rodParameters.delta_t = 1.e-3;

  compareAnalyticAndNumericMu(wrench, errorTolerance, rodParameters);
}

// singularity with a4 = 0 
BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Mu_CaseII_singular_a4_1)
{
  Eigen::Vector3d wrench(3., 0., -4.);
  static const double errorTolerance = 1.e-8; // tolerance on the error of dq(i) / da(j) between
                                              // analytical and numerically integrated expressions
  qserl::rod2d::Parameters rodParameters;
  rodParameters.radius = 1.;
  rodParameters.length = 1.;
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
  rodParameters.delta_t = 1.e-3;

  compareAnalyticAndNumericMu(wrench, errorTolerance, rodParameters);
}

// singularity with a5 = 0 _ not handled TODO
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Mu_CaseII_singular_a5_1)
//{
//  Eigen::Vector3d wrench(1., 4., 0.);
//  static const double errorTolerance = 1.e-8; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-3;
//  compareAnalyticAndNumericMu(wrench, errorTolerance, rodParameters);
//}

// singularity with a4 = a5 = 0 
BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Mu_CaseII_singular_a4_a5_1)
{
  Eigen::Vector3d wrench(1., 0., 0.);
  static const double errorTolerance = 1.e-8; // tolerance on the error of dq(i) / da(j) between
                                              // analytical and numerically integrated expressions
  qserl::rod2d::Parameters rodParameters;
	rodParameters.radius = 1.;
	rodParameters.length = 1.;
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
	rodParameters.delta_t = 1.e-3;
  compareAnalyticAndNumericMu(wrench, errorTolerance, rodParameters);
}

BOOST_AUTO_TEST_SUITE_END();

/* ------------------------------------------------------------------------- */
/* Analytic benchmarking of mu(t)                        									   */
/* ------------------------------------------------------------------------- */
BOOST_AUTO_TEST_SUITE(Analytic_Mu_Benchmarking)

BOOST_AUTO_TEST_CASE(Analytic_Mu_Benchmarking_Full_RandomSet)
{
  static const size_t numRuns = 100000;

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
  Eigen::Vector3d mu;
  qserl::rod2d::MotionConstantsMu motionConstants;
  int successfullMotionConstants = 0;
  int successfullMu = 0;

	qserl::util::TimePoint startBenchTime = qserl::util::getTimePoint();
  for (size_t idxRun = 0; idxRun < numRuns; ++idxRun)
  {
    for (int k = 0; k < 3; ++k)
    {
      wrench[k] = ((rand() % 10000) * ( aSpaceUpperBounds[k] - aSpaceLowerBounds[k]) ) / 1.e4 + 
        aSpaceUpperBounds[k];
    }

    if (qserl::rod2d::computeMotionConstantsMu(wrench, motionConstants))
    {
      ++successfullMotionConstants;
      static const double t = 1.; // we compute at t = 1 
      if (qserl::rod2d::computeMuAtPositionT(t, motionConstants, mu))
      {
        ++successfullMu;
      }
    }
  }
	double benchTimeMs = static_cast<double>(qserl::util::getElapsedTimeMsec(startBenchTime).count());
	BOOST_TEST_MESSAGE( "Num runs: " << numRuns << " / success Motion Constants:" << successfullMotionConstants 
    << " / success mu:" << successfullMu );
	BOOST_TEST_MESSAGE( "Benchmarking total time: " << benchTimeMs << "ms for "
    << numRuns << " analytic mu computations" );
	double benchTimePerMuUs = benchTimeMs * 1.e3 / static_cast<double>(numRuns);
	BOOST_TEST_MESSAGE( "  Avg. computation time per mu = " << benchTimePerMuUs << "us" );
	
}

BOOST_AUTO_TEST_SUITE_END();
