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

#include <boost/test/unit_test.hpp>

#include "qserl/rod2d/analytic_q.h"
#include "qserl/rod2d/workspace_integrated_state.h"
#include "qserl/rod2d/rod.h"
#include "qserl/util/timer.h"
#include "util/lie_algebra_utils.h"

/* ------------------------------------------------------------------------- */
/* Analytic vs. Numeric Tests for q values              									   */
/* ------------------------------------------------------------------------- */
void
compareAnalyticAndNumericQ(const Eigen::Vector3d& i_wrench,
                           double i_errorTolerance,
                           const qserl::rod2d::Parameters& i_rodParameters,
                           const qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions& i_integrationOptions)
{
  // Integrate numerically mu(t)
  // ----------------------------------

  static const qserl::rod2d::Displacement2D identityDisp = qserl::rod2d::Displacement2D::Zero();

  qserl::rod2d::Wrench2D wrench2D;
  wrench2D[0] = i_wrench[1];
  wrench2D[1] = i_wrench[2];
  wrench2D[2] = i_wrench[0];

  qserl::rod2d::WorkspaceIntegratedStateShPtr rodState = qserl::rod2d::WorkspaceIntegratedState::create(wrench2D,
                                                                                                        identityDisp,
                                                                                                        i_rodParameters);
  rodState->integrationOptions(i_integrationOptions);
  qserl::rod2d::WorkspaceIntegratedState::IntegrationResultT integrationStatus = rodState->integrate();
  BOOST_CHECK(integrationStatus == qserl::rod2d::WorkspaceIntegratedState::IR_VALID);

  // Compute analytically q(t)
  // ----------------------------------
  qserl::rod2d::MotionConstantsQ motionConstants;
  bool motionConstantsSuccess = qserl::rod2d::computeMotionConstantsQ(i_wrench, motionConstants);

  BOOST_CHECK(motionConstantsSuccess);

  double maxError = 0.;
  if(motionConstantsSuccess)
  {
    if(i_integrationOptions.integrator == qserl::rod2d::WorkspaceIntegratedState::IN_RK4)
    {
      const size_t numNodes = i_rodParameters.numberOfNodes();
      for(size_t idxNode = 0; idxNode < numNodes; ++idxNode)
      {
        const double t = static_cast<double>(idxNode) / static_cast<double>(numNodes - 1);
        Eigen::Vector3d q_dot_al, q_al;;
        bool qSuccess = qserl::rod2d::computeQAtPositionT(t, i_wrench, motionConstants, q_dot_al, q_al);

        BOOST_CHECK(qSuccess);

        // Compare analytic vs. numerically integrated q(t)
        // TODO also check q_dot(t)
        // ----------------------------------

        qserl::rod2d::Displacement2D q_disp2D_num = rodState->nodes()[idxNode];

        for(size_t j = 0; j < 3; ++j)
        {
          size_t jNum = j == 0 ? 2 : j - 1;
          // if testing angle theta, analytic solutions in the case II can give the total amount
          // of twist (i.e. integration of theta_dot) instead of using normalized theta to [-pi;pi]
          // as used in numerical integration thanks to the use of atan2
          // Note that in the case I we have to use atan2 in analytic form of theta, see notes for more details.
          // So to compare theta between analytic and numerical form, we normalize it to [-pi;pi]
          // using the (costly but simple) formula: theta_norm = atan2(sin(theta), cos(theta))
          if(j == 0)
          {
            q_al[j] = atan2(sin(q_al[j]), cos(q_al[j]));
            q_disp2D_num[jNum] = atan2(sin(q_disp2D_num[jNum]), cos(q_disp2D_num[jNum]));
          }
          double err_q_j = abs(q_al[j] - q_disp2D_num[jNum]);
          if(j == 0)
          {
            err_q_j = atan2(sin(err_q_j), cos(err_q_j));
          }

          maxError = std::max(maxError, err_q_j);
          BOOST_CHECK_SMALL(err_q_j, i_errorTolerance);
        }
      }
    }
    //else if (i_integrationOptions.integrator == qserl::rod2d::WorkspaceIntegratedState::IN_RK45)
    //{
    //  const double t =  1.;
    //  Eigen::Vector3d q_dot_al, q_al;;
    //  bool qSuccess = qserl::rod2d::computeQAtPositionT(t, i_wrench, motionConstants, q_dot_al, q_al);

    //  BOOST_CHECK( qSuccess );	

    //  // Compare analytic vs. numerically integrated q(t)
    //  // TODO also check q_dot(t)
    //  // ----------------------------------

    //  BOOST_CHECK( rodState->nodes().size() == 2 );	

    //  qserl::rod2d::Displacement2D q_disp2D_num = rodState->nodes()[1];

    //  for (size_t j = 0; j < 3; ++j)
    //  {
    //    size_t jNum = j == 0 ? 2 : j-1; 
    //    // if testing angle theta, analytic solutions in the case II can give the total amount
    //    // of twist (i.e. integration of theta_dot) instead of using normalized theta to [-pi;pi]
    //    // as used in numerical integration thanks to the use of atan2
    //    // Note that in the case I we have to use atan2 in analytic form of theta, see notes for more details.
    //    // So to compare theta between analytic and numerical form, we normalize it to [-pi;pi]
    //    // using the (costly but simple) formula: theta_norm = atan2(sin(theta), cos(theta))
    //    if (j == 0)
    //    {
    //      q_al[j] = atan2(sin(q_al[j]), cos(q_al[j]));
    //      q_disp2D_num[jNum] = atan2(sin(q_disp2D_num[jNum]), cos(q_disp2D_num[jNum]));
    //    }
    //    double err_q_j = abs(q_al[j] - q_disp2D_num[jNum]);
    //    if (j == 0)
    //      err_q_j = atan2(sin(err_q_j), cos(err_q_j));

    //    maxError = std::max(maxError, err_q_j);
    //    BOOST_CHECK_SMALL( err_q_j, i_errorTolerance );
    //  }
    //}else{
    //  BOOST_ERROR( "unhandled integrator type" );
    //}

    BOOST_TEST_MESSAGE("  Max error = " << maxError << " (must be less than " << i_errorTolerance << ")");
  }
}

BOOST_AUTO_TEST_SUITE(AnalyticVsNumericTests_Q_RK4)

BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseI_set1_RK4)
{
  Eigen::Vector3d wrench(2.3777, -49.6303, -9.8917);
  static const double errorTolerance = 1.e-3; // tolerance on the error of dq(i) / da(j) between
  // analytical and numerically integrated expressions
  qserl::rod2d::Parameters rodParameters;
  rodParameters.radius = 1.;
  rodParameters.length = 1.;
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
  rodParameters.delta_t = 1.e-4;            // integration resolution will impact on divergence with
  // analytical forms

  // set integration options
  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
  integrationOptions.stop_if_unstable = false;
  integrationOptions.keepMuValues = true;
  integrationOptions.keepJdet = true;
  integrationOptions.keepMMatrices = false;
  integrationOptions.keepJMatrices = true;
  integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK4;

  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters, integrationOptions);
}

BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseI_set2_RK4)
{
  Eigen::Vector3d wrench(-1.2339, -21.8067, -12.0168);
  static const double errorTolerance = 1.e-3; // tolerance on the error of dq(i) / da(j) between
  // analytical and numerically integrated expressions
  qserl::rod2d::Parameters rodParameters;
  rodParameters.radius = 1.;
  rodParameters.length = 1.;
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
  rodParameters.delta_t = 1.e-4;              // integration resolution will impact on divergence with
  // analytical forms

  // set integration options
  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
  integrationOptions.stop_if_unstable = false;
  integrationOptions.keepMuValues = true;
  integrationOptions.keepJdet = true;
  integrationOptions.keepMMatrices = false;
  integrationOptions.keepJMatrices = true;
  integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK4;

  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters, integrationOptions);
}

BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseI_set3_RK4)
{
  Eigen::Vector3d wrench(2.5, 8., -16.);
  static const double errorTolerance = 1.e-3; // tolerance on the error of dq(i) / da(j) between
  // analytical and numerically integrated expressions
  qserl::rod2d::Parameters rodParameters;
  rodParameters.radius = 1.;
  rodParameters.length = 1.;
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
  rodParameters.delta_t = 1.e-4;              // integration resolution will impact on divergence with
  // analytical forms

  // set integration options
  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
  integrationOptions.stop_if_unstable = false;
  integrationOptions.keepMuValues = true;
  integrationOptions.keepJdet = true;
  integrationOptions.keepMMatrices = false;
  integrationOptions.keepJMatrices = true;
  integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK4;

  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters, integrationOptions);
}

// singularity with a3 = 0 _ not handled TODO
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseI_singular_a3_1_RK4)
//{
//  Eigen::Vector3d wrench(0., 4., 8.);
//  static const double errorTolerance = 1.e-2; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-3;

//// set integration options
//qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
//integrationOptions.stop_if_unstable = false;
//integrationOptions.keepMuValues = true;
//integrationOptions.keepJdet = true;
//integrationOptions.keepMMatrices = false;
//integrationOptions.keepJMatrices = true;
//integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK4;

//compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters, integrationOptions);
//}

// singularity with a4 = 0
BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseI_singular_a4_1_RK4)
{
  Eigen::Vector3d wrench(2., 0., 6.);
  // XXX high deviation on dq1 / da here from numerical solutions 
  static const double errorTolerance = 1.e-3; // tolerance on the error of dq(i) / da(j) between
  // analytical and numerically integrated expressions
  qserl::rod2d::Parameters rodParameters;
  rodParameters.radius = 1.;
  rodParameters.length = 1.;
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
  rodParameters.delta_t = 1.e-4;              // integration resolution will impact on divergence with
  // analytical forms

  // set integration options
  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
  integrationOptions.stop_if_unstable = false;
  integrationOptions.keepMuValues = true;
  integrationOptions.keepJdet = true;
  integrationOptions.keepMMatrices = false;
  integrationOptions.keepJMatrices = true;
  integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK4;

  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters, integrationOptions);
}

// singularity with a5 = 0 
BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseI_singular_a5_1_RK4)
{
  Eigen::Vector3d wrench(-1., -8., 0.);
  static const double errorTolerance = 1.e-3; // tolerance on the error of dq(i) / da(j) between
  // analytical and numerically integrated expressions
  qserl::rod2d::Parameters rodParameters;
  rodParameters.radius = 1.;
  rodParameters.length = 1.;
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
  rodParameters.delta_t = 1.e-4;

  // set integration options
  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
  integrationOptions.stop_if_unstable = false;
  integrationOptions.keepMuValues = true;
  integrationOptions.keepJdet = true;
  integrationOptions.keepMMatrices = false;
  integrationOptions.keepJMatrices = true;
  integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK4;

  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters, integrationOptions);
}

// singularity with a3 = a4 = 0 _ not handled TODO
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseI_singular_a3_a4_1_RK4)
//{
//  Eigen::Vector3d wrench(0., 0., -8.);
//  static const double errorTolerance = 1.e-2; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-3;
//  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters);
//}

// singularity with a3 = a5 = 0 _ SINGULAR case (abnormal) not handled TODO
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseI_singular_a3_a5_1_RK4)
//{
//  Eigen::Vector3d wrench(0., 5., 0.);
//  static const double errorTolerance = 1.e-2; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-3;
//  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters);
//}

BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseI_close_singular_a3_1_RK4)
{
  Eigen::Vector3d wrench(-1.e-4, 1.5, 17.);
  static const double errorTolerance = 1.e-3; // tolerance on the error of dq(i) / da(j) between
  // analytical and numerically integrated expressions
  qserl::rod2d::Parameters rodParameters;
  rodParameters.radius = 1.;
  rodParameters.length = 1.;
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
  rodParameters.delta_t = 1.e-4;

  // set integration options
  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
  integrationOptions.stop_if_unstable = false;
  integrationOptions.keepMuValues = true;
  integrationOptions.keepJdet = true;
  integrationOptions.keepMMatrices = false;
  integrationOptions.keepJMatrices = true;
  integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK4;

  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters, integrationOptions);
}

BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseII_set1_RK4)
{
  Eigen::Vector3d wrench(-4.1337, 87.7116, 18.0966);
  static const double errorTolerance = 1.e-3;

  qserl::rod2d::Parameters rodParameters;
  rodParameters.radius = 1.;
  rodParameters.length = 1.;
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
  rodParameters.delta_t = 1.e-4;

  // set integration options
  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
  integrationOptions.stop_if_unstable = false;
  integrationOptions.keepMuValues = true;
  integrationOptions.keepJdet = true;
  integrationOptions.keepMMatrices = false;
  integrationOptions.keepJMatrices = true;
  integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK4;

  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters, integrationOptions);
}

BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseII_set2_RK4)
{
  Eigen::Vector3d wrench(4.1748, 23.4780, 4.0259);
  static const double errorTolerance = 1.e-3;

  qserl::rod2d::Parameters rodParameters;
  rodParameters.radius = 1.;
  rodParameters.length = 1.;
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
  rodParameters.delta_t = 1.e-4;

  // set integration options
  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
  integrationOptions.stop_if_unstable = false;
  integrationOptions.keepMuValues = true;
  integrationOptions.keepJdet = true;
  integrationOptions.keepMMatrices = false;
  integrationOptions.keepJMatrices = true;
  integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK4;

  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters, integrationOptions);
}

BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseII_set3_RK4)
{
  Eigen::Vector3d wrench(1., 4., 0.5);
  static const double errorTolerance = 1.e-3;

  qserl::rod2d::Parameters rodParameters;
  rodParameters.radius = 1.;
  rodParameters.length = 1.;
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
  rodParameters.delta_t = 1.e-4;

  // set integration options
  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
  integrationOptions.stop_if_unstable = false;
  integrationOptions.keepMuValues = true;
  integrationOptions.keepJdet = true;
  integrationOptions.keepMMatrices = false;
  integrationOptions.keepJMatrices = true;
  integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK4;

  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters, integrationOptions);
}

// singularity with a4 = 0 
BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseII_singular_a4_1_RK4)
{
  Eigen::Vector3d wrench(3., 0., -4.);
  static const double errorTolerance = 1.e-3; // tolerance on the error of dq(i) / da(j) between
  // analytical and numerically integrated expressions
  qserl::rod2d::Parameters rodParameters;
  rodParameters.radius = 1.;
  rodParameters.length = 1.;
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
  rodParameters.delta_t = 1.e-4;

  // set integration options
  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
  integrationOptions.stop_if_unstable = false;
  integrationOptions.keepMuValues = true;
  integrationOptions.keepJdet = true;
  integrationOptions.keepMMatrices = false;
  integrationOptions.keepJMatrices = true;
  integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK4;

  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters, integrationOptions);
}

// singularity with a5 = 0 _ not handled TODO
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseII_singular_a5_1_RK4)
//{
//  Eigen::Vector3d wrench(1., 4., 0.);
//  static const double errorTolerance = 1.e-2; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-3;
//  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters);
//}

// singularity with a4 = a5 = 0 _ not handled TODO
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseII_singular_a4_a5_1_RK4)
//{
//  Eigen::Vector3d wrench(1., 0., 0.);
//  static const double errorTolerance = 1.e-2; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-3;
//  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters);
//}

BOOST_AUTO_TEST_SUITE_END();

/* ------------------------------------------------------------------------- */
/* Analytic benchmarking of q(t)                        									   */
/* ------------------------------------------------------------------------- */
BOOST_AUTO_TEST_SUITE(Analytic_Q_Benchmarking)

BOOST_AUTO_TEST_CASE(Analytic_Q_Benchmarking_Full_RandomSet)
{
  static const size_t numRuns = 10000;

  // set base A-space bounds
  Eigen::Matrix<double, 3, 1> aSpaceUpperBounds, aSpaceLowerBounds;
  static const double maxTorque = 6.28;
  static const double maxForce = 100;

  aSpaceUpperBounds[0] = maxTorque;
  aSpaceLowerBounds[0] = -maxTorque;
  for(int k = 1; k < 3; ++k)
  {
    aSpaceUpperBounds[k] = maxForce;
    aSpaceLowerBounds[k] = -maxForce;
  }

  Eigen::Vector3d wrench;
  Eigen::Vector3d dq, q;
  qserl::rod2d::MotionConstantsQ motionConstants;
  int successfullMotionConstants = 0;
  int successfullQ = 0;

  qserl::util::TimePoint startBenchTime = qserl::util::getTimePoint();
  for(size_t idxRun = 0; idxRun < numRuns; ++idxRun)
  {
    for(int k = 0; k < 3; ++k)
    {
      wrench[k] = ((rand() % 10000) * (aSpaceUpperBounds[k] - aSpaceLowerBounds[k])) / 1.e4 +
                  aSpaceUpperBounds[k];
    }

    if(qserl::rod2d::computeMotionConstantsQ(wrench, motionConstants))
    {
      ++successfullMotionConstants;
      static const double t = 1.; // we compute at t = 1 
      if(qserl::rod2d::computeQAtPositionT(t, wrench, motionConstants, dq, q))
      {
        ++successfullQ;
      }
    }
  }
  double benchTimeMs = static_cast<double>(qserl::util::getElapsedTimeMsec(startBenchTime).count());
  BOOST_TEST_MESSAGE("Num runs: " << numRuns << " / success Motion Constants:" << successfullMotionConstants
                                  << " / success q:" << successfullQ);
  BOOST_TEST_MESSAGE("Benchmarking total time: " << benchTimeMs << "ms for "
                                                 << numRuns << " analytic (q_dot, q) computations");
  double benchTimePerQUs = benchTimeMs * 1.e3 / static_cast<double>(numRuns);
  BOOST_TEST_MESSAGE("  Avg. computation time per (q_dot, q) = " << benchTimePerQUs << "us");

}

BOOST_AUTO_TEST_SUITE_END();


/* ------------------------------------------------------------------------- */
/* Numeric benchmarking of q(1)                        									   */
/* ------------------------------------------------------------------------- */
BOOST_AUTO_TEST_SUITE(Numeric_Q1_Benchmarking)

BOOST_AUTO_TEST_CASE(Numeric_Q1_Benchmarking_Full_RandomSet)
{
  static const size_t numRuns = 10000;

  // set base A-space bounds
  Eigen::Matrix<double, 3, 1> aSpaceUpperBounds, aSpaceLowerBounds;
  static const double maxTorque = 6.28;
  static const double maxForce = 100;

  for(int k = 0; k < 2; ++k)
  {
    aSpaceUpperBounds[k] = maxForce;
    aSpaceLowerBounds[k] = -maxForce;
  }
  aSpaceUpperBounds[2] = maxTorque;
  aSpaceLowerBounds[2] = -maxTorque;

  qserl::rod2d::Wrench2D wrench_XYT;
  //qserl::rod2d::Displacement2D q1;
  int successfullQ = 0;

  qserl::rod2d::Parameters rodParameters;
  // set appropriate elasticity parameters
  rodParameters.radius = 0.01;
  rodParameters.length = 1.;
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
  rodParameters.delta_t = 1.e-3;  // with RK4: for a 1.e-3 max error on q(t) w.r.t analytic forms, should take at least 1.e-4 for delta_t

  // set integration options
  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
  integrationOptions.stop_if_unstable = false;
  integrationOptions.keepMuValues = false;
  integrationOptions.keepJdet = false;
  integrationOptions.keepMMatrices = false;
  integrationOptions.keepJMatrices = false;
  integrationOptions.computeJacobians = false;
  integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK4;

  static const qserl::rod2d::Displacement2D identityDisp = qserl::rod2d::Displacement2D::Zero();

  qserl::util::TimePoint startBenchTime = qserl::util::getTimePoint();
  for(size_t idxRun = 0; idxRun < numRuns; ++idxRun)
  {
    for(int k = 0; k < 3; ++k)
    {
      wrench_XYT[k] = ((rand() % 10000) * (aSpaceUpperBounds[k] - aSpaceLowerBounds[k])) / 1.e4 +
                      aSpaceUpperBounds[k];
    }

    qserl::rod2d::WorkspaceIntegratedStateShPtr rodState = qserl::rod2d::WorkspaceIntegratedState::create(wrench_XYT,
                                                                                                          identityDisp,
                                                                                                          rodParameters);
    BOOST_CHECK(rodState);
    rodState->integrationOptions(integrationOptions);
    // not singular (zero volume, so should never happen by random sampling
    qserl::rod2d::WorkspaceIntegratedState::IntegrationResultT integrationStatus = rodState->integrate();
    //BOOST_CHECK ( integrationStatus != qserl::rod2d::WorkspaceIntegratedState::IR_SINGULAR );
    if(integrationStatus == qserl::rod2d::WorkspaceIntegratedState::IR_VALID)
    {
      ++successfullQ;
    }
  }
  double benchTimeMs = static_cast<double>(qserl::util::getElapsedTimeMsec(startBenchTime).count());
  BOOST_TEST_MESSAGE("Num runs: " << numRuns << " / success q1:" << successfullQ);
  BOOST_TEST_MESSAGE("Benchmarking total time: " << benchTimeMs << "ms for "
                                                 << numRuns << " numeric (q) computations");
  double benchTimePerQUs = benchTimeMs * 1.e3 / static_cast<double>(numRuns);
  BOOST_TEST_MESSAGE("  Avg. computation time per (q) = " << benchTimePerQUs << "us");

}

BOOST_AUTO_TEST_SUITE_END();


/** --------------Testing with Runge Kutta 45 ---------------  */
//
//BOOST_AUTO_TEST_SUITE(AnalyticVsNumericTests_Q_RK45)
//
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseI_set1_RK45)
//{
//  Eigen::Vector3d wrench(2.3777, -49.6303, -9.8917);
//  static const double errorTolerance = 1.e-3; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-2;            // integration resolution will impact on divergence with 
//                                            // analytical forms
//
//  // set integration options
//  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
//  integrationOptions.stop_if_unstable = false;
//  integrationOptions.keepMuValues = true;
//  integrationOptions.keepJdet = true;
//  integrationOptions.keepMMatrices = false;
//  integrationOptions.keepJMatrices = true;
//  integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK45;
//
//  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters, integrationOptions);
//}
//
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseI_set2_RK45)
//{
//  Eigen::Vector3d wrench(-1.2339, -21.8067, -12.0168);
//  static const double errorTolerance = 1.e-3; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-2;              // integration resolution will impact on divergence with 
//                                              // analytical forms 
//
//  // set integration options
//  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
//  integrationOptions.stop_if_unstable = false;
//  integrationOptions.keepMuValues = true;
//  integrationOptions.keepJdet = true;
//  integrationOptions.keepMMatrices = false;
//  integrationOptions.keepJMatrices = true;
//  integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK45;
//
//  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters, integrationOptions);
//}
//
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseI_set3_RK45)
//{
//  Eigen::Vector3d wrench(2.5, 8., -16.);
//  static const double errorTolerance = 1.e-3; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-2;              // integration resolution will impact on divergence with 
//                                              // analytical forms 
//
//  // set integration options
//  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
//  integrationOptions.stop_if_unstable = false;
//  integrationOptions.keepMuValues = true;
//  integrationOptions.keepJdet = true;
//  integrationOptions.keepMMatrices = false;
//  integrationOptions.keepJMatrices = true;
//  integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK45;
//
//  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters, integrationOptions);
//}
//
//// singularity with a4 = 0
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseI_singular_a4_1_RK45)
//{
//  Eigen::Vector3d wrench(2., 0., 6.);
//  // XXX high deviation on dq1 / da here from numerical solutions 
//  static const double errorTolerance = 1.e-3; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-2;              // integration resolution will impact on divergence with 
//                                              // analytical forms
//
//  // set integration options
//  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
//  integrationOptions.stop_if_unstable = false;
//  integrationOptions.keepMuValues = true;
//  integrationOptions.keepJdet = true;
//  integrationOptions.keepMMatrices = false;
//  integrationOptions.keepJMatrices = true;
//  integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK45;
//
//  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters, integrationOptions);
//}
//
//// singularity with a5 = 0 
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseI_singular_a5_1_RK45)
//{
//  Eigen::Vector3d wrench(-1., -8., 0.);
//  static const double errorTolerance = 1.e-3; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-2;
//
//  // set integration options
//  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
//  integrationOptions.stop_if_unstable = false;
//  integrationOptions.keepMuValues = true;
//  integrationOptions.keepJdet = true;
//  integrationOptions.keepMMatrices = false;
//  integrationOptions.keepJMatrices = true;
//  integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK45;
//
//  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters, integrationOptions);
//}
//
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseI_close_singular_a3_1_RK45)
//{
//  Eigen::Vector3d wrench(-1.e-4, 1.5, 17.);
//  static const double errorTolerance = 1.e-3; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-2;
//
//  // set integration options
//  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
//  integrationOptions.stop_if_unstable = false;
//  integrationOptions.keepMuValues = true;
//  integrationOptions.keepJdet = true;
//  integrationOptions.keepMMatrices = false;
//  integrationOptions.keepJMatrices = true;
//  integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK45;
//
//  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters, integrationOptions);
//}
//
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseII_set1_RK45)
//{
//  Eigen::Vector3d wrench(-4.1337, 87.7116, 18.0966);
//  static const double errorTolerance = 1.e-3; 
//
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-2;
//
//  // set integration options
//  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
//  integrationOptions.stop_if_unstable = false;
//  integrationOptions.keepMuValues = true;
//  integrationOptions.keepJdet = true;
//  integrationOptions.keepMMatrices = false;
//  integrationOptions.keepJMatrices = true;
//  integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK45;
//
//  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters, integrationOptions);
//}
//
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseII_set2_RK45)
//{
//  Eigen::Vector3d wrench(4.1748, 23.4780, 4.0259);
//  static const double errorTolerance = 1.e-3; 
//
//  qserl::rod2d::Parameters rodParameters;
//  rodParameters.radius = 1.;
//  rodParameters.length = 1.;
//  rodParameters.integrationTime = 1.;
//  rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//  rodParameters.delta_t = 1.e-2;
//
//  // set integration options
//  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
//  integrationOptions.stop_if_unstable = false;
//  integrationOptions.keepMuValues = true;
//  integrationOptions.keepJdet = true;
//  integrationOptions.keepMMatrices = false;
//  integrationOptions.keepJMatrices = true;
//  integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK45;
//
//  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters, integrationOptions);
//}
//
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseII_set3_RK45)
//{
//  Eigen::Vector3d wrench(1., 4., 0.5);
//  static const double errorTolerance = 1.e-3; 
//
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-2;
//
//  // set integration options
//  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
//  integrationOptions.stop_if_unstable = false;
//  integrationOptions.keepMuValues = true;
//  integrationOptions.keepJdet = true;
//  integrationOptions.keepMMatrices = false;
//  integrationOptions.keepJMatrices = true;
//  integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK45;
//
//  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters, integrationOptions);
//}
//
//// singularity with a4 = 0 
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_Q_CaseII_singular_a4_1_RK45)
//{
//  Eigen::Vector3d wrench(3., 0., -4.);
//  static const double errorTolerance = 1.e-3; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//  rodParameters.radius = 1.;
//  rodParameters.length = 1.;
//  rodParameters.integrationTime = 1.;
//  rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//  rodParameters.delta_t = 1.e-2;
//
//  // set integration options
//  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
//  integrationOptions.stop_if_unstable = false;
//  integrationOptions.keepMuValues = true;
//  integrationOptions.keepJdet = true;
//  integrationOptions.keepMMatrices = false;
//  integrationOptions.keepJMatrices = true;
//  integrationOptions.integrator = qserl::rod2d::WorkspaceIntegratedState::IN_RK45;
//
//  compareAnalyticAndNumericQ(wrench, errorTolerance, rodParameters, integrationOptions);
//}
//
//BOOST_AUTO_TEST_SUITE_END();

//#endif
