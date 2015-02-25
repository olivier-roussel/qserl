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
#include "qserl/rod2d/workspace_integrated_state.h"
#include "qserl/rod2d/rod.h"
#include "util/lie_algebra_utils.h"
#include "util/timer.h"

/* ------------------------------------------------------------------------- */
/* Analytic vs. Numeric Tests for the Jacobian dq / da  									   */
/* ------------------------------------------------------------------------- */
void compareAnalyticAndNumericDqDa(const Eigen::Vector3d& i_wrench, double i_errorTolerance,
  const qserl::rod2d::Parameters& i_rodParameters)
{
  // Integrate numerically the jacobian dq / da
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

  // Compute analytically the jacobian dq / da
  // ----------------------------------
  qserl::rod2d::MotionConstantsDqDa motionConstants;
  bool motionConstantsSuccess = qserl::rod2d::computeMotionConstantsDqDa(i_wrench, motionConstants);

  BOOST_CHECK( motionConstantsSuccess );	

  if (motionConstantsSuccess)
  {
    const size_t numNodes = i_rodParameters.numberOfNodes();
    double maxError = 0.;
    for (size_t idxNode = 0; idxNode < numNodes; ++idxNode)
    {
      const double t =  static_cast<double>(idxNode) / static_cast<double>(numNodes - 1);
      Eigen::Matrix3d dqda_al;
      bool dqdaSuccess = qserl::rod2d::computeDqDaAtPositionT(t, motionConstants, dqda_al);

      BOOST_CHECK( dqdaSuccess );	

      // Compare analytic vs. numerically integrated jacobians dq / da
      // ----------------------------------

      // as the jacobian from numerical integration is mapped to some velocities of se(2)
      // in the body frame by TqLq^-1, we have to bring them back to Tq(SE2) by TqLq
      const Eigen::Matrix3d& J_num_Te_Body = rodState->getJMatrix(idxNode);

      // As elements of Tq(SE2) show the same skew-symetry than elements of Te(SE2),
      // we can represent the as a vector 3x1
      std::vector<Eigen::Matrix3d> J_num_Tq(3, Eigen::Matrix3d::Zero());
      std::vector<Eigen::Matrix3d> J_num_Te_Body_Skew_SE2(3, Eigen::Matrix3d::Zero());
      const qserl::rod2d::Displacement2D& q_disp2D = rodState->nodes()[idxNode];
      const Eigen::Displacementd q_disp3D = qserl::rod2d::toDisplacement3D(q_disp2D);
      Eigen::Matrix3d q_mat;
      qserl::rod2d::toHomogeneousMatrix(q_disp2D, q_mat);
      Eigen::Vector3d dtheta_da; // dtheta_da(j) for j=1..3
      for (size_t j = 0; j < 3; ++j)
      {
        qserl::util::hat_SE2(J_num_Te_Body.col(j), J_num_Te_Body_Skew_SE2[j]);
        J_num_Tq[j] = q_mat * J_num_Te_Body_Skew_SE2[j];
        const double dsin_theta_da = J_num_Tq[j](1,0);
        const double sin_theta = q_mat(1,0);
        const double dcos_theta_da = J_num_Tq[j](0,0);
        const double cos_theta = q_mat(0,0);
        dtheta_da[j] = dsin_theta_da * cos_theta - dcos_theta_da * sin_theta;
      }

      for (size_t j = 0; j < 3; ++j)
      {
        size_t jNum = j == 0 ? 2 : j-1; 
        // compare dq1 / da (i.e. dtheta / da)
        const double err_dq1 = abs(dqda_al(0,j) - dtheta_da[jNum]);
        maxError = std::max(maxError, err_dq1);
        BOOST_CHECK_SMALL( err_dq1, i_errorTolerance );
        // compare dq2 / da (i.e. dx / da)
        const double err_dq2 = abs(dqda_al(1,j) - J_num_Tq[jNum](0,2));
        maxError = std::max(maxError, err_dq2);
        BOOST_CHECK_SMALL( err_dq2, i_errorTolerance );
        // compare dq3 / da (i.e. dy / da)
        const double err_dq3 = abs(dqda_al(2,j) - J_num_Tq[jNum](1,2));
        maxError = std::max(maxError, err_dq3);
        BOOST_CHECK_SMALL( err_dq3, i_errorTolerance );
      }
    }	

    BOOST_TEST_MESSAGE( "  Max error = " << maxError << " (must be less than " << i_errorTolerance << ")" );
  }
}

BOOST_AUTO_TEST_SUITE(AnalyticVsNumericTests_DqDa)

BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_DqDa_CaseI_set1)
{
  Eigen::Vector3d wrench(2.3777, -49.6303, -9.8917);
  static const double errorTolerance = 1.e-4; // tolerance on the error of dq(i) / da(j) between
                                              // analytical and numerically integrated expressions
  qserl::rod2d::Parameters rodParameters;
	rodParameters.radius = 1.;
	rodParameters.length = 1.;
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
	rodParameters.delta_t = 1.e-4;            // integration resolution will impact on divergence with 
                                            // analytical forms

  compareAnalyticAndNumericDqDa(wrench, errorTolerance, rodParameters);
}

BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_DqDa_CaseI_set2)
{
  Eigen::Vector3d wrench(-1.2339, -21.8067, -12.0168);
  static const double errorTolerance = 1.e-4; // tolerance on the error of dq(i) / da(j) between
                                              // analytical and numerically integrated expressions
  qserl::rod2d::Parameters rodParameters;
	rodParameters.radius = 1.;
	rodParameters.length = 1.;
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
	rodParameters.delta_t = 1.e-4;              // integration resolution will impact on divergence with 
                                              // analytical forms 

  compareAnalyticAndNumericDqDa(wrench, errorTolerance, rodParameters);
}

// singularity with a3 = 0 _ not handled TODO
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_DqDa_CaseI_singular_a3_1)
//{
//  Eigen::Vector3d wrench(0., 4., 8.);
//  static const double errorTolerance = 1.e-3; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-4;
//  compareAnalyticAndNumericDqDa(wrench, errorTolerance, rodParameters);
//}

// singularity with a4 = 0
BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_DqDa_CaseI_singular_a4_1)
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

  compareAnalyticAndNumericDqDa(wrench, errorTolerance, rodParameters);
}

// singularity with a5 = 0 _ not handled TODO
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_DqDa_CaseI_singular_a5_1)
//{
//  Eigen::Vector3d wrench(-1., -8., 0.);
//  static const double errorTolerance = 1.e-3; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-4;
//  compareAnalyticAndNumericDqDa(wrench, errorTolerance, rodParameters);
//}

// singularity with a3 = a4 = 0 _ not handled TODO
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_DqDa_CaseI_singular_a3_a4_1)
//{
//  Eigen::Vector3d wrench(0., 0., -8.);
//  static const double errorTolerance = 1.e-3; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-4;
//  compareAnalyticAndNumericDqDa(wrench, errorTolerance, rodParameters);
//}

// singularity with a3 = a5 = 0 _ not handled TODO
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_DqDa_CaseI_singular_a3_a5_1)
//{
//  Eigen::Vector3d wrench(0., 5., 0.);
//  static const double errorTolerance = 1.e-3; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-4;
//  compareAnalyticAndNumericDqDa(wrench, errorTolerance, rodParameters);
//}

BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_DqDa_CaseI_close_singular_a3_1)
{
  Eigen::Vector3d wrench(-1.e-4, 1.5, 17.);
  // warning _ as a[0] (i.e. a3) is close to zero, we are getting closer to singularity 
  // and both solutions diverge _ so the error tolerance is here lowered to 0.01
  static const double errorTolerance = 1.e-4; // tolerance on the error of dq(i) / da(j) between
                                              // analytical and numerically integrated expressions
  qserl::rod2d::Parameters rodParameters;
	rodParameters.radius = 1.;
	rodParameters.length = 1.;
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
	rodParameters.delta_t = 1.e-4;

  compareAnalyticAndNumericDqDa(wrench, errorTolerance, rodParameters);
}

BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_DqDa_CaseII_set1)
{
  Eigen::Vector3d wrench(-4.1337, 87.7116, 18.0966);
  static const double errorTolerance = 1.e-3; 

  qserl::rod2d::Parameters rodParameters;
	rodParameters.radius = 1.;
	rodParameters.length = 1.;
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
	rodParameters.delta_t = 1.e-4;

  compareAnalyticAndNumericDqDa(wrench, errorTolerance, rodParameters);
}

BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_DqDa_CaseII_set2)
{
  Eigen::Vector3d wrench(4.1748, 23.4780, 4.0259);
  static const double errorTolerance = 1.e-3; 

  qserl::rod2d::Parameters rodParameters;
  rodParameters.radius = 1.;
  rodParameters.length = 1.;
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
  rodParameters.delta_t = 1.e-4;

  compareAnalyticAndNumericDqDa(wrench, errorTolerance, rodParameters);
}

BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_DqDa_CaseII_set3)
{
  Eigen::Vector3d wrench(1., 4., 0.5);
  static const double errorTolerance = 1.e-4; 

  qserl::rod2d::Parameters rodParameters;
	rodParameters.radius = 1.;
	rodParameters.length = 1.;
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
	rodParameters.delta_t = 1.e-4;

  compareAnalyticAndNumericDqDa(wrench, errorTolerance, rodParameters);
}

// singularity with a4 = 0 
BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_DqDa_CaseII_singular_a4_1)
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

  compareAnalyticAndNumericDqDa(wrench, errorTolerance, rodParameters);
}

// singularity with a5 = 0 _ not handled TODO
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_DqDa_CaseII_singular_a5_1)
//{
//  Eigen::Vector3d wrench(1., 4., 0.);
//  static const double errorTolerance = 1.e-3; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-4;
//  compareAnalyticAndNumericDqDa(wrench, errorTolerance, rodParameters);
//}

// singularity with a4 = a5 = 0 _ not handled TODO
//BOOST_AUTO_TEST_CASE(AnalyticVsNumericTest_DqDa_CaseII_singular_a4_a5_1)
//{
//  Eigen::Vector3d wrench(1., 0., 0.);
//  static const double errorTolerance = 1.e-3; // tolerance on the error of dq(i) / da(j) between
//                                              // analytical and numerically integrated expressions
//  qserl::rod2d::Parameters rodParameters;
//	rodParameters.radius = 1.;
//	rodParameters.length = 1.;
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.delta_t = 1.e-4;
//  compareAnalyticAndNumericDqDa(wrench, errorTolerance, rodParameters);
//}

BOOST_AUTO_TEST_SUITE_END();

/* ------------------------------------------------------------------------- */
/* Analytic benchmarking of the Jacobian dq / da         									   */
/* ------------------------------------------------------------------------- */
BOOST_AUTO_TEST_SUITE(Analytic_DqDa_Benchmarking)

BOOST_AUTO_TEST_CASE(Analytic_DqDa_Benchmarking_Full_RandomSet)
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
      if (qserl::rod2d::computeDqDaAtPositionT(t, motionConstants, dqda))
      {
        ++successfullDqDa;
      }
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