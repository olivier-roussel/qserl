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

#include "qserl/rod2d/workspace_integrated_state.h"
#include "util/timer.h"

/* ------------------------------------------------------------------------- */
/* SingularConfigurations2DTests  																						 */
/* ------------------------------------------------------------------------- */
BOOST_AUTO_TEST_SUITE(SingularConfigurations2DTests)

BOOST_AUTO_TEST_CASE(SingularConfigurations2DTest_1)
{
	qserl::rod2d::Parameters rodDefaultParameters;

	// 1. ckeck null base wrench 
	static const qserl::rod2d::Displacement2D identityDisp = { { 0. } };
	qserl::rod2d::Wrench2D singularWrench1 = { { 0. } };
	qserl::rod2d::WorkspaceIntegratedStateShPtr rodIntegratedState1 = qserl::rod2d::WorkspaceIntegratedState::create(singularWrench1,
		/*rodDefaultParameters.numNodes, */identityDisp, rodDefaultParameters);
	BOOST_CHECK( rodIntegratedState1 );	
	// singular 
	BOOST_CHECK( !rodIntegratedState1->integrate() );	
}

BOOST_AUTO_TEST_CASE(SingularConfigurations2DTest_2)
{
	qserl::rod2d::Parameters rodDefaultParameters;

	// 2. ckeck signular base wrench (a[1] = a[2] = 0)
	static const qserl::rod2d::Displacement2D identityDisp = { { 0. } };
	qserl::rod2d::Wrench2D singularWrench2 = { { 0. } };
	singularWrench2[0] = -1.5;
	qserl::rod2d::WorkspaceIntegratedStateShPtr rodIntegratedState2 = qserl::rod2d::WorkspaceIntegratedState::create(singularWrench2,
		/*rodDefaultParameters.numNodes, */identityDisp, rodDefaultParameters);
	BOOST_CHECK( rodIntegratedState2  );	
	// singular 
	BOOST_CHECK( !rodIntegratedState2->integrate()  );	
}

BOOST_AUTO_TEST_SUITE_END();

/* ------------------------------------------------------------------------- */
/* InextensibleRodStability2DTests																						 */
/* ------------------------------------------------------------------------- */
BOOST_AUTO_TEST_SUITE(InextensibleRodStability2DTests)

BOOST_AUTO_TEST_CASE(InextensibleRodStability2DTest_stable1)
{
	qserl::rod2d::Parameters rodParameters;
	// set appropriate elasticity parameters
	rodParameters.radius = 0.01;
	rodParameters.length = 1.;
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
	rodParameters.delta_t = 0.01;

	// set integration options
	qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
	integrationOptions.stop_if_unstable = false;
	integrationOptions.keepMuValues = true;
	integrationOptions.keepJdet = true;
	integrationOptions.keepMMatrices = true;
	integrationOptions.keepJMatrices = true;

	// stable configuration
	static const qserl::rod2d::Displacement2D identityDisp = { { 0. } };
	qserl::rod2d::Wrench2D stableConf1;
	stableConf1[0] = 0.;
	stableConf1[1] = 0.;
	stableConf1[2] = 1.;
	qserl::rod2d::WorkspaceIntegratedStateShPtr rodStableState1 = qserl::rod2d::WorkspaceIntegratedState::create(stableConf1,
		/*rodParameters.numNodes, */identityDisp, rodParameters);
	BOOST_CHECK( rodStableState1 );	
	rodStableState1->integrationOptions(integrationOptions);
	// not singular 
	BOOST_CHECK( rodStableState1->integrate() );	
	// stable
	BOOST_CHECK( rodStableState1->isStable() );

	// check co-state mu terminal values (sould be constant so equal to base wrench)
	qserl::rod2d::Wrench2D mu_last = rodStableState1->wrench(rodStableState1->numNodes()-1);
	BOOST_CHECK_CLOSE( mu_last[0], stableConf1[0], 1.e-6 );
	BOOST_CHECK_CLOSE( mu_last[1], stableConf1[1], 1.e-6 );
	BOOST_CHECK_CLOSE( mu_last[2], stableConf1[2], 1.e-6 );

	// check state q terminal values
	qserl::rod2d::Displacement2D q_last = rodStableState1->nodes().back();
	BOOST_CHECK_CLOSE( q_last[0], sin(1.), 1.e-3 );
	BOOST_CHECK_CLOSE( q_last[1], 1. - cos(1.), 1.e-3 );
	BOOST_CHECK_CLOSE( q_last[2], 1. , 1.e-6 );
}

BOOST_AUTO_TEST_CASE(InextensibleRodStability2DTest_iterativeIntegration_stable1)
{
	qserl::rod2d::Parameters rodParameters;
	// set appropriate elasticity parameters
	rodParameters.radius = 0.01;
	rodParameters.length = 1.;
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
	rodParameters.delta_t = 0.01;

	// set integration options
	qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
	integrationOptions.stop_if_unstable = false;
	integrationOptions.keepMuValues = true;
	integrationOptions.keepJdet = true;
	integrationOptions.keepMMatrices = true;
	integrationOptions.keepJMatrices = true;

	// stable configuration
	static const qserl::rod2d::Displacement2D identityDisp = { { 0. } };
	qserl::rod2d::Wrench2D stableConf1;
	stableConf1[0] = 0.;
	stableConf1[1] = 0.;
	stableConf1[2] = 1.;
	qserl::rod2d::WorkspaceIntegratedStateShPtr rodStableState1 = qserl::rod2d::WorkspaceIntegratedState::create(stableConf1,
		identityDisp, rodParameters);
	BOOST_CHECK( rodStableState1 );	
	rodStableState1->integrationOptions(integrationOptions);
	// not singular 
	double tinv = 0.;
	static const qserl::rod2d::Wrench2D maxWrench = { { std::numeric_limits<double>::max() } };
	qserl::rod2d::WorkspaceIntegratedState::IntegrationResultT integrationResult = rodStableState1->integrateWhileValid(maxWrench, tinv);
	BOOST_CHECK( integrationResult != qserl::rod2d::WorkspaceIntegratedState::IR_SINGULAR );	
	// stable
	BOOST_CHECK( rodStableState1->isStable() );
	BOOST_CHECK( integrationResult != qserl::rod2d::WorkspaceIntegratedState::IR_UNSTABLE );
	// as configuraiton is stable and infinite upper wrench bounds, invalid point is not found if we stop at this max number of nodes,
	// so should return < 0 tinv value.
	BOOST_CHECK( tinv < 0. );	


	// check co-state mu terminal values (sould be constant so equal to base wrench)
	qserl::rod2d::Wrench2D mu_last = rodStableState1->wrench(rodStableState1->numNodes()-1);
	BOOST_CHECK_CLOSE( mu_last[0], stableConf1[0], 1.e-6 );
	BOOST_CHECK_CLOSE( mu_last[1], stableConf1[1], 1.e-6 );
	BOOST_CHECK_CLOSE( mu_last[2], stableConf1[2], 1.e-6 );

	// check state q terminal values
	qserl::rod2d::Displacement2D q_last = rodStableState1->nodes().back();
	BOOST_CHECK_CLOSE( q_last[0], sin(1.), 1.e-3 );
	BOOST_CHECK_CLOSE( q_last[1], 1. - cos(1.), 1.e-3 );
	BOOST_CHECK_CLOSE( q_last[2], 1. , 1.e-6 );
}

BOOST_AUTO_TEST_SUITE_END();


/* ------------------------------------------------------------------------- */
/* Inextensible2DBencnhmarks																										 */
/* ------------------------------------------------------------------------- */
//#ifndef _DEBUG
BOOST_AUTO_TEST_SUITE(Inextensible2DBencnhmarks)

BOOST_AUTO_TEST_CASE(Inextensible2DBencnhmark_1)
{
	qserl::rod2d::Parameters rodParameters;
	// set appropriate elasticity parameters
	rodParameters.radius = 0.01;
	rodParameters.length = 1.;
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
	rodParameters.delta_t = 0.005;
	
	// set base A-space bounds
	Eigen::Matrix<double, 6, 1> aSpaceUpperBounds, aSpaceLowerBounds;
	static const double maxTorque = 6.28;
	static const double maxForce = 100;
	for (int k = 0 ; k < 3 ; ++k)
	{
		aSpaceUpperBounds[k] = maxTorque;
		aSpaceLowerBounds[k] = -maxTorque;
	}
	for (int k = 3 ; k < 6 ; ++k)
	{
		aSpaceUpperBounds[k] = maxForce;
		aSpaceLowerBounds[k] = -maxForce;
	}

	// set integration options
	qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
	integrationOptions.stop_if_unstable = false;
	integrationOptions.keepMuValues = true;
	integrationOptions.keepJdet = true;
	integrationOptions.keepMMatrices = false;
	integrationOptions.keepJMatrices = true;

	static const qserl::rod2d::Displacement2D identityDisp = { { 0. } };
	util::TimePoint startBenchTime = util::getTimePoint();
	qserl::rod2d::Wrench2D wrench;
	static const int numSamples = 500;
	int validSamples = 0;
	srand(1);
	for (int i = 0 ; i < numSamples ; ++i)
	{
		wrench[0] = ((rand() % 10000) * ( aSpaceUpperBounds[3] - aSpaceLowerBounds[3]) ) / 1.e4 + aSpaceUpperBounds[3];
		wrench[1] = ((rand() % 10000) * ( aSpaceUpperBounds[4] - aSpaceLowerBounds[4]) ) / 1.e4 + aSpaceUpperBounds[4];
		wrench[2] = ((rand() % 10000) * ( aSpaceUpperBounds[0] - aSpaceLowerBounds[0]) ) / 1.e4 + aSpaceUpperBounds[0];
		bool isStable = false;
		qserl::rod2d::WorkspaceIntegratedStateShPtr rodState = qserl::rod2d::WorkspaceIntegratedState::create(wrench,
			/*rodParameters.numNodes, */identityDisp, rodParameters);
		rodState->integrationOptions(integrationOptions);
		BOOST_CHECK( rodState );	
		// not singular (zero volume, so should never happen by random sampling
		bool isNotSingular = rodState->integrate();
		BOOST_CHECK( isNotSingular );	
		if (isNotSingular && rodState->isStable())
			++validSamples;
	}
	const int numNodes = rodParameters.numberOfNodes();
	double benchTimeMs = util::getElapsedTimeMsec(startBenchTime).count();
	BOOST_TEST_MESSAGE( "Benchmarking inextensible 2D rods total time: " << benchTimeMs << "ms for " << numSamples << " integrated rods" );
	double benchTimePerRodUs = benchTimeMs * 1.e3 / static_cast<double>(numSamples);
	BOOST_TEST_MESSAGE( "  Avg. integration time per rod = " << benchTimePerRodUs << "us for " << numNodes << " nodes" );
	BOOST_TEST_MESSAGE( "  Avg. integration time per rod node = " << benchTimePerRodUs / static_cast<double>(numNodes) << "us" );
	double stabilityRatio = static_cast<double>(validSamples) / static_cast<double>(numSamples);
	// stability ratio should be arround 99%
	static const double minStabilityRatio = 0.9;
	BOOST_TEST_MESSAGE( "  Stability ratio = " << stabilityRatio << " (should be above " << minStabilityRatio << ")" );
	BOOST_CHECK( stabilityRatio > minStabilityRatio );	

}

BOOST_AUTO_TEST_CASE(Inextensible2DBencnhmark_iterativeIntegration_1)
{
	qserl::rod2d::Parameters rodParameters;
	// set appropriate elasticity parameters
	rodParameters.radius = 0.01;
	rodParameters.length = 1.;
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
	rodParameters.delta_t = 0.005;
	
	// set A-space bounds
	Eigen::Matrix<double, 6, 1> aSpaceUpperBounds, aSpaceLowerBounds;
	static const double maxTorque = 6.28;
	static const double maxForce = 100;
	for (int k = 0 ; k < 3 ; ++k)
	{
		aSpaceUpperBounds[k] = maxTorque;
		aSpaceLowerBounds[k] = -maxTorque;
	}
	for (int k = 3 ; k < 6 ; ++k)
	{
		aSpaceUpperBounds[k] = maxForce;
		aSpaceLowerBounds[k] = -maxForce;
	}

	// set integration options
	qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
	integrationOptions.stop_if_unstable = false;	/**< Ignored in this case. */
	integrationOptions.keepMuValues = true;
	integrationOptions.keepJdet = true;
	integrationOptions.keepMMatrices = false;
	integrationOptions.keepJMatrices = true;

	static const qserl::rod2d::Displacement2D identityDisp = { { 0. } };
	util::TimePoint startBenchTime = util::getTimePoint();
	qserl::rod2d::Wrench2D wrench;
	static const int numSamples = 500;
	int validSamples = 0;
	srand(1);
	for (int i = 0 ; i < numSamples ; ++i)
	{
		wrench[0] = ((rand() % 10000) * ( aSpaceUpperBounds[3] - aSpaceLowerBounds[3]) ) / 1.e4 + aSpaceUpperBounds[3];
		wrench[1] = ((rand() % 10000) * ( aSpaceUpperBounds[4] - aSpaceLowerBounds[4]) ) / 1.e4 + aSpaceUpperBounds[4];
		wrench[2] = ((rand() % 10000) * ( aSpaceUpperBounds[0] - aSpaceLowerBounds[0]) ) / 1.e4 + aSpaceUpperBounds[0];
		bool isStable = false;
		qserl::rod2d::WorkspaceIntegratedStateShPtr rodState = qserl::rod2d::WorkspaceIntegratedState::create(wrench,
			/*rodParameters.numNodes, */identityDisp, rodParameters);
		rodState->integrationOptions(integrationOptions);
		BOOST_CHECK( rodState );	
		// not singular (zero volume, so should never happen by random sampling
		double tconj = 0.;
		static const qserl::rod2d::Wrench2D maxWrench = { { std::numeric_limits<double>::max() } };
		qserl::rod2d::WorkspaceIntegratedState::IntegrationResultT integrationResult = rodState->integrateWhileValid(maxWrench, tconj);
		BOOST_CHECK( integrationResult != qserl::rod2d::WorkspaceIntegratedState::IR_SINGULAR );	// zero volume singular region, so should never happen
		BOOST_CHECK( integrationResult != qserl::rod2d::WorkspaceIntegratedState::IR_OUT_OF_WRENCH_BOUNDS );	// infinite bounds but finite integration time
		BOOST_CHECK( rodState->isStable() );		// should be always stable as we keep last stable configuration
		if (integrationResult != qserl::rod2d::WorkspaceIntegratedState::IR_SINGULAR && 
			integrationResult != qserl::rod2d::WorkspaceIntegratedState::IR_UNSTABLE)
			++validSamples;
	}
	const int numNodes = rodParameters.numberOfNodes();
	double benchTimeMs = util::getElapsedTimeMsec(startBenchTime).count();
	BOOST_TEST_MESSAGE( "Benchmarking inextensible 2D rods (iterative mode) total time: " << benchTimeMs << "ms for " << numSamples << " integrated rods" );
	double benchTimePerRodUs = benchTimeMs * 1.e3 / static_cast<double>(numSamples);
	BOOST_TEST_MESSAGE( "  Avg. integration time per rod = " << benchTimePerRodUs << "us for " << numNodes << " nodes" );
	BOOST_TEST_MESSAGE( "  Avg. integration time per rod node = " << benchTimePerRodUs / static_cast<double>(numNodes) << "us" );
	double stabilityRatio = static_cast<double>(validSamples) / static_cast<double>(numSamples);
	// stability ratio should be arround 99%
	static const double minStabilityRatio = 0.9;
	BOOST_TEST_MESSAGE( "  Stability ratio = " << stabilityRatio << " (should be above " << minStabilityRatio << ")" );
	BOOST_CHECK( stabilityRatio > minStabilityRatio );	
}

BOOST_AUTO_TEST_SUITE_END();
//#endif