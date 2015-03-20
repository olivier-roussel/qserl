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

#include "qserl/rod3d/workspace_integrated_state.h"
#include "util/timer.h"

/* ------------------------------------------------------------------------- */
/* SingularConfigurations3DTests  																						 */
/* ------------------------------------------------------------------------- */
BOOST_AUTO_TEST_SUITE(SingularConfigurations3DTests)

BOOST_AUTO_TEST_CASE(SingularConfigurations3DTest_1)
{
	qserl::rod3d::Parameters rodDefaultParameters;

	// 1. ckeck null base wrench 
	Eigen::Wrenchd singularWrench1(Eigen::Wrenchd::Zero());
	qserl::rod3d::WorkspaceIntegratedStateShPtr rodIntegratedState1 = qserl::rod3d::WorkspaceIntegratedState::create(singularWrench1,
		rodDefaultParameters.numNodes, Eigen::Displacementd::Identity(), rodDefaultParameters);
	BOOST_CHECK( rodIntegratedState1 );	
	// singular 
	BOOST_CHECK( !rodIntegratedState1->integrate() );	
}

BOOST_AUTO_TEST_CASE(SingularConfigurations3DTest_2)
{
	qserl::rod3d::Parameters rodDefaultParameters;

	// 2. ckeck signular base wrench (a[1] = a[2] = a[4] = a[5] = 0)
	Eigen::Wrenchd singularWrench2(Eigen::Wrenchd::Zero());
	singularWrench2[0] = -1.5;
	singularWrench2[3] = 2.;
	qserl::rod3d::WorkspaceIntegratedStateShPtr rodIntegratedState2 = qserl::rod3d::WorkspaceIntegratedState::create(singularWrench2,
		rodDefaultParameters.numNodes, Eigen::Displacementd::Identity(), rodDefaultParameters);
	BOOST_CHECK( rodIntegratedState2  );	
	// singular 
	BOOST_CHECK( !rodIntegratedState2->integrate()  );	
}

BOOST_AUTO_TEST_SUITE_END();

/* ------------------------------------------------------------------------- */
/* ExtensibleRodStability3DTests  																						 */
/* ------------------------------------------------------------------------- */
BOOST_AUTO_TEST_SUITE(ExtensibleRodStability3DTests)

BOOST_AUTO_TEST_CASE(ExtensibleRodStability3DTest_stable1)
{
	qserl::rod3d::Parameters rodParameters;
	// set appropriate elasticity parameters
	rodParameters.radius = 0.01;
	rodParameters.length = 1.;
	const double youngModulus = 15.4e6;  /** Default Young modulus of rubber: 15.4 MPa */
	const double shearModulus = 5.13e6;  /** Default Shear modulus of rubber: 5.13 MPa */
  rodParameters.setIsotropicStiffnessCoefficientsFromElasticityParameters(youngModulus, shearModulus);
	rodParameters.density = 1.1 * 10e3;   /** 1.10 kg/dm3 -> kg/m3, */ 
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod3d::Parameters::RM_EXTENSIBLE_SHEARABLE;
	rodParameters.numNodes = 100;

	// stable configuration
	Eigen::Wrenchd stableConf1;
	stableConf1[0] = -0.3967;
	stableConf1[1] = 0.2774;
	stableConf1[2] = 0.1067;
	stableConf1[3] = 0.54;
	stableConf1[4] = 1.501;
	stableConf1[5] = 0.2606;
	qserl::rod3d::WorkspaceIntegratedStateShPtr rodStableState1 = qserl::rod3d::WorkspaceIntegratedState::create(stableConf1,
		rodParameters.numNodes, Eigen::Displacementd::Identity(), rodParameters);
	BOOST_CHECK( rodStableState1 );	
	// not singular 
	BOOST_CHECK( rodStableState1->integrate() );	
	// stable
	BOOST_CHECK( rodStableState1->isStable() );
}

BOOST_AUTO_TEST_CASE(ExtensibleRodStability3DTest_stable2)
{
	qserl::rod3d::Parameters rodParameters;
	// set appropriate elasticity parameters
	rodParameters.radius = 0.01;
	rodParameters.length = 1.;
	const double youngModulus = 15.4e6;  /** Default Young modulus of rubber: 15.4 MPa */
	const double shearModulus = 5.13e6;  /** Default Shear modulus of rubber: 5.13 MPa */
  rodParameters.setIsotropicStiffnessCoefficientsFromElasticityParameters(youngModulus, shearModulus);
	rodParameters.density = 1.1 * 10e3;   /** 1.10 kg/dm3 -> kg/m3, */ 
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod3d::Parameters::RM_EXTENSIBLE_SHEARABLE;
	rodParameters.numNodes = 100;

	// stable configuration
	Eigen::Wrenchd stableConf2;
	stableConf2[0] = 0.5205;
	stableConf2[1] = 0.2989;
	stableConf2[2] = 0.0875;
	stableConf2[3] = 0.9518;
	stableConf2[4] = -0.8417;
	stableConf2[5] = -0.8075;
	qserl::rod3d::WorkspaceIntegratedStateShPtr rodStableState2 = qserl::rod3d::WorkspaceIntegratedState::create(stableConf2,
		rodParameters.numNodes, Eigen::Displacementd::Identity(), rodParameters);
	BOOST_CHECK( rodStableState2 );	
	// not singular 
	BOOST_CHECK( rodStableState2->integrate() );	
	// stable
	BOOST_CHECK( rodStableState2->isStable() );
}

BOOST_AUTO_TEST_CASE(ExtensibleRodStability3DTest_unstable1)
{
	qserl::rod3d::Parameters rodParameters;
	// set appropriate elasticity parameters
	rodParameters.radius = 0.01;
	rodParameters.length = 1.;
	const double youngModulus = 15.4e6;  /** Default Young modulus of rubber: 15.4 MPa */
	const double shearModulus = 5.13e6;  /** Default Shear modulus of rubber: 5.13 MPa */
  rodParameters.setIsotropicStiffnessCoefficientsFromElasticityParameters(youngModulus, shearModulus);
	rodParameters.density = 1.1 * 10e3;   /** 1.10 kg/dm3 -> kg/m3, */ 
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod3d::Parameters::RM_EXTENSIBLE_SHEARABLE;
	rodParameters.numNodes = 100;

	// unstable configuration
	Eigen::Wrenchd unstableConf1;
	unstableConf1[0] = -0.5885;
	unstableConf1[1] = -0.7467;
	unstableConf1[2] = 0.4277;
	unstableConf1[3] = -0.121;
	unstableConf1[4] = 0.0508;
	unstableConf1[5] = 0.9760;
	qserl::rod3d::WorkspaceIntegratedStateShPtr rodUnstableState1 = qserl::rod3d::WorkspaceIntegratedState::create(unstableConf1,
		rodParameters.numNodes, Eigen::Displacementd::Identity(), rodParameters);
	BOOST_CHECK( rodUnstableState1 );	
	// not singular 
	BOOST_CHECK( rodUnstableState1->integrate() );	
	// unstable
	BOOST_CHECK( !rodUnstableState1->isStable() );
}

BOOST_AUTO_TEST_CASE(ExtensibleRodStability3DTest_unstable2)
{
	qserl::rod3d::Parameters rodParameters;
	// set appropriate elasticity parameters
	rodParameters.radius = 0.01;
	rodParameters.length = 1.;
	const double youngModulus = 15.4e6;  /** Default Young modulus of rubber: 15.4 MPa */
	const double shearModulus = 5.13e6;  /** Default Shear modulus of rubber: 5.13 MPa */
  rodParameters.setIsotropicStiffnessCoefficientsFromElasticityParameters(youngModulus, shearModulus);
	rodParameters.density = 1.1 * 10e3;   /** 1.10 kg/dm3 -> kg/m3, */ 
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod3d::Parameters::RM_EXTENSIBLE_SHEARABLE;
	rodParameters.numNodes = 100;

	// unstable configuration
	Eigen::Wrenchd unstableConf2;
	unstableConf2[0] = -0.4294;
	unstableConf2[1] = -0.3144;
	unstableConf2[2] = -0.4496;
	unstableConf2[3] = -1.1369;
	unstableConf2[4] = -1.6963;
	unstableConf2[5] = -1.7618;
	qserl::rod3d::WorkspaceIntegratedStateShPtr rodUnstableState2 = qserl::rod3d::WorkspaceIntegratedState::create(unstableConf2,
		rodParameters.numNodes, Eigen::Displacementd::Identity(), rodParameters);
	BOOST_CHECK( rodUnstableState2 );	
	// not singular 
	BOOST_CHECK( rodUnstableState2->integrate() );	
	// unstable
	BOOST_CHECK( !rodUnstableState2->isStable() );
}

BOOST_AUTO_TEST_SUITE_END();

/* ------------------------------------------------------------------------- */
/* InextensibleRodStability3DTests																						 */
/* ------------------------------------------------------------------------- */
BOOST_AUTO_TEST_SUITE(InextensibleRodStability3DTests)

BOOST_AUTO_TEST_CASE(InextensibleRodStability3DTest_stable1)
{
	qserl::rod3d::Parameters rodParameters;
	// set appropriate elasticity parameters
	rodParameters.radius = 0.01;
	rodParameters.length = 1.;
	const double youngModulus = 15.4e6;  /** Default Young modulus of rubber: 15.4 MPa */
	const double shearModulus = 5.13e6;  /** Default Shear modulus of rubber: 5.13 MPa */
  rodParameters.setIsotropicStiffnessCoefficientsFromElasticityParameters(youngModulus, shearModulus);
	rodParameters.density = 1.1 * 10e3;   /** 1.10 kg/dm3 -> kg/m3, */ 
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod3d::Parameters::RM_INEXTENSIBLE;
	rodParameters.numNodes = 100;

	// stable configuration
	Eigen::Wrenchd stableConf1;
	stableConf1[0] = -0.3967;
	stableConf1[1] = 0.2774;
	stableConf1[2] = 0.1067;
	stableConf1[3] = 0.54;
	stableConf1[4] = 1.501;
	stableConf1[5] = 0.2606;
	qserl::rod3d::WorkspaceIntegratedStateShPtr rodStableState1 = qserl::rod3d::WorkspaceIntegratedState::create(stableConf1,
		rodParameters.numNodes, Eigen::Displacementd::Identity(), rodParameters);
	BOOST_CHECK( rodStableState1 );	
	// not singular 
	BOOST_CHECK( rodStableState1->integrate() );	
	// stable
	BOOST_CHECK( rodStableState1->isStable() );
}

BOOST_AUTO_TEST_CASE(InextensibleRodStability3DTest_stable2)
{
	qserl::rod3d::Parameters rodParameters;
	// set appropriate elasticity parameters
	rodParameters.radius = 0.01;
	rodParameters.length = 1.;
	const double youngModulus = 15.4e6;  /** Default Young modulus of rubber: 15.4 MPa */
	const double shearModulus = 5.13e6;  /** Default Shear modulus of rubber: 5.13 MPa */
  rodParameters.setIsotropicStiffnessCoefficientsFromElasticityParameters(youngModulus, shearModulus);
	rodParameters.density = 1.1 * 10e3;   /** 1.10 kg/dm3 -> kg/m3, */ 
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod3d::Parameters::RM_INEXTENSIBLE;
	rodParameters.numNodes = 100;

	// stable configuration 
	Eigen::Wrenchd stableConf2;
	stableConf2[0] = 0.5205;
	stableConf2[1] = 0.2989;
	stableConf2[2] = 0.0875;
	stableConf2[3] = 0.9518;
	stableConf2[4] = -0.8417;
	stableConf2[5] = -0.8075;
	qserl::rod3d::WorkspaceIntegratedStateShPtr rodStableState2 = qserl::rod3d::WorkspaceIntegratedState::create(stableConf2,
		rodParameters.numNodes, Eigen::Displacementd::Identity(), rodParameters);
	BOOST_CHECK( rodStableState2 );	
	// not singular 
	BOOST_CHECK( rodStableState2->integrate() );	
	// stable
	BOOST_CHECK( rodStableState2->isStable() );
}

BOOST_AUTO_TEST_CASE(InextensibleRodStability3DTest_unstable1)
{
	qserl::rod3d::Parameters rodParameters;
	// set appropriate elasticity parameters
	rodParameters.radius = 0.01;
	rodParameters.length = 1.;
	const double youngModulus = 15.4e6;  /** Default Young modulus of rubber: 15.4 MPa */
	const double shearModulus = 5.13e6;  /** Default Shear modulus of rubber: 5.13 MPa */
  rodParameters.setIsotropicStiffnessCoefficientsFromElasticityParameters(youngModulus, shearModulus);
	rodParameters.density = 1.1 * 10e3;   /** 1.10 kg/dm3 -> kg/m3, */ 
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod3d::Parameters::RM_INEXTENSIBLE;
	rodParameters.numNodes = 100;

	// unstable configuration
	Eigen::Wrenchd unstableConf1;
	unstableConf1[0] = -0.5885;
	unstableConf1[1] = -0.7467;
	unstableConf1[2] = 0.4277;
	unstableConf1[3] = -0.121;
	unstableConf1[4] = 0.0508;
	unstableConf1[5] = 0.9760;
	qserl::rod3d::WorkspaceIntegratedStateShPtr rodUnstableState1 = qserl::rod3d::WorkspaceIntegratedState::create(unstableConf1,
		rodParameters.numNodes, Eigen::Displacementd::Identity(), rodParameters);
	BOOST_CHECK( rodUnstableState1 );	
	// not singular 
	BOOST_CHECK( rodUnstableState1->integrate() );	
	// unstable
	BOOST_CHECK( !rodUnstableState1->isStable() );
}

BOOST_AUTO_TEST_CASE(InextensibleRodStability3DTest_unstable2)
{
	qserl::rod3d::Parameters rodParameters;
	// set appropriate elasticity parameters
	rodParameters.radius = 0.01;
	rodParameters.length = 1.;
	const double youngModulus = 15.4e6;  /** Default Young modulus of rubber: 15.4 MPa */
	const double shearModulus = 5.13e6;  /** Default Shear modulus of rubber: 5.13 MPa */
  rodParameters.setIsotropicStiffnessCoefficientsFromElasticityParameters(youngModulus, shearModulus);
	rodParameters.density = 1.1 * 10e3;   /** 1.10 kg/dm3 -> kg/m3, */ 
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod3d::Parameters::RM_INEXTENSIBLE;
	rodParameters.numNodes = 100;

	// unstable configuration
	Eigen::Wrenchd unstableConf2;
	unstableConf2[0] = -0.4294;
	unstableConf2[1] = -0.3144;
	unstableConf2[2] = -0.4496;
	unstableConf2[3] = -1.1369;
	unstableConf2[4] = -1.6963;
	unstableConf2[5] = -1.7618;
	qserl::rod3d::WorkspaceIntegratedStateShPtr rodUnstableState2 = qserl::rod3d::WorkspaceIntegratedState::create(unstableConf2,
		rodParameters.numNodes, Eigen::Displacementd::Identity(), rodParameters);
	BOOST_CHECK( rodUnstableState2 );	
	// not singular 
	BOOST_CHECK( rodUnstableState2->integrate() );	
	// unstable
	BOOST_CHECK( !rodUnstableState2->isStable() );
}

BOOST_AUTO_TEST_SUITE_END();

/* ------------------------------------------------------------------------- */
/* Extensible3DBencnhmarks																									*/
/* ------------------------------------------------------------------------- */
#ifndef _DEBUG

BOOST_AUTO_TEST_SUITE(Extensible3DBencnhmarks)

BOOST_AUTO_TEST_CASE(Extensible3DBencnhmark_1)
{
	qserl::rod3d::Parameters rodParameters;
	// set appropriate elasticity parameters
	rodParameters.radius = 0.01;
	rodParameters.length = 1.;
	const double youngModulus = 15.4e6;  /** Default Young modulus of rubber: 15.4 MPa */
	const double shearModulus = 5.13e6;  /** Default Shear modulus of rubber: 5.13 MPa */
  rodParameters.setIsotropicStiffnessCoefficientsFromElasticityParameters(youngModulus, shearModulus);
	rodParameters.density = 1.1 * 10e3;   /** 1.10 kg/dm3 -> kg/m3, */ 
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod3d::Parameters::RM_EXTENSIBLE_SHEARABLE;
	rodParameters.numNodes = 200;
	
	// set A-space bounds
	Eigen::Matrix<double, 6, 1> aSpaceUpperBounds, aSpaceLowerBounds;
	static const double maxTorque = 0.4;
	static const double maxForce = 1.8;
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
	qserl::rod3d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
	integrationOptions.computeJ_nu_sv = false;
	integrationOptions.stop_if_unstable = false;
	integrationOptions.keepMuValues = true;
	integrationOptions.keepJdet = true;
	integrationOptions.keepMMatrices = false;
	integrationOptions.keepJMatrices = true;

	qserl::util::TimePoint startBenchTime = qserl::util::getTimePoint();
	Eigen::Wrenchd wrench;
	static const int numSamples = 2000;
	int validSamples = 0;
	for (int i = 0 ; i < numSamples ; ++i)
	{
		for (int k = 0 ; k < 6 ; ++k)
			wrench[k] = ((rand() % 10000) * ( aSpaceUpperBounds[k] - aSpaceLowerBounds[k]) ) / 1.e4 + aSpaceLowerBounds[k];
		bool isStable = false;
		qserl::rod3d::WorkspaceIntegratedStateShPtr rodState = qserl::rod3d::WorkspaceIntegratedState::create(wrench,
			rodParameters.numNodes, Eigen::Displacementd::Identity(), rodParameters);
		rodState->integrationOptions(integrationOptions);
		BOOST_CHECK( rodState );	
		// not singular (zero volume, so should never happen by random sampling
		bool isNotSingular = rodState->integrate();
		BOOST_CHECK( isNotSingular );	
		if (isNotSingular && rodState->isStable())
			++validSamples;
	}
	double benchTimeMs = qserl::util::getElapsedTimeMsec(startBenchTime).count();
	BOOST_TEST_MESSAGE( "Benchmarking 3D extensible rods total time: " << benchTimeMs << "ms for " << numSamples << " integrated rods" );
	double benchTimePerRodUs = benchTimeMs * 1.e3 / static_cast<double>(numSamples);
	BOOST_TEST_MESSAGE( "  Avg. integration time per rod = " << benchTimePerRodUs << "us for " << rodParameters.numNodes << " nodes" );
	BOOST_TEST_MESSAGE( "  Avg. integration time per rod node = " << benchTimePerRodUs / static_cast<double>(rodParameters.numNodes) << "us" );
	double stabilityRatio = static_cast<double>(validSamples) / static_cast<double>(numSamples);
	// stability ratio should be arround 99%
	static const double minStabilityRatio = 0.9;
	BOOST_TEST_MESSAGE( "  Stability ratio = " << stabilityRatio << " (should be above " << minStabilityRatio << ")" );
	BOOST_CHECK( stabilityRatio > minStabilityRatio );	

}

BOOST_AUTO_TEST_SUITE_END();

#endif