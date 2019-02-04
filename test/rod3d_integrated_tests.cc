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

#include "qserl/rod3d/workspace_integrated_state.h"
#include "qserl/util/timer.h"
#include "util/lie_algebra_utils.h"


/* ------------------------------------------------------------------------- */
/* InextensibleRod3DRawDataTests																						 */
/* ------------------------------------------------------------------------- */
BOOST_AUTO_TEST_SUITE(InextensibleRod3DRawDataTests)

BOOST_AUTO_TEST_CASE(InextensibleRod3DRawDataTests_compare_external_data1)
{
  qserl::rod3d::Parameters rodParameters;
  // set appropriate elasticity parameters
  rodParameters.radius = 0.01;
  rodParameters.length = 1.;
  rodParameters.stiffnessCoefficients = Eigen::Matrix<double, 6, 1>::Ones();
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod3d::Parameters::RM_INEXTENSIBLE;
  rodParameters.numNodes = 200000;

  // stable configuration
  qserl::rod3d::Wrench stableConf1;
  stableConf1[0] = 5.7449;
  stableConf1[1] = -0.1838;
  stableConf1[2] = 3.7734;
  stableConf1[3] = -71.6227;
  stableConf1[4] = -15.6477;
  stableConf1[5] = 83.1471;
  qserl::rod3d::WorkspaceIntegratedStateShPtr rodStableState1 = qserl::rod3d::WorkspaceIntegratedState::create(
      stableConf1,
      rodParameters.numNodes,
      qserl::rod3d::Displacement::Identity(),
      rodParameters);
  BOOST_CHECK(rodStableState1);
  qserl::rod3d::WorkspaceIntegratedState::IntegrationResultT status = rodStableState1->integrate();
  // not singular
  BOOST_CHECK(status != qserl::rod3d::WorkspaceIntegratedState::IR_SINGULAR);
  // stable
  BOOST_CHECK(rodStableState1->isStable());
  BOOST_CHECK(status != qserl::rod3d::WorkspaceIntegratedState::IR_UNSTABLE);

  // compare q(1)
  Eigen::Matrix4d external_q1_mat;
  external_q1_mat << 0.131015450583688, -0.221694028519735, -0.966271452117711, 0.502570527667144,
      0.670919241188877, -0.697739138440250, 0.251057844719277, 0.224474097471111,
      -0.729858559768751, -0.681183316364205, 0.057314761488519, -0.515527857600866,
      0., 0., 0., 1.;
  const Eigen::Matrix4d impl_q1_mat = qserl::util::GetHomogenousMatrix(rodStableState1->nodes()[rodParameters.numNodes -
                                                                                                1]);
  for(int i = 0; i < 4; ++i)
  {
    for(int j = 0; j < 4; ++j)
    {
      BOOST_CHECK_CLOSE(external_q1_mat(i, j), impl_q1_mat(i, j), 0.1);
    }
  }

  // compare J(1)
  Eigen::Matrix<double, 6, 6> external_J1_mat;
  external_J1_mat
      << 0.863757430643127, 0.068239275411261, -0.311326745708524, -0.020469934592801, 0.008818885895480, -0.004128259697585,
      0.109312742029055, 0.060603348805207, -0.110555839878126, -0.016787411090804, 0.002089286520661, -0.010250408399195,
      -0.227912963349722, 0.005536330236083, 0.159710593966831, -0.016637254942726, 0.006340580390063, 0.024850086613713,
      -0.020469998324946, 0.008818707958365, -0.004128975395289, -0.001054045186435, -0.002919725344545, 0.001802353003749,
      -0.016788282449621, 0.002089215355955, -0.010250316008463, 0.001747674215122, 0.005553562171767, 0.002865488168595,
      -0.016635282303093, 0.006341004305033, 0.024849079504003, 0.004688562486302, -0.003222635977263, 0.006209331205126;

  const Eigen::Matrix<double, 6, 6>& impl_J1_mat = rodStableState1->getJMatrix(rodParameters.numNodes - 1);
  for(int i = 0; i < 6; ++i)
  {
    for(int j = 0; j < 6; ++j)
    {
      BOOST_CHECK_CLOSE(external_J1_mat(i, j), impl_J1_mat(i, j), 0.1);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END();

/* ------------------------------------------------------------------------- */
/* SingularConfigurations3DTests  																						 */
/* ------------------------------------------------------------------------- */
BOOST_AUTO_TEST_SUITE(SingularConfigurations3DTests)

BOOST_AUTO_TEST_CASE(SingularConfigurations3DTest_1)
{
  qserl::rod3d::Parameters rodDefaultParameters;

  // 1. ckeck null base wrench
  qserl::rod3d::Wrench
  singularWrench1(qserl::rod3d::Wrench::Zero());
  qserl::rod3d::WorkspaceIntegratedStateShPtr rodIntegratedState = qserl::rod3d::WorkspaceIntegratedState::create(
      singularWrench1,
      rodDefaultParameters.numNodes,
      qserl::rod3d::Displacement::Identity(),
      rodDefaultParameters);
  BOOST_CHECK(rodIntegratedState);
  qserl::rod3d::WorkspaceIntegratedState::IntegrationResultT status = rodIntegratedState->integrate();
  // singular
  BOOST_CHECK(status == qserl::rod3d::WorkspaceIntegratedState::IR_SINGULAR);
}

BOOST_AUTO_TEST_CASE(SingularConfigurations3DTest_2)
{
  qserl::rod3d::Parameters rodDefaultParameters;

  // 2. ckeck signular base wrench (a[1] = a[2] = a[4] = a[5] = 0)
  qserl::rod3d::Wrench
  singularWrench2(qserl::rod3d::Wrench::Zero());
  singularWrench2[0] = -1.5;
  singularWrench2[3] = 2.;
  qserl::rod3d::WorkspaceIntegratedStateShPtr rodIntegratedState = qserl::rod3d::WorkspaceIntegratedState::create(
      singularWrench2,
      rodDefaultParameters.numNodes,
      qserl::rod3d::Displacement::Identity(),
      rodDefaultParameters);
  BOOST_CHECK(rodIntegratedState);
  qserl::rod3d::WorkspaceIntegratedState::IntegrationResultT status = rodIntegratedState->integrate();
  // singular
  BOOST_CHECK(status == qserl::rod3d::WorkspaceIntegratedState::IR_SINGULAR);
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
  //rodParameters.density = 1.1 * 10e3;   /** 1.10 kg/dm3 -> kg/m3, */
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod3d::Parameters::RM_EXTENSIBLE_SHEARABLE;
  rodParameters.numNodes = 100;

  // stable configuration
  qserl::rod3d::Wrench stableConf1;
  stableConf1[0] = -0.3967;
  stableConf1[1] = 0.2774;
  stableConf1[2] = 0.1067;
  stableConf1[3] = 0.54;
  stableConf1[4] = 1.501;
  stableConf1[5] = 0.2606;
  qserl::rod3d::WorkspaceIntegratedStateShPtr rodStableState1 = qserl::rod3d::WorkspaceIntegratedState::create(
      stableConf1,
      rodParameters.numNodes,
      qserl::rod3d::Displacement::Identity(),
      rodParameters);
  BOOST_CHECK(rodStableState1);
  qserl::rod3d::WorkspaceIntegratedState::IntegrationResultT status = rodStableState1->integrate();
  // not singular
  BOOST_CHECK(status != qserl::rod3d::WorkspaceIntegratedState::IR_SINGULAR);
  // stable
  BOOST_CHECK(rodStableState1->isStable());
  BOOST_CHECK(status != qserl::rod3d::WorkspaceIntegratedState::IR_UNSTABLE);
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
  //rodParameters.density = 1.1 * 10e3;   /** 1.10 kg/dm3 -> kg/m3, */
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod3d::Parameters::RM_EXTENSIBLE_SHEARABLE;
  rodParameters.numNodes = 100;

  // stable configuration
  qserl::rod3d::Wrench stableConf2;
  stableConf2[0] = 0.5205;
  stableConf2[1] = 0.2989;
  stableConf2[2] = 0.0875;
  stableConf2[3] = 0.9518;
  stableConf2[4] = -0.8417;
  stableConf2[5] = -0.8075;
  qserl::rod3d::WorkspaceIntegratedStateShPtr rodStableState1 = qserl::rod3d::WorkspaceIntegratedState::create(
      stableConf2,
      rodParameters.numNodes,
      qserl::rod3d::Displacement::Identity(),
      rodParameters);
  BOOST_CHECK(rodStableState1);
  qserl::rod3d::WorkspaceIntegratedState::IntegrationResultT status = rodStableState1->integrate();
  // not singular
  BOOST_CHECK(status != qserl::rod3d::WorkspaceIntegratedState::IR_SINGULAR);
  // stable
  BOOST_CHECK(rodStableState1->isStable());
  BOOST_CHECK(status != qserl::rod3d::WorkspaceIntegratedState::IR_UNSTABLE);
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
  //rodParameters.density = 1.1 * 10e3;   /** 1.10 kg/dm3 -> kg/m3, */
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod3d::Parameters::RM_EXTENSIBLE_SHEARABLE;
  rodParameters.numNodes = 100;

  // unstable configuration
  qserl::rod3d::Wrench unstableConf1;
  unstableConf1[0] = -0.5885;
  unstableConf1[1] = -0.7467;
  unstableConf1[2] = 0.4277;
  unstableConf1[3] = -0.121;
  unstableConf1[4] = 0.0508;
  unstableConf1[5] = 0.9760;
  qserl::rod3d::WorkspaceIntegratedStateShPtr rodUnstableState1 = qserl::rod3d::WorkspaceIntegratedState::create(
      unstableConf1,
      rodParameters.numNodes,
      qserl::rod3d::Displacement::Identity(),
      rodParameters);
  BOOST_CHECK(rodUnstableState1);
  qserl::rod3d::WorkspaceIntegratedState::IntegrationResultT status = rodUnstableState1->integrate();
  // not singular
  BOOST_CHECK(status != qserl::rod3d::WorkspaceIntegratedState::IR_SINGULAR);
  // ubstable
  BOOST_CHECK(!rodUnstableState1->isStable());
  BOOST_CHECK(status == qserl::rod3d::WorkspaceIntegratedState::IR_UNSTABLE);
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
  //rodParameters.density = 1.1 * 10e3;   /** 1.10 kg/dm3 -> kg/m3, */
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod3d::Parameters::RM_EXTENSIBLE_SHEARABLE;
  rodParameters.numNodes = 100;

  // unstable configuration
  qserl::rod3d::Wrench unstableConf2;
  unstableConf2[0] = -0.4294;
  unstableConf2[1] = -0.3144;
  unstableConf2[2] = -0.4496;
  unstableConf2[3] = -1.1369;
  unstableConf2[4] = -1.6963;
  unstableConf2[5] = -1.7618;
  qserl::rod3d::WorkspaceIntegratedStateShPtr rodUnstableState2 = qserl::rod3d::WorkspaceIntegratedState::create(
      unstableConf2,
      rodParameters.numNodes,
      qserl::rod3d::Displacement::Identity(),
      rodParameters);
  BOOST_CHECK(rodUnstableState2);
  qserl::rod3d::WorkspaceIntegratedState::IntegrationResultT status = rodUnstableState2->integrate();
  // not singular
  BOOST_CHECK(status != qserl::rod3d::WorkspaceIntegratedState::IR_SINGULAR);
  // ubstable
  BOOST_CHECK(!rodUnstableState2->isStable());
  BOOST_CHECK(status == qserl::rod3d::WorkspaceIntegratedState::IR_UNSTABLE);
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
  //rodParameters.density = 1.1 * 10e3;   /** 1.10 kg/dm3 -> kg/m3, */
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod3d::Parameters::RM_INEXTENSIBLE;
  rodParameters.numNodes = 100;

  // stable configuration
  qserl::rod3d::Wrench stableConf1;
  stableConf1[0] = -0.3967;
  stableConf1[1] = 0.2774;
  stableConf1[2] = 0.1067;
  stableConf1[3] = 0.54;
  stableConf1[4] = 1.501;
  stableConf1[5] = 0.2606;
  qserl::rod3d::WorkspaceIntegratedStateShPtr rodStableState1 = qserl::rod3d::WorkspaceIntegratedState::create(
      stableConf1,
      rodParameters.numNodes,
      qserl::rod3d::Displacement::Identity(),
      rodParameters);
  BOOST_CHECK(rodStableState1);
  qserl::rod3d::WorkspaceIntegratedState::IntegrationResultT status = rodStableState1->integrate();
  // not singular
  BOOST_CHECK(status != qserl::rod3d::WorkspaceIntegratedState::IR_SINGULAR);
  // stable
  BOOST_CHECK(rodStableState1->isStable());
  BOOST_CHECK(status != qserl::rod3d::WorkspaceIntegratedState::IR_UNSTABLE);
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
  //rodParameters.density = 1.1 * 10e3;   /** 1.10 kg/dm3 -> kg/m3, */
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod3d::Parameters::RM_INEXTENSIBLE;
  rodParameters.numNodes = 100;

  // stable configuration
  qserl::rod3d::Wrench stableConf2;
  stableConf2[0] = 0.5205;
  stableConf2[1] = 0.2989;
  stableConf2[2] = 0.0875;
  stableConf2[3] = 0.9518;
  stableConf2[4] = -0.8417;
  stableConf2[5] = -0.8075;
  qserl::rod3d::WorkspaceIntegratedStateShPtr rodStableState2 = qserl::rod3d::WorkspaceIntegratedState::create(
      stableConf2,
      rodParameters.numNodes,
      qserl::rod3d::Displacement::Identity(),
      rodParameters);
  BOOST_CHECK(rodStableState2);
  qserl::rod3d::WorkspaceIntegratedState::IntegrationResultT status = rodStableState2->integrate();
  // not singular
  BOOST_CHECK(status != qserl::rod3d::WorkspaceIntegratedState::IR_SINGULAR);
  // stable
  BOOST_CHECK(rodStableState2->isStable());
  BOOST_CHECK(status != qserl::rod3d::WorkspaceIntegratedState::IR_UNSTABLE);
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
  //rodParameters.density = 1.1 * 10e3;   /** 1.10 kg/dm3 -> kg/m3, */
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod3d::Parameters::RM_INEXTENSIBLE;
  rodParameters.numNodes = 100;

  // unstable configuration
  qserl::rod3d::Wrench unstableConf1;
  unstableConf1[0] = -0.5885;
  unstableConf1[1] = -0.7467;
  unstableConf1[2] = 0.4277;
  unstableConf1[3] = -0.121;
  unstableConf1[4] = 0.0508;
  unstableConf1[5] = 0.9760;
  qserl::rod3d::WorkspaceIntegratedStateShPtr rodUnstableState1 = qserl::rod3d::WorkspaceIntegratedState::create(
      unstableConf1,
      rodParameters.numNodes,
      qserl::rod3d::Displacement::Identity(),
      rodParameters);
  BOOST_CHECK(rodUnstableState1);
  qserl::rod3d::WorkspaceIntegratedState::IntegrationResultT status = rodUnstableState1->integrate();
  // not singular
  BOOST_CHECK(status != qserl::rod3d::WorkspaceIntegratedState::IR_SINGULAR);
  // ubstable
  BOOST_CHECK(!rodUnstableState1->isStable());
  BOOST_CHECK(status == qserl::rod3d::WorkspaceIntegratedState::IR_UNSTABLE);
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
  //rodParameters.density = 1.1 * 10e3;   /** 1.10 kg/dm3 -> kg/m3, */
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod3d::Parameters::RM_INEXTENSIBLE;
  rodParameters.numNodes = 100;

  // unstable configuration
  qserl::rod3d::Wrench unstableConf2;
  unstableConf2[0] = -0.4294;
  unstableConf2[1] = -0.3144;
  unstableConf2[2] = -0.4496;
  unstableConf2[3] = -1.1369;
  unstableConf2[4] = -1.6963;
  unstableConf2[5] = -1.7618;
  qserl::rod3d::WorkspaceIntegratedStateShPtr rodUnstableState2 = qserl::rod3d::WorkspaceIntegratedState::create(
      unstableConf2,
      rodParameters.numNodes,
      qserl::rod3d::Displacement::Identity(),
      rodParameters);
  BOOST_CHECK(rodUnstableState2);
  qserl::rod3d::WorkspaceIntegratedState::IntegrationResultT status = rodUnstableState2->integrate();
  // not singular
  BOOST_CHECK(status != qserl::rod3d::WorkspaceIntegratedState::IR_SINGULAR);
  // ubstable
  BOOST_CHECK(!rodUnstableState2->isStable());
  BOOST_CHECK(status == qserl::rod3d::WorkspaceIntegratedState::IR_UNSTABLE);
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
  //rodParameters.density = 1.1 * 10e3;   /** 1.10 kg/dm3 -> kg/m3, */
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod3d::Parameters::RM_EXTENSIBLE_SHEARABLE;
  rodParameters.numNodes = 200;

  // set A-space bounds
  Eigen::Matrix<double, 6, 1> aSpaceUpperBounds, aSpaceLowerBounds;
  static const double maxTorque = 0.4;
  static const double maxForce = 1.8;
  for(int k = 0; k < 3; ++k)
  {
    aSpaceUpperBounds[k] = maxTorque;
    aSpaceLowerBounds[k] = -maxTorque;
  }
  for(int k = 3; k < 6; ++k)
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
  qserl::rod3d::Wrench wrench;
  static const int numSamples = 2000;
  int validSamples = 0;
  for(int i = 0; i < numSamples; ++i)
  {
    for(int k = 0; k < 6; ++k)
    {
      wrench[k] = ((rand() % 10000) * (aSpaceUpperBounds[k] - aSpaceLowerBounds[k])) / 1.e4 + aSpaceLowerBounds[k];
    }
    bool isStable = false;
    qserl::rod3d::WorkspaceIntegratedStateShPtr rodState = qserl::rod3d::WorkspaceIntegratedState::create(wrench,
                                                                                                          rodParameters.numNodes,
                                                                                                          qserl::rod3d::Displacement::Identity(),
                                                                                                          rodParameters);
    rodState->integrationOptions(integrationOptions);
    BOOST_CHECK(rodState);
    // not singular (zero volume, so should never happen by random sampling
    qserl::rod3d::WorkspaceIntegratedState::IntegrationResultT status = rodState->integrate();
    bool isNotSingular = status != qserl::rod3d::WorkspaceIntegratedState::IR_SINGULAR;
    BOOST_CHECK(isNotSingular);
    if(isNotSingular && rodState->isStable())
    {
      ++validSamples;
    }
  }
  double benchTimeMs = qserl::util::getElapsedTimeMsec(startBenchTime).count();
  BOOST_TEST_MESSAGE(
      "Benchmarking 3D extensible rods total time: " << benchTimeMs << "ms for " << numSamples << " integrated rods");
  double benchTimePerRodUs = benchTimeMs * 1.e3 / static_cast<double>(numSamples);
  BOOST_TEST_MESSAGE(
      "  Avg. integration time per rod = " << benchTimePerRodUs << "us for " << rodParameters.numNodes << " nodes");
  BOOST_TEST_MESSAGE(
      "  Avg. integration time per rod node = " << benchTimePerRodUs / static_cast<double>(rodParameters.numNodes)
                                                << "us");
  double stabilityRatio = static_cast<double>(validSamples) / static_cast<double>(numSamples);
  // stability ratio should be arround 99%
  static const double minStabilityRatio = 0.9;
  BOOST_TEST_MESSAGE("  Stability ratio = " << stabilityRatio << " (should be above " << minStabilityRatio << ")");
  BOOST_CHECK(stabilityRatio > minStabilityRatio);

}

BOOST_AUTO_TEST_SUITE_END();

#endif

/* ------------------------------------------------------------------------- */
/* InextensibleRod3DJacobiansTests																						 */
/* ------------------------------------------------------------------------- */
//BOOST_AUTO_TEST_SUITE(InextensibleRod3DJacobiansTests)
//
//BOOST_AUTO_TEST_CASE(InextensibleRod3DJacobiansTests_compare_external_data1)
//{
//	qserl::rod3d::Parameters rodParameters;
//	// set appropriate elasticity parameters
//	rodParameters.radius = 0.01;
//	rodParameters.length = 1.;
//	rodParameters.stiffnessCoefficients = Eigen::Matrix<double, 6, 1>::Ones();
//	rodParameters.integrationTime = 1.;
//	rodParameters.rodModel = qserl::rod3d::Parameters::RM_INEXTENSIBLE;
//	rodParameters.numNodes = 1000;
//
//	// stable configuration
//	Wrench stableConf1;
//	stableConf1[0] = 5.7449;
//	stableConf1[1] = -0.1838;
//	stableConf1[2] = 3.7734;
//	stableConf1[3] = -71.6227;
//	stableConf1[4] = -15.6477;
//	stableConf1[5] = 83.1471;
//	qserl::rod3d::WorkspaceIntegratedStateShPtr rodStableState1 = qserl::rod3d::WorkspaceIntegratedState::create(stableConf1,
//		rodParameters.numNodes, se3::SE3::Identity(), rodParameters);
//	BOOST_CHECK( rodStableState1 );	
//	qserl::rod3d::WorkspaceIntegratedState::IntegrationResultT status = rodStableState1->integrate();
//	// not singular 
//	BOOST_CHECK( status != WorkspaceIntegratedState::IR_SINGULAR );	
//	// stable
//	BOOST_CHECK( rodStableState1->isStable() );
//	BOOST_CHECK( status != WorkspaceIntegratedState::IR_UNSTABLE );	
//}
//
//BOOST_AUTO_TEST_SUITE_END();