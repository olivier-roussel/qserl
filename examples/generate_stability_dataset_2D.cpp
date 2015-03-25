// generate_stability_dataset_2D

#include <qserl/rod2d/rod.h>

#include "qserl/rod2d/analytic_q.h"
#include "qserl/rod2d/analytic_dqda.h"

#include <fstream>
#include <boost/math/constants/constants.hpp>

int main()
{
  // A-space bounds
  static const double maxTorque = 2 * boost::math::constants::pi<double>();
	static const double maxForce = 100;

  static const int numSamplesTorque = 10;
  static const int numSamplesForce = 10;

  static const size_t numNodes = 1001;		// discretization along the rod to check stability
  static const double kStabilityThreshold = 1.e-7;				/** Threshold for Jacobian determinant. */
  static const double kStabilityTolerance = 1.e-8;			/** Tolerance for which Jacobian determinant vanishes. */

	// set base A-space bounds
	Eigen::Matrix<double, 3, 1> aSpaceUpperBounds, aSpaceLowerBounds;

  aSpaceUpperBounds[0] = maxTorque;
  aSpaceLowerBounds[0] = -maxTorque;
	for (int k = 1 ; k < 3 ; ++k)
	{
		aSpaceUpperBounds[k] = maxForce;
		aSpaceLowerBounds[k] = -maxForce;
	}

  static const std::string stabilityDatasetFilename("stability_dataset_2D.txt");

  // write headers to files
  std::ofstream ssamples(stabilityDatasetFilename, std::ofstream::out);
  ssamples << "# Planar rod case" << std::endl;
	for (int k = 0 ; k < 3 ; ++k)
	{
		ssamples << "# A_low[" << k+3 << "]=" << aSpaceLowerBounds[k] << " / A_upper[" << k+3 << "]=" << aSpaceUpperBounds[k] << std::endl;
	}

	ssamples << "# Number of samples (torque) = " << numSamplesTorque << std::endl;
	ssamples << "# Number of samples (force) = " << numSamplesForce << std::endl;
  ssamples << "# a3 a4 a5 stable energy" << std::endl;

  
  Eigen::Vector3d wrench;
  qserl::rod2d::MotionConstantsQ motionConstants;

  qserl::rod2d::Parameters rodParameters;
	rodParameters.radius = 1.;
	rodParameters.length = 1.;
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
	rodParameters.delta_t = 1. / static_cast<double>(numNodes);

  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
  integrationOptions.stop_if_unstable = false;
	integrationOptions.keepMuValues = true;
	integrationOptions.keepJdet = true;
	integrationOptions.keepMMatrices = false;
	integrationOptions.keepJMatrices = true;

	//qserl::util::TimePoint startBenchTime = qserl::util::getTimePoint();

  // proceed to samples of q1
  std::cout << "Starting generation of " << numSamples << " for the rod tip position...";
	for (int idxSample = 0 ; idxSample < numSamples ; 
	{
    for (int k = 0; k < 3; ++k)
    {
      wrench[k] = ((rand() % 10000) * ( aSpaceUpperBounds[k] - aSpaceLowerBounds[k]) ) / 1.e4 + 
        aSpaceUpperBounds[k];
    }

    if (qserl::rod2d::computeMotionConstantsQ(wrench, motionConstants))
    {
      //++successfullMotionConstants;
      if (qserl::rod2d::computeQAtPositionT(1., wrench, motionConstants, q1_dot, q1))
      {
        //static const qserl::rod2d::Displacement2D identityDisp = { { 0., 0., 0. } };
        //qserl::rod2d::Wrench2D wrench2D;
        //wrench2D[0] = wrench[1];
        //wrench2D[1] = wrench[2];
        //wrench2D[2] = wrench[0];
        //qserl::rod2d::WorkspaceIntegratedStateShPtr rodState = qserl::rod2d::WorkspaceIntegratedState::create(wrench2D,
        //  identityDisp, rodParameters);
        //rodState->integrationOptions(integrationOptions);
        //qserl::rod2d::WorkspaceIntegratedState::IntegrationResultT integrationStatus = rodState->integrate();

        // check stability 
        bool isStable = true;
        bool isThresholdOn = false;
        double prevDet = 0.;
        for (size_t idxNode = 0 ; idxNode < numNodes && isStable; ++idxNode)
        {
          const double t = static_cast<double>(idxNode) / static_cast<double>(numNodes-1);
          if (qserl::rod2d::computeDqDaAtPositionT(t, motionConstants, dqda))
          {
            const double detJ = dqda.determinant();
            if (abs(detJ) > kStabilityThreshold)
              isThresholdOn = true;
            if (isThresholdOn && ( abs(detJ) < kStabilityTolerance || 
              (idxNode > 0 && detJ * prevDet < 0.) ) )
              isStable = false;
            else
              prevDet = detJ;
          }
        }
        //if (isStable && integrationStatus != qserl::rod2d::WorkspaceIntegratedState::IR_VALID)
        //  std::cout << "[DEBUG] found stable analytically but not with numerical" << std::endl;

        //if (!isStable && integrationStatus!= qserl::rod2d::WorkspaceIntegratedState::IR_UNSTABLE)
        //  std::cout << "[DEBUG] found unstable analytically but not with numerical" << std::endl;

        if (isStable)
        {
          ++successSamples;

          // normalize the rotation
          const double theta_init = q1[0];
          q1[0] = atan2(sin(theta_init), cos(theta_init));

          ssamples << wrench[0] << '\t' << wrench[1] << '\t' << wrench[2] << '\t';
          ssamples << q1[0] << '\t' << q1[1] << '\t' << q1[2] << std::endl;
        }
      }
    }
  }
  std::cout << "DONE" << std::endl;

  // proceed to initial guesses
  std::cout << "Starting generation of " << numInitialGuesses << " initial guesses...";
  int successGuesses = 0;
  while (successGuesses < numInitialGuesses)
  {
    for (int k = 0; k < 3; ++k)
    {
      wrench[k] = ((rand() % 10000) * ( aSpaceUpperBounds[k] - aSpaceLowerBounds[k]) ) / 1.e4 + 
        aSpaceUpperBounds[k];
    }

    if (qserl::rod2d::computeMotionConstantsDqDa(wrench, motionConstants))
    {
      // check stability 
      bool isStable = true;
      bool isThresholdOn = false;
      double prevDet = 0.;
      for (size_t idxNode = 0 ; idxNode < numNodes && isStable; ++idxNode)
      {
        const double t = static_cast<double>(idxNode) / static_cast<double>(numNodes-1);
        if (qserl::rod2d::computeDqDaAtPositionT(t, motionConstants, dqda))
        {
          const double detJ = dqda.determinant();
          if (abs(detJ) > kStabilityThreshold)
            isThresholdOn = true;
          if (isThresholdOn && ( abs(detJ) < kStabilityTolerance || 
            (idxNode > 0 && detJ * prevDet < 0.) ) )
            isStable = false;
          else
            prevDet = detJ;
        }
      }
      if (isStable)
      {
        ++successGuesses;
        sguesses << wrench[0] << '\t' << wrench[1] << '\t' << wrench[2] << std::endl;
      }
    }
  }
  std::cout << "DONE" << std::endl;
}
