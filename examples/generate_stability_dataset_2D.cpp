// generate_stability_dataset_2D

#include <qserl/rod2d/rod.h>

#include "qserl/rod2d/analytic_q.h"
#include "qserl/rod2d/analytic_energy.h"

#include <fstream>
#include <boost/math/constants/constants.hpp>
#include <boost/chrono.hpp>

int
main()
{
#ifdef _OPENMP
  static const int numThreads = 24;
  omp_set_num_threads(numThreads);
  std::cout << "Using  OpenMP "<< _OPENMP << " with " << numThreads << " treads" << std::endl;

#endif

  // A-space bounds
  static const double maxTorque = 3 * boost::math::constants::pi<double>();
  static const double maxForce = 150;

  static const int numSamplesTorque = 200;
  static const int numSamplesForce = 200;

  static const size_t numNodes = 1001;    // discretization along the rod to check stability
  static const double kStabilityThreshold = 1.e-7;        /** Threshold for Jacobian determinant. */
  static const double kStabilityTolerance = 1.e-8;      /** Tolerance for which Jacobian determinant vanishes. */

  bool filterNonSurfacePoints = true;

  // set base A-space bounds
  static const int numSamplesTotal = numSamplesTorque * numSamplesForce * numSamplesForce;

  Eigen::Matrix<double, 3, 1> aSpaceUpperBounds, aSpaceLowerBounds;

  aSpaceUpperBounds[0] = maxTorque;
  aSpaceLowerBounds[0] = -maxTorque;
  for(int k = 1; k < 3; ++k)
  {
    aSpaceUpperBounds[k] = maxForce;
    aSpaceLowerBounds[k] = -maxForce;
  }
  static const double da_force = 2 * maxForce / static_cast<double>(numSamplesForce);
  static const double da_torque = 2 * maxTorque / static_cast<double>(numSamplesTorque);

  static const char* stabilityDatasetFilename = "stability_dataset_2D.txt";

  qserl::rod2d::Parameters rodParameters;
  rodParameters.radius = 1.;
  rodParameters.length = 1.;
  rodParameters.integrationTime = 1.;
  rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
  rodParameters.delta_t = 1. / static_cast<double>(numNodes);

  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
  integrationOptions.stop_if_unstable = true;
  integrationOptions.keepMuValues = false;
  integrationOptions.keepJdet = false;
  integrationOptions.keepMMatrices = false;
  integrationOptions.keepJMatrices = false;

  int successfullMotionConstants = 0;
  int done = 0;
  std::vector<Eigen::Vector3d> dataset_a_TXY(numSamplesTotal);
  std::vector<int> dataset_stability(numSamplesTotal, -1);
  std::vector<double> dataset_energy(numSamplesTotal, -1.);

  boost::chrono::high_resolution_clock::time_point startBenchTime = boost::chrono::high_resolution_clock::now();

  // proceed to samples of q1
  std::cout << "Starting generation of " << numSamplesTotal << " rod configurations for stability checking..";
#pragma omp parallel for schedule(dynamic)
  for(int idxSample = 0; idxSample < numSamplesTotal; ++idxSample)
  {
    const int i1 = (idxSample % (numSamplesTorque * numSamplesForce)) % numSamplesTorque;
    const int i2 = (idxSample % (numSamplesTorque * numSamplesForce)) / numSamplesTorque;
    const int i3 = idxSample / (numSamplesTorque * numSamplesForce);

    Eigen::Vector3d wrench_TXY;
    wrench_TXY[0] = i1 * da_torque;
    wrench_TXY[1] = i2 * da_force;
    wrench_TXY[2] = i3 * da_force;

    qserl::rod2d::MotionConstantsQ motionConstants;
    double energy = -1.;
    qserl::rod2d::WorkspaceIntegratedState::IntegrationResultT integrationStatus = qserl::rod2d::WorkspaceIntegratedState::IR_NUMBER_OF_INTEGRATION_RESULTS;
    if(qserl::rod2d::computeMotionConstantsQ(wrench_TXY, motionConstants))
    {
      ++successfullMotionConstants;
      qserl::rod2d::computeTotalElasticEnergy(motionConstants, energy);

      static const qserl::rod2d::Displacement2D identityDisp = {{0., 0., 0.}};
      qserl::rod2d::Wrench2D wrench_XYT;
      wrench_XYT[0] = wrench_TXY[1];
      wrench_XYT[1] = wrench_TXY[2];
      wrench_XYT[2] = wrench_TXY[0];
      qserl::rod2d::WorkspaceIntegratedStateShPtr rodState = qserl::rod2d::WorkspaceIntegratedState::create(wrench_XYT,
                                                                                                            identityDisp,
                                                                                                            rodParameters);
      rodState->integrationOptions(integrationOptions);
      integrationStatus = rodState->integrate();
    }
#pragma omp critical
    {
      // store to data array
      dataset_a_TXY[idxSample] = wrench_TXY;
      dataset_stability[idxSample] = static_cast<int>(integrationStatus);
      dataset_energy[idxSample] = energy;
      ++done;
      if(done % (numSamplesTotal / 100) == 0)
      {
        std::cout << "[PROGRESS] Done :" << static_cast<double>(done) * 100 / static_cast<double>(numSamplesTotal)
                  << std::endl;
      }
    }
  }
  std::cout << "DONE" << std::endl;
  double benchTimeMs = boost::chrono::duration_cast<boost::chrono::milliseconds>(
      boost::chrono::high_resolution_clock::now() - startBenchTime).count();

  std::cout << "Total time: " << benchTimeMs << "ms for " << numSamplesTotal << " samples" << std::endl;
  std::cout << "Per sample time =" << (benchTimeMs * 1.e3 / static_cast<double>(numSamplesTotal)) << "us" << std::endl;
  std::cout << "Singular samples =" << numSamplesTotal - successfullMotionConstants << " / " << numSamplesTotal
            << std::endl;

  if(filterNonSurfacePoints)
  {
    // brute force method complexity in O(n^2)
    static const int kNumNeighb = 27;
    static const double forceThresholdDistance = da_force * 1.1; // 10% error margin
    static const double torqueThresholdDistance = da_torque * 1.1; // 10% error margin

    std::cout << "[PROGRESS] Filtering non surface points..." << std::endl;
    std::vector<Eigen::Vector3d> dataset_a_TXY_filtered(numSamplesTotal);
    std::vector<int> dataset_stability_filtered(numSamplesTotal);
    std::vector<double> dataset_energy_filtered(numSamplesTotal);

    for(size_t i = 0; i < numSamplesTotal; ++i)
    {
      const int stabilityI = dataset_stability[i];
      const Eigen::Vector3d& nodeI = dataset_a_TXY[i];
      int numStableNeighb = 0;
      // if stable point, check if has less than n^3 - 1 stable neighbours
      // if so, then it is a surface point
      for(size_t j = 0; j < numSamplesTotal; ++j)
      {
        if(i != j && dataset_stability[j] == stabilityI)
        {
          const Eigen::Vector3d& nodeJ = dataset_a_TXY[j];
          if(abs(nodeI[0] - nodeJ[0]) < (torqueThresholdDistance &&
                                         abs(nodeI[1] - nodeJ[1]) < forceThresholdDistance &&
                                         abs(nodeI[2] - nodeJ[2]) < forceThresholdDistance))
          {
            ++numStableNeighb;
          }
        }
      } // end for neighb j
      if(numStableNeighb < kNumNeighb)
      {
        dataset_a_TXY_filtered[i] = nodeI;
        dataset_stability_filtered[i] = stabilityI;
        dataset_energy_filtered[i] = dataset_energy[i];
      }
      if(i % (numSamplesTotal / 100) == 0)
      {
        std::cout << "[PROGRESS] Done :" << static_cast<double>(i) * 100 / static_cast<double>(numSamplesTotal)
                  << std::endl;
      }

    }
    std::cout << "DONE" << std::endl;
  }

  // output to file
  std::cout << "Writing data to file..." << std::endl;

  // write headers to files
  std::ofstream ssamples(stabilityDatasetFilename, std::ofstream::out);
  ssamples << "# Planar rod case" << std::endl;
  ssamples << "# TXY convention" << std::endl;
  for(int k = 0; k < 3; ++k)
  {
    ssamples << "# A_low[" << k + 3 << "]=" << aSpaceLowerBounds[k] << " / A_upper[" << k + 3 << "]="
             << aSpaceUpperBounds[k] << std::endl;
  }

  ssamples << "# Number of samples (torque) = " << numSamplesTorque << std::endl;
  ssamples << "# Number of samples (force) = " << numSamplesForce << std::endl;
  ssamples << "#  Total number of samples = " << numSamplesTotal << std::endl;
  ssamples << "# delta_a3 (torque) = " << da_torque << std::endl;
  ssamples << "# delta_a[4,5] (force) = " << da_force << std::endl;
  ssamples << "# a3(T) a4(Fx) a5(Fy) stability energy" << std::endl;

  for(int idxSample = 0; idxSample < numSamplesTotal; ++idxSample)
  {
    ssamples << dataset_a_TXY[idxSample][0] << "\t"
             << dataset_a_TXY[idxSample][1] << "\t"
             << dataset_a_TXY[idxSample][2] << "\t"
             << dataset_stability[idxSample] << "\t"
             << dataset_energy[idxSample] << std::endl;
  }
  std::cout << "DONE" << std::endl;

  // pause
  char ch;
  std::cin >> ch;
}
