// Requires the header "qserl/rod2d/rod.h"
#include <qserl/rod2d/rod.h>


//
#include "qserl/rod2d/analytic_dqda.h"
#include "qserl/rod2d/workspace_integrated_state.h"
#include "qserl/rod2d/rod.h"
#include "util/lie_algebra_utils.h"
//
int main()
{
  Eigen::Vector3d wrench(-4.1337, 87.7116, 18.0966);
  static const double errorTolerance = 1.e-3; // tolerance on the error of dq(i) / da(j) between
                                              // analytical and numerically integrated expressions

  // Integrate numerically the jacobian dq / da
  // ----------------------------------
  qserl::rod2d::Parameters rodParameters;
	rodParameters.radius = 1.;
	rodParameters.length = 1.;
	rodParameters.integrationTime = 1.;
	rodParameters.rodModel = qserl::rod2d::Parameters::RM_INEXTENSIBLE;
	rodParameters.delta_t = 0.005;

	static const qserl::rod2d::Displacement2D identityDisp = { { 0., 0., 0. } };
  // set integration options
  qserl::rod2d::WorkspaceIntegratedState::IntegrationOptions integrationOptions;
  integrationOptions.stop_if_unstable = false;
	integrationOptions.keepMuValues = true;
	integrationOptions.keepJdet = true;
	integrationOptions.keepMMatrices = false;
	integrationOptions.keepJMatrices = true;

  qserl::rod2d::Wrench2D wrench2D;
  wrench2D[0] = wrench[1];
  wrench2D[1] = wrench[2];
  wrench2D[2] = wrench[0];

  qserl::rod2d::WorkspaceIntegratedStateShPtr rodState = qserl::rod2d::WorkspaceIntegratedState::create(wrench2D,
    identityDisp, rodParameters);
		rodState->integrationOptions(integrationOptions);
  qserl::rod2d::WorkspaceIntegratedState::IntegrationResultT integrationStatus = rodState->integrate();
  assert( integrationStatus == qserl::rod2d::WorkspaceIntegratedState::IR_VALID );	

  // Compute analytically the jacobian dq / da
  // ----------------------------------
  qserl::rod2d::MotionConstantsDqDa motionConstants;
  bool motionConstantsSuccess = qserl::rod2d::computeMotionConstantsDqDa(wrench, motionConstants);

  assert( motionConstantsSuccess );	

  const size_t numNodes = rodParameters.numberOfNodes();
  //for (size_t idxNode = 0; idxNode < numNodes; ++idxNode)
  size_t idxNode = numNodes - 1;
  {
    const double t =  static_cast<double>(idxNode) / static_cast<double>(numNodes - 1);
    Eigen::Matrix3d dqda_al;
    bool dqdaSuccess = qserl::rod2d::computeDqDaAtPositionT(t, motionConstants, dqda_al);
    
    assert( dqdaSuccess );	

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
	    //BOOST_CHECK_CLOSE( dqda_al(0,j), dtheta_da[j], errorTolerance );
      const double err1 = abs(dqda_al(0,j) - dtheta_da[jNum]);
      if (err1 > errorTolerance)
        std::cout << "err1=" << err1 << std::endl;
      // compare dq2 / da (i.e. dx / da)
	    //BOOST_CHECK_CLOSE( dqda_al(1,j), J_num_Tq[j](0,0), errorTolerance );
      const double err2 = abs(dqda_al(1,j) - J_num_Tq[jNum](0,2));
      if (err2 > errorTolerance)
        std::cout << "err2=" << err2 << std::endl;
      // compare dq3 / da (i.e. dy / da)
	    //BOOST_CHECK_CLOSE( dqda_al(2,j), J_num_Tq[j](0,1), errorTolerance );
      const double err3 = abs(dqda_al(2,j) - J_num_Tq[jNum](1,2));
      if (err3 > errorTolerance)
        std::cout << "err3=" << err3 << std::endl;
    }
  }	
  int toto = 1;
}

#if 0
int main ()
{
  using namespace qserl;
  using namespace rod2d;

  // Define rod parameters
  Parameters rodParameters;
  rodParameters.radius = 0.01;                            // radius of the rod
  rodParameters.rodModel = Parameters::RM_INEXTENSIBLE;   // in the planar case, only the inextensible rod model is implemented
  rodParameters.delta_t = 0.01;                           // the discretization resolution dt
                                                          // The number of discretized rod nodes will be
                                                          // N = (1/dt) + 1. Default value of 1e-2 should be sufficient for
                                                          // most usages.
  // Create the rod
  RodShPtr planarRod = Rod::create(rodParameters);

  // Define coordinates in A-space for the rod (its parameterization)...
  Wrench2D baseWrench;
  baseWrench[0] = -8.;
  baseWrench[1] = 24.;
  baseWrench[2] = 1.57;

  // ... the rod base position ...
  const Displacement2D rodBasePosition = { { 1., 0., -3.1416 } };

  // ... some integration options if needed ...
  // for example, here we want in addition to the rod geometry its internal wrenches
  // and Jacobian matrices (dq / da).
  WorkspaceIntegratedState::IntegrationOptions integrationOptions;
  integrationOptions.keepMuValues = true;
  integrationOptions.keepJMatrices = true;

  // ... and compute corresponding configuration
  // note that the given coordinates must not be singular (otherwise the returned
  // boolean isNotSingular will be false). See Rod::isConfigurationSingular() for details.
	WorkspaceIntegratedState::IntegrationResultT intResult = planarRod->integrateStateFromBaseWrench(baseWrench,
    rodBasePosition, integrationOptions);

	if (intResult != WorkspaceIntegratedState::IR_SINGULAR)
  {
    WorkspaceIntegratedStateShPtr integratedState1 = planarRod->integratedState();

    // here we can access to rod geometry (and more as wrenches as we defined some specific
    // integration options

    // For example, we can express the rod geometry (node by node) in the global frame
    // instead of rod base frame this way
    const Eigen::Displacementd base3D = toDisplacement3D(rodBasePosition);
    for (std::vector<Displacement2D>::const_iterator itNode = integratedState1->nodes().begin();
      itNode != integratedState1->nodes().end(); ++itNode)
    {
                // we will tranform our 2D 'pseudo'-displacement to 3D displacements
                // which will make transformations and manipulation easier
                const Eigen::Displacementd node3D = toDisplacement3D(*itNode);
                const Eigen::Displacementd node3DGlobal = base3D * node3D;

                // ... do something with the current rod node position ...
    }

  }
} //main

#endif
