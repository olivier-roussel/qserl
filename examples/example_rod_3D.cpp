// Requires the header "qserl/rod3d/rod.h"
#include <qserl/rod3d/rod.h>

int
main()
{
  using namespace qserl;
  using namespace rod3d;

  // Define rod parameters
  Parameters rodParameters;
  rodParameters.radius = 0.01;                            // radius of the rod
  rodParameters.rodModel = Parameters::RM_INEXTENSIBLE;   // in 3D, could be inextensible rod model or extensible and shearable one.
  rodParameters.numNodes = 100;                          // The number of discretized rod nodes N where
  // N = (1/dt) + 1. Default value of 1e-2 should be sufficient for
  // most usages. As we are integrating numerically on a manifold,
  // keep in mind the error is not handled properly here.
  rodParameters.stiffnessCoefficients = Eigen::Matrix<double, 6, 1>::Ones();  // stiffness coefficients each deformation axis.
  // if inextensible model is used only the three first are relevant.
  // see documentation for more details.

  // Create the rod
  RodShPtr rod = Rod::create(rodParameters);

  // Define coordinates in A-space for the rod (its parameterization)...
  Wrench baseWrench;
  baseWrench[0] = 5.7449; // tx
  baseWrench[1] = -0.1838; // ty
  baseWrench[2] = 3.7734; // tz
  baseWrench[3] = -71.6227; // fx
  baseWrench[4] = -15.6477; // fy
  baseWrench[5] = 83.1471; // fz

  // ... the rod base position (here identity element of SE(3) ) ...
  const Displacement rodBasePosition(Displacement::Identity());

  // ... some integration options if needed ...
  // for example, here we want in addition to the rod geometry its internal wrenches
  // and Jacobian matrices (dq / da).
  WorkspaceIntegratedState::IntegrationOptions integrationOptions;
  integrationOptions.keepMuValues = true;
  integrationOptions.keepJMatrices = true;

  // ... and compute corresponding configuration
  // note that the given coordinates must not be singular (otherwise the returned
  // boolean isNotSingular will be false). See Rod::isConfigurationSingular() for details.
  WorkspaceIntegratedState::IntegrationResultT intResult = rod->integrateStateFromBaseWrench(baseWrench,
                                                                                             rodBasePosition,
                                                                                             integrationOptions);

  if(intResult != WorkspaceIntegratedState::IR_SINGULAR)
  {
    WorkspaceIntegratedStateShPtr integratedState1 = rod->integratedState();

    // here we can access to rod geometry (and more as wrenches as we defined some specific
    // integration options

    // For example, we can express the rod geometry (node by node) in the global frame
    // instead of rod base frame this way
    // Note this could be also directly given by the method WorkspaceState::nodesAbsolute6DPositions()
    for(Displacements::const_iterator itNode = integratedState1->nodes().begin();
        itNode != integratedState1->nodes().end(); ++itNode)
    {
      const Displacement nodeAbsolutePosition = rodBasePosition * (*itNode);

      // ... do something with the current rod node position ...
    }
  }
} 
