// Requires the header "qserl/rod2d/rod.h"
#include <qserl/rod2d/rod.h>


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
