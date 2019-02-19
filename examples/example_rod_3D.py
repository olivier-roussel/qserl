from qserl.rod3d import *
import numpy as np

rodParameters = Parameters()
rodParameters.radius = 0.01
rodParameters.rodModel = Parameters.RM_INEXTENSIBLE
rodParameters.numNodes = 100
rodParameters.integrationTime = 1
# N = (1/dt) + 1. Default value of 1e-2 should be sufficient for
# most usages. As we are integrating numerically on a manifold,
# keep in mind the error is not handled properly here.
rodParameters.stiffnessCoefficients = np.ones((6,1)) # stiffness coefficients each deformation axis.
# if inextensible model is used only the three first are relevant.
# see documentation for more details.

# Create the rod
rod = Rod.create(rodParameters)

# Define coordinates in A-space for the rod (its parameterization)...
baseWrench = np.ones ((6,1))
baseWrench[0] = 5.7449; # tx
baseWrench[1] = -0.1838; # ty
baseWrench[2] = 3.7734; # tz
baseWrench[3] = -71.6227; # fx
baseWrench[4] = -15.6477; # fy
baseWrench[5] = 83.1471; # fz

# ... the rod base position (here identity element of SE(3) ) ...
rodBasePosition = np.eye(4)

# ... some integration options if needed ...
# for example, here we want in addition to the rod geometry its internal wrenches
# and Jacobian matrices (dq / da).
integrationOptions = WorkspaceIntegratedState.IntegrationOptions()
integrationOptions.keepMuValues = True
integrationOptions.keepJMatrices = True

# ... and compute corresponding configuration
# note that the given coordinates must not be singular (otherwise the returned
# boolean isNotSingular will be false). See Rod::isConfigurationSingular() for details.
intResult = rod.integrateStateFromBaseWrench(baseWrench,
        rodBasePosition,
        integrationOptions)

if intResult != WorkspaceIntegratedState.IR_SINGULAR:
  integratedState1 = rod.integratedState()

  # here we can access to rod geometry (and more as wrenches as we defined some specific
  # integration options

  # For example, we can express the rod geometry (node by node) in the global frame
  # instead of rod base frame this way
  # Note this could be also directly given by the method WorkspaceState::nodesAbsolute6DPositions()
  for i in range (integratedState1.numNodes()):
    node = integratedState1.node(i)
    nodeAbsolutePosition = rodBasePosition * node;
    # ... do something with the current rod node position ...
