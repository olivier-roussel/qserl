from qserl.rod3d import *
import numpy as np

N = 100
rodParameters = Parameters()
rodParameters.radius = 0.01
rodParameters.rodModel = Parameters.RM_INEXTENSIBLE
rodParameters.numNodes = N
rodParameters.integrationTime = 1.

rodParameters.stiffnessCoefficients = np.ones((6,1))

rod = Rod.create(rodParameters)
# baseWrench = np.matrix ([ 5.7449, -0.1838, 3.7734, -71.6227, -15.6477, 83.1471,]).T
baseWrench = np.ones ((6,1))

rodBasePosition = np.eye(4)

integrationOptions = WorkspaceIntegratedState.IntegrationOptions ()
integrationOptions.keepJMatrices = True

intResult = rod.integrateStateFromBaseWrench(baseWrench, rodBasePosition, integrationOptions)
if intResult == WorkspaceIntegratedState.IR_VALID:
    state = rod.integratedState()

    oMd = state.node(N-1)
    oMd[:3,3] *= 0.5
    ik = InverseKinematics (rod)
    ik.maxIterations = 100
    ik.errorThreshold = 1e-3
    ik.verbosity = 1
    # Low scale improves converge success rate but increases
    # the number of iterations. (Between 0 and 1)
    ik.scale = 1.
    res = ik.compute (state, N-1, oMd)
