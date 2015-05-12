=====================================================
Quasi-Static Elastic Rods Library (qserl) quick-start
=====================================================

:Author: Olivier Roussel
  
:Last update: 12/05/2015

.. |arw| unicode:: U+02794

-------------------------

This quick-start guide introduces very briefly how to use the qserl library to manipulate
quasi-static elastic rods. Note that rods are assumed to have constant circular cross-section.

Dependencies
>>>>>>>>>>>>

Qserl depends on the following packages:

 - Boost (>= 1.55)
 - Eigen with the unsupported module Lgsm.
   Eigen v3.2.4 including Lgsm is accessible through:
   git@github.com:olivier-roussel/eigen3.git
	
-------------------------
	
Planar rods
>>>>>>>>>>>

Planar (i.e. 2-dimensional) elastic rods code is packaged entirely inside the ``rod2d`` namespace,
itself inside the global ``qserl`` namespace of the library.

As a reminder, a planar quasi-static elastic rod can be parameterized globally by 3-dimensional coordinates,
which correspond to its internal wrench applied at its base (see Bre13_ for details). 
The rod is assumed to be handled by both ends, with given boundary conditions in SE(2).

The geometry of the rod is discretized in a given number of nodes and expressed in the frame of its base.
(see ``qserl::rod2d::WorkspaceState`` class for details).

In this section, we will see how to create a planar rod parameterized by its internal body wrench 
at its base and how to compute the rod entire geometry and internal wrenches
from this 3-dimensional parameterization.

The following code can be found in examples/example_rod_2D.cpp.

.. warning::

	In the planar case, the following limitations hold so far:
	
	 - the elasticity parameter cannot be defined
	 - the length of the rod must be unitary

.. code:: cpp

	// Requires the header "qserl/rod2d/rod.h"
	
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
					Eigen::Displacementd node3D = toDisplacement3D(*itNode);
					Eigen::Displacementd node3DGlobal = base3D * node3D;

					// ... do something with the current rod node position ...
		} 

	}

.. _Bre13: http://bretl.csl.illinois.edu/s/Bretl2014.pdf

-------------------------

3-dimensional rods
>>>>>>>>>>>>>>>>>>

3-dimensional elastic rods code is packaged entirely inside the ``rod3d`` namespace,
itself inside the global ``qserl`` namespace of the library.

3D Quasi-static elastic rods can be parameterized globally by 6-dimensional coordinates,
which correspond to its internal wrench applied at its base (see Bre13_ for details). 
The rod is assumed to be handled by both ends, with given boundary conditions in SE(3).

All static parameters of the rod are described in the Parameters structure.
This includes:

	- the radius of the rod (i.e. of the cross-section)
	- length of the rod
	- the deformation model which can be:
		
		- inextensible: only model deformations around the three rotations, i.e. torsion (X axis) and bending (Y and Z axis).
		- extensible and shearable: model deformations around the three rotations, i.e. torsion (X axis) and bending (Y and Z axis), and along the three translations, i.e. elongation/compression (X axis) and shearing (Y and Z axis).
		
	- stiffness coefficients:
	
	Defines the elasticity stiffness coefficient in all deformations directions. This is a 6-dimensional vector with the following index / deformation correspondance:
	
		- 0: torsional stiffness (X axis)
		- 1: bending stiffness (Y axis)
		- 2: bending stiffness (Z axis)
		- 3: elongation / compression stiffness (X axis). Only for extensible and shearable deformation model.
		- 4: shearing stiffness (Y axis). Only for extensible and shearable deformation model.
		- 5: shearing stiffness (Z axis). Only for extensible and shearable deformation model.
	
	- number of discretized nodes
	
The geometry of the rod is discretized in nodes and each is expressed in the local frame of its base.

The following code creates a 3-dimensional rod and computes its entire geometry and internal wrenches from its
internal wrench at its base (i.e. its 6-dimensional parameterization):

.. code:: cpp

	// Requires the header "qserl/rod3d/rod.h"
	
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
	Eigen::Wrenchd baseWrench;
	baseWrench.tx() = 5.7449;
	baseWrench.ty() = -0.1838;
	baseWrench.tz() = 3.7734;
	baseWrench.fx() = -71.6227;
	baseWrench.fy() = -15.6477;
	baseWrench.fz() = 83.1471;
	
	// ... the rod base position (here identity element of SE(3) ) ...
	const Eigen::Displacementd rodBasePosition(Eigen::Displacementd::Identity());
	
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
	rodBasePosition, integrationOptions);
	
	if (intResult != WorkspaceIntegratedState::IR_SINGULAR)
	{
		WorkspaceIntegratedStateShPtr integratedState1 = rod->integratedState();
		
		// here we can access to rod geometry (and more as wrenches as we defined some specific
		// integration options
		
		// For example, we can express the rod geometry (node by node) in the global frame
		// instead of rod base frame this way
		// Note this could be also directly given by the method WorkspaceState::nodesAbsolute6DPositions()
		for (std::vector<Eigen::Displacementd>::const_iterator itNode = integratedState1->nodes().begin();
			itNode != integratedState1->nodes().end(); ++itNode)
		{
			const Eigen::Displacementd nodeAbsolutePosition = rodBasePosition * (*itNode);
		
			// ... do something with the current rod node position ...
		}
	}


.. _Bre13: http://bretl.csl.illinois.edu/s/Bretl2014.pdf
