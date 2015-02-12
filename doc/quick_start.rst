=====================================================
Quasi-Static Elastic Rods Library (qserl) quick-start
=====================================================

:Author: Olivier Roussel
  
:Date: 12/02/2015

.. |arw| unicode:: U+02794

-------------------------

This quick-start guide introduces very briefly how to use the qserl library to manipulate
quasi-static elastic rods.

Dependencies
>>>>>>>>>>>>

Qserl depends on the following packages:

 - Boost (>= 1.55)
 - Eigen with the unsupported module Lgsm.
   An Eigen v3.1 with Lgsm fork is accessible through:
   ssh://hg@bitbucket.org/oroussel/eigen
	
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
	rodParameters.radius = 0.01;				// radius of the rod
	rodParameters.rodModel = Parameters::RM_INEXTENSIBLE; 	// in the planar case, only the inextensible rod model is implemented 
	rodParameters.delta_t = 0.01;				// the discretization resolution dt
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
	bool isNotSingular = planarRod->integrateStateFromBaseWrench(baseWrench, rodBasePosition,
		integrationOptions);		

	if (isNotSingular)
	{
		WorkspaceIntegratedStateShPtr integratedState = rod->integratedState();
		
		// here we can access to rod geometry (and more as wrenches as we defined some specific
		// integration options
		
		// For example, we can express the rod geometry (node by node) in the global frame
		// instead of rod base frame this way
		for (std::vector<Displacement2D>::const_iterator itNode = integratedState->nodes().begin();
			itNode != integratedState->nodes().end(); ++itNode)
		{
			// we will tranform our 2D 'pseudo'-displacement to 3D displacements
			// which will make transformations and manipulation easier

		}
		
	}
	
.. _Bre13: http://bretl.csl.illinois.edu/s/Bretl2014.pdf

-------------------------

3-dimensional rods
>>>>>>>>>>>>>>>>>>

TODO
