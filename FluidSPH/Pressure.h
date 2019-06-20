/*
File Name: Pressure.cpp
Copyright © 2018
Original authors: Srinivasan Thiagarajan, Sanketh Bhat, Benjamin Evans
Refactored by Sanketh Bhat, Benjamin Evans
Written under the supervision of David I. Schwartz, Ph.D., and
supported by a professional development seed grant from the B. Thomas
Golisano College of Computing & Information Sciences
(https://www.rit.edu/gccis) at the Rochester Institute of Technology.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Description:
In this example we demonstrate the use of SPH to emulate fluid motion.
Fluid mechanics can be implemented in 2 ways: Eularian based or Lagragian
based. Eularian based fluid simulation are done using a grid. There are points
on the grid, and the particles contained within the grid follow a specific set
of rules. In Eularian approach, you need to account for conservation of mass
explicitly.
The lagragian approach accoutns for conservation of mass implicitly, Since each
cluster of particles interact with each other and are separatly accountable.
The lagragian approach, considers the forces caused by all the surrounding particles.
It interpolates between the position of the surrounding particles to get the overall
force acting on the selected particle. This process of interpolation is called SPH.

In SPH, we use smoothing kernels to inetrpolate based on the distance from the particle.
We use different kernels for different "aspects" of fluid properties. For pressure, we
implement a spike kernel, as the pressure sould increase drastically as the distance
gets smaller. But we use a poly6 smoothing kernel for density distributions and surface
tension. We use the gradient or laplacian of the kernel, based on whichever one is more applicable.
For more info see the referenced papers.

In SPH fluid simulation, the particles each consist of mass, velocity and acceleration.
The particles expereince density change, forces due to pressure, Viscosity, surface tension
and collision amongst themselves.

Press "SHIFT" to start simulation
Use "SPACE" to toggle gravity in x-axis, or use "W" to toggle gravity in y-axis.

References:
Nicholas Gallagher
Lagrangian Fluid Dynamics Using Smoothed Particles Hydrodynamics by Micky Kelager
Particle-Based Fluid Simulation for Interactive Applications by Matthias Muller, David Charypar and Markus Gross
*/
#ifndef _PRESSURE_H
#define _PRESSURE_H

#include "Definitons.h"

#pragma region
//================================================================
//						PRESSURE
//================================================================
//spike kernel for calculating the pressure force
glm::vec3 spikeKernelPoly6Gradient(glm::vec3 r)
{
	/*
	We are using spike kernel to smooth pressure. We are using the spike
	kernel because we need the pressure to increase alsmot exponentially as
	the distance between the two positons decreases.
	*/
	glm::vec3 grad(0.0f);;
	float R = glm::length(r);
	grad = glm::normalize(r) * (H - R) * (H - R) * (-45.0f);
	grad /= (PI * powf(H, 6));// *fmax(R, FLT_EPSILON));

	return grad;
}

//Calculate the force experienced by a particle due to another particle
glm::vec3 pressureForcePerParticle(Particle &r, Particle &p)
{
	/*
	PV = nRT

	n = mass/MolarMAss = 1000g /18 = 55.55555
	R = 0.0083144621(75) amu (km/s)2 K−1
	T= 293.15 K

	V = mass/density
	P = nRT * mass / density
	*/

	float P1 = K * 13.533444f  * r.density / (r.mass), P2 = K * 13.533444f * p.density / (p.mass);		// Calculating the pressure.

																										/*
																										Here we compute the force caused due to pressure difference between the particles.
																										It is done by calculating the pressure difference between the two positions,
																										and use the spike kernel to calculatethe force.
																										*/
	glm::vec3 fp = (P1 + P2) * p.mass * spikeKernelPoly6Gradient(r.position - p.position) / (2.0f * p.density);

	return fp;
}
//-----------------------------------------------------------------

#pragma endregion Pressure


#endif _PRESSURE_H

