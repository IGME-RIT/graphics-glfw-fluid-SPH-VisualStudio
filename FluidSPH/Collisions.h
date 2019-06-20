/*
File Name: Collisions.cpp
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
#ifndef _RESOLVE_COLLISION_H
#define _RESOLVE_COLLISION_H

#include "Definitons.h"
#include "Physics.h"


void boundVelocities()
{
	std::vector<Particle *>::iterator it;

	// This function goes through all the edge grids.
	// For all the particles there, it checks if they leave the bounding area. 
	// If they are outside the bounding area, then check if they continue to 
	// move outwards of the bounding volume, then change the component of velocity
	// which is along the surface normal.

	for (int i = 0; i < Grid_Size; i++)
	{
		for (int j = 0; j < Grid_Size; j++)
		{
			//X-axis
			for (it = grid[0][i][j].begin(); it != grid[0][i][j].end(); it++)
			{
				if ((*it)->position.x < 0 && ((*it)->velocity.x < 0 || (*it)->acceleration.x < 0))
				{
					(*it)->velocity.x *= DAMPENING_CONSTANT;
					(*it)->position.x = 0.0f;
				}
			}
			for (it = grid[Grid_Size - 1][i][j].begin(); it != grid[Grid_Size - 1][i][j].end(); it++)
			{
				if ((*it)->position.x > BoundarySizeX && ((*it)->velocity.x > 0 || (*it)->acceleration.x > 0))
				{
					(*it)->velocity.x *= DAMPENING_CONSTANT;
					(*it)->position.x = BoundarySizeX;
				}
			}
			//Y-Axis
			for (it = grid[i][0][j].begin(); it != grid[i][0][j].end(); it++)
			{
				if ((*it)->position.y < 0 && ((*it)->velocity.y < 0 || (*it)->acceleration.y < 0))
				{
					(*it)->velocity.y *= -0.1f;// DAMPENING_CONSTANT;
					(*it)->position.y = 0.0f;
				}
			}

			for (it = grid[i][Grid_Size - 1][j].begin(); it != grid[i][Grid_Size - 1][j].end(); it++)
			{
				if ((*it)->position.y > BoundarySizeY && ((*it)->velocity.y > 0 || (*it)->acceleration.y > 0))
				{
					(*it)->velocity.y *= DAMPENING_CONSTANT;
					(*it)->position.y = BoundarySizeY;
				}
			}

			//Z-axis
			for (it = grid[i][j][0].begin(); it != grid[i][j][0].end(); it++)
			{
				if ((*it)->position.z < 0 && ((*it)->velocity.z < 0 || (*it)->acceleration.z < 0))
				{
					(*it)->velocity.z *= DAMPENING_CONSTANT;
					(*it)->position.z = 0.0f;
				}
			}
			for (it = grid[i][j][Grid_Size - 1].begin(); it != grid[i][j][Grid_Size - 1].end(); it++)
			{
				if ((*it)->position.z > BoundarySizeZ && ((*it)->velocity.z > 0 || (*it)->acceleration.z > 0))
				{
					(*it)->velocity.z *= DAMPENING_CONSTANT;
					(*it)->position.z = BoundarySizeZ;
				}
			}

		}
	}


	//NOT needed

	//for (int i = 0; i < Number_of_particels; i++)
	//{
	//	if (particles[i].position.y < 0 && particles[i].velocity.y < 0)
	//	{
	//		particles[i].velocity.y *= -0.7f;
	//		particles[i].position.y = 0.0f;
	//		//particles[i].acceleration.y += CONTAINER_SPRING_CONSTANT;
	//		//particles[i].acceleration.y += particles[i].velocity.y * CONTAINER_DAMPENING;
	//	}

	//	if (particles[i].position.x > BoundarySize && particles[i].velocity.x > 0)
	//	{
	//		particles[i].velocity.x *= -0.7f;
	//		particles[i].position.x = BoundarySize;
	//		//particles[i].acceleration.x -= CONTAINER_SPRING_CONSTANT;
	//		//particles[i].acceleration.x -= particles[i].velocity.x * CONTAINER_DAMPENING;
	//	}

	//	if (particles[i].position.x < 0 && particles[i].velocity.x < 0)
	//	{
	//		particles[i].velocity.x *= -0.7f;
	//		particles[i].position.x = 0.0f;
	//		//particles[i].acceleration.x += CONTAINER_SPRING_CONSTANT;
	//		//particles[i].acceleration.x += particles[i].velocity.x * CONTAINER_DAMPENING;
	//	}
	//}
}

void updateVelocities()
{
	/*
	This function updates the acceleration of each particle.
	*/

	std::vector<Particle *>::iterator it;
	float k;													// Smoothed Color
	glm::vec3 Fpressure, Fviscosity, n, Fexternal, Fsurface, Finternal, Ftotal;
	float l;
	for (int i = 0; i < Number_of_particels; i++)
	{
		//Reset the values in forces
		Fpressure = glm::vec3(0.0f);
		Fviscosity = glm::vec3(0.0f);
		Fexternal = glm::vec3(0.0f);
		Finternal = glm::vec3(0.0f);
		Fsurface = glm::vec3(0.0f);
		k = 0;
		n = Fpressure;

		for (it = neighbors[i].begin(); it != neighbors[i].end(); it++)
		{
			l = glm::length(particles[i].position - (*it)->position);

			if (l <= H && l > 0.0f)
			{
				Fpressure += pressureForcePerParticle(particles[i], *(*it));
				Fviscosity += viscosityForcePerParticle(particles[i], *(*it));
				// n is the direction of surface tension force.
				// For particles which are not in the outside surface teh n will sum upto 0.
				// For all the particles in the surface, the value will be non zero.
				n += ((*it)->mass * smoothKernelPoly6Gradient(particles[i].position - (*it)->position)) / (*it)->density;
				k += ((*it)->mass * smoothKernelPoly6Laplacian(particles[i].position - (*it)->position)) / (*it)->density;
			}
		}

		Fpressure *= -1.0f;
		Fviscosity *= particles[i].viscosity;

		float CFNLength = glm::length(n);
		//Calculate the surface tension force
		if (CFNLength > COLOR_FIELD_THRESHOLD && CFNLength != 0.0f)
		{
			Fsurface = (-SIGMA) * k * (n / CFNLength);
		}

		Finternal = Fviscosity + Fpressure;

		Fexternal = (G * particles[i].density) + Fsurface;

		Ftotal = Finternal + Fexternal;

		particles[i].acceleration = Ftotal / particles[i].density;
	}

	boundVelocities();
}

glm::vec3 EulerIntegrator(glm::vec3 pos, float h, glm::vec3 &velocity, glm::vec3 acc)
{
	glm::vec3 P;

	velocity += h * acc;

	//Calculate the displacement in that time step with the current velocity.
	P = pos + (h * velocity);

	//return the position P
	return P;
}


void integrate(float dt)
{
	for (int i = 0; i < Number_of_particels; i++)
	{
		particles[i].position = EulerIntegrator(particles[i].position, dt, particles[i].velocity, particles[i].acceleration);
	}
}



#endif _RESOLVE_COLLISION_H
