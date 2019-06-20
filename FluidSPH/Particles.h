
/*
File Name: Particles.h
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

#ifndef _PARTICLES_H
#define _PARTICLES_H

#include "Definitons.h"

void setup()
{
	float divisionX = BoundarySizeX / Grid_Size, divisionY = BoundarySizeY / Grid_Size, divisionZ = BoundarySizeZ / Grid_Size;

	int i, j, k, l;
	for (i = 0; i < Number_of_particels; i++)
	{
		particles[i].position.x = (float)(i % (Grid_Size - 1))*divisionX + (0.1f);	//* (i%2)
		particles[i].position.y = (float)(i / (Grid_Size - 1))*divisionY + 1.0f;
		particles[i].position.z = (float)(i / ((Grid_Size - 1)* (Grid_Size - 1)))*divisionZ + 0.1f;

		particles[i].density = DENSITY;
		particles[i].mass = MASS;
		particles[i].viscosity = VISCOSITY;
		particles[i].velocity = glm::vec3(0.0f);
		particles[i].acceleration = glm::vec3(0.0f);
	}

	for (i = 0; i < Grid_Size; i++)
	{
		for (j = 0; j < Grid_Size; j++)
		{
			for (k = 0; k < Grid_Size; k++)
			{
				grid[i][j][k].clear();
			}
		}
	}
}

void setParticles()
{
	float divisionX = BoundarySizeX / Grid_Size, divisionY = BoundarySizeY / Grid_Size, divisionZ = BoundarySizeZ / Grid_Size;

	int i, j, k, l;
	for (i = 0; i < Number_of_particels; i++)
	{
		particles[i].position.x = (float)(i % (Grid_Size - 1))*divisionX + (0.1f);	//* (i%2)
		particles[i].position.y = (float)(i / (Grid_Size - 1))*divisionY + 1.0f;
		particles[i].position.z = (float)(i / ((Grid_Size - 1)* (Grid_Size - 1)))*divisionZ + 0.1f;

		particles[i].density = DENSITY;
		particles[i].mass = MASS;
		particles[i].viscosity = VISCOSITY;
		particles[i].velocity = glm::vec3(0.0f);
		particles[i].acceleration = glm::vec3(0.0f);
	}

	for (i = 0; i < Grid_Size; i++)
	{
		for (j = 0; j < Grid_Size; j++)
		{
			for (k = 0; k < Grid_Size; k++)
			{
				grid[i][j][k].clear();
			}
		}
	}
}


void clear_tree()
{
	// Clear the grid. ( matrix of vectors )
	int i, j, k;
	for (i = 0; i < Grid_Size; i++)
	{
		for (j = 0; j < Grid_Size; j++)
		{
			for (k = 0; k < Grid_Size; k++)
			{
				grid[i][j][k].clear();
			}
		}
	}
}

void catergorizeParticles()
{
	// This function categorizes each particle into the respectivel grid using its position's first value as the index for the grid.
	clear_tree();
	Particle p;
	float divisionX = BoundarySizeX / Grid_Size, divisionY = BoundarySizeY / Grid_Size, divisionZ = BoundarySizeZ / Grid_Size;
	int x, y, z;
	for (int i = 0; i < Number_of_particels; i++)
	{
		p = particles[i];
		x = (int)floor(p.position.x / (divisionX));
		y = (int)floor(p.position.y / (divisionY));
		z = (int)floor(p.position.z / (divisionZ));

		x = (x < 0) ? 0 : x;
		y = (y < 0) ? 0 : y;
		z = (z < 0) ? 0 : z;

		x = (x > Grid_Size - 1) ? Grid_Size - 1 : x;
		y = (y > Grid_Size - 1) ? Grid_Size - 1 : y;
		z = (z > Grid_Size - 1) ? Grid_Size - 1 : z;

		grid[x][y][z].push_back(&particles[i]);


	}
}

std::vector<Particle *>  getNeighborsforPoint(int x, int y, int z, Particle r)
{
	// Get neighbors for a particle in a specific grid. 

	std::vector<Particle *> p;
	bool xFlagL = false, yFlagL = false, zFlagL = false;
	bool xFlagM = false, yFlagM = false, zFlagM = false;
	std::vector<Particle *>::iterator it;

	p.clear();

	for (it = grid[x][y][z].begin(); it != grid[x][y][z].end(); it++)
	{
		p.push_back((*it));
	}

	// We have to get the particles within the support raius. So, we take the number of grids which encompass 1 support radius. (H/boundarySize)
	for (int i = 1; i <= (H*Grid_Size) / BoundarySizeX; i++)
	{
		xFlagL = false; yFlagL = false; zFlagL = false;
		xFlagM = false; yFlagM = false; zFlagM = false;
		// X-axis
		if (x - i >= 0)
		{
			for (it = grid[x - i][y][z].begin(); it != grid[x - i][y][z].end(); it++)
			{
				p.push_back((*it));
			}
			xFlagL = true;
		}
		if (x + i < Grid_Size)
		{
			for (it = grid[x + i][y][z].begin(); it != grid[x + i][y][z].end(); it++)
			{
				p.push_back((*it));
			}
			xFlagM = true;
		}

		//Y axis
		if (y - i >= 0)
		{
			for (it = grid[x][y - i][z].begin(); it != grid[x][y - i][z].end(); it++)
			{
				p.push_back((*it));
			}
			yFlagL = true;

		}
		if (y + i < Grid_Size)
		{
			for (it = grid[x][y + i][z].begin(); it != grid[x][y + i][z].end(); it++)
			{
				p.push_back((*it));
			}
			yFlagM = true;
		}

		// Z-axis
		if (z - i >= 0)
		{
			for (it = grid[x][y][z - i].begin(); it != grid[x][y][z - i].end(); it++)
			{
				p.push_back((*it));
			}
			zFlagL = true;
		}
		if (z + i < Grid_Size)
		{
			for (it = grid[x][y][z + i].begin(); it != grid[x][y][z + i].end(); it++)
			{
				p.push_back((*it));
			}
			zFlagM = true;
		}

		/*
		if (xFlagL && yFlagL)
		{
		for (it = grid[x - i][y - i][z].begin(); it != grid[x - i][y - i][z].end(); it++)
		{
		p.push_back((*it));
		}

		if (zFlagL)
		{
		for (it = grid[x - i][y - i][z - i].begin(); it != grid[x - i][y - i][z - i].end(); it++)
		{
		p.push_back((*it));
		}
		}
		if (zFlagM)
		{
		for (it = grid[x - i][y - i][z + i].begin(); it != grid[x - i][y - i][z + i].end(); it++)
		{
		p.push_back((*it));
		}
		}
		}

		if (xFlagM && yFlagM)
		{
		for (it = grid[x + i][y + i][z].begin(); it != grid[x + i][y + i][z].end(); it++)
		{
		p.push_back((*it));
		}

		if (zFlagL)
		{
		for (it = grid[x + i][y + i][z - i].begin(); it != grid[x + i][y + i][z - i].end(); it++)
		{
		p.push_back((*it));
		}
		}
		if (zFlagM)
		{
		for (it = grid[x + i][y + i][z + i].begin(); it != grid[x + i][y + i][z + i].end(); it++)
		{
		p.push_back((*it));
		}
		}
		}

		if (xFlagL && yFlagM)
		{
		for (it = grid[x - i][y + i][z].begin(); it != grid[x - i][y + i][z].end(); it++)
		{
		p.push_back((*it));
		}

		if (zFlagL)
		{
		for (it = grid[x - i][y + i][z - i].begin(); it != grid[x - i][y + i][z - i].end(); it++)
		{
		p.push_back((*it));
		}
		}
		if (zFlagM)
		{
		for (it = grid[x - i][y + i][z + i].begin(); it != grid[x - i][y + i][z + i].end(); it++)
		{
		p.push_back((*it));
		}
		}
		}

		if (xFlagM && yFlagL)
		{
		for (it = grid[x + i][y - i][z].begin(); it != grid[x + i][y - i][z].end(); it++)
		{
		p.push_back((*it));
		}

		if (zFlagL)
		{
		for (it = grid[x + i][y - i][z - i].begin(); it != grid[x + i][y - i][z - i].end(); it++)
		{
		p.push_back((*it));
		}
		}
		if (zFlagM)
		{
		for (it = grid[x + i][y - i][z + i].begin(); it != grid[x + i][y - i][z + i].end(); it++)
		{
		p.push_back((*it));
		}
		}
		}*/
	}
	return p;
}


//get all the neighbors for respective particle
void getNeighbors()
{
	Particle p;
	float divisionX = BoundarySizeX / Grid_Size, divisionY = BoundarySizeY / Grid_Size, divisionZ = BoundarySizeZ / Grid_Size;
	int x, y, z;

	for (int i = 0; i < Number_of_particels; i++)
	{
		p = particles[i];

		x = (int)floor(p.position.x / divisionX);
		y = (int)floor(p.position.y / divisionY);
		z = (int)floor(p.position.z / divisionZ);

		if (x == 1)
			x = 1;
		/*if (!((x < 0) || (x > Grid_Size - 1) || (y < 0) || (y > Grid_Size - 1) || (z < 0) || (z > Grid_Size - 1)))
		neighbors[i] = getNeighborsforPoint(x, y, z, particles[i]);
		else
		neighbors[i].clear();*/

		x = (x < 0) ? 0 : x;
		y = (y < 0) ? 0 : y;
		z = (z < 0) ? 0 : z;

		x = (x > Grid_Size - 1) ? Grid_Size - 1 : x;
		y = (y > Grid_Size - 1) ? Grid_Size - 1 : y;
		z = (z > Grid_Size - 1) ? Grid_Size - 1 : z;

		neighbors[i] = getNeighborsforPoint(x, y, z, particles[i]);
	}
}



#endif _PARTICLES_H