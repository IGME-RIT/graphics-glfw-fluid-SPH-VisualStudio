#ifndef _DENSITY_H
#define _DENSITY_H

#include "Definitons.h"


#pragma region 
//===============================================================
//						DENSITY
//===============================================================
//The smoothing kernel for density calculation
float smoothKernelPoly6(glm::vec3 r)
{
	//This smoothing kernel is used to compute the density of the particle. 
	//This kernel forms a sort of a bell-curve, which is what we want for density. 
	//We need the density to be MAX value and not INFINITY and decrease as the distance increases from 0.
	// For more information and graphs, please refer to the papers listed above.
	/*
	W(r,h) = 315 * (h^2 - |r|^2)/64 * pi * h^9				 When 0<= |r| <= h
	W(r,h) = 0												 When |r| > h

	h = support radius
	*/
	float R = glm::length(r);
	float result = (315.0f * powf((H * H) - (R*R), 3)) / (64.0f * PI * powf(H, 9));

	return result;
}

//Calculate the change in density
float densityChange(Particle &r, Particle &p)
{
	return p.mass * smoothKernelPoly6(r.position - p.position);
}

// Update the density of the a single particle
void updateParticleDensity(Particle &p, int counterValue)
{
	std::vector<Particle *> neighbor;

	//Get the neighbouring particles of the selected particle.
	neighbor = neighbors[counterValue];

	float density = 0.0f;

	// Update density of the current particle
	for (int i = 0; i < neighbor.size(); i++)
	{
		//For each particle in the selected particle's vicinity, computeteh effect on density
		if (glm::length(p.position - neighbor[i]->position) < H)
			density += densityChange(p, *neighbor[i]);
	}

	p.density = density;

	//Since a lot forces are inversly proportional to density, we set to a value slightly higher than 0, 
	//if it is 0. This prevent divide by 0 errors. 
	if (density == 0.0f)
	{
		p.density = FLT_EPSILON;
	}
}

// Update the densities of all the particles
void updateDensities()
{
	//Call the function for each particle.
	for (int i = 0; i < Number_of_particels; i++)
	{
		updateParticleDensity(particles[i], i);
	}
}
//----------------------------------------------------------------
#pragma endregion Density

#pragma region
	//This function computes the gradient of the smooth Kernel
glm::vec3 smoothKernelPoly6Gradient(glm::vec3 r)
{
	/*
	gradient W(r,h) = -945 * r * (h^2 - |r|^2)/ 32 * pi * h
	*/
	float R = glm::length(r);
	glm::vec3 result = r * (-954.0f * powf((H*H) - (R*R), 2));

	result /= (32 * PI * powf(H, 9));

	return result;
}

//This function computes the laplacian of the smooth Kernel
float smoothKernelPoly6Laplacian(glm::vec3 r)
{
	/*
	Laplacian W(r,h) = -945 * (h^2 - |r|^2) * (3h^2 - 7|r|^2) / 32 * pi * h^9
	*/
	float R = glm::length(r);
	float L = -945 / (32 * PI * powf(H, 9));

	L *= ((H * H) - (R * R)) * ((3.0f * H * H) - (7.0f * R * R));

	return L;
}

float smoothColorField(Particle &r, Particle &p)
{
	return p.mass * smoothKernelPoly6(r.position - p.position) / p.density;
}
#pragma endregion Surface Tension





#endif _DENSITY_H





