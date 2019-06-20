/*
File Name: main.cpp
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

#include "GLRender.h"
#include "Physics.h"
#include "Collisions.h"



#pragma region Global Data member
// Global data members
// This is your reference to your shader program.
// This will be assigned with glCreateProgram().
// This program will run on your GPU.


double time = 0.0;
double timebase = 0.0;
double accumulator = 0.0;
double physicsStep = 0.022; // This is the number of milliseconds we intend for the physics to update.

bool start = false;
float radius = 0.25f;

							

#pragma endregion




#pragma region util_functions
// Functions called between every frame. game logic

// This runs once every physics timestep.
void update(float t)
{
	//Catergorize the particles into their respective grids.
	catergorizeParticles();

	//Each particles collects info on the particles surrounding them
	getNeighbors();

	//Update the densities at each particle location
	updateDensities();

	//update the acceleration of each particle
	if (start)
		updateVelocities();

	//Resolve collisions
	//findAndResolveCollisions();

	//Integrate the particle (update the 
	integrate(t);
}

// This runs once every frame to determine the FPS and how often to call update based on the physics step.
void checkTime()
{
	// Get the current time.
	time = glfwGetTime();

	// Get the time since we last ran an update.
	double dt = time - timebase;

	// If more time has passed than our physics timestep.
	if (dt > physicsStep)
	{

		timebase = time; // set new last updated time

						 // Limit dt
		if (dt > 0.25)
		{
			dt = 0.25;
		}
		accumulator += dt;

		// Update physics necessary amount
		while (accumulator >= physicsStep)
		{
			update(physicsStep);
			accumulator -= physicsStep;
		}
	}
}

// This function runs every frame

void renderParticles()
{
	
	glLineWidth(1.0f);

	glm::vec3 p;
	glm::mat4 m;
	glEnable(GL_POINT_SMOOTH);
	glUniformMatrix4fv(uniMVP, 1, GL_FALSE, glm::value_ptr(PV));
	glColor3f(1.0f, 1.0f, 1.0f);
	glPointSize(POINTSIZE);

	glBegin(GL_POINTS);

	for (int i = 0; i < Number_of_particels; i++)
	{
		glVertex3fv((float*)&particles[i].position);
	}
	glEnd();
}


// This function is used to handle key inputs.
// It is a callback funciton. i.e. glfw takes the pointer to this function (via function pointer) and calls this function every time a key is pressed in the during event polling.
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_SPACE && (action == GLFW_PRESS))
	{
		if (G.x >= 0)
		{

			view = glm::lookAt(glm::vec3(BoundarySizeX / 2.0f, 0.5f, 3.0f), glm::vec3(BoundarySizeX / 2.0f, 0.5f, 0.0f), glm::vec3(1.0f, 1.0f, 0.0f));

			PV = proj * view;

			G.x = -14.8f;
		}
		else
		{
			view = glm::lookAt(glm::vec3(BoundarySizeX / 2.0f, 0.5f, 3.0f), glm::vec3(BoundarySizeX / 2.0f, 0.5f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));

			PV = proj * view;

			G.x = 0;
		}

	}


	if (key == GLFW_KEY_W && (action == GLFW_PRESS || action == GLFW_REPEAT))
	{
		if (G.y == 0)
		{
			G.y = -9.8f;
		}
		else
			G.y = 0.0f;
	}

	if (key == GLFW_KEY_LEFT_SHIFT && (action == GLFW_PRESS || action == GLFW_REPEAT))
	{
		start = true;
	}

}
#pragma endregion

void main()
{
	// Initializes most things needed before the main loop
	init();

	// Sends the funtion as a funtion pointer along with the window to which it should be applied to.
	glfwSetKeyCallback(window, key_callback);

	std::cout << "\n Press \"SHIFT\" to start simulation.";
	std::cout << "\n Use \"SPACE\" to toggle gravity in x - axis.";
	std::cout << "\n use \"W\" to toggle gravity in y - axis.";


	setup();
	// Enter the main loop.
	while (!glfwWindowShouldClose(window))
	{
		
		//checkTime();
		update(physicsStep);

		// Call the render function(s).
		renderScene();

		renderParticles();

		// Swaps the back buffer to the front buffer
		// Remember, you're rendering to the back buffer, then once rendering is complete, you're moving the back buffer to the front so it can be displayed.
		glfwSwapBuffers(window);

		// Checks to see if any events are pending and then processes them.
		glfwPollEvents();
	}

	//Cleans shaders and the program and frees up GLFW memory
	cleanup();

	return;
}