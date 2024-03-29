/*

Flockuda: A numerical model of predator-prey dynamics based on the Molecular Dynamics approach of Lee et al. (2006).

Copyright (C) 2019 Christian T. Jacobs

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/


#include "prey.h"
#include "forces.h"
#include <stdio.h>


__host__ void Prey::initialise(float x0, float x1)
{
    // Location.
    x[0] = x0;
    x[1] = x1;

    // Velocity.
    v[0] = 0.0;
    v[1] = 0.0;

    save();

    return;
}

__host__ void Prey::save()
{   
    // Save information from the previous timestep.
    xold[0] = x[0];
    xold[1] = x[1];
    vold[0] = v[0];
    vold[1] = v[1];

    return;
}

__host__ void initialise_prey(Prey *p, Configuration config, float *xrandom, float *yrandom)
{
    // Initialise prey at random locations throughout the domain.
    for(int i=0; i < config.nprey; ++i)
    {
        p[i].initialise(config.Lx*xrandom[i], config.Ly*yrandom[i]);
    }
    return;
}

__host__ void write_prey(H5PartFile *output, Prey *p, Configuration config, int it)
{
    // Record the timestep.
    H5PartSetStep(output, it);

    // Collect all prey location data into X and Y arrays.
    float x[config.nprey];
    float y[config.nprey];
    float z[config.nprey];
    for(int i=0; i < config.nprey; ++i)
    {
        x[i] = p[i].x[0];
        y[i] = p[i].x[1];
        z[i] = 0;
    }

    // Write to .h5part file.
    H5PartWriteDataFloat32(output, "PreyX", x);
    H5PartWriteDataFloat32(output, "PreyY", y);
    H5PartWriteDataFloat32(output, "PreyZ", z);
}

__host__ void save_prey(Prey *p, Configuration config)
{
    for(int i=0; i < config.nprey; ++i)
    {
        p[i].save();
    }

    return;
}

__global__ void prey_velocity(Prey *p, Configuration config, float xp0, float xp1)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float f[2];  // The force acting on the prey.

    if (i < config.nprey)
    {
        for(int d=0; d<2; ++d)
        {
            // Compute force terms.
            f[d] = prey_alignment(p, config, d) + prey_attraction(p, config, d) + prey_repulsion(p, config, d) - prey_friction(p, config, d) - prey_avoid(p, config, xp0, xp1, d);

            // Compute velocity using F = ma as per Equation 1 of Lee et al. (2006).
            p[i].v[d] = p[i].vold[d] + config.dt*(f[d]/config.mass);
        }
    }

    return;
}

__global__ void prey_location(Prey *p, Configuration config)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < config.nprey)
    {
        for(int d=0; d<2; ++d)
        {
            // Compute location solving dx/dt = v as per Equation 1 of Lee et al. (2006).
            p[i].x[d] = p[i].xold[d] + config.dt*p[i].v[d];

            // Apply periodic boundary condition.
            if(d == 0 && p[i].x[d] > config.Lx)
            {
                p[i].x[d] -= config.Lx;
                p[i].v[d] = 0;
            }
            else
            {
                if (d == 0 && p[i].x[d] < 0)
                {
                    p[i].x[d] += config.Lx;
                    p[i].v[d] = 0;
                }
                else
                {
                    if (d == 1 && p[i].x[d] > config.Ly)
                    {
                        p[i].x[d] -= config.Ly;
                        p[i].v[d] = 0;
                    }
                    else
                    {
                        if (d == 1 && p[i].x[d] < 0)
                        {
                            p[i].x[d] += config.Ly;
                            p[i].v[d] = 0;
                        }
                    }
                }
            }
        }
    }
    return;
}

__host__ void prey_centre(Prey *p, Configuration config, float *centre)
{
    // The centre of the prey flock based on computing the average of all prey locations.
    centre[0] = 0;
    centre[1] = 0;
    for(int i=0; i<config.nprey; ++i)
    {
        centre[0] += p[i].xold[0];
        centre[1] += p[i].xold[1];
    }
    centre[0] /= config.nprey;
    centre[1] /= config.nprey;
    return;
}
