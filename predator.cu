/*

Flockuda: A numerical model of predator-prey dynamics based on the Molecular Dynamics approach of Lee et al. (2006).

Copyright (C) 2019 Christian T. Jacobs

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/


#include "predator.h"
#include "configuration.h"


__host__ void Predator::initialise(float x0, float x1, float v0, float v1)
{
    // Location.
    x[0] = x0;
    x[1] = x1;

    // Velocity.
    v[0] = v0;
    v[1] = v1;

    save();

    return;
}

__host__ void Predator::save()
{   
    // Save information from the previous timestep.
    xold[0] = x[0];
    xold[1] = x[1];
    vold[0] = v[0];
    vold[1] = v[1];

    return;
}

__host__ void initialise_predator(Predator *p, Configuration config)
{
    p->initialise(0, 0, 0, 0);
    return;
}

__host__ void write_predator(H5PartFile *output, Predator *p, int it)
{
    float x[1] = {p->x[0]};
    float y[1] = {p->x[1]};
    float z[1] = {0};
    
    // Record the timestep.
    H5PartSetStep(output, it);

    // Write to .h5part file.
    H5PartWriteDataFloat32(output, "PredatorX", x);
    H5PartWriteDataFloat32(output, "PredatorY", y);
    H5PartWriteDataFloat32(output, "PredatorZ", z);
    return;
}

__host__ void save_predator(Predator *p)
{
    p->save();
    return;
}

__host__ void predator_velocity(Predator *p, Configuration config, float *c, float *xrandom)
{
    float xpc[2] = {c[0] - p->x[0], c[1] - p->x[1]};
    float magnitude = sqrt(pow(c[0] - p->x[0], 2) + pow(c[1] - p->x[1], 2));
    for(int d=0; d<2; ++d)
    {
        // Compute velocity as per Equation 2 of Lee et al. (2006).
        p->v[d] = (config.kappa/config.dt)*(xpc[d]/magnitude);
    }
    return;
}

__host__ void predator_location(Predator *p, Configuration config)
{
    for(int d=0; d<2; ++d)
    {
        // Compute location solving dx/dt = v as per Equation 2 of Lee et al. (2006).
        p->x[d] = p->xold[d] + config.dt*p->v[d];

        // Apply periodic boundary condition.
        if(d == 0 && p->x[d] > config.Lx)
        {
            p->x[d] -= config.Lx;
        }
        else
        {
            if (d == 0 && p->x[d] < 0)
            {
                p->x[d] += config.Lx;
            }
            else
            {
                if (d == 1 && p->x[d] > config.Ly)
                {
                    p->x[d] -= config.Ly;
                }
                else
                {
                    if (d == 1 && p->x[d] < 0)
                    {
                        p->x[d] += config.Ly;
                    }
                }
            }
        }
    }
    return;
}
