/*

Flockuda: A numerical model of predator-prey dynamics based on a Molecular Dynamics approach.

Copyright (C) 2019 Christian Thomas Jacobs

*/

#include "predator.h"

__host__ void Predator::initialise(float kappa, float x0, float x1, float v0, float v1)
{   
    // Strength of bias (i.e. speed of attack).
    k = kappa;

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

__host__ void initialise_predator(Predator *p)
{
    float kappa = 4.0;
    p->initialise(kappa, 0, 0, 0, 0);
    return;
}

__host__ void write_predator(H5PartFile *output, Predator *p, int it)
{
    float x[1] = {p->x[0]};
    float y[1] = {p->x[1]};
    
    // Record the timestep.
    H5PartSetStep(output, it);

    // Write to .h5part file.
    H5PartWriteDataFloat32(output, "PredatorX", x);
    H5PartWriteDataFloat32(output, "PredatorY", y);
    return;
}

__host__ void save_predator(Predator *p)
{
    p->save();
    return;
}

__host__ void predator_velocity(Predator *p, float *c, float *xrandom, float dt)
{
    float xpc[2] = {c[0] - p->x[0], c[1] - p->x[1]};
    float magnitude = sqrt(pow(c[0] - p->x[0], 2) + pow(c[1] - p->x[1], 2));
    for(int d=0; d<2; ++d)
    {
        // Compute velocity as per Equation 2 of Lee et al. (2006).
        p->v[d] = (p->k/dt)*(xpc[d]/magnitude);
    }
    return;
}

__host__ void predator_location(Predator *p, float dt)
{
    float Lx = 1000.0;
    float Ly = 1000.0;
    for(int d=0; d<2; ++d)
    {
        // Compute location solving dx/dt = v as per Equation 2 of Lee et al. (2006).
        p->x[d] = p->xold[d] + dt*p->v[d];

        // Apply periodic boundary condition.
        if(d == 0 && p->x[d] > Lx)
        {
            p->x[d] -= Lx;
        }
        else
        {
            if (d == 0 && p->x[d] < 0)
            {
                p->x[d] += Lx;
            }
            else
            {
                if (d == 1 && p->x[d] > Ly)
                {
                    p->x[d] -= Ly;
                }
                else
                {
                    if (d == 1 && p->x[d] < 0)
                    {
                        p->x[d] += Ly;
                    }
                }
            }
        }
    }
    return;
}
