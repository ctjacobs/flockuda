#include "predator.h"

__host__ void Predator::initialise(float kappa, float x0, float x1, float v0, float v1)
{   
    // Strength of bias.
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
    
    H5PartSetStep(output, it);
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
    float xpc[2];
    float magnitude = sqrt(pow(c[0] - p->x[0], 2) + pow(c[1] - p->x[1], 2));

    for(int d=0; d<2; ++d)
    {
        xpc[0] = c[0] - p->x[0];
        xpc[1] = c[1] - p->x[1];

        // Compute velocity using F = ma.
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
        // Compute location solving dx/dt = v.
        p->x[d] = p->xold[d] + dt*p->v[d];

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
