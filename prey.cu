#include "prey.h"
#include "forces.h"

__device__ void Prey::initialise(float mass, float x0, float x1)
{   
    // Mass.
    m = mass;

    // Location.
    x[0] = x0;
    x[1] = x1;

    // Velocity.
    v[0] = 0.0;
    v[1] = 0.0;

    save();

    return;
}

__device__ void Prey::save()
{   
    // Save information from the previous timestep.
    xold[0] = x[0];
    xold[1] = x[1];
    vold[0] = v[0];
    vold[1] = v[1];

    return;
}

__global__ void initialise_prey(Prey *p, float *xrandom, float *yrandom, int nprey, float Lx, float Ly)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float mass = 1.0;
    if (i < nprey)
    {
        p[i].initialise(mass, Lx*xrandom[i], Ly*yrandom[i]);
    }
    return;
}

__host__ void write_prey(H5PartFile *output, Prey *p, int nprey, int it)
{
    H5PartSetStep(output, it);

    float x[nprey];
    float y[nprey];
    for(int i=0; i < nprey; ++i)
    {
        x[i] = p[i].x[0];
        y[i] = p[i].x[1];
    }
    H5PartWriteDataFloat32(output, "PreyX", x);
    H5PartWriteDataFloat32(output, "PreyY", y);
}

__global__ void save_prey(Prey *p, int nprey)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nprey)
    {
        p[i].save();
    }

    return;
}

__global__ void prey_velocity(Prey *p, int nprey, float xp0, float xp1, float *xrandom, float dt)
{
    float f[2];

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nprey)
    {
        for(int d=0; d<2; ++d)
        {
            // Compute force terms.
            f[d] = prey_alignment(p, nprey, d) + prey_attraction(p, nprey, d) + prey_repulsion(p, nprey, d) - prey_friction(p, nprey, d) - prey_avoid(p, nprey, xp0, xp1, d);
            //printf("%f \n", prey_friction(p, nprey, d));

            // Compute velocity using F = ma.
            p[i].v[d] = p[i].vold[d] + dt*(f[d]/p[i].m);
        }
    }

    return;
}

__global__ void prey_location(Prey *p, int nprey, float dt)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float Lx = 1000.0;
    float Ly = 1000.0;
    if (i < nprey)
    {
        for(int d=0; d<2; ++d)
        {
            // Compute location solving dx/dt = v.
            p[i].x[d] = p[i].xold[d] + dt*p[i].v[d];

            if(d == 0 && p[i].x[d] > Lx)
            {
                p[i].x[d] -= Lx;
                p[i].v[d] = 0;
            }
            else
            {
                if (d == 0 && p[i].x[d] < 0)
                {
                    p[i].x[d] += Lx;
                    p[i].v[d] = 0;
                }
                else
                {
                    if (d == 1 && p[i].x[d] > Ly)
                    {
                        p[i].x[d] -= Ly;
                        p[i].v[d] = 0;
                    }
                    else
                    {
                        if (d == 1 && p[i].x[d] < 0)
                        {
                            p[i].x[d] += Ly;
                            p[i].v[d] = 0;
                        }
                    }
                }
            }
        }
    }
    return;
}

__host__ void prey_centre(Prey *p, int nprey, float *c)
{
    c[0] = 0;
    c[1] = 0;
    for(int i=0; i<nprey; ++i)
    {
        c[0] += p[i].xold[0];
        c[1] += p[i].xold[1];
    }
    c[0] /= nprey;
    c[1] /= nprey;
    return;
}
