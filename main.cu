/* Flockuda

Copyright (C) 2019 Christian Thomas Jacobs

*/

#include <iostream>
#include <H5Part.h>
#include <cuda.h>
#include <curand.h>
#include "predator.h"
#include "prey.h"
using namespace std;

__global__ void initialise(Prey *p, float *xrandom, float *yrandom, int nprey, float Lx, float Ly);
__global__ void prey_velocity(Prey *p, int nprey, float dt);
__device__ float prey_alignment(Prey *p, int nprey, int dimension);
__device__ float prey_attraction(Prey *p, int nprey, int dimension);
__device__ float prey_repulsion(Prey *p, int nprey, int dimension);
__device__ float prey_friction(Prey *p, int nprey, int dimension);
__global__ void prey_location(Prey *p, int nprey, float dt);
__global__ void save(Prey *p, int nprey);
__host__ void write_prey(H5PartFile *output, Prey *p, int nprey, int it);
__host__ void write_predator(H5PartFile *output, Predator *p, int it);

int main()
{
    cout << "Flockuda v1.0.0" << endl;

    // Domain
    float Lx = 1000.0;
    float Ly = 200.0;

    // Timestepping.
    float t = 0.0;
    float dt = 0.5;
    int it = 0;
    int nt = 1000;

    // Predator.
    Predator *predator;
    cudaMallocManaged(&predator, sizeof(Predator));

    // Prey.
    int nprey = 200;
    Prey *prey;
    cudaMallocManaged(&prey, nprey*sizeof(Prey));

    // Random number generator.
    curandGenerator_t generator;
    float *xrandom;
    float *yrandom;
    cudaMallocManaged(&xrandom, nprey*sizeof(float));
    cudaMallocManaged(&yrandom, nprey*sizeof(float));
    curandCreateGenerator(&generator, CURAND_RNG_PSEUDO_MTGP32);
    curandSetPseudoRandomGeneratorSeed(generator, unsigned(time(NULL)));
    curandGenerateUniform(generator, xrandom, nprey);
    curandGenerateUniform(generator, yrandom, nprey);

    // Files.
    H5PartFile *output_prey = H5PartOpenFile("prey.h5part", H5PART_WRITE);
    H5PartSetNumParticles(output_prey, nprey);
    H5PartFile *output_predator = H5PartOpenFile("predator.h5part", H5PART_WRITE);
    H5PartSetNumParticles(output_predator, 1);


    // Initialise.
    initialise<<<200, 200>>>(prey, xrandom, yrandom, nprey, Lx, Ly);
    cudaDeviceSynchronize();

    write_prey(output_prey, prey, nprey, it);
    write_predator(output_predator, predator, it);

    
    while(it < nt)
    {
        cout << it << "\t" << t << endl;

        H5PartSetStep(output_prey, it);
        H5PartSetStep(output_predator, it);

        // Compute prey velocities.
        prey_velocity<<<200, 200>>>(prey, nprey, dt);
        cudaDeviceSynchronize();

        cout << prey[0].v[0] << " " << prey[0].v[1] << endl;

        prey_location<<<200, 200>>>(prey, nprey, dt);
        cudaDeviceSynchronize();

        save<<<200, 200>>>(prey, nprey);
        cudaDeviceSynchronize();

        write_prey(output_prey, prey, nprey, it);
        write_predator(output_predator, predator, it);

        // Update time.
        it += 1;
        t += dt;
    }

    // Free unified memory.
    cudaFree(predator);
    cudaFree(prey);
    curandDestroyGenerator(generator);
    cudaFree(xrandom);
    cudaFree(yrandom);

    // Close output streams.
    H5PartCloseFile(output_prey);
    H5PartCloseFile(output_predator);

    return 0;
}

__global__ void prey_velocity(Prey *p, int nprey, float dt)
{
    float f[2];

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nprey)
    {
        for(int d=0; d<2; ++d)
        {
            // Compute force terms.
            f[d] = prey_alignment(p, nprey, d) + prey_attraction(p, nprey, d) + prey_repulsion(p, nprey, d) - prey_friction(p, nprey, d);

            // Compute velocity using F = ma.
            p[i].v[d] = dt*(p[i].vold[d] + (1.0/p[i].m)*(f[d]));
        }
    }

    return;
}

__device__ float prey_alignment(Prey *p, int nprey, int dimension)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float g = 0.5;
    float sum = 0.0;
    float magnitude = 0.0;
    for(int j=0; j<nprey; ++j)
    {
        if(i != j)
        {
            magnitude = sqrt(pow(p[i].x[0] - p[j].x[0], 2) + pow(p[i].x[1] - p[j].x[1], 2));
            if(fabs(p[j].v[dimension]) != 0)
            {
                sum += (g/magnitude)*(p[j].v[dimension]/fabs(p[j].v[dimension]));
            }
        }
    }

    return sum;
}

__device__ float prey_attraction(Prey *p, int nprey, int dimension)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float catt = 7.0;
    float latt = 100.0;
    float sum = 0.0;
    float magnitude = 0.0;
    float f = 0.0;
    for(int j=0; j<nprey; ++j)
    {
        if(i != j)
        {
            magnitude = sqrt(pow(p[i].x[0] - p[j].x[0], 2) + pow(p[i].x[1] - p[j].x[1], 2));
            f = (p[i].x[dimension] - p[j].x[dimension])/magnitude;
            sum += exp(-magnitude/latt)*f;
        }
    }

    return catt*sum;
}

__device__ float prey_repulsion(Prey *p, int nprey, int dimension)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float crep = 10.0;
    float lrep = 150.0;
    float sum = 0.0;
    float magnitude = 0.0;
    float f = 0.0;
    for(int j=0; j<nprey; ++j)
    {
        if(i != j)
        {
            magnitude = sqrt(pow(p[i].x[0] - p[j].x[0], 2) + pow(p[i].x[1] - p[j].x[1], 2));
            f = (p[i].x[dimension] - p[j].x[dimension])/magnitude;
            sum += exp(-magnitude/lrep)*f;
        }
    }

    return -crep*sum;
}

__device__ float prey_friction(Prey *p, int nprey, int dimension)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float gamma = 0.05;
    return gamma*p[i].v[dimension];
}

__global__ void prey_location(Prey *p, int nprey, float dt)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nprey)
    {
        for(int d=0; d<2; ++d)
        {
            // Compute location solving dx/dt = v.
            p[i].x[d] = p[i].xold[d] + dt*p[i].v[d];
        }
    }

    return;
}

__global__ void save(Prey *p, int nprey)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nprey)
    {
        p[i].save();
    }

    return;
}

__global__ void initialise(Prey *p, float *xrandom, float *yrandom, int nprey, float Lx, float Ly)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nprey)
    {
        p[i].initialise(10.0, Lx*xrandom[i], Ly*yrandom[i]);
    }
    return;
}

__host__ void write_prey(H5PartFile *output, Prey *p, int nprey, int it)
{
    H5PartSetStep(output, it);

    float x[nprey];
    float y[nprey];
    for(int i=0; i <= nprey; ++i)
    {
        x[i] = p[i].x[0];
        y[i] = p[i].x[1];
    }
    H5PartWriteDataFloat32(output, "PreyX", x);
    H5PartWriteDataFloat32(output, "PreyY", y);
}

__host__ void write_predator(H5PartFile *output, Predator *p, int it)
{
    float x[1] = {p->x[0]};
    float y[1] = {p->x[1]};
    
    H5PartSetStep(output, it);
    H5PartWriteDataFloat32(output, "PredatorX", x);
    H5PartWriteDataFloat32(output, "PredatorY", y);
}
