/* Flockuda

Copyright (C) 2019 Christian Thomas Jacobs

*/

#include <iostream>
#include <H5Part.h>
#include <cuda.h>
#include <curand.h>
#include "prey.h"
using namespace std;

__global__ void initialise(Prey *p, float *xrandom, float *yrandom, int nprey, float Lx, float Ly);
__global__ void prey_velocity(Prey *p, int nprey, float dt);
__global__ void prey_location(Prey *p, int nprey, float dt);
__global__ void save(Prey *p, int nprey);
__host__ void write(H5PartFile *output, Prey *p, int nprey, int it);

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
    int nt = 5;

    // Prey.
    int nprey = 2048;
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

    // File.
    H5PartFile *output = H5PartOpenFile("prey.h5part", H5PART_WRITE);
    H5PartSetNumParticles(output, nprey);

    // Initialise.
    initialise<<<1024, 1024>>>(prey, xrandom, yrandom, nprey, Lx, Ly);
    cudaDeviceSynchronize();

    write(output, prey, nprey, it);

    while(it < nt)
    {
        cout << it << "\t" << t << endl;

        H5PartSetStep(output, it);

        // Compute prey velocities.
        prey_velocity<<<1024, 1024>>>(prey, nprey, dt);
        cudaDeviceSynchronize();

        cout << prey[0].x[0] << " " << prey[0].x[1] << endl;

        prey_location<<<1024, 1024>>>(prey, nprey, dt);
        cudaDeviceSynchronize();

        save<<<1024, 1024>>>(prey, nprey);
        cudaDeviceSynchronize();

        write(output, prey, nprey, it);

        // Update time.
        it += 1;
        t += dt;
    }

    // Free unified memory.
    cudaFree(prey);
    curandDestroyGenerator(generator);
    cudaFree(xrandom);
    cudaFree(yrandom);

    // Close output stream.
    H5PartCloseFile(output);

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
            f[d] = 100.0;

            // Compute velocity using F = ma.
            p[i].v[d] = (1.0/dt)*p[i].vold[d] + (1.0/p[i].m)*(f[d]);
        }
    }

    return;
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

__host__ void write(H5PartFile *output, Prey *p, int nprey, int it)
{
    H5PartSetStep(output, it);

    float x[nprey];
    float y[nprey];
    for(int i=0; i <= nprey; ++i)
    {
        x[i] = p[i].x[0];
        y[i] = p[i].x[1];
    }
    H5PartWriteDataFloat32(output, "x", x);
    H5PartWriteDataFloat32(output, "y", y);
}
