/*

Flockuda: A numerical model of predator-prey dynamics based on a Molecular Dynamics approach.

Copyright (C) 2019 Christian Thomas Jacobs

*/

#include <iostream>
#include <H5Part.h>
#include <cuda.h>
#include <curand.h>
#include "predator.h"
#include "prey.h"
using namespace std;


int main()
{
    // All equations solved within are based on those described by Lee et al. (2006), "Prey-Flock Deformation under a Predator's Attack", Journal of the Korean Physical Society, 48:S236--S240.

    cout << "Flockuda v1.0.0" << endl;

    // Domain specification.
    float Lx = 1000.0;
    float Ly = 1000.0;

    // Timestepping parameters.
    float t = 0.0;
    float dt = 0.2;
    int it = 0;
    int nt = 1000;

    // Predator.
    Predator *predator;
    cudaMallocManaged(&predator, sizeof(Predator));

    // Prey.
    int nprey = 200;
    Prey *prey;
    cudaMallocManaged(&prey, nprey*sizeof(Prey));

    // Centre of flock.
    float centre[2];

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

    // Output file streams.
    H5PartFile *output_prey = H5PartOpenFile("prey.h5part", H5PART_WRITE);
    H5PartSetNumParticles(output_prey, nprey);
    H5PartFile *output_predator = H5PartOpenFile("predator.h5part", H5PART_WRITE);
    H5PartSetNumParticles(output_predator, 1);

    // Initialise.
    initialise_prey<<<1, nprey>>>(prey, xrandom, yrandom, nprey, Lx, Ly);
    cudaDeviceSynchronize();
    initialise_predator(predator);
    cudaDeviceSynchronize();

    // Write initial condition.
    write_prey(output_prey, prey, nprey, it);
    write_predator(output_predator, predator, it);

    // Timestepping loop.
    while(it < nt)
    {
        cout << "Iteration " << it << "\t Time: " << t << endl;        

        // Compute the centre of the flock.
        prey_centre(prey, nprey, centre);

        // Compute predator velocity.
        predator_velocity(predator, centre, xrandom, dt);
        predator_location(predator, dt);
        save_predator(predator);

        // Compute prey velocities.
        prey_velocity<<<1, nprey>>>(prey, nprey, predator->x[0], predator->x[1], dt);
        cudaDeviceSynchronize();
        prey_location<<<1, nprey>>>(prey, nprey, dt);
        cudaDeviceSynchronize();
        save_prey<<<1, nprey>>>(prey, nprey);
        cudaDeviceSynchronize();

        // Write prey and predator positions to file.
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

