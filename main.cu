/*

Flockuda: A numerical model of predator-prey dynamics based on the Molecular Dynamics approach of Lee et al. (2006).

Copyright (C) 2019 Christian T. Jacobs

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/


#include <iostream>
#include <fstream>
#include <H5Part.h>
#include <cuda.h>
#include <curand.h>
#include <cassert>
#include "configuration.h"
#include "predator.h"
#include "prey.h"
using namespace std;


int main(int argc, char **argv)
{
    // All equations solved within are based on those described by Lee et al. (2006), "Prey-Flock Deformation under a Predator's Attack", Journal of the Korean Physical Society, 48:S236--S240.

    cout << "Flockuda v1.0.0" << endl;

    // Read in the simulation's configuration from a file.
    assert(argc == 2);
    Configuration config;
    #if DEBUG
        config.read("tests/test.cfg");
        assert(config.Lx == 10.0);
        assert(config.Ly == 20.0);
        assert(config.mass == 1.0);
        assert(config.nprey == 200);
        assert(config.R == 100.0);
    #else
        config.read(argv[1]);
    #endif

    // Initialise timestepping variables.
    float t = 0.0;
    int it = 0;

    // Predator.
    Predator *predator;
    cudaMallocManaged(&predator, sizeof(Predator));

    // Prey.
    Prey *prey;
    cudaMallocManaged(&prey, config.nprey*sizeof(Prey));

    // Centre of flock.
    float centre[2];

    // Random number generator.
    curandGenerator_t generator;
    float *xrandom;
    float *yrandom;
    cudaMallocManaged(&xrandom, config.nprey*sizeof(float));
    cudaMallocManaged(&yrandom, config.nprey*sizeof(float));
    curandCreateGenerator(&generator, CURAND_RNG_PSEUDO_MTGP32);
    curandSetPseudoRandomGeneratorSeed(generator, unsigned(time(NULL)));
    curandGenerateUniform(generator, xrandom, config.nprey);
    curandGenerateUniform(generator, yrandom, config.nprey);

    // Output file streams in .h5part (HDF5 Particle) format.
    H5PartFile *output_prey = H5PartOpenFile("prey.h5part", H5PART_WRITE);
    H5PartSetNumParticles(output_prey, config.nprey);
    H5PartFile *output_predator = H5PartOpenFile("predator.h5part", H5PART_WRITE);
    H5PartSetNumParticles(output_predator, 1);

    // Initialise.
    cudaDeviceSynchronize();
    initialise_prey(prey, config, xrandom, yrandom);
    cudaDeviceSynchronize();
    initialise_predator(predator, config);
    cudaDeviceSynchronize();

    // Write initial condition.
    write_prey(output_prey, prey, config, it);
    write_predator(output_predator, predator, it);

    // Timestepping loop.
    while(it < config.nt)
    {
        cout << "Iteration " << it << "\t Time: " << t << endl;        

        // Compute the centre of the flock.
        prey_centre(prey, config, centre);

        // Compute predator velocity.
        predator_velocity(predator, config, centre, xrandom);
        predator_location(predator, config);
        save_predator(predator);

        // Compute prey velocities on the CUDA-enabled graphics processing unit (GPU).
        prey_velocity<<<1, config.nprey>>>(prey, config, predator->x[0], predator->x[1]);
        cudaDeviceSynchronize();
        prey_location<<<1, config.nprey>>>(prey, config);
        cudaDeviceSynchronize();
        save_prey(prey, config);
        cudaDeviceSynchronize();

        // Write prey and predator positions to file.
        write_prey(output_prey, prey, config, it);
        write_predator(output_predator, predator, it);

        // Update time.
        it += 1;
        t += config.dt;
    }

    #if DEBUG
        // Check that the predator's position has not changed.
        assert(predator->x[0] == 0);
        assert(predator->x[1] == 0);

        // Check finish time and iteration count.
        assert(it == 1);
        assert(fabs(t - 0.2) < 1e-6);

        // Check output files exist.
        ifstream f;
        f.open("predator.h5part");
        assert(f.good());
        f.close();
        f.open("prey.h5part");
        assert(f.good());
        f.close();
    #endif

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

