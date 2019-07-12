/*

Flockuda: A numerical model of predator-prey dynamics based on the Molecular Dynamics approach of Lee et al. (2006).

Copyright (C) 2019 Christian T. Jacobs

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/


#include "predator.h"
#include "prey.h"


__device__ float prey_alignment(Prey *p, Configuration config, int dimension)
{
    // The force of alignment based on the inverse square of distance between individual prey.

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float sum = 0.0;
    float magnitude_ij;
    float magnitude_j;
    for(int j=0; j<config.nprey; ++j)
    {
        if(i != j)
        {
            magnitude_ij = sqrt(pow(p[i].xold[0] - p[j].xold[0], 2) + pow(p[i].xold[1] - p[j].xold[1], 2));
            magnitude_j = sqrt(pow(p[j].vold[0], 2) + pow(p[j].vold[1], 2));
            if(magnitude_j > 1e-7 && magnitude_ij > 1e-7)
            {
                sum += (config.g/pow(magnitude_ij, 2))*(p[j].vold[dimension]/magnitude_j);
            }
        }
    }

    return sum;
}

__device__ float prey_attraction(Prey *p, Configuration config, int dimension)
{
    // The attracting force between individual prey. This term is used to retain flock formation.

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float sum = 0.0;
    float magnitude = 0.0;
    float fhat;  // Unit distance vector.
    for(int j=0; j<config.nprey; ++j)
    {
        if(i != j)
        {
            magnitude = sqrt(pow(p[i].xold[0] - p[j].xold[0], 2) + pow(p[i].xold[1] - p[j].xold[1], 2));
            fhat = (p[i].xold[dimension] - p[j].xold[dimension])/magnitude;
            if(magnitude > 1e-7)
                sum += exp(-magnitude/config.latt)*fhat;
        }
    }

    return config.catt*sum;
}

__device__ float prey_repulsion(Prey *p, Configuration config, int dimension)
{
    // The repelling force between individual prey. This term is used to avoid collisions between prey.

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float sum = 0.0;
    float magnitude = 0.0;
    float fhat = 0.0;  // Unit distance vector.
    for(int j=0; j<config.nprey; ++j)
    {
        if(i != j)
        {
            magnitude = sqrt(pow(p[i].xold[0] - p[j].xold[0], 2) + pow(p[i].xold[1] - p[j].xold[1], 2));
            fhat = (p[i].xold[dimension] - p[j].xold[dimension])/magnitude;
            if(magnitude > 1e-7)
                sum += exp(-magnitude/config.lrep)*fhat;
        }
    }

    return -config.crep*sum;
}

__device__ float prey_friction(Prey *p, Configuration config, int dimension)
{
    // A frictional drag force used to prevent the prey from moving too quickly.

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    return config.gamma*p[i].vold[dimension];
}

__device__ float prey_avoid(Prey *p, Configuration config, float xp0, float xp1, int dimension)
{
    // A force in response to predator attack. This term allows the prey to avoid the predator by moving away from it.

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float xp[2] = {xp0, xp1};  // Location of the predator.
    float magnitude = sqrt(pow(xp[0] - p[i].xold[0], 2) + pow(xp[1] - p[i].xold[1], 2));
    return config.cavoid*(1.0/(1.0 + exp(config.omega*(magnitude - config.R))))*(xp[dimension] - p[i].xold[dimension])/magnitude;
}
