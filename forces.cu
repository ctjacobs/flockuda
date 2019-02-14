/*

Flockuda: A numerical model of predator-prey dynamics based on a Molecular Dynamics approach.

Copyright (C) 2019 Christian Thomas Jacobs

*/

#include "predator.h"
#include "prey.h"

__device__ float prey_alignment(Prey *p, int nprey, int dimension)
{
    // The force of alignment based on the inverse square of distance between individual prey.

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float g = 0.5;
    float sum = 0.0;
    float magnitude_ij;
    float magnitude_j;
    for(int j=0; j<nprey; ++j)
    {
        if(i != j)
        {
            magnitude_ij = sqrt(pow(p[i].xold[0] - p[j].xold[0], 2) + pow(p[i].xold[1] - p[j].xold[1], 2));
            magnitude_j = sqrt(pow(p[j].vold[0], 2) + pow(p[j].vold[1], 2));
            if(magnitude_j > 1e-7 && magnitude_ij > 1e-7)
            {
                sum += (g/pow(magnitude_ij, 2))*(p[j].vold[dimension]/magnitude_j);
            }
        }
    }

    return sum;
}

__device__ float prey_attraction(Prey *p, int nprey, int dimension)
{
    // The attracting force between individual prey. This term is used to retain flock formation.

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float catt = 0.2;
    float latt = 100.0;  // Range of attractive forces.
    float sum = 0.0;
    float magnitude = 0.0;
    float fhat;  // Unit distance vector.
    for(int j=0; j<nprey; ++j)
    {
        if(i != j)
        {
            magnitude = sqrt(pow(p[i].xold[0] - p[j].xold[0], 2) + pow(p[i].xold[1] - p[j].xold[1], 2));
            fhat = (p[i].xold[dimension] - p[j].xold[dimension])/magnitude;
            if(magnitude > 1e-7)
                sum += exp(-magnitude/latt)*fhat;
        }
    }

    return catt*sum;
}

__device__ float prey_repulsion(Prey *p, int nprey, int dimension)
{
    // The repelling force between individual prey. This term is used to avoid collisions between prey.

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float crep = 0.5;
    float lrep = 150.0;  // Range of repulsive forces.
    float sum = 0.0;
    float magnitude = 0.0;
    float fhat = 0.0;  // Unit distance vector.
    for(int j=0; j<nprey; ++j)
    {
        if(i != j)
        {
            magnitude = sqrt(pow(p[i].xold[0] - p[j].xold[0], 2) + pow(p[i].xold[1] - p[j].xold[1], 2));
            fhat = (p[i].xold[dimension] - p[j].xold[dimension])/magnitude;
            if(magnitude > 1e-7)
                sum += exp(-magnitude/lrep)*fhat;
        }
    }

    return -crep*sum;
}

__device__ float prey_friction(Prey *p, int nprey, int dimension)
{
    // A frictional drag force used to prevent the prey from moving too quickly.

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float gamma = 0.05;
    return gamma*p[i].vold[dimension];
}

__device__ float prey_avoid(Prey *p, int nprey, float xp0, float xp1, int dimension)
{
    // A force in response to predator attack. This term allows the prey to avoid the predator by moving away from it.

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float cavoid = 50;
    float omega = 0.2;  // Scaling factor based on predator size.
    float R = 100.0;  // Risk range.
    float xp[2] = {xp0, xp1};  // Location of the predator.
    float magnitude = sqrt(pow(xp[0] - p[i].xold[0], 2) + pow(xp[1] - p[i].xold[1], 2));
    return cavoid*(1.0/(1.0 + exp(omega*(magnitude - R))))*(xp[dimension] - p[i].xold[dimension])/magnitude;
}
