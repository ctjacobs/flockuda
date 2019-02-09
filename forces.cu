#include "predator.h"
#include "prey.h"

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

