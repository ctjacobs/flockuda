#include <stdlib.h>
#include "predator.h"

__device__ void Predator::initialise(float kappa, float x0, float x1, float v0, float v1)
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

__device__ void Predator::save()
{   
    // Save information from the previous timestep.
    xold[0] = x[0];
    xold[1] = x[1];
    vold[0] = v[0];
    vold[1] = v[1];

    return;
}
