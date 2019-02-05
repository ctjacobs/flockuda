#include <stdlib.h>
#include "prey.h"

__device__ void Prey::initialise(double mass, double x0, double x1)
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
