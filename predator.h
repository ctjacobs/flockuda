/*

Flockuda: A numerical model of predator-prey dynamics based on a Molecular Dynamics approach.

Copyright (C) 2019 Christian Thomas Jacobs

*/

#include <H5Part.h>

class Predator
{
    public:

        float k;  // Strength of bias (i.e. speed of attack).

        float x[2];  // Location of the predator in two-dimensional space.
        float xold[2];  // Location of the predator in two-dimensional space at the previous timestep.

        float v[2];  // Velocity of the predator.
        float vold[2];  // Velocity of the predator at the previous timestep.

        __host__ void initialise(float kappa, float x0, float x1, float v0, float v1);
        __host__ void save();
};

__host__ void initialise_predator(Predator *p);
__host__ void write_predator(H5PartFile *output, Predator *p, int it);
__host__ void save_predator(Predator *p);
__host__ void predator_velocity(Predator *p, float *c, float *xrandom, float dt);
__host__ void predator_location(Predator *p, float dt);
