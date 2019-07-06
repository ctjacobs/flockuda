/*

Flockuda: A numerical model of predator-prey dynamics based on the Molecular Dynamics approach of Lee et al. (2006).

Copyright (C) 2019 Christian T. Jacobs

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

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
__host__ void predator_location(Predator *p, float dt, float Lx, float Ly);
