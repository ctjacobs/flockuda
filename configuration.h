/*

Flockuda: A numerical model of predator-prey dynamics based on the Molecular Dynamics approach of Lee et al. (2006).

Copyright (C) 2019 Christian T. Jacobs

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/


#pragma once


class Configuration
{
    public:

        // Domain specification.
        float Lx;
        float Ly;

        // Timestepping parameters.
        float dt;
        int nt;

        // Number of prey.
        int nprey;

        // The mass of each prey.
        float mass;

        // Strength of bias (i.e. speed of attack).
        float kappa;

        float g;

        // Attraction coefficients.
        float catt;
        float latt;  // Range of attractive forces.

        // Repulsion coefficients.
        float crep;
        float lrep;  // Range of repulsive forces.

        float gamma;

        float cavoid;
        float omega;  // Scaling factor based on predator size.
        float R;  // Risk radius.

        __host__ void read(char *path);
};
