#include <H5Part.h>

class Prey
{
    public:

        float m;  // Mass of the prey.

        float x[2];  // Location of the prey in two-dimensional space.
        float xold[2];  // Location of the prey in two-dimensional space at the previous timestep.

        float v[2];  // Velocity of the prey.
        float vold[2];  // Velocity of the prey at the previous timestep.

        __device__ void initialise(float mass, float x0, float x1);
        __device__ void save();
};

__global__ void initialise_prey(Prey *p, float *xrandom, float *yrandom, int nprey, float Lx, float Ly);
__host__ void write_prey(H5PartFile *output, Prey *p, int nprey, int it);
__global__ void save_prey(Prey *p, int nprey);
__global__ void prey_velocity(Prey *p, int nprey, float dt);
__global__ void prey_location(Prey *p, int nprey, float dt);
__host__ void centre(Prey *p, int nprey, float *c);
