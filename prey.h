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
