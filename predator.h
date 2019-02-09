class Predator
{
    public:

        float k;  // Strength of bias.

        float x[2];  // Location of the predator in two-dimensional space.
        float xold[2];  // Location of the predator in two-dimensional space at the previous timestep.

        float v[2];  // Velocity of the predator.
        float vold[2];  // Velocity of the predator at the previous timestep.

        __device__ void initialise(float kappa, float x0, float x1, float v0, float v1);
        __device__ void save();
};
