class Prey
{
    public:

        double m;  // Mass of the prey.

        double x[2];  // Location of the prey in two-dimensional space.
        double xold[2];  // Location of the prey in two-dimensional space at the previous timestep.

        double v[2];  // Velocity of the prey.
        double vold[2];  // Velocity of the prey at the previous timestep.

        __device__ void initialise(double mass, double x0, double x1);
        __device__ void save();
};
