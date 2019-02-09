#include "predator.h"
#include "forces.h"

__host__ void Predator::initialise(float kappa, float x0, float x1, float v0, float v1)
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

__host__ void Predator::save()
{   
    // Save information from the previous timestep.
    xold[0] = x[0];
    xold[1] = x[1];
    vold[0] = v[0];
    vold[1] = v[1];

    return;
}

__host__ void initialise_predator(Predator *p)
{
    float kappa = 0.5;
    p->initialise(kappa, 0, 0, 0.5, 0.5);
    return;
}

__host__ void write_predator(H5PartFile *output, Predator *p, int it)
{
    float x[1] = {p->x[0]};
    float y[1] = {p->x[1]};
    
    H5PartSetStep(output, it);
    H5PartWriteDataFloat32(output, "PredatorX", x);
    H5PartWriteDataFloat32(output, "PredatorY", y);
    return;
}

__host__ void save_predator(Predator *p)
{
    p->save();
    return;
}
