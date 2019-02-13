__device__ float prey_alignment(Prey *p, int nprey, int dimension);
__device__ float prey_attraction(Prey *p, int nprey, int dimension);
__device__ float prey_repulsion(Prey *p, int nprey, int dimension);
__device__ float prey_friction(Prey *p, int nprey, int dimension);
__device__ float prey_avoid(Prey *p, int nprey, float xp0, float xp1, int dimension);
