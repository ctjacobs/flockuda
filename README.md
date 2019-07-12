# Flockuda

Flockuda is a predator-prey model that runs on a CUDA-enabled NVIDIA Graphics Processing Unit (GPU). The model is based on the Molecular Dynamics approach described by Lee et al. (2006), "Prey-Flock Deformation under a Predator's Attack", Journal of the Korean Physical Society, 48:S236--S240. However, no stochastic forces are currently included. Furthermore, the model implements a direct approach to force calculation, such that the computational complexity is O(N^2) where N is the number of prey.

## Dependencies

The [NVIDIA CUDA C compiler](https://developer.nvidia.com/cuda-toolkit) is required in order to compile Flockuda. The [H5Part library](https://code.lbl.gov/projects/h5part/) is also required, since the prey and predator solution fields are written to an .h5part file. Flockuda is known to successfully compile with v9.1.85 of the NVIDIA CUDA C compiler and v1.6.6 of H5Part. Users may need to adjust the paths to the H5Part library in the `Makefile`.

## Usage
A simulation can be set up by editing the parameters in `example.cfg`, then compiling and executing the code by running the following commands at the command line:

```
make
./flockuda example.cfg
```

## Authors

* [Christian T. Jacobs](http://christianjacobs.uk/)

## Copyright statement

Copyright (C) 2019 Christian T. Jacobs
