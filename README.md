# Biotissue
multiscale biological tissue implementation using amsi

## Installation
The easiest way to build mumfim is through [spack](https://github.com/spack/spack). See the
instructions [here](https://github.com/jacobmerson/mumfim-spack)

## Notes on Problem Setup
1. In simmodeler you can use `val*$t` in the dirichlet or neumann boundary condition. This will
   apply `val` load linearly accross the number of timesteps you selected. E.g. If you apply a
   displacement `1*$t`, the final displacement will be `1`, and if you have `20` timesteps,
   each step will apply `0.05`.
2. Currently you cannot output the load on a face with a dirichlet boundary condition
   it will cause a crash.
