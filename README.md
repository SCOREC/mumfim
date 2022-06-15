# MuMFiM
multiscale biological tissue implementation using amsi
![build](https://github.com/SCOREC/mumfim/actions/workflows/run_tests.yml/badge.svg?branch=develop)

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

# Microscale Solver
This solver is designed to be a truss/beam solver which can operate on its own, or as part
of a biotissue multiscale problem.

## Current Limitations
- Only Truss Networks are implemented
- Implicit only solver
- Uses single precicion linear algebra solver

## Input Files

### Network Input File
The network input file describes the geometry of the fiber network to be read in.
It has the following description:
```
num-nodes num-fibers 0
node-1-x node-1-y node-1-z
...
node-n-x node-n-y node-n-z
fiber-1-node-1 fiber-1-node-2
...
fiber-m-node-1 fiber-m-node-2
```
The node ids use zero based indexing. See [here](../test/microscale/fiber_networks/del_4450seedL5_new_1.txt) for an example.

### Parameter Input File
The parameter file describes the material parameters of each of the fibers. It should be named the
same as the network input file with `.params` appended to the end. E.g. if the network input file
is name `test_net.txt` the params file should be named `test_net.txt.params`.

```
num-parameters num-fibers
reaction-id  fiber-radius  youngs-modulus ... density
fiber-1-parameter
...
fiber-m-parameter
```
The `reaction-id` selects the type of reaction you would like to use. Currently there are the following options and parameters:

| Reaction Type | `reaction-id` | Parameters                                                                 |
|---------------|:-------------:|----------------------------------------------------------------------------|
| Linear        | 0             | `fiber-radius`, `youngs-modulus`, `density`                                |
| Linear        | 0             | `fiber-radius`, `youngs-modulus`, `B`, `transition-length-ratio`, `density`|

See [here](../test/microscale/fiber_networks/del_4450seedL5_new_1.txt.params) for an example.

The current version of this input file is limited to linear truss elements. A future version will
allow the user to select the element type and enter the appropriate number of parameters as needed.

### Microscale only input file
The microscale only input file is a relatively self explanatory yaml input file. See
[here](../test/microscale/fiber_only.yaml) for an example.

# Preparing the Model
1. Setup model and mesh in simmodeler. Attribute definition template in Simmodeler called `macro`. 
2. Convert mesh to pumi mesh (requires pumi +simmetrix)
   ```console
   $ convert SimmetrixModel.smd SimmetrixMesh.sms OutPumiMesh.smb [--native-model=native_model.x_t]
   ```
3. Convert model to PUMI model (requires pumi +simmetrix)
   ```console
   $ mdlConvert SimmetrixModel.smd PumiModel.dmg  
   ```
4. Convert model atributes to yaml (requires model-traits +simmetrix)
   ```console
   $ smd2yaml SimmetrixModel.smd > ModelTraitsModel.yaml
   ```
