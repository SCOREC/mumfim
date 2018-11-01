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
The node ids use zero based indexing. See [here](./test/fiber_networks/del_4450seedL5_new_1.txt) for an example.

### Parameter Input File
The parameter file describes the material parameters of each of the fibers. It should be named the
same as the network input file with `.params` appended to the end. E.g. if the network input file
is name `test_net.txt` the params file should be named `test_net.txt.params`.

```
num-parameters num-fibers
param-id  fiber-radius  youngs-modulus fiber-density
fiber-1-parameter
...
fiber-m-parameter
```
The parameter id must start at zero and be monotonically increasing.
See [here](./test/fiber_networks/del_4450seedL5_new_1.txt.params) for an example.

The current version of this input file is limited to linear truss elements. A future version will
allow the user to select the element type and enter the appropriate number of parameters as needed.

### Microscale only input file
The microscale only input file is a relatively self explanatory yaml input file. See 
[here](./test/fiber_only.yaml) for an example.
