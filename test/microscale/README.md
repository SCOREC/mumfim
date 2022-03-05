# Tests

This readme gives some information on what the tests do and how to prepare them
for adding new data.

## List of tests

- test_global_stiffness

    This test verifies that the global stiffness matrix matches with a "clean room"
    implementation written in python. Currently there are two networks that are tested.
    If you would like to add more.
    
    If you would like to clean the networks before use you can use the [fix_broken_networks.py](https://gitlab.com/jacob.merson/fiber_netgen/blob/master/examples/fix_broken_networks.py).
    script which is part of the [fiber_netgen](https://gitlab.com/jacob.merson/fiber_netgen) package.

    Generating the stiffness matrix can be done using the [compute_global_truss.py](./compute_global_truss.py) script. This script will also give you some statistics
    about your stiffness matrix to tell you if the geometry is reasonable or not.
    If you have more than 6 zero/negative eigenvalues, then there is a problem because this
    represents unconstrained modes in the system. If this is the case try using the network cleanup utility.

- test_elemental stiffness
    This test computest the elemental stiffness matrices for a given mesh. Currently no automated
    tests use this exe, but it is useful for debugging purposes.

## Relavent test scripts

- [compute_global_truss.py](./compute_global_truss.py)
    This scripts generates a global stiffness matrix from a biotissue input network.
- [read_abaqus_mtx.py](./read_abaqus_mtx.py)
    This script reads an abaqus stiffness matrix and performs the validity checks on that matrix.
    This is somewhat useful for debugging.
- [stiffness_validity.py](./stiffness_validity.py)
    This gives a utility for computing some stiffness matrix validity. E.g. Positive definiteness.
