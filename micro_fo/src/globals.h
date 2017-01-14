#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <ctype.h>
#include <cfloat>
#include <strings.h>
#include <time.h>

/* Dimensions of the RVE. The RVE is a cube of unit length. 
  The following six variables provide the distance of the 
  sides of the cube from its center */

// These values are used if not specified in file
#define DEFINE_CUBE_BOTTOM -0.5
#define DEFINE_CUBE_TOP     0.5
#define DEFINE_CUBE_LEFT   -0.5
#define DEFINE_CUBE_RIGHT   0.5
#define DEFINE_CUBE_BACK   -0.5
#define DEFINE_CUBE_FRONT   0.5

// Periodic boundary conditions (network files must also be periodic format)
// If non-periodic, then all dofs on boundary are dirichlet, 0 or specified displacement
#define PERIODIC false

// Timoshenko or Euler-Bernoulli beams
#define TIMOSHENKO false

// Include nonlinear geometric stiffness matrix
#define GEOMETRIC_STIFFNESS false

// true for truss - 3 dofs per node, axial forces only
// false for beam - 6 dofs per node, axial, bending, and torsion
#define TRUSS true

// true if there is additional column in RVE file that specifies fiber type.
// used to distinguish between fibers with different properties.
#define SPECIFY_FIBER_TYPE true

// If using truss
// true for nonlinear force: f = E A / B (exp(B/2 (lambda^2 - 1)) - 1)
// false for linear force: f = E A ( lambda - 1 )
// where lambda = current length / original length
#define NONLINEAR_FORCE true

// True to include several segments of code helpful for size effects tests
#define FIBER_ONLY_SIZE_EFFECT_TEST false

