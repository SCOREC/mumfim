#include "FiberNetwork.h"
#include "InitGuess.h"
#include "apfUtil.h"

#include <apfShape.h>
#include <apfMesh.h>

#include <cmath>
#include <numeric>

namespace bio
{
  FiberNetwork::FiberNetwork(apf::Mesh * f)
    : fn(f)
    , fn_u(NULL)
    , dim(f->getDimension())
  { }
}
