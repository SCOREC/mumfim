#ifndef MUMFIM_FIBER_NETWORK_IO_H_
#define MUMFIM_FIBER_NETWORK_IO_H_
#include <apf.h>
#include <apfMesh2.h>
#include <iostream>
#include "FiberReactions.h"
namespace mumfim
{
  apf::Mesh2 * loadFromStream(std::istream & strm);
  template <typename O>
  void loadParamsFromStream(apf::Mesh2 * msh, std::istream & strm, O rctns);
}  // namespace mumfim
#include "FiberNetworkIO_impl.h"
#endif
