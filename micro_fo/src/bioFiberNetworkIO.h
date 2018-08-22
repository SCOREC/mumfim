#ifndef BIO_FIBER_NETWORK_IO_H_
#define BIO_FIBER_NETWORK_IO_H_
#include <apf.h>
#include <apfMesh2.h>
#include <iostream>
#include "bioFiberReactions.h"
namespace bio
{
  apf::Mesh2 * loadFromStream(std::istream & strm);
  apf::Mesh2 * loadFromFile(const std::string & fnm);
  template <typename O>
  void loadParamsFromFile(apf::Mesh2 * msh, const std::string & fnm, O rctns);
}  // namespace bio
#include "bioFiberNetworkIO_impl.h"
#endif
