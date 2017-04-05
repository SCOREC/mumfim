#ifndef BIO_UTIL_H_
#define BIO_UTIL_H_
#include <apf.h>
namespace bio
{
  void mapGlobalToLocal(apf::Mesh * msh,
                        apf::MeshEntity * e,
                        apf::Vector3 const& global,
                        apf::Vector3 & local);
}
#endif
