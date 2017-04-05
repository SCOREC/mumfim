#include "bioUtil.h"
namespace bio
{
  void mapGlobalToLocal(apf::Mesh * msh,
                        apf::MeshEntity * ent,
                        apf::Vector3 const& global,
                        apf::Vector3 & local)
  {
    ma::Affine i = ma::getMap(msh,ent));
    local = i * global;
  }
}
