#include <bioApfPointers.h>
#include <apfMesh.h>
#include <apf.h>

namespace mumfim
{
  mesh_unique_ptr_type make_unique(apf::Mesh * mesh)
  {
    return mesh_unique_ptr_type{mesh, &deleteMesh};
  }
  void deleteMesh(apf::Mesh * mesh)
  {
    mesh->destroyNative();
    apf::destroyMesh(mesh);
  }
}
