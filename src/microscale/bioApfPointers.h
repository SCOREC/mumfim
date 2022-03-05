#ifndef BIO_APF_CUSTOM_DELETERS_H__
#define BIO_APF_CUSTOM_DELETERS_H__
#include <memory>
namespace apf {
  class Mesh;
}
namespace bio
{
  // auto deleteMesh = [](apf::Mesh* mesh)
  //{
  // mesh->destroyNative();
  // apf::destroyMesh(mesh);
  //};
  // we we can modify this to use a lambda with c++14? because we can
  // use auto type deduction. This is desireable since it is costly to
  // use a function pointer rather than a lambda
  void deleteMesh(apf::Mesh * mesh);
  using mesh_unique_ptr_type =
      std::unique_ptr<apf::Mesh, decltype(&deleteMesh)>;
  mesh_unique_ptr_type make_unique(apf::Mesh *);
}  // namespace bio
#endif
