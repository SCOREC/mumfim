#ifndef MUMFIM_UTIL_IMPL_H_
#define MUMFIM_UTIL_IMPL_H_
namespace mumfim
{
  template <typename O>
    void calcDimMeasures(apf::Mesh * msh, int dim, O msrs)
  {
    apf::MeshIterator * it = NULL;
    apf::MeshEntity * ent = NULL;
    for(it = msh->begin(dim); (ent = msh->iterate(it));)
    {
      apf::MeshElement * elm = apf::createMeshElement(msh,ent);
      *msrs++ = apf::measure(elm);
      apf::destroyMeshElement(elm);
    }
    msh->end(it);
  }
}
#endif
