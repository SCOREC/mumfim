#ifndef BIO_UTIL_IMPL_H_
#define BIO_UTIL_IMPL_H_
namespace bio
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
  template <typename I, typename O>
    void mapGlobalsToLocals(apf::Mesh * msh,
                            apf::MeshEntity * e,
                            I gbl_crds_bgn,
                            I gbl_crds_end,
                            O lcl_crds)
  {
    for(auto crd = gbl_crds_bgn; crd != gbl_crds_end; ++crd)
    {
      apf::Vector3 lcl;
      mapGlobalToLocal(msh,e,*crd,lcl);
      *lcl_crds++ = lcl;
    }
  }
}
#endif
