#include <simClassified.h>
#include <apf.h>
#include <cassert>
namespace bio
{
  template <typename O>
    void getFieldComponentsClassified(apf::Field * fld,
                                      apf::ModelEntity * mdl_ent,
                                      int nm_cmps,
                                      int * cmps,
                                      O out)
  {
    for(int ent_dim = 0; ent_dim < 3; ent_dim++)
    {
      amsi::SimClassifiedIter bgn = amsi::beginClassified(apf::getMesh(fld),
                                                          mdl_ent,
                                                          ent_dim);
      amsi::SimClassifiedIter end = amsi::endClassified(bgn);
      getFieldComponentsOn(fld,nm_cmps,cmps,bgn,end,out);
    }
  }
  template <typename I, typename O>
    void getFieldComponentsOn(apf::Field * fld,
                              int nm_cmps,
                              int * cmps,
                              I bgn,
                              I end,
                              O out)
  {
    for(I it = bgn; it != end; it++)
      getEntFieldComponents(fld,*it,nm_cmps,cmps,out);
  }
  template <typename O>
    void getFieldComponentsOnEnt(apf::Field * fld,
                                 apf::MeshEntity * ent,
                                 int nm_cmps,
                                 int * cmps,
                                 O out)
  {
    apf::FieldShape * fs = apf::getShape(fld);
    int cnt = fs->countNodesOn(apf::getMesh(fld)->getType(ent));
    for(int nde = 0; nde < cnt; nde++)
    {
      int fld_cmps = apf::countComponents(fld);
      assert(nm_cmps < fld_cmps);
      double fld_data[fld_cmps];
      apf::getComponents(fld,ent,nde,fld_data);
      for(int ii = 0; ii < nm_cmps; ++ii)
      {
        assert(cmps[ii] < fld_cmps);
        *out++ = fld_data[cmps[ii]];
      }
    }
  }
}
