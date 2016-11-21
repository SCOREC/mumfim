#include "bioVariableRecovery.h"
#include <SimFEA.h>
#include <apfField.h>
#include <apfSIM.h>
#include <cassert>
#include <list>
#include <numeric>
namespace bio
{
  void getEntFieldComponent(apf::Field * fld, int crd, std::vector<double>& crds, pEntity ent)
  {
    apf::MeshEntity * apf_ent = apf::castEntity(ent);
    int cnt = fld->countNodesOn(apf_ent);
    for(int nde = 0; nde < cnt; nde++)
    {
      int cmpnts = apf::countComponents(fld);
      assert(crd < cmpnts);
      double nd_crds[cmpnts];
      apf::getComponents(fld,apf_ent,nde,nd_crds);
      crds.push_back(nd_crds[crd]);
    }
  }
  void getFieldComponentOn(apf::Field * fld, int crd, std::vector<double>& crds, pMesh msh, pGEntity mdl_ent)
  {
    for(int ent_dim = 0; ent_dim < 3; ent_dim++)
    {
      std::list<pEntity> msh_ents;
      amsi::Model_GetClassifiedEntities(msh,mdl_ent,ent_dim,msh_ents);
      getEntsFieldComponent(fld,crd,crds,msh_ents.begin(),msh_ents.end());
    }
  }
}
