#ifndef BIO_VARIABLE_RECOVERY_H_
#define BIO_VARIABLE_RECOVERY_H_
#include <apf.h>
#include <MeshSim.h>
#include <vector>
namespace bio
{
  // calclate the average value of field coordinate cord for all nodes on all entities classified on mdl_ent
  void getFieldComponentOn(apf::Field * fld, int crd, std::vector<double>& crds, pMesh msh, pGEntity mdl_ent);
  template <typename T>
    void getEntsFieldComponent(apf::Field * fld, int crd, std::vector<double>& crds, T begin, T end);
  void getEntFieldComponent(apf::Field * fld, int crd, std::vector<double>& crds, pEntity ent);
}
#endif
#include "bioVariableRecovery_impl.h"
