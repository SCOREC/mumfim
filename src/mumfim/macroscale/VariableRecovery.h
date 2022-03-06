#ifndef MUMFIM_VARIABLE_RECOVERY_H_
#define MUMFIM_VARIABLE_RECOVERY_H_
#include <apf.h>
#include <apfMesh.h>
#include <vector>
namespace mumfim
{
  template <typename I, typename O>
    void getFieldComponentsOnEnts(apf::Field * fld,
                                  int nm_cmps,
                                  int * cmps,
                                  I begin,
                                  I end,
                                  O out);
  template <typename O>
    void getFieldComponentsOnEnt(apf::Field * fld,
                                 apf::MeshEntity * ent,
                                 int nm_cmps,
                                 int * cmps,
                                 O out);
}
#include "VariableRecovery_impl.h"
#endif
