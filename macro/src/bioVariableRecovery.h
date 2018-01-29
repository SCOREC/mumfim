#ifndef BIO_VARIABLE_RECOVERY_H_
#define BIO_VARIABLE_RECOVERY_H_
#include <apf.h>
#include <apfMesh.h>
#include <MeshSim.h>
#include <vector>
namespace bio
{
  template <typename O>
    void getFieldComponentsClassified(apf::Field * fld,
                                      apf::ModelEntity * mdl_ent,
                                      int nm_cmps,
                                      int * cmps,
                                      O out);
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
#include "bioVariableRecovery_impl.h"
#endif
