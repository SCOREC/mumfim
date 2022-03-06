#ifndef MUMFIM_MODELTRAITS_H
#define MUMFIM_MODELTRAITS_H
#include <apf.h>
#include <model_traits/ModelTraits.h>
#include <apfMesh.h>
namespace mumfim
{
  void GetModelTraitNodeGeometry(apf::Mesh * mesh,
                              const mt::ModelTraitNode * mtn,
                              std::vector<apf::ModelEntity *> & ents);
}
#endif  // BIOTISSUE_BIOMODELTRAITS_H
