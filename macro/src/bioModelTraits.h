#ifndef BIOTISSUE_BIOMODELTRAITS_H
#define BIOTISSUE_BIOMODELTRAITS_H
#include <apf.h>
#include <model_traits/ModelTraits.h>
#include <apfMesh.h>
namespace bio
{
  void GetModelTraitNodeGeometry(apf::Mesh * mesh,
                              const mt::ModelTraitNode * mtn,
                              std::vector<apf::ModelEntity *> & ents);
}
#endif  // BIOTISSUE_BIOMODELTRAITS_H
