#include "bioModelTraits.h"
namespace bio
{
  void GetModelTraitNodeGeometry(apf::Mesh * mesh,
                              const mt::ModelTraitNode * mtn,
                              std::vector<apf::ModelEntity *> & ents)
  {
    for (auto & trait : mtn->GetModelTraits())
    {
      const auto& geometry_set = std::dynamic_pointer_cast<mt::DimIdGeometrySet>(trait.first);
      if(geometry_set != nullptr) {
        for(const auto& geom : *geometry_set)
        {
          ents.push_back(mesh->findModelEntity(geom.GetDimension(),geom.GetID()));
        }
      }
    }
  }
}
