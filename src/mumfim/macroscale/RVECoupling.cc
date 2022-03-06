#include "RVECoupling.h"
#include <model_traits/AssociatedModelTraits.h>
namespace mumfim
{
  MicroscaleType getMicroscaleType(
      const mt::AssociatedCategoryNode * category_node)
  {
    if (category_node)
    {
      if (category_node->GetType() != "multiscale model")
      {
        std::cerr << "must have multiscale model type\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      if (category_node->FindCategoryByType("isotropic_neohookean") != nullptr)
      {
        return MicroscaleType::ISOTROPIC_NEOHOOKEAN;
      }
      else if (category_node->FindCategoryByType("fiber only") != nullptr)
      {
        return MicroscaleType::FIBER_ONLY;
      }
      else
      {
        std::cerr << "Multiscale model type is not a valid. Must be "
                     "\"isotropic_neohookean\" or \"fiber only\".\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
    }
    return MicroscaleType::NONE;
  }
}  // namespace mumfim
