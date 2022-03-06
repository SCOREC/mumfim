#include <amsiNonlinearConvergenceOperator.h>
#include <model_traits/AssociatedModelTraits.h>
#include <model_traits/ModelTraits.h>
#include "bioVolumeConvergence.h"
namespace mumfim
{
  template <typename O>
  void buildLASConvergenceOperators(const mt::CategoryNode * solution_strategy,
                                    amsi::MultiIteration * it,
                                    amsi::LAS * las,
                                    O out)
  {
    auto convergence_operators_cat =
        mt::GetCategoriesByType(solution_strategy, "convergence operator");
    for (const auto * convergence_operator_cat : convergence_operators_cat)
    {
      for (const auto & convergence_operator :
           convergence_operator_cat->GetCategoryNodes())
      {
        auto * cvg =
            amsi::createConvergenceOperator(&convergence_operator, it, las)
                .release();
        if (cvg != nullptr)
        {
          *out++ = cvg;
        }
      }
    }
  }
  /**
   * Extracts the volume convergence regions from a simmetrix case
   * @param ss simmetrix case
   * @param it amsi iterator
   * @param u displacement
   * @param vl_tks list of volumes tracked for volume convergence
   * @param out typically std::back_inserter to some sort of list type.
   */
  template <typename I, typename O>
  void buildVolConvergenceOperators(const mt::CategoryNode * solution_strategy,
                                    amsi::MultiIteration * it,
                                    apf::Field * u,
                                    I vl_tks,
                                    O out)
  {
    const auto * convergence_operator =
        mt::GetCategoryByType(solution_strategy, "convergence operator");
    auto volume_convergence_operators =
        mt::GetCategoriesByType(convergence_operator, "volume convergence");
    for (const auto * volume_convergence : volume_convergence_operators)
    {
      const auto * regions =
          mt::GetCategoryByType(volume_convergence, "regions");
      const auto * track_volume =
          mt::GetCategoryModelTraitNodeByType(regions, "track volume");
      if (track_volume == nullptr)
      {
        std::cerr << "track volume required on volume convergence.\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      const auto & model_traits = track_volume->GetModelTraits();
      if (model_traits.size() != 1 || track_volume->GetName().empty())
      {
        std::cerr << "track volume region should only have one model trait and "
                     "must have a name!\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      *out++ = buildVolConvergenceOperator(*volume_convergence, it,
                                           vl_tks[track_volume->GetName()], u)
                   .release();
    }
  }
}  // namespace mumfim
