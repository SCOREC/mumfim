#include "bioVolumeConvergence.h"
#include <amsiCasters.h>
#include <amsiNonlinearAnalysis.h>
#include <amsiNonlinearConvergenceOperator.h>
#include <gmi.h>
#include "bioModelTraits.h"
namespace mumfim
{
  struct CalcDV : public amsi::to_R1
  {
    VolCalc * v;
    double operator()() { return fabs(v->getV() - v->getPV()); }
  };
  struct CalcDV0 : public amsi::to_R1
  {
    VolCalc * v;
    double operator()() { return fabs(v->getV() - v->getV0()); }
  };
  struct CalcDVPS : public amsi::to_R1
  {
    VolCalc * v;
    double operator()() { return fabs(v->getV() - v->getVPS()); }
  };
  struct CalcPV : public amsi::to_R1
  {
    VolCalc * v;
    double operator()() { return v->getPV(); }
  };
  struct CalcV : public amsi::to_R1
  {
    VolCalc * v;
    double operator()() { return v->getV(); }
  };
  struct CalcV0 : public amsi::to_R1
  {
    VolCalc * v;
    double operator()() { return v->getV0(); }
  };
  struct CalcVPS : public amsi::to_R1
  {
    VolCalc * v;
    double operator()() { return v->getVPS(); }
  };
  std::unique_ptr<amsi::Convergence> buildVolConvergenceOperator(
      const mt::CategoryNode & category_node,
      amsi::MultiIteration * it,
      VolCalc * vl,
      apf::Field * fld)
  {
    if (category_node.GetType() != "volume convergence")
    {
      std::cerr << "volume convergence should be created from a convergence "
                   "operator with the volume convergence type.\n";
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    const auto * regions_nd = category_node.FindCategoryNodeByType("regions");
    if (regions_nd == nullptr)
    {
      std::cerr << "volume convergence operator must have regions.\n";
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    const auto * track_volume =
        mt::GetCategoryModelTraitNodeByType(regions_nd, "track volume");
    if (track_volume == nullptr)
    {
      std::cerr
          << "volume convergence operator should have a track volume region.\n";
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    std::vector<apf::ModelEntity *> mdl_ents;
    auto * mesh = apf::getMesh(fld);
    GetModelTraitNodeGeometry(mesh, track_volume, mdl_ents);
    std::cout
        << "Volume convergence operator discovered, effects model entitites : ";
    gmi_model * mdl = apf::getMesh(fld)->getModel();
    for (auto mdl_ent = mdl_ents.begin(); mdl_ent != mdl_ents.end(); ++mdl_ent)
      std::cout << mdl->ops->tag(mdl, (gmi_ent *)*mdl_ent) << " ";
    std::cout << std::endl;
    const auto * reference_mt = mt::GetCategoryModelTraitByType<mt::IntMT>(
        &category_node, "reference value");
    if (reference_mt == nullptr)
    {
      std::cerr << "volume convergence should have a reference value\n";
      MPI_Abort(AMSI_COMM_SCALE, 1);
    }
    int ref_tp = (*reference_mt)();
    amsi::to_R1 * dv = NULL;
    amsi::to_R1 * ref_v = NULL;
    if (ref_tp == 0)  //  initial
    {
      CalcDV0 * v = new CalcDV0;
      v->v = vl;
      dv = v;
      CalcV0 * r = new CalcV0;
      r->v = vl;
      ref_v = r;
    }
    else if (ref_tp == 1)  // load_step
    {
      CalcDVPS * v = new CalcDVPS;
      v->v = vl;
      dv = v;
      CalcVPS * r = new CalcVPS;
      r->v = vl;
      ref_v = r;
    }
    else if (ref_tp == 2)  // iteration
    {
      CalcDV * v = new CalcDV;
      v->v = vl;
      dv = v;
      CalcVPS * r = new CalcVPS;
      r->v = vl;
      ref_v = r;
    }
    const auto * cap_mt =
        mt::MTCast<mt::IntMT>(category_node.FindModelTraitNode("iteration cap")
                                  ->GetModelTraits()[0]
                                  .second.get());
    const auto * eps_mt = category_node.FindModelTraitNode("epsilon")
                              ->GetModelTraits()[0]
                              .second.get();
    if (cap_mt)
    {
      // lifetime of this iteration is linked to the lifetime of the associated
      // MultiIteration
      amsi::Iteration * stop_at_max_iters =
          new amsi::StopAtMaxIters((*cap_mt)());
      it->addIteration(stop_at_max_iters);
    }
    const auto * scalar_eps = mt::MTCast<mt::ScalarMT>(eps_mt);
    if (scalar_eps)
    {
      return createUpdatingConvergence(
          it, dv, new amsi::MTConstantEpsilon((*scalar_eps)()), ref_v);
    }
    const auto * function_eps = mt::MTCast<mt::ScalarFunctionMT<1>>(eps_mt);
    if (function_eps)
    {
      return createUpdatingConvergence(
          it, dv, new amsi::MTUpdatingEpsilon(*function_eps), ref_v);
    }
    std::cerr << "Epsilon must be function or scalar type.\n";
    MPI_Abort(AMSI_COMM_WORLD, 1);
    return nullptr;
  }
}  // namespace mumfim
