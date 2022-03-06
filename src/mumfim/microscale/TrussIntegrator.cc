#include "TrussIntegrator.h"
#include "Utility.h"
namespace mumfim
{
  void TrussIntegrator::inElement(apf::MeshElement * me)
  {
    elmt = apf::createElement(u, me);
    nen = apf::countNodes(elmt);
    // apf::getVectorNodes(elmt,N);
    lo = apf::measure(me);
    apf::MeshEntity * ent = apf::getMeshEntity(me);
    apf::MeshElement * ume = apf::createMeshElement(xu, ent);
    l = apf::measure(ume);
    apf::destroyMeshElement(ume);
    es = amsi::buildApfElementalSystem(elmt, nm);
    es->zero();
    msh->getIntTag(ent, rct_tg, &tg);
    dim = msh->getDimension();
    apf::MeshEntity * vs[2];
    msh->getDownward(ent, 0, &vs[0]);
    apf::Element * dsp_elm = apf::createElement(xu, me);
    apf::NewArray<apf::Vector3> crds;
    apf::getVectorNodes(dsp_elm, crds);
    int nd_cnt = apf::countNodes(dsp_elm);
    apf::destroyElement(dsp_elm);
    // unit vector of the fiber
    spans_l = (crds[nd_cnt - 1] - crds[0]) / l;
  }
  void TrussIntegrator::atPoint(const apf::Vector3 &, double, double)
  {
    auto f_dfdl = (*frs)[tg].forceReaction(lo, l);
    double f = f_dfdl.first;
    // the scalar version of the force derivative
    double dfdl = f_dfdl.second;
    // why are we dividing by l here?
    double fl = f / l;
    double dfdl_fl = dfdl - fl;
    double frc = 0.0;
    for (int ii = 0; ii < dim; ii++)
    {
      // parametric coordinate along length
      frc = spans_l[ii] * f;
      es->fe(ii) = -frc;
      es->fe(dim + ii) = frc;
    }
    // rctn is df/du
    apf::Matrix3x3 rctn =
        apf::tensorProduct(spans_l, spans_l * dfdl_fl) + eye() * fl;
    double op = -1.0;
    for (int ii = 0; ii < 2; ii++)
    {
      op *= -1.0;
      for (int jj = 0; jj < 2; jj++)
      {
        op *= -1.0;
        for (int kk = 0; kk < dim; kk++)
          for (int ll = 0; ll < dim; ll++)
            es->ke(ii * dim + kk, jj * dim + ll) += rctn[kk][ll] * op;
      }
    }
  }
  void TrussIntegrator::outElement()
  {
    auto ops = las::getLASOps<las::MICRO_BACKEND>();
    ops->assemble(k, es->nedof(), &es->dofs(0), es->nedof(), &es->dofs(0),
                  &es->ke(0, 0));
    ops->assemble(f, es->nedof(), &es->dofs(0), &es->fe(0));
    amsi::destroyApfElementalSystem(es);
    apf::destroyElement(elmt);
  }
}  // namespace mumfim
