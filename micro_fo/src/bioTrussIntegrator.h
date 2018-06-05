#ifndef BIO_TRUSS_INTEGRATOR_H_
#define BIO_TRUSS_INTEGRATOR_H_
#include "bioFiberNetwork.h"
#include "bioFiberReactions.h"
#include "bioUtil.h"
#include <amsiElementalSystem.h>
#include <las.h>
#include <lasConfig.h>
#include <apf.h>
#include <apfDynamicMatrix.h>
#include <apfElementalSystem.h> //amsi
#include <apfFunctions.h> //amsi
#include <apfMeasure.h> //amsi
#include <cassert>
namespace las
{
  class Mat;
  class Vec;
}
namespace bio
{
  class ElementalSystem;
  /**
   * An integrator class to produce elemental jacobian matrices
   *  and force vectors for trusses using the provided FiberReaction.
   *  The elemental systems are assembled into the global system formed
   *  by the skMat k and skVec f using the supplied numbering, which
   *  should already be numbered.
   * @todo Bill : refactor and simplify, it feels like there is too much
   *              happening in this class
   */
  class TrussIntegrator : public apf::Integrator
  {
  protected:
    apf::Field * u;
    apf::Field * xu;
    apf::Mesh * msh;
    apf::MeshTag * id_tg; //debug
    apf::Numbering * nm;
    apf::Element * elmt;
    int nen;
    apf::NewArray<apf::Vector3> N;
    amsi::ElementalSystem2 * es;
    double l;
    double lo;
    apf::Vector3 spans_l;
    int dim;
    FiberReaction ** frs;
    FiberReaction * fr;
    apf::MeshTag * rct_tg;
    las::Mat * k;
    las::Vec * f;
  public:
  TrussIntegrator(apf::Numbering * n,
                  apf::Field * u_,
                  apf::Field * xu_,
                  FiberReaction ** frs_,
                  las::Mat * k_,
                  las::Vec * f_,
                  int o)
    : apf::Integrator(o)
      , u(u_)
      , xu(xu_)
      , msh(apf::getMesh(u))
      , id_tg(msh->findTag("id"))
      , nm(n)
      , elmt(NULL)
      , nen(0)
      , N()
      , es(NULL)
      , l(0.0)
      , lo(0.0)
      , spans_l()
      , dim()
      , frs(frs_)
      , fr()
      , rct_tg(msh->findTag("fiber_reaction"))
      , k(k_)
      , f(f_)
    { }
    void inElement(apf::MeshElement * me)
    {
      elmt = apf::createElement(u,me);
      nen = apf::countNodes(elmt);
      apf::getVectorNodes(elmt,N);
      lo = apf::measure(me);
      apf::MeshEntity * ent = apf::getMeshEntity(me);
      apf::MeshElement * ume = apf::createMeshElement(xu,ent);
      l = apf::measure(ume);
      apf::destroyMeshElement(ume);
      es = amsi::buildApfElementalSystem(elmt,nm);
      es->zero();
      int tg = -1;
      msh->getIntTag(ent,rct_tg,&tg);
      assert(tg == 0 || tg == 1);
      fr = frs[tg];
      dim = msh->getDimension();
      apf::MeshEntity * vs[2];
      msh->getDownward(ent,0,&vs[0]);
      apf::Element * dsp_elm = apf::createElement(xu,me);
      apf::NewArray<apf::Vector3> crds;
      apf::getVectorNodes(dsp_elm,crds);
      int nd_cnt = apf::countNodes(dsp_elm);
      apf::destroyElement(dsp_elm);
      spans_l = (crds[nd_cnt-1] - crds[0]) / l;
      int id = -1;
      msh->getIntTag(ent,id_tg,&id);
    }
    void atPoint(const apf::Vector3 &, double, double)
    {
      auto f_dfdl = fr->forceReaction(lo,l);
      double f = f_dfdl.first;
      double dfdl = f_dfdl.second;
      double fl = f / l;
      double dfdl_fl = dfdl - fl;
      double frc = 0.0;
      for(int ii = 0; ii < dim; ii++)
      {
        // parametric coordinate along length
        frc = spans_l[ii] * f;
        es->fe(ii)     = -frc;
        es->fe(dim+ii) =  frc;
      }
      apf::Matrix3x3 rctn = apf::tensorProduct(spans_l,spans_l*dfdl_fl) + eye()*fl;
      double op = -1.0;
      for(int ii = 0; ii < 2; ii++)
      {
        op *= -1.0;
        for(int jj = 0; jj < 2; jj++)
        {
          op *= -1.0;
          for(int kk = 0; kk < dim; kk++)
            for(int ll = 0; ll < dim; ll++)
              es->ke(ii*dim + kk, jj*dim + ll) += rctn[kk][ll] * op;
        }
      }
      /*
      // have ke and fe, modify fe based on fixed dofs
      apf::NewArray<int> dofs;
      apf::getElementNumbers(nm,apf::getMeshEntity(elmt),dofs);
      apf::DynamicVector u_fxd(nen*dim);
      u_fxd.zero();
      for(int en = 0; en < nen; ++en)
      {
        for(int ii = 0; ii < dim; ++ii)
        {
          int ldof = en * dim + ii;
          if(dofs[ldof] < 0) // if fixed
            u_fxd[ldof] = N[en][ii];
        }
      }
      // would prefer matrix-vector multiplication, but ke is wrapped up in a generic interface so can't use it as a DynamicMatrix
      for(int ii = 0; ii < nen*dim; ++ii)
        for(int jj = 0; jj < nen*dim; ++jj)
          es->fe(ii) += es->ke(ii,jj) * u_fxd(jj);
      */
    }
    void outElement()
    {
      auto ops = las::getLASOps<las::sparskit>();
      ops->assemble(k,es->nedof(),&es->dofs(0),es->nedof(),&es->dofs(0),&es->ke(0,0));
      ops->assemble(f,es->nedof(),&es->dofs(0),&es->fe(0));
      amsi::destroyApfElementalSystem(es);
      apf::destroyElement(elmt);
    }
  };
}
#endif
