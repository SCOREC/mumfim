#ifndef BIO_TRUSS_INTEGRATOR_H_
#define BIO_TRUSS_INTEGRATOR_H_
#include <amsiElementalSystem.h>
#include <apf.h>
#include <apfDynamicMatrix.h>
#include <apfElementalSystem.h>  //amsi
#include <apfFunctions.h>        //amsi
#include <apfMeasure.h>          //amsi
#include <las.h>
#include <lasConfig.h>
#include <cassert>
#include <limits>
#include "bioFiberNetwork.h"
#include "bioFiberReactions.h"
#include "bioUtil.h"
namespace las
{
  class Mat;
  class Vec;
}  // namespace las
namespace bio
{
  // class ElementalSystem;
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
    apf::MeshTag * id_tg;  // debug
    apf::Numbering * nm;
    apf::Element * elmt;
    int nen;
    amsi::ElementalSystem2 * es;
    double l;
    double lo;
    apf::Vector3 spans_l;
    int dim;
    FiberNetwork::reaction_ptr_type frs;
    // the integer tag that corresponds to the fiber reaction
    int tg;
    apf::MeshTag * rct_tg;
    las::Mat * k;
    las::Vec * f;
    int id;
    public:
    TrussIntegrator(apf::Numbering * n,
                    apf::Field * u,
                    apf::Field * xu,
                    FiberNetwork::reaction_ptr_type frs,
                    las::Mat * k,
                    las::Vec * f,
                    int o)
        : apf::Integrator(o)
        , u(u)
        , xu(xu)
        , msh(apf::getMesh(u))
        , id_tg(msh->findTag("id"))
        , nm(n)
        , elmt(NULL)
        , nen(0)
        //, N()
        , es(NULL)
        , l(0.0)
        , lo(0.0)
        , spans_l()
        , dim()
        , frs(frs)
        , rct_tg(msh->findTag("fiber_reaction"))
        , k(k)
        , f(f)
        , id(-1)
    {
    }
    virtual void inElement(apf::MeshElement * me);
    virtual void atPoint(const apf::Vector3 &, double, double);
    virtual void outElement();
  };
}  // namespace bio
#endif
