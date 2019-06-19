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
    FiberReaction ** frs;
    FiberReaction * fr;
    apf::MeshTag * rct_tg;
    las::Mat * k;
    las::Vec * f;
    int id;
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
        //, N()
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
        , id(-1)
    {
    }
    virtual void inElement(apf::MeshElement * me);
    virtual void atPoint(const apf::Vector3 &, double, double);
    virtual void outElement();
  };
  class ExplicitTrussIntegrator : public TrussIntegrator
  {
    protected:
    double sound_speed;
    double delta_t_crit;
    double delta_t_crit_elmt;
    double stiffness_damping_factor;
    // velocity
    las::Vec * v;
    las::Vec * f_damp;
    double * vel_arr;
    apf::DynamicVector fe_damp;
    public:
    virtual void inElement(apf::MeshElement * me) override;
    virtual void atPoint(const apf::Vector3 &, double, double) override;
    virtual void outElement() override;
    void setStiffnessDampingFactor(double alpha) {stiffness_damping_factor = alpha; }
    // call this before processing a second time so we capture an increase in
    // the critical time step.
    void resetCritTimeStep()
    {
      delta_t_crit = std::numeric_limits<double>::max();
    }
    void updateF(las::Vec * vec) { f = vec; }
    void updateFDamp(las::Vec * vec) { f_damp = vec; }
    double getCriticalTimeStep() const { return delta_t_crit; }
    ExplicitTrussIntegrator(apf::Numbering * n,
                            apf::Field * u_,
                            apf::Field * xu_,
                            FiberReaction ** frs_,
                            las::Mat * k_,
                            las::Vec * f_,
                            las::Vec * f_damp_,
                            las::Vec * v_,
                            int o)
        : TrussIntegrator(n, u_, xu_, frs_, k_, f_, o)
        , delta_t_crit(std::numeric_limits<double>::max())
        , stiffness_damping_factor(0)
        , v(v_)
        , f_damp(f_damp_)
    {
      // std::cout<<"Explicit Truss Integrator"<<std::endl;
      // std::cout<<es->nedof()<<std::endl;
      fe_damp = apf::DynamicVector();
    }
  };
}  // namespace bio
#endif
