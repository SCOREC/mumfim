#ifndef BIO_FIBER_REACTIONS_H_
#define BIO_FIBER_REACTIONS_H_

#include <cmath>
#include <utility>
#include <vector>

namespace bio
{

  class FiberReaction
  {
  public:
    virtual double F(double l, double l0) = 0;
    virtual double dFdl(double l, double l0) = 0;
  };

  class NonlinearReaction
  {
  protected:
    double E;
    double A;
    double B;
  public:
    void setYoungs(double e) { E = e; }
    void setFiberArea(double a) { A = a; }
    void setNonlinearity(double b) { B = b; }

    double F(double l, double l0)
    {
      return ((E*A)/B) * (exp((0.5*B)*((l/l0)*(l/l0) - 1) - 1));
    }

    double dFdl(double l, double l0)
    {
      return ((E*A*l)/(l0*l0)) * (exp((0.5*B)*((l/l0)*(l/l0) - 1)));
    }
  };
/*
  class FiberReaction
  {
  protected:
    std::string nm;
  public:
    FiberReaction(const std::string & n) : nm(n) {}
    virtual double f(const Element & original, const Element & deformed) = 0; 
    virtual double df_dl(const Element & original, const Element & deformed) = 0;
    const std::string & getName() {return nm;}
  };

  class FiberReactionAssignment
  {
  protected:
    std::vector<FiberReaction*> reactions;
    std::vector<int> assignment;
  public:
    FiberReactionAssignment();
    FiberReaction * getReaction(int fiber) { return reactions[assignment[fiber]]; }
    int addReaction(FiberReaction * r) { reactions.push_back(r); }
    void assignReaction(int fiber, int reaction) {assignment[fiber] = reaction;}
  };
  
  class NonLinearReaction : public FiberReaction
  {
  protected:
    double fbr_ar;
    double B;
    double E;
    double lexp;
    NonLinearReaction();
  public:
    NonLinearReaction(double fa, double b, double e, double le)
      : FiberReaction("nonlinear"), fbr_ar(fa), B(b), E(e), lexp(le)
    {}
    virtual double f(const Element & original, const Element & deformed);
    virtual double df_dl(const Element & original, const Element & deformed);
  };
  
  // Note: essentially just a linear spring
  class LinearReaction : public FiberReaction
  {
  protected:
    double fbr_ar;
    double E;
  public:
    LinearReaction(double fa, double e)
      : FiberReaction("linear"), fbr_ar(fa), E(e)
    {}
    virtual double f(const Element & original, const Element & deformed);
    virtual double df_dl(const Element & original, const Element & deformed);
  };
*/
  /*
  struct BeamReaction : public FiberReaction
  {
    double * forceReaction(double,double) {return NULL;}
  };
  */
}
#endif
