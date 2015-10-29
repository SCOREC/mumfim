#ifndef BIO_FIBER_REACTIONS_H_
#define BIO_FIBER_REACTIONS_H_

#include <cmath>
#include <utility>

namespace bio
{
  struct Element;// fiber
  class FiberReaction
  {
  public:
    virtual double f(const Element & e) = 0; 
    virtual double df_dl(const Element & e) = 0;
    //virtual std::pair<double,double> forceReaction(double,double) const= 0;
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
      : fbr_ar(fa), B(b), E(e), lexp(le)
    {}
    virtual double f(const Element & e);
    virtual double df_dl(const Element & e);
  };
  
  // Note: essentially just a linear spring
  class LinearReaction : public FiberReaction
  {
  protected:
    double fbr_ar;
    double E;
  public:
    LinearReaction(double fa, double e)
      : fbr_ar(fa), E(e)
    {}
    virtual double f(const Element & e);
    virtual double df_dl(const Element & e);
  };

  /*
  struct BeamReaction : public FiberReaction
  {
    double * forceReaction(double,double) {return NULL;}
  };
  */
}
#endif
