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

  class NonlinearReaction : public FiberReaction
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

  class LinearReaction : public FiberReaction
  {
  protected:
    double E;
    double A;
  public:
    void setYoungs(double e) { E = e; }
    void setFiberArea(double a) { A = a; }

    double F(double l, double l0)
    {
      return ((E*A) / l0)*((l/l0) - 1.0);
    }

    double dFdl(double l, double l0)
    {
      return ((E*A) / (l0*l0));
    }
  };
}
#endif
