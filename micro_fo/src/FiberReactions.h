#ifndef FIBERREACTIONS_H_
#define FIBERREACTIONS_H_
#include <cmath>
#include <cstddef>
#include <utility>
enum FiberConstitutive
{
  LINEAR = 0,
  NONLINEAR = 1
};
struct FiberReaction
{
  virtual std::pair<double,double> forceReaction(double,double) const= 0;
};
struct NonlinearReaction : public FiberReaction
{
  double fiber_area;
  double B;
  double E;
  double length_ratio_trns;
  std::pair<double,double> forceReaction(double orig_length, double length) const
  {
    std::pair<double,double> result;
    // Reminder: Every time you change E you have to change the derivatives of the forces too. See gfunc.c
    double length_ratio = length / orig_length;
    double tension = E * fiber_area / B;
    if(length_ratio < length_ratio_trns)
//    if(false)
    {
      double green_strain = 0.5 * (length_ratio * length_ratio - 1.0);
      double expBeps = exp(B * green_strain);
      result.first = tension * (expBeps - 1.0);
      result.second = ((fiber_area * E * length) / (orig_length * orig_length)) * expBeps;
    }
    else
    {
      double length_trns = length_ratio_trns * orig_length;
      double green_strain_trns = 0.5 * (length_ratio_trns * length_ratio_trns - 1.0);
      double expBeps_trns = exp(B * green_strain_trns);
      double m = ( (fiber_area * E * length_trns ) / (orig_length * orig_length) )  * expBeps_trns; // slope
      double b = tension * (expBeps_trns - 1.0) - m * length_ratio_trns; // y-intercept 
      result.first = m * length_ratio + b;
      result.second = m / orig_length;      
    }
    return result;
  }
};
struct LinearReaction : public FiberReaction
{
  double fiber_area;
  double E;
  std::pair<double,double> forceReaction(double orig_length, double length) const
  {
    std::pair<double,double> result;
    double length_ratio = length / orig_length;
    result.first = ((E*fiber_area)/orig_length) * (length_ratio - 1.0);
    result.second = (E * fiber_area) / (orig_length * orig_length);
    return result;
  }
};
struct BeamReaction : public FiberReaction
{
  double * forceReaction(double,double) {return NULL;}
};
#endif
