#ifndef __AMPLITUDE_H__
#define __AMPLITUDE_H__
#include <cassert>
#include <iostream>
namespace bio
{
// FIXME Use CRTP here because we call this a lot, and we don't want class
// overhead!
#define diffEPS 1E-5
#define PI 3.141592653589793238462643383279502884197169399375105820974
  struct Amplitude
  {
    virtual double operator()(double time) const = 0;
    virtual double derivative(double t) const
    {
      return (operator()(t + diffEPS) - operator()(t - diffEPS)) /
             (2 * diffEPS);
    }
    virtual double secondDerivative(double t) const
    {
      return (operator()(t + diffEPS) -
              2 * operator()(t) + operator()(t - diffEPS)) /
             (diffEPS * diffEPS);
    }
    virtual ~Amplitude() {}
  };
  struct SmoothAmp : public Amplitude
  {
    double t_max;
    //SmoothAmp(double t_max) : t_max(t_max) { assert(t_max >= 1.0); }
    SmoothAmp(double t_max) : t_max(t_max) { }
    double operator()(double t) const override
    {
      double xi = t / t_max;
      return xi * xi * xi * (10 - 15 * xi + 6 * xi * xi);
    }
    virtual double derivative(double t) const override
    {
      double xi = t / t_max;
      double xisq = xi * xi;
      return 1.0 / t_max * (30 * xisq - 60 * xi * xisq + 30 * xisq * xisq);
    }
    virtual double secondDerivative(double t) const override
    {
      double xi = t / t_max;
      double xisq = xi * xi;
      return 1.0 / (t_max * t_max) * (60 * xi - 180 * xisq + 120 * xisq * xi);
    }
  };
  struct SmoothAmpHold : public Amplitude
  {
    double t_max;
    double t_hold;
    SmoothAmpHold(double t_max, double t_hold) : t_max(t_max), t_hold(t_hold)
    {
      //assert(t_max >= 1.0);
    }
    double operator()(double t) const override
    {
      if (t <= t_max)
      {
        double xi = t / t_max;
        return xi * xi * xi * (10 - 15 * xi + 6 * xi * xi);
      }
      else
        return 1.0;
    }
    virtual double derivative(double t) const override
    {
      if (t <= t_max)
      {
        double xi = t / t_max;
        double xisq = xi * xi;
        return 1.0 / t_max * (30 * xisq - 60 * xi * xisq + 30 * xisq * xisq);
      }
      else
        return 0;
    }
    virtual double secondDerivative(double t) const override
    {
      if (t <= t_max)
      {
        double xi = t / t_max;
        double xisq = xi * xi;
        return 1.0 / (t_max * t_max) * (60 * xi - 180 * xisq + 120 * xisq * xi);
      }
      else
        return 0;
    }
  };
}  // namespace bio
#endif
