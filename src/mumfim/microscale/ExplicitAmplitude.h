#ifndef MUMFIM_AMPLITUDE_H_
#define MUMFIM_AMPLITUDE_H_
#include <iostream>
namespace mumfim
{
// FIXME Use CRTP here because we call this a lot, and we don't want class
// overhead!
#define diffEPS 1E-5
#define PI 3.141592653589793238462643383279502884197169399375105820974
  struct Amplitude
  {
    virtual double operator()(double time) const = 0;
    [[nodiscard]] virtual double derivative(double t) const
    {
      return (operator()(t + diffEPS) - operator()(t - diffEPS)) /
             (2 * diffEPS);
    }
    [[nodiscard]] virtual double secondDerivative(double t) const
    {
      return (operator()(t + diffEPS) -
              2 * operator()(t) + operator()(t - diffEPS)) /
             (diffEPS * diffEPS);
    }
    virtual ~Amplitude() = default;
  };
  struct SmoothAmp : public Amplitude
  {
    double t_max;
    explicit SmoothAmp(double t_max) : t_max(t_max) { }
    double operator()(double t) const override
    {
      if ( t<= t_max)
      {
        double xi = t / t_max;
        return xi * xi * xi * (10 - 15 * xi + 6 * xi * xi);
      }
      else
        return 1.0;
    }
    virtual double derivative(double t) const override
    {
      if(t <= t_max)
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
      if(t <= t_max)
      {
        double xi = t / t_max;
        double xisq = xi * xi;
        return 1.0 / (t_max * t_max) * (60 * xi - 180 * xisq + 120 * xisq * xi);
      }
      else
        return 0;
    }
  };
  struct SmoothAmpHold : public Amplitude
  {
    double t_max;
    double t_hold;
    SmoothAmpHold(double t_max, double t_hold) : t_max(t_max), t_hold(t_hold)
    {
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
}  // namespace mumfim
#endif
