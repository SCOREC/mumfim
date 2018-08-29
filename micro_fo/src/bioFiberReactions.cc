#include "bioFiberReactions.h"
namespace bio
{
  std::pair<double, double> LinearReaction::forceReaction(double orig_length,
                                                          double length) const
  {
    std::pair<double, double> result;
    double length_ratio = length / orig_length;
    result.first = ((E * fiber_area) / orig_length) * (length_ratio - 1.0);
    result.second = (E * fiber_area) / (orig_length * orig_length);
    return result;
  }
  std::pair<double, double> NonlinearReaction::forceReaction(
      double orig_length,
      double length) const
  {
    std::pair<double, double> result;
    // Reminder: Every time you change E you have to change the derivatives of
    // the forces too. See gfunc.c
    double length_ratio = length / orig_length;
    double tension = E * fiber_area / B;
    double cmp_ratio_trns = 1.0;
    if (length_ratio > length_ratio_trns)
    {
      double green_strain_trns =
          0.5 * (length_ratio_trns * length_ratio_trns - 1.0);
      double expBeps_trns = exp(B * green_strain_trns);
      double m = fiber_area * E * length_ratio_trns * expBeps_trns;  // slope
      double b = tension * (expBeps_trns - 1.0) -
                 m * length_ratio_trns;  // y-intercept
      result.first = m * length_ratio + b;
      result.second = m / orig_length;
    }
    else if (length_ratio < cmp_ratio_trns)
    {
      double green_strain_trns = 0.5 * (cmp_ratio_trns * cmp_ratio_trns - 1.0);
      double expBeps_trns = exp(B * green_strain_trns);
      double m = fiber_area * E * cmp_ratio_trns * expBeps_trns;  // slope
      double b =
          tension * (expBeps_trns - 1.0) - m * cmp_ratio_trns;  // y-intercept
      result.first = m * length_ratio + b;
      result.second = m / orig_length;
    }
    else
    {
      double green_strain = 0.5 * (length_ratio * length_ratio - 1.0);
      double expBeps = exp(B * green_strain);
      result.first = tension * (expBeps - 1.0);
      result.second =
          ((fiber_area * E * length) / (orig_length * orig_length)) * expBeps;
    }
    return result;
  }
  std::pair<double, double> BeamReaction::forceReaction(double, double) const
  {
    return std::make_pair(0.0, 0.0);
  }
}  // namespace bio
