#include "bioFiberReactions.h"
#include <vector>
namespace mumfim
{
  std::pair<double, double> LinearReaction::forceReaction(double orig_length,
                                                          double length) const
  {
    std::pair<double, double> result;
    double length_ratio = length / orig_length;
    double green_strain = 1.0 / 2.0 * (length_ratio * length_ratio - 1);
    result.first = length_ratio * E * fiber_area * green_strain;
    result.second = E * fiber_area / orig_length *
                    (green_strain + length_ratio * length_ratio);
    return result;
  }
  static double computeNLForce(double alpha,
                               double B,
                               double length_ratio,
                               double green_strain)
  {
    return alpha * length_ratio * (exp(B * green_strain) - 1);
  }
  static double computeNLDeriv(double alpha,
                               double B,
                               double length_ratio,
                               double orig_length,
                               double green_strain)
  {
    return alpha / orig_length *
           ((1 + length_ratio * length_ratio * B) * exp(B * green_strain) - 1);
  }
  std::pair<double, double> NonlinearReaction::forceReaction(
      double orig_length,
      double length) const
  {
    std::pair<double, double> result;
    double length_ratio = length / orig_length;
    double alpha = E * fiber_area / B;
    if (length_ratio > length_ratio_trns)
    {
      double green_strain_trns =
          1.0 / 2.0 * (length_ratio_trns * length_ratio_trns - 1);
      // slope
      double m = computeNLDeriv(alpha, B, length_ratio_trns, orig_length,
                                green_strain_trns);
      // intercept
      double b = computeNLForce(alpha, B, length_ratio_trns, green_strain_trns);
      result.first = m * (length_ratio - length_ratio_trns) * orig_length + b;
      result.second = m;
    }
    else if (length_ratio < cmp_ratio_trns)
    {
      double green_strain_cmp =
          1.0 / 2.0 * (cmp_ratio_trns * cmp_ratio_trns - 1);
      // slope
      double m = computeNLDeriv(alpha, B, cmp_ratio_trns, orig_length,
                                green_strain_cmp);
      // intercept
      double b = computeNLForce(alpha, B, cmp_ratio_trns, green_strain_cmp);
      result.first = m * (length_ratio - cmp_ratio_trns) * orig_length + b;
      result.second = m;
    }
    else
    {
      double green_strain = 1.0 / 2.0 * (length_ratio * length_ratio - 1);
      result.first = computeNLForce(alpha, B, length_ratio, green_strain);
      result.second =
          computeNLDeriv(alpha, B, length_ratio, orig_length, green_strain);
    }
    return result;
  }
  std::pair<double, double> BeamReaction::forceReaction(double, double) const
  {
    return std::make_pair(0.0, 0.0);
  }
}  // namespace mumfim
