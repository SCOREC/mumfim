#ifndef BIO_FIBERREACTIONS_H_
#define BIO_FIBERREACTIONS_H_
#include <amsiEnumOps.h>
#include <cmath>
#include <cstddef>
#include <utility>
namespace bio
{
#define FBR_RCT_TYPES(OP) OP(linear), OP(nonlinear), OP(num_fbr_constitutive)
  enum FiberConstitutive
  {
    FBR_RCT_TYPES(MAKE_ENUM_OP)
  };
  const char * getFiberConstitutiveString(int ii);
  struct FiberReaction
  {
    virtual std::pair<double, double> forceReaction(double, double) const = 0;
  };
  struct NonlinearReaction : public FiberReaction
  {
    double fiber_area;
    double B;
    double E;
    // tension transition to compressive regime
    double length_ratio_trns;
    // compressive transition to linear regime
    double cmp_ratio_trns = 1.0;
    std::pair<double, double> forceReaction(double orig_length,
                                            double length) const;
  };
  struct LinearReaction : public FiberReaction
  {
    double fiber_area;
    double E;
    std::pair<double, double> forceReaction(double orig_length,
                                            double length) const;
  };
  struct BeamReaction : public FiberReaction
  {
    std::pair<double, double> forceReaction(double, double) const;
  };
}  // namespace bio
#endif
