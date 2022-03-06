#ifndef MUMFIM_FIBERREACTIONS_H_
#define MUMFIM_FIBERREACTIONS_H_
#include <amsiEnumOps.h>
#include <cmath>
#include <cstddef>
#include <utility>
namespace mumfim
{
#define FBR_RCT_TYPES(OP) OP(linear), OP(nonlinear), OP(num_fbr_constitutive)
  enum FiberConstitutive
  {
    FBR_RCT_TYPES(MAKE_ENUM_OP)
  };
  const char * getFiberConstitutiveString(int ii);
  struct FiberReaction
  {
    virtual ~FiberReaction(){}
    virtual std::pair<double, double> forceReaction(double, double) const = 0;
    [[nodiscard]] double getYoungModulus() const noexcept { return E; }
    [[nodiscard]] double getFiberArea() const noexcept { return fiber_area; }
    [[nodiscard]] double getFiberDensity() const noexcept
    {
      return fiber_density;
    }
    double fiber_area;
    double fiber_density;
    double E;
  };
  struct NonlinearReaction : public FiberReaction
  {
    double B;
    // tension transition to compressive regime
    double length_ratio_trns;
    // compressive transition to linear regime
    double cmp_ratio_trns;
    std::pair<double, double> forceReaction(double orig_length,
                                            double length) const;
  };
  struct LinearReaction : public FiberReaction
  {
    std::pair<double, double> forceReaction(double orig_length,
                                            double length) const;
  };
  struct BeamReaction : public FiberReaction
  {
    std::pair<double, double> forceReaction(double, double) const;
  };
}  // namespace mumfim
#endif
