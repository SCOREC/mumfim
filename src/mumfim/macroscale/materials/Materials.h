#ifndef MUMFIM_MACROSCALE_MATERIALS_H
#define MUMFIM_MACROSCALE_MATERIALS_H
#include <amsiDeformation.h>
#include <mumfim/exceptions.h>
#include <mumfim/macroscale/materials/MaterialResult.h>

// Note: If the function call overhead for these functions is too high from
// the hot loop of the integrator, we may put them in header to allow inline
namespace mumfim
{
  struct micro_fo_result;

  [[nodiscard]] MaterialResult NeohookeanMaterial(
      double shear_modulus,
      double poissons_ratio,
      const apf::Matrix3x3 & deformation_gradient);
  [[nodiscard]] MaterialResult MultiscaleMaterial(micro_fo_result * rslt);
}  // namespace mumfim
#ifdef MUMFIM_ENABLE_TORCH
#include <mumfim/macroscale/materials/TorchMaterial.h>
#endif

#endif
