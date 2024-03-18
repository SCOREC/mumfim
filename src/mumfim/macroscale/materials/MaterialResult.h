#ifndef MUMFIM_SRC_MUMFIM_MACROSCALE_MATERIALS_MATERIALRESULT_H
#define MUMFIM_SRC_MUMFIM_MACROSCALE_MATERIALS_MATERIALRESULT_H
#include <ElementalSystem.h>

struct MaterialResult
{
  apf::DynamicMatrix cauchy_stress;
  apf::DynamicMatrix material_stiffness;
};

#endif  // MUMFIM_SRC_MUMFIM_MACROSCALE_MATERIALS_MATERIALRESULT_H
