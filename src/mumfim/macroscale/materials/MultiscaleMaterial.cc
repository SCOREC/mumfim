#include "Materials.h"
#include "mumfim/macroscale/RVECoupling.h"

namespace mumfim {
  MaterialResult MultiscaleMaterial(micro_fo_result * rslt) {
    apf::DynamicMatrix cauchy_stress(3,3);
    apf::DynamicMatrix material_stiffness(6,6);
    cauchy_stress(0,0) = rslt->data[0];
    cauchy_stress(1,1) = rslt->data[1];
    cauchy_stress(2,2) = rslt->data[2];
    cauchy_stress(1,2) = cauchy_stress(2,1) = rslt->data[3];
    cauchy_stress(0,2) = cauchy_stress(2,0) = rslt->data[4];
    cauchy_stress(0,1) = cauchy_stress(1,0) = rslt->data[5];
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 6; ++j)
        material_stiffness(i, j) = rslt->data[9 + i * 6 + j];
    return {.cauchy_stress = cauchy_stress, .material_stiffness = material_stiffness};
  }
}
