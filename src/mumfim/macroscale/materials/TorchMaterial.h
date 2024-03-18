#ifndef MUMFIM_SRC_MUMFIM_MACROSCALE_MATERIALS_TORCHMATERIAL_H
#define MUMFIM_SRC_MUMFIM_MACROSCALE_MATERIALS_TORCHMATERIAL_H

#include <mumfim/macroscale/materials/MaterialResult.h>
#include <torch/script.h>

namespace mumfim
{

  // PyTorch Material Functor
  struct TorchMaterial
  {
    explicit TorchMaterial(const std::string & model_path);
    explicit TorchMaterial(const torch::jit::Module & Model);

    MaterialResult operator()(const apf::Matrix3x3 & deformation_gradient,
                              apf::MeshEntity * element,
                              int integration_point);

    private:
    // model that takes the right-cauchy-green-deformation tensor and returns
    // the energy
    torch::jit::script::Module model;
    // reuse memory
    torch::Tensor deformation_gradient_tensor;

    // inputs are the right cauchy_green deformation_tensor
    std::vector<torch::jit::IValue> inputs;
  };

}  // namespace mumfim

#endif  // MUMFIM_SRC_MUMFIM_MACROSCALE_MATERIALS_TORCHMATERIAL_H
