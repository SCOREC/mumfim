#include <apfDynamicMatrix.h>
#include <apfMatrixUtil.h>
#include <fmt/format.h>
#include <mumfim/exceptions.h>
#include <mumfim/macroscale/materials/TorchMaterial.h>
#include <torch/script.h>
#include <torch/torch.h>

#include <iostream>

namespace mumfim
{
  static void ComputeRightCauchyGreen(
      const torch::Tensor & deformation_gradient,
      torch::Tensor & right_cauchy_green)
  {
    torch::matmul_outf(torch::transpose(deformation_gradient, 1, 2),
                       deformation_gradient, right_cauchy_green);
  }

  [[nodiscard]] static at::Tensor ComputeRightCauchyGreen(
      at::Tensor & deformation_gradient)
  {
    return torch::einsum("kli,klj->kij",
                         {deformation_gradient, deformation_gradient});
  }

  [[nodiscard]] static int VoigtIndex(int i, int j)
  {
    // 00, 11, 22, 12, 02, 01
    if (i > j)
    {
      std::swap(i, j);
    }
    if (i == j)
    {
      return i;
    }
    if (i == 1 && j == 2)
    {
      return 3;
    }
    else if (i == 0 && j == 2)
    {
      return 4;
    }
    else if (i == 0 && j == 1)
    {
      return 5;
    }
    else
    {
      throw material_error{"Invalid Voigt index"};
    }
  }

  [[nodiscard]] static std::pair<int, int> Tensor4Index(int i)
  {
    if (i == 0)
    {
      return {0, 0};
    }
    else if (i == 1)
    {
      return {1, 1};
    }
    else if (i == 2)
    {
      return {2, 2};
    }
    else if (i == 3)
    {
      return {1, 2};
    }
    else if (i == 4)
    {
      return {0, 2};
    }
    else if (i == 5)
    {
      return {0, 1};
    }
    else
    {
      throw material_error{"Invalid Tensor4 index"};
    }
  }

  [[nodiscard]] static torch::Tensor Tensor4ToVoigt(
      const torch::Tensor & tensor)
  {
    using torch::indexing::None;
    using torch::indexing::Slice;
    if (tensor.ndimension() != 5)
      throw material_error{"Tensor4ToVoigt only works for 5D tensors"};

    const auto & tensor_shapes = tensor.sizes();
    for (int i = 1; i < 5; ++i)
    {
      if (tensor_shapes[i] != 3)
        throw material_error{
            "Tensor4ToVoigt only works for 5D tensors with nx3x3x3x3 shapes"};
    }
    // if (!torch::allclose(tensor, torch::transpose(tensor, 1, 2), 1E-4) ||
    //     !torch::allclose(tensor, torch::transpose(tensor, 2, 3), 1E-4))
    //{
    //   // throw material_error{"Tensor4ToVoigt only works for symmetric
    //   // tensors"};
    //   //std::cerr << "Warning: Tensor4ToVoigt called with non symmetric
    //   tensor\n";
    // }
    auto voigt_tensor = torch::zeros({tensor_shapes[0], 6, 6});
    for (int i = 0; i < 6; ++i)
    {
      for (int j = 0; j < 6; ++j)
      {
        voigt_tensor.index({Slice{None, None}, i, j}) = tensor.index(
            {Slice{None, None}, Tensor4Index(i).first, Tensor4Index(i).second,
             Tensor4Index(j).first, Tensor4Index(j).second});
      }
    }
    return voigt_tensor;
  }

  [[nodiscard]] static torch::Tensor VoigtToTensor4(const torch::Tensor & voigt)
  {
    using torch::indexing::None;
    using torch::indexing::Slice;
    if (voigt.ndimension() != 3)
      throw material_error{"VoigtToTensor4 only works for 3D tensors"};
    const auto & voigt_shapes = voigt.sizes();
    for (int i = 1; i < 3; ++i)
    {
      if (voigt_shapes[i] != 6)
        throw material_error{
            "VoigtToTensor4 only works for 3D tensors with nx6x6 shapes"};
    }
    auto tensor4 = torch::zeros({voigt_shapes[0], 3, 3, 3, 3});
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        auto idx1 = VoigtIndex(i, j);
        for (int k = 0; k < 3; ++k)
        {
          for (int l = 0; l < 3; ++l)
          {
            auto idx2 = VoigtIndex(k, l);
            tensor4.index({Slice{None, None}, i, j, k, l}) =
                voigt.index({Slice{None, None}, idx1, idx2});
          }
        }
      }
    }
    return tensor4;
  }

  [[nodiscard]] static torch::Tensor TotalLagrangianStiffnessToULStiffness(
      const torch::Tensor & material_stiffness,
      const torch::Tensor & deformation_gradient)
  {
    using torch::indexing::None;
    using torch::indexing::Slice;
    auto Jinv = 1.0 / torch::linalg::det(deformation_gradient);
    auto tensor = torch::einsum(
        "kmi, knj, kijrs, kpr, kqs -> kmnpq",
        {deformation_gradient, deformation_gradient, material_stiffness,
         deformation_gradient, deformation_gradient});
    auto ul_material_stiffness =
        Jinv.index({Slice{None, None}, None, None}) * Tensor4ToVoigt(tensor);
    return ul_material_stiffness;
  }

  [[nodiscard]] static at::Tensor MandelProbeInv(int batch_size)
  {
    auto T = torch::zeros({batch_size, 6, 6});
    auto sr2 = -sqrt(2.0) / 2.0;
    for (int i = 0; i < batch_size; ++i)
    {
      T.index({i, 0, 0}) = 1;
      T.index({i, 1, 1}) = 1;
      T.index({i, 2, 2}) = 1;
      T.index({i, 3, 3}) = 1. / 2.;
      T.index({i, 4, 4}) = 1. / 2.;
      T.index({i, 5, 5}) = 1. / 2.;
      T.index({i, 0, 4}) = sr2;
      T.index({i, 0, 5}) = sr2;
      T.index({i, 1, 3}) = sr2;
      T.index({i, 1, 5}) = sr2;
      T.index({i, 2, 3}) = sr2;
      T.index({i, 2, 4}) = sr2;
    }
    return T;
  }

  [[nodiscard]] static at::Tensor ProbeDirections(int batch_size)
  {
    auto directions = torch::zeros({6, batch_size, 3, 3});
    for (int i = 0; i < batch_size; ++i)
    {
      directions.index({0, i, 0, 0}) = 1;
      directions.index({1, i, 1, 1}) = 1;
      directions.index({2, i, 2, 2}) = 1;
      directions.index({3, i, 1, 1}) = 1;
      directions.index({3, i, 1, 2}) = 1;
      directions.index({3, i, 2, 1}) = 1;
      directions.index({3, i, 2, 2}) = 1;
      directions.index({4, i, 0, 0}) = 1;
      directions.index({4, i, 0, 2}) = 1;
      directions.index({4, i, 2, 0}) = 1;
      directions.index({4, i, 2, 2}) = 1;
      directions.index({5, i, 0, 0}) = 1;
      directions.index({5, i, 0, 1}) = 1;
      directions.index({5, i, 1, 0}) = 1;
      directions.index({5, i, 1, 1}) = 1;
    }
    return directions;
  }

  static void MandelToVoigt(torch::Tensor & tensor)
  {
    using torch::indexing::None;
    using torch::indexing::Slice;

    tensor.index({Slice{None, None}, Slice{3, None}, Slice{0, 3}}) /= sqrt(2);
    tensor.index({Slice{None, None}, Slice{0, 3}, Slice{3, None}}) /= sqrt(2);
    tensor.index({Slice{None, None}, Slice{3, None}, Slice{3, None}}) /= 2;
  }

  static void VoigtToMandel(torch::Tensor & tensor)
  {
    using torch::indexing::None;
    using torch::indexing::Slice;

    tensor.index({Slice{None, None}, Slice{3, None}, Slice{0, 3}}) *= sqrt(2);
    tensor.index({Slice{None, None}, Slice{0, 3}, Slice{3, None}}) *= sqrt(2);
    tensor.index({Slice{None, None}, Slice{3, None}, Slice{3, None}}) *= 2;
  }

  [[nodiscard]] static torch::Tensor ComputeMaterialStiffness(
      const torch::Tensor & PK2,
      const torch::Tensor & right_cauchy_green)
  {
    using torch::indexing::None;
    using torch::indexing::Slice;
    int batch_size = (int)PK2.size(0);
    auto material_stiffness = torch::empty({batch_size, 6, 6});
    auto probe_directions = ProbeDirections(batch_size);
    auto mandel_probe_inv = MandelProbeInv(batch_size);
    std::vector<int> dir1_vec{0, 1, 2, 1, 0, 0};
    std::vector<int> dir2_vec{0, 1, 2, 2, 2, 1};
    for (int i = 0; i < 6; ++i)
    {
      auto grad = torch::autograd::grad({PK2}, {right_cauchy_green},
                                        {probe_directions[i]}, true, true)[0];
      assert(torch::allclose(grad, torch::transpose(grad, 1, 2), 1E-4));
      for (int j = 0; j < 6; ++j)
      {
        material_stiffness.index({Slice{None, None}, j, i}) =
            grad.index({Slice{None, None}, dir1_vec[j], dir2_vec[j]});
      }
    }
    VoigtToMandel(material_stiffness);
    //// use probe inverse to get true derivatives (in Mandel form)
    material_stiffness = torch::matmul(material_stiffness, mandel_probe_inv);
    MandelToVoigt(material_stiffness);
    // at this point we have the Voigt form of dPK2/dC
    //  dPK2/dE = 2*dPK2/dC (d^2u/dEdE = 4*d^2u/dCdC)...take derivative w.r.t.
    //  PK2
    material_stiffness *= 2.0;
    return material_stiffness;
  }

  [[nodiscard]] static torch::Tensor PK2ToCauchy(const torch::Tensor & PK2,
                                                 const torch::Tensor & F)
  {
    using torch::indexing::None;
    using torch::indexing::Slice;
    auto Jinv = 1.0 / torch::linalg::det(F);
    auto sigma = torch::einsum("k, kip,kpl,kjl->kij", {Jinv, F, PK2, F});
    auto sigma_reduced = torch::empty({PK2.sizes()[0], 6});
    for (int i = 0; i < 6; ++i)
    {
      auto idx = Tensor4Index(i);
      sigma_reduced.index({Slice{None, None}, i}) =
          sigma.index({Slice{None, None}, idx.first, idx.second});
    }
    return sigma_reduced;
  }

  struct DerivativeResult
  {
    at::Tensor cauchy_stress;
    at::Tensor material_stiffness;
  };

  [[nodiscard]] static DerivativeResult ComputeCauchyStressStiffness(
      const at::Tensor & energy,
      const at::Tensor & deformation_gradient,
      const at::Tensor & right_cauchy_green)
  {
    DerivativeResult result;
    auto energy_sum = torch::sum(energy, 0);
    auto PK2 = 2 * torch::autograd::grad({energy_sum}, {right_cauchy_green}, {},
                                         true, true)[0];
    result.cauchy_stress = PK2ToCauchy(PK2, deformation_gradient);
    auto material_stiffness = ComputeMaterialStiffness(PK2, right_cauchy_green);
    auto val_tensor = VoigtToTensor4(material_stiffness);
    result.material_stiffness =
        TotalLagrangianStiffnessToULStiffness(val_tensor, deformation_gradient);
    return result;
  }

  static void ApfMatrixToTorchTensor(const apf::Matrix3x3 & apf_matrix,
                                     torch::Tensor & tensor)
  {
    // accessor
    auto tensor_a = tensor.accessor<float, 3>();
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; j++)
      {
        // tensor.index({i, j}) = apf_matrix[i][j];
        tensor_a[0][i][j] = static_cast<float>(apf_matrix[i][j]);
      }
    }
  }

  [[nodiscard]] static apf::DynamicMatrix StressToApfDynamicMatrix(
      const torch::Tensor & cauchy_stress)
  {
    auto cauchy_stress_i = cauchy_stress.squeeze();
    if (cauchy_stress_i.ndimension() != 1 || cauchy_stress_i.sizes()[0] != 6)
    {
      throw material_error(
          "cauchy stress tensor should be 1 dimensional with 6 elements");
    }
    auto cauchy_stress_a = cauchy_stress_i.accessor<float, 1>();
    apf::DynamicMatrix apf_matrix(3, 3);
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        apf_matrix(i, j) =
            static_cast<double>(cauchy_stress_a[VoigtIndex(i, j)]);
      }
    }
    return apf_matrix;
  }

  [[nodiscard]] static apf::DynamicMatrix TorchTensorToApfDynamicMatrix(
      const torch::Tensor & tensor_in)
  {
    auto tensor = tensor_in.squeeze();
    if (tensor.ndimension() != 2)
    {
      throw material_error("Tensor must be 2 dimensional");
    }
    auto nrow = tensor.sizes()[0];
    auto ncol = tensor.sizes()[0];
    apf::DynamicMatrix apf_matrix(nrow, ncol);
    // We don't directly do a memcopy because the tensor may have a float type
    // whereas the apf matrix is a double type.
    auto tensor_a = tensor.accessor<float, 2>();
    for (int i = 0; i < nrow; ++i)
    {
      for (int j = 0; j < ncol; ++j)
      {
        apf_matrix(i, j) = static_cast<double>(tensor_a[i][j]);
      }
    }
    return apf_matrix;
  }

  TorchMaterial::TorchMaterial(const torch::jit::Module & Model) : model(Model)
  {
    try
    {
      deformation_gradient_tensor = torch::empty({1, 3, 3});
      inputs.emplace_back(torch::empty({1, 3, 3}));
    }
    // forward torch exception as
    catch (const torch::Error & e)
    {
      throw material_error(fmt::format("mumfim TorchMaterial: {}", e.what()));
    }
  }

  TorchMaterial::TorchMaterial(const std::string & model_path)
  {
    try
    {
      model = torch::jit::load(model_path);
      deformation_gradient_tensor = torch::empty({1, 3, 3});
      inputs.emplace_back(torch::empty({1, 3, 3}));
    }
    // forward torch exception as
    catch (const torch::Error & e)
    {
      throw material_error(fmt::format("mumfim TorchMaterial: {}", e.what()));
    }
  }

  MaterialResult TorchMaterial::operator()(
      const apf::Matrix3x3 & deformation_gradient,
      apf::MeshEntity * /* element */,
      int /* integration_point */)
  {
    try
    {
      ApfMatrixToTorchTensor(deformation_gradient, deformation_gradient_tensor);
      // we set requires_grad to false because we do in place tensor operations
      inputs[0].toTensor().set_requires_grad(false);
      ComputeRightCauchyGreen(deformation_gradient_tensor,
                              inputs[0].toTensor());
      // set requires_grad to true so we can take derivatives w.r.t. right
      // cauchy green deformation tensor
      inputs[0].toTensor().set_requires_grad(true);
      auto energy = model.forward(inputs).toTensor();
      auto [cauchy_stress, stiffness] = ComputeCauchyStressStiffness(
          energy, deformation_gradient_tensor, inputs[0].toTensor());

      MaterialResult result;
      // stress is 6x1 for symmetric components, but material result
      // expects the full 3x3 matrix
      result.cauchy_stress = StressToApfDynamicMatrix(cauchy_stress);
      result.material_stiffness = TorchTensorToApfDynamicMatrix(stiffness);
      return result;
    }
    // forward torch exception as
    catch (const torch::Error & e)
    {
      throw material_error(fmt::format("mumfim TorchMaterial: {}", e.what()));
    }
  }

}  // namespace mumfim
