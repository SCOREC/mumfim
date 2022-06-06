#ifndef MUMFIM_SRC_MUMFIM_MACROSCALE_UPDATEDLAGRANGIANMATERIALINTEGRATOR_H
#define MUMFIM_SRC_MUMFIM_MACROSCALE_UPDATEDLAGRANGIANMATERIALINTEGRATOR_H
#include <ElementalSystem.h>
#include <amsiDeformation.h>
#include <mumfim/exceptions.h>
#include "materials/Materials.h"
#include "ApfMatrixUtils.h"
#include <sstream>
namespace mumfim
{
  template <typename MaterialConstitutiveFunction>
  class UpdatedLagrangianMaterialIntegrator : public amsi::ElementalSystem
  {
    static_assert(std::is_invocable_r_v<MaterialResult,
                                        MaterialConstitutiveFunction,
                                        apf::Matrix3x3, apf::MeshEntity*, int>,
                  "Material must take deformation gradient and return stress "
                  "and material stiffness.");

    public:
    UpdatedLagrangianMaterialIntegrator(MaterialConstitutiveFunction constitutive_function,
                                        apf::Field * strn,
                                        apf::Field * strs,
                                        apf::Field * u,
                                        apf::Field * dfm_grd,
                                        int o)
        : ElementalSystem(u, o)
        , current_integration_point(0)
        , strain_field(strn)
        , stress_field(strs)
        , dim(-1)
        , ref_lmnt(nullptr)
        , du_lmnt(nullptr)
        , dfm_grd_fld(dfm_grd)
        , constitutive_function_(std::move(constitutive_function))
    {
    }
    void inElement(apf::MeshElement * me) override
    {
      ElementalSystem::inElement(me);
      ref_lmnt =
          apf::createMeshElement(apf::getMesh(f), apf::getMeshEntity(me));
      du_lmnt = apf::createElement(f, ref_lmnt);
      dim = apf::getDimension(me);
      current_integration_point = 0;
    }
    void outElement() override
    {
      apf::destroyElement(du_lmnt);
      apf::destroyMeshElement(ref_lmnt);
      ElementalSystem::outElement();
    }
    bool includesBodyForces() override { return true; }
    // after the refactor, the micro-scale and macro-scale stress tesnors are
    // indexed identically Note that the system is processed using the current
    // coordinate field. I.E. the volume, and gradients of the default element
    // are updated to the current iteration
    void atPoint(apf::Vector3 const & p, double w, double dV) override
    {
      if(dV <= 0) {
        std::stringstream ss;
        ss <<"Elemement has negative Jacobian! (" << dV <<")";
        throw mumfim_error(ss.str());
      }
      apf::MeshEntity * mesh_entity = apf::getMeshEntity(me);
      // Compute the deformation gradient and strain tensor
      apf::Matrix3x3 F;
      amsi::deformationGradient(du_lmnt, p, F);
      auto detF = apf::getDeterminant(F);
      if(detF<=0) {
        std::stringstream ss;
        ss <<"Elemement has negative deformation gradient! (" << detF <<")\n";
        ss<<"F=\n"<<F;
        throw mumfim_error(ss.str());
      }
      const MaterialResult consitutive_response =
          constitutive_function_(F, mesh_entity, current_integration_point);
      const auto & material_stiffness = consitutive_response.material_stiffness;
      const auto & cauchy_stress = consitutive_response.cauchy_stress;
      int & nen = nenodes;   // = 4 (tets)
      int & nedof = nedofs;  // = 12 (tets)
      // fill the matrix form of the stress tensor
      apf::NewArray<apf::Vector3> grads;
      apf::getShapeGrads(e, p, grads);
      apf::NewArray<double> N;
      apf::getShapeValues(e, p, N);
      apf::DynamicMatrix BL(6, nedof);  // linear strain disp
      BL.zero();
      for (int ii = 0; ii < nen; ii++)
      {
        BL(0, dim * ii) = grads[ii][0];      // N_(ii,1)
        BL(1, dim * ii + 1) = grads[ii][1];  // N_(ii,2)
        BL(2, dim * ii + 2) = grads[ii][2];  // N_(ii,3)
        BL(5, dim * ii) = grads[ii][1];      // N_(ii,2)
        BL(5, dim * ii + 1) = grads[ii][0];  // N_(ii,1)
        BL(3, dim * ii + 1) = grads[ii][2];  // N_(ii,3)
        BL(3, dim * ii + 2) = grads[ii][1];  // N_(ii,2)
        BL(4, dim * ii) = grads[ii][2];      // N_(ii,3)
        BL(4, dim * ii + 2) = grads[ii][0];  // N_(ii,1)
      }
      apf::DynamicMatrix K0(nedof, nedof);
      apf::DynamicMatrix tmp(6, nedof);
      apf::DynamicMatrix BLT(nedof, 6);
      apf::transpose(BL, BLT);
      apf::multiply(material_stiffness, BL, tmp);
      apf::multiply(BLT, tmp, K0);
      // Nonlinear terms (geometric stiffness)
      /*=======================================================================
        Initial Stress (or Geometric) component of Tangent Stiffness Matrix: K1
        - See Bathe PP.557, Table 6.6 (Updated Lagrangian)
        ========================================================================*/
      apf::DynamicMatrix BNL(9, nedof);  // nonlinear strain disp
      BNL.zero();
      for (int i = 0; i < nenodes; ++i)
      {
        BNL(0, i * dim) = grads[i][0];
        BNL(1, i * dim) = grads[i][1];
        BNL(2, i * dim) = grads[i][2];
        BNL(3, i * dim + 1) = grads[i][0];
        BNL(4, i * dim + 1) = grads[i][1];
        BNL(5, i * dim + 1) = grads[i][2];
        BNL(6, i * dim + 2) = grads[i][0];
        BNL(7, i * dim + 2) = grads[i][1];
        BNL(8, i * dim + 2) = grads[i][2];
      }
      apf::DynamicMatrix tau(9, 9);
      tau.zero();
      for (int i = 0; i < dim; ++i)
      {
        for (int j = 0; j < dim; ++j)
        {
          tau(i, j) = cauchy_stress(i, j);
          tau(i + dim, j + dim) = cauchy_stress(i, j);
          tau(i + 2 * dim, j + 2 * dim) = cauchy_stress(i, j);
        }
      }
      apf::DynamicMatrix K1(nedof, nedof);
      K1.zero();
      apf::DynamicMatrix BNLT(nedof, 9);
      apf::DynamicMatrix BNLTxtau(9, nedof);
      apf::transpose(BNL, BNLT);
      apf::multiply(BNLT, tau, BNLTxtau);
      apf::multiply(BNLTxtau, BNL, K1);
      // for force vector calculation
      apf::DynamicMatrix BLTxSV(nedof, 1);
      BLTxSV.zero();
      // multiply linear portion of B with the Cauchy stress
      // apf::DynamicVector cauchyVoigt(num_symmetric_components(dim));
      auto cauchyVoigt = dynamic_matrix_to_voigt(cauchy_stress);
      apf::DynamicVector BLTxCauchyVoigt(nedof);
      apf::multiply(BLT, cauchyVoigt, BLTxCauchyVoigt);
      for (int ii = 0; ii < nedof; ii++)
      {
        fe[ii] -= BLTxCauchyVoigt(ii) * w * dV;  // P - F
        for (int jj = 0; jj < nedof; jj++)
        {
          Ke(ii, jj) += (K0(ii, jj) + K1(ii, jj)) * w * dV;
        }
      }
      // Calculate rightCauchyGreen tensor.
      auto green_strain = computeGreenLagrangeStrain(F);
      apf::setMatrix(dfm_grd_fld, mesh_entity, current_integration_point, F);
      apf::setMatrix(strain_field, mesh_entity, current_integration_point, green_strain);
      apf::Matrix3x3 cauchy_stress_matrix;
      for(size_t i=0; i<cauchy_stress.getRows(); ++i) {
        for(size_t j=0; j<cauchy_stress.getColumns(); ++j) {
          cauchy_stress_matrix[i][j] = cauchy_stress(i,j);
        }
      }
      apf::setMatrix(stress_field, mesh_entity, current_integration_point, cauchy_stress_matrix);
      current_integration_point++;
    }
    int current_integration_point;

    private:
    apf::Field * strain_field;
    apf::Field * stress_field;
    int dim;
    apf::MeshElement * ref_lmnt;
    apf::Element * du_lmnt;
    apf::Field * dfm_grd_fld;
    MaterialConstitutiveFunction constitutive_function_;
  };
}  // namespace mumfim
#endif  // MUMFIM_SRC_MUMFIM_MACROSCALE_UPDATEDLAGRANGIANMATERIALINTEGRATOR_H
