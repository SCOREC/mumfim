#include "NonlinearTissue.h"
#include <amsiControlService.h>
#include <apfFunctions.h>
#include <apfLabelRegions.h>
#include <array>
#include <cstdlib>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include "NeoHookeanIntegrator.h"
#include "TrnsIsoNeoHookeanIntegrator.h"
#include "gmi.h"
namespace mumfim
{
  NonlinearTissue::NonlinearTissue(apf::Mesh * mesh,
                                   const mt::CategoryNode & analysis_case,
                                   MPI_Comm cm)
      : TissueBase(mesh, analysis_case, {}, {}, "macro", cm)
      , constitutives()
      , dv_prev(0.0)
      , load_step(0)
      , iteration(0)
  {
    apf_primary_field =
        apf::createLagrangeField(apf_mesh, "displacement", apf::VECTOR, 1);
    apf::zeroField(apf_primary_field);
    delta_u = apf::createLagrangeField(apf_mesh, "displacement_delta",
                                       apf::VECTOR, 1);
    apf::zeroField(delta_u);
    apf_primary_delta_field = delta_u;
    // create a current coordinate field from the CurrentCoordFunc (x = X+u)
    xpyfnc =
        new amsi::XpYFunc(apf_mesh->getCoordinateField(), apf_primary_field);
    // user field is not zeroed because it is the sum of two other fields
    current_coords =
        apf::createUserField(apf_mesh, "current_coordinates", apf::VECTOR,
                             apf::getLagrange(1), xpyfnc);
    // take the current coordinate field, and subtract increment in displacement
    // to get the coords at the previous increment (this is sorta hacky but it
    // should work)
    prv_crd_fnc = new amsi::XpYFunc(current_coords, delta_u, 1, -1);
    prev_coords =
        apf::createUserField(apf_mesh, "previous_coordinates", apf::VECTOR,
                             apf::getLagrange(1), prv_crd_fnc);
    apf_primary_numbering = apf::createNumbering(apf_primary_field);
    strs = apf::createIPField(apf_mesh, "stress", apf::MATRIX, 1);
    apf::zeroField(strs);
    strn = apf::createIPField(apf_mesh, "strain", apf::MATRIX, 1);
    apf::zeroField(strn);
    dfm_grd = apf::createIPField(apf_mesh, "F", apf::MATRIX, 1);
    apf::zeroField(dfm_grd);
    stf_vrtn =
        apf::createIPField(apf_mesh, "stiffness_variation", apf::SCALAR, 1);
    apf::zeroField(stf_vrtn);
    axl_yngs_mod =
        apf::createIPField(apf_mesh, "axial_youngs_modulus", apf::SCALAR, 1);
    apf::zeroField(axl_yngs_mod);
    amsi::applyUniqueRegionTags(apf_mesh);
    // zero stiffness_variation and axial_youngs_modulus fields
    const auto * stiffness_gradient =
        problem_definition.unassociated->FindCategoryNodeByType(
            "stiffness gradient");
    if (stiffness_gradient != nullptr)
    {
      for (const auto & grad_nd : stiffness_gradient->GetCategoryNodes())
      {
        stf_vrtn_cnst.push_back(
            buildStiffnessVariation(grad_nd, stf_vrtn));
      }
    }
    for (auto & cnst : stf_vrtn_cnst)
    {
      cnst->populate_stf_vrtn_fld();
    }
    const auto * incompressible =
        problem_definition.unassociated->FindCategoryNodeByType(
            "incompressible");
    if (incompressible != nullptr)
    {
      for (const auto & incompressible_constraint :
           incompressible->GetCategoryNodes())
      {
        vol_cnst.push_back(buildVolumeConstraint(incompressible_constraint,
                                                 apf_primary_numbering));
      }
    }
    static constexpr int dimension = 3;
    auto * gmodel = mesh->getModel();
    struct gmi_ent * gent;
    auto * it = gmi_begin(gmodel, dimension);
    while ((gent = gmi_next(gmodel, it)))
    {
      int tag = gmi_tag(gmodel, gent);
      const auto * model_traits =
          problem_definition.associated->Find({dimension, tag});
      const auto * material_model =
          mt::GetCategoryByType(model_traits, "material model");
      const auto * continuum_model =
          mt::GetPrimaryCategoryByType(material_model, "continuum model");
      const auto * youngs_modulus =
          mt::GetCategoryModelTraitByType<mt::ScalarMT>(continuum_model,
                                                        "youngs modulus");
      const auto * poisson_ratio =
          mt::GetCategoryModelTraitByType<mt::ScalarMT>(continuum_model,
                                                        "poisson ratio");
      if (youngs_modulus == nullptr || poisson_ratio == nullptr)
      {
        std::cerr << " \"poisson ratio\" and \"youngs modulus\" are required "
                     "types for the continuum model.\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      if (continuum_model->GetType() == "isotropic_neohookean")
      {
        // should check to make sure the continuum model is iso lin ela for init
        // solve?
        constitutives[reinterpret_cast<apf::ModelEntity *>(gent)] =
            std::make_unique<NeoHookeanIntegrator>(
                this, apf_primary_field, dfm_grd, current_coords,
                (*youngs_modulus)(), (*poisson_ratio)(), 1);
      }
      else if (continuum_model->GetType() == "transverse_isotropic")
      {
        const auto * axis = mt::GetCategoryModelTraitByType<mt::VectorMT>(
            continuum_model, "axis");
        const auto * axial_shear_modulus =
            mt::GetCategoryModelTraitByType<mt::ScalarMT>(
                continuum_model, "axial shear modulus");
        const auto * axial_youngs_modulus =
            mt::GetCategoryModelTraitByType<mt::ScalarMT>(
                continuum_model, "axial youngs modulus");
        std::array<double, 3> axs = {(*axis)(0), (*axis)(1), (*axis)(2)};
        constitutives[reinterpret_cast<apf::ModelEntity *>(gent)] =
            std::make_unique<TrnsIsoNeoHookeanIntegrator>(
                this, apf_primary_field, stf_vrtn, dfm_grd, axl_yngs_mod,
                current_coords, (*youngs_modulus)(), (*poisson_ratio)(),
                &axs[0], (*axial_shear_modulus)(), (*axial_youngs_modulus)(),
                1);
      }
    }
    gmi_end(gmodel, it);
    dirichlet_bcs.push_back(
        amsi::DirichletBCEntry{.categories = {"displacement"}, .mt_name = "x"});
    dirichlet_bcs.push_back(
        amsi::DirichletBCEntry{.categories = {"displacement"}, .mt_name = "y"});
    dirichlet_bcs.push_back(
        amsi::DirichletBCEntry{.categories = {"displacement"}, .mt_name = "z"});
    neumann_bcs.push_back(
        amsi::NeumannBCEntry{.categories = {"pressure"},
                             .mt_name = "magnitude",
                             .mt_type = amsi::NeumannBCType::pressure});
    neumann_bcs.push_back(
        amsi::NeumannBCEntry{.categories = {"traction"},
                             .mt_name = "direction",
                             .mt_type = amsi::NeumannBCType::traction});
  }
  NonlinearTissue::~NonlinearTissue()
  {
    delete xpyfnc;
    apf::destroyField(current_coords);
    delete prv_crd_fnc;
    apf::destroyField(prev_coords);
    // destroying this field currently throws an error
    // most likely because it is zeroed and core does not
    // properly deal with this...
    // apf::destroyField(delta_u);
    apf::destroyField(strs);
    apf::destroyField(strn);
    apf::destroyField(stf_vrtn);
    apf::destroyField(axl_yngs_mod);
  }
  void NonlinearTissue::computeInitGuess(amsi::LAS * las)
  {
    // LinearTissue lt(model, mesh, prob_def, solution_strategy, analysis_comm);
    LinearTissue lt(apf_mesh, problem_definition, solution_strategy, output,
                    analysis_comm);
    lt.setSimulationTime(T);
    LinearSolver(&lt, las);
    las->iter();
    apf::copyData(delta_u, lt.getField());
    apf::copyData(apf_primary_field, lt.getField());
  }
  void NonlinearTissue::step()
  {
    for (auto cnst = vol_cnst.begin(); cnst != vol_cnst.end(); cnst++)
      (*cnst)->step();
    iteration = 0;
    load_step++;
  }
  void NonlinearTissue::iter()
  {
    for (auto cnst = vol_cnst.begin(); cnst != vol_cnst.end(); cnst++)
      (*cnst)->iter();
    iteration++;
  }
  void NonlinearTissue::Assemble(amsi::LAS * las)
  {
    ApplyBC_Neumann(las);
    // custom iterator would be perfect for switching for multiscale version
    // FIXME !!! here createMeshELement is using the mesh (Original coords),
    // BUT Constitutives Should be set to Use the deformed coords for consistency
    // with the multiscale analysis
    AssembleIntegratorIntoLAS(las, [this](apf::MeshEntity * me, int)
    { return constitutives[apf_mesh->toModel(me)].get(); }, current_coords);
    double nrm = 0.0;
    las->GetVectorNorm(nrm);
    // process constraints
    for (auto cnst = vol_cnst.begin(); cnst != vol_cnst.end(); cnst++)
      (*cnst)->apply(las);
    // las->GetVectorNorm(nrm);
  }
  void NonlinearTissue::UpdateDOFs(const double * sol)
  {
    // accumulate displacement deltas into primary field
    amsi::AccumOp ac_op;
    amsi::FreeApplyOp frac_op(apf_primary_numbering, &ac_op);
    amsi::ApplyVector(apf_primary_numbering, apf_primary_field, sol,
                      first_local_dof, &frac_op)
        .run();
    // amsi::PrintField(apf_primary_field,std::cout).run();
    amsi::WriteOp wr_op;
    amsi::FreeApplyOp frwr_op(apf_primary_numbering, &wr_op);
    amsi::ApplyVector(apf_primary_numbering, delta_u, sol, first_local_dof,
                      &frwr_op)
        .run();
    // amsi::PrintField(delta_u,std::cout).run();
    apf::synchronize(apf_primary_field);
    apf::synchronize(delta_u);
  }
  /**
   * get the Neumann boundary condition value on the specified entity
   * @param ent model entity to query
   * @param frc will be filled with the value of the NeumannBCs
   */
  void NonlinearTissue::getLoadOn(apf::ModelEntity * ent, double * frc)
  {
    auto model_dimension = apf_mesh->getModelType(ent);
    auto model_tag = apf_mesh->getModelTag(ent);
    auto * associated_traits =
        problem_definition.associated->Find({model_dimension, model_tag});
    if (associated_traits == nullptr)
    {
      std::cerr << "Attempting to get Neumann BC on [" << model_dimension << ","
                << model_tag << "].";
      std::cerr << " No boundary conditions associated with that geometry.\n";
      exit(1);
    }
    for (const auto & neumann_bc : neumann_bcs)
    {
      const mt::AssociatedCategoryNode * category_node = associated_traits;
      for (const auto & category : neumann_bc.categories)
      {
        category_node = category_node->FindCategoryByType(category);
        if (category_node == nullptr)
        {
          break;
        }
      }
      if (category_node == nullptr)
      {
        std::cerr << "Invalid Neumann BC provided\n";
        exit(1);
      }
      const auto * load_trait =
          mt::GetCategoryModelTraitByType(category_node, neumann_bc.mt_name);
      const auto * const_load_trait = mt::MTCast<mt::VectorMT>(load_trait);
      const auto * time_function_load =
          mt::MTCast<mt::VectorFunctionMT<1>>(load_trait);
      if (const_load_trait != nullptr)
      {
        for (int i = 0; i < 3; ++i)
        {
          frc[i] = (*const_load_trait)(i);
        }
      }
      else if (time_function_load != nullptr)
      {
        for (int i = 0; i < 3; ++i)
        {
          frc[i] = (*time_function_load)(i, T);
        }
      }
      else
      {
        std::cerr << "Load must be either a scalar value, or function of time "
                     "only.\n";
        exit(1);
      }
    }
  }
  void NonlinearTissue::recoverSecondaryVariables(int /* unused load_step */ )
  {
    /*
    //#ifdef SCOREC
    int rnk = -1;
    MPI_Comm_rank(AMSI_COMM_SCALE, &rnk);
    std::stringstream fnm;
    fnm << amsi::fs->getResultsDir() << "/qlty.stp_" << load_step << ".rnk_"
        << rnk << ".dat";
    // analyze and print the quality of the elements
    apf::Field* qfld = amsi::analyzeMeshQuality(apf_mesh, apf_primary_field);
    std::ofstream file(fnm.str().c_str(), std::ofstream::out);
    amsi::PrintField(qfld, file).run();
    apf::destroyField(qfld);
    */
    //#endif
  }
  void NonlinearTissue::storeStrain(apf::MeshElement * me, double * strain)
  {
    apf::MeshEntity * m_ent = apf::getMeshEntity(me);
    apf::Matrix3x3 eps(strain[0], strain[3], strain[5], strain[3], strain[1],
                       strain[4], strain[5], strain[4], strain[2]);
    apf::setMatrix(strn, m_ent, 0, eps);
  }
  void NonlinearTissue::storeStrain(apf::MeshElement * me, apf::Matrix3x3 eps)
  {
    apf::MeshEntity * m_ent = apf::getMeshEntity(me);
    apf::setMatrix(strn, m_ent, 0, eps);
  }
  void NonlinearTissue::storeStress(apf::MeshElement * me, double * stress)
  {
    apf::MeshEntity * m_ent = apf::getMeshEntity(me);
    apf::Matrix3x3 sigma(stress[0], stress[3], stress[5], stress[3], stress[1],
                         stress[4], stress[5], stress[4], stress[2]);
    apf::setMatrix(strs, m_ent, 0, sigma);
  }
  void NonlinearTissue::storeStress(apf::MeshElement * me, apf::Matrix3x3 eps)
  {
    apf::MeshEntity * m_ent = apf::getMeshEntity(me);
    apf::setMatrix(strs, m_ent, 0, eps);
  }
}  // namespace mumfim
