#include "MultiscaleTissue.h"
#include <amsiControlService.h>  // amsi
#include <amsiDetectOscillation.h>
#include <model_traits/AssociatedModelTraits.h>
#include <map>
#include <string>
#include "mumfim/microscale/MicroFOParams.h"
#include "mumfim/macroscale/ULMultiscaleHydrostaticPressureIntegrator.h"
#include "mumfim/macroscale/ULMultiscaleIntegrator.h"
namespace mumfim
{
  struct StochasticFieldData
  {
    int orientation_num_bins;
    std::string orientation_filename;
    int alignment_num_bins;
    std::string alignment_filename;
  };
  static StochasticFieldData read_stochastic_field(
      const mt::AssociatedCategoryNode & stochastic_field)
  {
    const auto * alignment_field =
        stochastic_field.FindCategoryByType("alignment field");
    const auto * orientation_field =
        stochastic_field.FindCategoryByType("orientation field");
    if (alignment_field == nullptr || orientation_field == nullptr)
    {
      std::cerr << "\"alignment field\" and \"orientation field\" are "
                   "required for stochastic field definition.\n";
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    const auto * alignment_filename_cat =
        alignment_field->FindCategoryByType("filename");
    const auto * alignment_num_bins_cat =
        alignment_field->FindCategoryByType("number of bins");
    const auto * orientation_filename_cat =
        orientation_field->FindCategoryByType("filename");
    const auto * orientation_num_bins_cat =
        orientation_field->FindCategoryByType("number of bins");
    if (alignment_filename_cat == nullptr ||
        alignment_num_bins_cat == nullptr ||
        orientation_filename_cat == nullptr ||
        orientation_num_bins_cat == nullptr)
    {
      std::cerr
          << "\"filename\" and \"number of bins\" categories are required for "
             "\"orientation field\" and \"alignment field\".\n";
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    const auto * alignment_filename =
        mt::MTCast<mt::StringMT>(alignment_filename_cat->GetModelTrait());
    const auto * orientation_filename =
        mt::MTCast<mt::StringMT>(orientation_filename_cat->GetModelTrait());
    const auto * alignment_num_bins =
        mt::MTCast<mt::IntMT>(alignment_num_bins_cat->GetModelTrait());
    const auto * orientation_num_bins =
        mt::MTCast<mt::IntMT>(orientation_num_bins_cat->GetModelTrait());
    if (alignment_filename == nullptr || alignment_num_bins == nullptr ||
        orientation_filename == nullptr || orientation_num_bins == nullptr)
    {
      std::cerr
          << "\"filename\" and \"number of bins\" require model traits for "
             "\"orientation field\" and \"alignment field\".\n";
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    return StochasticFieldData{
        .orientation_num_bins = (*orientation_num_bins)(),
        .orientation_filename = (*orientation_filename)(),
        .alignment_num_bins = (*alignment_num_bins)(),
        .alignment_filename = (*alignment_filename)(),
    };
  }
  MultiscaleTissue::MultiscaleTissue(apf::Mesh * mesh,
                                     const mt::CategoryNode & analysis_case,
                                     MPI_Comm cm,
                                     const amsi::Multiscale & amsi_multiscale)
      : NonlinearTissue(mesh, analysis_case, cm)
      , mltscl(NULL)
      , crt_rve(apf::createIPField(apf_mesh, "micro_rve_type", apf::SCALAR, 1))
      , prv_rve(apf::createIPField(apf_mesh, "micro_old_type", apf::SCALAR, 1))
      , compute_ornt_3D(false)
      , compute_ornt_2D(false)
      , ornt_3D(NULL)
      , ornt_2D(NULL)
      , fo_cplg(apf_mesh,
                crt_rve,
                prv_rve,
                apf::getShape(apf_primary_field)->getOrder(),
                amsi_multiscale)
      , nm_rves(0)
      , rve_dirs()
      , multiscale_(amsi_multiscale)
  {
    ornt_3D = apf::createIPField(apf_mesh, "ornt_tens_3D", apf::MATRIX, 1);
    ornt_2D = apf::createIPField(apf_mesh, "ornt_tens_2D", apf::MATRIX, 1);
    mltscl = new ULMultiscaleIntegrator(&fo_cplg, strn, strs, apf_primary_field,
                                        dfm_grd, 1);
    M2m_id =
        amsi::getRelationID(multiscale_.getMultiscaleManager(),
                            multiscale_.getScaleManager(), "macro", "micro_fo");
    m2M_id =
        amsi::getRelationID(multiscale_.getMultiscaleManager(),
                            multiscale_.getScaleManager(), "micro_fo", "macro");
    apf::zeroField(ornt_3D);
    apf::zeroField(ornt_2D);
    apf::zeroField(crt_rve);
    apf::zeroField(prv_rve);
    // DEBUG
    test_inc_dfm = apf::createIPField(apf_mesh, "inc_F", apf::MATRIX, 1);
    // END DEBUG
    // get properties for output orientation tensor
    const auto * whole_model_traits = output.associated->Find({9, 1});
    if (whole_model_traits)
    {
      const auto * ornt_tens =
          whole_model_traits->FindCategoryByType("output orientation tensor");
      if (ornt_tens)
      {
        const auto * ornt_tens_3D =
            mt::GetCategoryModelTraitByType(ornt_tens, "3D Orientation Tensor");
        const auto * ornt_tens_2D =
            mt::GetCategoryByType(ornt_tens, "2D Orientation Tensor");
        compute_ornt_3D = (ornt_tens_3D != nullptr);
        compute_ornt_2D = (ornt_tens_2D != nullptr);
        if (ornt_tens_2D != nullptr)
        {
          const auto * axis = mt::GetCategoryModelTraitByType<mt::VectorMT>(
              ornt_tens_2D, "axis");
          if (axis == nullptr)
          {
            std::cerr << "2D orientation requires axis.\n";
            MPI_Abort(AMSI_COMM_WORLD, 1);
          }
          ornt_2D_axis[0] = (*axis)(0);
          ornt_2D_axis[1] = (*axis)(1);
          ornt_2D_axis[2] = (*axis)(2);
        }
      }
    }
  }
  MultiscaleTissue::~MultiscaleTissue()
  {
    delete mltscl;
    // destroying these fields causes segfult...needs investigation
    apf::destroyField(crt_rve);
    apf::destroyField(prv_rve);
    apf::destroyField(ornt_3D);
    apf::destroyField(ornt_2D);
  }
  void MultiscaleTissue::Assemble(amsi::LAS * las)
  {
    computeRVEs();
#ifdef LOGRUN
    amsi::Log state = amsi::activateLog("tissue_efficiency");
    amsi::log(state) << load_step << ", " << iteration << ", " << MPI_Wtime()
                     << ", "
                     << "start_fea" << std::endl;
#endif
    apfFEA::ApplyBC_Neumann(las);
    apf::MeshIterator * it = apf_mesh->begin(analysis_dim);
    apf::MeshEntity * me = NULL;
    while ((me = apf_mesh->iterate(it)))
    {
      apf::MeshElement * mlm = apf::createMeshElement(current_coords, me);
      amsi::ElementalSystem * sys =
          getIntegrator(me, 0);  // ERROR: assumes 1 type per ent
      apf::Element * elm = apf::createElement(sys->getField(), mlm);
      sys->process(mlm);
      apf::NewArray<apf::Vector3> dofs;
      apf::getVectorNodes(elm, dofs);
      apf::NewArray<int> ids;
      apf::getElementNumbers(apf_primary_numbering, me, ids);
      AssembleDOFs(las, sys->numElementalDOFs(), &ids[0], &dofs[0],
                   &sys->getKe()(0, 0), &sys->getfe()(0),
                   sys->includesBodyForces());
      apf::destroyElement(elm);
      apf::destroyMeshElement(mlm);
    }
    apf_mesh->end(it);
    for (auto cnst = vol_cnst.begin(); cnst != vol_cnst.end(); cnst++)
      (*cnst)->apply(las);
#ifdef LOGRUN
    amsi::log(state) << load_step << ", " << iteration << ", " << MPI_Wtime()
                     << ", "
                     << "end_fea" << std::endl;
    amsi::log(state) << load_step << ", " << iteration << ", " << MPI_Wtime()
                     << ", "
                     << "start_solve" << std::endl;
#endif
  }
  void MultiscaleTissue::computeRVEs()
  {
#ifdef LOGRUN
    amsi::Log state = amsi::activateLog("tissue_efficiency");
    amsi::log(state) << load_step << ", " << iteration << ", " << MPI_Wtime()
                     << ", "
                     << "start_rves" << std::endl;
#endif
    std::vector<micro_fo_data> fo_data;
    serializeRVEData(
        std::back_inserter(fo_data));  // serialize data for fiber_only
    fo_cplg.sendRVEData(fo_data);
    fo_cplg.recvRVEData();
#ifdef LOGRUN
    amsi::log(state) << load_step << ", " << iteration << ", " << MPI_Wtime()
                     << ", "
                     << "end_rves" << std::endl;
#endif
  }
  void MultiscaleTissue::initMicro()
  {
    amsi::ControlService * cs = amsi::ControlService::Instance();
    fo_cplg.initCoupling();
    loadRVELibraryInfo();
    int num_rve_tps = rve_dirs.size();
    cs->scaleBroadcast(M2m_id, &num_rve_tps);
    std::vector<int> rve_cnts(num_rve_tps);
    int ii = 0;
    std::vector<MPI_Request> rqsts;
    for (auto tp = rve_dirs.begin(); tp != rve_dirs.end(); ++tp)
    {
      cs->aSendBroadcast(std::back_inserter(rqsts), M2m_id, tp->c_str(),
                         tp->size() + 1);
      rve_cnts[ii++] = rve_dir_cnts[std::distance(rve_dirs.begin(), tp)];
    }
    cs->aSendBroadcast(std::back_inserter(rqsts), M2m_id, &rve_cnts[0],
                       num_rve_tps);
    // wait for the sends to complete if we don't do this it can cause
    // a crash on systems without hardware MPI buffers
    MPI_Waitall(rqsts.size(), &rqsts[0], MPI_STATUSES_IGNORE);
  }
  void MultiscaleTissue::updateMicro()
  {
    updateRVETypes();
    updateRVEExistence();
  }
  void MultiscaleTissue::updateRVETypes()
  {
    nm_rves = 0;
    apf::MeshIterator * it = apf_mesh->begin(analysis_dim);
    for (apf::MeshEntity * me = NULL; (me = apf_mesh->iterate(it));)
    {
      apf::MeshElement * mlm = apf::createMeshElement(apf_mesh, me);
      int ip = apf::countIntPoints(mlm, 1);
      for (int ii = 0; ii < ip; ++ii)
      {
        apf::setScalar(prv_rve, me, ii, apf::getScalar(crt_rve, me, ii));
        MicroscaleType nw_tp = updateRVEType(me);
        apf::setScalar(crt_rve, me, ii, static_cast<int>(nw_tp));
        if (nw_tp == MicroscaleType::FIBER_ONLY) nm_rves++;
      }
      apf::destroyMeshElement(mlm);
    }
    apf_mesh->end(it);
  }
  MicroscaleType MultiscaleTissue::updateRVEType(apf::MeshEntity * me)
  {
    // TODO: Use error estimate as determination of microscale TYPE
    // apf::MeshElement * me = apf::createMeshElement(apf_mesh,m_ent);
    // apf::Element * err_elmt = apf::createElement(apf_size_field,me);
    // get value of size field on element
    // if value is ?less? than some threshold, FIBER_ONLY, otherwise
    // constitutive
    //
    // multiscale based purely off of simmodeler specification, no adaptivity,
    // need differentiate between initialization and updating
    auto * gent = apf_mesh->toModel(me);
    const auto * pd = problem_definition.associated->Find(
        {apf_mesh->getModelType(gent), apf_mesh->getModelTag(gent)});
    const auto * material_model = pd->FindCategoryByType("material model");
    const auto * multiscale_model =
        material_model->FindCategoryByType("multiscale model");
    return getMicroscaleType(multiscale_model);
  }
  void MultiscaleTissue::updateRVEExistence()
  {
    std::vector<micro_fo_header> nw_hdrs;
    std::vector<micro_fo_params> nw_prms;
    std::vector<micro_fo_init_data> nw_data;
    std::vector<int> to_dlt;
    std::vector<micro_fo_solver> slvr_prms;
    std::vector<micro_fo_int_solver> slvr_int_prms;
    fo_cplg.updateRVEDeletion(std::back_inserter(to_dlt));
    serializeNewRVEData(
        std::back_inserter(nw_hdrs), std::back_inserter(nw_prms),
        std::back_inserter(nw_data), std::back_inserter(slvr_prms),
        std::back_inserter(slvr_int_prms));
    fo_cplg.deleteRVEs(to_dlt.begin(), to_dlt.end());
    size_t add_ptrn = fo_cplg.addRVEs(nw_hdrs.size());
    fo_cplg.sendNewRVEs(add_ptrn, nw_hdrs, nw_prms, nw_data, slvr_prms,
                        slvr_int_prms);
    fo_cplg.updateRecv();
    nm_rves += nw_hdrs.size();
  }
  void MultiscaleTissue::recoverSecondaryVariables(int step)
  {
    NonlinearTissue::recoverSecondaryVariables(step);
    fo_cplg.recvRVEStepData();
    // loop over all mesh regions
    apf::MeshIterator * it = apf_mesh->begin(analysis_dim);
    while (apf::MeshEntity * meshEnt = apf_mesh->iterate(it))
    {
      if (apf_mesh->isOwned(meshEnt))
      {
        apf::MeshElement * mlm = apf::createMeshElement(apf_mesh, meshEnt);
        int numIP = apf::countIntPoints(mlm, 1);
        for (int ip = 0; ip < numIP; ++ip)
        {
          micro_fo_step_result * stp_rslt =
              fo_cplg.getRVEStepResult(meshEnt, 0);
          // this is kindof ugly due to repetitions, but it takes out some
          // logic from the tight loop
          if (compute_ornt_2D && compute_ornt_3D)
          {
            apf::Matrix3x3 orn_tens_3D;
            apf::Matrix3x3 orn_tens_2D;
            for (int i = 0; i < 3; ++i)
            {
              for (int j = 0; j < 3; ++j)
              {
                orn_tens_3D[i][j] = stp_rslt->data[i * 3 + j];
                orn_tens_2D[i][j] = stp_rslt->data[9 + i * 3 + j];
              }
            }
            apf::setMatrix(ornt_3D, meshEnt, ip, orn_tens_3D);
            apf::setMatrix(ornt_2D, meshEnt, ip, orn_tens_2D);
          }
          else if (compute_ornt_3D)
          {
            apf::Matrix3x3 orn_tens_3D;
            for (int i = 0; i < 3; ++i)
            {
              for (int j = 0; j < 3; ++j)
              {
                orn_tens_3D[i][j] = stp_rslt->data[i * 3 + j];
              }
            }
            apf::setMatrix(ornt_3D, meshEnt, ip, orn_tens_3D);
          }
          else if (compute_ornt_2D)
          {
            apf::Matrix3x3 orn_tens_2D;
            for (int i = 0; i < 3; ++i)
            {
              for (int j = 0; j < 3; ++j)
              {
                orn_tens_2D[i][j] = stp_rslt->data[9 + i * 3 + j];
              }
            }
            apf::setMatrix(ornt_2D, meshEnt, ip, orn_tens_2D);
          }
        }
        apf::destroyMeshElement(mlm);
      }
    }
    apf_mesh->end(it);
    if (compute_ornt_2D) apf::synchronize(ornt_2D);
    if (compute_ornt_3D) apf::synchronize(ornt_3D);
  }
  amsi::ElementalSystem * MultiscaleTissue::getIntegrator(apf::MeshEntity * me,
                                                          int ip)
  {
    MicroscaleType tp =
        static_cast<MicroscaleType>(apf::getScalar(crt_rve, me, ip));
    switch (tp)
    {
      case MicroscaleType::NONE:
        return constitutives[apf_mesh->toModel(me)].get();
      case MicroscaleType::FIBER_ONLY:
        return mltscl;
      case MicroscaleType::ISOTROPIC_NEOHOOKEAN:
        return mltscl;
      default:
        std::cerr << "Applying a macroscale integrator. We should not expect "
                     "to get here.\n";
        return constitutives[apf_mesh->toModel(me)].get();
    }
  }
  void MultiscaleTissue::loadRVELibraryInfo()
  {
    apf::MeshEntity * ent;
    auto * it = apf_mesh->begin(3);
    while ((ent = apf_mesh->iterate(it)))
    {
      auto * gent = apf_mesh->toModel(ent);
      auto gdim = apf_mesh->getModelType(gent);
      auto gtag = apf_mesh->getModelTag(gent);
      const auto * pd = problem_definition.associated->Find({gdim, gtag});
      const auto * material_model = pd->FindCategoryByType("material model");
      if (material_model == nullptr)
      {
        std::cerr << "material model required on (" << gdim << "," << gtag
                  << ").\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      const auto * multiscale_model =
          material_model->FindCategoryByType("multiscale model");
      if (multiscale_model == nullptr)
      {
        std::cerr << "multiscale model required on (" << gdim << "," << gtag
                  << ").\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      auto micro_tp = getMicroscaleType(multiscale_model);
      multiscale_model = &multiscale_model->GetCategories()[0];
      const auto * stochastic_field =
          multiscale_model->FindCategoryByType("stochastic field");
      if (micro_tp == MicroscaleType::FIBER_ONLY)
      {
        const auto * directory_mt =
            mt::GetCategoryModelTraitByType<mt::StringMT>(multiscale_model,
                                                          "directory");
        const auto * prefix_mt = mt::GetCategoryModelTraitByType<mt::StringMT>(
            multiscale_model, "prefix");
        const auto * count_mt = mt::GetCategoryModelTraitByType<mt::IntMT>(
            multiscale_model, "count");
        if (directory_mt == nullptr || prefix_mt == nullptr ||
            count_mt == nullptr)
        {
          std::cerr << "\"directory\", \"prefix\", and \"count\" categories "
                       "are required components of the multiscale model.\n";
          MPI_Abort(AMSI_COMM_WORLD, 1);
        }
        std::string directory = (*directory_mt)();
        std::string prefix = (*prefix_mt)();
        auto tp = directory + "/";
        int count = (*count_mt)();
        if (stochastic_field)
        {
          auto stochastic_field_data = read_stochastic_field(*stochastic_field);
          auto search_alignment = stochastic_field_map.find(
              stochastic_field_data.alignment_filename);
          if (search_alignment == stochastic_field_map.end())
          {
            stochastic_field_map[stochastic_field_data.alignment_filename] =
                std::make_shared<GridData>(
                    stochastic_field_data.alignment_filename.c_str());
          }
          auto search_orientation = stochastic_field_map.find(
              stochastic_field_data.orientation_filename);
          if (search_orientation == stochastic_field_map.end())
          {
            stochastic_field_map[stochastic_field_data.orientation_filename] =
                std::make_shared<GridData>(
                    stochastic_field_data.orientation_filename.c_str());
          }
          auto orientation_field =
              stochastic_field_map[stochastic_field_data.orientation_filename];
          auto alignment_field =
              stochastic_field_map[stochastic_field_data.alignment_filename];
          auto centroid = apf::getLinearCentroid(apf_mesh, ent);
          tp += getNetworkSuffix(*alignment_field, *orientation_field, prefix,
                                 centroid[0], centroid[1], centroid[2],
                                 stochastic_field_data.alignment_num_bins,
                                 stochastic_field_data.orientation_num_bins);
        }
        else
        {
          tp += prefix;
        }
        if (std::find(rve_dirs.begin(), rve_dirs.end(), tp) == rve_dirs.end())
        {
          rve_dirs.push_back(tp);
          rve_dir_cnts.push_back(count);
        }
      }
    }
    apf_mesh->end(it);
  }
  void MultiscaleTissue::getExternalRVEData(apf::MeshEntity * ent,
                                            micro_fo_header & hdr,
                                            micro_fo_params & prm,
                                            micro_fo_solver & slvr,
                                            micro_fo_int_solver & int_slvr)
  {
    auto * model_entity = apf_mesh->toModel(ent);
    const auto * pd = problem_definition.associated->Find(
        {apf_mesh->getModelType(model_entity),
         apf_mesh->getModelTag(model_entity)});
    const auto * material_model = pd->FindCategoryByType("material model");
    const auto * multiscale_model =
        material_model->FindCategoryByType("multiscale model");
    // pGEntity smdl_ent = EN_whatIn(reinterpret_cast<pEntity>(ent));
    // pAttribute mm = GEN_attrib(smdl_ent, "material model");
    // pAttribute sm = Attribute_childByType(mm, "multiscale model");
    MicroscaleType micro_tp = getMicroscaleType(multiscale_model);
    // FIXME this should be static_cast to underlying type...add utility
    // to_underlying function
    hdr.data[RVE_TYPE] = static_cast<int>(micro_tp);
    if (micro_tp == MicroscaleType::FIBER_ONLY)
    {
      multiscale_model = multiscale_model->FindCategoryByType("fiber only");
      const auto * radius_cat = multiscale_model->FindCategoryByType("radius");
      const auto * volume_fraction_cat =
          multiscale_model->FindCategoryByType("volume fraction");
      if (radius_cat == nullptr || volume_fraction_cat == nullptr)
      {
        std::cerr
            << "\"radius\" and \"volume fraction\" are required categories.\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      // either linear or nonlinear
      const auto * force_reaction_cat =
          multiscale_model->FindCategoryByType("force reaction");
      if (force_reaction_cat == nullptr)
      {
        std::cerr << "\"force reaction\" category required in the multiscale "
                     "model.\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      force_reaction_cat = &force_reaction_cat->GetCategories()[0];
      if (force_reaction_cat == nullptr)
      {
        std::cerr << "\"force reaction\" must contain the reaction type "
                     "category (linear,nonlinear)\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      const auto * youngs_modulus_cat =
          force_reaction_cat->FindCategoryByType("youngs modulus");
      if (youngs_modulus_cat == nullptr)
      {
        std::cerr << "\"youngs modulus\" category is required for all "
                     "multiscale model types.\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      const auto * linear_transition_cat =
          force_reaction_cat->FindCategoryByType("linear transition");
      const auto * nonlinearity_parameter_cat =
          force_reaction_cat->FindCategoryByType("nonlinearity parameter");
      if (force_reaction_cat->GetType() == "nonlinear" &&
          (linear_transition_cat == nullptr ||
           nonlinearity_parameter_cat == nullptr))
      {
        std::cerr << "\"linear transition\" and \"nonlinearity parameter\" are "
                     "required on a nonlinear material model.\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      const auto * radius =
          mt::MTCast<mt::ScalarMT>(radius_cat->GetModelTrait());
      const auto * volume_fraction =
          mt::MTCast<mt::ScalarMT>(volume_fraction_cat->GetModelTrait());
      const auto * youngs_modulus =
          mt::MTCast<mt::ScalarMT>(youngs_modulus_cat->GetModelTrait());
      const auto * nonlinearity_parameter =
          nonlinearity_parameter_cat
              ? mt::MTCast<mt::ScalarMT>(
                    nonlinearity_parameter_cat->GetModelTrait())
              : nullptr;
      const auto * linear_transition =
          linear_transition_cat
              ? mt::MTCast<mt::ScalarMT>(linear_transition_cat->GetModelTrait())
              : nullptr;
      // get properties for output orienation tensor
      hdr.data[COMPUTE_ORIENTATION_3D] = compute_ornt_3D;
      hdr.data[COMPUTE_ORIENTATION_2D] = compute_ornt_2D;
      prm.data[ORIENTATION_AXIS_X] = ornt_2D_axis[0];
      prm.data[ORIENTATION_AXIS_Y] = ornt_2D_axis[1];
      prm.data[ORIENTATION_AXIS_Z] = ornt_2D_axis[2];
      prm.data[FIBER_RADIUS] = (*radius)();
      prm.data[VOLUME_FRACTION] = (*volume_fraction)();
      prm.data[YOUNGS_MODULUS] = (*youngs_modulus)();
      prm.data[NONLINEAR_PARAM] =
          nonlinearity_parameter ? (*nonlinearity_parameter)() : 0.0;
      prm.data[LINEAR_TRANSITION] =
          linear_transition ? (*linear_transition)() : 0.0;
    }
    else if (micro_tp == MicroscaleType::ISOTROPIC_NEOHOOKEAN)
    {
      multiscale_model =
          mt::GetCategoryByType(multiscale_model, "isotropic_neohookean");
      const auto * youngs_modulus =
          mt::GetCategoryModelTraitByType<mt::ScalarMT>(multiscale_model,
                                                        "youngs modulus");
      const auto * poisson_ratio =
          mt::GetCategoryModelTraitByType<mt::ScalarMT>(multiscale_model,
                                                        "poisson ratio");
      if (youngs_modulus == nullptr || poisson_ratio == nullptr)
      {
        std::cerr << "\"youngs modulus\" and \"poisson ratio\" are required "
                     "categories in an isotropic neohookean material.\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      prm.data[YOUNGS_MODULUS] = (*youngs_modulus)();
      // to avoid adding extra communication we use the nonlinear parameters
      // slot in the parameter pack
      prm.data[NONLINEAR_PARAM] = (*poisson_ratio)();
    }
    const auto * ss = solution_strategy.associated->GetNullGeometry();
    const auto * microscale_convergence =
        mt::GetPrimaryCategoryByType(ss, "microscale convergence operator");
    if (microscale_convergence == nullptr)
    {
      std::cerr << "\"microscale convergence operator\" must exist in the "
                   "solution strategy.\n";
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    const auto * solver_tolerance =
        mt::GetCategoryModelTraitByType<mt::ScalarMT>(microscale_convergence,
                                                      "micro solver tolerance");
    slvr.data[MICRO_SOLVER_TOL] =
        solver_tolerance ? (*solver_tolerance)() : 1E-6;
    if (microscale_convergence->GetType() == "explicit timestep")
    {
      // get explicit params
      int_slvr.data[MICRO_SOLVER_TYPE] = static_cast<int>(SolverType::Explicit);
      const auto * amplitude_cat =
          mt::GetPrimaryCategoryByType(microscale_convergence, "amplitude");
      const auto * viscous_damping_factor =
          mt::GetCategoryModelTraitByType<mt::ScalarMT>(
              microscale_convergence, "viscous damping factor");
      const auto * critical_time_factor =
          mt::GetCategoryModelTraitByType<mt::ScalarMT>(
              microscale_convergence, "critical time scale factor");
      const auto * energy_check_epsilon =
          mt::GetCategoryModelTraitByType<mt::ScalarMT>(microscale_convergence,
                                                        "energy check epsilon");
      const auto * serial_gpu_cutoff =
          mt::GetCategoryModelTraitByType<mt::IntMT>(
              microscale_convergence, "serial to GPU dof cutoff");
      if (amplitude_cat == nullptr || viscous_damping_factor == nullptr ||
          critical_time_factor == nullptr || energy_check_epsilon == nullptr ||
          serial_gpu_cutoff == nullptr)
      {
        std::cerr
            << "invalid explicit timestep microscale convergence operator.\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      int_slvr.data[SERIAL_GPU_CUTOFF] = (*serial_gpu_cutoff)();
      slvr.data[CRITICAL_TIME_SCALE_FACTOR] = (*critical_time_factor)();
      slvr.data[ENERGY_CHECK_EPSILON] = (*energy_check_epsilon)();
      slvr.data[VISCOUS_DAMPING_FACTOR] = (*viscous_damping_factor)();
      // load data
      if (amplitude_cat->GetType() == "Smooth Step")
      {
        int_slvr.data[AMPLITUDE_TYPE] =
            static_cast<int>(AmplitudeType::SmoothStep);
      }
      else if (amplitude_cat->GetType() == "Smooth Step and Hold")
      {
        int_slvr.data[AMPLITUDE_TYPE] =
            static_cast<int>(AmplitudeType::SmoothStepHold);
        const auto * hold_time = mt::GetCategoryModelTraitByType<mt::ScalarMT>(
            amplitude_cat, "hold time");
        if (hold_time == nullptr)
        {
          std::cerr << "\"hold time\" missing from \"Smooth Step and Hold\" "
                       "amplitude.\n";
          MPI_Abort(AMSI_COMM_WORLD, 1);
        }
        slvr.data[HOLD_TIME] = (*hold_time)();
      }
      else
      {
        std::cerr << amplitude_cat->GetType()
                  << " is not a valid amplitude type. attDefs and "
                     "MultiscaleTissue are out of sync."
                  << std::endl;
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      const auto * load_time = mt::GetCategoryModelTraitByType<mt::ScalarMT>(
          amplitude_cat, "load time");
      if (load_time == nullptr)
      {
        std::cerr << "\"hold time\" missing from amplitude definition.\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      slvr.data[LOAD_TIME] = (*load_time)();
      const auto * history_output = mt::GetPrimaryCategoryByType(
          microscale_convergence, "history output");
      const auto * field_output =
          mt::GetPrimaryCategoryByType(microscale_convergence, "field output");
      if (history_output)
      {
        const auto * num_iterations =
            mt::GetCategoryModelTraitByType<mt::IntMT>(history_output,
                                                       "number of iterations");
        if (num_iterations == nullptr)
        {
          std::cerr << "history output must have \"number of iterations\"\n";
          MPI_Abort(AMSI_COMM_WORLD, 1);
        }
        int_slvr.data[PRINT_HISTORY_FREQUENCY] = (*num_iterations)();
      }
      else
      {
        int_slvr.data[PRINT_HISTORY_FREQUENCY] = 0;
      }
      if (field_output)
      {
        if (field_output->GetType() == "by iterations")
        {
          const auto * num_iterations =
              mt::GetCategoryModelTraitByType<mt::IntMT>(
                  field_output, "number of iterations");
          if (num_iterations == nullptr)
          {
            std::cerr << "\"by iterations\" field output must have \"number of "
                         "iterations\".\n";
            MPI_Abort(AMSI_COMM_WORLD, 1);
          }
          int_slvr.data[PRINT_FIELD_BY_NUM_FRAMES] = 0;
          int_slvr.data[PRINT_FIELD_FREQUENCY] = (*num_iterations)();
        }
        else if (field_output->GetType() == "by frames")
        {
          int_slvr.data[PRINT_FIELD_BY_NUM_FRAMES] = 1;
          const auto * num_frames = mt::GetCategoryModelTraitByType<mt::IntMT>(
              field_output, "number of frames");
          if (num_frames == nullptr)
          {
            std::cerr << "\"by fames\" field output must have \"number of "
                         "frames\".\n";
            MPI_Abort(AMSI_COMM_WORLD, 1);
          }
          int_slvr.data[PRINT_FIELD_FREQUENCY] = (*num_frames)();
        }
        else
        {
          std::cerr << field_output->GetType()
                    << " is not a valid field output type. MultiscaleTissue "
                       "and attDefs are out of sync."
                    << std::endl;
          MPI_Abort(AMSI_COMM_WORLD, 1);
        }
      }
      else
      {
        int_slvr.data[PRINT_FIELD_FREQUENCY] = 0;
        int_slvr.data[PRINT_FIELD_BY_NUM_FRAMES] = 0;
      }
    }
    else if (microscale_convergence->GetType() ==
             "implicit nonlinear iteration")
    {
      int_slvr.data[MICRO_SOLVER_TYPE] = static_cast<int>(SolverType::Implicit);
      const auto * detect_oscillation = mt::GetPrimaryCategoryByType(
          microscale_convergence, "oscillation detection");
      if (detect_oscillation == nullptr)
      {
        std::cerr << "\"oscillation detection\" is required on implicit "
                     "microscale solver.\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      const auto * num_attempts = mt::GetCategoryModelTraitByType<mt::IntMT>(
          detect_oscillation, "number of attempts");
      const auto * cut_factor = mt::GetCategoryModelTraitByType<mt::IntMT>(
          detect_oscillation, "attempt cut factor");
      const auto * iteration_cap = mt::GetCategoryModelTraitByType<mt::IntMT>(
          detect_oscillation, "micro iteration cap");
      const auto * prev_iteration_factor =
          mt::GetCategoryModelTraitByType<mt::ScalarMT>(detect_oscillation,
                                                        "previous itr factor");
      if (detect_oscillation->GetType() == "iteration only")
      {
        int_slvr.data[DETECT_OSCILLATION_TYPE] =
            static_cast<int>(amsi::DetectOscillationType::IterationOnly);
        if (iteration_cap == nullptr)
        {
          std::cerr << "\"iteration cap\" required for iteration only "
                       "oscillation detection.\n";
          MPI_Abort(AMSI_COMM_WORLD, 1);
        }
      }
      else if (detect_oscillation->GetType() == "previous residual")
      {
        int_slvr.data[DETECT_OSCILLATION_TYPE] =
            static_cast<int>(amsi::DetectOscillationType::PrevNorm);
        if (prev_iteration_factor == nullptr)
        {
          std::cerr << "\"previous itr factor\" required for previous norm "
                       "oscillation detection.\n";
          MPI_Abort(AMSI_COMM_WORLD, 1);
        }
      }
      else if (detect_oscillation->GetType() == "combined")
      {
        int_slvr.data[DETECT_OSCILLATION_TYPE] =
            static_cast<int>(amsi::DetectOscillationType::IterationPrevNorm);
        if (iteration_cap == nullptr || prev_iteration_factor == nullptr)
        {
          std::cerr << "\"iteration cap\" and \"previous itr factor\" required "
                       "for iteration only oscillation detection.\n";
          MPI_Abort(AMSI_COMM_WORLD, 1);
        }
      }
      else
      {
        std::cerr << detect_oscillation->GetType()
                  << " is not a valid oscillation detection type. attDefs and "
                     "MultiscaleTissue are out of sync.\n";
        std::abort();
      }
      if (num_attempts == nullptr || cut_factor == nullptr)
      {
        std::cerr << "\"number of attempts\" and \"attempt cut factor\" are "
                     "required for oscillation detection.\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      slvr.data[PREV_ITER_FACTOR] =
          prev_iteration_factor ? (*prev_iteration_factor)() : 0;
      int_slvr.data[MAX_MICRO_CUT_ATTEMPT] = (*num_attempts)();
      int_slvr.data[MICRO_ATTEMPT_CUT_FACTOR] = (*cut_factor)();
      int_slvr.data[MAX_MICRO_ITERS] = iteration_cap ? (*iteration_cap)() : 0;
      const auto * convergence_tolerance =
          mt::GetCategoryModelTraitByType<mt::ScalarMT>(
              microscale_convergence, "micro convergence tolerance");
      if (convergence_tolerance == nullptr)
      {
        std::cerr << "\"micro convergence tolerance\" is required for implicit "
                     "analysis.\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      slvr.data[MICRO_CONVERGENCE_TOL] = (*convergence_tolerance)();
      // should not choose number less than 1E-6 for now due to LAS using single
      // precision
      assert(slvr.data[MICRO_SOLVER_TOL] <= slvr.data[MICRO_CONVERGENCE_TOL]);
    }
    else
    {
      std::cerr << microscale_convergence->GetType()
                << " is not a valid microscale convergence type. attDefs and "
                   " MultiscaleTissue are out of sync.\n";
      std::abort();
    }
  }
  void MultiscaleTissue::getInternalRVEData(apf::MeshEntity * rgn,
                                            micro_fo_header & hdr,
                                            micro_fo_params &,
                                            micro_fo_init_data & dat)
  {
    int idx = getRVEDirectoryIndex(rgn);
    hdr.data[RVE_DIR_TYPE] = idx;
    hdr.data[FIELD_ORDER] = apf::getShape(apf_primary_field)->getOrder();
    hdr.data[ELEMENT_TYPE] = apf_mesh->getType(rgn);
    hdr.data[GAUSS_ID] = -1;  // needs to be set for each IP
    apf::MeshElement * mlm = apf::createMeshElement(apf_mesh, rgn);
    apf::Element * ce = apf::createElement(apf_mesh->getCoordinateField(), mlm);
    apf::NewArray<apf::Vector3> Ni;
    apf::getVectorNodes(ce, Ni);
    int nds = apf::countNodes(ce);
    for (int nd = 0; nd < nds; ++nd)
      Ni[nd].toArray(&dat.init_data[nd * 3]);
    apf::destroyElement(ce);
    apf::destroyMeshElement(mlm);
  }
  int MultiscaleTissue::getRVEDirectoryIndex(apf::MeshEntity * ent)
  {
    apf::ModelEntity * gEnt = apf_mesh->toModel(ent);
    int ii = -1;
    const auto * model_traits = problem_definition.associated->Find(
        {apf_mesh->getModelType(gEnt), apf_mesh->getModelTag(gEnt)});
    const auto * material_model =
        model_traits->FindCategoryByType("material model");
    if (material_model == nullptr)
    {
      std::cerr << "\"material model\" must exist in \"problem definition\".\n";
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    const auto * multiscale_model =
        mt::GetCategoryByType(material_model, "multiscale model");
    if (multiscale_model == nullptr)
    {
      std::cerr
          << "\"multiscale model\" material must exist in multiscale region.\n";
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    MicroscaleType micro_tp = getMicroscaleType(multiscale_model);
    multiscale_model = &multiscale_model->GetCategories()[0];
    const auto * stochastic_field =
        mt::GetCategoryByType(multiscale_model, "stochastic field");
    std::string tp;
    if (micro_tp == MicroscaleType::FIBER_ONLY)
    {
      if (multiscale_model == nullptr)
      {
        std::cerr << "fiber only type should exist in the multiscale model";
        std::abort();
      }
      const auto * directory_mt = mt::GetCategoryModelTraitByType<mt::StringMT>(
          multiscale_model, "directory");
      const auto * prefix_mt = mt::GetCategoryModelTraitByType<mt::StringMT>(
          multiscale_model, "prefix");
      if (directory_mt == nullptr || prefix_mt == nullptr)
      {
        std::cerr
            << "\"directory\" and \"prefix\" categories are components of "
               "the multiscale model.\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
      std::string directory = (*directory_mt)();
      std::string prefix = (*prefix_mt)();
      if (stochastic_field != nullptr)
      {
        auto stochastic_field_data = read_stochastic_field(*stochastic_field);
        auto search_alignment =
            stochastic_field_map.find(stochastic_field_data.alignment_filename);
        if (search_alignment == stochastic_field_map.end())
        {
          std::cerr << "Alignment field "
                    << stochastic_field_data.alignment_filename
                    << " should be loaded in the loadRVELibraryInfo function"
                    << std::endl;
          MPI_Abort(AMSI_COMM_WORLD, 1);
        }
        auto search_orientation = stochastic_field_map.find(
            stochastic_field_data.orientation_filename);
        if (search_orientation == stochastic_field_map.end())
        {
          std::cerr << "Alignment field "
                    << stochastic_field_data.orientation_filename
                    << " should be loaded in the loadRVELibraryInfo function"
                    << std::endl;
          MPI_Abort(AMSI_COMM_WORLD, 1);
        }
        auto glbPnt = apf::getLinearCentroid(apf_mesh, ent);
        tp = directory + "/";
        tp += getNetworkSuffix(*(search_alignment->second),
                               *(search_orientation->second), prefix, glbPnt[0],
                               glbPnt[1], glbPnt[2],
                               stochastic_field_data.alignment_num_bins,
                               stochastic_field_data.orientation_num_bins);
      }
      else
      {
        tp = directory + "/" + prefix;
      }
      auto fnd = std::find(rve_dirs.begin(), rve_dirs.end(), tp);
      if (fnd != rve_dirs.end())
      {
        ii = std::distance(rve_dirs.begin(), fnd);
      }
      else
      {
        std::cerr << "network directory name " << tp
                  << " not found in library.\n";
        MPI_Abort(AMSI_COMM_WORLD, 1);
      }
    }
    return ii;
  }
}  // namespace mumfim
