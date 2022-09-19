#include "MultiscaleRVEAnalysis.h"
#include <amsiCommunicationManager.h>
#include <amsiControlService.h>
#include <amsiDetectOscillation.h>
#include <amsiNonlinearAnalysis.h>
#include <amsiTaskManager.h>
#include <apfMDS.h>
#include <apfMeshIterator.h>  // amsi
#include <gmi.h>
#include <mumfim/microscale/Verbosity.h>
#include <cassert>
#include <sstream>
//#include "BatchedNeohookeanAnalysis.h"
#include "FiberNetworkLibrary.h"
#include "FiberRVEAnalysis.h"
#include "MultiscaleCoupling.h"
#include "MultiscaleMicroFOParams.h"
#include "NeoHookeanRVEAnalysis.h"
#include "mumfim/exceptions.h"
namespace mumfim
{
  MultiscaleRVEAnalysis::~MultiscaleRVEAnalysis() = default;
  MultiscaleRVEAnalysis::MultiscaleRVEAnalysis(
      const amsi::Multiscale & amsi_multiscale)
      : recv_ptrn()
      , send_ptrn()
      , rve_dd(nullptr)
      , M2m_id()
      , m2M_id()
      , multiscale_(amsi_multiscale)
      , rve_tp_cnt(0)
      , hdrs()
      , prms()
      , slvr_prms()
      , slvr_int_prms()
      , rve_tp_dirs()
      , network_library()
      , macro_iter(0)
      , macro_step(0)
      , initial_update(true)
  {
    if (multiscale_.getMultiscaleManager() == nullptr)
    {
      throw mumfim_error{
          "Multiscale RVE Analysis cannot run without a multiscale manager"};
    }
    if (multiscale_.getScaleManager() == nullptr)
    {
      throw mumfim_error{
          "Multiscale RVE Analysis cannot run without a scale manager"};
    }
    M2m_id =
        amsi::getRelationID(multiscale_.getMultiscaleManager(),
                            multiscale_.getScaleManager(), "macro", "micro_fo");
    m2M_id =
        amsi::getRelationID(multiscale_.getMultiscaleManager(),
                            multiscale_.getScaleManager(), "micro_fo", "macro");
  }
  void MultiscaleRVEAnalysis::init()
  {
    initCoupling();
    initAnalysis();
  }
  void MultiscaleRVEAnalysis::initCoupling()
  {
    auto * cs = multiscale_.getControlService();
    rve_dd = amsi::createDataDistribution(
        multiscale_.getScaleManager()->getLocalTask(), "macro_fo_data");
    recv_ptrn = cs->RecvCommPattern("micro_fo_data", "macro", "macro_fo_data",
                                    "micro_fo");
    cs->CommPattern_Reconcile(recv_ptrn);
    send_ptrn = cs->CommPattern_UseInverted(recv_ptrn, "macro_fo_data",
                                            "micro_fo", "macro");
    cs->CommPattern_Assemble(send_ptrn);
    cs->CommPattern_Reconcile(send_ptrn);
  }
  void MultiscaleRVEAnalysis::initAnalysis()
  {
    int num_rve_tps = 0;
    amsi::ControlService * cs = amsi::ControlService::Instance();
    cs->scaleBroadcast(M2m_id, &num_rve_tps);
    rve_tp_dirs.resize(num_rve_tps);
    std::vector<MPI_Request> rqsts;
    // the order of receipt might be non-deterministic. need to handle that
    for (int ii = 0; ii < num_rve_tps; ++ii)
    {
      int cnt = 0;
      while ((cnt = cs->aRecvBroadcastSize<char>(M2m_id)) == 0)
      {
      }
      rve_tp_dirs[ii] = new char[cnt];
      // don't have to block to wait since we know the message was available for
      // size info
      cs->aRecvBroadcast(std::back_inserter(rqsts), M2m_id, &rve_tp_dirs[ii][0],
                         cnt);
    }
    rve_tp_cnt.resize(num_rve_tps);
    cs->aRecvBroadcast(std::back_inserter(rqsts), M2m_id, &rve_tp_cnt[0],
                       num_rve_tps);
    // assert(rqsts.size() == num_rve_tps);
    //  note rather than waiting for all of the communication, we can interleave
    //  creating some of the new structures with some of the communication
    MPI_Waitall(rqsts.size(), &rqsts[0], MPI_STATUSES_IGNORE);
  }
  void MultiscaleRVEAnalysis::updateCoupling()
  {
    int rank;
    std::vector<micro_fo_init_data> inis;
    std::vector<int> to_delete;
    amsi::ControlService * cs = amsi::ControlService::Instance();
    cs->RemoveData(recv_ptrn, to_delete);
    if (to_delete.size() > 0)
    {
      std::cerr << "Deleting RVEs is not currently supported!" << std::endl;
      MPI_Abort(AMSI_COMM_WORLD, EXIT_FAILURE);
    }
    std::vector<int> to_add;
    std::vector<int> empty;
    size_t recv_init_ptrn = cs->AddData(recv_ptrn, empty, to_add);
    if (!initial_update && to_add.size() > 0)
    {
      std::cerr << "Adding RVEs on the fly is not currently supported!"
                << std::endl;
      MPI_Abort(AMSI_COMM_WORLD, EXIT_FAILURE);
    }
    initial_update = false;
    cs->Communicate(recv_init_ptrn, hdrs,
                    amsi::mpi_type<mumfim::micro_fo_header>());
    cs->Communicate(recv_init_ptrn, prms,
                    amsi::mpi_type<mumfim::micro_fo_params>());
    cs->Communicate(recv_init_ptrn, inis,
                    amsi::mpi_type<mumfim::micro_fo_init_data>());
    cs->Communicate(recv_init_ptrn, slvr_prms,
                    amsi::mpi_type<mumfim::micro_fo_solver>());
    cs->Communicate(recv_init_ptrn, slvr_int_prms,
                    amsi::mpi_type<mumfim::micro_fo_int_solver>());
    // store the current communicator which should be AMSI_COMM_SCALE
    MPI_Comm comm = PCU_Get_Comm();
    // All of the APF mesh operations need to have the PCU communicator
    // set to self, otherwise it will try to do things in parallel
    PCU_Switch_Comm(MPI_COMM_SELF);
    std::vector<std::shared_ptr<const FiberNetwork>> fiber_networks;
    std::vector<std::shared_ptr<const MicroSolutionStrategy>>
        solution_strategies;
    for (std::size_t i = 0; i < to_add.size(); ++i)
    {
      micro_fo_header & hdr = hdrs[i];
      micro_fo_params & prm = prms[i];
      micro_fo_solver & slvr_prm = slvr_prms[i];
      micro_fo_int_solver & slvr_int_prm = slvr_int_prms[i];
      MicroscaleType micro_tp = static_cast<MicroscaleType>(hdr.data[RVE_TYPE]);
      if (micro_tp == MicroscaleType::FIBER_ONLY)
      {
        int tp = hdr.data[RVE_DIR_TYPE];
        int rnd = rand() % rve_tp_cnt[tp];
        std::string fiber_network_file =
            rve_tp_dirs[tp] + std::to_string(rnd + 1) + ".txt";
        auto fiber_network = network_library.load(
            fiber_network_file, fiber_network_file + ".params", tp, rnd);
        fiber_networks.push_back(fiber_network);
        solution_strategies.push_back(
            serializeSolutionStrategy(slvr_prm, slvr_int_prm));
        // FIXME This only needs to be done once per Fiber network, and
        // FIBER_RADIUS should come from network properties, not the multiscale
        // anlysis...FIBER_RADIUS needs to be removed from the communicated data
        double pi = 4 * atan(1);
        double fbr_area = pi * prm.data[FIBER_RADIUS] * prm.data[FIBER_RADIUS];
        double fbr_vol_frc = prm.data[VOLUME_FRACTION];
        double scale_factor = calcScaleConversion(
            fiber_network->getNetworkMesh(), fbr_area, fbr_vol_frc);
        fiber_network->setScaleConversion(scale_factor);
      }
      else if (micro_tp == MicroscaleType::ISOTROPIC_NEOHOOKEAN)
      {
        std::cerr
            << "Currently Neohookean Analysis is not supported in batched mode"
            << std::endl;
        //*rve = initNeoHookeanRVEAnalysisFromMultiscale(prm);
        MPI_Abort(AMSI_COMM_WORLD, EXIT_FAILURE);
      }
      else
      {
        std::cerr << "The Microscale/RVE type is not valid" << std::endl;
        MPI_Abort(AMSI_COMM_WORLD, EXIT_FAILURE);
      }
    }
    if (to_add.size() > 0)
    {
       batched_analysis = BatchedAnalysisType{
           new BatchedFiberRVEAnalysisExplicit<Scalar, LocalOrdinal,
                                               Kokkos::DefaultExecutionSpace>{
               std::move(fiber_networks), std::move(solution_strategies)}};
      //batched_analysis = BatchedAnalysisType{
      //    new BatchedNeohookeanAnalysis<Scalar, LocalOrdinal,
      //                                  Kokkos::DefaultExecutionSpace>(
      //        static_cast<int>(fiber_networks.size()),10000,0.3)};
    }
      // reset the PCU communictor back to its original state so that
      // our scale wide communications can proceed
      PCU_Switch_Comm(comm);
      cs->CommPattern_UpdateInverted(recv_ptrn, send_ptrn);
      cs->CommPattern_Assemble(send_ptrn);
      cs->CommPattern_Reconcile(send_ptrn);
    }
    void MultiscaleRVEAnalysis::run()
    {
      amsi::ControlService * cs = amsi::ControlService::Instance();
      if (macro_step == 0 && macro_iter == 0)
      {
        updateCoupling();
      }
      else
      {
        std::cerr << "Very surprised that I am here!\n";
        std::abort();
      }
      // send the initial step result data to the macroscale to output
      // get the size of the step results vector
      std::vector<micro_fo_step_result> step_results(hdrs.size());
      Kokkos::DualView<Scalar * [3][3]> deformation_gradient(
          "deformation gradient", hdrs.size());
      Kokkos::DualView<Scalar * [6]> stress("stress", hdrs.size());
      Kokkos::DualView<Scalar * [6][6]> material_stiffness("material stiffness",
                                                           hdrs.size());
      Kokkos::DualView<Scalar * [3][3]> orientation_tensor("orientation tensor",
                                                           hdrs.size());
      Kokkos::DualView<Scalar * [3]> orientation_tensor_normal(
          "orientation tensor normal", hdrs.size());
      // communicate the step results back to the macro scale
      if (macro_step == 0 && macro_iter == 0)
      {
        recoverMultiscaleStepResults(
            orientation_tensor, orientation_tensor_normal,
            batched_analysis.get(), hdrs, prms, step_results);
        cs->Communicate(send_ptrn, step_results,
                        amsi::mpi_type<micro_fo_step_result>());
      }
      bool sim_complete{false};
      while (!sim_complete)
      {
        bool step_complete{false};
        int step_accepted{0};
        while (!step_complete)
        {
          // migration
          // if (macro_iter == 0) updateCoupling();
          // send the initial microscale rve states back to the macroscale
          std::vector<micro_fo_data> data;
          cs->Communicate(recv_ptrn, data, amsi::mpi_type<micro_fo_data>());
          std::vector<micro_fo_result> results(data.size());
          PCU_Switch_Comm(MPI_COMM_SELF);
          int rank = -1;
          MPI_Comm_rank(AMSI_COMM_WORLD, &rank);
          // int ii = 0;
          MUMFIM_V1(double t0 = MPI_Wtime();)
          auto deformation_gradient_h = deformation_gradient.h_view;
          auto stress_h = stress.h_view;
          auto material_stiffness_h = material_stiffness.h_view;
          deformation_gradient.modify<Kokkos::HostSpace>();
          // fill the deformation gradient data
          for (std::size_t i = 0; i < results.size(); ++i)
          {
            for (int ei = 0; ei < 3; ++ei)
            {
              for (int ej = 0; ej < 3; ++ej)
              {
                deformation_gradient_h(i, ei, ej) = data[i].data[ei * 3 + ej];
              }
            }
          }
          batched_analysis->run(deformation_gradient, stress);
          batched_analysis->computeMaterialStiffness(material_stiffness);
          stress.sync<Kokkos::HostSpace>();
          material_stiffness.sync<Kokkos::HostSpace>();
          // fill the results data
          for (std::size_t i = 0; i < results.size(); ++i)
          {
            double * sigma = &(results[i].data[0]);
            double * avgVolStress = &(results[i].data[6]);
            double * matStiffMatrix = &(results[i].data[9]);
            for (int j = 0; j < 6; ++j)
            {
              sigma[j] = stress_h(i, j);
            }
            // avg vol stress which we don't currently use
            for (int j = 0; j < 3; ++j)
            {
              avgVolStress[j] = 0;
            }
            for (int j = 0; j < 36; ++j)
            {
              int row = (j) / 6;
              int col = (j) % 6;
              matStiffMatrix[j] = material_stiffness_h(i, row, col);
            }
          }
          MUMFIM_V1(double t1 = MPI_Wtime();)
          MUMFIM_V1(std::cout << "Computed " << results.size() << " RVEs in "
                              << t1 - t0 << " seconds. On rank " << rank << "."
                              << std::endl;)
          // PCU_Switch_Comm(AMSI_COMM_SCALE);
          cs->Communicate(send_ptrn, results,
                          amsi::mpi_type<micro_fo_result>());
          macro_iter++;

          cs->scaleBroadcast(M2m_id, &step_accepted);
          bool step_complete = (step_accepted > 0);
          if(step_accepted) {
            batched_analysis->accept();
          }
          //batched_analysis->accept();
        }
        //batched_analysis->accept();
        // get the size of the step results vector
        std::vector<micro_fo_step_result> step_results(hdrs.size());
        // recover step results and set the step results vector
        recoverMultiscaleStepResults(
            orientation_tensor, orientation_tensor_normal,
            batched_analysis.get(), hdrs, prms, step_results);
        // communicate the step results back to the macro scale
        cs->Communicate(send_ptrn, step_results,
                        amsi::mpi_type<micro_fo_step_result>());
        macro_iter = 0;
        macro_step++;
        cs->scaleBroadcast(M2m_id, &sim_complete);
      }
    }
    std::unique_ptr<MicroSolutionStrategy> serializeSolutionStrategy(
        micro_fo_solver & slvr, micro_fo_int_solver & slvr_int)
    {
      auto solver_type =
          static_cast<SolverType>(slvr_int.data[MICRO_SOLVER_TYPE]);
      auto solution_strategy = std::unique_ptr<MicroSolutionStrategy>{nullptr};
      if (solver_type == SolverType::Explicit)
      {
        solution_strategy.reset(new MicroSolutionStrategyExplicit);
        auto sse = static_cast<MicroSolutionStrategyExplicit *>(
            solution_strategy.get());
        sse->total_time = slvr.data[LOAD_TIME] + slvr.data[HOLD_TIME];
        sse->load_time = slvr.data[LOAD_TIME];
        sse->crit_time_scale_factor = slvr.data[CRITICAL_TIME_SCALE_FACTOR];
        sse->visc_damp_coeff = slvr.data[VISCOUS_DAMPING_FACTOR];
        sse->energy_check_eps = slvr.data[ENERGY_CHECK_EPSILON];
        sse->ampType =
            static_cast<AmplitudeType>(slvr_int.data[AMPLITUDE_TYPE]);
        sse->print_history_frequency = slvr_int.data[PRINT_HISTORY_FREQUENCY];
        sse->print_field_frequency = slvr_int.data[PRINT_FIELD_FREQUENCY];
        sse->print_field_by_num_frames =
            slvr_int.data[PRINT_FIELD_BY_NUM_FRAMES];
        sse->serial_gpu_cutoff = slvr_int.data[SERIAL_GPU_CUTOFF];
      }
      else if (solver_type == SolverType::Implicit)
      {
        solution_strategy.reset(new MicroSolutionStrategy);
      }
      else
      {
        std::cerr << "Attempting to use a type of solver which has not been "
                     "implemented yet.\n";
        std::cerr << "This is most likely configuration error, but it could "
                     "also happen if there";
        std::cerr << " is an MPI communication error" << std::endl;
        std::abort();
      }
      // set the combined solver parameters
      if (solution_strategy != nullptr)
      {
        solution_strategy->cnvgTolerance = slvr.data[MICRO_CONVERGENCE_TOL];
        solution_strategy->slvrTolerance = slvr.data[MICRO_SOLVER_TOL];
        solution_strategy->oscPrms.maxIterations =
            slvr_int.data[MAX_MICRO_ITERS];
        solution_strategy->oscPrms.maxMicroCutAttempts =
            slvr_int.data[MAX_MICRO_CUT_ATTEMPT];
        solution_strategy->oscPrms.microAttemptCutFactor =
            slvr_int.data[MICRO_ATTEMPT_CUT_FACTOR];
        solution_strategy->oscPrms.oscType =
            static_cast<amsi::DetectOscillationType>(
                slvr_int.data[DETECT_OSCILLATION_TYPE]);
        solution_strategy->oscPrms.prevNormFactor = slvr.data[PREV_ITER_FACTOR];
      }
      return solution_strategy;
    }
  }  // namespace mumfim
