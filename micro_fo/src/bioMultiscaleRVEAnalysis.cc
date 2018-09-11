#include "bioMultiscaleCoupling.h"
#include "bioMultiscaleRVEAnalysis.h"
#include <amsiDetectOscillation.h>
#include <amsiNonlinearAnalysis.h>
#include <apfFEA.h>  // amsi
#include <apfMDS.h>
//#include <apfMatrixUtil.h>    //amsi
#include <apfMeshIterator.h>  // amsi
#include <bioVerbosity.h>
#include <gmi.h>
#include <lasCSRCore.h>
#include <cassert>
#include "bioFiberNetworkIO.h"
#include "bioFiberRVEAnalysis.h"
#include "bioMultiscaleMicroFOParams.h"
//#include "bioRVEVolumeTerms.h"
namespace bio
{
  MultiscaleRVEAnalysis::~MultiscaleRVEAnalysis()
  {
    delete bfrs;
    for (auto v = vecs.begin(); v != vecs.end(); ++v)
    {
      delete (*v);
    }
    for (auto rve = ans.begin(); rve != ans.end(); ++rve)
    {
      destroyAnalysis(*rve);
      (*rve) = NULL;
    }
    assert(fns.size() == sprs.size() && fns.size() == dofs_cnt.size());
    // need to delete any fns that were not used as part of an analysis
    for (std::size_t i = 0; i < fns.size(); ++i)
    {
      for (int j = 0; j < rve_tp_cnt[i]; ++j)
      {
        delete fns[i][j];
        las::destroySparsity<las::CSR *>(sprs[i][j]);
        if (meshes[i][j])
        {
          meshes[i][j]->destroyNative();
          apf::destroyMesh(meshes[i][j]);
          meshes[i][j] = NULL;
        }
      }
      delete[] fns[i];
      delete[] sprs[i];
      delete[] dofs_cnt[i];
      delete[] meshes[i];
    }
    fns.clear();
    sprs.clear();
    dofs_cnt.clear();
    vecs.clear();
    meshes.clear();
  }
  MultiscaleRVEAnalysis::MultiscaleRVEAnalysis()
      : eff()
      , wgt()
      , tmg()
      , recv_ptrn()
      , send_ptrn()
      , rve_dd(NULL)
      , M2m_id()
      , m2M_id()
      //    , dat_tp()
      , rve_tp_cnt(0)
      , fns()
      , hdrs()
      , prms()
      , slvr_prms()
      , slvr_int_prms()
      , sprs()
      , bfrs(NULL)
      , vecs()
      , dofs_cnt()
      , macro_iter(0)
      , macro_step(0)
  {
    M2m_id = amsi::getRelationID(amsi::getMultiscaleManager(),
                                 amsi::getScaleManager(),
                                 "macro",
                                 "micro_fo");
    m2M_id = amsi::getRelationID(amsi::getMultiscaleManager(),
                                 amsi::getScaleManager(),
                                 "micro_fo",
                                 "macro");
  }
  void MultiscaleRVEAnalysis::init()
  {
    initLogging();
    initCoupling();
    initAnalysis();
  }
  void MultiscaleRVEAnalysis::initLogging()
  {
    eff = amsi::activateLog("micro_fo_efficiency");
    wgt = amsi::activateLog("micro_fo_weights");
    tmg = amsi::activateLog("micro_fo_timing");
    rve_tp_lg = amsi::activateLog("rve_type_log");
  }
  void MultiscaleRVEAnalysis::initCoupling()
  {
    amsi::ControlService * cs = amsi::ControlService::Instance();
    rve_dd = amsi::createDataDistribution(amsi::getLocal(), "macro_fo_data");
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
    char ** rve_tp_dirs = new char *[num_rve_tps];
    MPI_Request rqsts[num_rve_tps];
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
      cs->aRecvBroadcast(&rqsts[ii], M2m_id, &rve_tp_dirs[ii][0], cnt);
    }
    MPI_Status stss[num_rve_tps];
    MPI_Waitall(num_rve_tps, &rqsts[0], &stss[0]);
    MPI_Request hdr_rqst;
    rve_tp_cnt.resize(num_rve_tps);
    cs->aRecvBroadcast(&hdr_rqst, M2m_id, &rve_tp_cnt[0], num_rve_tps);
    MPI_Status hdr_sts;
    MPI_Waitall(1, &hdr_rqst, &hdr_sts);
    // Read in all the fiber network meshes and reactions
    for (int ii = 0; ii < num_rve_tps; ++ii)
    {
      fns.push_back(new FiberNetworkReactions *[rve_tp_cnt[ii]]);
      sprs.push_back(new las::Sparsity *[rve_tp_cnt[ii]]);
      dofs_cnt.push_back(new int[rve_tp_cnt[ii]]);
      meshes.push_back(new apf::Mesh2 *[rve_tp_cnt[ii]]);
    }
    int dof_max = -1;
    int nnz_max = -1;
    PCU_Switch_Comm(MPI_COMM_SELF);
    for (int ii = 0; ii < num_rve_tps; ii++)
    {
      for (int jj = 0; jj < rve_tp_cnt[ii]; ++jj)
      {
        std::stringstream fl;
        fl << rve_tp_dirs[ii] << jj + 1 << ".txt";
        FiberNetworkReactions * fn_rctns = new FiberNetworkReactions;
        meshes[ii][jj] = loadFromFile(fl.str());
        apf::Mesh2 * fn = meshes[ii][jj];
        fn_rctns->fileName = fl.str();
        fl << ".params";
        loadParamsFromFile(fn, fl.str(), std::back_inserter(fn_rctns->rctns));
        fns[ii][jj] = fn_rctns;
        apf::Field * u = apf::createLagrangeField(fn, "u", apf::VECTOR, 1);
        apf::zeroField(u);
        apf::Numbering * n = apf::createNumbering(u);
        int dofs = apf::NaiveOrder(n);
        sprs[ii][jj] = las::createCSR(n, dofs);
        int nnz = ((las::CSR*)sprs[ii][jj])->getNumNonzero();
        dofs_cnt[ii][jj] = dofs;
        apf::destroyNumbering(n);
        apf::destroyField(u);
        dof_max = dofs > dof_max ? dofs : dof_max;
        nnz_max = nnz > nnz_max ? nnz : nnz_max;
      }
    }
    assert(nnz_max > 0);
    assert(dof_max > 0);
    PCU_Switch_Comm(AMSI_COMM_SCALE);
    bfrs = new las::SparskitBuffers(dof_max, nnz_max);
  }
  void MultiscaleRVEAnalysis::updateCoupling()
  {
    std::vector<micro_fo_init_data> inis;
    std::vector<int> to_delete;
    amsi::ControlService * cs = amsi::ControlService::Instance();
    cs->RemoveData(recv_ptrn, to_delete);
    for (auto idx = to_delete.rbegin(); idx != to_delete.rend(); ++idx)
    {
      // FIXME need to clear out LinearStructs,
      // sparsity, dofs_cnt, etc here.
      // We also cannot just destroy the analysis w/o removing the corresponding
      // terms from the ans vector see vector erase
      // destroyAnalysis(ans[*idx]);
    }
    std::vector<int> to_add;
    std::vector<int> empty;
    size_t recv_init_ptrn = cs->AddData(recv_ptrn, empty, to_add);
    ans.resize(ans.size() + to_add.size());
    cs->Communicate(recv_init_ptrn, hdrs,
                    amsi::mpi_type<bio::micro_fo_header>());
    cs->Communicate(recv_init_ptrn, prms,
                    amsi::mpi_type<bio::micro_fo_params>());
    cs->Communicate(recv_init_ptrn, inis,
                    amsi::mpi_type<bio::micro_fo_init_data>());
    cs->Communicate(recv_init_ptrn, slvr_prms,
                    amsi::mpi_type<bio::micro_fo_solver>());
    cs->Communicate(recv_init_ptrn, slvr_int_prms,
                    amsi::mpi_type<bio::micro_fo_int_solver>());
    PCU_Switch_Comm(MPI_COMM_SELF);
    int ii = 0;
    for (auto rve = ans.begin(); rve != ans.end(); ++rve)
    {
      if (*rve == NULL)
      {
        micro_fo_header & hdr = hdrs[ii];
        micro_fo_params & prm = prms[ii];
        micro_fo_init_data & dat = inis[ii];
        micro_fo_solver & slvr_prm = slvr_prms[ii];
        micro_fo_int_solver & slvr_int_prm = slvr_int_prms[ii];
        int tp = hdr.data[RVE_TYPE];
        int rnd = rand() % rve_tp_cnt[tp];
        apf::Mesh * msh_cpy =
            apf::createMdsMesh(gmi_load(".null"), meshes[tp][rnd]);
        FiberNetwork * fn = new FiberNetwork(msh_cpy);
        fn->setFiberReactions(fns[tp][rnd]->rctns);
        vecs.push_back(
            createLinearStructs(dofs_cnt[tp][rnd],slvr_prm.data[MICRO_SOLVER_TOL] ,sprs[tp][rnd], bfrs));
        *rve = initFromMultiscale(
            fn, vecs[ii], hdr, prm, dat, slvr_prm, slvr_int_prm);
        fn->setRVEType(ii);
        BIO_V2(
            // print the list of fiber network names to file
            amsi::log(rve_tp_lg)
                << ii << " " << fns[tp][rnd]->fileName << std::endl;)
        ++ii;
      }
    }
    PCU_Switch_Comm(AMSI_COMM_SCALE);
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::ofstream rve_tp_lg_fs(
        amsi::fs->getResultsDir() + "/rve_tp." + std::to_string(rank) + ".log",
        std::ios::out | std::ios::app);
    amsi::flushToStream(rve_tp_lg, rve_tp_lg_fs);
    cs->CommPattern_UpdateInverted(recv_ptrn, send_ptrn);
    cs->CommPattern_Assemble(send_ptrn);
    cs->CommPattern_Reconcile(send_ptrn);
  }
  void MultiscaleRVEAnalysis::run()
  {
    amsi::ControlService * cs = amsi::ControlService::Instance();
    bool sim_complete = false;
    while (!sim_complete)
    {
      bool step_complete = false;
      while (!step_complete)
      {
        // migration
        if (macro_iter == 0) updateCoupling();
        std::vector<micro_fo_data> data;
        cs->Communicate(recv_ptrn, data, amsi::mpi_type<micro_fo_data>());
        std::vector<micro_fo_result> results(data.size());
        PCU_Switch_Comm(MPI_COMM_SELF);
        int rank = -1;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        int ii = 0;
        BIO_V1(double t0 = MPI_Wtime();)
        for (auto rve = ans.begin(); rve != ans.end(); ++rve)
        {
          // TODO this is somewhat inefficient and will not be needed when we
          // directly communicate the DeformationGradient type
          DeformationGradient dfmGrd;
          for (int jj = 0; jj < 9; ++jj)
          {
            dfmGrd[jj] = data[ii].data[jj];
          }
          bool result = (*rve)->run(dfmGrd);
          if (!result)
          {
            std::cerr << "Failed during  Macro Step " << macro_step
                      << " and macro iteration " << macro_iter << "."
                      << std::endl;
            // write the fiber network out if we fail, so we can test externally
            std::stringstream sout;
            sout << amsi::fs->getResultsDir() << "/"
                 << "rnk_" << rank << "_fn_" << (*rve)->getFn()->getRVEType()
                 << "_step_" << macro_step << "_iter_" << macro_iter;
            apf::writeVtkFiles(
                sout.str().c_str(), (*rve)->getFn()->getNetworkMesh(), 1);
            // if the rve didn't work, crash the analysis
            std::abort();  // should I use MPI_Abort() here?
          }
          // we've converged and have not reset the state of the vectors,
          // matrices, and buffers the inversion of the tangent stiffness matrix
          // should be available in the buffers?
          recoverMultiscaleResults(*rve, &results[ii]);
          ii++;
#ifdef WRITE_MICRO_PER_ITER
          std::stringstream sout;
          sout << "rnk_" << rank << "_fn_" << (*rve)->getFn()->getRVEType()
               << "_step_" << macro_step << "_iter_" << macro_iter;
          apf::writeVtkFiles(
              sout.str().c_str(), (*rve)->getFn()->getNetworkMesh(), 1);
#endif
        }
        BIO_V1(double t1 = MPI_Wtime();)
        BIO_V1(std::cout << "Computed " << ans.size() << " RVEs in " << t0 - t1
                         << " seconds." << std::endl;)
        PCU_Switch_Comm(AMSI_COMM_SCALE);
        cs->Communicate(send_ptrn, results, amsi::mpi_type<micro_fo_result>());
        macro_iter++;
        cs->scaleBroadcast(M2m_id, &step_complete);
      }
      // get the size of the step results vector
      std::vector<micro_fo_step_result> step_results(hdrs.size());
      // recover step results and set the step results vector
      int i = 0;
      PCU_Switch_Comm(MPI_COMM_SELF);
      for (auto rve = ans.begin(); rve != ans.end(); ++rve)
      {
        micro_fo_header & hdr = hdrs[i];
        micro_fo_params & prm = prms[i];
        recoverMultiscaleStepResults(*rve, hdr, prm, &step_results[i]);
        ++i;
      }
      PCU_Switch_Comm(AMSI_COMM_SCALE);
      // communicate the step results back to the macro scale
      amsi::ControlService * cs = amsi::ControlService::Instance();
      cs->Communicate(
          send_ptrn, step_results, amsi::mpi_type<micro_fo_step_result>());
#ifdef WRITE_MICRO_PER_STEP
      for (auto rve = ans.begin(); rve != ans.end(); ++rve)
      {
        std::stringstream sout;
        int rnk = -1;
        MPI_Comm_rank(MPI_COMM_WORLD, &rnk);
        int ii = 0;
        sout << "rnk_" << rnk << "_fn_" << (*rve)->getFn()->getRVEType()
             << "_step_" << macro_step << "_iter_" << macro_iter;
        apf::writeVtkFiles(sout.str().c_str(), (*rve)->fn->getNetworkMesh(), 1);
        sout.str("");
        ii++;
      }
#endif
      macro_iter = 0;
      macro_step++;
      cs->scaleBroadcast(M2m_id, &sim_complete);
    }
  }
}  // namespace bio
