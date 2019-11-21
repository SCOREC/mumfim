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
#include <lasCorePETSc.h>
#include <cassert>
#include "bioFiberNetworkIO.h"
#include "bioFiberRVEAnalysis.h"
#include "bioMultiscaleMicroFOParams.h"
#include <fstream>
#include <sstream>
#include "bioNeoHookeanRVEAnalysis.h"
namespace bio
{
  void MultiscaleRVEAnalysis::loadNetworkIntoLibrary(const std::string & network_name,
      size_t net_type, size_t net_id, int & num_dofs, int & nnz)
  {
      MPI_Comm comm = PCU_Get_Comm();
      PCU_Switch_Comm(MPI_COMM_SELF);
      FiberNetworkReactions * fn_rctns = new FiberNetworkReactions;
      meshes[net_type][net_id] = loadFromFile(network_name);
      apf::Mesh2 * fn = meshes[net_type][net_id];
      fn_rctns->fileName = network_name;
      loadParamsFromFile(fn, network_name+".params", std::back_inserter(fn_rctns->rctns));
      fns[net_type][net_id] = fn_rctns;
      apf::Field * u = apf::createLagrangeField(fn, "u", apf::VECTOR, 1);
      apf::zeroField(u);
      apf::Numbering * n = apf::createNumbering(u);
      num_dofs = apf::NaiveOrder(n);
      dofs_cnt[net_type][net_id] = num_dofs;
      //apf::NaiveOrder(n);
#if defined MICRO_USING_SPARSKIT
      sprs[net_type][net_id] = las::createCSR(n, num_dofs);
      nnz = ((las::CSR*)sprs[net_type][net_id])->getNumNonzero();
#elif defined MICRO_USING_PETSC
      sprs[net_type][net_id] = las::createPetscSparsity(n, num_dofs, PCU_Get_Comm());
#endif
      apf::destroyNumbering(n);
      apf::destroyField(u);
      PCU_Switch_Comm(comm);
  }
  MultiscaleRVEAnalysis::~MultiscaleRVEAnalysis()
  {
    delete bfrs;
    for (auto v = vecs.begin(); v != vecs.end(); ++v)
    {
      delete (*v);
    }
    for (auto rve = ans.begin(); rve != ans.end(); ++rve)
    {
      delete *rve;
      (*rve) = NULL;
    }
    assert(fns.size() == sprs.size() && fns.size() == dofs_cnt.size());
    // need to delete any fns that were not used as part of an analysis
    for (std::size_t i = 0; i < fns.size(); ++i)
    {
      for (int j = 0; j < rve_tp_cnt[i]; ++j)
      {
        delete fns[i][j];
        las::destroySparsity<las::MICRO_BACKEND>(sprs[i][j]);
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
      delete[] rve_tp_dirs[i];
      delete[] meshes[i];
    }
    fns.clear();
    sprs.clear();
    dofs_cnt.clear();
    rve_tp_dirs.clear();
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
      , rve_tp_dirs()
      , macro_iter(0)
      , macro_step(0)
      , nnz_max(-1)
  {
    M2m_id = amsi::getRelationID(amsi::getMultiscaleManager(),
                                 amsi::getScaleManager(),
                                 "macro",
                                 "micro_fo");
    m2M_id = amsi::getRelationID(amsi::getMultiscaleManager(),
                                 amsi::getScaleManager(),
                                 "micro_fo",
                                 "macro");
    //PCU_Switch_Comm(MPI_COMM_SELF);
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
      cs->aRecvBroadcast(std::back_inserter(rqsts), M2m_id, &rve_tp_dirs[ii][0], cnt);
    }
    rve_tp_cnt.resize(num_rve_tps);
    cs->aRecvBroadcast(std::back_inserter(rqsts), M2m_id, &rve_tp_cnt[0], num_rve_tps);
    //assert(rqsts.size() == num_rve_tps);
    // note rather than waiting for all of the communication, we can interleave
    // creating some of the new structures with some of the communication
    MPI_Waitall(rqsts.size(), &rqsts[0], MPI_STATUSES_IGNORE);
    // Read in all the fiber network meshes and reactions
    for (int ii = 0; ii < num_rve_tps; ++ii)
    {
      fns.push_back(new FiberNetworkReactions *[rve_tp_cnt[ii]]);
      sprs.push_back(new las::Sparsity *[rve_tp_cnt[ii]]);
      dofs_cnt.push_back(new int[rve_tp_cnt[ii]]);
      meshes.push_back(new apf::Mesh2 *[rve_tp_cnt[ii]]);
    }
    for (int ii = 0; ii < num_rve_tps; ii++)
    {
      for (int jj = 0; jj < rve_tp_cnt[ii]; ++jj)
      {
        // initialize the library to all 0/NULL
        fns[ii][jj] = NULL;
        sprs[ii][jj] = NULL;
        dofs_cnt[ii][jj] = 0;
        meshes[ii][jj] = NULL;
      }
    }
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
#ifdef MICRO_USING_SPARSKIT
    // resize the sparskit buffers
    if(rve_tp_cnt.size() >0)
    {
      // make minimum size buffer so that we can get a pointer.
      // this will be resized after all of the relevant meshes are loaded
      bfrs = new las::SparskitBuffers(1);
    }
#endif
    for (auto rve = ans.begin(); rve != ans.end(); ++rve)
    {
      if (*rve == NULL)
      {
        micro_fo_header & hdr = hdrs[ii];
        micro_fo_params & prm = prms[ii];
        micro_fo_init_data & dat = inis[ii];
        micro_fo_solver & slvr_prm = slvr_prms[ii];
        micro_fo_int_solver & slvr_int_prm = slvr_int_prms[ii];
        MicroscaleType micro_tp = static_cast<MicroscaleType>(hdr.data[RVE_TYPE]);
        if(micro_tp == MicroscaleType::FIBER_ONLY)
        {
          int tp = hdr.data[RVE_DIR_TYPE];
          //int rnd = 0;//rand() % rve_tp_cnt[tp];
          int rnd = rand() % rve_tp_cnt[tp];
          // load the mesh into the library if it hasn't already been loaded
          if(meshes[tp][rnd] == NULL)
          {
            std::stringstream fl;
            fl << rve_tp_dirs[tp] << rnd + 1 << ".txt";
            int num_dofs, nnz = 0;
            loadNetworkIntoLibrary(fl.str(), tp, rnd, num_dofs, nnz);
            // if the network is already loaded into the library, then comparing
            // its nnz against nnz_max cannot increase nnz_max, so only setting nnz_max here
            // should be safe
            nnz_max = (nnz>nnz_max) ? nnz : nnz_max;
          }
          apf::Mesh * msh_cpy =
              apf::createMdsMesh(gmi_load(".null"), meshes[tp][rnd]);
          // Fiber network takes ownership of the mesh copy
          FiberNetwork * fn = new FiberNetwork(msh_cpy);
          fn->setFiberReactions(fns[tp][rnd]->rctns);
          vecs.push_back(
              createLinearStructs(dofs_cnt[tp][rnd],slvr_prm.data[MICRO_SOLVER_TOL] ,sprs[tp][rnd], bfrs));
          *rve = initFiberRVEAnalysisFromMultiscale(
              fn, vecs[ii], hdr, prm, dat, slvr_prm, slvr_int_prm);
          fn->setRVEType(ii);
          BIO_V2(
              // print the list of fiber network names to file
              amsi::log(rve_tp_lg)
                  << ii << " " << fns[tp][rnd]->fileName << std::endl;)
          ++ii;
        }
        else if (micro_tp == MicroscaleType::ISOTROPIC_NEOHOOKEAN)
        {
          *rve = initNeoHookeanRVEAnalysisFromMultiscale(prm);
          assert(*rve);
        }
        else
        {
          std::cerr<<"The Microscale/RVE type is not valid"<<std::endl;
          MPI_Abort(AMSI_COMM_WORLD, 1);
        }
      }
    }
#ifdef MICRO_USING_SPARSKIT
    // resize the sparskit buffers
    // the max nnz will be negative 1 on any subsequent load steps
    if(rve_tp_cnt.size() > 0)
    {
      assert(nnz_max > 0);
      // magic number that doesn't seem to require much rescaling
      // if this buffer size doesn't work we can try nnz_max*10*log(nnz_max)
      bfrs->resizeMatrixBuffer(nnz_max*100);
      //bfrs = new las::SparskitBuffers(dof_max, nnz_max*100);
    }
#endif
    // resize sparskit buffers to max size here.
    PCU_Switch_Comm(AMSI_COMM_SCALE);
    int rank = -1;
    MPI_Comm_rank(AMSI_COMM_WORLD, &rank);
    std::stringstream ss;
    ss << amsi::fs->getResultsDir() << "/rve_tp." << rank << ".log";
    std::ofstream rve_tp_lg_fs(ss.str().c_str(), (std::ios::out | std::ios::app));
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
        MPI_Comm_rank(AMSI_COMM_WORLD, &rank);
        int ii = 0;
        BIO_V1(double t0 = MPI_Wtime();)
        for (auto rve = ans.begin(); rve != ans.end(); ++rve)
        {
          double * sigma = &(results[ii].data[0]);
          double * avgVolStress = &(results[ii].data[6]);
          double * matStiffMatrix = &(results[ii].data[9]);
          // TODO this is somewhat inefficient and will not be needed when we
          // directly communicate the DeformationGradient type
          DeformationGradient dfmGrd;
          for (int jj = 0; jj < 9; ++jj)
          {
            dfmGrd[jj] = data[ii].data[jj];
          }
          bool result = (*rve)->run(dfmGrd, sigma, true);
          // if the rve didn't work, crash the analysis
          if (!result)
          {
            std::cerr << "Failed during  Macro Step " << macro_step
                      << " and macro iteration " << macro_iter << ".\n"
                      << "Applied deformation gradient was: "
                      << "F=";
                      for(int i=0; i<9;++i) {
                        std::cerr<<dfmGrd[i]<<" ";
                      }
                      std::cerr<<std::endl;
            // use MPI abort, so that we don't wait for networks to complete if one
            // has failed.
            MPI_Abort(AMSI_COMM_WORLD, 1);
          }
          //recoverMultiscaleResults(*rve, &results[ii]);
          (*rve)->computeAvgVolStress(avgVolStress);
          (*rve)->computeMaterialStiffness(matStiffMatrix);
          FiberRVEAnalysis* FRveAns = dynamic_cast<FiberRVEAnalysis *>(*rve);
          if(FRveAns)
            convertStressQuantities(FRveAns, sigma, matStiffMatrix);
          ii++;
#ifdef WRITE_MICRO_PER_ITER
          FiberRVEAnalysis * FNRve = dynamic_cast<FiberRVEAnalysis *>(*rve);
          if(FNRve)
          {
            std::stringstream sout;
            sout << "rnk_" << rank << "_fn_" << (*FNRve)->getFn()->getFNRveType()
                 << "_step_" << macro_step << "_iter_" << macro_iter;
            apf::writeVtkFiles(
                sout.str().c_str(), (*FNRve)->getFn()->getNetworkMesh(), 1);
          }
          else
          {
            std::cerr<<"Not Writing microscale to output per iteration since this RVE type has no mesh representation"<<std::endl;
          }
#endif
        }
        BIO_V1(double t1 = MPI_Wtime();)
        BIO_V1(std::cout << "Computed " << ans.size() << " RVEs in " << t1 - t0
                         << " seconds. On rank "<<rank<<"."<< std::endl;)
        //PCU_Switch_Comm(AMSI_COMM_SCALE);
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
        FiberRVEAnalysis * FNRve2 = dynamic_cast<FiberRVEAnalysis *>(*rve);
        if(FNRve2)
        {
          std::stringstream sout;
          int rnk = -1;
          MPI_Comm_rank(AMSI_COMM_WORLD, &rnk);
          int ii = 0;
          sout << "rnk_" << rnk << "_fn_" << (*FNRve2)->getFn()->getRVEType()
               << "_step_" << macro_step << "_iter_" << macro_iter;
          apf::writeVtkFiles(sout.str().c_str(), (*FNRve2)->fn->getNetworkMesh(), 1);
          sout.str("");
          ii++;
        }
        else
        {
          std::cerr<<"Not Writing microscale to output per iteration since this RVE type has no mesh representation"<<std::endl;
        }
      }
#endif
      macro_iter = 0;
      macro_step++;
      cs->scaleBroadcast(M2m_id, &sim_complete);
    }
  }
}  // namespace bio
