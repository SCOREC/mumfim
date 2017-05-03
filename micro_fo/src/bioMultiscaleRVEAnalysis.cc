#include "bioMultiscaleRVEAnalysis.h"
namespace bio
{
  void applyMultiscaleCoupling(FiberRVEAnalysis * ans, micro_fo_data * data)
  {

  }
  void recoverMultiscaleResults(FiberRVEAnalysis * ans, micro_fo_results * data)
  {

  }
  MultiscaleRVEAnalysis::MultiscaleRVEAnaysis()
    : eff()
    , wgt()
    , tmg()
    , recv_ptrn()
    , send_ptrn()
    , rve_dd(NULL)
    , M2m_id()
    , m2M_id()
    , hdr_tp()
    , prm_tp()
    , ini_tp()
    , dat_tp()
    , rst_tp()
    , fns()
    , sprs()
    , bfrs(NULL)
    , macro_iter(0)
    , macro_step(0)
  { }
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
  }
  void MultiscaleRVEAnalysis::initCoupling()
  {
    amsi::ControlService * cs = amsi::ControlService::Instance();
    rve_dd = amsi::createDataDistribution(amsi::getLocal(),"macro_fo_data");
    recv_ptrn = cs->RecvCommPattern("micro_fo_data","macro",
                                    "macro_fo_data","micro_fo");
    cs->CommPattern_Reconcile(recv_ptrn);
    send_ptrn = cs->CommPattern_UseInverted(recv_ptrn,"macro_fo_data",
                                            "micro_fo","macro");
    mdt.MultiscaleDataTypesMPICommit();
    hdr_tp = mdt.micro_fo_header_data_type;
    prm_tp = mdt.micro_fo_parameter_data_type;
    ini_tp = mdt.micro_fo_init_data_type;
    dat_tp = mdt.micro_fo_data_type;
    rst_tp = mdt.micro_fo_result_type;
  }
  void MultiscaleRVEAnalysis::initAnalysis()
  {
    int num_rve_tps = 1;
    cs->scaleBroadcast(M2m_id,&num_rve_tps);
    char ** rve_tp_dirs = new char * [num_rve_tps];
    MPI_Request rqsts[num_rve_tps];
    // the order of receipt might be non-deterministic. need to handle that
    for(int ii = 0; ii < num_rve_tps; ++ii)
    {
      int cnt = 0;
      while((cnt = cs->aRecvBroadcastSize<char>(M2m_id)) == 0)
      { }
      rve_tp_dirs[ii] = new char [cnt];
      // don't have to block to wait since we know the message was available for size info
      cs->aRecvBroadcast(&rqsts[ii],M2m_id,&rve_tp_dirs[ii][0],cnt);
    }
    MPI_Status stss[num_rve_tps];
    MPI_Waitall(num_rve_tps,&rqsts[0],&stss[0]);
    MPI_Request hdr_rqst;
    rve_tp_cnts.resize(num_rve_tps];
    cs->aRecvBroadcast(&hdr_rqst,M2m_id,&rve_tp_cnts[0],num_rve_tps);
    MPI_Status hdr_sts;
    MPI_Waitall(1,&hdr_rqst,&hdr_sts);
    // Read in all the fiber networks
    for(int ii = 0; ii < num_rve_tps; ++ii)
    {
      fnd.push_back(new FiberNetwork* [rve_tp_cnts[ii]]);
      sprs.push_back(new CSR* [rve_tp_cnts[ii]]);
    }
    int dof_max = -1;
    NetworkLoader ldr;
    for(int ii = 0; ii < num_rve_tps; ii++)
    {
      for(int jj = 0; jj < rve_tp_cnts[ii]; ++jj)
      {
        std::stringstream fl;
        fl << rve_tp_dirs[ii] << jj+1 << ".txt";
        std::ifstream strm(fl.str());
        FiberNetwork * fn = fbr_ntwrks[ii][jj] = ldr.fromStream(strm);
        int dofs = fn->getDofCount();
        sprs_strcts[ii][jj] = createCSR(fn->getUNumbering(),dofs)
        dof_max = dofs > dof_max ? dofs : dof_max;
      }
    }
    bfrs = new SparskitBuffers(dof_max);
  }
  void MultiscaleRVEAnalysis::updateCoupling()
  {
    std::vector<micro_fo_header> hdrs;
    std::vector<micro_fo_params> prms;
    std::vector<micro_fo_init_data> inis;
    std::vector<int> to_delete;
    amsi::ControlService * cs = amsi::ControlService::Instance();
    cs->RemoveData(recv_ptrn,to_delete);
    for(auto idx = to_delete.rbegin(); idx != to_delete.rend(); ++idx)
      destroyAnalysis(ans[*idx]);
    std::vector<int> to_add;
    std::vector<int> empty;
    size_t recv_init_ptrn = cs->addData(recv_ptrn,empty,to_add);
    int ii = 0;
    ans.resize(ans.size()+to_add.size());
    cs->Communicate(recv_init_ptrn, hdrs, hdr_tp);
    cs->Communicate(recv_init_ptrn, prms, prm_tp);
    cs->Communicate(recv_init_ptrn, inis, ini_tp);
    for(auto rve = ans.begin(); rve != ans.end(); ++rve)
    {
      if(*rve == NULL)
      {
        double pt[3];
        micro_fo_header & hdr = hdrs[ii];
        micro_fo_params & prm = prms[ii];
        micro_fo_init_data & dat = dats[ii];
        apf::Vector3 gp;
        apf::getGaussPoint(hdr.data[ELEMENT_TYPE],1,hdr.data[GAUSS_ID],gp);
        gp.toArray(pt);
        int nenodes = NumElementNodes(hdr.data[ELEMENT_TYPE]);
        int tp = hdr.data[RVE_TYPE];
        int rnd = rand() % rve_tp_cnt[tp];
        *rve = makeAnalysis(fns[tp][rnd]);
      }
    }
  }
  void MultiscaleRVEAnalysis::run()
  {
    amsi::ControlService * cs = amsi::ControlService::Instance();
    int sim_complete = 0;
    while(!sim_complete)
    {
      int step_complete = 0;
      while(!step_complete)
      {
        // migration
        std::vector<micro_fo_data> data;
        cs->Communicate(recv_ptrn,data,micro_fo_data_type);
        std::vector<micro_fo_result> results(data.size());
        int ii = 0;
        for(auto rve = ans.begin(); rve != ans.end(); ++rve)
        {
          applyMultiscaleCoupling(*rve,&data[ii]);
          FiberRVEIteration itr(*rve);
          FiberRVEConvergence cnvrg(*rve);
          amsi::nonlinearSolve(&itr,&cnvrg);
          recoverMultiscaleResults(*rve,&results[ii]);
          ++ii;
        }
        cs->Communicate(send_ptrn,results,micro_fo_result_type);
        macro_iter++;
        cs->scaleBroadcaste(M2m_id,&step_complete);
      }
      macro_iter = 0;
      macro_step++;
    }
  }
}
