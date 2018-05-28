#ifndef BIO_MULTISCALE_RVE_ANALYSIS_H_
#define BIO_MULTISCALE_RVE_ANALYSIS_H_
#include "bioFiberRVEAnalysis.h"
#include "bioMicroFOMultiscale.h"
#include "lasSparskit.h"
#include <amsiReporter.h>
#include <amsiMultiscale.h>
namespace bio
{
  void applyMultiscaleCoupling(FiberRVEAnalysis * ans, micro_fo_data * data);
  void recoverMultiscaleResults(FiberRVEAnalysis * ans, micro_fo_result * data);
  class MultiscaleRVEAnalysis
  {
  protected:
    // logging
    amsi::Log eff;
    amsi::Log wgt;
    amsi::Log tmg;
    // multiscale
    size_t recv_ptrn;
    size_t send_ptrn;
    amsi::DataDistribution * rve_dd;
    size_t M2m_id;
    size_t m2M_id;
    //MicroFODatatypes dat_tp;
    // analysis
    std::vector<int> rve_tp_cnt;
    std::vector<FiberNetworkReactions **> fns;
    std::vector<las::CSR **> sprs;
    std::vector<FiberRVEAnalysis*> ans;
    las::SparskitBuffers * bfrs;
    int macro_iter;
    int macro_step;
    // funcs
    void initLogging();
    void initCoupling();
    void initAnalysis();
    void updateCoupling();
  public:
    MultiscaleRVEAnalysis();
    virtual void init();
    virtual void run();
  };
}
#include "bioMultiscaleRVEAnalysis_impl.h"
#endif
