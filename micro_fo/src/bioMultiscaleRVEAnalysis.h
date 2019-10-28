#ifndef BIO_MULTISCALE_RVE_ANALYSIS_H_
#define BIO_MULTISCALE_RVE_ANALYSIS_H_
#include <amsiMultiscale.h>
#include <amsiReporter.h>
#include <lasSparskit.h>
#include <vector>
#include "bioFiberRVEAnalysis.h"
#include "bioMultiscaleMicroFOParams.h"
namespace bio
{
  class MultiscaleRVEAnalysis
  {
    protected:
    // logging
    amsi::Log eff;
    amsi::Log wgt;
    amsi::Log tmg;
    amsi::Log rve_tp_lg;
    // multiscale
    size_t recv_ptrn;
    size_t send_ptrn;
    amsi::DataDistribution * rve_dd;
    size_t M2m_id;
    size_t m2M_id;
    // MicroFODatatypes dat_tp;
    // analysis
    std::vector<int> rve_tp_cnt;
    std::vector<FiberNetworkReactions **> fns;
    std::vector<apf::Mesh2 **> meshes;
    std::vector<micro_fo_header> hdrs;
    std::vector<micro_fo_params> prms;
    std::vector<micro_fo_solver> slvr_prms;
    std::vector<micro_fo_int_solver> slvr_int_prms;
    std::vector<las::Sparsity **> sprs;
    std::vector<RVEAnalysis *> ans;
    las::SparskitBuffers * bfrs;
    std::vector<LinearStructs<las::MICRO_BACKEND> *> vecs;
    // a vector that stores the dof number for each type of network
    std::vector<int *> dofs_cnt;
    int macro_iter;
    int macro_step;
    // funcs
    void initLogging();
    void initCoupling();
    void initAnalysis();
    void updateCoupling();

    public:
    MultiscaleRVEAnalysis();
    ~MultiscaleRVEAnalysis();
    virtual void init();
    virtual void run();
  };
}  // namespace bio
#endif
