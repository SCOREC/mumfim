#ifndef BIO_MULTISCALE_RVE_ANALYSIS_H_
#define BIO_MULTISCALE_RVE_ANALYSIS_H_
#include <amsiMultiscale.h>
#include <amsiReporter.h>
#include <vector>
#include "bioFiberRVEAnalysis.h"
#include "bioMultiscaleMicroFOParams.h"
#include "bioFiberNetworkLibrary.h"
namespace bio
{
  class MultiscaleRVEAnalysis
  {
    private:
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
   
    // number of different RVE realizations of a given type
    // each different type (physical set of parameters) is
    // stored in a different folder.
    std::vector<int> rve_tp_cnt;

    std::vector<micro_fo_header> hdrs;
    std::vector<micro_fo_params> prms;
    std::vector<micro_fo_solver> slvr_prms;
    std::vector<micro_fo_int_solver> slvr_int_prms;

    std::vector<std::unique_ptr<RVEAnalysis>> ans;
    // this is a workspace that the microscale solves can use.
    // Note the current design is that all microscale solves
    // share the same workspace
    std::shared_ptr<void> mWorkspace;
    //  an vector that holds the names of the rve physical
    //  types. Here a type would refer to a set of networks
    //  with the same set of statistical properties i.e. density
    std::vector<char*> rve_tp_dirs;
    FiberNetworkLibrary network_library;
   
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

  std::unique_ptr<MicroSolutionStrategy> serializeSolutionStrategy(micro_fo_solver & slvr, micro_fo_int_solver &slvr_int);

}  // namespace bio
#endif
