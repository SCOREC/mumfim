#include "bioFiberRVEAnalysis.h"
#include "bioFiberRVEAnalysisExplicit.h"
#include "bioFiberRVEAnalysisStaticImplicit.h"
#include <memory>

namespace bio {
  std::unique_ptr<FiberRVEAnalysis> createFiberRVEAnalysis(
      std::unique_ptr<FiberNetwork> fn,
      const MicroSolutionStrategy & ss, void * bfrs) {

    LinearStructs<las::MICRO_BACKEND> * vecs,
    FiberRVEAnalysis * an = nullptr; 
    
    if (ss.slvrType == SolverType::Implicit) {
       an = new FiberRVEAnalysisSImplicit(fn, vecs, ss);
    }
    else if (ss.SolverType == SolverType::Implicit) {
      an = new FiberRVEAnalysisExplicit(fn, vecs, static_cast<const MicroSolutionStrategyExplicit &>(ss));
    }
    else {
      std::cerr<<"Incorrect Analysis type chosen!"<<std::endl;
      std::abort();
    }
    return an;
  }
}
