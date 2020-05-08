#include "bioFiberRVEAnalysis.h"
#include "bioFiberRVEAnalysisExplicit.h"
#include "bioFiberRVEAnalysisStaticImplicit.h"
#include <memory>

namespace bio {
  std::unique_ptr<FiberRVEAnalysis> createFiberRVEAnalysis(
      std::unique_ptr<FiberNetwork> fn,
      std::unique_ptr<MicroSolutionStrategy> ss) {

    auto analysis = std::unique_ptr<FiberRVEAnalysis>{nullptr};
    if (ss->slvrType == SolverType::Implicit) {
       analysis.reset(new FiberRVEAnalysisSImplicit(std::move(fn), std::move(ss)));
    }
    else if (ss->slvrType == SolverType::Implicit) {
      analysis.reset(new FiberRVEAnalysisExplicit(std::move(fn), std::move(ss)));
    }
    else {
      std::cerr<<"Incorrect Analysis type chosen!"<<std::endl;
      std::abort();
    }
    return analysis;
  }
}
