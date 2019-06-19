#include "bioFiberRVEAnalysis.h"
#include "bioFiberRVEAnalysisQuasiStaticExplicit.h"
namespace bio {
  FiberRVEAnalysis * createFiberRVEAnalysis(
      FiberNetwork * fn,
      LinearStructs<las::MICRO_BACKEND> * vecs,
      const MicroSolutionStrategy & ss,
      FiberRVEAnalysisType type) {
    FiberRVEAnalysis * an = NULL; 
    if (type == FiberRVEAnalysisType::StaticImplicit) {
       an = new FiberRVEAnalysisSImplicit(fn, vecs, ss);
    }
    else if (type == FiberRVEAnalysisType::QuasiStaticExplicit) {
      an = new FiberRVEAnalysisQSExplicit(fn, vecs, ss);
    }
    else {
      std::cerr<<"Incorrect Analysis type chosen!"<<std::endl;
      std::abort();
    }
    return an;
  }
}
