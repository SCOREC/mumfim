#ifndef BIO_FIBER_RVE_ANALYSIS_STATIC_IMPLICIT_H__
#define BIO_FIBER_RVE_ANALYSIS_STATIC_IMPLICIT_H__

#include "bioFiberRVEAnalysis.h"

namespace bio {
  class FiberRVEAnalysisSImplicit : public FiberRVEAnalysis
  {
    protected:
    virtual void computeCauchyStress(double sigma[6]) final;

    public:
    // constructors
    //explicit FiberRVEAnalysisSImplicit(const FiberRVEAnalysisSImplicit & an);
    FiberRVEAnalysisSImplicit(std::unique_ptr<FiberNetwork> fn,
                              std::unique_ptr<MicroSolutionStrategy> ss,
                              las::SparskitBuffers* sparksit_workspace);
    FiberRVEAnalysisSImplicit(FiberRVEAnalysisSImplicit && an) = default;
    virtual bool run(const DeformationGradient & dfmGrd, double sigma[6], bool update_coords=true) final;
    virtual SolverType getAnalysisType()
    {
      return SolverType::Implicit;
    }
    virtual void computeStiffnessMatrix()
    {
      auto ops = las::getLASOps<las::MICRO_BACKEND>();
      if(getK() == nullptr || getF() == nullptr)
      {
        std::cerr<<"ERROR! The stiffness matrix, and force vector should not be null!"<<std::endl;
        std::abort();
      }
      ops->zero(getK());
      ops->zero(getF());
      es->process(getFn()->getNetworkMesh(), 1);
      // finalize the vectors so we can set boundary condition
      // values
      las::finalizeMatrix<las::MICRO_BACKEND>(getK());
      las::finalizeVector<las::MICRO_BACKEND>(getF());
    }
  };
  class FiberRVEIterationSImplicit : public amsi::Iteration
  {
    protected:
    FiberRVEAnalysisSImplicit * an;

    public:
    FiberRVEIterationSImplicit(FiberRVEAnalysisSImplicit * a);
    void iterate();
  };
}
#endif
