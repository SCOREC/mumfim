#include "bioFiberRVEAnalysis.h"
namespace bio
{
  FiberRVEAnalysisQSExplicit::FiberRVEAnalysisQSExplicit(
      const FiberRVEAnalysisSImplicit & an)
      : FiberRVEAnalysis(an)
  {
  }
  FiberRVEAnalysisQSExplicit::FiberRVEAnalysisQSExplicit(FiberNetwork * fn,
                              LinearStructs<las::MICRO_BACKEND> * vecs,
                              const MicroSolutionStrategy & ss) : FiberRVEAnalysis(fn, vecs, ss) {
  }
  FiberRVEIterationQSExplicit::FiberRVEIterationQSExplicit(FiberRVEAnalysisQSExplicit * a)
      : amsi::Iteration(), an(a)
  {
  }
  // see belytschko box 6.1
  void FiberRVEIterationQSExplicit::iterate() {
    if(this->iteration() == 0) {
      // ONLY FIRST ITERATION
      // 1. initial conditions and initialization
      // 2. getforce
      // 3. compute accelerations a(n) = inv(M)(f-Cv(n-1/2))
    }
    // ITERATIONS START HERE!
    // 4. update time t(n+1) = t(n)+delta_t(n+1/2), t(n+1/2)=1/2(t(n)+t(n+1))
    // 5. partial update of nodal velocities v(n+1/2)=v(n)+(t(n+1/2)-t(n))a(n)
    // 6. enforce dirichlet boundary conditions
    // 7. update nodal displacements u(n+1) = delta_t(n+1/2)*v(n+1/2)
    // 8. getforce
    // 9. compute a(n+1)
    // 10. partial update of nodal velocities v(n+1)=v(n+1/2)+(t(n+1)-t(n+1/2))a(n+1)
    // 11. check energy balance
    Iteration::iterate();
  }
  bool FiberRVEAnalysisQSExplicit::run(
      const DeformationGradient & dfmGrd)
  {
    // While "time" < 1
    //     run iteration
    //     check convergence
    //         is balance of energy preserved?
    //         is the kinetic energy < threshold*total energy? (are we
    //         quasistatic)
    return false;
  }
}  // namespace bio
