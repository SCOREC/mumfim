#include "bioFiberRVEAnalysis.h"
#include "bioVerbosity.h"
namespace bio
{
  struct val_gen
  {
    val_gen(FiberRVEAnalysis * a) : an(a), prv_nrm(1.0) {}
    FiberRVEAnalysis * an;
    double prv_nrm;
    double operator()()
    {
      auto ops = las::getLASOps<las::MICRO_BACKEND>();
      double nrm = ops->norm(an->getF());
      double val = fabs(prv_nrm - nrm);
      prv_nrm = nrm;
      return val;
    }
  };
  struct eps_gen
  {
    eps_gen(double eps) : eps(eps) {}
    double operator()(int) { return eps; }

    protected:
    double eps;
  };
  struct ref_gen
  {
    double operator()() { return 1.0; }
  };
  FiberRVEAnalysisSImplicit::FiberRVEAnalysisSImplicit(
      const FiberRVEAnalysisSImplicit & an)
      : FiberRVEAnalysis(an)
  {
    es = createMicroElementalSystem(fn, getK(), getF(), FiberRVEAnalysisType::StaticImplicit);
  }
  FiberRVEAnalysisSImplicit::FiberRVEAnalysisSImplicit(FiberNetwork * fn,
                              LinearStructs<las::MICRO_BACKEND> * vecs,
                              const MicroSolutionStrategy & ss) : FiberRVEAnalysis(fn, vecs, ss) {
    es = createMicroElementalSystem(fn, getK(), getF(), FiberRVEAnalysisType::StaticImplicit);
  }
  FiberRVEIterationSImplicit::FiberRVEIterationSImplicit(FiberRVEAnalysisSImplicit * a)
      : amsi::Iteration(), an(a)
  {
  }
  void FiberRVEIterationSImplicit::iterate()
  {
    auto ops = las::getLASOps<las::MICRO_BACKEND>();
    ops->zero(an->getK());
    ops->zero(an->getU());
    ops->zero(an->getF());
    apf::Mesh * fn = an->getFn()->getNetworkMesh();
    an->es->process(fn, 1);
    // finalize the vectors so we can set boundary condition
    // values
    las::finalizeMatrix<las::MICRO_BACKEND>(an->vecs->getK());
    las::finalizeVector<las::MICRO_BACKEND>(an->vecs->getF());
    applyRVEBC(an->bnd_nds[RVE::all].begin(),
               an->bnd_nds[RVE::all].end(),
               an->getFn()->getUNumbering(),
               an->getK(),
               an->getF());
    an->getSlv()->solve(an->getK(), an->getU(), an->getF());
    amsi::WriteOp wrt;
    amsi::AccumOp acm;
    amsi::FreeApplyOp fr_wrt(an->getFn()->getUNumbering(), &wrt);
    amsi::FreeApplyOp fr_acm(an->getFn()->getUNumbering(), &acm);
    double * s = NULL;
    ops->get(an->getU(), s);
    amsi::ApplyVector(
        an->getFn()->getUNumbering(), an->getFn()->getdUField(), s, 0, &fr_wrt)
        .run();
    amsi::ApplyVector(
        an->getFn()->getUNumbering(), an->getFn()->getUField(), s, 0, &fr_acm)
        .run();
    ops->restore(an->getU(), s);
    BIO_V2(int rank = -1; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
           std::cout << "RVE: " << an->getFn()->getRVEType() << " Iter: "
                     << this->iteration() << " Rank: " << rank << "\n";)
    Iteration::iterate();
  }
  bool FiberRVEAnalysisSImplicit::run(const DeformationGradient & dfmGrd)
  {
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    unsigned int maxMicroAttempts = 0;   // parameter
    unsigned int microAttemptCutFactor;  // parameter
    bool solveSuccess = false;
    unsigned int microAttemptCount = 1;
    unsigned int attemptCutFactor;
    do
    {
      // create a deep copy of the analysis
      // Note the current implementation of copy does not deep copy
      // the sparskit matrices, vectors, or solver
      auto tmpRVE = new FiberRVEAnalysisSImplicit(*this);
      val_gen vg(tmpRVE);
      eps_gen eg(tmpRVE->solver_eps);
      ref_gen rg;
      maxMicroAttempts = tmpRVE->max_cut_attempt;
      microAttemptCutFactor = tmpRVE->attempt_cut_factor;
      attemptCutFactor = std::pow(microAttemptCutFactor, microAttemptCount - 1);
      BIO_V1(if (attemptCutFactor > 1) std::cout
                 << "Micro Attempt: " << microAttemptCount
                 << " cutting the original applied displacement by: "
                 << attemptCutFactor << " on rank: " << rank << "\n";)
      assert(maxMicroAttempts > 0);
      DeformationGradient appliedDefm;
      bool microIterSolveSuccess = true;
      for (unsigned int microAttemptIter = 1;
           microAttemptIter <= attemptCutFactor; ++microAttemptIter)
      {
        BIO_V3(std::cout << "Rank: " << rank << " F=";)
        for (int j = 0; j < 9; ++j)
        {
          double I = j % 4 == 0 ? 1 : 0;
          // cut the displacement gradient (not the deformation gradient)
          appliedDefm[j] =
              (((dfmGrd[j] - I) * microAttemptIter) / attemptCutFactor) + I;
          BIO_V3(std::cout << appliedDefm[j] << " ";)
        }
        BIO_V3(std::cout << "\n";)
        applyGuessSolution(tmpRVE, appliedDefm);
        FiberRVEIterationSImplicit rveItr(tmpRVE);
        std::vector<amsi::Iteration *> itr_stps;
        itr_stps.push_back(&rveItr);
        amsi::MultiIteration itr(itr_stps.begin(), itr_stps.end());
        amsi::UpdatingConvergence<decltype(&vg), decltype(&eg), decltype(&rg)>
            resid_cnvrg(&itr, &vg, &eg, &rg);
        std::vector<amsi::Convergence *> cnvrg_stps;
        cnvrg_stps.push_back(&resid_cnvrg);
        amsi::MultiConvergence cnvrg(cnvrg_stps.begin(), cnvrg_stps.end());
        amsi::Iteration * osc_itr =
            amsi::createOscillationDetection<decltype(&resid_cnvrg)>(
                tmpRVE->detect_osc_type, &resid_cnvrg, &rveItr,
                tmpRVE->max_itrs, tmpRVE->prev_itr_factor);
        itr.addIteration(osc_itr);
        // solve is successful if the current solve and all previous
        // cutIterations are successful
        microIterSolveSuccess =
            (amsi::numericalSolve(&itr, &cnvrg) && microIterSolveSuccess);
        // cleanup the oscillation detection memory
        delete osc_itr;
        // don't bother computing the rest of the attempt if any
        // subiteration fails, for our current use case we don't care what
        // made us fail, we will try to reduce the load and try again.
        if (!microIterSolveSuccess) break;
      }
      // if the attempt was completely successful then the overall solve
      // was successful
      if (microIterSolveSuccess)
      {
        solveSuccess = true;
        // note that the destructor for *this should get called automatically
        *this = *tmpRVE;
        tmpRVE = NULL;
      }
      ++microAttemptCount;
    } while (solveSuccess == false && (microAttemptCount <= maxMicroAttempts));
    if (!solveSuccess)
    {
      std::cerr << "RVE: " << this->getFn()->getRVEType()
                << " failed to converge in " << microAttemptCount - 1
                << " attempts on processor " << rank << ".\n";
      return false;
    }
    return true;
  }
}  // namespace bio
