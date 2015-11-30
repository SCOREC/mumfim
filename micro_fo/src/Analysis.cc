#include "Analysis.h"
namespace bio
{
  struct FiberRVEAnalysis
  {
    FiberNetwork * fn;
    RVE * rve;
    apf::Integrator * itgr;
    apf::FieldOp * writesol;
    apf::FieldOp * accumsol;
    las::skMat * k;
    las::skVec * u;
    las::skVec * f;
  };
  FiberRVEAnalysis * makeAnalysis(const std::string & fnm)
  {
    FiberRVEAnalysis * an = new FiberRVEAnalysis;
    an->fn = loadFromFile(fnm);
    an->rve = new RVE;
    int ndofs = apf::NaiveOrder(fn->getNumbering());
    an->f = las::makeVec(ndofs);
    an->u = las::makeVec(ndofs);
    an->k = new skMat(createCSR(dofs,ndofs));
    an->writesol = new ApplySolution(dofs,u,0,true);
    an->accumsol = new ApplySolution(dofs,u,0,false);
    an->itgr = new TrussIntegrator(1,dofs,new LinearReaction,an->k,an->f);
  }
  void destroyAnalysis(FiberRVEAnalysis * fa)
  {
    delete fa->itgr;
    delete fa->accumsol;
    delete fa->writesol;
    delete fa->k;
    las::destroyVec(fa->u);
    las::destroyVec(fa->f);
    delete fa->rve;
    delete fa->fn;
    delete fa;
    fa = NULL;
  }
  FiberRVEIteration::FiberRVEIteration(FiberRVEAnalysis * a)
    : num::Iteration()
    , an(a)
  {}
  void FiberRVEIteration::iterate()
  {
    an->itgr->process(an->fn->getNetworkMesh());
    applyRVEForceBC(an->f,an->rve,an->fn);
    //solve(an->k,an->u,an->f);
    an->accumsol->apply(an->fn->getIncrementalDispField());
    an->writesol->apply(an->fn->getDisplacementField());
  }
  FiberRVEConvergence::FiberRVEConvergence(FiberRVEAnalysis * a, double e)
    : num::Convergence()
    , an(a)
    , eps(e)
    , resid_im(0.0)
  {}
  bool FiberRVEConvergece::converged()
  {
    bool result = false;
    double resid = norm(*(an->u));
    if(resid < eps)
      result = true;
    return reslt;
  }
};