#ifndef BIO_SOLVER_H_
#define BIO_SOLVER_H_

#include "FiberNetwork.h"
#include "FiberReactions.h"
#include "RVE.h"
#include "TrussIntegrator.h"

namespace bio
{

  class NewtonRaphson
  {
  protected:
    int load_step;
    int iteration;
    double epsilon;
  public:
    virtual void run(int cp)
    {
      for(load_step = 0; load_step < cp; load_step++)
	step();
    }
    virtual void step()
    {
      while(!iter())
      {
	load_step++;
      }
    }
    virtual bool iter()
    {
      iteration++;
      return true;
    }
  };

  class FiberRVEAnalysis : public NewtonRaphson
  {
  protected:
    FiberNetwork * fn;
    RVE * rve;
    apf::Integrator * itgr;

    apf::FieldOp * writesol;
    apf::FieldOp * accumsol;
    
    skVec f;
    skVec u;
    skMat * k;
  public:
    FiberRVEAnalysis(FiberNetwork * f, RVE * r)
      : fn(f)
      , rve(r)
      , itgr(NULL)
      , writesol(NULL)
      , accumsol(NULL)
      , f()
      , u()
      , k()
    { }
    void init()
    {
      apf::Numbering * dofs = fn->getNumbering();
      // need to use naive order... adjreorder tries
      // to query faces which don't exist for the fiber network
      int ndofs = apf::NaiveOrder(dofs);
      f = makeVec(ndofs);
      u = makeVec(ndofs);
      k = new skMat(createCSR(dofs,ndofs));
      writesol = new ApplySolution(dofs,u,0,true);
      accumsol = new ApplySolution(dofs,u,0,false);
      itgr = new TrussIntegrator(1,dofs,new LinearReaction,k,&f);
    }
    void step()
    {
      NewtonRaphson::step();
    }
    bool iter()
    {
      itgr->process(fn->getNetworkMesh());
      applyRVEForceBC(&f,rve,fn);
      //solve(k,u,f);
      accumsol->apply(fn->getIncrementalDispField());
      writesol->apply(fn->getDisplacementField());
      /*
      if(norm(u) < epsilon)
	return true;
      else
      {
	NewtonRaphson::iter();
	return false;
      }
      */
    }
  };
  void RVEanalysis(RVE * rve, FiberNetwork * fn, double epsilon)
  {
    /*
    apf::Numbering * dofs = fn->getNumbering();
    int ndofs = apf::AdjReorder(dofs);
    skVec f = makeVec(ndofs);
    skMat k(makeStructure(dofs,ndofs));
    TrussIntegrator t(1,num,new LinearReaction,k,f);

    t.process(fn->getMesh());
    applyRVEForceBC(f,rve,fn);
    skVec u = solve(k,f);
    ApplySolution(dof,u).apply(fn->getIncrementalDispField());
    */
  }
  /*
  void multiscaleRVEAnalysis(MacroInfo * macro,
			     RVE * rve,
			     FiberNetwork * fn)
  {
    // use macroscale information to setup initial displacements

    ///RVEAnalysis(rve,fn,1e-12);

    // compute the macroscale stress terms
  }
  */
}

#endif
