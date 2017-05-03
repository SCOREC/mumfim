#include "bioFiberRVEAnalysis.h"
#include "bioFiber.h"
#include "bioFiberNetwork.h"
#include "bioFiberNetworkIO2.h"
#include "bioRVE2.h"
#include "bioTrussIntegrator.h"
#include "lasSparskit.h"
namespace bio
{
  amsi::ElementalSystem * createMicroElementalSystem(FiberNetwork * fn)
  {
    amsi::ElementalSystem * es = NULL;
    FiberMemeber tp = fn->getFiberMember();
    if(tp == FiberMemeber::truss)
      es = new TrussIntegrator();
  }
  FiberRVEAnalysis * makeFiberRVEAnalysis(FiberNetwork * fn)
  {
    FiberRVEAnalysis * an = new FiberRVEAnalysis;
    an->fn = fn;
    an->rve = new RVE(fn->getDim());
    getBoundaryVerts(an->rve,an->fn,std::back_inserter(bnd_nds));
    apf::Numbering * udofs = an->fn->getUNumbering();
    int nudofs = apf::NaiveOrder(udofs);
    apf::Numbering * wdofs = an->fn->getdWNumbering();
    int nwdofs = apf::NaiveOrder(wdofs);
    apf::SetNumberingOffset(wdof,nudofs);
    int ndofs = nudofs + nwdofs;
    an->es = createMicroElementalSystem(fn);
    an->f = las::createSparskitVec(udofs);
    an->u = las::createSparskitVec(udofs);
    // todo : won't work when we have a field for moments.
    an->k = las::createSparskitMat(las::createCSR(udofs,nudofs));
    an->ops = las::initSparskitOps();
  }
  void destroyAnalysis(FiberRVEAnalysis * fa)
  {
    delete fa->k;
    las::destroyVec(fa->u);
    las::destroyVec(fa->f);
    fa->rve = NULL;
    fa->fn = NULL;
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
    applyRVEForceBC(an->bnd_nds.begin(),
                    an->bnd_nds.end(),
                    an->fn->getUNumbering(),
                    an->ops,
                    an->f);
    an->ops->solve(an->k,&an->u,&an->f);
    amsi::WriteOp wrt;
    amsi::AccumOp acm;
    amsi::FreeApplyOp fr_wrt(an->fn->getUNumbering(),&wrt);
    amsi::FreeApplyOp fr_acm(an->fn->getUNumbering(),&acm);
    double * s = NULL;
    an->ops->get(&an->u,s);
    amsi::ApplyVector(an->fn->getUNumbering(),
                      an->fn->getdUField(),
                      s,0,&fr_wrt).run();
    amsi::ApplyVector(an->fn->getUNumbering(),
                      an->fn->getUField(),
                      s,0,fr_acm).run();
    an->ops->restore(&an->u,s);
  }
  FiberRVEConvergence::FiberRVEConvergence(FiberRVEAnalysis * a, double e)
    : num::Convergence()
    , an(a)
    , eps(e)
    , resid_im(0.0)
  {}
  bool FiberRVEConvergence::converged()
  {
    bool result = false;
    double resid = 0.0;
    //double resid = norm(*(an->u));
    if(resid < eps)
      result = true;
    return result;
  }
  double calcStiffness(FiberRVEAnalysis *)
  {
    return 0.0;
  }
  void calcStress(FiberRVEAnalysis * fra, apf::Matrix3 & sigma)
  {
    sigma.zero();
    std::vector<apf::MeshEntity*> bnd;
    getBoundaryVerts(fra->rve,fra->fn,RVE::side::all,std::back_inserter(bnd));
    apf::Numbering * nm = fra->fn->getNumbering();
    apf::Field * u = fra->fn->getDisplacementField();
    apf::Mesh * fn = fra->fn->getNetworkMesh();
    for(auto vrt = bnd.begin(); vrt != bnd.end(); ++vrt)
    {
      int dof[3] = {};
      dof[0] = apf::getNumber(nm,*vrt,0,0);
      dof[1] = apf::getNumber(nm,*vrt,0,1);
      dof[2] = apf::getNumber(nm,*vrt,0,2);
      apf::Vector3 x;
      fn->getPoint(*vrt,0,x);
      apf::Vector3 u;
      apf::getVector(u,*vrt,0,u);
      x += u;
      for(int ii = 0; ii < 3; ++ii)
        for(int jj = 0; jj < 3; ++jj)
          sigma(ii,jj) += fra->f(dof[ii]) * x(jj);
    }
    // this is just the symmetric part of the matrix, there should be a standalone operation for this...
    sigma(0,1) = sigma(1,0) = 0.5 * (sigma(0,1) + sigma(1,0));
    sigma(0,2) = sigma(2,0) = 0.5 * (sigma(0,2) + sigma(2,0));
    sigma(1,2) = sigma(2,1) = 0.5 * (sigma(1,2) + sigma(2,1));
  }
  void tdYdXr(FiberRVEAnalysis * fra)
  {

  }
  void calcFEMJacobian(FiberRVEAnalysis * fra)
  {

  }
  void solver(FiberRVEAnalysis * fra)
  {
    
  }
};
