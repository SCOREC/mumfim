#include "bioFiberRVEAnalysis.h"
#include "bioFiber.h"
#include "bioFiberNetwork.h"
#include "bioFiberNetworkIO2.h"
#include "bioRVE2.h"
#include "bioTrussIntegrator.h"
#include "lasSparskit.h"
namespace bio
{
  apf::Integrator * createMicroElementalSystem(FiberNetwork * fn, las::LasOps * ops, las::Mat * k, las::Vec * f)
  {
    apf::Integrator * es = NULL;
    FiberMember tp = fn->getFiberMember();
    if(tp == FiberMember::truss)
      es = new TrussIntegrator(fn->getUNumbering(),fn->getFiberReactions(),ops,k,f,1);
    return es;
  }
  FiberRVEAnalysis * makeFiberRVEAnalysis(FiberNetwork * fn, las::SparskitBuffers * b)
  {
    FiberRVEAnalysis * an = new FiberRVEAnalysis;
    an->fn = fn;
    an->rve = new RVE(fn->getDim());
    getBoundaryVerts(an->rve,an->fn,RVE::side::all,std::back_inserter(an->bnd_nds));
    apf::Numbering * udofs = an->fn->getUNumbering();
    int nudofs = apf::NaiveOrder(udofs);
    apf::Numbering * wdofs = an->fn->getdWNumbering();
    int nwdofs = apf::NaiveOrder(wdofs);
    apf::SetNumberingOffset(wdofs,nudofs);
    int ndofs = nudofs + nwdofs;
    an->f = las::createSparskitVector(ndofs);
    an->u = las::createSparskitVector(ndofs);
    // todo : won't work when we have a field for moments.
    an->k = las::createSparskitMatrix(las::createCSR(udofs,nudofs));
    if(b == NULL)
      b = new las::SparskitBuffers(ndofs); // TODO memory leak (won't be hit in multi-scale)
    an->ops = las::initSparskitOps(b);
    an->es = createMicroElementalSystem(fn,an->ops,an->k,an->f);
    return an;
  }
  void destroyAnalysis(FiberRVEAnalysis * fa)
  {
    las::deleteSparskitMatrix(fa->k);
    las::deleteSparskitVector(fa->u);
    las::deleteSparskitVector(fa->f);
    delete fa->rve;
    fa->fn = NULL;
    delete fa;
  }
  FiberRVEIteration::FiberRVEIteration(FiberRVEAnalysis * a)
    : amsi::Iteration()
    , an(a)
  {}
  void FiberRVEIteration::iterate()
  {
    // need to apply dirichlet bcs possibly?
    an->es->process(an->fn->getNetworkMesh());
    applyRVEForceBC(an->bnd_nds.begin(),
                    an->bnd_nds.end(),
                    an->fn->getUNumbering(),
                    an->ops,
                    an->f,
                    an->k);
    an->ops->solve(an->k,an->u,an->f);
    amsi::WriteOp wrt;
    amsi::AccumOp acm;
    amsi::FreeApplyOp fr_wrt(an->fn->getUNumbering(),&wrt);
    amsi::FreeApplyOp fr_acm(an->fn->getUNumbering(),&acm);
    double * s = NULL;
    an->ops->get(an->u,s);
    amsi::ApplyVector(an->fn->getUNumbering(),
                      an->fn->getdUField(),
                      s,0,&fr_wrt).run();
    amsi::ApplyVector(an->fn->getUNumbering(),
                      an->fn->getUField(),
                      s,0,&fr_acm).run();
    an->ops->restore(an->u,s);
  }
  FiberRVEConvergence::FiberRVEConvergence(FiberRVEAnalysis * a, double e)
    : amsi::Convergence()
    , an(a)
    , eps(e)
    , resid_im(0.0)
  {}
  bool FiberRVEConvergence::converged()
  {
    // don't bother recalculating force vector, just overconverge to get more accuracy
    bool result = false;
    double resid = an->ops->norm(an->f);
    if(resid < eps)
      result = true;
    return result;
  }
  void calcStress(FiberRVEAnalysis * fra, apf::Matrix3x3 & sigma)
  {
    for(int ii = 0; ii < 3; ++ii)
      for(int jj = 0; jj < 3; ++jj)
        sigma[ii][jj] = 0.0;
    double * f = NULL;
    fra->ops->get(fra->f,f);
    std::vector<apf::MeshEntity*> bnd;
    getBoundaryVerts(fra->rve,fra->fn,RVE::side::all,std::back_inserter(bnd));
    apf::Numbering * nm = fra->fn->getUNumbering();
    apf::Field * uf = fra->fn->getUField();
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
      apf::getVector(uf,*vrt,0,u);
      x += u;
      for(int ii = 0; ii < 3; ++ii)
        for(int jj = 0; jj < 3; ++jj)
          sigma[ii][jj] += f[dof[ii]] * x[jj];
    }
    fra->ops->restore(fra->f,f);
    // this is just the symmetric part of the matrix, there should be a standalone operation for this...
    sigma[0][1] = sigma[1][0] = 0.5 * (sigma[0][1] + sigma[1][0]);
    sigma[0][2] = sigma[2][0] = 0.5 * (sigma[0][2] + sigma[2][0]);
    sigma[1][2] = sigma[2][1] = 0.5 * (sigma[1][2] + sigma[2][1]);
  }
};
