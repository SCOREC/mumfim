#include "bioFiberRVEAnalysis.h"
#include "bioFiber.h"
#include "bioFiberNetwork.h"
#include "bioFiberNetworkIO.h"
#include "lasSparskit.h"
#include <apfMeshIterator.h>
#include <PCU.h>
namespace bio
{
  // todo: rename (shouldn't have reference to micro in a single-scale file)
  apf::Integrator * createMicroElementalSystem(FiberNetwork * fn,
                                               las::LasOps * ops,
                                               las::Mat * k,
                                               las::Vec * f)
  {
    apf::Integrator * es = NULL;
    FiberMember tp = fn->getFiberMember();
    if(tp == FiberMember::truss)
      es = new TrussIntegrator(fn->getUNumbering(),
                               fn->getUField(),
                               fn->getXpUField(),
                               &fn->getFiberReactions()[0],
                               ops,
                               k,
                               f,
                               1);
    return es;
  }
  FiberRVEAnalysis * makeFiberRVEAnalysis(FiberNetwork * fn,
                                          las::CSR * csr,
                                          las::SparskitBuffers * b)
  {
    FiberRVEAnalysis * an = new FiberRVEAnalysis;
    an->fn = fn;
    // todo determine rve size from input
    an->rve = new RVE(0.5,fn->getDim());
    auto bgn = amsi::apfMeshIterator(fn->getNetworkMesh(),0);
    decltype(bgn) end = amsi::apfEndIterator(fn->getNetworkMesh());
    getBoundaryVerts(an->rve,
                     an->fn->getNetworkMesh(),
                     bgn,
                     end,
                     RVE::side::all,
                     std::back_inserter(an->bnd_nds));
    apf::Numbering * udofs = an->fn->getUNumbering();
    int ndofs = apf::NaiveOrder(udofs);
    an->f0 = las::createSparskitVector(ndofs);
    an->f = las::createSparskitVector(ndofs);
    an->u = las::createSparskitVector(ndofs);
    an->k = las::createSparskitMatrix(csr); // no direct connection between csr and u field internal to this function, might need to refactor for that
    if(b == NULL)
      b = new las::SparskitBuffers(ndofs); // TODO memory leak (won't be hit in multi-scale)
    an->ops = las::initSparskitOps();
    an->ops->zero(an->f0);
    an->slv = las::createSparskitLUSolve(b);
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
    // free fixed bcs so all dofs get numbers
    freeRVEBC(an->bnd_nds.begin(),
              an->bnd_nds.end(),
              an->fn->getUNumbering());
    // number all dofs
    apf::NaiveOrder(an->fn->getUNumbering());
    an->ops->zero(an->k);
    an->ops->zero(an->u);
    an->ops->zero(an->f);
    // set the rows/cols for each boundary node
    //  to identity rows, and fix the dofs
    //  themselves so all assembly routines
    //  don't write to the fixed rows, leaving
    //  them as identity rows
    applyRVEBC(an->bnd_nds.begin(),
               an->bnd_nds.end(),
               an->fn->getUNumbering(),
               an->ops,
               an->k,
               an->f);
    apf::Mesh * fn = an->fn->getNetworkMesh();
    apf::MeshEntity * me = NULL;
    apf::MeshIterator * itr = fn->begin(1);
    while((me = fn->iterate(itr)))
    {
      apf::MeshElement * mlm = apf::createMeshElement(fn,me);
      an->es->process(mlm);
      apf::destroyMeshElement(mlm);
    }
    fn->end(itr);
    if(!PCU_Comm_Self())
    {
      std::ofstream fout("micro_fo_matrix_0");
      printSparskitMat(fout,an->k);
    }
    if(this->iteration() == 0)
      an->ops->axpy(1.0,an->f,an->f0);
    an->slv->solve(an->k,an->u,an->f);
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
    Iteration::iterate();
  }
  void calcStress(FiberRVEAnalysis * fra, apf::Matrix3x3 & sigma)
  {
    for(int ii = 0; ii < 3; ++ii)
      for(int jj = 0; jj < 3; ++jj)
        sigma[ii][jj] = 0.0;
    double * f = NULL;
    fra->ops->get(fra->f,f);
    auto & bnd = fra->bnd_nds;
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
