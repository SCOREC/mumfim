#include "bioFiberRVEAnalysis.h"
#include "bioFiber.h"
#include "bioFiberNetwork.h"
#include "bioFiberNetworkIO.h"
#include <las.h>
#include <lasSys.h>
#include <lasConfig.h>
#include <apfMeshIterator.h>
#include <PCU.h>
namespace bio
{
  // todo: rename (shouldn't have reference to micro in a single-scale file)
  apf::Integrator * createMicroElementalSystem(FiberNetwork * fn,
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
                               k,
                               f,
                               1);
    return es;
  }
  FiberRVEAnalysis * makeFiberRVEAnalysis(FiberNetwork * fn,
                                          las::Sparsity * csr,
                                          las::SparskitBuffers * b)
  {
    FiberRVEAnalysis * an = new FiberRVEAnalysis;
    an->fn = fn;
    // todo determine rve size from input
    an->rve = new RVE(0.5,fn->getDim());
    auto bgn = amsi::apfMeshIterator(fn->getNetworkMesh(),0);
    decltype(bgn) end = amsi::apfEndIterator(fn->getNetworkMesh());
    for(int sd = RVE::side::bot; sd <= RVE::side::all; ++sd)
    {
      getBoundaryVerts(an->rve,
                       an->fn->getNetworkMesh(),
                       bgn,
                       end,
                       static_cast<RVE::side>(sd),
                       std::back_inserter(an->bnd_nds[sd]));
    }
    apf::Numbering * udofs = an->fn->getUNumbering();
    int ndofs = apf::NaiveOrder(udofs);
    las::LasCreateVec * vb = las::getVecBuilder<las::sparskit>(0);
    las::LasCreateMat * mb = las::getMatBuilder<las::sparskit>(0);
    an->f0 = vb->create(ndofs,LAS_IGNORE,MPI_COMM_SELF);
    an->f = vb->create(ndofs,LAS_IGNORE,MPI_COMM_SELF);
    an->u = vb->create(ndofs,LAS_IGNORE,MPI_COMM_SELF);
    an->k = mb->create(ndofs,LAS_IGNORE,csr,MPI_COMM_SELF);
    if(b == NULL)
      b = new las::SparskitBuffers(ndofs); // TODO memory leak (won't be hit in multi-scale)
    auto ops = las::getLASOps<las::sparskit>();
    ops->zero(an->f0);
    an->slv = las::createSparskitLUSolve(b,1e-6);
    an->es = createMicroElementalSystem(fn,an->k,an->f);
    an->dx_fn_dx_rve_set = false;
    return an;
  }
  void destroyAnalysis(FiberRVEAnalysis * fa)
  {
    auto md = las::getMatBuilder<las::sparskit>(0);
    auto vd = las::getVecBuilder<las::sparskit>(0);
    md->destroy(fa->k);
    vd->destroy(fa->u);
    vd->destroy(fa->f);
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
    /*
    PCU_Switch_Comm(MPI_COMM_SELF);
    if(!PCU_Comm_Self())
      apf::writeVtkFiles("rve",an->rve->getMesh());
    PCU_Switch_Comm(AMSI_COMM_SCALE);
    */
    auto ops = las::getLASOps<las::sparskit>();
    ops->zero(an->k);
    ops->zero(an->u);
    ops->zero(an->f);
    apf::Mesh * fn = an->fn->getNetworkMesh();
    apf::MeshEntity * me = NULL;
    apf::MeshIterator * itr = fn->begin(1);
    int ii = 0;
    std::stringstream sout;
    while((me = fn->iterate(itr)))
    {
      apf::MeshElement * mlm = apf::createMeshElement(fn,me);
      an->es->process(mlm);
      apf::destroyMeshElement(mlm);
      ++ii;
    }
    fn->end(itr);
    applyRVEBC(an->bnd_nds[RVE::all].begin(),
               an->bnd_nds[RVE::all].end(),
               an->fn->getUNumbering(),
               an->k,
               an->f);
    if(this->iteration() == 0)
      ops->axpy(1.0,an->f,an->f0);
    an->slv->solve(an->k,an->u,an->f);
    amsi::WriteOp wrt;
    amsi::AccumOp acm;
    amsi::FreeApplyOp fr_wrt(an->fn->getUNumbering(),&wrt);
    amsi::FreeApplyOp fr_acm(an->fn->getUNumbering(),&acm);
    double * s = NULL;
    ops->get(an->u,s);
    amsi::ApplyVector(an->fn->getUNumbering(),
                      an->fn->getdUField(),
                      s,0,&fr_wrt).run();
    amsi::ApplyVector(an->fn->getUNumbering(),
                      an->fn->getUField(),
                      s,0,&fr_acm).run();
    ops->restore(an->u,s);
    Iteration::iterate();
  }
  void calcStress(FiberRVEAnalysis * fra, apf::Matrix3x3 & sigma)
  {
    auto ops = las::getLASOps<las::sparskit>();
    for(int ii = 0; ii < 3; ++ii)
      for(int jj = 0; jj < 3; ++jj)
        sigma[ii][jj] = 0.0;
    double * f = NULL;
    ops->get(fra->f,f);
    auto & bnd = fra->bnd_nds[RVE::all];
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
    ops->restore(fra->f,f);
    // this is just the symmetric part of the matrix, there should be a standalone operation for this...
    sigma[0][1] = sigma[1][0] = 0.5 * (sigma[0][1] + sigma[1][0]);
    sigma[0][2] = sigma[2][0] = 0.5 * (sigma[0][2] + sigma[2][0]);
    sigma[1][2] = sigma[2][1] = 0.5 * (sigma[1][2] + sigma[2][1]);
  }
};
