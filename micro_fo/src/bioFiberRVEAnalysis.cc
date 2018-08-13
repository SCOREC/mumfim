#include "bioFiberRVEAnalysis.h"
#include "bioFiber.h"
#include "bioFiberNetwork.h"
#include "bioFiberNetworkIO.h"
#include <las.h>
#include <lasSys.h>
#include <lasConfig.h>
#include <apfMeshIterator.h>
#include <PCU.h>
#include <bioVerbosity.h>
namespace bio
{
  // todo: rename (shouldn't have reference to micro in a single-scale file)
  static apf::Integrator * createMicroElementalSystem(FiberNetwork * fn,
                                               las::Mat * k,
                                               las::Vec * f)
  {
    apf::Integrator * es = NULL;
    FiberMember tp = fn->getFiberMember();
    if(tp == FiberMember::truss)
      es = new TrussIntegrator(fn->getUNumbering(),
                               fn->getUField(),
                               fn->getXpUField(),
                               &(fn->getFiberReactions()[0]),
                               k,
                               f,
                               1);
    return es;
  }
  FiberRVEAnalysisVecs::FiberRVEAnalysisVecs(int ndofs, las::Sparsity *csr,
                                             las::SparskitBuffers *b)
  {
    las::LasCreateVec *vb = las::getVecBuilder<las::sparskit>(0);
    las::LasCreateMat *mb = las::getMatBuilder<las::sparskit>(0);
    this->f = vb->create(ndofs, LAS_IGNORE, MPI_COMM_SELF);
    this->u = vb->create(ndofs, LAS_IGNORE, MPI_COMM_SELF);
    this->k = mb->create(ndofs, LAS_IGNORE, csr, MPI_COMM_SELF);
    // FIXME memory leak (won't be hit in multi-scale)
    // because we provide buffers
    if (b == NULL)
      b = new las::SparskitBuffers(
          ndofs);  
    this->slv = las::createSparskitLUSolve(b, 1e-6);
  }
  FiberRVEAnalysisVecs::~FiberRVEAnalysisVecs() {
      auto md = las::getMatBuilder<las::sparskit>(0);
      auto vd = las::getVecBuilder<las::sparskit>(0);
      md->destroy(k);
      vd->destroy(u);
      vd->destroy(f);
      delete slv;
  }
  FiberRVEAnalysis::FiberRVEAnalysis(const FiberRVEAnalysis &an)
  {
    fn = new FiberNetwork(*an.fn);
    multi = new MultiscaleRVE(*an.multi);
    // need to keep the same rve in the MultiscaleRVE and in the
    // FiberRVEAnalysis
    rve = multi->getRVE();
    // you must recompute the boundary nodes to get the correct mesh entities
    auto bgn = amsi::apfMeshIterator(fn->getNetworkMesh(), 0);
    decltype(bgn) end = amsi::apfEndIterator(fn->getNetworkMesh());
    for (int sd = RVE::side::bot; sd <= RVE::side::all; ++sd) {
      getBoundaryVerts(this->rve, this->fn->getNetworkMesh(), bgn, end,
                       static_cast<RVE::side>(sd),
                       std::back_inserter(this->bnd_nds[sd]));
    }
    // We share the vecs across copies. k,u,f are zeroed in the
    // FiberRVEIteration, so in current usage this is OK. The solver is just a
    // set of buffers and a tolerance, so we should ok by copying the pointer
    vecs = an.vecs;
    es = createMicroElementalSystem(fn, getK(), getF());
    dx_fn_dx_rve = an.dx_fn_dx_rve;
    dx_fn_dx_rve_set = an.dx_fn_dx_rve_set;
    solver_eps = an.solver_eps;
    prev_itr_factor = an.prev_itr_factor;
    max_cut_attempt = an.max_cut_attempt;
    attempt_cut_factor = an.attempt_cut_factor;
    max_itrs = an.max_itrs;
    detect_osc_type = an.detect_osc_type;
  }
  // TODO determine rve size from input
  FiberRVEAnalysis::FiberRVEAnalysis(FiberNetwork *fn,
                                     FiberRVEAnalysisVecs *vecs,
                                     micro_fo_solver &slvr,
                                     micro_fo_int_solver &slvr_int)
      : fn(fn)
      , rve(new RVE(0.5, fn->getDim()))
      , vecs(vecs)
      , solver_eps(slvr.data[MICRO_SOLVER_EPS])
      , prev_itr_factor(slvr.data[PREV_ITER_FACTOR])
      , max_cut_attempt(slvr_int.data[MAX_MICRO_CUT_ATTEMPT])
      , attempt_cut_factor(slvr_int.data[MICRO_ATTEMPT_CUT_FACTOR])
      , max_itrs(slvr_int.data[MAX_MICRO_ITERS])
      , detect_osc_type(static_cast<amsi::DetectOscillationType>(
            slvr_int.data[DETECT_OSCILLATION_TYPE]))
  {
    auto bgn = amsi::apfMeshIterator(fn->getNetworkMesh(),0);
    decltype(bgn) end = amsi::apfEndIterator(fn->getNetworkMesh());
    for(int sd = RVE::side::bot; sd <= RVE::side::all; ++sd)
    {
      getBoundaryVerts(this->rve,
                       this->fn->getNetworkMesh(),
                       bgn,
                       end,
                       static_cast<RVE::side>(sd),
                       std::back_inserter(this->bnd_nds[sd]));
    }
    apf::Numbering* udofs = fn->getUNumbering(); 
    apf::NaiveOrder(udofs);
    this->es = createMicroElementalSystem(fn,getK(),getF());
    this->dx_fn_dx_rve_set = false;
  }
  FiberRVEAnalysis::~FiberRVEAnalysis() {

    delete fn;
    fn = NULL;
    delete rve;
    rve = NULL;
    delete es;
    es = NULL;
    // don't delete vecs because they are manage externally
    vecs = NULL;
    delete multi;
    multi = NULL;
    for(int i=0; i<RVE::side::all+1; ++i) {
      bnd_nds[i].clear();
    }
  }
  FiberRVEAnalysis *createFiberRVEAnalysis(FiberNetwork *fn,
                                         FiberRVEAnalysisVecs *vecs,
                                         micro_fo_solver &slvr,
                                         micro_fo_int_solver &slvr_int)
  {
    FiberRVEAnalysis *an = new FiberRVEAnalysis(fn, vecs, slvr, slvr_int);
    return an;
  }
  FiberRVEAnalysis *initFromMultiscale(
      FiberNetwork *fn, FiberRVEAnalysisVecs* vecs,
      micro_fo_header &hdr, micro_fo_params &prm, micro_fo_init_data &ini,
      micro_fo_solver &slvr, micro_fo_int_solver &slvr_int)
  {
    FiberRVEAnalysis *rve = createFiberRVEAnalysis(fn, vecs, slvr, slvr_int);
    rve->multi = new MultiscaleRVE(rve->rve, fn, hdr, prm, ini);
    return rve;
  }
  void destroyAnalysis(FiberRVEAnalysis *fa)
  {
    delete fa;
    fa = NULL;
  }
  FiberRVEAnalysisVecs *createFiberRVEAnalysisVecs(int ndofs,
                                                   las::Sparsity *csr,
                                                   las::SparskitBuffers *bfrs)
  {
    return new FiberRVEAnalysisVecs(ndofs, csr, bfrs);
  }
  void destroyFiberRVEAnalysisVecs(FiberRVEAnalysisVecs *vecs)
  {
    delete vecs;
    vecs = NULL;
  }
  FiberRVEAnalysis *copyAnalysis(FiberRVEAnalysis *an)
  {
    return new FiberRVEAnalysis(*an);
  }
  FiberRVEIteration::FiberRVEIteration(FiberRVEAnalysis *a)
      : amsi::Iteration(), an(a)
  {
  }
  void FiberRVEIteration::iterate()
  {
    auto ops = las::getLASOps<las::sparskit>();
    ops->zero(an->getK());
    ops->zero(an->getU());
    ops->zero(an->getF());
    apf::Mesh * fn = an->getFn()->getNetworkMesh();
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
               an->getFn()->getUNumbering(),
               an->getK(),
               an->getF());
    an->getSlv()->solve(an->getK(),an->getU(),an->getF());
    amsi::WriteOp wrt;
    amsi::AccumOp acm;
    amsi::FreeApplyOp fr_wrt(an->getFn()->getUNumbering(),&wrt);
    amsi::FreeApplyOp fr_acm(an->getFn()->getUNumbering(),&acm);
    double * s = NULL;
    ops->get(an->getU(),s);
    amsi::ApplyVector(an->getFn()->getUNumbering(),
                      an->getFn()->getdUField(),
                      s,0,&fr_wrt).run();
    amsi::ApplyVector(an->getFn()->getUNumbering(),
                      an->getFn()->getUField(),
                      s,0,&fr_acm).run();
    ops->restore(an->getU(),s);
    BIO_V2(int rank = -1; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
           std::cout << "RVE: " << an->getFn()->getRVEType() << " Iter: "
                     << this->iteration() << " Rank: " << rank << "\n";)
    Iteration::iterate();
  }
  void calcStress(FiberRVEAnalysis * fra, apf::Matrix3x3 & sigma)
  {
    auto ops = las::getLASOps<las::sparskit>();
    for(int ii = 0; ii < 3; ++ii)
      for(int jj = 0; jj < 3; ++jj)
        sigma[ii][jj] = 0.0;
    double * f = NULL;
    ops->get(fra->getF(),f);
    auto & bnd = fra->bnd_nds[RVE::all];
    apf::Numbering * nm = fra->getFn()->getUNumbering();
    apf::Field * uf = fra->getFn()->getUField();
    apf::Mesh * fn = fra->getFn()->getNetworkMesh();
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
    ops->restore(fra->getF(),f);
    // this is just the symmetric part of the matrix, there should be a standalone operation for this...
    sigma[0][1] = sigma[1][0] = 0.5 * (sigma[0][1] + sigma[1][0]);
    sigma[0][2] = sigma[2][0] = 0.5 * (sigma[0][2] + sigma[2][0]);
    sigma[1][2] = sigma[2][1] = 0.5 * (sigma[1][2] + sigma[2][1]);
  }
};
