#include "bioFiberRVEAnalysis.h"
#include <PCU.h>
#include <apfMeshIterator.h>
#include <bioVerbosity.h>
#include <las.h>
#include <lasConfig.h>
#include <lasSys.h>
#include "bioFiber.h"
#include "bioFiberNetwork.h"
#include "bioFiberNetworkIO.h"
#include "bioMultiscaleCoupling.h"
#include "bioVerbosity.h"
#include "bioMicroFOParams.h"
namespace bio
{
  // todo: rename (shouldn't have reference to micro in a single-scale file)
  static apf::Integrator * createMicroElementalSystem(FiberNetwork * fn,
                                                      las::Mat * k,
                                                      las::Vec * f)
  {
    apf::Integrator * es = NULL;
    FiberMember tp = fn->getFiberMember();
    if (tp == FiberMember::truss)
      es = new TrussIntegrator(fn->getUNumbering(),
                               fn->getUField(),
                               fn->getXpUField(),
                               &(fn->getFiberReactions()[0]),
                               k,
                               f,
                               1);
    return es;
  }
  LinearStructs::LinearStructs(int ndofs,
                                             las::Sparsity * csr,
                                             las::SparskitBuffers * b)
  {
    las::LasCreateVec * vb = las::getVecBuilder<las::sparskit>(0);
    las::LasCreateMat * mb = las::getMatBuilder<las::sparskit>(0);
    this->f = vb->create(ndofs, LAS_IGNORE, MPI_COMM_SELF);
    this->u = vb->create(ndofs, LAS_IGNORE, MPI_COMM_SELF);
    this->k = mb->create(ndofs, LAS_IGNORE, csr, MPI_COMM_SELF);
    // FIXME memory leak (won't be hit in multi-scale)
    // because we provide buffers
    if (b == NULL) b = new las::SparskitBuffers(ndofs);
    this->slv = las::createSparskitLUSolve(b, 1e-6);
  }
  LinearStructs::~LinearStructs()
  {
    auto md = las::getMatBuilder<las::sparskit>(0);
    auto vd = las::getVecBuilder<las::sparskit>(0);
    md->destroy(k);
    vd->destroy(u);
    vd->destroy(f);
    delete slv;
  }
  FiberRVEAnalysis::FiberRVEAnalysis(const FiberRVEAnalysis & an)
  {
    fn = new FiberNetwork(*an.fn);
    rve = new RVE(*an.rve);
    // note if we have a singlescale run, the MultiscaleRVE pointer will be null
    if (an.multi)
    {
      multi = new MultiscaleRVE(*an.multi);
      multi->setRVE(rve);
    }
    else
    {
      multi = NULL;
    }
    // you must recompute the boundary nodes to get the correct mesh entities
    auto bgn = amsi::apfMeshIterator(fn->getNetworkMesh(), 0);
    decltype(bgn) end = amsi::apfEndIterator(fn->getNetworkMesh());
    for (int sd = RVE::side::bot; sd <= RVE::side::all; ++sd)
    {
      getBoundaryVerts(this->rve,
                       this->fn->getNetworkMesh(),
                       bgn,
                       end,
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
  FiberRVEAnalysis::FiberRVEAnalysis(FiberNetwork * fn,
                                     LinearStructs * vecs,
                                     micro_fo_solver & slvr,
                                     micro_fo_int_solver & slvr_int)
      : fn(fn)
      , multi(NULL)
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
    auto bgn = amsi::apfMeshIterator(fn->getNetworkMesh(), 0);
    decltype(bgn) end = amsi::apfEndIterator(fn->getNetworkMesh());
    for (int sd = RVE::side::bot; sd <= RVE::side::all; ++sd)
    {
      getBoundaryVerts(this->rve,
                       this->fn->getNetworkMesh(),
                       bgn,
                       end,
                       static_cast<RVE::side>(sd),
                       std::back_inserter(this->bnd_nds[sd]));
    }
    apf::Numbering * udofs = fn->getUNumbering();
    apf::NaiveOrder(udofs);
    this->es = createMicroElementalSystem(fn, getK(), getF());
    this->dx_fn_dx_rve_set = false;
  }
  FiberRVEAnalysis::~FiberRVEAnalysis()
  {
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
    for (int i = 0; i < RVE::side::all + 1; ++i)
    {
      bnd_nds[i].clear();
    }
  }
  struct val_gen
  {
    val_gen(FiberRVEAnalysis * a) : an(a), prv_nrm(1.0) {}
    FiberRVEAnalysis * an;
    double prv_nrm;
    double operator()()
    {
      auto ops = las::getLASOps<las::sparskit>();
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
  bool FiberRVEAnalysis::run(const DeformationGradient & dfmGrd)
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
      FiberRVEAnalysis * tmpRVE = copyAnalysis(this);
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
           microAttemptIter <= attemptCutFactor;
           ++microAttemptIter)
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
        applyDeformation(tmpRVE, appliedDefm);
        FiberRVEIteration rveItr(tmpRVE);
        std::vector<amsi::Iteration *> itr_stps = {&rveItr};
        amsi::MultiIteration itr(itr_stps.begin(), itr_stps.end());
        amsi::UpdatingConvergence<decltype(&vg), decltype(&eg), decltype(&rg)>
            resid_cnvrg(&itr, &vg, &eg, &rg);
        std::vector<amsi::Convergence *> cnvrg_stps = {&resid_cnvrg};
        amsi::MultiConvergence cnvrg(cnvrg_stps.begin(), cnvrg_stps.end());
        amsi::Iteration * osc_itr =
            amsi::createOscillationDetection<decltype(&resid_cnvrg)>(
                tmpRVE->detect_osc_type,
                &resid_cnvrg,
                &rveItr,
                tmpRVE->max_itrs,
                tmpRVE->prev_itr_factor);
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
  FiberRVEAnalysis * createFiberRVEAnalysis(FiberNetwork * fn,
                                            LinearStructs * vecs,
                                            micro_fo_solver & slvr,
                                            micro_fo_int_solver & slvr_int)
  {
    FiberRVEAnalysis * an = new FiberRVEAnalysis(fn, vecs, slvr, slvr_int);
    return an;
  }
  FiberRVEAnalysis * initFromMultiscale(FiberNetwork * fn,
                                        LinearStructs * vecs,
                                        micro_fo_header & hdr,
                                        micro_fo_params & prm,
                                        micro_fo_init_data & ini,
                                        micro_fo_solver & slvr,
                                        micro_fo_int_solver & slvr_int)
  {
    FiberRVEAnalysis * rve = createFiberRVEAnalysis(fn, vecs, slvr, slvr_int);
    rve->multi = new MultiscaleRVE(rve->rve, fn, hdr, prm, ini);
    return rve;
  }
  void destroyAnalysis(FiberRVEAnalysis * fa)
  {
    delete fa;
    fa = NULL;
  }
  LinearStructs * createLinearStructs(int ndofs,
                                                    las::Sparsity * csr,
                                                    las::SparskitBuffers * bfrs)
  {
    return new LinearStructs(ndofs, csr, bfrs);
  }
  void destroyLinearStructs(LinearStructs * vecs)
  {
    delete vecs;
    vecs = NULL;
  }
  FiberRVEAnalysis * copyAnalysis(FiberRVEAnalysis * an)
  {
    return new FiberRVEAnalysis(*an);
  }
  FiberRVEIteration::FiberRVEIteration(FiberRVEAnalysis * a)
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
    while ((me = fn->iterate(itr)))
    {
      apf::MeshElement * mlm = apf::createMeshElement(fn, me);
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
  void calcStress(FiberRVEAnalysis * fra, apf::Matrix3x3 & sigma)
  {
    auto ops = las::getLASOps<las::sparskit>();
    for (int ii = 0; ii < 3; ++ii)
      for (int jj = 0; jj < 3; ++jj)
        sigma[ii][jj] = 0.0;
    double * f = NULL;
    ops->get(fra->getF(), f);
    auto & bnd = fra->bnd_nds[RVE::all];
    apf::Numbering * nm = fra->getFn()->getUNumbering();
    apf::Field * uf = fra->getFn()->getUField();
    apf::Mesh * fn = fra->getFn()->getNetworkMesh();
    for (auto vrt = bnd.begin(); vrt != bnd.end(); ++vrt)
    {
      int dof[3] = {};
      dof[0] = apf::getNumber(nm, *vrt, 0, 0);
      dof[1] = apf::getNumber(nm, *vrt, 0, 1);
      dof[2] = apf::getNumber(nm, *vrt, 0, 2);
      apf::Vector3 x;
      fn->getPoint(*vrt, 0, x);
      apf::Vector3 u;
      apf::getVector(uf, *vrt, 0, u);
      x += u;
      for (int ii = 0; ii < 3; ++ii)
        for (int jj = 0; jj < 3; ++jj)
          sigma[ii][jj] += f[dof[ii]] * x[jj];
    }
    ops->restore(fra->getF(), f);
    // this is just the symmetric part of the matrix, there should be a
    // standalone operation for this...
    sigma[0][1] = sigma[1][0] = 0.5 * (sigma[0][1] + sigma[1][0]);
    sigma[0][2] = sigma[2][0] = 0.5 * (sigma[0][2] + sigma[2][0]);
    sigma[1][2] = sigma[2][1] = 0.5 * (sigma[1][2] + sigma[2][1]);
  }
  void applyDeformation(FiberRVEAnalysis * ans, const DeformationGradient & dfmGrd)
  {
    int d = ans->rve->getDim();
    assert(d == 3 || d == 2);
    apf::Matrix3x3 F;
    for (int ei = 0; ei < d; ++ei)
      for (int ej = 0; ej < d; ++ej)
        F[ei][ej] = dfmGrd[ei * d + ej];
    apf::Mesh * rve_msh = ans->rve->getMesh();
    apf::Field * rve_du = ans->rve->getdUField();
    apf::Field * rve_u = ans->rve->getUField();
    // apply deformation gradient to the rve cube
    ApplyDeformationGradient(F, rve_msh, rve_du, rve_u).run();
    apf::Mesh * fn_msh = ans->getFn()->getNetworkMesh();
    apf::Field * fn_du = ans->getFn()->getdUField();
    apf::Field * fn_u = ans->getFn()->getUField();
    // only apply first order continuation if we have derivative information
    // from the last iteration
    if (ans->dx_fn_dx_rve_set)
    {
      int dim = ans->rve->getDim();
      int nnd = ans->rve->numNodes();
      apf::DynamicVector rve_duv(nnd * dim);
      // array of the displacement field change since last itr of the rve cube
      apf::NewArray<apf::Vector3> rve_dus;
      apf::Element * du_elmt =
          apf::createElement(rve_du, ans->rve->getMeshEnt());
      apf::getVectorNodes(du_elmt, rve_dus);
      apf::destroyElement(du_elmt);
      apf::NewArray<int> rve_dofs;
      apf::getElementNumbers(
          ans->rve->getNumbering(), ans->rve->getMeshEnt(), rve_dofs);
      // reorder rve dofs with dofids since dx_fn_dx_rve uses them
      for (int nd = 0; nd < nnd; ++nd)
      {
        for (int dm = 0; dm < dim; ++dm)
        {
          rve_duv[rve_dofs[nd * dim + dm]] = rve_dus[nd][dm];
        }
      }
      int fn_dofs = ans->getFn()->getDofCount();
      apf::DynamicVector du_fn(fn_dofs);
      apf::multiply(ans->dx_fn_dx_rve, rve_duv, du_fn);
      ApplyDeformationGradient app_F(F, fn_msh, fn_du, fn_u);
      // Apply the deformation gradient to the boundary nodes and do not add
      // extra term from first order continuation e.g. set du_fn to zero on
      // boundary
      for (auto vrt = ans->bnd_nds[RVE::side::all].begin();
           vrt != ans->bnd_nds[RVE::side::all].end();
           ++vrt)
      {
        for (int dm = 0; dm < dim; ++dm)
        {
          int dof = apf::getNumber(ans->getFn()->getUNumbering(), *vrt, 0, dm);
          du_fn(dof) = 0.0;
        }
        app_F.inEntity(*vrt);
        app_F.atNode(0);
        app_F.outEntity();
      }
      // note that we want to accumulate the du*fx_fn_dx_rve onto the u field
      ApplySolution(fn_u, ans->getFn()->getUNumbering(), &du_fn[0], 0, true)
          .run();
    }
    else
    {
      // apply deformation gradient to the fiber network mesh
      ApplyDeformationGradient(F, fn_msh, fn_du, fn_u).run();
    }
  }
};  // namespace bio
