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
  apf::Integrator * createImplicitMicroElementalSystem(FiberNetwork * fn,
                                               las::Mat * k,
                                               las::Vec * f)
  {
    apf::Integrator * es = NULL;
    FiberMember tp = fn->getFiberMember();
    if (tp == FiberMember::truss) {
        es = new TrussIntegrator(fn->getUNumbering(), fn->getUField(),
                                 fn->getXpUField(), &(fn->getFiberReactions()[0]),
                                 k, f, 1);
    }
    return es;
  }
  template <>
  LinearStructs<las::sparskit>::LinearStructs(int ndofs,
                               double solver_tol,
                               las::Sparsity * csr,
                               void * bfrs,
                               las::Sparsity * massCsr)
  {
    if(massCsr == NULL)
      massCsr = csr;
    las::SparskitBuffers * b = static_cast<las::SparskitBuffers *>(bfrs);
    las::LasCreateVec * vb = las::getVecBuilder<las::MICRO_BACKEND>(0);
    las::LasCreateMat * mb = las::getMatBuilder<las::MICRO_BACKEND>(0);
    this->k = mb->create(ndofs, LAS_IGNORE, csr, MPI_COMM_SELF);
    this->f = vb->createRHS(this->k);
    this->u = vb->createLHS(this->k);
    // FIXME memory leak (won't be hit in multi-scale)
    // because we provide buffers
    if (b == NULL) b = new las::SparskitBuffers(ndofs);
    this->slv = las::createSparskitLUSolve(b, 1e-6);
  }
  template <>
  LinearStructs<las::petsc>::LinearStructs(int ndofs,
                               double solver_tol,
                               las::Sparsity * sprs,
                               void * bfrs,
                               las::Sparsity * massSprs)
  {
    if(massSprs == NULL)
      massSprs = sprs;
    las::LasCreateVec * vb = las::getVecBuilder<las::MICRO_BACKEND>(0);
    las::LasCreateMat * mb = las::getMatBuilder<las::MICRO_BACKEND>(0);
    this->k = mb->create(ndofs, 1, sprs, MPI_COMM_SELF);
    this->f = vb->createRHS(this->k);
    this->u = vb->createLHS(this->k);
    // FIXME memory leak (won't be hit in multi-scale)
    // because we provide buffers
    this->slv = las::createPetscLUSolve(MPI_COMM_SELF);
  }
  template <typename T>
  LinearStructs<T>::~LinearStructs()
  {
    auto md = las::getMatBuilder<las::MICRO_BACKEND>(0);
    auto vd = las::getVecBuilder<las::MICRO_BACKEND>(0);
    md->destroy(k);
    vd->destroy(u);
    vd->destroy(f);
    delete slv;
  }
  // specifically instantiate linear structs for our
  // solver backends
  template class LinearStructs<las::sparskit>;
  template class LinearStructs<las::petsc>;

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
                                     LinearStructs<las::MICRO_BACKEND> * vecs,
                                     const MicroSolutionStrategy & ss)
      : fn(fn)
      , multi(NULL)
      , rve(new RVE(0.5, fn->getDim()))
      , vecs(vecs)
      , solver_eps(ss.cnvgTolerance)
      , prev_itr_factor(ss.oscPrms.prevNormFactor)
      , max_cut_attempt(ss.oscPrms.maxMicroCutAttempts)
      , attempt_cut_factor(ss.oscPrms.microAttemptCutFactor)
      , max_itrs(ss.oscPrms.maxIterations)
      , detect_osc_type(ss.oscPrms.oscType)
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
  FiberRVEAnalysis * createFiberRVEAnalysis(FiberNetwork * fn,
                                            LinearStructs<las::MICRO_BACKEND> * vecs,
                                            micro_fo_solver & slvr,
                                            micro_fo_int_solver & slvr_int)
  {
    // TODO this somewhat inefficient, and will be fixed when we directly
    // communicate the solution strategy struct
    MicroSolutionStrategyExplicit ss;
    ss.cnvgTolerance = slvr.data[MICRO_CONVERGENCE_TOL];
    ss.slvrTolerance = slvr.data[MICRO_SOLVER_TOL];
    ss.total_time = slvr.data[LOAD_TIME]+slvr.data[HOLD_TIME];
    ss.load_time = slvr.data[LOAD_TIME];
    ss.crit_time_scale_factor = slvr.data[CRITICAL_TIME_SCALE_FACTOR];
    ss.visc_damp_coeff = slvr.data[VISCOUS_DAMPING_FACTOR];
    ss.energy_check_eps = slvr.data[ENERGY_CHECK_EPSILON];
    ss.slvrType = static_cast<SolverType>(slvr_int.data[MICRO_SOLVER_TYPE]);
    ss.ampType = static_cast<AmplitudeType>(slvr_int.data[AMPLITUDE_TYPE]);
    ss.print_history_frequency = slvr_int.data[PRINT_HISTORY_FREQUENCY];
    ss.print_field_frequency = slvr_int.data[PRINT_FIELD_FREQUENCY];
    ss.print_field_by_num_frames = slvr_int.data[PRINT_FIELD_BY_NUM_FRAMES];
    ss.serial_gpu_cutoff = slvr_int.data[SERIAL_GPU_CUTOFF];
    ss.oscPrms.maxIterations = slvr_int.data[MAX_MICRO_ITERS];
    ss.oscPrms.maxMicroCutAttempts = slvr_int.data[MAX_MICRO_CUT_ATTEMPT];
    ss.oscPrms.microAttemptCutFactor = slvr_int.data[MICRO_ATTEMPT_CUT_FACTOR];
    ss.oscPrms.oscType = static_cast<amsi::DetectOscillationType>(
        slvr_int.data[DETECT_OSCILLATION_TYPE]);
    ss.oscPrms.prevNormFactor = slvr.data[PREV_ITER_FACTOR];
    if(ss.slvrType == SolverType::Implicit)
      return createFiberRVEAnalysis(fn, vecs, ss, FiberRVEAnalysisType::StaticImplicit);
    else if(ss.slvrType == SolverType::Explicit)
      return createFiberRVEAnalysis(fn, vecs, ss, FiberRVEAnalysisType::Explicit);
    else
    {
      std::cerr<<"Incorrect solver type selected when initializing rve analysis"<<std::endl;
      abort();
    }
  }
  FiberRVEAnalysis * initFromMultiscale(FiberNetwork * fn,
                                        LinearStructs<las::MICRO_BACKEND> * vecs,
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
  LinearStructs<las::MICRO_BACKEND> * createLinearStructs(int ndofs,double solver_tol,
                                      las::Sparsity * csr,
                                      void * bfrs, las::Sparsity * massCsr
                                      )
  {
    return new LinearStructs<las::MICRO_BACKEND>(ndofs, solver_tol, csr, bfrs, massCsr);
  }
  void destroyLinearStructs(LinearStructs<las::MICRO_BACKEND> * vecs)
  {
    delete vecs;
    vecs = NULL;
  }
  /*
  FiberRVEAnalysis * copyAnalysis(FiberRVEAnalysis * an)
  {
    return new FiberRVEAnalysis(*an);
  }
  */
  // TODO: Change this to use the frozen field from the mesh database
  void calcStress(FiberRVEAnalysis * fra, apf::Matrix3x3 & sigma)
  {
    auto ops = las::getLASOps<las::MICRO_BACKEND>();
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
  void applyGuessSolution(FiberRVEAnalysis * ans, const DeformationGradient & dfmGrd)
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
