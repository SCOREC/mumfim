#include "FiberRVEAnalysis.h"
#include <PCU.h>
#include <apfMeshIterator.h>
#include <las.h>
#include <lasConfig.h>
#include <lasSys.h>
#include "Fiber.h"
#include "FiberNetwork.h"
#include "MicroFOParams.h"
#include "MultiscaleMicroFOParams.h"
#include "Utility.h"
#include <apfMatrixUtil.h>
#include <lasCSRCore.h>
#include <lasSparskitExterns.h>
namespace mumfim
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
                                 fn->getXpUField(), fn->getFiberReactions(),
                                 k, f, 1);
    }
    return es;
  }
  template <>
  LinearStructs<las::sparskit>::LinearStructs(FiberNetwork * fiber_network,
                                              int ndofs,
                                              double /*solver_tol*/,
                                              std::shared_ptr<void> bfrs)
      : mSparsity(las::createCSR(fiber_network->getUNumbering(), ndofs),&las::destroySparsity<las::MICRO_BACKEND>)
  {
    if (bfrs == nullptr)
    {
      mBuffers = std::shared_ptr<void>{
        new las::SparskitBuffers(ndofs)};
    }
    else
    {
      mBuffers = bfrs;
    }
    auto buffer = static_cast<las::SparskitBuffers*>(mBuffers.get());

    las::LasCreateVec * vb = las::getVecBuilder<las::MICRO_BACKEND>(0);
    las::LasCreateMat * mb = las::getMatBuilder<las::MICRO_BACKEND>(0);
    this->k = mb->create(ndofs, LAS_IGNORE, mSparsity.get(), MPI_COMM_SELF);
    this->f = vb->createRHS(this->k);
    this->u = vb->createLHS(this->k);
    buffer->resizeMatrixBuffer(ndofs * 100);
    this->slv = std::unique_ptr<las::Solve>(las::createSparskitLUSolve(buffer, 1e-6));
  }
  /*
  template <>
  LinearStructs<las::petsc>::LinearStructs(FiberNetwork* fiber_network, int ndofs,
                               double solver_tol,
                               std::shared_ptr<void> bfrs)
                              : mSparsity(las::createPetscSparsity(fiber_network->getUNumbering(), ndofs, PCU_Get_Comm()),
                                  &las::destroySparsity<las::MICRO_BACKEND>)
  {
    las::LasCreateVec * vb = las::getVecBuilder<las::MICRO_BACKEND>(0);
    las::LasCreateMat * mb = las::getMatBuilder<las::MICRO_BACKEND>(0);
    this->k = mb->create(ndofs, 1, mSparsity.get(), MPI_COMM_SELF);
    this->f = vb->createRHS(this->k);
    this->u = vb->createLHS(this->k);
    // FIXME memory leak (won't be hit in multi-scale)
    // because we provide buffers
    this->slv = std::unique_ptr<las::Solve>(las::createPetscLUSolve(MPI_COMM_SELF));
  }
  */
  template <typename T>
  LinearStructs<T>::~LinearStructs()
  {
    auto md = las::getMatBuilder<las::MICRO_BACKEND>(0);
    auto vd = las::getVecBuilder<las::MICRO_BACKEND>(0);
    md->destroy(k);
    vd->destroy(u);
    vd->destroy(f);
  }
  // specifically instantiate linear structs for our
  // solver backends
  template class LinearStructs<las::sparskit>;
  template class LinearStructs<las::petsc>;

  FiberRVEAnalysis::FiberRVEAnalysis(const FiberRVEAnalysis & an)
    : RVEAnalysis(an)
  {
    mFiberNetwork = std::unique_ptr<FiberNetwork>(new FiberNetwork(*an.mFiberNetwork));
    rve = std::unique_ptr<RVE>(new RVE(*an.rve));
    // you must recompute the boundary nodes to get the correct mesh entities
    auto bgn = amsi::apfMeshIterator(mFiberNetwork->getNetworkMesh(), 0);
    decltype(bgn) end = amsi::apfEndIterator(mFiberNetwork->getNetworkMesh());
    for (int sd = RVE::side::bot; sd <= RVE::side::all; ++sd)
    {
      getBoundaryVerts(this->rve.get(),
                       this->mFiberNetwork->getNetworkMesh(),
                       bgn,
                       end,
                       static_cast<RVE::side>(sd),
                       std::back_inserter(this->bnd_nds[sd]));
    }
    // We share the vecs across copies. k,u,f are zeroed in the
    // FiberRVEIteration, so in current usage this is OK. The solver is just a
    // set of buffers and a tolerance, so we should ok by copying the pointer
    vecs = an.vecs;
    solver_eps = an.solver_eps;
    prev_itr_factor = an.prev_itr_factor;
    max_cut_attempt = an.max_cut_attempt;
    attempt_cut_factor = an.attempt_cut_factor;
    max_itrs = an.max_itrs;
    detect_osc_type = an.detect_osc_type;
  }
  // TODO determine rve size from input
  FiberRVEAnalysis::FiberRVEAnalysis(std::unique_ptr<FiberNetwork> fn,
                                     std::unique_ptr<MicroSolutionStrategy> ss,
                                     std::shared_ptr<void> workspace)
      : mFiberNetwork(std::move(fn)), mSolutionStrategy(std::move(ss))
      , rve(new RVE(0.5, mFiberNetwork->getDim()))
      , solver_eps(mSolutionStrategy->cnvgTolerance)
      , prev_itr_factor(mSolutionStrategy->oscPrms.prevNormFactor)
      , max_cut_attempt(mSolutionStrategy->oscPrms.maxMicroCutAttempts)
      , attempt_cut_factor(mSolutionStrategy->oscPrms.microAttemptCutFactor)
      , max_itrs(mSolutionStrategy->oscPrms.maxIterations)
      , detect_osc_type(mSolutionStrategy->oscPrms.oscType)
  {
    auto bgn = amsi::apfMeshIterator(mFiberNetwork->getNetworkMesh(), 0);
    decltype(bgn) end = amsi::apfEndIterator(mFiberNetwork->getNetworkMesh());
    for (int sd = RVE::side::bot; sd <= RVE::side::all; ++sd)
    {
      getBoundaryVerts(this->rve.get(),
                       this->mFiberNetwork->getNetworkMesh(),
                       bgn,
                       end,
                       static_cast<RVE::side>(sd),
                       std::back_inserter(this->bnd_nds[sd]));
    }
    auto num_dofs = mFiberNetwork->getDofCount();
    //vecs = std::make_shared<LinearStructs<las::MICRO_BACKEND>>(createLinearStructs(mFiberNetwork.get(), num_dofs,solver_eps, workspace));
    vecs = createLinearStructs(getFn(), num_dofs, solver_eps, workspace);
  }
  FiberRVEAnalysis::~FiberRVEAnalysis()
  {
    //for (int i = 0; i < RVE::side::all + 1; ++i)
    //{
    //  bnd_nds[i].clear();
    //}
  }
  FiberRVEAnalysis& FiberRVEAnalysis::operator=(FiberRVEAnalysis&& other)
  {
    RVEAnalysis::operator=(other);
    mFiberNetwork = std::move(other.mFiberNetwork);
    other.mFiberNetwork = nullptr;
    rve = std::move(other.rve);
    for (int sd = RVE::side::bot; sd <= RVE::side::all; ++sd)
    {
      bnd_nds[sd] = std::move(other.bnd_nds[sd]);
    }
    // We share the vecs across copies. k,u,f are zeroed in the
    // FiberRVEIteration, so in current usage this is OK. The solver is just a
    // set of buffers and a tolerance, so we should ok by copying the pointer
    vecs = std::move(other.vecs);
    solver_eps = other.solver_eps;
    prev_itr_factor = other.prev_itr_factor;
    max_cut_attempt = other.max_cut_attempt;
    attempt_cut_factor = other.attempt_cut_factor;
    max_itrs = other.max_itrs;
    detect_osc_type = other.detect_osc_type;
    return *this; 
  }
  // TODO This can be turned into a kokkos kernel
  void FiberRVEAnalysis::computeCauchyStress(double sigma[6])
  {
    auto ops = las::getLASOps<las::MICRO_BACKEND>();
    int dim = this->getFn()->getNetworkMesh()->getDimension();
    apf::Field * xpu = this->getFn()->getXpUField();
    apf::Vector3 xpu_val;
    apf::Vector3 f_val;
    apf::Numbering * num = this->getFn()->getUNumbering();
    apf::Vector3 coord, U;
    apf::Field* coord_field = this->getFn()->getNetworkMesh()->getCoordinateField();
    apf::Field* u_field = this->getFn()->getUField();
    int dofs[3];
    double * f = NULL;
    ops->get(this->getF(), f);
    apf::Matrix3x3 strs(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for (auto nd = this->bnd_nds[RVE::all].begin();
         nd != this->bnd_nds[RVE::all].end();
         ++nd)
    {
      apf::getVector(u_field, *nd, 0, U);
      apf::getVector(coord_field, *nd, 0, coord);
      apf::getVector(xpu, *nd, 0, xpu_val);
      for (int ii = 0; ii < this->getFn()->getDim(); ++ii)
        dofs[ii] = apf::getNumber(num, *nd, 0, ii);
      strs[0][0] += xpu_val[0] * f[dofs[0]];
      strs[0][1] += xpu_val[1] * f[dofs[0]];
      strs[0][2] += xpu_val[2] * f[dofs[0]];
      strs[1][0] += xpu_val[0] * f[dofs[1]];
      strs[1][1] += xpu_val[1] * f[dofs[1]];
      strs[1][2] += xpu_val[2] * f[dofs[1]];
      strs[2][0] += xpu_val[0] * f[dofs[2]];
      strs[2][1] += xpu_val[1] * f[dofs[2]];
      strs[2][2] += xpu_val[2] * f[dofs[2]];
    }
    // take the symmetric part of the current stress matrix
    apf::Matrix3x3 sym_strs = amsi::symmetricPart(strs);
    amsi::mat2VoigtVec(dim, sym_strs, &sigma[0]);
  }
  std::unique_ptr<LinearStructs<las::MICRO_BACKEND>> createLinearStructs(FiberNetwork* fiber_network, int ndofs,double solver_tol,
                                      std::shared_ptr<void> bfrs)
  {
    return std::unique_ptr<LinearStructs<las::MICRO_BACKEND>>(new LinearStructs<las::MICRO_BACKEND>(fiber_network, ndofs, solver_tol, bfrs));
  }
  // since we no longer compute the corner derivatives we need to thing about how
  // to do first order continuation...
  void applyGuessSolution(FiberRVEAnalysis * ans, const DeformationGradient & dfmGrd)
  {
    int d = ans->getRVE()->getDim();
    assert(d == 3 || d == 2);
    apf::Matrix3x3 F;
    for (int ei = 0; ei < d; ++ei)
      for (int ej = 0; ej < d; ++ej)
        F[ei][ej] = dfmGrd[ei * d + ej];
    apf::Mesh * rve_msh = ans->getRVE()->getMesh();
    apf::Field * rve_du = ans->getRVE()->getdUField();
    apf::Field * rve_u = ans->getRVE()->getUField();
    // apply deformation gradient to the rve cube
    ApplyIncrementalDeformationGradient(F, rve_msh, rve_du, rve_u).run();
    apf::Mesh * fn_msh = ans->getFn()->getNetworkMesh();
    apf::Field * fn_du = ans->getFn()->getdUField();
    apf::Field * fn_u = ans->getFn()->getUField();
    // apply deformation gradient to the fiber network mesh
    //ApplyDeformationGradient(F, fn_msh, fn_du, fn_u).run();
    ApplyIncrementalDeformationGradient(F, fn_msh, fn_du, fn_u).run();
  }
};  // namespace mumfim
