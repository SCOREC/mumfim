#include <pcu_util.h>
#include "bioFiberRVEAnalysis.h"
#include "bioMassIntegrator.h"
namespace bio
{
  /*
  FiberRVEAnalysisQSExplicit::FiberRVEAnalysisQSExplicit(
      const FiberRVEAnalysisSImplicit & an)
      : FiberRVEAnalysis(an)
  {
  }
  */
  // static void applyBndNdsToVec(an->getFn()->getANumbering(),
  //                             an->getFn()->getAField(),
  //                             an->getA(),
  //                             an->bnd_nds[RVE::all].begin(),
  //                             an->bnd_nds[RVE::all].end());
  template <typename I>
  static void applyBndNdsToVec(apf::Numbering * num,
                               apf::Field * field,
                               las::Vec * vec,
                               I bnd_bgn,
                               I bnd_end)
  {
    apf::Mesh * mesh = apf::getMesh(field);
    apf::FieldShape * shp = apf::getShape(field);
    int ncomp = apf::countComponents(field);
    double * cmps = new double[ncomp];
    auto ops = las::getLASOps<las::MICRO_BACKEND>();
    int dof;
    for (I ent = bnd_bgn; ent != bnd_end; ++ent)
    {
      int nnodes = shp->countNodesOn(mesh->getType(*ent));
      for (int i = 0; i < nnodes; ++i)
      {
        memset(&cmps[0], 0.0, sizeof(double) * ncomp);
        apf::getComponents(field, *ent, i, cmps);
        for (int cmp = 0; cmp < ncomp; ++cmp)
        {
          if (apf::isNumbered(num, *ent, i, cmp))
          {
            dof = apf::getNumber(num, *ent, i, cmp);
            ops->set(vec, 1, &dof, &(cmps[i]));
          }
        }
      }
    }
    delete [] cmps;
  }
  class ApplyDeformationGradientExplicit : public ApplyDeformationGradient
  {
    protected:
    apf::Matrix3x3 dVdX;
    apf::Matrix3x3 dAdX;
    apf::Field * vField;
    apf::Field * aField;
    int numComp;
    public:
    template <class IteratorType>
    ApplyDeformationGradientExplicit(IteratorType ent_bgn,
                                     IteratorType ent_end,
                                     apf::Matrix3x3 F,
                                     Amplitude * amp,
                                     double time,
                                     apf::Mesh * mesh,
                                     apf::Field * uField,
                                     apf::Field * duField,
                                     apf::Field * vField,
                                     apf::Field * aField
                                     )
        : ApplyDeformationGradient(ent_bgn, ent_end, F, mesh, duField, uField)
        , vField(vField)
        , aField(aField)
    {
      // only apply the portion of the deformation that we expect at time t
      FmI = FmI * (*amp)(time);
      dVdX = FmI * amp->derivative(time);
      dAdX = FmI * amp->secondDerivative(time);
    }
    virtual void atNode(int nd)
    {
      ApplyDeformationGradient::atNode(nd);
      // apply velocity and acceleration fields
      apf::Vector3 nd_v = dVdX * nd_xyz;
      apf::Vector3 nd_a = dAdX * nd_xyz;
      apf::setVector(vField, ent, nd, nd_v);
      apf::setVector(aField, ent, nd, nd_a);
    }
  };
  FiberRVEAnalysisQSExplicit::FiberRVEAnalysisQSExplicit(
      FiberNetwork * fn,
      LinearStructs<las::MICRO_BACKEND> * vecs,
      const MicroSolutionStrategy & ss)
      : FiberRVEAnalysis(fn, vecs, ss)
      , internal_energy(0)
      , kinetic_energy(0)
      , external_energy(0)
      , time(0)
  {
    apf::Integrator * massIntegrator = new MassIntegrator(
        fn->getUNumbering(), fn->getUField(), this->getM(),
        &(fn->getFiberReactions()[0]), 2, MassLumpType::RowSum);
    massIntegrator->process(this->getFn()->getNetworkMesh(), 1);
    delete massIntegrator;
    es = createExplicitMicroElementalSystem(fn, getK(), getF(), getFInt());
  }
  FiberRVEIterationQSExplicit::FiberRVEIterationQSExplicit(
      FiberRVEAnalysisQSExplicit * a,
      DeformationGradient appliedDefm,
      Amplitude * amplitude,
      double loadTime,
      double timeStepFactor,
      double massDampingFactor,
      double stiffnessDampingFactor,
      unsigned int printSteps)
      : amsi::Iteration()
      , an(a)
      , appliedDefm(appliedDefm)
      , amplitude(amplitude)
      , last_iter(false)
      , print_steps(printSteps)
      , load_time(loadTime)
      , time_step_factor(timeStepFactor)
      , mass_damping_factor(massDampingFactor)
      , stiffness_damping_factor(stiffnessDampingFactor)
  {
    assert(time_step_factor <= 1.0);
    auto vb = las::getVecBuilder<las::MICRO_BACKEND>(0);
    tmp = vb->createRHS(an->getC());
  }
  // see belytschko box 6.1
  void FiberRVEIterationQSExplicit::iterate()
  {
    // TODO check deformation wave speed
    auto ops = las::getLASOps<las::MICRO_BACKEND>();
    ops->zero(an->getC());  // don't bother zeroing because we reset it anyways
    las::MatMatAdd * smsma = las::getMatMatAdd<las::MICRO_BACKEND>();
    las::Mat * c;
    las::MatVecMult * mvm = las::getMatVecMult<las::MICRO_BACKEND>();
    amsi::WriteOp wrt;
    amsi::AccumOp accum;
    if (__builtin_expect(this->iteration() == 0, false))
    {
      ops->zero(an->getU());
      ops->zero(an->getV());
      ops->zero(an->getA());
      ops->zero(an->getK());
      ops->zero(an->getF());
      ops->zero(an->getFInt());
      ops->zero(an->getFExt());
      ops->zero(an->getFDamp());
      ops->zero(an->getPrevFInt());
      ops->zero(an->getPrevFExt());
      ops->zero(an->getPrevFDamp());
      an->es->process(an->getFn()->getNetworkMesh(), 1);
    }
    // ITERATIONS START HERE!
    // 4. update time t(n+1) = t(n)+delta_t(n+1/2), t(n+1/2)=1/2(t(n)+t(n+1))
    double delta_t =
        dynamic_cast<ExplicitTrussIntegrator *>(an->es)->getCriticalTimeStep() *
        time_step_factor;
    if (__builtin_expect(this->iteration() == 0, false))
    {
      std::cout << "Assuming a constant timestep (dt=" << delta_t
                << "),\nthis should take approx: "
                << static_cast<unsigned long>(load_time / delta_t)
                << " timesteps";
      std::cout << " to acheive T_max=" << load_time << std::endl;
    }
    double t_npone = an->getTime() + delta_t;
    if (__builtin_expect(t_npone >= load_time, false))
    {
      delta_t = t_npone - load_time;
      t_npone = load_time;
      last_iter = true;
    }
    double t_nphalf = 0.5 * (an->getTime() + t_npone);
    assert(t_nphalf);
    // 5. partial update of nodal velocities v(n+1/2)=v(n)+(t(n+1/2)-t(n))a(n)
    ops->axpy((t_nphalf - an->getTime()), an->getA(), an->getV());
    // 6. enforce dirichlet boundary conditions
    // FIXME we should directly apply the deformation gradient to
    // the array, so we don't wase compute time copying data around!
    apf::Matrix3x3 F(
        appliedDefm.data[0], appliedDefm.data[1], appliedDefm.data[2],
        appliedDefm.data[3], appliedDefm.data[4], appliedDefm.data[5],
        appliedDefm.data[6], appliedDefm.data[7], appliedDefm.data[8]);
    ApplyDeformationGradientExplicit(
        an->bnd_nds[RVE::all].begin(), an->bnd_nds[RVE::all].end(), F,
        amplitude, t_npone, an->getFn()->getNetworkMesh(),
        an->getFn()->getUField(), an->getFn()->getdUField(),
        an->getFn()->getVField(), an->getFn()->getAField())
        .run();
    applyBndNdsToVec(an->getFn()->getVNumbering(), an->getFn()->getVField(),
                     an->getV(), an->bnd_nds[RVE::all].begin(),
                     an->bnd_nds[RVE::all].end());
    // 7. update nodal displacements u(n+1) = u(n)+delta_t(n+1/2)*v(n+1/2)
    ops->axpy(delta_t, an->getV(), an->getU());
    // apply u boundary conditions
    applyBndNdsToVec(an->getFn()->getUNumbering(), an->getFn()->getUField(),
                     an->getU(), an->bnd_nds[RVE::all].begin(),
                     an->bnd_nds[RVE::all].end());
    // 8. getforce
    ops->zero(an->getK());
    ops->zero(an->getF());
    ops->zero(an->getFInt());
    an->es->process(an->getFn()->getNetworkMesh(), 1);
    // TODO we can make this more efficient by computing these energies in the
    // force integrator
    // Wint_n+delta_t_nphalf/2*v_nphalf^T(f_int_n+f_int_npone)
    // ops->axpy(1, an->getF(), an->getPrevFInt());
    las::VecVecAdd * vva = las::getVecVecAdd<las::MICRO_BACKEND>();
    vva->exec(0.5, an->getFInt(), 0.5, an->getPrevFInt(), tmp);
    an->setInternalEnergy(an->getInternalEnergy() +
                          delta_t / 2.0 * ops->dot(an->getV(), tmp));
    // Wext = Wext_n+delta_t_nphalf/2*v_nphalf^T(f_ext_n+f_ext_npone)
    vva->exec(0.5, an->getFExt(), 0.5, an->getPrevFExt(), tmp);
    an->setExternalEnergy(an->getExternalEnergy() +
                          delta_t / 2.0 * ops->dot(an->getV(), tmp));
    // compute the rayleigh damping matrix
    // FIXME there is some memory shit going on here...shouldn't need to set c
    // to nullptr
    auto mb = las::getMatBuilder<las::MICRO_BACKEND>(0);
    mb->destroy(an->getC());
    c = 0;
    smsma->exec(mass_damping_factor, an->getM(), stiffness_damping_factor,
                an->getK(), &c);
    an->setC(c);
    // add the damping term into the force vector
    mvm->exec(an->getC(), an->getV(), an->getFDamp());
    // Wdmp = Wdmp_n+delta_t_nphalf/2*v_nphalf^T(f_dmp_n+f_dmp_npone)
    vva->exec(0.5, an->getFDamp(), 0.5, an->getPrevFDamp(), tmp);
    an->setDampingEnergy(an->getDampingEnergy() +
                         delta_t / 2.0 * ops->dot(an->getV(), tmp));
    ops->axpy(-1, an->getFDamp(), an->getF());
    // 9. compute a(n+1)
    an->getSlv()->solve(an->getM(), an->getA(), an->getF());
    applyBndNdsToVec(an->getFn()->getANumbering(), an->getFn()->getAField(),
                     an->getA(), an->bnd_nds[RVE::all].begin(),
                     an->bnd_nds[RVE::all].end());
    // 10. partial update of nodal velocities
    // v(n+1)=v(n+1/2)+(t(n+1)-t(n+1/2))a(n+1)
    ops->axpy((t_npone - t_nphalf), an->getA(), an->getV());
    applyBndNdsToVec(an->getFn()->getVNumbering(), an->getFn()->getVField(),
                     an->getV(), an->bnd_nds[RVE::all].begin(),
                     an->bnd_nds[RVE::all].end());
    // 11. check energy balance (actual checking in convergence operators ...
    // just compute it here) W_kin = 1/2*v^TMv
    // Wkin_npone = 1/2*(v_npone)^T.M.v_npone
    mvm->exec(an->getM(), an->getV(), tmp);
    an->setKineticEnergy(0.5 * ops->dot(an->getV(), tmp));
    an->setTime(t_npone);
    // swap the force vectors
    if (__builtin_expect(
            ((this->iteration() % print_steps) == 0) || (last_iter == true),
            false))
    {
      std::cout << "Iteration: " << this->iteration()
                << " Simulation time: " << an->getTime() << "\n";
      std::cout << "Anorm: " << ops->norm(an->getA()) << std::endl;
      std::cout << "Vnorm: " << ops->norm(an->getV()) << std::endl;
      std::cout << "Unorm: " << ops->norm(an->getU()) << std::endl;
      std::cout << "Fnorm: " << ops->norm(an->getF()) << std::endl;
      // std::cout << "FIntNorm: " << "curr: " << ops->norm(an->getFInt()) <<
      //             "  prev: "<< ops->norm(an->getPrevFInt())<<std::endl;
      // std::cout << "FExtNorm: " << "curr: " << ops->norm(an->getFExt()) <<
      //             "  prev: "<< ops->norm(an->getPrevFExt())<<std::endl;
      // std::cout << "FDampNorm: " << "curr: " << ops->norm(an->getFDamp()) <<
      //             "  prev: "<< ops->norm(an->getPrevFDamp())<<std::endl;
      std::cout << "Kinetic Energy: " << an->getKineticEnergy() << std::endl;
      std::cout << "Internal Energy: " << an->getInternalEnergy() << std::endl;
      std::cout << "External Energy: " << an->getExternalEnergy() << std::endl;
      std::cout << "Damping Energy: " << an->getDampingEnergy() << std::endl;
      // only write fields when outputing data
      double * s = 0;
      // write a, v, u data to field
      // TODO all of these field writes could be combined to reduce the number
      // of mesh loops to 1 acceleration field
      ops->get(an->getA(), s);
      amsi::ApplyVector(an->getFn()->getANumbering(), an->getFn()->getAField(),
                        s, 0, &wrt)
          .run();
      ops->restore(an->getA(), s);
      // velocity field
      ops->get(an->getV(), s);
      amsi::ApplyVector(an->getFn()->getVNumbering(), an->getFn()->getVField(),
                        s, 0, &wrt)
          .run();
      ops->restore(an->getV(), s);
      // displacement field
      ops->get(an->getU(), s);
      amsi::ApplyVector(an->getFn()->getUNumbering(), an->getFn()->getUField(),
                        s, 0, &wrt)
          .run();
      ops->restore(an->getU(), s);
      ops->get(an->getF(), s);
      amsi::ApplyVector(an->getFn()->getFNumbering(), an->getFn()->getFField(),
                        s, 0, &wrt)
          .run();
      ops->restore(an->getU(), s);
      std::stringstream sout;
       std::string folder = "test_explicit/";
      //std::string folder = "te/";
      sout << "iter_" << this->iteration();
      apf::writeVtkFiles((folder + sout.str()).c_str(),
                         an->getFn()->getNetworkMesh(), 1);
      pvd_data.push_back(amsi::PvdData(sout.str(), t_npone, 0));
      std::string fname = "test_explicit.pvd";
      amsi::writePvdFile(folder + fname, pvd_data);
    }
    an->updateFExt();
    an->updateFInt();
    an->updateFDamp();
    Iteration::iterate();
  }
  struct EpsGen
  {
    EpsGen(double eps) : eps(eps) {}
    double operator()(int) { return eps; }
    protected:
    double eps;
  };
  struct EnergyBalValGen
  {
    EnergyBalValGen(FiberRVEAnalysisQSExplicit * a) : an(a) {}
    FiberRVEAnalysisQSExplicit * an;
    double operator()()
    {
      double val = fabs(an->getKineticEnergy() + an->getInternalEnergy() -
                        an->getExternalEnergy());
      return val;
    }
  };
  struct EnergyBalRefGen
  {
    EnergyBalRefGen(FiberRVEAnalysisQSExplicit * a) : an(a) {}
    FiberRVEAnalysisQSExplicit * an;
    double operator()()
    {
      return std::max(
          an->getExternalEnergy(),
          std::max(an->getInternalEnergy(), an->getKineticEnergy()));
    }
  };
  struct KineticEnergyValGen
  {
    KineticEnergyValGen(FiberRVEAnalysisQSExplicit * a) : an(a) {}
    FiberRVEAnalysisQSExplicit * an;
    double operator()() { return an->getKineticEnergy(); }
  };
  struct KineticEnergyRefGen
  {
    KineticEnergyRefGen(FiberRVEAnalysisQSExplicit * a) : an(a) {}
    FiberRVEAnalysisQSExplicit * an;
    double operator()() { return an->getTotalEnergy(); }
  };
  struct KineticEnergyEpsGen
  {
    KineticEnergyEpsGen(double eps, double startup_eps, int startup)
        : eps(eps), startup_eps(startup_eps), startup(startup)
    {
    }
    double operator()(int itr) { return (itr > startup) ? eps : startup_eps; }

    protected:
    double eps;
    double startup_eps;
    int startup;
  };
  class TimeConvergence : public amsi::Convergence
  {
    FiberRVEIterationQSExplicit * itr;

    public:
    TimeConvergence(FiberRVEIterationQSExplicit * itr) : itr(itr) {}
    virtual bool converged() { return itr->getCompleted(); }
  };
  struct LinearAmp : public Amplitude
  {
    LinearAmp(double t_max) : t_max(t_max) {}
    double t_max;
    double operator()(double t) { return t / t_max; }
    double derivative(double t) { return 1.0 / t_max; }
    double secondDerivative(double t) { return 0; }
    double integral(double t) { return (t * t) / (2 * t_max); }
    virtual ~LinearAmp() {}
  };
  bool FiberRVEAnalysisQSExplicit::run(const DeformationGradient & dfmGrd)
  {
    double energyBalEps = 0.1;
    double kinEnergyPercent = 10;
    double kinEnergyPercentStartup = 500;
    int kinEnergyStartupSteps = 1000;
    double loadTime = 10;
    double timeStepFactor = 0.95;
    double massDampingFactor = 0;
    double stiffnessDampingFactor = -0.1;
    unsigned int printSteps = 1000;
    EnergyBalRefGen energy_bal_rg(this);
    EnergyBalValGen energy_bal_vg(this);
    EpsGen energy_bal_eg(energyBalEps);
    KineticEnergyRefGen kin_energy_rg(this);
    KineticEnergyValGen kin_energy_vg(this);
    KineticEnergyEpsGen kin_energy_eg(kinEnergyPercent / 100.0,
                                      kinEnergyPercentStartup / 100.0,
                                      kinEnergyStartupSteps);
    LinearAmp linAmp(loadTime);
    FiberRVEIterationQSExplicit itr(this, dfmGrd, &linAmp, loadTime,
                                    timeStepFactor, massDampingFactor,
                                    stiffnessDampingFactor, printSteps);
    TimeConvergence time_cnvrg(&itr);
    amsi::UpdatingConvergence<decltype(&energy_bal_vg),
                              decltype(&energy_bal_eg),
                              decltype(&energy_bal_rg)>
        cnvrg_energ_bal(&itr, &energy_bal_vg, &energy_bal_eg, &energy_bal_rg,
                        true);
    amsi::UpdatingConvergence<decltype(&kin_energy_vg),
                              decltype(&kin_energy_eg),
                              decltype(&kin_energy_rg)>
        cnvrg_kinetic_energy(&itr, &kin_energy_vg, &kin_energy_eg,
                             &kin_energy_rg, true);
    std::vector<amsi::Convergence *> cnvrg_stps;
    // cnvrg_stps.push_back(&cnvrg_energ_bal);
    cnvrg_stps.push_back(&cnvrg_kinetic_energy);
    // this needs to be after the other convergence checks...so they we properly
    // run and fail if the energy balance isn't correct.
    cnvrg_stps.push_back(&time_cnvrg);
    amsi::MultiConvergence cnvrg(cnvrg_stps.begin(), cnvrg_stps.end());
    bool rslt = amsi::numericalSolve(&itr, &cnvrg);
    if (cnvrg_stps[0]->failed())
    {
      std::cerr << "The percentage of kinetic energy was above "
                << kinEnergyPercent
                << " percent of the total energy in iteration "
                << itr.iteration() << std::endl;
      std::cerr << "Kinetic Energy: " << getKineticEnergy() << std::endl;
      std::cerr << "Total Energy: " << getTotalEnergy() << std::endl;
      std::cerr << "Ratio: " << getKineticEnergy() / getTotalEnergy()
                << std::endl;
    }
    return rslt;
  }
}  // namespace bio
