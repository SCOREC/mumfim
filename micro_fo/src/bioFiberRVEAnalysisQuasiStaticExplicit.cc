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
  class ApplyDeformationGradientExplicit : public ApplyDeformationGradient
  {
    protected:
    apf::Matrix3x3 dVdX;
    apf::Matrix3x3 dAdX;
    apf::Field * vField;
    apf::Field * aField;
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
                                     apf::Field * aField)
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
      //apf::setVector(vField, ent, nd, nd_v);
      //apf::setVector(aField, ent, nd, nd_a);
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
    es = createMicroElementalSystem(fn, getK(), getF(),
                                    FiberRVEAnalysisType::QuasiStaticExplicit);
    // apf::zeroField(fn->getVField());
    // apf::zeroField(fn->getAField());
  }
  FiberRVEIterationQSExplicit::FiberRVEIterationQSExplicit(
      FiberRVEAnalysisQSExplicit * a,
      DeformationGradient appliedDefm,
      Amplitude * amplitude, double loadTime)
      : amsi::Iteration()
      , an(a)
      , appliedDefm(appliedDefm)
      , amplitude(amplitude)
      , last_iter(false)
      , print_steps(1000)
      , load_time(loadTime)
  {
    u_arr = new double[an->getFn()->getDofCount()];
    v_arr = new double[an->getFn()->getDofCount()];
    a_arr = new double[an->getFn()->getDofCount()];
    f_arr = new double[an->getFn()->getDofCount()];
  }
  // see belytschko box 6.1
  void FiberRVEIterationQSExplicit::iterate()
  {
    // TODO check deformation wave speed
    auto ops = las::getLASOps<las::MICRO_BACKEND>();
    ops->zero(an->getC());  // don't bother zeroing because we reset it anyways
    // should be between 0.8 and 0.95 (Belytscheko)
    double alpha = 0.8;
    //// scale terms for Raleigh Damping
    double mC = 0;  // mass scale term
    //double kC = -0.1;  // stiffness scale term
    double kC = 0;  // stiffness scale term
    las::ScalarMatScalarMatAdd * smsma =
        las::getSparskitScalarMatScalarMatAdd();
    las::Mat * c;
    las::MatVecMult * mvm = las::getSparskitMatVecMult();
    las::LasCreateVec * vb = las::getVecBuilder<las::MICRO_BACKEND>(0);
    las::Vec * tmp = vb->createRHS(an->getC());
    amsi::WriteOp wrt;
    amsi::AccumOp accum;
    if (__builtin_expect(this->iteration() == 0, false))
    {
      ops->zero(an->getU());
      ops->zero(an->getV());
      ops->zero(an->getA());

      ops->zero(an->getK());
      ops->zero(an->getF());
      an->es->process(an->getFn()->getNetworkMesh(), 1);
      ops->set(an->getU(), u_arr);
      ops->set(an->getV(), v_arr);
      ops->set(an->getA(), a_arr);
      /*
      // ONLY FIRST ITERATION
      // 1. initial conditions and initialization
      // TODO use ToArray to put a,v,u into the linear structs, might not always
      // be necessary if maintain linear structs for each fiber network
      // individually
      amsi::ToArray(an->getFn()->getUNumbering(), an->getFn()->getUField(),
                    u_arr, 0, &wrt)
          .run();
      amsi::ToArray(an->getFn()->getVNumbering(), an->getFn()->getVField(),
                    v_arr, 0, &wrt)
          .run();
      amsi::ToArray(an->getFn()->getANumbering(), an->getFn()->getAField(),
                    a_arr, 0, &wrt)
          .run();
      // 2. getforce
      ops->zero(an->getF());
      an->es->process(an->getFn()->getNetworkMesh(), 1);
      // 3. compute accelerations a(n) = inv(M)(f-Cv(n-1/2))
      mvm->exec(an->getC(), an->getV(), tmp);
      ops->axpy(-1, tmp, an->getF());
      // note the force vector here is f_ext-f_int-C*v
      an->getSlv()->solve(an->getM(), an->getA(), an->getF());
      */
    }
    //std::cout<<"("<<iteration()<<") F_norm_1: " << ops->norm(an->getF())<<std::endl;
    // ITERATIONS START HERE!
    // 4. update time t(n+1) = t(n)+delta_t(n+1/2), t(n+1/2)=1/2(t(n)+t(n+1))
    double delta_t =
        dynamic_cast<ExplicitTrussIntegrator *>(an->es)->getCriticalTimeStep() *
        alpha;
    if (__builtin_expect(this->iteration() == 0, false))
    {
      std::cout<<"Assuming a constant timestep (dt="<<delta_t<<"),\nthis should take approx: "<<static_cast<unsigned long>(load_time/delta_t)<<" timesteps";
      std::cout<<" to acheive T_max="<<load_time<<std::endl;
    }
    double t_npone = an->getTime() + delta_t;
    if (__builtin_expect(t_npone >= load_time,false))
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
    //amsi::ToArray(an->getFn()->getVNumbering(), an->getFn()->getVField(), v_arr,
    //              0, &wrt)
    //    .run(an->bnd_nds[RVE::all].begin(), an->bnd_nds[RVE::all].end());
    // 7. update nodal displacements u(n+1) = u(n)+delta_t(n+1/2)*v(n+1/2)
    ops->axpy(delta_t, an->getV(), an->getU());
    // apply u boundary conditions
    amsi::ToArray(an->getFn()->getUNumbering(), an->getFn()->getUField(), u_arr,
                  0, &wrt)
        .run(an->bnd_nds[RVE::all].begin(), an->bnd_nds[RVE::all].end());
    // 8. getforce
    ops->zero(an->getK());
    ops->zero(an->getF());
    an->es->process(an->getFn()->getNetworkMesh(), 1);
    // compute the rayleigh damping matrix
    smsma->exec(mC, an->getM(), kC, an->getK(), &c);
    an->setC(c);
    // add the damping term into the force vector
    mvm->exec(an->getC(), an->getV(), tmp);
    ops->axpy(-1, tmp, an->getF());
    // 9. compute a(n+1)
    an->getSlv()->solve(an->getM(), an->getA(), an->getF());
    amsi::ToArray(an->getFn()->getANumbering(), an->getFn()->getAField(), a_arr,
                  0, &wrt)
        .run(an->bnd_nds[RVE::all].begin(), an->bnd_nds[RVE::all].end());
    // 10. partial update of nodal velocities
    // v(n+1)=v(n+1/2)+(t(n+1)-t(n+1/2))a(n+1)
    ops->axpy((t_npone - t_nphalf), an->getA(), an->getV());
    amsi::ToArray(an->getFn()->getVNumbering(), an->getFn()->getVField(), v_arr,
                  0, &wrt)
        .run(an->bnd_nds[RVE::all].begin(), an->bnd_nds[RVE::all].end());
    // 11. check energy balance (actual checking in convergence operators ...
    // just compute it here) W_kin = 1/2*v^TMv
    an->setTime(t_npone);
    if (__builtin_expect(((this->iteration() % print_steps) == 0) || (last_iter == true),false))
    {
      apf::Matrix3x3 I(1,0,0,0,1,0,0,0,1);
      std::cout << "Iteration: " << this->iteration()<<" Simulation time: " << an->getTime() << "\n";
      std::cout<<"Anorm: "<<ops->norm(an->getA())<<std::endl;
      std::cout<<"Vnorm: "<<ops->norm(an->getV())<<std::endl;
      std::cout<<"Unorm: "<<ops->norm(an->getU())<<std::endl;
      std::cout<<"Fnorm: "<<ops->norm(an->getF())<<std::endl;
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
      std::string folder = "test_explicit_no_damp/";
      sout << "iter_" << this->iteration();
      apf::writeVtkFiles((folder + sout.str()).c_str(),
                         an->getFn()->getNetworkMesh(), 1);
      pvd_data.push_back(amsi::PvdData(sout.str(), t_npone, 0));
      std::string fname = "test_explicit_no_damp.pvd";
      amsi::writePvdFile(folder + fname, pvd_data);
    }
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
    double operator()() { std::cout<<"Time convergence"<<std::endl;return an->getTotalEnergy(); }
  };
  class TimeConvergence : public amsi::Convergence
  {
    FiberRVEIterationQSExplicit* itr;

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
    double kinEnergyPercent = 5.0;
    double loadTime = 10.0;
    EnergyBalRefGen energy_bal_rg(this);
    EnergyBalValGen energy_bal_vg(this);
    EpsGen energy_bal_eg(energyBalEps);
    KineticEnergyRefGen kin_energy_rg(this);
    KineticEnergyValGen kin_energy_vg(this);
    EpsGen kin_energy_eg(kinEnergyPercent / 100.0);
    LinearAmp linAmp(loadTime);
    FiberRVEIterationQSExplicit itr(this, dfmGrd, &linAmp, loadTime);
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
    //cnvrg_stps.push_back(&cnvrg_energ_bal);
    //cnvrg_stps.push_back(&cnvrg_kinetic_energy);
    // this needs to be after the other convergence checks...so they we properly
    // run and fail if the energy balance isn't correct.
    cnvrg_stps.push_back(&time_cnvrg);
    amsi::MultiConvergence cnvrg(cnvrg_stps.begin(), cnvrg_stps.end());
    return amsi::numericalSolve(&itr, &cnvrg);
  }
}  // namespace bio
