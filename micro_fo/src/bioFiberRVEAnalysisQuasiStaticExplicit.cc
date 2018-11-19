#include <pcu_util.h>
#include "bioFiberRVEAnalysis.h"
#include "bioMassIntegrator.h"
#include <fstream>
#include <sys/stat.h>
#include "bioFiberRVEAnalysisQuasiStaticExplicit.h"
//#include <bioDynamicRelaxation.h>
namespace bio
{
  FiberRVEAnalysisQSExplicit::FiberRVEAnalysisQSExplicit(
      FiberNetwork * fn,
      LinearStructs<las::MICRO_BACKEND> * vecs,
      const MicroSolutionStrategy & ss)
      : FiberRVEAnalysis(fn, vecs, ss)
      , internal_energy(0)
      , kinetic_energy(0)
      , external_energy(0)
      , damping_energy(0)
      , time(0)
      , m_inv(NULL)
  {
    f_int_field = apf::createLagrangeField(fn->getNetworkMesh(), "f_int_field", apf::VECTOR, 1);
    apf::zeroField(f_int_field);
    f_ext_field = apf::createLagrangeField(fn->getNetworkMesh(), "f_ext_field", apf::VECTOR, 1);
    apf::zeroField(f_ext_field);
    f_damp_field = apf::createLagrangeField(fn->getNetworkMesh(), "f_damp_field", apf::VECTOR, 1);
    apf::zeroField(f_damp_field);
    f_inertial_field = apf::createLagrangeField(fn->getNetworkMesh(), "f_inertial_field", apf::VECTOR, 1);
    apf::zeroField(f_inertial_field);
    auto vb = las::getVecBuilder<las::MICRO_BACKEND>(0);
    f_inertial = vb->createRHS(getC());
    apf::Integrator * massIntegrator = new MassIntegrator(
        fn->getUNumbering(), fn->getUField(), this->getM(),
        &(fn->getFiberReactions()[0]), 2, MassLumpType::RowSum);
    massIntegrator->process(this->getFn()->getNetworkMesh(), 1);
    delete massIntegrator;
    // note this requires that we are usign lumped mass
    las::MatDiagonalInverse * mdiag = las::getMatDiagonalInverse<las::MICRO_BACKEND>();
    mdiag->exec(1.0, this->getM(), m_inv);
    // note that the integrator computes the internal forces!
    es = createExplicitMicroElementalSystem(fn, getK(), getFInt(), getFDamp(), getV());
  }
  FiberRVEAnalysisQSExplicit::~FiberRVEAnalysisQSExplicit()
  {
    las::destroyVector(m_inv);
  }
  ExplicitOutputWriter::ExplicitOutputWriter(std::string folder, std::string pvdName, FiberRVEAnalysisQSExplicit * an)
      : outputFrame(0)
      , pvdData(std::vector<amsi::PvdData>())
      , folder(folder)
      , pvdName(pvdName)
      , an(an)
      , ops(las::getLASOps<las::MICRO_BACKEND>())
      {
        // if the folder doesn't exist create it
        std::string fname = folder+"/";
        struct stat sb;
        if(stat(fname.c_str(),&sb) != 0)
        {
          mkdir(fname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        }
        // this clears the micro csv file.
        fname = folder+"/"+"micro_results.csv";
        std::ofstream strm(fname);
      }
      void ExplicitOutputWriter::write(int iteration)
      {
        std::string fname = folder+"/"+"micro_results.csv";
        std::ofstream strm;
        strm.open(fname , std::ios::out | std::ios::app);
        assert(strm.is_open());
        if(__builtin_expect(outputFrame == 0, 0))
          strm<<"Frame, Iteration, Simulation_Time, Kinetic_Energy, Internal_Energy, External_Energy, Damping_Energy, "<<
                     " A_Norm, V_norm, U_norm, F_int_norm, F_ext_norm, F_damp_norm, F_norm\n";
        strm<<outputFrame<<", ";
        strm<<iteration<<", ";
        strm<<an->getTime()<<", ";
        strm<<an->getKineticEnergy()<<", ";
        strm<<an->getInternalEnergy()<<", ";
        strm<<an->getExternalEnergy()<<", ";
        strm<<an->getDampingEnergy()<<", ";
        strm<<ops->norm(an->getA())<<", ";
        strm<<ops->norm(an->getV())<<", ";
        strm<<ops->norm(an->getU())<<", ";
        strm<<ops->norm(an->getFInt())<<", ";
        strm<<ops->norm(an->getFExt())<<", ";
        strm<<ops->norm(an->getFDamp())<<", ";
        strm<<ops->norm(an->getF())<<"\n";
        // only write fields when outputing data
        // FIXME the get functions allocate data which must be freed after use.
        // This is currently a memory leak!
        double * s = 0;
        // write a, v, u data to field
        // TODO all of these field writes could be combined to reduce the number
        // of mesh loops to 1 acceleration field
        amsi::WriteOp wrt;
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
        ops->get(an->getF(), s);
        amsi::ApplyVector(an->getFn()->getFNumbering(), an->getFn()->getFField(),
                          s, 0, &wrt)
            .run();
        // test fields
        ops->get(an->getFInt(), s);
        amsi::ApplyVector(an->getFn()->getFNumbering(), an->getFIntField(),
                          s, 0, &wrt)
            .run();
        ops->get(an->getFExt(), s);
        amsi::ApplyVector(an->getFn()->getFNumbering(), an->getFExtField(),
                          s, 0, &wrt)
            .run();
        ops->get(an->getFDamp(), s);
        amsi::ApplyVector(an->getFn()->getFNumbering(), an->getFDampField(),
                          s, 0, &wrt)
            .run();
        ops->get(an->getFInertial(), s);
        amsi::ApplyVector(an->getFn()->getFNumbering(), an->getFInertialField(),
                          s, 0, &wrt)
            .run();
        ops->get(an->getDeltaU(), s);
        amsi::ApplyVector(an->getFn()->getUNumbering(), an->getFn()->getdUField(),
                          s, 0, &wrt)
            .run();
        std::stringstream sout;
        sout << "frame_" << outputFrame++;
        apf::writeVtkFiles((folder+ "/" + sout.str()).c_str(),
                           an->getFn()->getNetworkMesh(), 1);
        pvdData.push_back(amsi::PvdData(sout.str(), an->getTime(), 0));
        amsi::writePvdFile(folder + "/" + pvdName, pvdData);
      }
  FiberRVEIterationQSExplicit::FiberRVEIterationQSExplicit(
      FiberRVEAnalysisQSExplicit * a,
      DeformationGradient appliedDefm,
      Amplitude * amplitude,
      double loadTime,
      double timeStepFactor,
      double massDampingFactor,
      double stiffnessDampingFactor,
      unsigned int printSteps,
      ExplicitOutputWriter * outputWriter)
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
      , outputWriter(outputWriter)
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
    //ops->zero(an->getC());  // don't bother zeroing because we reset it anyways
    //las::MatMatAdd * smsma = las::getMatMatAdd<las::MICRO_BACKEND>();
    //las::Mat * c;
    las::MatVecMult * mvm = las::getMatVecMult<las::MICRO_BACKEND>();
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
      static_cast<ExplicitTrussIntegrator *>(an->es)->setStiffnessDampingFactor(stiffness_damping_factor);
      outputWriter->write(0);
      an->es->process(an->getFn()->getNetworkMesh(), 1);
    }
    // ITERATIONS START HERE!
    // 4. update time t(n+1) = t(n)+delta_t(n+1/2), t(n+1/2)=1/2(t(n)+t(n+1))
    double delta_t =
        static_cast<ExplicitTrussIntegrator *>(an->es)->getCriticalTimeStep() *
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
    // get delta_u
    ops->zero(an->getDeltaU());
    ops->axpy(delta_t, an->getV(), an->getDeltaU());
    applyBndNdsToVec(an->getFn()->getUNumbering(), an->getFn()->getdUField(),
                     an->getDeltaU(), an->bnd_nds[RVE::all].begin(),
                     an->bnd_nds[RVE::all].end());
    // apply u boundary conditions
    // FIXME...these boundary conditions also need to get applied to DeltaU
    applyBndNdsToVec(an->getFn()->getUNumbering(), an->getFn()->getUField(),
                     an->getU(), an->bnd_nds[RVE::all].begin(),
                     an->bnd_nds[RVE::all].end());
    // 8. getforce
    //ops->zero(an->getK()); // no longer explicitly compute K
    // FIXME add mass damping...need to multiply by scalar...
    //mvm->exec(an->getM(), an->getV(), an->getFDamp());
    ops->zero(an->getFDamp());
    ops->zero(an->getFInt());
    an->es->process(an->getFn()->getNetworkMesh(), 1);
    // compute the rayleigh damping matrix
    //auto mb = las::getMatBuilder<las::MICRO_BACKEND>(0);
    //mb->destroy(an->getC());
    //c = 0;
    //smsma->exec(mass_damping_factor, an->getM(), stiffness_damping_factor,
    //            an->getK(), &c);
    //an->setC(c);
    // compute the damping force
    //mvm->exec(an->getC(), an->getV(), an->getFDamp());
    ops->zero(an->getF());
    ops->axpy(-1, an->getFInt(), an->getF());
    ops->axpy(-1, an->getFDamp(), an->getF());
    // we don't allow the use to apply a force, so this is not needed here
    //ops->axpy(1, an->getFExt(), an->getF());
    // sum the forces to the total force vector (note -FInt is summed in the force integrator)
    // 9. compute a(n+1)
    //an->getSlv()->solve(an->getM(), an->getA(), an->getF());
    las::HadamardProduct * hp = las::getHadamardProduct<las::MICRO_BACKEND>();
    hp->exec(an->getMInv(), an->getF(), an->getA());
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
    ops->zero(an->getFExt()); 
    ops->axpy(1, an->getFInt(), an->getFExt());
    ops->axpy(1, an->getFDamp(), an->getFExt());
    // compute inertial forces and add into external work
    mvm->exec(an->getM(), an->getA(), an->getFInertial());
    ops->axpy(1, an->getFInertial(), an->getFExt());
    //mvm->exec(an->getM(), an->getA(), tmp);
    //ops->axpy(1, tmp, an->getFExt());
    // Wint_n+delta_t_nphalf/2*v_nphalf^T(f_int_n+f_int_npone)
    las::VecVecAdd * vva = las::getVecVecAdd<las::MICRO_BACKEND>();
    vva->exec(0.5, an->getFInt(), 0.5, an->getPrevFInt(), tmp);
    an->setInternalEnergy(an->getInternalEnergy() +
                          0.5 * ops->dot(an->getDeltaU(), tmp));
    // Wext = Wext_n+delta_t_nphalf/2*v_nphalf^T(f_ext_n+f_ext_npone)
    vva->exec(0.5, an->getFExt(), 0.5, an->getPrevFExt(), tmp);
    an->setExternalEnergy(an->getExternalEnergy() +
                          0.5 * ops->dot(an->getDeltaU(), tmp));
    // Wdmp = Wdmp_n+delta_t_nphalf/2*v_nphalf^T(f_dmp_n+f_dmp_npone)
    vva->exec(0.5, an->getFDamp(), 0.5, an->getPrevFDamp(), tmp);
    an->setDampingEnergy(an->getDampingEnergy() +
                         0.5 * ops->dot(an->getDeltaU(), tmp));
    // just compute it here) W_kin = 1/2*v^TMv
    // Wkin_npone = 1/2*(v_npone)^T.M.v_npone
    mvm->exec(an->getM(), an->getV(), tmp);
    an->setKineticEnergy(0.5 * ops->dot(an->getV(), tmp));
    an->setTime(t_npone);
    if (__builtin_expect(
            ((this->iteration() % print_steps) == 0) || (last_iter == true),
            false))
    {
      outputWriter->write(iteration()+1);
    }
    // swap the force vectors
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
    double loadTime = 0.5E-2;
    //double loadTime = 0.5E-1;
    //double loadTime = 10;
    double timeStepFactor = 0.95;
    double massDampingFactor = 0;
    double stiffnessDampingFactor = 0.1;
    unsigned int printSteps = 10000;
    ExplicitOutputWriter writer("aba/", "test_explicit.pvd", this);
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
                                    stiffnessDampingFactor, printSteps, &writer);
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
    //cnvrg_stps.push_back(&cnvrg_kinetic_energy);
    // this needs to be after the other convergence checks...so they we properly
    // run and fail if the energy balance isn't correct.
    cnvrg_stps.push_back(&time_cnvrg);
    amsi::MultiConvergence cnvrg(cnvrg_stps.begin(), cnvrg_stps.end());
    bool rslt = amsi::numericalSolve(&itr, &cnvrg);
    if (cnvrg_kinetic_energy.failed())
    {
      std::cerr << "The percentage of kinetic energy was above "
                << kinEnergyPercent
                << " percent of the total energy in iteration "
                << itr.iteration() << std::endl;
      std::cerr << "Kinetic Energy: " << getKineticEnergy() << std::endl;
      std::cerr << "Total Energy: " << getTotalEnergy() << std::endl;
      std::cerr << "Ratio: " << getKineticEnergy() / getTotalEnergy()
                << std::endl;
      writer.write(itr.iteration()+1);
    }
    return rslt;
  }
  apf::Integrator * createExplicitMicroElementalSystem(FiberNetwork * fn,
                                               las::Mat * k,
                                               las::Vec * f,
                                               las::Vec * f_damp,
                                               las::Vec * v)
  {
    apf::Integrator * es = NULL;
    FiberMember tp = fn->getFiberMember();
    if (tp == FiberMember::truss) {
        es = new ExplicitTrussIntegrator(fn->getUNumbering(), fn->getUField(),
                                 fn->getXpUField(), &(fn->getFiberReactions()[0]),
                                 k, f, f_damp, v, 1);
      }
    return es;
  }
}  // namespace bio
