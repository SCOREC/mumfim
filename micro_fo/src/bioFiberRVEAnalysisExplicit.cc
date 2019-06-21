#include "bioFiberRVEAnalysisExplicit.h"
#include <PCU.h>
#include <apf.h>
#include <apfConvert.h>
#include <apfMesh2.h>
#include <las.h>
#include <lionPrint.h>
#include <mpi.h>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include "bioExplicitAmplitude.h"
#include "bioExplicitOutputWriter.h"
#include "bioFiberNetworkIO.h"
#include "bioFiberRVEAnalysisExplicit_impl.h"
#include "bioMassIntegrator.h"
#include "bioMicroFOConfig.h"
#ifdef ENABLE_KOKKOS
#include <Kokkos_Core.hpp>
#endif
namespace bio
{
  FiberRVEAnalysisExplicit::FiberRVEAnalysisExplicit(FiberNetwork * fn,
                             LinearStructs<las::MICRO_BACKEND> * vecs,
                             const MicroSolutionStrategyExplicit & ss)
        : FiberRVEAnalysis(fn, vecs, static_cast<MicroSolutionStrategy>(ss))
        , visc_damp_coeff(ss.visc_damp_coeff)
        , print_history_frequency(ss.print_history_frequency)
        , print_field_frequency(ss.print_field_frequency)
        , print_field_by_num_frames(ss.print_field_by_num_frames)
        , total_time(ss.total_time)
        , load_time(ss.load_time)
        , serial_gpu_cutoff(ss.serial_gpu_cutoff)
        , crit_time_scale_factor(ss.crit_time_scale_factor)
        , energy_check_eps(ss.energy_check_eps)
        , disp_bound_nfixed(0)
        , disp_bound_dof(NULL)
        , disp_bound_vals(NULL)
        , fiber_elastic_modulus(fn->getFiberReactions()[0]->E)
        , fiber_area(fn->getFiberReactions()[0]->fiber_area)
        , fiber_density(fn->getFiberReactions()[0]->fiber_density)
        , amp(NULL)
        , system_initialized(false)
        , itr_prev(0)
    {
      es = createImplicitMicroElementalSystem(fn, getK(), getF());
      std::stringstream sout;
      int rnk = -1;
      MPI_Comm_rank(MPI_COMM_WORLD, &rnk);
      sout << "rnk_" << rnk << "_fn_" << getFn()->getRVEType()<<"_explicit";
      analysis_name = sout.str();
      std::string dir = amsi::fs ? amsi::fs->getResultsDir() : ".";
      writer = new ExplicitOutputWriter(fn->getNetworkMesh(),
                 (dir + "/" + analysis_name).c_str(),
                 (analysis_name + ".pvd").c_str());
      if(ss.ampType == AmplitudeType::SmoothStep)
      {
        amp = new SmoothAmp(total_time);
      }
      else if(ss.ampType == AmplitudeType::SmoothStepHold)
      {
        amp = new SmoothAmpHold(load_time,total_time);
      }
      else
      {
        std::cerr<<"Incorrect Amplitude type selected ("<<static_cast<int>(ss.ampType)<<")"<<std::endl;
        std::abort();
      }
      //amp = new SmoothAmpHold(1.0, total_time);
      assert(fiber_density > 0);
      assert(fiber_area > 0);
      assert(fiber_elastic_modulus > 0);
    }
  // FIXME Note this is the same as a single implicit iteration...
  // we should probably use the same function call that is there if this works
  void FiberRVEAnalysisExplicit::relaxSystem() {
    auto ops = las::getLASOps<las::MICRO_BACKEND>();
    ops->zero(this->getK());
    ops->zero(this->getU());
    ops->zero(this->getF());
    apf::Mesh * fn = this->getFn()->getNetworkMesh();
    this->es->process(fn, 1);
    // finalize the vectors so we cthis set boundary condition
    // values
    las::finalizeMatrix<las::MICRO_BACKEND>(this->vecs->getK());
    las::finalizeVector<las::MICRO_BACKEND>(this->vecs->getF());
    applyRVEBC(this->bnd_nds[RVE::all].begin(),
               this->bnd_nds[RVE::all].end(),
               this->getFn()->getUNumbering(),
               this->getK(),
               this->getF());
    this->getSlv()->solve(this->getK(), this->getU(), this->getF());
    //amsi::SubtractOp acm;
    //amsi::WriteScalarMultOp mlt(-1.0);
    amsi::WriteOp wrt;
    amsi::AccumOp acm;
    amsi::FreeApplyOp fr_mlt(this->getFn()->getUNumbering(), &wrt);
    amsi::FreeApplyOp fr_acm(this->getFn()->getUNumbering(), &acm);
    double * s = NULL;
    ops->get(this->getU(), s);
    amsi::ApplyVector(
        this->getFn()->getUNumbering(), this->getFn()->getdUField(), s, 0, &fr_mlt)
        .run();
    amsi::ApplyVector(
        this->getFn()->getUNumbering(), this->getFn()->getUField(), s, 0, &fr_acm)
        .run();
    ops->restore(this->getU(), s);
  }
  bool FiberRVEAnalysisExplicit::run(const DeformationGradient & dfmGrd)
  {
    BIO_V3(
    std::cout << "F=[";
    std::cout << dfmGrd[0] << " " << dfmGrd[1] << " " << dfmGrd[2] << "; ";
    std::cout << dfmGrd[3] << " " << dfmGrd[4] << " " << dfmGrd[5] << "; ";
    std::cout << dfmGrd[6] << " " << dfmGrd[7] << " " << dfmGrd[8];
    std::cout << "]\n";
    );
    computeDisplcamentBC(dfmGrd);
    apf::Mesh * mesh = getFn()->getNetworkMesh();
    // we do look for the field because if a previous analysis
    // was run, then we need to let the ETFEM code know
    apf::Field * massField = mesh->findField("nodalMass");
    // TODO check if we set the u field to something different somewhere...
    // possibly the displacements aren't the same when we leave and return to the solver?
    // just for testing purposes, we are using the du field since it isn't
    // used otherwise
    apf::Field * f_ext_field = getFn()->getdUField();
    apf::zeroField(getFn()->getVField());
    apf::zeroField(getFn()->getAField());
    //apf::zeroField(getFn()->getFField());
    //apf::zeroField(getFn()->getdUField());
    //apf::zeroField(getFn()->getUField());
    //apf::freeze(getFn()->getUField());
    apf::Field * coords = getFn()->getNetworkMesh()->getCoordinateField(); 
    apf::freeze(getFn()->getUField());
    double * u_arr = apf::getArrayData(getFn()->getUField());
    bool rtn;
#ifdef ENABLE_KOKKOS 
    if (this->getFn()->getDofCount() < serial_gpu_cutoff)
    {
#endif
      // run serial analysis
      ExplicitAnalysisSerial analysis(
          static_cast<apf::Mesh2 *>(mesh),
          total_time, fiber_elastic_modulus, fiber_area, fiber_density, amp,
          analysis_name, visc_damp_coeff, print_history_frequency,
          print_field_frequency, print_field_by_num_frames,
          crit_time_scale_factor, energy_check_eps,coords,
          getFn()->getUField(), getFn()->getVField(), getFn()->getAField(),
          getFn()->getFField(), f_ext_field, massField, writer);
      assert(disp_bound_dof && disp_bound_vals);
      analysis.setDispBC(disp_bound_nfixed, disp_bound_dof, disp_bound_init_vals, disp_bound_vals);
      rtn = analysis.run(itr_prev);
      system_initialized = true;
#ifdef ENABLE_KOKKOS 
    }
    else
    {
      // run gpu analysis
      ExplicitAnalysisKokkos analysis(
          static_cast<apf::Mesh2 *>(mesh),
          total_time, fiber_elastic_modulus, fiber_area, fiber_density, amp,
          analysis_name, visc_damp_coeff, print_history_frequency,
          print_field_frequency, print_field_by_num_frames,
          crit_time_scale_factor, energy_check_eps, coords,
          getFn()->getUField(), getFn()->getVField(), getFn()->getAField(),
          getFn()->getFField(), f_ext_field, massField, writer);
      assert(disp_bound_dof && disp_bound_vals);
      analysis.setDispBC(disp_bound_nfixed, disp_bound_dof, disp_bound_init_vals, disp_bound_vals);
      rtn = analysis.run(itr_prev);
      system_initialized = true;
    }
#endif
    // delete the boundary condition data
    // TODO fix this so we aren't newing and deleting these arrays every iteration
    delete [] disp_bound_dof;
    disp_bound_dof=NULL;
    delete [] disp_bound_vals;
    disp_bound_vals=NULL;
    delete [] disp_bound_init_vals;
    disp_bound_init_vals = NULL;
    // accumulate the du from the explicit analysis into the u field
    //apf::freeze(getFn()->getdUField());
    //double * du = apf::getArrayData(getFn()->getdUField());
    //amsi::AccumOp acm;
    //amsi::ApplyVector(getFn()->getUNumbering(), getFn()->getUField(),  du, 0, &acm).run();
    //relaxSystem();
    return rtn;
  }
  void FiberRVEAnalysisExplicit::computeDisplcamentBC(
      const DeformationGradient & dfmGrd)
  {
    disp_bound_nfixed = 3 * bnd_nds[RVE::all].size();
    disp_bound_dof = new int[disp_bound_nfixed];
    disp_bound_vals = new double[disp_bound_nfixed];
    disp_bound_init_vals = new double[disp_bound_nfixed];
    //apf::Field * current_coords_field = getFn()->getXpUField();
    apf::Field * current_coords_field = getFn()->getNetworkMesh()->getCoordinateField();
    apf::Field * disp_field = getFn()->getUField();
    apf::Vector3 coords, disp;
    int node_num[3];
    for (std::size_t i = 0; i < bnd_nds[RVE::all].size(); ++i)
    {
      apf::MeshEntity * nd = bnd_nds[RVE::all][i];
      apf::getVector(current_coords_field, nd, 0, coords);
      apf::getVector(disp_field, nd, 0, disp);
      node_num[0] = apf::getNumber(getFn()->getUNumbering(), nd, 0, 0);
      node_num[1] = apf::getNumber(getFn()->getUNumbering(), nd, 0, 1);
      node_num[2] = apf::getNumber(getFn()->getUNumbering(), nd, 0, 2);
      disp_bound_dof[i * 3] = node_num[0];
      disp_bound_dof[i * 3 + 1] = node_num[1];
      disp_bound_dof[i * 3 + 2] = node_num[2];
      disp_bound_vals[i * 3] = (dfmGrd[0] - 1) * coords[0] +
                               dfmGrd[1] * coords[1] + dfmGrd[2] * coords[2]-disp[0];
      disp_bound_vals[i * 3 + 1] = dfmGrd[3] * coords[0] +
                                   (dfmGrd[4] - 1) * coords[1] +
                                   dfmGrd[5] * coords[2]-disp[1];
      disp_bound_vals[i * 3 + 2] = dfmGrd[6] * coords[0] +
                                   dfmGrd[7] * coords[1] +
                                   (dfmGrd[8] - 1) * coords[2]-disp[2];
      disp_bound_init_vals[i*3] = disp[0];
      disp_bound_init_vals[i*3+1] = disp[1];
      disp_bound_init_vals[i*3+2] = disp[2];
    }
  }
  void FiberRVEAnalysisExplicit::copyForceDataToForceVec()
  {
    auto ops = las::getLASOps<las::MICRO_BACKEND>();
    apf::freeze(getFn()->getFField());
    double * f_arr = apf::getArrayData(getFn()->getFField());
    ops->set(getF(), f_arr);
  }
  void FiberRVEAnalysisExplicit::copyDispDataToDispVec()
  {
    auto ops = las::getLASOps<las::MICRO_BACKEND>();
    apf::freeze(getFn()->getUField());
    double * u_arr = apf::getArrayData(getFn()->getUField());
    ops->set(getU(), u_arr);
  }
  void FiberRVEAnalysisExplicit::computeStiffnessMatrix()
  {
    auto ops = las::getLASOps<las::MICRO_BACKEND>();
    ops->zero(getK());
    ops->zero(getF());
    apf::Integrator * es =
        createImplicitMicroElementalSystem(getFn(), getK(), getF());
    es->process(getFn()->getNetworkMesh(), 1);
    delete es;
  }
}  // namespace bio
