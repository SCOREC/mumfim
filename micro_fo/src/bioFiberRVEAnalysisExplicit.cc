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
#include <Kokkos_Core.hpp>
namespace bio
{
  FiberRVEAnalysisExplicit::FiberRVEAnalysisExplicit(std::unique_ptr<FiberNetwork> fn,
                             std::unique_ptr<MicroSolutionStrategy> ss,
                             std::shared_ptr<void> workspace)
        : FiberRVEAnalysis(std::move(fn),std::move(ss), workspace)
        , serial_gpu_cutoff(static_cast<MicroSolutionStrategyExplicit*>(mSolutionStrategy.get())->serial_gpu_cutoff)
        , total_time(static_cast<MicroSolutionStrategyExplicit*>(mSolutionStrategy.get())->total_time)
        , load_time(static_cast<MicroSolutionStrategyExplicit*>(mSolutionStrategy.get())->load_time)
        , fiber_elastic_modulus(mFiberNetwork->getFiberReaction(0).E)
        , fiber_area(mFiberNetwork->getFiberReaction(0).fiber_area)
        , fiber_density(mFiberNetwork->getFiberReaction(0).fiber_density)
        , amp(NULL)
        , visc_damp_coeff(static_cast<MicroSolutionStrategyExplicit*>(mSolutionStrategy.get())->visc_damp_coeff)
        , print_history_frequency(static_cast<MicroSolutionStrategyExplicit*>(mSolutionStrategy.get())->print_history_frequency)
        , print_field_frequency(static_cast<MicroSolutionStrategyExplicit*>(mSolutionStrategy.get())->print_field_frequency)
        , print_field_by_num_frames(static_cast<MicroSolutionStrategyExplicit*>(mSolutionStrategy.get())->print_field_by_num_frames)
        , crit_time_scale_factor(static_cast<MicroSolutionStrategyExplicit*>(mSolutionStrategy.get())->crit_time_scale_factor)
        , energy_check_eps(static_cast<MicroSolutionStrategyExplicit*>(mSolutionStrategy.get())->energy_check_eps)
        , disp_bound_nfixed(0)
        , disp_bound_dof(NULL)
        , disp_bound_vals(NULL)
        , system_initialized(false)
        , itr_prev(0)
    {
      std::stringstream sout;
      int rnk = -1;
      MPI_Comm_rank(AMSI_COMM_WORLD, &rnk);
      sout << "rnk_" << rnk << "_fn_" << getFn()->getRVEType()<<"_explicit";
      analysis_name = sout.str();
      std::string dir = amsi::fs ? amsi::fs->getResultsDir() : ".";
      writer = new ExplicitOutputWriter(mFiberNetwork->getNetworkMesh(),
                 (dir + "/" + analysis_name).c_str(),
                 (analysis_name + ".pvd").c_str());
      if(static_cast<MicroSolutionStrategyExplicit*>(mSolutionStrategy.get())->ampType == AmplitudeType::SmoothStep)
      {
        amp = new SmoothAmp(total_time);
      }
      else if(static_cast<MicroSolutionStrategyExplicit*>(mSolutionStrategy.get())->ampType == AmplitudeType::SmoothStepHold)
      {
        amp = new SmoothAmpHold(load_time,total_time);
      }
      else
      {
        std::cerr<<"Incorrect Amplitude type selected ("<<static_cast<int>(static_cast<MicroSolutionStrategyExplicit*>(mSolutionStrategy.get())->ampType)<<")"<<std::endl;
        std::abort();
      }
      //amp = new SmoothAmpHold(1.0, total_time);
      assert(fiber_density > 0);
      assert(fiber_area > 0);
      assert(fiber_elastic_modulus > 0);
      old_u = std::vector<double>();
    }
  bool FiberRVEAnalysisExplicit::run(const DeformationGradient & dfmGrd, double sigma[6], bool update_coords)
  {
    BIO_V3(
    std::cout << "F=[";
    std::cout << dfmGrd[0] << " " << dfmGrd[1] << " " << dfmGrd[2] << "; ";
    std::cout << dfmGrd[3] << " " << dfmGrd[4] << " " << dfmGrd[5] << "; ";
    std::cout << dfmGrd[6] << " " << dfmGrd[7] << " " << dfmGrd[8];
    std::cout << "]\n";
    );
    if(!update_coords)
    {
      old_u.resize(getFn()->getDofCount());
      old_rve_u.resize(getRVE()->getDofCount());
      old_rve_du.resize(getRVE()->getDofCount());
      amsi::WriteOp wrt;
      // for the fiber network only bother save the U field since we don't currently use the
      // dU field (other than for debugging purposes)
      amsi::ToArray(getFn()->getUNumbering(), getFn()->getUField(), &old_u[0], 0, &wrt).run();
      amsi::ToArray(getRVE()->getNumbering(), getRVE()->getUField(), &old_rve_u[0], 0, &wrt).run();
      amsi::ToArray(getRVE()->getNumbering(), getRVE()->getdUField(), &old_rve_du[0], 0, &wrt).run();
    }
    // see what happens when we initially apply affine solution
    applyGuessSolution(this, dfmGrd);
    computeDisplacementBC(dfmGrd);
    // since for debugging use we are putting the external forces into the du field
    // we need to zero that field after apply the guess solution...otherwise we apply
    // external forces at every point in the field which is bad...
    apf::zeroField(getFn()->getdUField());
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
    apf::Field * coords = getFn()->getNetworkMesh()->getCoordinateField(); 
    bool rtn;
    if (this->getFn()->getDofCount() < serial_gpu_cutoff)
    {
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
    computeCauchyStress(sigma);
    if(!update_coords)
    {
      amsi::WriteOp wrt;
      amsi::ApplyVector(getFn()->getUNumbering(), getFn()->getUField(), &old_u[0], 0, &wrt).run();
      amsi::ApplyVector(getRVE()->getNumbering(), getRVE()->getUField(), &old_rve_u[0], 0, &wrt).run();
      amsi::ApplyVector(getRVE()->getNumbering(), getRVE()->getdUField(), &old_rve_du[0], 0, &wrt).run();
    }
    else
    {
      for(int i=0; i<6; ++i)
        curStress[i] = sigma[i];
    }
    // delete the boundary condition data
    // TODO fix this so we aren't newing and deleting these arrays every iteration
    delete [] disp_bound_dof;
    disp_bound_dof=NULL;
    delete [] disp_bound_vals;
    disp_bound_vals=NULL;
    delete [] disp_bound_init_vals;
    disp_bound_init_vals = NULL;
    return rtn;
  }
  /*
  void FiberRVEAnalysisExplicit::computeDisplacementBC(
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
     // version 1
      disp_bound_vals[i * 3] = (dfmGrd[0] - 1) * coords[0] +
                               dfmGrd[1] * coords[1] + dfmGrd[2] * coords[2]-disp[0];
      disp_bound_vals[i * 3 + 1] = dfmGrd[3] * coords[0] +
                                   (dfmGrd[4] - 1) * coords[1] +
                                   dfmGrd[5] * coords[2]-disp[1];
      disp_bound_vals[i * 3 + 2] = dfmGrd[6] * coords[0] +
                                   dfmGrd[7] * coords[1] +
                                   (dfmGrd[8] - 1) * coords[2]-disp[2];
     // Version 2
      //disp_bound_vals[i * 3] = (dfmGrd[0] - 1) * coords[0] +
      //                         dfmGrd[1] * coords[1] + dfmGrd[2] * coords[2];
      //disp_bound_vals[i * 3 + 1] = dfmGrd[3] * coords[0] +
      //                             (dfmGrd[4] - 1) * coords[1] +
      //                             dfmGrd[5] * coords[2];
      //disp_bound_vals[i * 3 + 2] = dfmGrd[6] * coords[0] +
      //                             dfmGrd[7] * coords[1] +
      //                             (dfmGrd[8] - 1) * coords[2];
      disp_bound_init_vals[i*3] = disp[0];
      disp_bound_init_vals[i*3+1] = disp[1];
      disp_bound_init_vals[i*3+2] = disp[2];
    }
  }
  */
  void FiberRVEAnalysisExplicit::computeDisplacementBC(
      const DeformationGradient & dfmGrd)
  {
    disp_bound_nfixed = 3 * bnd_nds[RVE::all].size();
    disp_bound_dof = new int[disp_bound_nfixed];
    disp_bound_vals = new double[disp_bound_nfixed];
    disp_bound_init_vals = new double[disp_bound_nfixed];
    apf::Field * du_field = getFn()->getdUField();
    apf::Field * disp_field = getFn()->getUField();
    apf::Vector3 du, disp;
    int node_num[3];
    for (std::size_t i = 0; i < bnd_nds[RVE::all].size(); ++i)
    {
      apf::MeshEntity * nd = bnd_nds[RVE::all][i];
      apf::getVector(du_field, nd, 0, du);
      apf::getVector(disp_field, nd, 0, disp);
      node_num[0] = apf::getNumber(getFn()->getUNumbering(), nd, 0, 0);
      node_num[1] = apf::getNumber(getFn()->getUNumbering(), nd, 0, 1);
      node_num[2] = apf::getNumber(getFn()->getUNumbering(), nd, 0, 2);
      disp_bound_dof[i * 3] = node_num[0];
      disp_bound_dof[i * 3 + 1] = node_num[1];
      disp_bound_dof[i * 3 + 2] = node_num[2];
      disp_bound_vals[i * 3] = du[0];
      disp_bound_vals[i * 3 + 1] = du[1];
      disp_bound_vals[i * 3 + 2] = du[2];
      disp_bound_init_vals[i*3] = disp[0]-du[0];
      disp_bound_init_vals[i*3+1] = disp[1]-du[1];
      disp_bound_init_vals[i*3+2] = disp[2]-du[2];
    }
  }
  void FiberRVEAnalysisExplicit::copyForceDataToForceVec()
  {
    auto ops = las::getLASOps<las::MICRO_BACKEND>();
    unsigned int num_dofs = getFn()->getDofCount();
    double * s = new double[num_dofs];
    amsi::WriteOp wrt;
    amsi::ToArray(getFn()->getFNumbering(), getFn()->getFField(), &s[0], 0, &wrt).run();
    ops->set(getF(), s);
  }
  void FiberRVEAnalysisExplicit::copyDispDataToDispVec()
  {
    auto ops = las::getLASOps<las::MICRO_BACKEND>();
    unsigned int num_dofs = getFn()->getDofCount();
    double * s = new double[num_dofs];
    amsi::WriteOp wrt;
    amsi::ToArray(getFn()->getUNumbering(), getFn()->getUField(), &s[0], 0, &wrt).run();
    ops->set(getU(), s);
  }
  void FiberRVEAnalysisExplicit::computeCauchyStress(double sigma[6])
  {
    copyForceDataToForceVec();
    FiberRVEAnalysis::computeCauchyStress(sigma);
  }
  void FiberRVEAnalysisExplicit::computeMaterialStiffness(double C[36])
  {
    //double diff_load_time = 2.0;
    //double diff_total_time = 2.5;
    ////double diff visc_damp_coeff = 0.6;
    //Amplitude * old_amp = amp;
    //double old_load_time = load_time;
    //double old_total_time = total_time;
    //// set values before computing finite difference
    //load_time = diff_load_time;
    //total_time = diff_total_time;
    //amp = new SmoothAmpHold(diff_load_time,diff_total_time);
    // compute material stiffness with default finite difference formulation
    RVEAnalysis::computeMaterialStiffness(C);
    // reset variables for the explicit class, so we don't "feel" the effects
    // from changing them when we do our full runs.
    //delete amp;
    //amp = old_amp;
    //total_time = old_total_time;
    //load_time = old_total_time;
  }
}  // namespace bio
