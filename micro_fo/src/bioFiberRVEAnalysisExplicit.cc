#include "bioFiberRVEAnalysisExplicit.h"
#include <PCU.h>
#include "bioExplicitAmplitude.h"
#include <apf.h>
#include <apfConvert.h>
#include <apfMesh2.h>
#include "bioFiberNetworkIO.h"
#include "bioMassIntegrator.h"
#include <lionPrint.h>
#include <mpi.h>
#include "bioExplicitOutputWriter.h"
#include "bioFiberRVEAnalysisExplicit_impl.h"
#include <Kokkos_Core.hpp>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <las.h>
namespace bio
{
  bool FiberRVEAnalysisExplicit::run(const DeformationGradient & dfmGrd)
  {
    computeDisplcamentBC(dfmGrd);
    if (this->getFn()->getDofCount() < serial_gpu_cutoff)
    {
      // run serial analysis
      ExplicitAnalysisSerial analysis(
          static_cast<apf::Mesh2 *>(this->getFn()->getNetworkMesh()),
          total_time, fiber_elastic_modulus, fiber_area, fiber_density, amp,
          analysis_name, visc_damp_coeff, print_history_frequency,
          print_field_frequency, print_field_by_num_frames,
          crit_time_scale_factor, energy_check_eps, getFn()->getUField(),
          getFn()->getVField(), getFn()->getAField(), getFn()->getFField());
      analysis.setDispBC(disp_bound_nfixed, disp_bound_dof, disp_bound_vals);
      return analysis.run();
    }
    else
    {
      // run gpu analysis
      ExplicitAnalysisKokkos analysis(
          static_cast<apf::Mesh2 *>(this->getFn()->getNetworkMesh()),
          total_time, fiber_elastic_modulus, fiber_area, fiber_density, amp,
          analysis_name, visc_damp_coeff, print_history_frequency,
          print_field_frequency, print_field_by_num_frames,
          crit_time_scale_factor, energy_check_eps, getFn()->getUField(),
          getFn()->getVField(), getFn()->getAField(), getFn()->getFField());
      analysis.setDispBC(disp_bound_nfixed, disp_bound_dof, disp_bound_vals);
      return analysis.run();
    }
  }
  void FiberRVEAnalysisExplicit::computeDisplcamentBC(
      const DeformationGradient & dfmGrd)
  {
    disp_bound_nfixed = 3 * bnd_nds[RVE::all].size();
    disp_bound_dof = new int[disp_bound_nfixed];
    disp_bound_vals = new double[disp_bound_nfixed];
    apf::Field * current_coords_field = getFn()->getXpUField();
    apf::Vector3 coords;
    int node_num[3];
    for (std::size_t i = 0; i < bnd_nds[RVE::all].size(); ++i)
    {
      apf::MeshEntity * nd = bnd_nds[RVE::all][i];
      apf::getVector(current_coords_field, nd, 0, coords);
      node_num[0] = apf::getNumber(getFn()->getUNumbering(), nd, 0, 0);
      node_num[1] = apf::getNumber(getFn()->getUNumbering(), nd, 0, 1);
      node_num[2] = apf::getNumber(getFn()->getUNumbering(), nd, 0, 2);
      disp_bound_dof[i * 3] = node_num[0];
      disp_bound_dof[i * 3 + 1] = node_num[1];
      disp_bound_dof[i * 3 + 2] = node_num[2];
      disp_bound_vals[i * 3] = (dfmGrd[0] - 1) * coords[0] +
                               dfmGrd[1] * coords[1] + dfmGrd[2] * coords[2];
      disp_bound_vals[i * 3 + 1] = dfmGrd[3] * coords[0] +
                                   (dfmGrd[4] - 1) * coords[1] +
                                   dfmGrd[5] * coords[2];
      disp_bound_vals[i * 3 + 2] = dfmGrd[6] * coords[0] +
                                   dfmGrd[7] * coords[1] +
                                   (dfmGrd[8] - 1) * coords[2];
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
    apf::Integrator * es = createImplicitMicroElementalSystem(getFn(), getK(),getF());
    es->process(getFn()->getNetworkMesh(), 1);
    delete es;
  }
}  // namespace bio
