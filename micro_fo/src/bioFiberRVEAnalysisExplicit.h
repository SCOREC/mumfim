#ifndef BIO_FIBER_RVE_ANALYSIS_EXPLICIT
#define BIO_FIBER_RVE_ANALYSIS_EXPLICIT
#include "bioExplicitAmplitude.h"
#include "bioFiberNetwork.h"
#include "bioFiberRVEAnalysis.h"
#include <string>
namespace bio
{
  // TODO FIXME deal with network with multiple material properties
  class FiberRVEAnalysisExplicit : public FiberRVEAnalysis
  {
    public:
    FiberRVEAnalysisExplicit(FiberNetwork * fn,
                             LinearStructs<las::MICRO_BACKEND> * vecs,
                             const MicroSolutionStrategyExplicit & ss)
        : FiberRVEAnalysis(fn, vecs, static_cast<MicroSolutionStrategy>(ss))
        , visc_damp_coeff(ss.visc_damp_coeff)
        , print_history_frequency(ss.print_history_frequency)
        , print_field_frequency(ss.print_field_frequency)
        , print_field_by_num_frames(ss.print_field_by_num_frames)
        , total_time(ss.total_time)
        , serial_gpu_cutoff(ss.serial_gpu_cutoff)
        , crit_time_scale_factor(ss.crit_time_scale_factor)
        , energy_check_eps(ss.energy_check_eps)
        , disp_bound_nfixed(0)
        , disp_bound_nodes(NULL)
        , disp_bound_dof(NULL)
        , disp_bound_vals(NULL)
        , fiber_elastic_modulus(fn->getFiberReactions()[0]->E)
        , fiber_area(fn->getFiberReactions()[0]->fiber_area)
        , fiber_density(fn->getFiberReactions()[0]->fiber_density)
        , amp(NULL)
    {
      std::stringstream sout;
      int rnk = -1;
      MPI_Comm_rank(MPI_COMM_WORLD, &rnk);
      sout << "rnk_" << rnk << "_fn_" << getFn()->getRVEType()<<"_explicit";
      analysis_name = sout.str();
      if(ss.ampType == AmplitudeType::SmoothStep)
      {
        amp = new SmoothAmp(total_time);
      }
      else
      {
        std::abort();
      }
    }
    ~FiberRVEAnalysisExplicit()
    {
      delete amp;
      delete [] disp_bound_nodes;
      delete [] disp_bound_dof;
      delete [] disp_bound_vals;
    }
    virtual bool run(const DeformationGradient & dfmGrd);
    virtual FiberRVEAnalysisType getAnalysisType()
    {
      return FiberRVEAnalysisType::Explicit;
    }
    protected:
    // This parameters gives the number of fibers under which the serial
    // implementation will be used.
    int serial_gpu_cutoff;
    double total_time;
    double fiber_elastic_modulus;
    double fiber_area;
    double fiber_density;
    Amplitude * amp;
    std::string analysis_name;
    double visc_damp_coeff;
    unsigned long print_history_frequency;
    unsigned long print_field_frequency;
    bool print_field_by_num_frames;
    double crit_time_scale_factor;
    double energy_check_eps;
    int disp_bound_nfixed;
    int * disp_bound_nodes;
    int * disp_bound_dof;
    double * disp_bound_vals;
    void computeDisplcamentBC(const DeformationGradient & dfmGrd);
  };
}  // namespace bio
#endif
