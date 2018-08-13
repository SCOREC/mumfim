#ifndef BIO_MULTISCALE_RVE_H_
#define BIO_MULTISCALE_RVE_H_
#include "bioRVE.h"
#include "bioMicroFOMultiscale.h"
#include <amsiDetectOscillation.h>
namespace bio
{
  /*
   * \brief computes RVE scaling factor
   * Computes the scaling factor. This correlates to the RVE side length in "physical space"
   * \param fn a pointer to the fiber network data structure
   * \param fbr_area fiber cross-sectional area (from experiment)
   * \fbr_vol_frc volume fraction of the fiber network (from experiment)
   * \warning This function assumes that all of the fibers have the same cross-sectional area
   */
  double calcRVEDimensionality(FiberNetwork * fn, double fbr_area, double fbr_vol_frc);
  // coupling terms that involve terms from micro and macro
  class MultiscaleRVE
  {
  private:
    RVE * rve;
    int gss_id;
    int dim;
    apf::Vector3 lcl_gss;
    apf::Mesh * macro;
    apf::Field * macro_u;
    apf::MeshEntity * macro_ent;
    apf::MeshElement * macro_melmnt;
    apf::Element * macro_elmnt;
    int nnd; // num nodes effecting the macroscale element
    double fbr_area;
    double fbr_vl_frc;
    double rve_dim;
    double scale_conversion;
  public:
  MultiscaleRVE(RVE *rve, FiberNetwork *fn, micro_fo_header &hdr,
                micro_fo_params &prms, micro_fo_init_data &data);
  MultiscaleRVE(const MultiscaleRVE &mrve);
  RVE *getRVE() const { return rve; }
  ~MultiscaleRVE();
  void calcdRVEdFE(apf::DynamicMatrix &dRVEdFE);
  double getScaleConversion() const { return scale_conversion; }
  };
}
#endif
