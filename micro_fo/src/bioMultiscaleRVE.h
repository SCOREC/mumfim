#ifndef BIO_MULTISCALE_RVE_H_
#define BIO_MULTISCALE_RVE_H_
#include "bioRVE2.h"
#include "bioMicroFOMultiscale.h"
namespace bio
{
  double calcRVEDimensionality(FiberNetwork * fn, double fbr_area, double fbr_vol_frc);
  class MultiscaleRVE
  {
  private:
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
    void dCidFE(apf::DynamicMatrix & dRVEdFE, const int ii, const apf::Vector3 & ci, const double rve_dim);
  public:
    MultiscaleRVE(RVE * rve,
                  FiberNetwork * fn,
                  micro_fo_header & hdr,
                  micro_fo_params & prms,
                  micro_fo_init_data & dat);
    ~MultiscaleRVE();
    void calcdRVEdFE(apf::DynamicMatrix & dRVEdFE, const RVE * rve);
  };
}
#endif
