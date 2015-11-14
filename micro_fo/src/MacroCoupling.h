#ifndef BIO_MACRO_COUPLING_H
#define BIO_MACRO_COUPLING_H

#include "apfUtil.h"
#include "RVE.h"

#include <apf.h>
#include <apfDynamicMatrix.h>

namespace bio
{
  /**
   * A class to manage the multi-scale coupling information relating a single fiber network
   *  to a single macro-scale integration point and mesh element.
   */
  class MacroInfo
  {
  protected:
    int gss_id;
    int dim;

    apf::Vector3 lcl_gss;
    apf::Mesh * macro_msh;
    apf::MeshEntity * macro_ent;
    apf::MeshElement * macro_melmnt;
    apf::Element * macro_elmnt;
    int nnd; // num nodes effecting the macro element

    double fbr_area;
    double fbr_vl_frc;
    double rve_dim;
    
    void dCidFE(apf::DynamicMatrix&,const int,const apf::Vector3 &, double);
  public:
    MacroInfo();
    
    /**
     * Calculate the term relating the macro-scale nodal displacements to the
     *  micro-scale cube corner displacements.
     * @param drve_dfe The \f$ N_{dofs_{rve}} \times N_{dofs_{fe}} \f$ matrix approximating
     *                 the displacement caused on the RVE nodes due to the displacement
     *                 of the macro-scale element nodes.
     * @param rve The RVE effected by the displacement.
     */
    void calcdRVEdFE(apf::DynamicMatrix & drve_dfe, const RVE * rve);
  };

  /**
   * Calculate the macro-scale dimension of the RVE cube (xyz).
   * @param fn The FiberNetwork to determine the dimensionality of.
   * @param fbr_area The cross-sectional area of a single fiber in the macro-scale coordinate system.
   * @param fbr_vl_frc The fiber volume fraction (% fibers per unit volume)
   * @return The dimension of the RVE (in the macro-scale coordinate system).
   */
  double calcRVEDimensionality(const FiberNetwork * fn,double fbr_area, double fbr_vl_frc);
}

#endif
