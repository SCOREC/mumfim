#ifndef BIO_RVE_H_
#define BIO_RVE_H_
#include "bioUtil.h"
#include <apfFunctions.h> //amsi
#include <apfFieldOp.h>
#include <apf.h>
#include <apfDynamicVector.h>
#include <apfMesh.h>
#include <cassert>
namespace bio
{
  class FiberNetwork;
  // TODO : pull calculation of dV_dx_rve in here since the term is derived entirely from RVE values
  class RVE
  {
  protected:
    int dim;
    double crd;
    apf::Mesh * cbe;
    apf::MeshEntity * cbe_e;
    apf::Element * cbe_u_e;
    apf::Field * cbe_u;
    apf::Field * cbe_du;
    amsi::XpYFunc * xpufnc;
    apf::Field * cbe_xpu;
    apf::Numbering * cbe_dof;
    int cbe_dof_cnt;
  public:
    RVE(double cr = 0.5, int d = 3);
    RVE(const RVE & rve);
    ~RVE();
    /// This is an enum of the sides of the RVE
    enum side // [ -y -z x z -x y]
    {
      bot = 0,  ///< bottom face
      frt = 1,  ///< front face
      rgt = 2,  ///< right face
      bck = 3,  ///< back face
      lft = 4,  ///< left face
      top = 5,  ///< top face
      all = 6   ///< all faces
    };
    /*
     * Returns all mesh entities on a specific side.
     */
    apf::MeshEntity * getSide(side sd) const;
    /**
     * Get the number of nodes for the RVE
     * @return The number of nodes on the RVE (4 for 2d, 8 for 3d, -1 for failure)
     */
    int numNodes() const { return dim == 2 ? 4 : dim == 3 ? 8 : -1; }
    /*
     * Return the number of degrees of freedom in the rve cube
     */
    int getDofCount() const { return cbe_dof_cnt; }
    /*
     * Get the dimensionality of the RVE
     */
    int getDim() const {return dim;}
    /*
     * Measure the RVE based on dimensionality
     *  gives area for a 2d RVE and volume for a 3d
     *  rve.
     * @note This is the displaced RVE measure,
     *  not the reference measure.
     */
    double measureDu() const
    {
      apf::MeshElement * mlmt = apf::createMeshElement(cbe_xpu,cbe_e);
      double msr = apf::measure(mlmt);
      apf::destroyMeshElement(mlmt);
      return msr;
    }
    /*
     * Measure the RVE based on dimensionality
     *  gives area for a 2d RVE and volume for a 3d
     *  rve.
     * @note This is the reference RVE measure,
     *  not the displaced measure.
     */
    double measure() const
    {
      apf::MeshElement * mlmt = apf::createMeshElement(cbe,cbe_e);
      double msr = apf::measure(mlmt);
      apf::destroyMeshElement(mlmt);
      return msr;
    }
    /**
     * Retrieve the coordinate related to a side of the RVE. This operates on the
     *  reference configuration of the RVE, so all faces of the cube are axis-aligned
     *  and only a single coordinate is of importance.
     * @param sd The side for which to get the primary coordinate, ordered as [-y -z x z -x y]
     *           in keeping with the standard ordering of faces on a hexahedral element.
     * @return The primary coordinate for the side specified by the parameter.
     */
    double sideCoord(side sd) const
    {
      static const double op[6] = {-1.0,-1.0,1.0,1.0,-1.0,1.0};
      return crd * op[sd];
    }
    /**
     * Determine whether a given coordinate in the dimensionless RVE space is on
     * (within double-precision \f$\epsilon\f$ of) the boundary.
     * @note This works on the undeformed RVE, so any displacements of the RVE
     *       cube will not effect the result of this function
     * @param crd The coordinate in the RVE space to check
     * @return Whether any coordinate of the parameter is within \f$ \epsilon \f$
     *         of an initial boundary of the RVE.
     */
    bool onBoundary(const apf::Vector3 & crd, side sd, double eps = 1e-8) const
    {
      if(sd == side::all)
        return fabs(sideCoord(side::rgt) - crd[0]) < eps ||
               fabs(sideCoord(side::lft) - crd[0]) < eps ||
               fabs(sideCoord(side::bot) - crd[1]) < eps ||
               fabs(sideCoord(side::top) - crd[1]) < eps ||
               fabs(sideCoord(side::frt) - crd[2]) < eps ||
               fabs(sideCoord(side::bck) - crd[2]) < eps;
      int idx = -1;
      if(sd == side::lft || sd == side::rgt)
        idx = 0;
      else if (sd == side::bot || sd == side::top)
        idx = 1;
      else if (sd == side::frt || sd == side::bck)
        idx = 2;
      return fabs(sideCoord(sd) - crd[idx]) < eps;
     }
    apf::Mesh * getMesh() const { return cbe; }
    apf::MeshEntity * getMeshEnt() const { return cbe_e; }
    apf::Element * getElement()
    {
      if(cbe_u_e)
        apf::destroyElement(cbe_u_e);
      cbe_u_e = apf::createElement(cbe_u,cbe_e);
      return cbe_u_e;
    }
    apf::Numbering * getNumbering() const { return cbe_dof; }
    apf::Field * getUField() const { return cbe_u; }
    apf::Field * getdUField() const { return cbe_du; }
    apf::Field * getXpUField() const { return cbe_xpu; }
  };
  /**
   * Compile a list of MeshEntities in the fiber network which have at least a single
   * node classified 'on' the boundary of the RVE.
   */
  template <typename I, typename O>
    void getBoundaryVerts(const RVE * rve, apf::Mesh * msh,
                          I bgn_vrts, I end_vrts,
                          RVE::side sd, O nds);
  /**
   * Apply a vector of displacements to the coordinates of the rve mesh
   * @param rve The rve to displace
   * @param du A vector of the correct length (24 for 3d, 8 for 2d) containing
   *           displacement values for each rve node in the typically finite element
   *           ordering for a cube
   */
  void displaceRVE(RVE * rve, const apf::DynamicVector & du);
  void getRVEDisplacement(RVE * rve, apf::DynamicVector & u);
  void getRVEReferenceCoords(RVE * rve, apf::DynamicVector & xyz_0);
  void getRVECoords(RVE * rve, apf::DynamicVector & xyz);
  //
  void updateRVEBounds(RVE * rve, FiberNetwork * fn, const double disp[6]);
}
#include <bioRVE_impl.h>
#endif
