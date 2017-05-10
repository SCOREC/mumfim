#ifndef BIO_RVE2_H_
#define BIO_RVE2_H_
//#include "MicroFOUtil.h"
#include "bioUtil.h"
#include <amsiLAS2.h>
#include <apfFieldOp.h>
#include <apf.h>
#include <apfDynamicVector.h>
#include <apfMesh.h>
#include <cassert>
namespace bio
{
  class FiberNetwork;
  class RVE
  {
  public:
    enum side // [ -y -z x z -x y]
    {
      bot = 0,
      frt = 1,
      rgt = 2,
      bck = 3,
      lft = 4,
      top = 5,
      all = 6
    };
  protected:
    int dim;
    double crd;
    apf::Mesh * cbe;
    apf::MeshEntity * cbe_e;
    apf::Element * cbe_u_e;
    apf::Field * cbe_u;
    apf::Numbering * cbe_dof;
  protected:
  public:
    RVE(int d = 3);
    /**
     * Get the number of nodes for the RVE
     * @return The number of nodes on the RVE (4 for 2d, 8 for 3d, -1 for failure)
     */
    int numNodes() const { return dim == 2 ? 4 : dim == 3 ? 8 : -1; }
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
    apf::MeshEntity * getMeshEnt() { return cbe_e; }
    apf::Element * getElement() { return cbe_u_e; }
    apf::Numbering * getNumbering() { return cbe_dof; }
    apf::Field * getField() { return cbe_u; }
  };
  /**
   * A FieldOp designed to take the shape function values from an 'enclosing' element
   *  and apply them as displacements on 'contained' mesh elements. Used by the MacroCoupling
   *  to interpolate the displacements on the cube contatining the rve fiber network onto
   *  the individual nodes in the fiber network.
   */
  class InterpOnto : public amsi::FieldOp
  {
  protected:
    apf::Element * src_elmnt;
    apf::MeshEntity * cur_ent;
    apf::Field * dest_fld;
    apf::Field * crd_fld;
  public:
    InterpOnto(apf::Element * se, apf::Field * d) :
      src_elmnt(se),
      cur_ent(NULL),
      dest_fld(d),
      crd_fld(NULL)
    {
      crd_fld = apf::getMesh(dest_fld)->getCoordinateField();
    }
    bool inEntity(apf::MeshEntity * me)
    {
      cur_ent = me;
      return true;
    }
    void outEntity() {}
    // todo (h) : fix this (this might be fixed... i dunno anymore, lol)
    void atNode(int nde)
    {
      apf::Vector3 p;
      apf::Vector3 u;
      apf::getVector(crd_fld,cur_ent,nde,p); // get xyz coord of fiber node
      //mapGlobalToLocal(apf::getMesh(crd_fld),apf::getMeshEntity(src_elmnt),p,u);
      apf::getVector(src_elmnt,p,u);
      apf::setVector(dest_fld,cur_ent,nde,u);
    }
  };
  /**
   * Use the displacement of the RVE to displace the fiber network nodes
   *  'inside' of the RVE
   */
  void forwardRVEDisplacement(RVE * rve, FiberNetwork * fn);
  /**
   * Calculate the position of the corner nodes of the square/cube enclosing the RVE
   *  in the cartesian space of the problem.
   * @param rve_crds The coordinates of the 8 vertices of the RVE (using the standard
   *                 finite element ordering) in the macro-scale coordinate system.
   * @param rve_dim The dimensionality of the RVE in the macro-scale coordinate
   *                system, using calculated by calling calcRVEDimensionality().
   * @param gbl_gss The macro-scale coordinate at the center of the RVE
   */
  void calcGlobalRVECoords(apf::DynamicArray<apf::Vector3> & rve_crds,
                           double rve_dim,
                           const apf::Vector3 & gbl_gss);
  /**
   * Compile a list of MeshEntities in the fiber network which have at least a single
   * node classified 'on' the boundary of the RVE.
   */
  template <typename O>
    void getBoundaryVerts(const RVE * rve, const FiberNetwork * fn, RVE::side sd, O nds);
  /**
   * Apply a vector of displacements to the coordinates of the rve mesh
   * @param rve The rve to displace
   * @param du A vector of the correct length (24 for 3d, 8 for 2d) containing
   *           displacement values for each rve node in the typically finite element
   *           ordering for a cube
   */
  void displaceRVE(RVE * rve, const apf::DynamicVector & du);
  //
  void alignFiberNetwork(RVE * rve, FiberNetwork * fn, const double align_vec[3]);
  void affineDeformation(RVE * rve, FiberNetwork * fn, const double disp[6]);
  void updateRVEBounds(RVE * rve, FiberNetwork * fn, const double disp[6]);
  double calcFiberDensity(RVE * rve,FiberNetwork * fn);
}
#include <bioRVE2_impl.h>
#endif
