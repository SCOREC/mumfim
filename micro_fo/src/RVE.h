#ifndef BIO_RVE_H_
#define BIO_RVE_H_
#include "apfUtil.h"
#include "MicroFOUtil.h"
#include "lasSparskit.h"
#include <apf.h>
#include <apfDynamicVector.h>
#include <apfElement.h>
#include <apfField.h>
#include <apfMesh.h>
#include <cassert>
namespace bio
{
  class FiberNetwork;
  class RVE
  {
  protected:
    int dim;
    double hd;
    apf::Mesh * cbe;
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
    double sideCoord(int sd) const
    {
      assert(sd >= 0 && sd < 6);
      static double op[6] = {-1.0,-1.0,1.0,1.0,-1.0,1.0};
      return hd * op[sd];
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
    bool onBoundary(const apf::Vector3 & crd) const
    {
      return close(sideCoord(2),crd[0]) || close(sideCoord(4),crd[0])
          || close(sideCoord(0),crd[1]) || close(sideCoord(5),crd[1])
          || close(sideCoord(1),crd[2]) || close(sideCoord(3),crd[2]);
    }
    apf::Element * getElement() { return cbe_u_e;}
    apf::Numbering * getNumbering() { return cbe_dof; }
    apf::Field * getField() { return cbe_u; }
  };
  /**
   * A FieldOp designed to take the shape function values from an 'enclosing' element
   *  and apply them as displacements on 'contained' mesh elements. Used by the MacroCoupling
   *  to interpolate the displacements on the cube contatining the rve fiber network onto
   *  the individual nodes in the fiber network.
   */
  class InterpOnto : public apf::FieldOp
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
    }
    void outEntity() {}

    // todo (h) : fix this (this might be fixed... i dunno anymore, lol)
    void atNode(int nde)
    {
      apf::Vector3 p;
      apf::Vector3 u;
      apf::getVector(crd_fld,cur_ent,nde,p); // get xyz coord of fiber node
      calcLocalCoord(src_elmnt->getMesh(),apf::getMeshEntity(src_elmnt),p);
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
   * @param rve The RVE to check against.
   * @param fn The FiberNetwork to check all nodes for being of the RVE boundary.
   * @param bnds The MeshEntities determined to lie on the boundary of the RVE.
   * @note Only checks against the reference configuration of the RVE
   */
  void calcBoundaryNodes(const RVE * rve,
			 const FiberNetwork * fn,
			 std::vector<apf::MeshEntity*> & bnds);
  /**
   * Set the force vector value associated with dofs on nodes lying on the RVE
   *  boundaries to zero.
   * @param f The force vector
   * @param rve The rve which defines which nodes are on the boundary
   * @param fn The fiber network to check for boundary nodes
   */
  void applyRVEForceBC(las::skVec * f,
		                   RVE * rve,
		                   FiberNetwork * fn);
  /**
   * Apply a vector of displacements to the coordinates of the rve mesh
   * @param rve The rve to displace
   * @param du A vector of the correct length (24 for 3d, 8 for 2d) containing
   *           displacement values for each rve node in the typically finite element
   *           ordering for a cube
   */
  void displaceRVE(RVE * rve,
                   const apf::DynamicVector & du);
}
#endif
