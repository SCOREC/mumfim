#ifndef BIO_RVE_H_
#define BIO_RVE_H_

#include "MicroFOUtil.h"

#include <apf.h>
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
    
  protected:
    apf::Element * getElement() { return cbe_u_e; }
  public:
    RVE();
    
    /**
     * Get the number of nodes for the RVE
     * @return The number of nodes on the RVE (4 for 2d, 8 for 3d, -1 for failure)
     */
    int numNodes() const { return dim == 2 ? 4 : dim == 3 ? 8 : -1; }

    double sideCoord(int sd)
    {
      assert(sd >= 0 && sd < 6);
      static double op[6] = {-1.0,-1.0,1.0,1.0,-1.0,1.0};
      return hd * op[sd];
    }

    bool onBoundary(const apf::Vector3 & crd)
    {
      return close(sideCoord(2),crd[0]) || close(sideCoord(4),crd[0])
          || close(sideCoord(0),crd[1]) || close(sideCoord(5),crd[1])
          || close(sideCoord(1),crd[2]) || close(sideCoord(3),crd[2]);
    }
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

    // todo (h) : fix this
    void atNode(int nde)
    {
      apf::Vector3 p;
      apf::Vector3 u;
      apf::getVector(crd_fld,cur_ent,nde,p); // get xyz coord of fiber node
      //p = calcLocalCoord(src_elmnt->getMesh(),apf::getMeshEntity(src_elmnt));
      
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
   * @note This is specific to the scale-coupling and should maybe go in MacroCoupling.h
   */
  void calcGlobalRVECoords(apf::DynamicArray<apf::Vector3> & rve_crds,
			   double rve_dim,
			   const apf::Vector3 & gbl_gss);

  void calcBoundaryNodes(const RVE * rve,
			 const FiberNetwork * fn,
			 std::vector<apf::MeshEntity*> & bnds);
}

#endif
