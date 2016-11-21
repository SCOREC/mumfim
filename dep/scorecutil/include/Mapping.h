#ifndef H_Mapping
#define H_Mapping

#include <vector>
#include <list>

#include <MeshSim.h>

#include "mTensor2.h"
#include "mPoint.h"

enum Trellis_mType {VERTEX,EDGE,TRI,QUAD,HEX,PRISM,PYRAMID,TET};
typedef Trellis_mType TRELLIS_MTYPE;

/**
   Mapping for mesh entities. We suppose that we have 2 systems of
   coordinates : (u,v,w) and (x,y,z). The mapping deals with these
   2 systems of coordinates.
 */
  typedef double (*LevelSetFunction)(mPoint &);

  class mVector;
  class mPoint;
  class Mapping;

  class Mapping
  {
    protected:

/******************************************************************
  Three new member variables and corresponding member functiones 
  are added by P.Hu  July, 2002
******************************************************************/
    TRELLIS_MTYPE entityType;
    int entityDim;
    Mapping* rootMapping; 

    pEntity ent;
    
    std::vector<mPoint> knots;
    int m_order; // the geometric approximation shape order of the mesh entity
                 // m_order = 1, linear; m_order = 2, quadratic

    
    public:
      /**
	 Constructor takes a mesh entity as input	 
       */
    Mapping(pEntity e);
      /**
	 Constructor takes a mesh entity, vector of all vertecies and entity type, entity dim  as input	 
       */
/****************************************************************
  One Question: Member variable 'knots' is not in the base class
  Mapping, therefore, this constructor can't really work here, if 
  move the variable 'knots' in the base calss, it will be better
****************************************************************/

      Mapping(Mapping *theMapping, 
	      std::vector<mPoint> *_input_knots, 
	      TRELLIS_MTYPE _type, int dim);
      
      virtual void buildMappings(Mapping * theMapping, LevelSetFunction _levelSetFunc, 
				 int flag, std::list<Mapping *> &out_list );
      
      /**
	 virtual destructor does nothing
       */
      virtual ~Mapping(){}
      /**
	 get the mesh entity
       */
      inline pEntity getEntity () const {return ent;}
      /**
	 set the entity type
      */
      inline void setEntityType(TRELLIS_MTYPE _type) {entityType = _type;}
      /**
	 get the entity type
      */
      inline TRELLIS_MTYPE getEntityType() const {return entityType;}
      /**
	 set the entity dimension
      */
      inline void setEntityDim(int _dim) {entityDim = _dim;}
      /**
	 get the entity dimension
      */
      inline int getEntityDim() const {return entityDim;}
      /**      
	 set the root mapping
      */
      inline void setRootMapping(Mapping *_theMapping) {rootMapping = _theMapping;}
      /**
	 get the root mapping
      */
      inline Mapping * getRootMapping() const {return rootMapping;}

      void local2rootlocal(double u, double v, double w, double &r_u, double &r_v, double &r_w);

      /**
	 Inversion of the mapping by a newton algorithm
       */      
      virtual bool invert(double x, double y, double z,
			  double &u, double &v, double &w) const;
      /**
	 Jacobian = d(x,y,z)/d(u,v,w) is a 3x3 matrix
	 computation of the inverse of the jacobian at (u,v,w), 
	 returns the determinant of the jacobian (not of the inverse).
       */
      virtual double jacInverse(double x, double y, double z, 
				mTensor2&) const;
      /**
	 returns determinant of jacobian at (u,v,w).
       */
      virtual double detJac(double u, double v, double w) const;
      /**
	 tells if a point (u,v,w) is inside the reference element or not
       */
      virtual bool inReferenceElement(double u, double v, double w) const;
      /**
	 checks if a point in REAL coordinates is inside the element ent
	 if it's the case, returns local coordinates
       */
      virtual bool interiorCheck (pEntity , const mPoint &p, double &u, double &v, double &w) const;   
      /**
	 returns the center of gravity of the element in local coordinates
       */
      virtual void COG (double &u, double &v, double &w) const;
      /**
	 Computes gradients in real world knowing them in refernce world
       */
      virtual double PushBack(double u, double v, double w, int vsize, 
			      mVector *vec ) const;
      virtual double PushBack( double u, double v, double w, int vsize, 
			       mTensor2 *vec ) const;		       

      inline double PushBack (double u, double v, double w, int vsize, 
			       std::vector<mVector> &vec) const
      {
		return PushBack (u,v,w,vsize,&(*vec.begin()));
      }      
      
      int simplexation ( int &dim , int &nbPnt, double *u, double *v, double *w, int combn [][4] ) const ; 


      //---------------------------------------------//
      //  P u r e   V i r t u a l   M e m b e r s    //
      //---------------------------------------------//


      /**
	 eval and deval functions have to be re-written for different mesh
	 mappings (curvilinear, bezier,...).	 
	 eval computes x(u,v,w) , y(u,v,w) and z(u,v,w) 
	 deval computes derivatives dx/du , dy/du ...
      */

      
      virtual void eval (double u, double v, double w,
			 double &x, double &y, double &z) const = 0;

      virtual void deval(double u, double v, double w,
			 double &dxdu, double &dydu, double &dzdu,
			 double &dxdv, double &dydv, double &dzdv,
			 double &dxdw, double &dydw, double &dzdw) const = 0;

      /**
	 compute the exterior normal to ent at a point u,v,w of the reference
	 face face
       */
      virtual void normalVector(pEntity face, double u , double v , double w, mVector &n) const = 0;
      
      /**
	 computes the normal vector to a face !!! If it's an edge, then
	 we suppose that the normal vector is t_edge x ez 
       */      
      virtual void normalVector(double u , double v , double w, mVector &n) const;

      /**
	 computes the bounding box of a mesh entity
       */
      virtual void boundingBox (  mPoint &min, mPoint &max ) const = 0;

      virtual int order() const = 0;
      virtual int geomOrder() const = 0;

  };

   // Added for compilation with FMDB --CWS 3/6/2008
  inline void  get_edges(pEntity pE, std::vector<pEdge> &v)
    {
      switch(EN_type(pE))
        {
        case 1:
          {
        if(E_numPoints((pEdge)pE))
          v.push_back((pEdge)pE);
          }
          break;
        case 2:
          {
            void *iter = 0;
            pEntity iadj;
            pPList fedgs = F_vertices((pFace)pE, 1);
            pPList edgs = F_edges((pFace)pE,1, (pVertex)PList_item(fedgs, 0));
            while ((iadj = (pEntity)PList_next(edgs,&iter))) {
              if(E_numPoints((pEdge)iadj))
                v.push_back((pEdge)iadj);
        }
            PList_delete(edgs);
            PList_delete(fedgs);
          }
          break;
        case 3:
          {
            void *iter = 0;
            pEntity iadj;
            pPList edgs = R_edges((pRegion)pE,1);
            while ((iadj = (pEntity)PList_next(edgs,&iter))){
          if(E_numPoints((pEdge)iadj))
            v.push_back((pEdge)iadj);
            }
            PList_delete(edgs);
          }
          break;
        default:
          throw 1;
        }
    }
  inline void M_GetEdges(pEntity ent, std::vector<pEdge> &e)
    { return get_edges(ent, e); }

#endif
