#ifndef H_MSMMAPPING
#define H_MSMMAPPING

#include "Mapping.h"
#include <MeshSim.h>
#include <vector>

typedef TRELLIS_MTYPE MTYPE;

inline int  dimension(pMesh);
inline MTYPE element_type(pEntity pe); 
inline void get_vertices(pEntity ent, std::vector<pVertex> &verts); 

inline int* build_representation(pEntity); 
inline int  edge_dir(pRegion pf, int n);
inline PList* faces(pEdge pe);

// Old C-style names
inline TRELLIS_MTYPE M_GetElementType(pEntity pe)
{ return element_type(pe); }

inline void M_GetVertices(pEntity ent, std::vector<pVertex> &verts) 
{ return get_vertices(ent, verts); }

inline int* M_BuildRepresentation(pEntity e) 
{ return build_representation(e); }

inline int  R_edgeDir(pRegion pf, int n) 
{ return edge_dir(pf, n); }

inline int  M_dimension(pMesh m) 
{ return dimension(m); }

inline PList* E_faces(pEdge pe) 
{ return faces(pe); }

// Implementations ///////////////////

inline int dimension(pMesh pm) {
  if(M_numRegions(pm) != 0)
    return 3;
  else if(M_numFaces(pm) != 0)
    return 2;
  else if(M_numEdges(pm) != 0)
    return 1;
  else 
    throw 1;
  return -1;
}


inline MTYPE element_type(pEntity pE)
{
  switch(EN_type(pE)) {
  case 0: 
    return VERTEX;
  case 1: 
    return EDGE;
  case 2:
    {
      int numEdges = F_numEdges(static_cast<pFace>(pE));
      if(numEdges == 3)      return TRI;
      else if(numEdges == 4) return QUAD;
      else throw 1;	
    }
    break;
  case 3:
    {
      /// There is a difference here btwn Luo's and my MCTK --CWS 3/6/2008
      int numFaces = R_numFaces(static_cast<pRegion>(pE));
      if(numFaces == 4)      return TET;
      else if(numFaces == 6) return HEX;
      else if(numFaces == 5) return PRISM;
      else throw 1;	
    }
    break;
  default:
    throw 1;
  }
}

inline int edge_dir(pRegion pf, int n)
{
  const static int 
    TEdVt[6][2] = {{0,1},
		   {1,2},
		   {0,2},
		   {0,3},
		   {1,3},
		   {2,3}};

  int x;
  pPList edgelist = R_edges(pf, 1);
  pEdge edge = static_cast<pEdge>(PList_item(edgelist, n));
  pVertex v1 = E_vertex(edge, 0);
  pPList vertexlist = R_vertices(pf, 1);
  int indx = TEdVt[n][0];
  pVertex v1_bar = static_cast<pVertex>(PList_item(vertexlist, indx));
  if(v1 == v1_bar) x = 1;
  else x = -1;
  PList_delete(edgelist);
  PList_delete(vertexlist);
  return x;
}


inline void get_vertices(pEntity pE, std::vector<pVertex>& v)
{
  switch(EN_type(pE))
    {
    case 1:
      v.push_back(E_vertex((pEdge)pE,0));
      v.push_back(E_vertex((pEdge)pE,1));
      break;
    case 2:
      { 
	void *iter = 0;
	pEntity iadj; 
	pPList vtxs = F_vertices((pFace)pE,1); 
	while (iadj = static_cast<pEntity>(PList_next(vtxs,&iter)))
	  v.push_back(static_cast<pVertex>(iadj));
	PList_delete(vtxs);
      }
      break;
    case 3:      
      { 
	void *iter = 0;
	pEntity iadj; 
	pPList vtxs = R_vertices(static_cast<pRegion>(pE),1);
	while (iadj = static_cast<pEntity>(PList_next(vtxs,&iter)))
	  v.push_back(static_cast<pVertex>(iadj));
	PList_delete(vtxs);
      }
      break;
    default:
      throw 1;
    }
}

inline int* build_representation(pEntity pE)
{
  std::vector<pVertex> verts;
  get_vertices(pE, verts);
  int* x = new int[2 + verts.size()];
  x[0] = 2 + verts.size();
  x[1] = element_type(pE);
  for(size_t i = 0; i != verts.size(); ++i)
    x[2+i] = EN_id(static_cast<pEntity>(verts[i]));
  return x;
}

inline PList* faces(pEdge pe) {
  pPList pl = PList_new();
  int num = E_numFaces(pe);
  int i;
  pFace pf;
  for(i=0;i<num;i++) {
    pf = E_face(pe, i);
    pl = PList_appUnique(pl, pf);
  }
  if(PList_size(pl) != num)
    throw 1;

  return pl;
}

#endif
