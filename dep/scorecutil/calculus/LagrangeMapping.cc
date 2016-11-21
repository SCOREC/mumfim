/* 
   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of Trellis written and maintained by the 
   Scientific Computation Research Center (SCOREC) at Rensselaer Polytechnic
   Intitute, Troy, NY, USA.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the Rensselaer SCOREC Public License.

   This program is distributed in the hope that it will be useful, 
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.
   
   You should have received a copy of the Rensselaer SCOREC Public License
   along with this program; if not, write to Rensselaer Polytechnic Institure,
   110 8th Street, SCOREC, Troy, NY  12180, USA
*/
#include <iostream>
#include <algorithm>
#include "LagrangeMapping.h"
#ifndef SIM
#include "AOMD_Internals.h"
#endif
#include "MSMapping.h"
#include "mPoint.h"

#include <cmath>
#include <cstdio>
#include <cassert>

  using std::cout;
  using std::vector;
  using std::abs;

  LagrangeMapping::LagrangeMapping (pEntity e)
    :Mapping(e)
  {
    vector<pVertex> verts;
    M_GetVertices (ent,verts);
    for(size_t i=0; i < verts.size(); ++i)
      {
	pPoint p =  V_point(verts[i]); 
	knots.push_back(mPoint(P_x(p),P_y(p),P_z(p)));
      }    
    entityType = M_GetElementType(ent);
    m_order = 1;

    //edges
    std::vector<pEdge> edges;
    M_GetEdges(ent, edges);
    if(edges.size() > 0)
      m_order = 2;
    for(size_t i=0; i<edges.size(); i++)
      {
	if(E_numPoints(edges[i])) {
	  pPoint p = E_point(edges[i], 0);
	  knots.push_back(mPoint(P_x(p),P_y(p),P_z(p)));
	}
      }
  }

  LagrangeMapping::LagrangeMapping(pEntity me, vector<mPoint> * input_knots) :
    Mapping(me)
  {
    knots.reserve(input_knots->size());
    size_t i = 0;
    while ( i < input_knots->size())
      {
	knots.push_back((*input_knots)[i]);
	i++;
      }

    if(me)
      entityType = M_GetElementType(me);
    else {
      // in case the mapping is independent of mesh entity
      //Luo's -- CWS 3/11/2008
      m_order = 1;
      switch(input_knots->size()) {
      case 8:
        entityType = HEX;
        entityDim = 3;
        break;
      case 4:
        entityType = TET;
        entityDim = 3;
        break;
      case 6:
        entityType = PRISM;
        entityDim = 3;
        break;
      case 5:
        entityType = PYRAMID;
        entityDim = 3;
        break;
      default:
        break;
      }
    }

    // determine the appropriate m_order
    switch(entityType) 
      {
      case VERTEX:
	break;
	
      case EDGE:
	if(i==3)
	  m_order = 2;
	break;

      case TRI:
	if(i==6)
	  m_order = 2;
	break;

      case QUAD:
	if(i==8)
	  m_order = 2;
	break;

      case TET:
	if(i==10)
	  m_order = 2;
	break;

      case HEX:
	if(i==20)
	  m_order = 2;
	break;

      case PRISM:
	if(i==15)
	  m_order = 2;
	break;

      case PYRAMID:
	if(i==13)
	  m_order = 2;
	break;
      }
    
  }

  
  int LagrangeMapping::order() const
  { 
    return m_order; 
  }

  int LagrangeMapping::geomOrder() const
  { 
    return m_order; 
  }


  void LagrangeMapping::GeomShapeFunctionQuad (double u, double v, double dum, double *x) const
  {
    if(m_order == 1) {
      x[0] =  0.25 * (1.-u) * (1.-v);
      x[1] =  0.25 * (1.+u) * (1.-v);
      x[2] =  0.25 * (1.+u) * (1.+v);
      x[3] =  0.25 * (1.-u) * (1.+v);
    }
    else if(m_order == 2) {
      x[0] = (1.0-u)*(1.0-v)/4.0-(1.0-u*u)*(1.0-v)/4.0-(1.0-u)*(1.0-v*v)/4.0;  
      x[1] = (1.0+u)*(1.0-v)/4.0-(1.0-u*u)*(1.0-v)/4.0-(1.0+u)*(1.0-v*v)/4.0;
      x[2] = (1.0+u)*(1.0+v)/4.0-(1.0+u)*(1.0-v*v)/4.0-(1.0-u*u)*(1.0+v)/4.0;
      x[3] = (1.0-u)*(1.0+v)/4.0-(1.0-u*u)*(1.0+v)/4.0-(1.0-u)*(1.0-v*v)/4.0;
      x[4] = (1.0-u*u)*(1.0-v)/2.0;
      x[5] = (1.0+u)*(1.0-v*v)/2.0;
      x[6] = (1.0-u*u)*(1.0+v)/2.0;
      x[7] = (1.0-u)*(1.0-v*v)/2.0;
     
    }
    else {
      std::cout<<"The lagrange mapping does not support cubic Quad"<<std::endl;
    }

    return;
  } 

  void LagrangeMapping::GeomShapeFunctionTet (double r, double s, double t, double *x) const
  {
    double u = 1.-r-s-t;
    
    if(m_order == 1) {
      x[0] =  u;
      x[1] =  r;
      x[2] =  s;
      x[3] =  t;
    }
    else if(m_order == 2) {
     
      x[0] = u*(2.0*u-1.0);
      x[1] = r*(2.0*r-1.0);
      x[2] = s*(2.0*s-1.0);
      x[3] = t*(2.0*t-1.0);
      x[4] = 4.0*r*u;
      x[5] = 4.0*s*r;
      x[6] = 4.0*s*u;
      x[7] = 4.0*t*u;
      x[8] = 4.0*r*t;
      x[9] = 4.0*s*t;
    }
    else {
      std::cout << "The lagrange mapping does not support cubic Tet"<<std::endl;
    }
  } 
  
  void LagrangeMapping::GeomShapeFunctionHex (double u, double v, double w, double *x) const
  {
    if( m_order == 1 ) {
      x[0] = 0.125 * (1.-u) * (1.-v) * (1.-w);
      x[1] = 0.125 * (1.+u) * (1.-v) * (1.-w);
      x[2] = 0.125 * (1.+u) * (1.+v) * (1.-w);
      x[3] = 0.125 * (1.-u) * (1.+v) * (1.-w);
      x[4] = 0.125 * (1.-u) * (1.-v) * (1.+w);
      x[5] = 0.125 * (1.+u) * (1.-v) * (1.+w);
      x[6] = 0.125 * (1.+u) * (1.+v) * (1.+w);
      x[7] = 0.125 * (1.-u) * (1.+v) * (1.+w);
    }
    else if ( m_order == 2 ) {
      x[0] = (1.0-u)*(1.0-v)*(1.0-w)/8.0-(1.0-u*u)*(1.0-v)*(1.0-w)/8.0-(1.0-u)*(1.0-v*v)*(1.0-w)/8.0-(1.0-u)*(1.0-v)*(1.0-w*w)/8.0;
      x[1] = (1.0+u)*(1.0-v)*(1.0-w)/8.0-(1.0-u*u)*(1.0-v)*(1.0-w)/8.0-(1.0+u)*(1.0-v*v)*(1.0-w)/8.0-(1.0+u)*(1.0-v)*(1.0-w*w)/8.0;
      x[2] = (1.0+u)*(1.0+v)*(1.0-w)/8.0-(1.0+u)*(1.0-v*v)*(1.0-w)/8.0-(1.0-u*u)*(1.0+v)*(1.0-w)/8.0-(1.0+u)*(1.0+v)*(1.0-w*w)/8.0;
      x[3] = (1.0-u)*(1.0+v)*(1.0-w)/8.0-(1.0-u*u)*(1.0+v)*(1.0-w)/8.0-(1.0-u)*(1.0-v*v)*(1.0-w)/8.0-(1.0-u)*(1.0+v)*(1.0-w*w)/8.0;


      x[4] = (1.0-u)*(1.0-v)*(1.0+w)/8.0-(1.0-u)*(1.0-v)*(1.0-w*w)/8.0-(1.0-u*u)*(1.0-v)*(1.0+w)/8.0-(1.0-u)*(1.0-v*v)*(1.0+w)/8.0;
      x[5] = (1.0+u)*(1.0-v)*(1.0+w)/8.0-(1.0-u*u)*(1.0-v)*(1.0+w)/8.0-(1.0+u)*(1.0-v*v)*(1.0+w)/8.0-(1.0+u)*(1.0-v)*(1.0-w*w)/8.0;
      x[6] = (1.0+u)*(1.0+v)*(1.0+w)/8.0-(1.0+u)*(1.0-v*v)*(1.0+w)/8.0-(1.0-u*u)*(1.0+v)*(1.0+w)/8.0-(1.0+u)*(1.0+v)*(1.0-w*w)/8.0;
      x[7] = (1.0-u)*(1.0+v)*(1.0+w)/8.0-(1.0-u)*(1.0-v*v)*(1.0+w)/8.0-(1.0-u*u)*(1.0+v)*(1.0+w)/8.0-(1.0-u)*(1.0+v)*(1.0-w*w)/8.0;

      x[8] = (1.0-u*u)*(1.0-v)*(1.0-w)/4.0;     
      x[9] = (1.0+u)*(1.0-v*v)*(1.0-w)/4.0;
      x[10] = (1.0-u*u)*(1.0+v)*(1.0-w)/4.0;
      x[11] = (1.0-u)*(1.0-v*v)*(1.0-w)/4.0;
      
      x[12] = (1.0-u)*(1.0-v)*(1.0-w*w)/4.0;
      x[13] = (1.0+u)*(1.0-v)*(1.0-w*w)/4.0;
      x[14] = (1.0+u)*(1.0+v)*(1.0-w*w)/4.0;
      x[15] = (1.0-u)*(1.0+v)*(1.0-w*w)/4.0;
      
      x[16] = (1.0-u*u)*(1.0-v)*(1.0+w)/4.0;
      x[17] = (1.0+u)*(1.0-v*v)*(1.0+w)/4.0;
      x[18] = (1.0-u*u)*(1.0+v)*(1.0+w)/4.0;
      x[19] = (1.0-u)*(1.0-v*v)*(1.0+w)/4.0;
      
    }
    else {
      std::cout << "The lagrange mapping does not support cubic Hex"<<std::endl;
    }
  } 

  
  void LagrangeMapping::GeomShapeFunctionPyr (double u, double v, double w, double *x) const
  {
    double ui, vi;
    
    ui = 2.*u/(1.-w); vi = 2.*v/(1.-w);

    if(m_order == 1) {
      
      x[0] = 0.125*(1.-ui)*(1.-vi)*(1.-w);
      x[1] = 0.125*(1.+ui)*(1.-vi)*(1.-w);
      x[2] = 0.125*(1.+ui)*(1.+vi)*(1.-w);
      x[3] = 0.125*(1.-ui)*(1.+vi)*(1.-w);
      
      x[4] = 0.5*(1.+w);
    
    }
    else if( m_order == 2 ) {
      x[0] = -(2.0+4.0*w*u*v+2.0*v*w+2.0*u*w+2.0*v*w*w+2.0*u*w*w-w*w*w*w+4.0*u*u*w+8.0*u*u*v+4.0*w*v*v+8.0*u*v*v-2.0*v*w*w*w-2.0*u*w*w*w+w*w*w-2.0*v-2.0*u-5.0*w-4.0*u*u-4.0*v*v+3.0*w*w-4.0*u*v*w*w)/pow(w-1.0,2.0)/8.0;
      x[1] = (-2.0+4.0*w*u*v-2.0*v*w+2.0*u*w-2.0*v*w*w+2.0*u*w*w+w*w*w*w-4.0*u*u*w-8.0*u*u*v-4.0*w*v*v+8.0*u*v*v+2.0*v*w*w*w-2.0*u*w*w*w-w*w*w+2.0*v-2.0*u+5.0*w+4.0*u*u+4.0*v*v-3.0*w*w-4.0*u*v*w*w)/pow(w-1.0,2.0)/8.0;
      x[2] = (-2.0-4.0*w*u*v+2.0*v*w+2.0*u*w+2.0*v*w*w+2.0*u*w*w+w*w*w*w-4.0*u*u*w+8.0*u*u*v-4.0*w*v*v+8.0*u*v*v-2.0*v*w*w*w-2.0*u*w*w*w-w*w*w-2.0*v-2.0*u+5.0*w+4.0*u*u+4.0*v*v-3.0*w*w+4.0*u*v*w*w)/pow(w-1.0,2.0)/8.0;
      x[3] = -(2.0-4.0*w*u*v-2.0*v*w+2.0*u*w-2.0*v*w*w+2.0*u*w*w-w*w*w*w+4.0*u*u*w-8.0*u*u*v+4.0*w*v*v+8.0*u*v*v+2.0*v*w*w*w-2.0*u*w*w*w+w*w*w+2.0*v-2.0*u-5.0*w-4.0*u*u-4.0*v*v+3.0*w*w+4.0*u*v*w*w)/pow(w-1.0,2.0)/8.0;
      x[4] = w/2.0+w*w/2.0;
      x[5] = (-w*w+2.0*w-1.0+4.0*u*u)*(w-1.0+2.0*v)/pow(w-1.0,2.0)/4.0;      
      x[6] = -(2.0*u-w+1.0)*(-w*w+2.0*w-1.0+4.0*v*v)/pow(w-1.0,2.0)/4.0;
      x[7] = -(-w*w+2.0*w-1.0+4.0*u*u)*(2.0*v-w+1.0)/pow(w-1.0,2.0)/4.0;
      x[8] = (w-1.0+2.0*u)*(-w*w+2.0*w-1.0+4.0*v*v)/pow(w-1.0,2.0)/4.0;
      x[9] = -(w-1.0+2.0*u)*(w-1.0+2.0*v)*(w+1.0)/(w-1.0)/4.0;
      x[10] = (2.0*u-w+1.0)*(w-1.0+2.0*v)*(w+1.0)/(w-1.0)/4.0;
      x[11] = -(2.0*u-w+1.0)*(2.0*v-w+1.0)*(w+1.0)/(w-1.0)/4.0;
      x[12] = (w-1.0+2.0*u)*(2.0*v-w+1.0)*(w+1.0)/(w-1.0)/4.0;
    }else{
      std::cout << "The lagrange mapping does not support high order Pyr"<<std::endl;
    }
  }
  
void LagrangeMapping::GradGeomShapeFunctionPyr (double u, double v, double w , mVector *result) const
{
    
  double ui, vi;

  ui = 2.*u/(1.-w); vi = 2.*v/(1.-w);
  
  if(m_order == 1) {
    result[0](0) = -0.25*(1.-vi); 
    result[0](1) = -0.25*(1.-ui);
    result[0](2) = 0.125*(ui*vi - 1.);
    
    result[1](0) = 0.25*(1.-vi); 
    result[1](1) = -0.25*(1.+ui);
    result[1](2) = -0.125*(ui*vi + 1.);
    
    result[2](0) = 0.25*(1.+vi); 
    result[2](1) = 0.25*(1.+ui);
    result[2](2) = 0.125*(ui*vi - 1.);
    
    result[3](0) = -0.25*(1.+vi); 
    result[3](1) = 0.25*(1.-ui);
    result[3](2) = -0.125*(ui*vi + 1.);
      
    result[4](0) = 0.0; 
    result[4](1) = 0.0;
    result[4](2) = 0.5;
  }
  else if(m_order == 2 ) {
    result[0](0) = -(4.0*v*w+2.0*w+2.0*w*w+8.0*u*w+16.0*u*v+8.0*v*v-2.0*w*w*w-2.0-8.0*u-4.0*v*w*w)/pow(w-1.0,2.0)/8.0;
    result[0](1) = -(4.0*u*w+2.0*w+2.0*w*w+8.0*u*u+8.0*v*w+16.0*u*v-2.0*w*w*w-2.0-8.0*v-4.0*u*w*w)/pow(w-1.0,2.0)/8.0;
    result[0](2) = -(4.0*u*v+2.0*v+2.0*u+4.0*v*w+4.0*u*w-4.0*w*w*w+4.0*u*u+4.0*v*v-6.0*v*w*w-6.0*u*w*w+3.0*w*w-5.0+6.0*w-8.0*w*u*v)/pow(w-1.0,2.0)/8.0+(2.0+4.0*w*u*v+2.0*v*w+2.0*u*w+2.0*v*w*w+2.0*u*w*w-w*w*w*w+4.0*u*u*w+8.0*u*u*v+4.0*w*v*v+8.0*u*v*v-2.0*v*w*w*w-2.0*u*w*w*w+w*w*w-2.0*v-2.0*u-5.0*w-4.0*u*u-4.0*v*v+3.0*w*w-4.0*u*v*w*w)/pow(w-1.0,3.0)/4.0;
    
    result[1](0) = (4.0*v*w+2.0*w+2.0*w*w-8.0*u*w-16.0*u*v+8.0*v*v-2.0*w*w*w-2.0+8.0*u-4.0*v*w*w)/pow(w-1.0,2.0)/8.0;
      result[1](1) = (4.0*u*w-2.0*w-2.0*w*w-8.0*u*u-8.0*v*w+16.0*u*v+2.0*w*w*w+2.0+8.0*v
		      -4.0*u*w*w)/pow(w-1.0,2.0)/8.0;
      result[1](2) = (4.0*u*v-2.0*v+2.0*u-4.0*v*w+4.0*u*w+4.0*w*w*w-4.0*u*u-4.0*v*v+6.0*v*w*w-6.0*u*w*w-3.0*w*w+5.0-6.0*w-8.0*w*u*v)/pow(w-1.0,2.0)/8.0-(-2.0+4.0*w*u*v-2.0*v*w+2.0*u*w-2.0*v*w*w+2.0*u*w*w+w*w*w*w-4.0*u*u*w-8.0*u*u*v-4.0*w*v*v+8.0*u*v*v+2.0*v*w*w*w-2.0*u*w*w*w-w*w*w+2.0*v-2.0*u+5.0*w+4.0*u*u+4.0*v*v-3.0*w*w-4.0*u*v*w*w)/pow(w-1.0,3.0)/4.0;
      
      result[2](0) = (-4.0*v*w+2.0*w+2.0*w*w-8.0*u*w+16.0*u*v+8.0*v*v-2.0*w*w*w-2.0+8.0*u+4.0*v*w*w)/pow(w-1.0,2.0)/8.0;
      result[2](1) = (-4.0*u*w+2.0*w+2.0*w*w+8.0*u*u-8.0*v*w+16.0*u*v-2.0*w*w*w-2.0+8.0*v+4.0*u*w*w)/pow(w-1.0,2.0)/8.0;
      result[2](2) = (-4.0*u*v+2.0*v+2.0*u+4.0*v*w+4.0*u*w+4.0*w*w*w-4.0*u*u-4.0*v*v-6.0*v*w*w-6.0*u*w*w-3.0*w*w+5.0-6.0*w+8.0*w*u*v)/pow(w-1.0,2.0)/8.0-(-2.0-4.0*w*u*v+2.0*v*w+2.0*u*w+2.0*v*w*w+2.0*u*w*w+w*w*w*w-4.0*u*u*w+8.0*u*u*v-4.0*w*v*v+8.0*u*v*v-2.0*v*w*w*w-2.0*u*w*w*w-w*w*w-2.0*v-2.0*u+5.0*w+4.0*u*u+4.0*v*v-3.0*w*w+4.0*u*v*w*w)/pow(w-1.0,3.0)/4.0;
      
      result[3](0) = -(-4.0*v*w+2.0*w+2.0*w*w+8.0*u*w-16.0*u*v+8.0*v*v-2.0*w*w*w-2.0-8.0*u+4.0*v*w*w)/pow(w-1.0,2.0)/8.0;
      result[3](1) = -(-4.0*u*w-2.0*w-2.0*w*w-8.0*u*u+8.0*v*w+16.0*u*v+2.0*w*w*w+2.0-8.0*v+4.0*u*w*w)/pow(w-1.0,2.0)/8.0;
      result[3](2) = -(-4.0*u*v-2.0*v+2.0*u-4.0*v*w+4.0*u*w-4.0*w*w*w+4.0*u*u+4.0*v*v+6.0*v*w*w-6.0*u*w*w+3.0*w*w-5.0+6.0*w+8.0*w*u*v)/pow(w-1.0,2.0)/8.0+(2.0-4.0*w*u*v-2.0*v*w+2.0*u*w-2.0*v*w*w+2.0*u*w*w-w*w*w*w+4.0*u*u*w-8.0*u*u*v+4.0*w*v*v+8.0*u*v*v+2.0*v*w*w*w-2.0*u*w*w*w+w*w*w+2.0*v-2.0*u-5.0*w-4.0*u*u-4.0*v*v+3.0*w*w+4.0*u*v*w*w)/pow(w-1.0,3.0)/4.0;
      
      result[4](0) = 0.0;      
      result[4](1) = 0.0;
      result[4](2) = 1.0/2.0+w;

      result[5](0) = 2.0*u*(w-1.0+2.0*v)/pow(w-1.0,2.0);      
      result[5](1) = (-w*w+2.0*w-1.0+4.0*u*u)/pow(w-1.0,2.0)/2.0;
      result[5](2) = (-2.0*w+2.0)*(w-1.0+2.0*v)/pow(w-1.0,2.0)/4.0+(-w*w+2.0*w-1.0+4.0*u*u)/pow(w-1.0,2.0)/4.0-(-w*w+2.0*w-1.0+4.0*u*u)*(w-1.0+2.0*v)/pow(w-1.0,3.0)/2.0;
      
      result[6](0) = -(-w*w+2.0*w-1.0+4.0*v*v)/pow(w-1.0,2.0)/2.0;      
      result[6](1) = -2.0*(2.0*u-w+1.0)*v/pow(w-1.0,2.0);
      result[6](2) = (-w*w+2.0*w-1.0+4.0*v*v)/pow(w-1.0,2.0)/4.0-(2.0*u-w+1.0)*(-2.0*w+2.0)/pow(w-1.0,2.0)/4.0+(2.0*u-w+1.0)*(-w*w+2.0*w-1.0+4.0*v*v)/pow(w-1.0,3.0)/2.0;
      
      result[7](0) = -2.0*u*(2.0*v-w+1.0)/pow(w-1.0,2.0);      
      result[7](1) = -(-w*w+2.0*w-1.0+4.0*u*u)/pow(w-1.0,2.0)/2.0;
      result[7](2) = -(-2.0*w+2.0)*(2.0*v-w+1.0)/pow(w-1.0,2.0)/4.0+(-w*w+2.0*w-1.0+4.0*u*u)/pow(w-1.0,2.0)/4.0+(-w*w+2.0*w-1.0+4.0*u*u)*(2.0*v-w+1.0)/pow(w-1.0,3.0)/2.0;

      result[8](0) = (-w*w+2.0*w-1.0+4.0*v*v)/pow(w-1.0,2.0)/2.0;      
      result[8](1) = 2.0*(w-1.0+2.0*u)*v/pow(w-1.0,2.0);
      result[8](2) = (-w*w+2.0*w-1.0+4.0*v*v)/pow(w-1.0,2.0)/4.0+(w-1.0+2.0*u)*(-2.0*w+2.0)/pow(w-1.0,2.0)/4.0-(w-1.0+2.0*u)*(-w*w+2.0*w-1.0+4.0*v*v)/pow(w-1.0,3.0)/2.0;

      result[9](0) = -(w-1.0+2.0*v)*(w+1.0)/(w-1.0)/2.0;      
      result[9](1) = -(w-1.0+2.0*u)*(w+1.0)/(w-1.0)/2.0;
      result[9](2) = -(w-1.0+2.0*v)*(w+1.0)/(w-1.0)/4.0-(w-1.0+2.0*u)*(w+1.0)/(w-1.0)/4.0-(w-1.0+2.0*u)*(w-1.0+2.0*v)/(w-1.0)/4.0+(w-1.0+2.0*u)*(w-1.0+2.0*v)*(w+1.0)/pow(w-1.0,2.0)/4.0;
      
      result[10](0) = (w-1.0+2.0*v)*(w+1.0)/(w-1.0)/2.0;
      result[10](1) = (2.0*u-w+1.0)*(w+1.0)/(w-1.0)/2.0;      
      result[10](2) = -(w-1.0+2.0*v)*(w+1.0)/(w-1.0)/4.0+(2.0*u-w+1.0)*(w+1.0)/(w-1.0)/4.0+(2.0*u-w+1.0)*(w-1.0+2.0*v)/(w-1.0)/4.0-(2.0*u-w+1.0)*(w-1.0+2.0*v)*(w+1.0)/pow(w-1.0,2.0)/4.0;

      result[11](0) = -(2.0*v-w+1.0)*(w+1.0)/(w-1.0)/2.0;      
      result[11](1) = -(2.0*u-w+1.0)*(w+1.0)/(w-1.0)/2.0;
      result[11](2) = (2.0*v-w+1.0)*(w+1.0)/(w-1.0)/4.0+(2.0*u-w+1.0)*(w+1.0)/(w-1.0)/4.0-(2.0*u-w+1.0)*(2.0*v-w+1.0)/(w-1.0)/4.0+(2.0*u-w+1.0)*(2.0*v-w+1.0)*(w+1.0)/pow(w-1.0,2.0)/4.0;
      
      result[12](0) = (2.0*v-w+1.0)*(w+1.0)/(w-1.0)/2.0;      
      result[12](1) = (w-1.0+2.0*u)*(w+1.0)/(w-1.0)/2.0;
      result[12](2) = (2.0*v-w+1.0)*(w+1.0)/(w-1.0)/4.0-(w-1.0+2.0*u)*(w+1.0)/(w-1.0)/4.0+(w-1.0+2.0*u)*(2.0*v-w+1.0)/(w-1.0)/4.0-(w-1.0+2.0*u)*(2.0*v-w+1.0)*(w+1.0)/pow(w-1.0,2.0)/4.0;

    }
    else {
      std::cout << "The lagrange mapping does not support high order Pyr"<<std::endl;
    }
  }


  void LagrangeMapping::GeomShapeFunctionTri ( double r, double s, double dum, double *x) const
  {
    double t = 1.-r-s;
    if(m_order == 1) {
      x[0] = t;
      x[1] = r;
      x[2] = s;
    }
    else if(m_order == 2) {
      x[0] = t * (2.*t - 1.);
      x[1] = r * (2.*r - 1.);
      x[2] = s * (2.*s - 1.);

      x[3] = 4. * t * r;
      x[4] = 4. * r * s;
      x[5] = 4. * t * s;
    }
    else {
      std::cout << "The lagrange mapping does not support cubic Tri"<<std::endl;
    }
  }

  void LagrangeMapping::GradGeomShapeFunctionTri (double u, double v, double dum , mVector *result) const
  {
   if (m_order == 1) {
      result[0] (0)  = -1.0;
      result[1] (0)  = 1.0;
      result[2] (0)  = 0.0;
      result[0] (1)  = -1.0;
      result[1] (1)  = 0.0;
      result[2] (1)  = 1.0;    
      result[0] (2) = 0.0; 
      result[1] (2) = 0.0; 
      result[2] (2) = 0.0; 
    } 
    else if (m_order == 2) {
    
      double w = 1.-u-v;
      double t2 = 1.0-4.0*w;
      double t3 = 4.0*u;
      double t6 = 4.0*v;
      
      result[0] (0)  = t2;
      result[1] (0)  = t3-1.0;
      result[2] (0)  = 0.0;
      result[3] (0)  = 4.0*w-4.0*u;
      result[4] (0)  = t6;
      result[5] (0)  = -t6;
      result[0] (1)  = t2;
      result[1] (1)  = 0.0;
      result[2] (1)  = t6-1.0;
      result[3] (1)  = -t3;
      result[4] (1)  = t3;
      result[5] (1)  = 4.0*w-4.0*v;
      
      result[0] (2) = 0.0; 
      result[1] (2) = 0.0; 
      result[2] (2) = 0.0; 
      result[3] (2) = 0.0; 
      result[4] (2) = 0.0; 
      result[5] (2) = 0.0; 
    }
    else {
      std::cout << "The lagrange mapping does not support cubic Tri"<<std::endl;
    }
   
  }

  void LagrangeMapping::GradGeomShapeFunctionQuad (double u, double v, double dum , mVector *result) const
  {
    if(m_order == 1) {
      double t1 = v-1.0;
      double t2 = v+1.0;
      double t3 = u-1.0;
      double t4 = -u-1.0;
      result[0] (0)  = t1/4.0;
      result[1] (0)  = -t1/4.0;
      result[2] (0)  = t2/4.0;
      result[3] (0)  = -t2/4.0;
      result[0] (1)  = t3/4.0;
      result[1] (1)  = t4/4.0;
      result[2] (1)  = -t4/4.0;
      result[3] (1)  = -t3/4.0;    
      result[0] (2) = 0.0; 
      result[1] (2) = 0.0; 
      result[2] (2) = 0.0; 
      result[3] (2) = 0.0; 
    }
    else if (m_order == 2) {
      
      result[0](0) = v/4.0+u*(1.0-v)/2.0-v*v/4.0;
      result[0](1) = u/4.0-u*u/4.0+(1.0-u)*v/2.0;
      
      result[1](0) = -v/4.0+u*(1.0-v)/2.0+v*v/4.0;
      result[1](1) = -u/4.0-u*u/4.0+(1.0+u)*v/2.0;

      result[2](0) = v/4.0+v*v/4.0+u*(1.0+v)/2.0;
      result[2](1) = u/4.0+(1.0+u)*v/2.0+u*u/4.0;

      result[3](0) = -v/4.0+u*(1.0+v)/2.0-v*v/4.0;
      result[3](1) = -u/4.0+u*u/4.0+(1.0-u)*v/2.0;

      result[4](0) = -u*(1.0-v);
      result[4](1) = -1.0/2.0+u*u/2.0;

      result[5](0) = 1.0/2.0-v*v/2.0;
      result[5](1) = -(1.0+u)*v;

      result[6](0) = -u*(1.0+v);
      result[6](1) = 1.0/2.0-u*u/2.0;

      result[7](0) = -1.0/2.0+v*v/2.0;
      result[7](1) = -(1.0-u)*v;

      
      result[0](2) = result[1](2) = result[2](2) = result[3](2) = result[4](2) = result[5](2) = result[6](2) = result[7](2) = 0.0;

    }
    else {
      std::cout << "The lagrange mapping does not support cubic Quad"<<std::endl; 
    }
  } 
  
  void LagrangeMapping::GradGeomShapeFunctionTet (double r, double s, double t, mVector *result) const
  {
   if(m_order == 1) {
      result[0] (0)  = -1.0;
      result[1] (0)  = 1.0;
      result[2] (0)  = 0.0;
      result[3] (0)  = 0.0;
      result[0] (1)  = -1.0;
      result[1] (1)  = 0.0;
      result[2] (1)  = 1.0;
      result[3] (1)  = 0.0;
      result[0] (2)  = -1.0;
      result[1] (2)  = 0.0;
      result[2] (2)  = 0.0;
      result[3] (2)  = 1.0;
    }
    else if(m_order == 2) {
      double u=1.-r-s-t;
      result[0] (0)  = 1.-4.*u;
      result[0] (1)  = 1.-4.*u;
      result[0] (2)  = 1.-4.*u;

      result[1] (0)  = 4.*r-1.;
      result[1] (1)  = 0.;
      result[1] (2)  = 0.;
      
      result[2] (0)  = 0.0;
      result[2] (1)  = 4.*s-1.;
      result[2] (2)  = 0.0;

      result[3] (0)  = 0.0;
      result[3] (1)  = 0.0;
      result[3] (2)  = 4.*t-1.;

      result[4] (0)  = 4.0*u-4.0*r;
      result[4] (1)  = -4.*r;
      result[4] (2)  = -4.0*r;

      result[5] (0)  = 4.0*s;
      result[5] (1)  = 4.*r;
      result[5] (2)  = 0.;

      result[6] (0)  = -4.*s;
      result[6] (1)  = 4.*u-4.*s;
      result[6] (2)  = -4.*s;

      result[7] (0)  = -4.0*t;
      result[7] (1)  = -4.*t;
      result[7] (2)  = 4.0*u-4.*t;

      result[8] (0)  = 4.0*t;
      result[8] (1)  = 0.;
      result[8] (2)  = 4.0*r;

      result[9] (0)  = 0.;
      result[9] (1)  = 4.*t;
      result[9] (2)  = 4.0*s;

    }
    else {
      std::cout << "The lagrange mapping does not support cubic Quad"<<std::endl; 
    }
  } 
  
  void LagrangeMapping::GradGeomShapeFunctionHex (double u, double v, double w , mVector *grad) const
  {
    if(m_order == 1) {
      grad[0][0] = -0.125*(1.-v)*(1.-w); 
      grad[0][1]= -0.125*(1.-u)*(1.-w); 
      grad[0][2] = -0.125*(1.-u)*(1.-v);
      
      grad[1][0] =  0.125*(1.-v)*(1.-w); 
      grad[1][1]= -0.125*(1.+u)*(1.-w); 
      grad[1][2] = -0.125*(1.+u)*(1.-v);
      
      grad[2][0] =  0.125*(1.+v)*(1.-w); 
      grad[2][1]=  0.125*(1.+u)*(1.-w); 
      grad[2][2] = -0.125*(1.+u)*(1.+v);
      
      grad[3][0] = -0.125*(1.+v)*(1.-w); 
      grad[3][1]=  0.125*(1.-u)*(1.-w); 
      grad[3][2] = -0.125*(1.-u)*(1.+v);
      
      grad[4][0] = -0.125*(1.-v)*(1.+w); 
      grad[4][1]= -0.125*(1.-u)*(1.+w); 
      grad[4][2] =  0.125*(1.-u)*(1.-v);
      
      grad[5][0] =  0.125*(1.-v)*(1.+w); 
      grad[5][1]= -0.125*(1.+u)*(1.+w); 
      grad[5][2] =  0.125*(1.+u)*(1.-v);
      
      grad[6][0] =  0.125*(1.+v)*(1.+w); 
      grad[6][1]=  0.125*(1.+u)*(1.+w); 
      grad[6][2] =  0.125*(1.+u)*(1.+v);
    
      grad[7][0] = -0.125*(1.+v)*(1.+w); 
      grad[7][1]=  0.125*(1.-u)*(1.+w); 
      grad[7][2] =  0.125*(1.-u)*(1.+v);
    }
    else if(m_order == 2) {
      //*********0
      grad[0][0] = -(1.0-v)*(1.0-w)/8.0+u*(1.0-v)*(1.0-w)/4.0+(1.0-v*v)*(1.0-w)/8.0+(1.0-v)*(1.0-w*w)/8.0;
      grad[0][1] = -(1.0-u)*(1.0-w)/8.0+(1.0-u*u)*(1.0-w)/8.0+(1.0-u)*v*(1.0-w)/4.0+(1.0-u)*(1.0-w*w)/8.0;
      grad[0][2] = -(1.0-u)*(1.0-v)/8.0+(1.0-u*u)*(1.0-v)/8.0+(1.0-u)*(1.0-v*v)/8.0+(1.0-u)*(1.0-v)*w/4.0;

      //*********1
      grad[1][0] = (1.0-v)*(1.0-w)/8.0+u*(1.0-v)*(1.0-w)/4.0-(1.0-v*v)*(1.0-w)/8.0-(1.0-v)*(1.0-w*w)/8.0;
      grad[1][1] = -(1.0+u)*(1.0-w)/8.0+(1.0-u*u)*(1.0-w)/8.0+(1.0+u)*v*(1.0-w)/4.0+(1.0+u)*(1.0-w*w)/8.0;
      grad[1][2] = -(1.0+u)*(1.0-v)/8.0+(1.0-u*u)*(1.0-v)/8.0+(1.0+u)*(1.0-v*v)/8.0+(1.0+u)*(1.0-v)*w/4.0;


      //*********2
      grad[2][0] = (1.0+v)*(1.0-w)/8.0-(1.0-v*v)*(1.0-w)/8.0+u*(1.0+v)*(1.0-w)/4.0-(1.0+v)*(1.0-w*w)/8.0;
      grad[2][1] = (1.0+u)*(1.0-w)/8.0+(1.0+u)*v*(1.0-w)/4.0-(1.0-u*u)*(1.0-w)/8.0-(1.0+u)*(1.0-w*w)/8.0;
      grad[2][2] = -(1.0+u)*(1.0+v)/8.0+(1.0+u)*(1.0-v*v)/8.0+(1.0-u*u)*(1.0+v)/8.0+(1.0+u)*(1.0+v)*w/4.0;

      //*********3
      grad[3][0] = -(1.0+v)*(1.0-w)/8.0+u*(1.0+v)*(1.0-w)/4.0+(1.0-v*v)*(1.0-w)/8.0+(1.0+v)*(1.0-w*w)/8.0;
      grad[3][1] = (1.0-u)*(1.0-w)/8.0-(1.0-u*u)*(1.0-w)/8.0+(1.0-u)*v*(1.0-w)/4.0-(1.0-u)*(1.0-w*w)/8.0;
      grad[3][2] = -(1.0-u)*(1.0+v)/8.0+(1.0-u*u)*(1.0+v)/8.0+(1.0-u)*(1.0-v*v)/8.0+(1.0-u)*(1.0+v)*w/4.0;

      //*********4
      grad[4][0] = -(1.0-v)*(1.0+w)/8.0+(1.0-v)*(1.0-w*w)/8.0+u*(1.0-v)*(1.0+w)/4.0+(1.0-v*v)*(1.0+w)/8.0;
      grad[4][1] = -(1.0-u)*(1.0+w)/8.0+(1.0-u)*(1.0-w*w)/8.0+(1.0-u*u)*(1.0+w)/8.0+(1.0-u)*v*(1.0+w)/4.0;
      grad[4][2] = (1.0-u)*(1.0-v)/8.0+(1.0-u)*(1.0-v)*w/4.0-(1.0-u*u)*(1.0-v)/8.0-(1.0-u)*(1.0-v*v)/8.0;

      //*********5
      grad[5][0] = (1.0-v)*(1.0+w)/8.0+u*(1.0-v)*(1.0+w)/4.0-(1.0-v*v)*(1.0+w)/8.0-(1.0-v)*(1.0-w*w)/8.0;
      grad[5][1] = -(1.0+u)*(1.0+w)/8.0+(1.0-u*u)*(1.0+w)/8.0+(1.0+u)*v*(1.0+w)/4.0+(1.0+u)*(1.0-w*w)/8.0;
      grad[5][2] = (1.0+u)*(1.0-v)/8.0-(1.0-u*u)*(1.0-v)/8.0-(1.0+u)*(1.0-v*v)/8.0+(1.0+u)*(1.0-v)*w/4.0;


      //*********6
      grad[6][0] = (1.0+v)*(1.0+w)/8.0-(1.0-v*v)*(1.0+w)/8.0+u*(1.0+v)*(1.0+w)/4.0-(1.0+v)*(1.0-w*w)/8.0;
      grad[6][1] = (1.0+u)*(1.0+w)/8.0+(1.0+u)*v*(1.0+w)/4.0-(1.0-u*u)*(1.0+w)/8.0-(1.0+u)*(1.0-w*w)/8.0;
      grad[6][2] = (1.0+u)*(1.0+v)/8.0-(1.0+u)*(1.0-v*v)/8.0-(1.0-u*u)*(1.0+v)/8.0+(1.0+u)*(1.0+v)*w/4.0;

      //*********7
      grad[7][0] = -(1.0+v)*(1.0+w)/8.0+(1.0-v*v)*(1.0+w)/8.0+u*(1.0+v)*(1.0+w)/4.0+(1.0+v)*(1.0-w*w)/8.0;
      grad[7][1] = (1.0-u)*(1.0+w)/8.0+(1.0-u)*v*(1.0+w)/4.0-(1.0-u*u)*(1.0+w)/8.0-(1.0-u)*(1.0-w*w)/8.0;
      grad[7][2] = (1.0-u)*(1.0+v)/8.0-(1.0-u)*(1.0-v*v)/8.0-(1.0-u*u)*(1.0+v)/8.0+(1.0-u)*(1.0+v)*w/4.0;

      //*********8
      grad[8][0] = -u*(1.0-v)*(1.0-w)/2.;
      grad[8][1] = -(1.0-u*u)*(1.0-w)/4.0;
      grad[8][2] = -(1.0-u*u)*(1.0-v)/4.0;

      //*********9
      grad[9][0] = (1.0-v*v)*(1.0-w)/4.0;
      grad[9][1] = -(1.0+u)*v*(1.0-w)/2.;
      grad[9][2] = -(1.0+u)*(1.0-v*v)/4.0;

      //*********10
      grad[10][0] = -u*(1.0+v)*(1.0-w)/2.;
      grad[10][1] = (1.0-u*u)*(1.0-w)/4.0;
      grad[10][2] = -(1.0-u*u)*(1.0+v)/4.0;

      //*********11
      grad[11][0] = -(1.0-v*v)*(1.0-w)/4.0;
      grad[11][1] = -(1.0-u)*v*(1.0-w)/2.;
      grad[11][2] = -(1.0-u)*(1.0-v*v)/4.0;

      //*********12
      grad[12][0] = -(1.0-v)*(1.0-w*w)/4.0;
      grad[12][1] = -(1.0-u)*(1.0-w*w)/4.0;
      grad[12][2] = -(1.0-u)*(1.0-v)*w/2.;

      //*********13
      grad[13][0] = (1.0-v)*(1.0-w*w)/4.0;
      grad[13][1] = -(1.0+u)*(1.0-w*w)/4.0;
      grad[13][2] = -(1.0+u)*(1.0-v)*w/2.;

      //*********14
      grad[14][0] = (1.0+v)*(1.0-w*w)/4.0;
      grad[14][1] = (1.0+u)*(1.0-w*w)/4.0;
      grad[14][2] = -(1.0+u)*(1.0+v)*w/2.;

      //*********15
      grad[15][0] = -(1.0+v)*(1.0-w*w)/4.0;
      grad[15][1] = (1.0-u)*(1.0-w*w)/4.0;
      grad[15][2] = -(1.0-u)*(1.0+v)*w/2.;

      //*********16
      grad[16][0] = -u*(1.0-v)*(1.0+w)/2.;
      grad[16][1] = -(1.0-u*u)*(1.0+w)/4.0;
      grad[16][2] = (1.0-u*u)*(1.0-v)/4.0;

      //*********17
      grad[17][0] = (1.0-v*v)*(1.0+w)/4.0;
      grad[17][1] = -(1.0+u)*v*(1.0+w)/2.;
      grad[17][2] = (1.0+u)*(1.0-v*v)/4.0;

      //*********18
      grad[18][0] = -u*(1.0+v)*(1.0+w)/2.;
      grad[18][1] = (1.0-u*u)*(1.0+w)/4.0;
      grad[18][2] = (1.0-u*u)*(1.0+v)/4.0;

      //*********19
      grad[19][0] = -(1.0-v*v)*(1.0+w)/4.0;
      grad[19][1] = -(1.0-u)*v*(1.0+w)/2.;
      grad[19][2] = (1.0-u)*(1.0-v*v)/4.0;
    }

    else {
      std::cout << "The lagrange mapping does not support cubic Hex"<<std::endl; 
    }

  } 

  void LagrangeMapping::GeomShapeFunctionLine (double u, double v, double w, double *result) const
  {
    if(m_order == 1) {
      result[0] = -u/2.0+1.0/2.0;
      result[1] = u/2.0+1.0/2.0;
    }
    else if(m_order == 2) {
      result[0] = u*(u-1.)/2.;
      result[1] = u*(u+1.)/2.;
      result[2] = 1.-u*u;
    }
    else {
      std::cout << "The lagrange mapping does not support cubic Line"<<std::endl;
    }
  }

  void LagrangeMapping::GradGeomShapeFunctionLine (double u, double v, double w, mVector *grad) const
  {
     if(m_order == 1) {
      grad[0](1) = grad[0](2) = 0.0;
      grad[1](1) = grad[1](2) = 0.0;
      grad[0](0) = (-0.5);
      grad[1](0) = (+0.5);
    }
    else if(m_order == 2) {
      grad[0](1) = grad[0](2) = 0.0;
      grad[1](1) = grad[1](2) = 0.0;
      grad[2](1) = grad[2](2) = 0.0;

      grad[0][0] = u-0.5;
      grad[1][0] = u+0.5;
      grad[2][0]= -2.*u;
    }
    else {
      std::cout << "The lagrange mapping does not support cubic Line"<<std::endl;
    }
  }
  
  void LagrangeMapping::GeomShapeFunctionPrism (double u, double v, double w, double *x) const
  {
    double k = 1.-u-v;
    
    if(m_order == 1) {
      x[0] =  (k) * (0.5*(1.-w));
      x[1] =  (u) * (0.5*(1.-w));
      x[2] =  (v) * (0.5*(1.-w));
      x[3] =  (k) * (0.5*(1.+w));
      x[4] =  (u) * (0.5*(1.+w));
      x[5] =  (v) * (0.5*(1.+w));
    }
    else if(m_order == 2) {
      
      x[0] = 0.5e0 * (0.1e1 - w) * k - 0.10e1 * (0.1e1 - w) * k * u - 0.10e1 * (0.1e1 - w) * k * v - 0.5e0 * k * (0.1e1 - w * w);
      x[1] = 0.5e0 * (0.1e1 - w) * u - 0.10e1 * (0.1e1 - w) * k * u - 0.10e1 * (0.1e1 - w) * u * v - 0.5e0 * u * (0.1e1 - w * w);
      x[2] = 0.5e0 * (0.1e1 - w) * v - 0.10e1 * (0.1e1 - w) * u * v - 0.10e1 * (0.1e1 - w) * k * v - 0.5e0 * v * (0.1e1 - w * w);
      x[3] = 0.5e0 * (0.1e1 + w) * k - 0.5e0 * k * (0.1e1 - w * w) - 0.10e1 * (0.1e1 + w) * k * u - 0.10e1 * (0.1e1 + w) * k * v;
      x[4] = 0.5e0 * (0.1e1 + w) * u - 0.5e0 * u * (0.1e1 - w * w) - 0.10e1 * (0.1e1 + w) * k * u - 0.10e1 * (0.1e1 + w) * u * v;
      x[5] = 0.5e0 * (0.1e1 + w) * v - 0.5e0 * v * (0.1e1 - w * w) - 0.10e1 * (0.1e1 + w) * u * v - 0.10e1 * (0.1e1 + w) * k * v;
      x[6] = 0.2e1 * (0.1e1 - w) * k * u;
      x[7] = 0.2e1 * (0.1e1 - w) * u * v;
      x[8] = 0.2e1 * (0.1e1 - w) * k * v;
      x[9] = k * (0.1e1 - w * w);
      x[10] = u * (0.1e1 - w * w);
      x[11] = v * (0.1e1 - w * w);
      x[12] = 0.2e1 * (0.1e1 + w) * k * u;
      x[13] = 0.2e1 * (0.1e1 + w) * u * v;
      x[14] = 0.2e1 * (0.1e1 + w) * k * v;
	}
    else {
      std::cout << "The lagrange mapping does not support cubic Prism"<<std::endl;
    }
  } 

  void LagrangeMapping::GradGeomShapeFunctionPrism (double u, double v, double w , mVector *grad) const
  {
    double k = 1.-u-v;

    if(m_order == 1) {
      grad[0][0] = -0.25*(1.-w); 
      grad[0][1]= -0.25*(1.-w); 
      grad[0][2] = -0.5*(1.-u-v);
      grad[1][0] =  0.25*(1.-w); 
      grad[1][1]= 0.0; 
      grad[1][2] = -0.5*(u);
      grad[2][0] =  0.0; 
      grad[2][1]= 0.25*(1.-w); 
      grad[2][2] = -0.5*(v);
      grad[3][0] = -0.25*(1.+w); 
      grad[3][1]= -0.25*(1.+w); 
      grad[3][2] = 0.5*(1.-u-v);
      grad[4][0] =  0.25*(1.+w); 
      grad[4][1]= 0.0; 
      grad[4][2] = 0.5*(u);
      grad[5][0] =  0.0; 
      grad[5][1]= 0.25*(1.+w); 
      grad[5][2] = 0.5*(v);
    }
    else if (m_order == 2) {
      grad[0][0] = 0.25e0 * w + 0.50e0 * (0.1e1 - w) * u + 0.50e0 * (0.1e1 - w) * v - 0.25e0 * w * w - 0.50e0 * (0.1e1 - w) * k;
      grad[0][1] = 0.25e0 * w + 0.50e0 * (0.1e1 - w) * u + 0.50e0 * (0.1e1 - w) * v - 0.25e0 * w * w - 0.50e0 * (0.1e1 - w) * k;
      grad[0][2] = -0.5e0 * k + 0.10e1 * k * u + 0.10e1 * k * v + 0.10e1 * k * w;
      grad[1][0] = 0.50e0 * (0.1e1 - w) * u - 0.25e0 * w - 0.50e0 * (0.1e1 - w) * k - 0.50e0 * (0.1e1 - w) * v + 0.25e0 * w * w;
      grad[1][1] = 0.0e0;
      grad[1][2] = -0.5e0 * u + 0.10e1 * k * u + 0.10e1 * u * v + 0.10e1 * u * w;
      grad[2][0] = 0.0e0;
      grad[2][1] = 0.50e0 * (0.1e1 - w) * v - 0.25e0 * w - 0.50e0 * (0.1e1 - w) * u - 0.50e0 * (0.1e1 - w) * k + 0.25e0 * w * w;
      grad[2][2] = -0.5e0 * v + 0.10e1 * u * v + 0.10e1 * k * v + 0.10e1 * v * w;
      grad[3][0] = -0.25e0 * w - 0.25e0 * w * w + 0.50e0 * (0.1e1 + w) * u + 0.50e0 * (0.1e1 + w) * v - 0.50e0 * (0.1e1 + w) * k;
      grad[3][1] = -0.25e0 * w - 0.25e0 * w * w + 0.50e0 * (0.1e1 + w) * u + 0.50e0 * (0.1e1 + w) * v - 0.50e0 * (0.1e1 + w) * k;

      grad[3][2] = 0.5e0 * k + 0.10e1 * k * w - 0.10e1 * k * u - 0.10e1 * k * v;
      grad[4][0] = 0.50e0 * (0.1e1 + w) * u + 0.25e0 * w + 0.25e0 * w * w - 0.50e0 * (0.1e1 + w) * k - 0.50e0 * (0.1e1 + w) * v;
      grad[4][1] = 0.0e0;
      grad[4][2] = 0.5e0 * u + 0.10e1 * u * w - 0.10e1 * k * u - 0.10e1 * u * v;
      grad[5][0] = 0.0e0;
      grad[5][1] = 0.50e0 * (0.1e1 + w) * v + 0.25e0 * w + 0.25e0 * w * w - 0.50e0 * (0.1e1 + w) * u - 0.50e0 * (0.1e1 + w) * k;
      grad[5][2] = 0.5e0 * v + 0.10e1 * v * w - 0.10e1 * u * v - 0.10e1 * k * v;
      grad[6][0] = -0.10e1 * (0.1e1 - w) * u + 0.10e1 * (0.1e1 - w) * k;
      grad[6][1] = -0.10e1 * (0.1e1 - w) * u;
      grad[6][2] = -0.2e1 * k * u;
      grad[7][0] = 0.10e1 * (0.1e1 - w) * v;
      grad[7][1] = 0.10e1 * (0.1e1 - w) * u;
      grad[7][2] = -0.2e1 * u * v;

      grad[8][0] = -0.10e1 * (0.1e1 - w) * v;
      grad[8][1] = -0.10e1 * (0.1e1 - w) * v + 0.10e1 * (0.1e1 - w) * k;
      grad[8][2] = -0.2e1 * k * v;
      grad[9][0] = -0.5e0 + 0.5e0 * w * w;
      grad[9][1] = -0.5e0 + 0.5e0 * w * w;
      grad[9][2] = -0.2e1 * k * w;
      grad[10][0] = 0.5e0 - 0.5e0 * w * w;
      grad[10][1] = 0.0e0;
      grad[10][2] = -0.2e1 * u * w;
      grad[11][0] = 0.0e0;

      grad[11][1] = 0.5e0 - 0.5e0 * w * w;
      grad[11][2] = -0.2e1 * v * w;
      grad[12][0] = -0.10e1 * (0.1e1 + w) * u + 0.10e1 * (0.1e1 + w) * k;
      grad[12][1] = -0.10e1 * (0.1e1 + w) * u;
      grad[12][2] = 0.2e1 * k * u;
      grad[13][0] = 0.10e1 * (0.1e1 + w) * v;
      grad[13][1] = 0.10e1 * (0.1e1 + w) * u;
      grad[13][2] = 0.2e1 * u * v;
      grad[14][0] = -0.10e1 * (0.1e1 + w) * v;
      grad[14][1] = -0.10e1 * (0.1e1 + w) * v + 0.10e1 * (0.1e1 + w) * k;
      grad[14][2] = 0.2e1 * k * v;

    }
    else {
      std::cout << "The lagrange mapping does not support cubic Prism"<<std::endl;
    }
  } 


  void LagrangeMapping::eval(double u, double v, double w, 
			     double &x, double &y, double &z) const
  {

    x = 0.0, y = 0.0, z = 0.0;
    double f[256];
    GeomShapeFunction (u,v,w,f);
    for(size_t i=0;i<knots.size();i++)
      {
	mPoint p =  knots[i];
	double fct = f[i];
	x += p(0) * fct;
	y += p(1) * fct;
	z += p(2) * fct;
      }
  }

  void LagrangeMapping::deval(double u, double v, double w, 
			      double &dxdu, double &dydu, double &dzdu,
			      double &dxdv, double &dydv, double &dzdv,
			      double &dxdw, double &dydw, double &dzdw) const
  {
    dxdu = 0.0, dydu = 0.0, dzdu = 0.0;
    dxdv = 0.0, dydv = 0.0, dzdv = 0.0;
    dxdw = 0.0, dydw = 0.0, dzdw = 0.0;
    mVector dus[256];
    GradGeomShapeFunction  (u ,  v ,  w, dus);
    for (size_t i = 0; i != knots.size(); ++i) {
	mVector du(dus[i]);
	mPoint p =  knots[i];

	const double xx = p(0);
	const double yy = p(1);
	const double zz = p(2);

	dxdu += xx * du(0);
	dydu += yy * du(0);
	dzdu += zz * du(0);

	dxdv += xx * du(1);
	dydv += yy * du(1);
	dzdv += zz * du(1);

	dxdw += xx * du(2);
	dydw += yy * du(2);
	dzdw += zz * du(2);
      }
  }

  void LagrangeMapping::normalVector(pEntity border, double u, double v, double w, mVector& n) const
  {

    // In case of a mesh mapping, the normal to an entity which is on
    // a border = sum of grad geom shape functions which node are not
    // on the border entity

    n[0] = n[1] = n[2] = 0.0;
  
    mTensor2 jInv;

    std::vector<pVertex> verts;
    M_GetVertices (ent,verts);
    std::vector<pEdge> edges;
    M_GetEdges (ent, edges);

 
    std::vector<pVertex> vb;
    M_GetVertices (border,vb);
    std::vector<pEdge> eb;
    M_GetEdges (border, eb);
    
    mVector dus[256];
    GradGeomShapeFunction(u ,  v ,  w, dus);
    size_t i, j;
    for(i = 0; i < verts.size(); ++i) {
      vector<pVertex>::iterator it = 
	std::find(vb.begin(),vb.end(),verts[i]);
      if(it == vb.end()) {
	n -= jInv * dus[i];
      }
    }

    for(j=0;j<edges.size(); j++, i++)
      {
	std::vector<pEdge>::iterator it = std::find (eb.begin(),eb.end(),edges[j]);
	if(it == eb.end())
	  {
	    n -= jInv * dus[i];
	  }
      }


    n.norm();
  }

  double LagrangeMapping::detJacDeriv(double &u, double &v, double &w, double *jDeriv, mVector *dus)
  {
    double jac[3][3];
    deval (u,v,w,jac[0][0],jac[1][0],jac[2][0],
	   jac[0][1],jac[1][1],jac[2][1],
	   jac[0][2],jac[1][2],jac[2][2]);
    GradGeomShapeFunction  (u ,  v ,  w, dus);
    size_t i;
    
    for( i=0; i<knots.size(); i++) {
    
      jDeriv[3*i] = dus[i](0)*(jac[1][1]*jac[2][2] - jac[1][2]*jac[2][1]) 
	- jac[1][0]*(dus[i](1)*jac[2][2] - dus[i](2)*jac[2][1]) 
	+ jac[2][0]*(dus[i](1)*jac[1][2] - jac[1][1]*dus[i](2));

      jDeriv[3*i+1] = jac[0][0]*(dus[i](1)*jac[2][2] - dus[i](2)*jac[2][1])
	- dus[i](0)*(jac[0][1]*jac[2][2] - jac[0][2]*jac[2][1])
	+ jac[2][0]*(jac[0][1]*dus[i](2) - dus[i](1)*jac[0][2]);

      jDeriv[3*i+2] = jac[0][0]*(jac[1][1]*dus[i](2) - jac[1][2]*dus[i](1))
	- jac[1][0]*(jac[0][1]*dus[i](2) - jac[0][2]*dus[i](1))
	+ dus[i](0)*(jac[0][1]*jac[1][2] - jac[0][2]*jac[1][1]);
    }
    
    double detJac;
    detJac = jac[0][0]*jac[1][1]*jac[2][2] + jac[0][2] *jac[1][0]*jac[2][1] +
	    jac[0][1]*jac[1][2]*jac[2][0] - jac[0][2] *jac[1][1]*jac[2][0] -
	    jac[0][0]*jac[1][2]*jac[2][1] - jac[0][1] *jac[1][0]*jac[2][2];
    
    return detJac;
    
  } 


  ConstantLagrangeMapping::ConstantLagrangeMapping(pEntity m)
    : LagrangeMapping(m)
  {
    det = LagrangeMapping::jacInverse(0,0,0,jac);
    if(det < 0)det = -det;
  }

  CylindricalCoordinatesLagrangeMapping::CylindricalCoordinatesLagrangeMapping(pEntity m)
    : LagrangeMapping(m)
  {  
  }

  double ConstantLagrangeMapping::detJac(double,double,double) const
  {
    return det;
  }

  double ConstantLagrangeMapping::jacInverse(double,double,double,mTensor2 &j) const
  {
    j = jac;
    return det;
  }

  double CylindricalCoordinatesLagrangeMapping::detJac(double u,double v,double w) const
  {
    mTensor2 j;
    double r,z,teta;
    eval(u,v,w,r,z,teta);
    //  printf("r = %12.5E det = %12.5E\n",r,LagrangeMapping::jacInverse(u,v,w,j) );
    return (r * LagrangeMapping::jacInverse(u,v,w,j));
  }

  double CylindricalCoordinatesLagrangeMapping::jacInverse(double u,double v,double w,mTensor2 &j) const
  {
    double r,z,teta;
    eval(u,v,w,r,z,teta);
    double det = LagrangeMapping::jacInverse (u,v,w,j);
    return (r * det);
  }

  int CylindricalCoordinatesLagrangeMapping::order() const
  { 
    return 2; 
  }


  RegularCubeLagrangeMapping::RegularCubeLagrangeMapping(pEntity e)
    : LagrangeMapping(e)
  {
    std::vector<pVertex> v;
    M_GetVertices (e,v);

    double x1,y1,z1;
    for (size_t i=0; i<v.size(); ++i) {
	pPoint pt = V_point(v[i]); 
	mPoint p = mPoint(P_x(pt),P_y(pt),P_z(pt)); 
	if(!i)
	  {
	    x0 = x1 = p(0);
	    y0 = y1 = p(1);
	    z0 = z1 = p(2);
	  }
	else
	  {
	    if(x0 > p(0))x0 = p(0);
	    if(y0 > p(1))y0 = p(1);
	    if(z0 > p(2))z0 = p(2);
	    if(x1 < p(0))x1 = p(0);
	    if(y1 < p(1))y1 = p(1);
	    if(z1 < p(2))z1 = p(2);
	  }
      }
    dx = x1-x0;
    dy = y1-y0;

    //  printf("%f %f %f %f\n",x0,y0,dx,dy);

    if(entityType == QUAD)
      dz = 1.0;
    else
      dz = z1-z0;
  }

  void RegularCubeLagrangeMapping::eval(double u, double v, double w,
					double &x, double &y, double &z) const
  {
    x = x0 + .5 * (1.+u) * dx;
    y = y0 + .5 * (1.+v) * dy;
    if(entityType == QUAD)
      z = 0.0;
    else
      z = z0 + .5 * (1.+w) * dz;
  }

  void RegularCubeLagrangeMapping::deval(double u, double v, double w,
					 double &dxdu, double &dydu, double &dzdu,
					 double &dxdv, double &dydv, double &dzdv,
					 double &dxdw, double &dydw, double &dzdw) const
  {
    dxdu = 0.5 * dx;
    dydv = 0.5 * dy;
    if(entityType == QUAD)
      dzdw = 0.0;
    else
      dzdw = .5 * dz;
    dxdv=dxdw=dydu=dydw=dzdu=dzdv = 0.0;

  }

  bool RegularCubeLagrangeMapping::invert(double x, double y, double z,
					  double &u, double &v, double &w) const
  {
    u = 2.* (x-x0)/dx -1.0;
    v = 2.* (y-y0)/dy -1.0;
    if(entityType == QUAD)
      w = 0.0;
    else
      w = 2.* (z-z0)/dz -1.0;

    //  printf("%f %f -> %f %f\n",x,y,u,v);

    return true;
  }

  double RegularCubeLagrangeMapping::detJac(double u, double v, double w) const
  {
    if(entityType == QUAD)
      return 0.25*(dx*dy);
    else
      return 0.125*(dx*dy*dz);
  }

  double RegularCubeLagrangeMapping::jacInverse(double u, double v, double w, 
						mTensor2& jInv) const
  {    
    jInv(0,0) = 2./dx;
    jInv(1,1) = 2./dy;
    if(entityType == QUAD)
      jInv(2,2) = 1.0;
    else
      jInv(2,2) = 2./dz;

    jInv(1,2) = jInv(2,1) 
      = jInv(2,0) = jInv(0,2) 
      = jInv(0,1) = jInv(1,0) = 0.0;
    return detJac(0,0,0);
  }

  void RegularCubeLagrangeMapping::normalVector (pEntity border , double u, double v, double w, mVector &n) const
  {
    double u1,v1,w1;
    mPoint p1;
    LagrangeMapping lm(border);
    lm.COG(u1,v1,w1);
    lm.eval(u1,v1,w1,p1(0),p1(1),p1(2));

    if( abs(p1(0) - x0) < 1.e-6 * dx ) n = mVector(-1,0,0);
    else if( abs(p1(0) - x0 - dx) < 1.e-6 * dx ) n = mVector(1,0,0);
    else if( abs(p1(1) - y0) < 1.e-6 * dy ) n = mVector(0,-1,0);
    else if( abs(p1(1) - y0 - dy) < 1.e-6 * dy ) n = mVector(0,1,0);
    else if( abs(p1(2) - z0) < 1.e-6 * dz ) n = mVector(0,0,-1);
    else if( abs(p1(2) - z0 - dz) < 1.e-6 * dz ) n = mVector(0,0,1);
    else assert (1==0);
  }

  double RegularCubeLagrangeMapping::PushBack (double u, double v, double w, int vsize, vector<mVector> &gr) const
  {
    double detJ = detJac(u,v,w);
  
    double xx = 2./dx;
    double yy = 2./dy;
    double zz = 2./dz;

    for(int i=0;i<vsize;i++)
      {
	gr[i](0) *= (xx);
	gr[i](1) *= (yy);
	if(entityType == HEX)
	  gr[i](2) *= (zz);	
      }
    return detJ;
  }

  void LagrangeMapping::boundingBox (  mPoint &pMin, mPoint &pMax ) const
  {
    std::vector<pVertex> v;
    M_GetVertices (ent,v);
    pPoint pt = V_point(v[0]); 
    pMin = pMax = mPoint(P_x(pt),P_y(pt),P_z(pt));    
    int size = v.size();;
    for(int i=1;i<size;i++)
      { 
	pt = V_point(v[i]); 
	mPoint p1 = mPoint(P_x(pt),P_y(pt),P_z(pt)); 
	if (p1(0) < pMin(0)) pMin(0) = p1(0);
	if (p1(1) < pMin(1)) pMin(1) = p1(1);
	if (p1(2) < pMin(2)) pMin(2) = p1(2);
	if (p1(0) > pMax(0)) pMax(0) = p1(0);
	if (p1(1) > pMax(1)) pMax(1) = p1(1);
	if (p1(2) > pMax(2)) pMax(2) = p1(2);
      }
  }
  
  void LagrangeMapping::GeomShapeFunction (double u, double v, double w, double *x) const
  {
    switch(entityType)
      {
      case VERTEX : x[0] = 1.0;break;
      case EDGE   : GeomShapeFunctionLine(u,v,w,x);break;
      case TRI    : GeomShapeFunctionTri(u,v,w,x);break;
      case QUAD   : GeomShapeFunctionQuad(u,v,w,x);break;
      case TET    : GeomShapeFunctionTet(u,v,w,x);break;
      case HEX    : GeomShapeFunctionHex(u,v,w,x);break;
      case PRISM  : GeomShapeFunctionPrism(u,v,w,x);break;
      case PYRAMID  : GeomShapeFunctionPyr(u,v,w,x);break;	
      default: 
	cout << "GeomshapeFunction: Entity type not supported\n"; 
      }
  }
  void LagrangeMapping::GradGeomShapeFunction (double u, double v, double w,
					       mVector* grad) const
  {
    switch(entityType)
      {
      case VERTEX : grad[0][0] = grad[0][1] = grad[0][2] = 0.0;break;
      case EDGE   : GradGeomShapeFunctionLine (u,v,w,grad);break;
      case TRI    : GradGeomShapeFunctionTri  (u,v,w,grad);break;
      case QUAD   : GradGeomShapeFunctionQuad (u,v,w,grad);break;
      case TET    : GradGeomShapeFunctionTet (u,v,w,grad);break;
      case HEX    : GradGeomShapeFunctionHex (u,v,w,grad);break;
      case PRISM  : GradGeomShapeFunctionPrism (u,v,w,grad);break;
      case PYRAMID  : GradGeomShapeFunctionPyr (u,v,w,grad);break;
      default: break; // need to throw an exception here
      }
  }
