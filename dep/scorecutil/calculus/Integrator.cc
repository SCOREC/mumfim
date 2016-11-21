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
#include "Integrator.h"
#ifndef SIM
#include "AOMD_Internals.h"
#endif
#include "IntPt.h"
#include "Mapping.h"
#include "MSMapping.h"
#include "GaussQuadrature.h"

int getNGQTPts(int);
int getNGQQPts(int);
int getNGQTetPts(int);
int getNGQHPts(int);
int getNGQPrismPts(int);
int getNGQPyramidPts(int);
int GaussLegendre1D(int,double **, double **);

GaussIntegrator::GaussIntegrator(pEntity e)
  : ent(e)
{
  entityType = M_GetElementType(e);
}

GaussIntegrator::GaussIntegrator(Mapping * theMapping)
{
  ent = theMapping->getEntity();
  entityType = theMapping->getEntityType();
}

int GaussIntegrator::nbIntegrationPoints(int order) const
{
  switch(M_GetElementType(ent))
    {
    case VERTEX : return 1;
    case EDGE   : return order/2+1;//(order<1)?1:order;
    case TRI    : return (order<1)?1:getNGQTPts(order);
    case TET    : return (order<1)?1:getNGQTetPts(order);
    case QUAD   : return getNGQQPts(order);
    case HEX    : return getNGQHPts(order);
    case PRISM  : return getNGQPrismPts(order);
    case PYRAMID : return getNGQPyramidPts(order);
    default: throw 1;
    }
}

void GaussIntegrator::iPoint(int i, int order , double &u, double &v, double &w, double &weight)const
{
  switch(entityType)
    {
    case VERTEX :
      u = v = w = 0.0;
      weight = 1.0;
      break;
    case EDGE :
      {
	double *pt,*wt;
	GaussLegendre1D(nbIntegrationPoints(order),&pt,&wt);
	u = pt[i]; 
	v = w = 0.0;
	weight = wt[i];
      }
      break;
    case TRI :
      {
	IntPt2d *pts = getGQTPts(order);
	u = pts[i].pt[0];
	v = pts[i].pt[1];
	w = 0.0;
	weight = 0.5 * pts[i].weight;
      }
      break;
    case QUAD :
      {
	IntPt2d *pts = getGQQPts(order);
	u = pts[i].pt[0];
	v = pts[i].pt[1];
	w = 0.0;
	weight = pts[i].weight;
      }
      break;
    case HEX :
      {
	IntPt3d *pts = getGQHPts(order);
	u = pts[i].pt[0];
	v = pts[i].pt[1];
	w = pts[i].pt[2];
	weight = pts[i].weight;
      }
      break;
    case TET:
      {
	const double sixth = 1.0/6.0;
	IntPt3d *pts = getGQTetPts(order);
	u = pts[i].pt[0];
	v = pts[i].pt[1];
	w = pts[i].pt[2];
	weight = sixth*pts[i].weight;
      }
      break;
     /// added from Luo's scorecutil --CWS 3/10/2008
    case PRISM :
      {
        IntPt3d *pts = getGQPrismPts(order);
        u = pts[i].pt[0];
        v = pts[i].pt[1];
        w = pts[i].pt[2];
        weight = pts[i].weight;
      }
      break;
    case PYRAMID:
      {
        IntPt3d *pts = getGQPyramidPts(order);
        u = pts[i].pt[0];
        v = pts[i].pt[1];
        w = pts[i].pt[2];
        weight = pts[i].weight;
      }
      break;
    default:
      assert(0);
      break;
    }
}

