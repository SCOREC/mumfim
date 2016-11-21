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
#ifndef SIM
#include "AOMDfwd.h"
#include "AOMD_Defs.h"
#include "BezierMapping.h"
#include "mEntity.h"
#include "mEdge.h"
#include "mFace.h"
#include "mPoint.h"
#include "mVertex.h"
#include "mAOMD.h"

#include <math.h>
#include <stdio.h>

using std::vector;

BezierMapping::BezierMapping (AOMD::mEntity *e, int order)
  :Mapping(e)
{
  entGetKnots();
  //printf("bezier mapping with %d knots\n",knots.size());
  for(unsigned int i=0; i<knots.size(); i++){
    //double f = BezierBase(i, u, v, w);
    mPoint p = knots[i];
    // printf("point (%d) = %f %f %f\n",i,p(0),p(1), p(2));
  }
  Bezierorder= order;
  //printf("The Bezierorder is: %d : BezierMapping \n", Bezierorder);
}
  
int BezierMapping::order() const
{
  //printf("The Bezierorder is: %d : Order()\n", Bezierorder);
  return 2;
}
  
int BezierMapping::geomOrder() const
{
  //printf("The order is: %d  : geomOrder\n", Bezierorder);
  return 2;
}
  
  
void BezierMapping::eval(double u, double v, double w, 
			   double &x, double &y, double &z) const
{
  x = 0.0, y = 0.0, z = 0.0;

  //  printf("evaluation of a bezier mapping with %d knots at %f %f\n",knots.size(),u,v);
  
  for(unsigned int i=0; i<knots.size(); i++){
    double f = BezierBase(i, u, v, w);
    // printf("%d, u=%f, v=%f, w=%f, f=%f", i, u, v, w,f);
    mPoint p = knots[i];
    // printf("The point is: x=%f, y=%f, z=%f", p(0), p(1), p(2));
    x += p(0)*f;
    y += p(1)*f;
    z += p(2)*f;
    //printf("x=%f, y=%f, z=%f\n", x, y, z);
  }

  //printf("\n\n");
}

void BezierMapping::boundingBox(mPoint &pMin, mPoint &pMax) const
{
  pMin = pMax = knots[0];
  for(unsigned int i=1; i<knots.size(); i++){
    mPoint p1 = knots[i];
    if (p1(0) < pMin(0)) pMin(0) = p1(0);
    if (p1(1) < pMin(1)) pMin(1) = p1(1);
    if (p1(2) < pMin(2)) pMin(2) = p1(2);
    if (p1(0) > pMax(0)) pMax(0) = p1(0);
    if (p1(1) > pMax(1)) pMax(1) = p1(1);
    if (p1(2) > pMax(2)) pMax(2) = p1(2);
  }
}

void BezierMapping::deval(double u, double v, double w,
			  double &dxdu, double &dydu, double &dzdu,
			  double &dxdv, double &dydv, double &dzdv,
			  double &dxdw, double &dydw, double &dzdw) const
{
  
  dxdu = 0.0, dydu = 0.0, dzdu = 0.0;
  dxdv = 0.0, dydv = 0.0, dzdv = 0.0;
  dxdw = 0.0, dydw = 0.0, dzdw = 0.0;
  //double t = 1.-u-v-w;

  mVector du;
  for(unsigned int i=0; i<knots.size(); i++)
    {
      //printf("The knots location is: x=%f, y=%f, z=%f\n", knots[i](0), knots[i](1), knots[i](2));
      GradBezierShapeFunction  (i , u ,  v ,  w, du);
      //printf("%d, u=%f, v=%f, w=%f, t=%f\n du[0]=%f, du[1]=%f, du[2]=%f\n", i, u, v, w, t, du[0], du[1], du[2]);
      dxdu += knots[i](0) * du[0];
      dydu += knots[i](1) * du[0];
      dzdu += knots[i](2) * du[0];

      dxdv += knots[i](0) * du[1];
      dydv += knots[i](1) * du[1];
      dzdv += knots[i](2) * du[1];

      dxdw += knots[i](0) * du[2];
      dydw += knots[i](1) * du[2];
      dzdw += knots[i](2) * du[2];

      //  printf("dxdu = %f,   dxdv = %f,   dxdw = %f   \n", dxdu, dxdv, dxdw);
      //        printf("dydu = %f,   dydv = %f,   dydw = %f   \n", dydu, dydv, dydw);
      //        printf("dzdu = %f,   dzdv = %f,   dzdw = %f   \n", dzdu, dzdv, dzdw);
    }

}

void BezierMapping::entGetKnots()
{
  switch(ent->getType())
    {
    case AOMD::mEntity::VERTEX : knots.push_back(((AOMD::mVertex*)ent)->point()); break;
    case AOMD::mEntity::EDGE   : getEdgeKnots(); break;
    case AOMD::mEntity::TRI    : getTriFaceKnots(); break;
    case AOMD::mEntity::TET    : getTetKnots(); break;
    default:
      throw 1;
      break; // need to throw an exception here
    } 
}

void BezierMapping::getEdgeKnots() 
{
  AOMD::mEdge *e = (AOMD::mEdge*)ent; 
  knots.push_back(e->vertex(0)->point());
  
  AOMD::mAttachablePointVector *ap = (AOMD::mAttachablePointVector*)ent->getData(AOMD::AOMD_Util::Instance()->getEmod());
  if(!ap){
    mPoint p0 = e->vertex(0)->point();
    mPoint p1 = e->vertex(1)->point();
    int n = order(); 
    for(int i=0; i<n-1; i++){
      mPoint p;
      for(int j=0; j<3; j++)
	p(j) = ((n-1-i)*p0(j) + (i+1)*p1(j))/(double)(n);
      knots.push_back(p);
    }
  }else{
    std::vector<mPoint>::iterator viter, last;
    viter = (ap->v).begin(); last = (ap->v).end();
    for(; viter!= last; viter++)
      knots.push_back(*viter);
  }
  knots.push_back(e->vertex(1)->point());
}

void BezierMapping::getTriFaceKnots()
{
  AOMD::mAttachablePointVector *ap;
  AOMD::mVertex* vert[3];
  vert[0] = (AOMD::mVertex*)ent->get(0, 0);
  knots.push_back(vert[0]->point());
  vert[1] = (AOMD::mVertex*)ent->get(0, 1);
  knots.push_back(vert[1]->point());
  vert[2] = (AOMD::mVertex*)ent->get(0, 2);
  knots.push_back(vert[2]->point());
  
  AOMD::mEdge* e[3] = {0,0,0};
  AOMD::mEdge temp1(vert[0],vert[1],0);
  AOMD::mEdge temp2(vert[0],vert[2],0);
  AOMD::mEdge temp3(vert[1],vert[2],0);


  for(int i=0;i<3;i++)
    {
      AOMD::mEdge *ee = (AOMD::mEdge*)ent->get(1,i);
      if(ee->equal(&temp1))e[0] = ee;
      if(ee->equal(&temp2))e[1] = ee;
      if(ee->equal(&temp3))e[2] = ee;
    }
  if(!e[0] || !e[1] || !e[2])throw 1;
{  
  for(int i=0; i<3; i++)
    {
      //printf("booh 2\n");
      ap = (AOMD::mAttachablePointVector*)e[i]->getData(AOMD::AOMD_Util::Instance()->getEmod());
      if(!ap)
	{
	  // printf("can not get the emode\n");
	  mPoint p0, p1;
	  int nn = order();
	  // printf("nn = %d : getTriKnots()\n",nn);
	  if(i==0)
	    {
	      p0 = vert[0]->point(); p1 = vert[1]->point();
	    }
	  else if(i==1)
	    {
	      p0 = vert[0]->point(); p1 = vert[2]->point();
	    }
	  else if(i==2)
	    {
	      p0 = vert[1]->point(); p1 = vert[2]->point();
	    }

	  for(int k=0; k<nn-1; k++){
	    mPoint p;
	    for(int j=0; j<3; j++)
	      p(j) = ((nn-1-k)*p0(j) + (k+1)*p1(j))/(double)(nn);
	    knots.push_back(p);
	  }
	  // printf("booh 4\n");
	}
      else{
	//printf("booh 3\n");
	int sign = 0;
	if(i==0 || i==1)
	  {
	    if((e[i]->get(0, 0))->equal(vert[0]))
	      sign = 1;
	  }
	else if(i==2)
	  {
	    if((e[i]->get(0, 0))->equal(vert[1]))
	      sign = 1;
	  }
	//printf("The sign is:   %d\n", sign);
	//printf("The size of the knots vector is : %d\n", (ap->v).size());
	if(sign == 1)
	  {
	    std::vector<mPoint>::iterator viter, last;
	    viter = (ap->v).begin(); last = (ap->v).end();
	    for(; viter!= last; viter++)
	      knots.push_back(*viter);
	  }
	else if(sign==0)
	  {
	    std::vector<mPoint>::reverse_iterator viter, last;
	    viter = (ap->v).rbegin(); last = (ap->v).rend();
	    for(; viter!=last; viter++)
	      knots.push_back(*viter);
	  }

      }
    }
}
  ap = (AOMD::mAttachablePointVector*)ent->getData(AOMD::AOMD_Util::Instance()->getFmod());
  if(ap){
    std::vector<mPoint>::iterator viter, last;
    viter = (ap->v).begin(); last = (ap->v).end();
    for(; viter!= last; viter++)
      knots.push_back(*viter);
  } 
}
  

void BezierMapping::getTetKnots() 
{
  AOMD::mAttachablePointVector *ap;
  AOMD::mVertex *vert0 = (AOMD::mVertex*)ent->get(0, 0);
  AOMD::mVertex *vert1 = (AOMD::mVertex*)ent->get(0, 1);
  AOMD::mVertex *vert2 = (AOMD::mVertex*)ent->get(0, 2);
  AOMD::mVertex *vert3 = (AOMD::mVertex*)ent->get(0, 3);

  
  AOMD::mEdge *e0=0,*e1=0,*e2=0,*e3=0,*e4=0,*e5=0;
  AOMD::mEdge temp0(vert0,vert1,0);
  AOMD::mEdge temp1(vert1,vert2,0);
  AOMD::mEdge temp2(vert2,vert0,0);
  AOMD::mEdge temp3(vert0,vert3,0);
  AOMD::mEdge temp4(vert1,vert3,0);
  AOMD::mEdge temp5(vert2,vert3,0);

  
  for(int i=0;i<6;i++)
    {
      AOMD::mEdge *ee = (AOMD::mEdge*)ent->get(1,i);
      if(ee->equal(&temp0))e0 = ee;
      if(ee->equal(&temp1))e1 = ee;
      if(ee->equal(&temp2))e2 = ee;
      if(ee->equal(&temp3))e3 = ee;
      if(ee->equal(&temp4))e4 = ee;
      if(ee->equal(&temp5))e5 = ee;
    }
  
  knots.push_back(vert0->point());

  ap = (AOMD::mAttachablePointVector*)e0->getData(AOMD::AOMD_Util::Instance()->getEmod());
  if(!ap){
    // printf("can not get the attach emod\n");
    mPoint p;
    for(int i=0; i<3; i++)
      p(i) = 0.5*((vert0->point())(i) + (vert1->point())(i));
    knots.push_back(p);
  }else{
     std::vector<mPoint>::iterator viter, last;
     viter = (ap->v).begin(); last = (ap->v).end();
     for(; viter!= last; viter++)
       knots.push_back(*viter);
  }

  knots.push_back(vert1->point());

  ap = (AOMD::mAttachablePointVector*)e2->getData(AOMD::AOMD_Util::Instance()->getEmod());
  if(!ap){
    mPoint p;
    for(int i=0; i<3; i++)
      p(i) = 0.5*((vert0->point())(i) + (vert2->point())(i));
    knots.push_back(p);
  }else{
    std::vector<mPoint>::iterator viter, last;
    viter = (ap->v).begin(); last = (ap->v).end();
    for(; viter!= last; viter++)
      knots.push_back(*viter);
  }    

  ap = (AOMD::mAttachablePointVector*)e1->getData(AOMD::AOMD_Util::Instance()->getEmod());
  if(!ap){
    mPoint p;
    for(int i=0; i<3; i++)
      p(i) = 0.5*((vert1->point())(i) + (vert2->point())(i));
    knots.push_back(p);
  }else{
    std::vector<mPoint>::iterator viter, last;
    viter = (ap->v).begin(); last = (ap->v).end();
    for(; viter!= last; viter++)
      knots.push_back(*viter);
  }
  
  
  knots.push_back(vert2->point());

  ap = (AOMD::mAttachablePointVector*)e3->getData(AOMD::AOMD_Util::Instance()->getEmod());
  if(!ap){
    mPoint p;
    for(int i=0; i<3; i++)
      p(i) = 0.5*((vert0->point())(i) + (vert3->point())(i));
    knots.push_back(p);
  }else{
    std::vector<mPoint>::iterator viter, last;
    viter = (ap->v).begin(); last = (ap->v).end();
    for(; viter!= last; viter++)
      knots.push_back(*viter);
  }
  
  
  ap = (AOMD::mAttachablePointVector*)e4->getData(AOMD::AOMD_Util::Instance()->getEmod());
  if(!ap){
    mPoint p;
    for(int i=0; i<3; i++)
      p(i) = 0.5*((vert1->point())(i) + (vert3->point())(i));
    knots.push_back(p);
  }else{
    std::vector<mPoint>::iterator viter, last;
    viter = (ap->v).begin(); last = (ap->v).end();
    for(; viter!= last; viter++)
      knots.push_back(*viter);
  }
  

  ap = (AOMD::mAttachablePointVector*)e5->getData(AOMD::AOMD_Util::Instance()->getEmod());
  if(!ap){
    mPoint p;
    for(int i=0; i<3; i++)
      p(i) = 0.5*((vert2->point())(i) + (vert3->point())(i));
    knots.push_back(p);
  }else{
    std::vector<mPoint>::iterator viter, last;
    viter = (ap->v).begin(); last = (ap->v).end();
    for(; viter!= last; viter++)
      knots.push_back(*viter);
  }
  

  knots.push_back(vert3->point());
  
}


double BezierMapping::BezierBase(int i, double u, double v, double w) const
{
  
  switch(ent->getType())
    {
    case AOMD::mEntity::VERTEX  : return 1.0;
    case AOMD::mEntity::EDGE    : return edgeBezierBase(i, u, v, w);
    case AOMD::mEntity::TRI     : return triFaceBezierBase(i, u, v, w);
    case AOMD::mEntity::TET     : return tetRegionBezierBase(i, u, v, w);
    default:
      throw 1;
      break;// need to throw an exception here
    }
}

double BezierMapping::edgeBezierBase(int i, double u, double v, double w) const
{
  int degree = order();
  
 
  if(degree == 2){
    switch(i)
      {
      case 0   : return 0.25*(u-1)*(u-1);
      case 1   : return 0.5*(1-u*u);
      case 2   : return 0.25*(u+1)*(u+1);
      default:
	throw 1;
	return 0.0;
	break;// need to throw an exception here
      }
  }else if(degree == 3){
    switch(i)
      {
      case 0   : return (1-u)*(1-u)*(1-u)/8.;
      case 1   : return 3.*(u+1)*(1-u)*(1-u)/8.;
      case 2   : return 3.*(u+1)*(u+1)*(1-u)/8.;
      case 3   : return (u+1)*(u+1)*(u+1)/8.;
      default:
	throw 1;
	return 0.0;
	break;// need to throw an exception here
      }
  }else if (degree == 4){
    switch(i)
      {
      case 0   : return (1-u)*(1-u)*(1-u)*(1-u)/16.;
      case 1   : return 4.*(1+u)*(1-u)*(1-u)*(1-u)/16.;
      case 2   : return 6.*(1+u)*(1+u)*(1-u)*(1-u)/16.;
      case 3   : return 4.*(1+u)*(1+u)*(1+u)*(1-u)/16.;
      case 4   : return (1+u)*(1+u)*(1+u)*(1+u)/16.;
      default:
	throw 1;
	return 0.0;
	break;// need to throw an exception here
      }
  }
  
  return 0.0;
}

double BezierMapping::triFaceBezierBase(int i, double u, double v, double t) const
{
  int degree = order();
  double w = 1.-u-v;

  if(degree == 2)
    {
      switch(i)
	{
	case 0  : return w*w;
	case 1  : return u*u;
	case 2  : return v*v;
	case 3  : return 2.*w*u;
	case 4  : return 2.*w*v;
	case 5  : return 2.*u*v;
	default:
	  throw 1;
	  return 0.0;
	  break;// need to throw an exception here
	}
    }
  else if(degree == 3)
    {
      switch(i)
	{
	case 0  :   return w*w*w;
	case 1  :   return u*u*u;
	case 2  :   return v*v*v;
	case 3  :   return 3.*w*w*u;
	case 4  :   return 3.*w*u*u;
	case 5  :   return 3.*w*w*v;
	case 6  :   return 3.*w*v*v;
	case 7  :   return 3.*u*u*v;
	case 8  :   return 3.*u*v*v;
	case 9  :   return 6.*v*u*w;
	}
    }
  else if(degree == 4)
    {
      switch(i)
	{
	case 0  :   return w*w*w*w;
	case 1  :   return u*u*u*u;
	case 2  :   return v*v*v*v;
	case 3  :   return 4.*w*w*w*u;
	case 4  :   return 6.*w*w*u*u;
	case 5  :   return 4.*w*u*u*u;
	case 6  :   return 4.*w*w*w*v;
	case 7  :   return 6.*w*w*v*v;
	case 8  :   return 4.*w*v*v*v;
	case 9  :   return 4.*u*u*u*v;
	case 10  :   return 6.*u*u*v*v;
	case 11  :   return 4.*u*v*v*v;
	case 12  :   return 12.*w*w*u*v;
	case 13  :   return 12.*w*u*u*v;
	case 14  :   return 12.*w*u*v*v;
	}
    }

  return 0.0;
}


double BezierMapping::tetRegionBezierBase(int i, double u, double v, double w) const
{
  double t = 1.-u-v-w;
  
  switch(i)
    {
    case 0  : return t*t;
    case 1  : return 2.*t*u;
    case 2  : return u*u;
    case 3  : return 2.*t*v;
    case 4  : return 2.*u*v;
    case 5  : return v*v;
    case 6  : return 2*t*w;
    case 7  : return 2*u*w;
    case 8  : return 2*v*w;
    case 9  : return w*w;
    default:
      throw 1;
      return 0.0;
      break;// need to throw an exception here
    }

  return 0.0;
}


void BezierMapping::GradBezierShapeFunction  (int iNod, double u,  double v,  double w, mVector &grad) const
{
  switch(ent->getType())
    {
    case AOMD::mEntity::VERTEX  : grad[0]=grad[1]=grad[2]=0.0; break;
    case AOMD::mEntity::EDGE    : GradBezierShapeFunctionLine(iNod,u,v,w,grad); break;
    case AOMD::mEntity::TRI     : GradBezierShapeFunctionTri(iNod,u,v,w,grad); break;
    case AOMD::mEntity::TET     : GradBezierShapeFunctionTet(iNod,u,v,w,grad); break;
    default:break; // need to throw an exception here
    }
}

void BezierMapping::GradBezierShapeFunctionLine(int iNod, double u, double v, double w, mVector &grad) const
{
  grad[1] = grad[2] = 0.0;
  int degree = order();

  if(degree == 2)
    {
      switch(iNod)
	{
	case 0: grad[0] = 0.5*(u-1); break;
	case 1: grad[0] = -u; break;
	case 2: grad[0] = 0.5*(u+1); break;
	default: break; //need to throw an exception here
	}
    }
  else if(degree == 3)
    {
      switch(iNod)
	{
	case 0 : grad[0] = -3.*(1.-u)*(1.-u)/8.;    break;
	case 1 : grad[0] = 3.*(-1.-2.*u+3.*u*u)/8.; break;
	case 2 : grad[0] = 3.*(-3.*u*u+1.-2.*u)/8.; break;
	case 3 : grad[0] = 3.*(1+u)*(1+u)/8.;        break;
	default:break;
	}
    }
  else if(degree == 4)
    {
      switch(iNod)
	{
	case 0 : grad[0] = 4.*(u-1.)*(u-1.)*(u-1.)/16.;   break;
	case 1 : grad[0] = 4.*(-2.+6.*u*u-4.*u*u*u)/16.;  break;
	case 2 : grad[0] = 6.*(-4.*u+4.*u*u*u)/16.;       break;
	case 3 : grad[0] = 4.*(-6.*u*u-4.*u*u*u+2.)/16.;  break;
	case 4 : grad[0] = 4.*(1.+u)*(1.+u)*(1.+u)/16.;   break;
	}
    }
}

void BezierMapping::GradBezierShapeFunctionTri(int iNod, double u, double v, double dum, mVector &grad) const
{
  double w=1.-u-v;
  grad(2) = 0.0;
  int degree = order();

  if(degree == 2)
    {
      switch(iNod)
	{
	case 0 :  grad[0]=-2.*w;    grad[1]=-2.*w;    break;
	case 1 :  grad[0]=2.*u;     grad[1]=0.0;      break;
	case 2 :  grad[0]=0.0;      grad[1]=2.*v;     break;
	case 3 :  grad[0]=2.*(w-u); grad[1]=-2.*u;    break;
	case 4 :  grad[0]=-2.*v;    grad[1]=2.*(w-v); break;
	case 5 :  grad[0]=2.*v;     grad[1]=2.*u;     break;
	default: break; // need to throw an exception here    
	}
    }
  else if(degree == 3)
    {
      switch(iNod)
	{
	case 0 : grad[0]=-3.*w*w;          grad[1]=-3.*w*w;         break;
	case 1 : grad[0]=3.*u*u;           grad[1]=0.0;             break;
	case 2 : grad[0]=0.0;              grad[1]=3.*v*v;            break;
	case 3 : grad[0]=3.*w*w-6.*w*u;    grad[1]=-6.*w*u;         break;
	case 4 : grad[0]=6.*w*u-3.*u*u;    grad[1]=-3.*u*u;         break;
	case 5 : grad[0]=-6.*w*v;          grad[1]=3.*w*w-6.*w*v;   break;
	case 6 : grad[0]=-3.*v*v;          grad[1]=6.*w*v-3.*v*v;   break;
	case 7 : grad[0]=6.*u*v;           grad[1]=3.*u*u;          break;
	case 8 : grad[0]=3.*v*v;           grad[1]=6.*u*v;          break;
	case 9 : grad[0]=6.*w*v-6.*u*v;    grad[1]=6.*w*u-6.*u*v;   break; 
	}
      
    }
  else if(degree == 4)
    {
      switch(iNod)
	{
	case 0 : grad[0]=-4.*w*w*w;            grad[1]=-4.*w*w*w;             break;
	case 1 : grad[0]=4.*u*u*u;             grad[1]=0.;                    break; 
	case 2 : grad[0]=0.0;                  grad[1]=4.*v*v*v;              break;
	case 3 : grad[0]=4.*w*w*w-12.*w*w*u;   grad[1]=-12.*w*w*u;            break; 
	case 4 : grad[0]=12.*w*w*u-12.*w*u*u;  grad[1]=-12.*w*u*u;            break;
	case 5 : grad[0]=12.*w*u*u-4.*u*u*u;   grad[1]=-4.*u*u*u;             break; 
	case 6 : grad[0]=-12.*w*w*v;           grad[1]=4.*w*w*w-12.*w*w*v;    break;
	case 7 : grad[0]=-12.*w*v*v;           grad[1]=12.*w*w*v-12.*w*v*v;   break;
	case 8 : grad[0]=-4.*v*v*v;            grad[1]=12.*w*v*v-4.*v*v*v;    break;
	case 9 : grad[0]=12.*u*u*v;            grad[1]=4.*u*u*u;              break; 
	case 10 : grad[0]=12.*u*v*v;           grad[1]=12.*u*u*v;             break;
	case 11 : grad[0]=4.*v*v*v;            grad[1]=12.*u*v*v;             break; 
	case 12 : grad[0]=12.*w*w*v-24.*u*v*w; grad[1]=12.*w*w*u-24.*u*v*w;   break;
	case 13 : grad[0]=24.*u*v*w-12.*u*u*v; grad[1]=12.*w*u*u-12.*u*u*v;   break; 
	case 14 : grad[0]=12.*w*v*v-12.*u*v*v; grad[1]=24.*u*v*w-12.*u*v*v;   break;
	}
    }
}




void BezierMapping::GradBezierShapeFunctionTet(int iNod, double u, double v, double w, mVector &grad) const
{
  double t=1.-u-v-w;

  switch(iNod)
    {
    case 0: grad[0]=-2.*t;    grad[1]=-2.*t;     grad[2]=-2.*t;     break; 
    case 1: grad[0]=2.*(t-u); grad[1]=-2.*u;     grad[2]=-2.*u;     break;
    case 2: grad[0]=2.*u;     grad[1]=0.0;       grad[2]=0.0;       break;
    case 3: grad[0]=-2.*v;    grad[1]=2.*(t-v);  grad[2]=-2.*v;     break;
    case 4: grad[0]=2.*v;     grad[1]=2.*u;      grad[2]=0.0;       break;
    case 5: grad[0]=0.0;      grad[1]=2.*v;      grad[2]=0.0;       break;
    case 6: grad[0]=-2.*w;    grad[1]=-2.*w;     grad[2]=2.*(t-w);  break;
    case 7: grad[0]=2.*w;     grad[1]=0.0;       grad[2]=2.*u;      break;
    case 8: grad[0]=0.0;      grad[1]=2.*w;      grad[2]=2.*v;      break;
    case 9: grad[0]=0.0;      grad[1]=0.0;       grad[2]=2.*w;      break;  
    default: break; // need to throw an exception here    
    }
}


void BezierMapping::normalVector(AOMD::mEntity *border,double u, double v, double w, mVector &n) const
{
  switch(ent->getType()){
  case AOMD::mEntity::TRI:
    triNormalVector(border, u, v, w, n);
    break;
  case AOMD::mEntity::TET:
    if(border->getType() == AOMD::mEntity::EDGE)
      tetEdgeNormalVector(border, u, v, w, n);
    else if(border->getType() == AOMD::mEntity::TRI)
      tetFaceNormalVector(border, u, v, w, n);
    break;
  default:
    throw 1;
    break;
  }
}

void BezierMapping::triNormalVector(AOMD::mEntity *border,double u, double v, double w, mVector &n) const
{
  n[0] = n[1] = n[2] = 0.0;
  
  mTensor2 jInv;
  jacInverse(u, v, w, jInv);
  AOMD::mVertex *vert1 = (AOMD::mVertex*)ent->get(0, 0);
  AOMD::mVertex *vert2 = (AOMD::mVertex*)ent->get(0, 1);
  AOMD::mVertex *vert3 = (AOMD::mVertex*)ent->get(0, 2);
  
  AOMD::mEdge temp1(vert1,vert2,0);
  AOMD::mEdge temp2(vert1,vert3,0);
  AOMD::mEdge temp3(vert2,vert3,0);
  
  int *index, degree=order(), count;
  count = degree+1;
  index = new int[count];
  if(temp1.equal((AOMD::mEdge*)border))
    {
      if(degree == 2)
	{index[0] = 0; index[1]=3; index[2]=1;}
      else if(degree == 3)
	{index[0] = 0; index[1]=3; index[2]=4; index[3]=1;}
      else if(degree == 4)
	{index[0] = 0; index[1]=3; index[2]=4; index[3]=5; index[4]=1;}
    }
  else if(temp2.equal((AOMD::mEdge*)border))
    {
      if(degree == 2)
	{index[0] = 0; index[1]=4; index[2]=2;}
      else if(degree == 3)
	{index[0] = 0; index[1]=5; index[2]=6; index[3]=2;}
      else if(degree == 4)
	{index[0] = 0; index[1]=6; index[2]=7; index[3]=8; index[4]=2;}
    }
  else if(temp3.equal((AOMD::mEdge*)border))
    {
      if(degree == 2)
	{index[0] = 1; index[1]=5; index[2]=2;}
      else if(degree == 3)
	{index[0] = 1; index[1]=7; index[2]=8; index[3]=2;}
      else if(degree == 4)
	{index[0] = 1; index[1]=9; index[2]=10; index[3]=11; index[4]=2;}
    }
  
  for(int i=0; i<count; i++){
    mVector gr;
    GradBezierShapeFunction(index[i], u, v, w, gr);
    gr = jInv * gr;
    n += gr;
  }
  
  delete []index;
  n.norm();  
}

void BezierMapping::tetEdgeNormalVector(AOMD::mEntity *border,double u, double v, double w, mVector &n) const
{

  n[0] = n[1] = n[2] = 0.0;
  
  mTensor2 jInv;
  jacInverse(u, v, w, jInv);

  AOMD::mVertex *vert0 = (AOMD::mVertex*)ent->get(0, 0);
  AOMD::mVertex *vert1 = (AOMD::mVertex*)ent->get(0, 1);
  AOMD::mVertex *vert2 = (AOMD::mVertex*)ent->get(0, 2);
  AOMD::mVertex *vert3 = (AOMD::mVertex*)ent->get(0, 3);

  AOMD::mEdge temp0(vert0,vert1,0);
  AOMD::mEdge temp1(vert1,vert2,0);
  AOMD::mEdge temp2(vert0,vert2,0);
  AOMD::mEdge temp3(vert0,vert3,0);
  AOMD::mEdge temp4(vert1,vert3,0);
  AOMD::mEdge temp5(vert2,vert3,0);

  int index[3];
  
  if(temp0.equal((AOMD::mEdge*)border)){
    index[0]=0; index[1]=1; index[2]=2;
  }else if(temp1.equal((AOMD::mEdge*)border)){
    index[0]=2; index[1]=4; index[2]=5;
  }else if(temp2.equal((AOMD::mEdge*)border)){
    index[0]=0; index[1]=3; index[2]=5;
  }else if(temp3.equal((AOMD::mEdge*)border)){
    index[0]=0; index[1]=6; index[2]=9;
  }else if(temp4.equal((AOMD::mEdge*)border)){
    index[0]=2; index[1]=7; index[2]=9;
  }else if(temp5.equal((AOMD::mEdge*)border)){
    index[0]=5; index[1]=8; index[2]=9;
  }
  
  for(int i=0; i<3; i++){
    mVector gr;
    GradBezierShapeFunction(index[i], u, v, w, gr);
    gr = jInv * gr;
    n += gr;
  }
  
  n.norm();
}


void BezierMapping::tetFaceNormalVector(AOMD::mEntity *border,double u, double v, double w, mVector &n) const
{

  n[0] = n[1] = n[2] = 0.0;
  
  mTensor2 jInv;
  jacInverse(u, v, w, jInv);

  AOMD::mVertex *vert0 = (AOMD::mVertex*)ent->get(0, 0);
  AOMD::mVertex *vert1 = (AOMD::mVertex*)ent->get(0, 1);
  AOMD::mVertex *vert2 = (AOMD::mVertex*)ent->get(0, 2);
  AOMD::mVertex *vert3 = (AOMD::mVertex*)ent->get(0, 3);

  AOMD::mFace temp0(vert0,vert1,vert2,0);
  AOMD::mFace temp1(vert0,vert1,vert3,0);
  AOMD::mFace temp2(vert1,vert2,vert3,0);
  AOMD::mFace temp3(vert0,vert2,vert3,0);

  int index[6];
  if(temp0.equal((AOMD::mFace*)border)){
    index[0]=0; index[1]=1; index[2]=2; 
    index[3]=3; index[4]=4; index[5]=5;
  }else if(temp1.equal((AOMD::mFace*)border)){
    index[0]=0; index[1]=1; index[2]=2; 
    index[3]=6; index[4]=7; index[5]=9;
  }else if(temp2.equal((AOMD::mFace*)border)){
    index[0]=2; index[1]=4; index[2]=5; 
    index[3]=7; index[4]=8; index[5]=9;
  }else if(temp3.equal((AOMD::mFace*)border)){
    index[0]=0; index[1]=3; index[2]=5; 
    index[3]=6; index[4]=8; index[5]=9;
  }

  for(int i=0; i<6; i++){
    mVector gr;
    GradBezierShapeFunction(index[i], u, v, w, gr);
    gr = jInv * gr;
    n += gr;
  }

  // printf("n : %f %f %f\n", n[0], n[1], n[2]);
  
  n.norm();
}
#endif
