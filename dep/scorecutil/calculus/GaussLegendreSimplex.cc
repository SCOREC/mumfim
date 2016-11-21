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
#include <math.h>
#include "IntPt.h"
#include "GaussLegendreSimplex.h"

void brickToTet(double xi,double eta, double zeta,
                double *r, double *s, double *t, double *J);

void quadToTri(double xi,double eta,
               double *r, double *s,double *J);
double quadToTriJac(double xi,double eta);
int GaussLegendre1D(int, double **,double **);

int GaussLegendreTet(int n1,int n2,int n3,IntPt3d *pts) 
{
  /* get degenerate n1Xn2Xn3 Gauss-Legendre scheme to integrate over a tet */
  int i,j,k,index=0;
  //int npnt = n1 * n2 * n3 ;
  double *pt1,*pt2,*pt3,*wt1,*wt2,*wt3,dJ;
  //const double six=6.000000000000000;

  GaussLegendre1D(n1,&pt1,&wt1);
  GaussLegendre1D(n2,&pt2,&wt2);
  GaussLegendre1D(n3,&pt3,&wt3);
  for(i=0; i < n1; i++) {
    for(j=0; j < n2; j++) {
      for(k=0; k < n3; k++) {
        brickToTet(pt1[i],pt2[j],pt3[k],&(pts[index].pt[0]),
                   &(pts[index].pt[1]),&(pts[index].pt[2]),&dJ);
        //GLrstw[index][3] = 1.0e0-GLrstw[index][0]-GLrstw[index][1]
	  // -GLrstw[index][2];
        //GLwt[index++] = dJ*wt1[i]*wt2[j]*wt3[k]*6.0;
	pts[index++].weight = dJ*wt1[i]*wt2[j]*wt3[k]*6.0;
      }
    }
  }
  return index;
}

//int GaussLegendreTri(int n1,int n2, double GLr[][3],double *GLwt) 
int GaussLegendreTri(int n1,int n2,IntPt2d *pts) 
{
  /* get degenerate n1Xn2 Gauss-Legendre scheme to integrate over a tri */
  int i,j,index=0;
  //int npnt = n1 * n2 ;
  double *pt1,*pt2,*wt1,*wt2,dJ;
  //const double two = 2.0000000000000000;

  GaussLegendre1D(n1,&pt1,&wt1);
  GaussLegendre1D(n2,&pt2,&wt2);
  for(i=0; i < n1; i++) {
    for(j=0; j < n2; j++) {
      //quadToTri(pt1[i],pt2[j],&GLr[index][0],&GLr[index][1],&dJ);
      quadToTri(pt1[i],pt2[j],&(pts[index].pt[0]),&(pts[index].pt[1]),&dJ);
      //GLwt[index++] = dJ*wt1[i]*wt2[j]*two;
      pts[index++].weight = dJ*wt1[i]*wt2[j]*2.0;
    }
  }
  return index;
}

void brickToTet(double xi,double eta, double zeta,
                      double *r, double *s, double *t, double *J) {
  double r1,rs1;
  *r = 0.5e0*(1.0e0+xi);
  r1 = 1.0e0-(*r);
  *s = 0.5e0*(1.0e0+eta)*r1;
  rs1 = 1.0e0-(*r)-(*s);
  *t = 0.5e0*(1.0e0+zeta)*rs1;
  *J = 0.125e0*r1*rs1;
}

void quadToTri(double xi,double eta,double *r, double *s, double *J) {
  double r1;
  *r = 0.5e0*(1.0e0+xi);
  r1 = 1.0e0-(*r);
  *s = 0.5e0*(1.0e0+eta)*r1;
  *J = 0.25e0*r1;  
}

double quadToTriJac(double,double eta) {
  return 0.125e0*(1.0e0-eta);
}



