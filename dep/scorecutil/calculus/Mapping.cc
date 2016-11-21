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
#include "Mapping.h"
#ifdef AOMD
#include "AOMD_Internals.h"
#include "mVertex.h"
#endif
#include "MSMapping.h"
#include "mPoint.h"

#include <cmath>
#include <cstdio>
#include <cassert>


using std::cout;
using std::vector;
using std::abs;



  Mapping::Mapping(pEntity e) : knots()
   /*
    : entityType(M_GetElementType(e)), 
      entityDim(EN_type(e)), 
      rootMapping(this), 
      ent(e)
  {
  }
   */
   // Luo's --CWS 3/10/2008
 {
    if(e) {
      entityType = M_GetElementType(e);
      entityDim = EN_type(e);
      rootMapping = this;
      ent = e;
      m_order = 1;
    }
  }

  
  Mapping::Mapping(Mapping * theMapping, 
		   std::vector<mPoint> *_input_knots, 
		   TRELLIS_MTYPE _type, int _dim)
    : entityType(_type), entityDim(_dim), ent(theMapping->getEntity())
  {
  }

  void Mapping::buildMappings(Mapping * theMapping, 
			      LevelSetFunction _levelSetFunc,
			      int flag,
			      std::list<Mapping*> &out_list )
  {}

  void Mapping::local2rootlocal(double u, double v, double w, double &r_u, double &r_v, double &r_w) {
    double x, y, z;
    if(rootMapping != this) {    
      eval(u,v, w, x, y, z);
      rootMapping->invert(x, y, z, r_u, r_v, r_w);
    } else {
      r_u = u;
      r_v = v;
      r_w = w;
    }
  }

  double Mapping::PushBack (double u, double v, double w, int vsize, mVector *gr) const {
    mTensor2 jInv;
    double detJac = jacInverse(u,v,w,jInv);
    for(int i = 0; i < vsize; i++) {
      gr[i] = jInv * gr[i];
    }
    return detJac;
  }
  
  double Mapping::PushBack (double u, double v, double w, int vsize, mTensor2 *gr) const
  {
    mTensor2 jInv;
    double detJac = jacInverse(u,v,w,jInv);

    for(int i=0; i<vsize; i++) {
/*
		std::cout << "gr of edge i=" << i;
		for (int j=0; j<3; ++j)
			std::cout << ' ' << gr[i](0,j);
		std::cout << "\ngs of edge i=" << i;
		for (j=0; j<3; ++j)
			std::cout << ' ' << gr[i](1,j);
		std::cout << "\ngt of edge i=" << i;
		for (j=0; j<3; ++j)
			std::cout << ' ' << gr[i](2,j);
		std::cout << std::endl;
*/
//		gr[i][0] *= jInv; gr[i][1] *= jInv; gr[i][2] *= jInv;

	//	gr[i] is a mTensor2 whose jk element is kth component of mVector Ni,j
	//	gr[i][0] = Ni,r,	gr[i][1] = Ni,s,	gr[i][2] = Ni,t
	//	jInv is the inverse of the Jacobian
	//		|r,x s,x t,x|
	//		|r,y s,y t,y|
	//		|r,x s,z t,z|
		gr[i] = jInv * gr[i];
	//	After above matrix multiplication:
	//	gr[i][0] = Ni,x, gr[i][1] = Ni,y, gr[i][2] = Ni,z

/*		std::cout << "gx of edge i=" << i;
		for (j=0; j<3; ++j)
			std::cout << ' ' << gr[i](0,j);
		std::cout << "\ngy of edge i=" << i;
		for (j=0; j<3; ++j)
			std::cout << ' ' << gr[i](1,j);
		std::cout << "\ngz of edge i=" << i;
		for (j=0; j<3; ++j)
			std::cout << ' ' << gr[i](2,j);
		std::cout << std::endl << std::endl;
*/
    }
    return detJac;
  }

#if 0
  bool Mapping::invert(double xp, double yp, double zp, double &Upos, double &Vpos, double &Wpos) const
  {
#define NR_PRECISION       1.e-6
#define NR_MAX_ITER        50
	
    double   x_est, y_est, z_est;
    double   u_new, v_new, w_new;
    double   Error = 1.0 ;
    int      iter = 1 ;

    COG(Upos,Vpos,Wpos);
	
    mTensor2 InvJacMatrix;

    while (Error > NR_PRECISION && iter < NR_MAX_ITER){
		
      iter++ ;

      jacInverse(Upos,Vpos,Wpos,InvJacMatrix);
      eval(Upos,Vpos,Wpos,x_est,y_est,z_est);

      //printf("%f %f %f %f %f %f\n",Upos,Vpos,Wpos,x_est,y_est,z_est);
      
      u_new = Upos + InvJacMatrix(0,0) * (xp-x_est) + InvJacMatrix(1,0) * (yp-y_est) +
	InvJacMatrix(2,0) * (zp-z_est) ;
      v_new = Vpos + InvJacMatrix(0,1) * (xp-x_est) + InvJacMatrix(1,1) * (yp-y_est) +
	InvJacMatrix(2,1) * (zp-z_est) ;
      w_new = Wpos + InvJacMatrix(0,2) * (xp-x_est) + InvJacMatrix(1,2) * (yp-y_est) +
	InvJacMatrix(2,2) * (zp-z_est) ;
    
      Error = (u_new - Upos) * (u_new - Upos) + 
	(v_new - Vpos) * (v_new - Vpos) + 
	(w_new - Wpos) * (w_new - Wpos) ;
    
      Upos = u_new;
      Vpos = v_new;
      Wpos = w_new;
    }
    if(Error > NR_PRECISION) {
	//printf("impossible to find %f %f %f \n",xp,yp,zp);
	 return false;
    }
    return true;
  }

#endif
  
  //invert mapping for mixed elements up to quadratic order
  void ludcmp_(double **a, int n, int *indx,double *d)
{
  int i, imax, j,k;
  double big,dum,sum,temp;
  vector<double> vv(n);
    
  *d=1.0;
  for(i=0; i<n; i++) {
    big=0.0;
    for(j=0;j<n; j++)
      if((temp = abs(a[i][j])) > big) big = temp;
    if(big==0.0) std::cerr<<"Singular matrix in routine LUDCMP"<<std::endl;
    vv[i]=1.0/big;
  }
  for(j=0; j<n;j++){
    for(i=0;i<j;i++){
      sum=a[i][j];
      for(k=0;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
    }
    big=0.0;
    for(i=j; i<n;i++){
      sum=a[i][j];
      for(k=0;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if((dum = vv[i] * abs(sum)) >= big){
	big=dum;
	imax=i;
      }
    }
    if(j != imax){
      for(k=0; k <n; k++){
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if( a[j][j] == 0.0) a[j][j] = 1.0e-20;
    if(j !=n){
      dum=1.0/(a[j][j]);
      for(i=j+1;i<n;i++) a[i][j] *= dum;
    }
  }
}


void lubksb_(double **a, int n, int *indx, double b[])
{
  int i,ii=-1,ip,j;
  double sum;
  
  for(i=0; i< n; i++){
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if(ii != -1)
      for(j=ii; j <= i-1; j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  for(i=n-1; i>=0; i--){
    sum=b[i];
    for(j=i+1;j<n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}

bool Mapping::invert(double xp, double yp, double zp, double &Upos, double &Vpos, double &Wpos) const
{
    
   
    double edgeA[3][3] = {
      {1.,  -1.,  1.},	{1.,   1.,  1.},{1.,   0.,  0.}};
    
    double triA[6][6] = {
      {1.,   0.,   0.,   0.,   0.,    0.},
      {1.,   1.,   0.,   1.,   0.,    0.},
      {1.,   0.,   1.,   0.,   1.,    0.},
      {1.,  0.5,   0., 0.25,   0.,    0.},
      {1.,  0.5,  0.5, 0.25, 0.25,  0.25},
      {1.,   0.,  0.5,   0., 0.25,    0.}};
    double quadA[8][8] = {
      {1.,  -1.,  -1.,  1.,  1.,  1.,  -1.,  -1.},
      {1.,   1.,  -1., -1.,  1.,  1.,  -1.,   1.},
      {1.,   1.,   1.,  1.,  1.,  1.,   1.,   1.},
      {1.,  -1.,   1., -1.,  1.,  1.,   1.,  -1.},
      {1.,   0.,  -1.,  0.,  0.,  1.,   0.,   0.},
      {1.,   1.,   0.,  0.,  1.,  0.,   0.,   0.},
      {1.,   0.,   1.,  0.,  0.,  1.,   0.,   0.},
      {1.,  -1.,   0.,  0.,  1.,  0.,   0.,   0.}};
    double tetA[10][10] = {
    //0     1     2     3     4     5     6     7     8     9
      {1.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.},
      {1.,   1.,   0.,   0.,   1.,   0.,   0.,   0.,   0.,   0.},
      {1.,   0.,   1.,   0.,   0.,   1.,   0.,   0.,   0.,   0.},
      {1.,   0.,   0.,   1.,   0.,   0.,   1.,   0.,   0.,   0.},
      {1.,  0.5,   0.,   0., 0.25,   0.,   0.,   0.,   0.,   0.},
      {1.,  0.5,  0.5,   0., 0.25, 0.25,   0., 0.25,   0.,   0.},
      {1.,   0.,  0.5,   0.,   0., 0.25,   0.,   0.,   0.,   0.},
      {1.,   0.,   0.,  0.5,   0.,   0., 0.25,   0.,   0.,   0.},
      {1.,  0.5,   0.,  0.5, 0.25,   0., 0.25,   0., 0.25,   0.},
      {1.,   0.,  0.5,  0.5,   0., 0.25, 0.25,   0.,   0., 0.25}};
    double hexA[20][20] = {
      //0    1    2    3    4    5    6    7    8   9  10   11   12   13   14   15   16   17   18  19
      {1., -1., -1., -1.,  1.,  1.,  1., -1.,  1., 1., 1., -1., -1., -1., -1., -1., -1.,  1.,  1.,  1.},
      {1.,  1., -1., -1., -1., -1.,  1.,  1.,  1., 1., 1., -1., -1.,  1., -1.,  1., -1.,  1., -1., -1.},
      {1.,  1.,  1., -1.,  1., -1., -1., -1.,  1., 1., 1.,  1., -1.,  1., -1.,  1.,  1., -1., -1.,  1.},
      {1., -1.,  1., -1., -1.,  1., -1.,  1.,  1., 1., 1.,  1., -1., -1., -1., -1.,  1., -1.,  1., -1.},
      {1., -1., -1.,  1.,  1., -1., -1.,  1.,  1., 1., 1., -1.,  1., -1.,  1., -1., -1., -1., -1.,  1.},
      {1.,  1., -1.,  1., -1.,  1., -1., -1.,  1., 1., 1., -1.,  1.,  1.,  1.,  1., -1., -1.,  1., -1.},
      {1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1., 1., 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.},
      {1., -1.,  1.,  1., -1., -1.,  1., -1.,  1., 1., 1.,  1.,  1., -1.,  1., -1.,  1.,  1., -1., -1.},
      {1.,  0., -1., -1.,  0.,  0.,  1.,  0.,  0., 1., 1.,  0.,  0.,  0., -1.,  0., -1.,  0.,  0., 0.},
      {1.,  1.,  0., -1.,  0., -1.,  0.,  0.,  1., 0., 1.,  0., -1.,  0.,  0.,  1.,  0.,  0.,  0., 0.},
      {1.,  0.,  1., -1.,  0.,  0., -1.,  0.,  0., 1., 1.,  0.,  0.,  0., -1.,  0.,  1.,  0.,  0., 0.},
      {1., -1.,  0., -1.,  0.,  1.,  0.,  0.,  1., 0., 1.,  0., -1.,  0.,  0., -1.,  0.,  0.,  0., 0.},
      {1., -1., -1.,  0.,  1.,  0.,  0.,  0.,  1., 1., 0., -1.,  0., -1.,  0.,  0.,  0.,  0.,  0., 0.},
      {1.,  1., -1.,  0., -1.,  0.,  0.,  0.,  1., 1., 0., -1.,  0.,  1.,  0.,  0.,  0.,  0.,  0., 0.},
      {1.,  1.,  1.,  0.,  1.,  0.,  0.,  0.,  1., 1., 0.,  1.,  0.,  1.,  0.,  0.,  0.,  0.,  0., 0.},
      {1., -1.,  1.,  0., -1.,  0.,  0.,  0.,  1., 1., 0.,  1.,  0., -1.,  0.,  0.,  0.,  0.,  0., 0.},
      {1.,  0., -1.,  1.,  0.,  0., -1.,  0.,  0., 1., 1.,  0.,  0.,  0.,  1.,  0., -1.,  0.,  0., 0.},
      {1.,  1.,  0.,  1.,  0.,  1.,  0.,  0.,  1., 0., 1.,  0.,  1.,  0.,  0.,  1.,  0.,  0.,  0., 0.},
      {1.,  0.,  1.,  1.,  0.,  0.,  1.,  0.,  0., 1., 1.,  0.,  0.,  0.,  1.,  0.,  1.,  0.,  0., 0.},
      {1., -1.,  0.,  1.,  0., -1.,  0.,  0.,  1., 0., 1.,  0.,  1.,  0.,  0., -1.,  0.,  0.,  0., 0.}};
    double prismA[15][15] = 
    {//0   1   2    3   4   5   6   7   8   9  10  11  12  13  14
      {1., 0., 0., -1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.},
      {1., 1., 0., -1.,-1., 0,  0., 0., 1., 0., 1., 1., 0., 0.,-1.},
      {1., 0., 1., -1., 0.,-1., 0., 0., 0., 1., 1., 0., 1.,-1., 0.},    
      {1., 0., 0.,  1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.},
      {1., 1., 0.,  1., 1., 0,  0., 0., 1., 0., 1., 1., 0., 0., 1.},
      {1., 0., 1.,  1., 0., 1., 0., 0., 0., 1., 1., 0., 1., 1., 0.},
      {1., 0.5, 0.,  -1., -0.5,  0.,   0.,   0., 0.25,   0.,  1.,  0.5,   0.,  0. , -0.25},
      {1., 0.5, 0.5, -1., -0.5, -0.5,0.25,-0.25, 0.25, 0.25,  1.,  0.5,  0.5,-0.25, -0.25},
      {1., 0.,  0.5, -1.,  0.0,-0.5,  0.0,  0.0, 0.,   0.25,  1.,  0.0,  0.5,-0.25,  0.0},
      {1., 0., 0.,  0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
      {1., 1., 0.,  0,  0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.},
      {1., 0., 1.,  0,  0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.},
      {1., 0.5, 0.,   1.,  0.5,  0.,   0.,   0., 0.25,   0.,  1.,  0.5,   0.,  0. ,  0.25},
      {1., 0.5, 0.5,  1.,  0.5,  0.5,0.25, 0.25, 0.25, 0.25,  1.,  0.5,  0.5, 0.25,  0.25},
      {1., 0.,  0.5,  1.,  0.0, 0.5,  0.0,  0.0, 0.,   0.25,  1.,  0.0,  0.5, 0.25,  0.0}}; 
    double pyrA[13][13] ={ 
      //0    1    2     3    4    5    6    7    8    9    10    11    12
      {1., -1., -1.,  -1.,  1.,  1.,  1.,  1.,  1.,  1.,  -1.,  -1.,  -1.},
      {1.,  1., -1.,  -1., -1.,  1.,  1.,  1., -1.,  1.,  -1.,   1.,   1.},
      {1.,  1.,  1.,  -1.,  1.,  1.,  1.,  1., -1., -1.,   1.,   1.,  -1.},
      {1., -1.,  1.,  -1., -1.,  1.,  1.,  1.,  1., -1.,   1.,  -1.,   1.},
      {1.,  0.,  0.,   1.,  0.,  0.,  0.,  1.,  0.,  0.,   0.,   0.,   0.},
      {1.,  0., -1.,  -1.,  0.,  0.,  1.,  1.,  0.,  1.,   0.,   0.,   0.},
      {1.,  1.,  0.,  -1.,  0.,  1.,  0.,  1., -1.,  0.,   0.,   0.,   0.},
      {1.,  0.,  1.,  -1.,  0.,  0.,  1.,  1.,  0., -1.,   0.,   0.,   0.},
      {1., -1.,  0.,  -1.,  0.,  1.,  0.,  1.,  1.,  0.,   0.,   0.,   0.},
      {1.,-0.5,-0.5,   0.,0.25,0.25,0.25,  0.,  0.,  0.,-0.125,-0.125, 0.},
      {1., 0.5,-0.5,   0.,-0.25,0.25,0.25, 0.,  0.,  0.,-0.125, 0.125, 0.}, 
      {1., 0.5, 0.5,   0.,0.25,0.25, 0.25, 0.,  0.,  0., 0.125, 0.125, 0.},
      {1.,-0.5, 0.5,   0.,-0.25,0.25,0.25, 0.,  0.,  0., 0.125,-0.125, 0.}};
  
    int *indx;
    size_t nshl;
    int dim,indx2[3];
    double xl[3],dxl[3];
    
    switch(entityType) {
    case VERTEX:
      break;
    case EDGE:
      nshl = (m_order == 1 ? 2 : 3);
      dim = 1;
      break;

    case TRI:
      nshl = (m_order == 1 ? 3 : 6);
      dim = 2;
      break;

    case QUAD:
      nshl = (m_order == 1 ? 4 : 8);
      dim = 2;
      break;
      
    case TET:
      nshl = (m_order == 1 ? 4 : 10);
      dim = 3;
      break;
      
    case HEX:
      nshl = (m_order == 1 ? 8 : 20);
      dim = 3;
      break;

    case PRISM:
      nshl = (m_order == 1 ? 6 : 15);
      dim = 3;
      break;

    case PYRAMID:
      nshl = (m_order == 1 ? 5 : 13);
      dim = 3;
      break;
      
    } 
    
    indx = new int[nshl];
    
    double *xel, *yel, *zel, xisol[3], fnumber, **A;
    xel = new double [nshl]; yel = new double [nshl]; zel = new double [nshl];
    A = new double*[3];
    for(int i=0; i<3; i++) A[i] = new double [nshl];
    double **M, **MS;
    size_t i, j;
    
    for(i=0; i<knots.size(); i++) {
      mPoint pt = knots[i];
      xel[i] = pt(0);
      yel[i] = pt(1);
      zel[i] = pt(2);
    }
      
    for(i=0; i<nshl;i++){
      A[0][i]=xel[i];
      A[1][i]=yel[i];
      A[2][i]=zel[i];
    }

    MS = new double*[3];
    for(i=0; i<3; i++) MS[i] = new double[3];
    M = new double*[nshl];
    for(i=0; i<nshl; i++) M[i] = new double[nshl];
   
    
    switch(entityType) {
    case VERTEX:
      break;
    case EDGE:
      for(i=0; i<nshl; i++)
	for(j=0; j<nshl; j++)
	  M[i][j] = edgeA[i][j];
      break;
      
    case TRI:
      for(i=0; i<nshl; i++)
	for(j=0; j<nshl; j++)
	  M[i][j] = triA[i][j];
      break;
      
    case QUAD:
      for(i=0; i<nshl; i++)
	for(j=0; j<nshl; j++)
	  M[i][j] = quadA[i][j];
      break;
      
    case TET:
      for(i=0; i<nshl; i++)
	for(j=0; j<nshl; j++)
	  M[i][j] = tetA[i][j];
      break;
      
    case HEX:
      for(i=0; i<nshl; i++)
	for(j=0; j<nshl; j++)
	  M[i][j] = hexA[i][j];
      break;
      
    case PRISM:
      for(i=0; i<nshl; i++)
	for(j=0; j<nshl; j++)
	  M[i][j] = prismA[i][j];
      break;
      
    case PYRAMID:

      for(i=0; i<nshl; i++)
	for(j=0; j<nshl; j++)
	  M[i][j] = pyrA[i][j];
      break;
    }
      

    ludcmp_(M, nshl, indx, &fnumber);
    lubksb_(M, nshl, indx, A[0] );
    lubksb_(M, nshl, indx, A[1] );
    lubksb_(M, nshl, indx, A[2] );
    // Now back substituting to get back the correct set of constants
    // for this element.
      
    for(i =0;i<3;i++){
      for(j=1; j<=3;j++){
	MS[i][j-1]= A[i][j];
      }
    }
    
    
    xisol[0] =  xp - A[0][0];
    xisol[1] =  yp - A[1][0];
    xisol[2] =  zp - A[2][0];
    
    xl[0] = xp;
    xl[1] = yp;
    xl[2] = zp;
    
    // LU decompsing the coeff matrix and solving for xisol
    ludcmp_( MS, dim, indx2, &fnumber);
    lubksb_( MS, dim, indx2, xisol);
    
      
    
    int flag = 1, loop = 0; double tol = 0.000001;
    
    
    while(flag){
      // do a specific number of newton corrections to correct the linear
      // approximation we have in xisol[3];
      
      // update the co-ordinates of the given point
      switch(entityType) {
      case VERTEX:
	break;
	
      case EDGE:
	dxl[0] = xl[0] - (A[0][0] + A[0][1]*xisol[0]);
	MS[0][0] = A[0][1];
	if(m_order == 2) {
	  dxl[0] -= A[0][2]*xisol[0]*xisol[0];
	  MS[0][0] += 2*A[0][2]*xisol[0];
	}
	
	dxl[1] = 0.; dxl[2] = 0.; xisol[1] = 0.; xisol[2] = 0.;
	break;
	
	
      case TRI:
	for(i=0; i<2; i++) {
	  dxl[i] = xl[i] - (A[i][0] + A[i][1]*xisol[0]+A[i][2]*xisol[1]);
	  MS[i][0] = A[i][1]; 
	  MS[i][1] = A[i][2];
	}
	//Quadratic 
	if(m_order == 2) {
	  for(i=0; i<2; i++) {
	    dxl[i] -= (A[i][3]*xisol[0]*xisol[0]+A[i][4]*xisol[1]*xisol[1]+ A[i][5]*xisol[0]*xisol[1]);
	    MS[i][0] += 2.*A[i][3]*xisol[0]+ A[i][5]*xisol[1];
	    MS[i][1] += 2.*A[i][4]*xisol[1]+ A[i][5]*xisol[0]; 
	  }
	}
	dxl[2] = 0.; xisol[2] = 0.;
	break;
	
      case QUAD:
	for(i=0; i<2; i++) {
	  dxl[i] = xl[i] - (A[i][0] + A[i][1]*xisol[0]+A[i][2]*xisol[1]+A[i][3]*xisol[0]*xisol[1]);
	  MS[i][0] = A[i][1] + A[i][3]*xisol[1]; 
	  MS[i][1] = A[i][2] + A[i][3]*xisol[0];
	}
	//Quadratic 
	if(m_order == 2) {
	  for(i=0; i<2; i++) {
	    dxl[i] -= (A[i][4]*xisol[0]*xisol[0]+A[i][4]*xisol[1]*xisol[1]+ 
		       A[i][6]*xisol[0]*xisol[0]*xisol[1]+A[i][7]*xisol[0]*xisol[1]*xisol[1]);
	    
	    MS[i][0] += 2.*A[i][4]*xisol[0]+ 2.*A[i][6]*xisol[0]*xisol[1]+A[i][7]*xisol[1]*xisol[1]; 
	    MS[i][1] += 2.*A[i][5]*xisol[1]+ 2.*A[i][7]*xisol[0]*xisol[1]+A[i][6]*xisol[0]*xisol[0]; 
	  }
	}
	dxl[2] = 0.; xisol[2] = 0.;
	break;
	
      case TET:
	for(i=0; i<3; i++) {
	  dxl[i] = xl[i] - (A[i][0] + A[i][1]*xisol[0] + A[i][2]*xisol[1] + A[i][3]*xisol[2]);
	  MS[i][0] = A[i][1];
	  MS[i][1] = A[i][2];
	  MS[i][2] = A[i][3];
	}
	//Quadratic
	if(m_order == 2) {
	  for(i=0; i<3; i++) {
	    dxl[i] -= (A[i][4]*xisol[0]*xisol[0] + A[i][5]*xisol[1]*xisol[1] + A[i][6]*xisol[2]*xisol[2] +
	               A[i][7]*xisol[0]*xisol[1] + A[i][8]*xisol[0]*xisol[2] + A[i][9]*xisol[1]*xisol[2]);
	    MS[i][0] += 2.*A[i][4]*xisol[0] + A[i][7]*xisol[1] + A[i][8]*xisol[2];
	    MS[i][1] += 2.*A[i][5]*xisol[1] + A[i][7]*xisol[0] + A[i][9]*xisol[2];
	    MS[i][2] += 2.*A[i][6]*xisol[2] + A[i][8]*xisol[0] + A[i][9]*xisol[1];
	  }
	}
	
	break;
	
      case HEX:
	for(i=0; i<3; i++) {
	  dxl[i] = xl[i] - (A[i][0] + A[i][1]*xisol[0] + A[i][2]*xisol[1] + A[i][3]*xisol[2] + 
			    A[i][4]*xisol[0]*xisol[1] + A[i][5]*xisol[0]*xisol[2] + A[i][6]*xisol[1]*xisol[2]
			    + A[i][7]*xisol[0]*xisol[1]*xisol[2]);
	  MS[i][0] = A[i][1] + A[i][4]*xisol[1] + A[i][5]*xisol[2] + A[i][7]*xisol[1]*xisol[2];
	  MS[i][1] = A[i][2] + A[i][4]*xisol[0] + A[i][6]*xisol[2] + A[i][7]*xisol[0]*xisol[2];
	  MS[i][2] = A[i][3] + A[i][5]*xisol[0] + A[i][6]*xisol[1] + A[i][7]*xisol[0]*xisol[1];
	}
	//Quadratic
	if(m_order == 2) {
	  for(i=0; i<3; i++) {
	    dxl[i] -= (A[i][8]*xisol[0]*xisol[0] + A[i][9]*xisol[1]*xisol[1] + A[i][10]*xisol[2]*xisol[2] +
		       A[i][11]*xisol[0]*xisol[0]*xisol[1] + A[i][12]*xisol[0]*xisol[0]*xisol[2] +
		       A[i][13]*xisol[1]*xisol[1]*xisol[0] + A[i][14]*xisol[1]*xisol[1]*xisol[2] +
		       A[i][15]*xisol[2]*xisol[2]*xisol[0] + A[i][16]*xisol[2]*xisol[2]*xisol[1] +
		       A[i][17]*xisol[0]*xisol[0]*xisol[1]*xisol[2] + 
		       A[i][18]*xisol[1]*xisol[1]*xisol[0]*xisol[2] +
		       A[i][19]*xisol[2]*xisol[2]*xisol[0]*xisol[1]);
	    MS[i][0] += (2.*A[i][8]*xisol[0] + 2.*A[i][11]*xisol[0]*xisol[1] + 2.*A[i][12]*xisol[0]*xisol[2] +
			 A[i][13]*xisol[1]*xisol[1] + A[i][15]*xisol[2]*xisol[2] + 
			 2.*A[i][17]*xisol[0]*xisol[1]*xisol[2] +
			 A[i][18]*xisol[1]*xisol[1]*xisol[2] + A[i][19]*xisol[2]*xisol[2]*xisol[1]);
	    MS[i][1] += (2.*A[i][9]*xisol[1] + A[i][11]*xisol[0]*xisol[0] + 2.*A[i][13]*xisol[0]*xisol[1] +
			 2.*A[i][14]*xisol[1]*xisol[2] + A[i][16]*xisol[2]*xisol[2] + 
			 A[i][17]*xisol[0]*xisol[0]*xisol[2] +
			 2.*A[i][18]*xisol[0]*xisol[1]*xisol[2] + A[i][19]*xisol[2]*xisol[2]*xisol[0]);
	    MS[i][2] += (2.*A[i][10]*xisol[2] + A[i][12]*xisol[0]*xisol[0] + A[i][14]*xisol[2]*xisol[2] +
			 2.*A[i][15]*xisol[2]*xisol[0] + 2.*A[i][16]*xisol[2]*xisol[1] + 
			 A[i][17]*xisol[0]*xisol[0]*xisol[1] +
			 A[i][18]*xisol[1]*xisol[1]*xisol[0] + 2.*A[i][19]*xisol[0]*xisol[1]*xisol[2]);
	  }
	}
	break;
	
      case PRISM:
	for(i=0; i<3; i++) {
	  dxl[i] = xl[i] - (A[i][0] + A[i][1]*xisol[0] + A[i][2]*xisol[1] + 
			    A[i][3]*xisol[2] +  A[i][4]*xisol[0]*xisol[2] + A[i][5]*xisol[1]*xisol[2]);
	  MS[i][0] = A[i][1] + A[i][4]*xisol[2];
	  MS[i][1] = A[i][2] + A[i][5]*xisol[2];
	  MS[i][2] = A[i][3] + A[i][4]*xisol[0] + A[i][5]*xisol[1];
	}
	//Quadratic
	if(m_order == 2) {
	  for(i=0; i<3; i++) {
	    dxl[i] -=  (A[i][6]*xisol[0]*xisol[1] + A[i][7]*xisol[0]*xisol[1]*xisol[2] +
			A[i][8]*xisol[0]*xisol[0] + A[i][9]*xisol[1]*xisol[1] + A[i][10]*xisol[2]*xisol[2] + 
			A[i][11]*xisol[0]*xisol[2]*xisol[2] + A[i][12]*xisol[1]*xisol[2]*xisol[2] + 
			A[i][13]*xisol[1]*xisol[1]*xisol[2] + A[i][14]*xisol[0]*xisol[0]*xisol[2] );
	
	    MS[i][0] += (A[i][6]*xisol[1] + A[i][7]*xisol[1]*xisol[2] + 2.*A[i][8]*xisol[0] + 
			 A[i][11]*xisol[2]*xisol[2] + 2.*A[i][14]*xisol[0]*xisol[2] );
	    MS[i][1] += (A[i][6]*xisol[0] + A[i][7]*xisol[1]*xisol[2] + 2.*A[i][9]*xisol[1] + 
			 A[i][12]*xisol[2]*xisol[2] + 2.*A[i][13]*xisol[1]*xisol[2] );
	    MS[i][2] += (A[i][7]*xisol[1]*xisol[2] + 2.*A[i][10]*xisol[2] + 2.*A[i][11]*xisol[0]*xisol[2] + 
			 2.*A[0][12]*xisol[1]*xisol[2] + A[i][13]*xisol[1]*xisol[1] + A[i][14]*xisol[0]*xisol[0]);
	  }
	}
	break;
	
      case PYRAMID:
	for(i=0; i<3; i++) {
	  dxl[i] = xl[i] - (A[i][0] + A[i][1]*xisol[0] + A[i][2]*xisol[1] + 
			    A[i][3]*xisol[2] +  A[i][4]*xisol[0]*xisol[1] );
	  MS[i][0] = A[i][1] + A[i][4]*xisol[1];
	  MS[i][1] = A[i][2] + A[i][4]*xisol[0];
	  MS[i][2] = A[i][3];
	}
	//Quadratic

	if(m_order == 2) {
	  for(i=0; i<3; i++) {
	    dxl[i] -= (A[i][5]*xisol[0]*xisol[0] + A[i][6]*xisol[1]*xisol[1] + A[i][7]*xisol[2]*xisol[2]
                       + A[i][8]*xisol[0]*xisol[2] + A[i][9]*xisol[1]*xisol[2] 
		       + A[i][10]*xisol[0]*xisol[0]*xisol[1] + A[i][11]*xisol[1]*xisol[1]*xisol[0]
		       + A[i][12]*xisol[0]*xisol[1]*xisol[2]);
	    MS[i][0] += (2.*A[i][5]*xisol[0] + A[i][8]*xisol[2] + 2.*A[i][10]*xisol[0]*xisol[1] 
			 + A[i][11]*xisol[1]*xisol[1] + A[i][12]*xisol[1]*xisol[2]);
	    MS[i][1] += (2.*A[i][6]*xisol[1] + A[i][9]*xisol[2] + A[i][10]*xisol[0]*xisol[0] 
			 + 2.*A[i][11]*xisol[0]*xisol[1] + A[i][12]*xisol[0]*xisol[2]);
	    MS[i][2] += (2.*A[i][7]*xisol[2] + A[i][8]*xisol[0] + A[i][9]*xisol[1] 
			 +  A[i][12]*xisol[0]*xisol[1]);
	  }
	}    
	  
	  
	break;

      }
	
	
      
      ludcmp_( MS, dim, indx2, &fnumber);
      lubksb_( MS, dim, indx2, dxl);
      
      
      for(i=0;i<3;i++) xisol[i]=xisol[i]+dxl[i];
      if( (fabs(dxl[0]) < tol) && (fabs(dxl[1]) < tol) && (fabs(dxl[2]) < tol) )
	flag = 0;
      loop ++;
      
      //maximum recursive correction is 10 steps for linear and 50 for quadratic
      if(loop == 10 && m_order == 1)
	flag = 0;
      else if(m_order == 2 && loop == 50)
	flag = 0;
    }
    
    int truth =1;                        
    if(loop == 10 && m_order == 1)
      truth = 0;
    else if(loop == 50 && m_order == 2)
      truth = 0;
    
    if (truth){
      Upos = xisol[0];                 
      Vpos = xisol[1];                 
      Wpos = xisol[2];  
    }else {
      std::cout<<"Can not get the inverse mapping params!...."<<std::endl;
    }
    
    // Fixed incorrect new/delete pairing --CWS 12/16/07
    for(i=0; i<3; i++) {delete [] MS[i]; delete [] A[i];}
    for(i=0; i<nshl; i++) delete [] M[i]; 
    delete [] MS;
    delete [] M;
    delete [] A;
    delete [] xel;
    delete [] yel;
    delete [] zel;
    delete [] indx;

    return truth;
    

  }

  // certainly not the fastest way !!
  double Mapping::detJac(double u, double v, double w) const 
  {
    mTensor2 t;
    return jacInverse(u,v,w,t);
  }

// P. Hu added for Level Set Integration
  /*
  double Mapping::detJac_way2(double u, double v, double w) const
  {
    mTensor2 t;
    return jacInverse_way2(u,v,w,t);
  }
  */

  /**
     This function is valid for any mesh mapping
  */
  double Mapping::jacInverse(double u, double v, double w, mTensor2 &Invjac) const 
  {

    double jac[3][3];
    double DetJac;

    deval (u,v,w,jac[0][0],jac[0][1],jac[0][2],
	   jac[1][0],jac[1][1],jac[1][2],
	   jac[2][0],jac[2][1],jac[2][2]);

    switch(entityDim)
      {
	// may take into account non x-y elements and curved elements, curved lines,...
      case 1:
	{
	
	  DetJac = sqrt ( jac[0][0]*jac[0][0] + jac[0][1]*jac[0][1]
			  + jac[0][2]*jac[0][2] );
	
	  // vrai is true in french
	  int vrai=0;
	  for(int i=0;i<3;i++){
	    if(jac[0][i]==0){   // if this component of the normal vector is zero
	      vrai =1;
	      for(int j=0;j<3;j++){
		if (j==i)
		  jac[1][j] = 1; // then the component of the second normal must be one,
		else
		  jac[1][j] = 0; // and the other components of the second normal are zero.
	      }
	      continue;
	    }
	  }                           
	  // surface equation with normal vector n : n1*X + n2*Y + n3*Z = 0
	  if(!vrai){ // looking for the second normal vector in the plane z = 0
	    double temp = sqrt ( jac[0][0]*jac[0][0] + jac[0][1]*jac[0][1] );
	    jac[1][0] = -jac[0][1]/temp;
	    jac[1][1] =  jac[0][0]/temp;
	    jac[1][2] =  0;
	  }
	  // The third normal vector
	  jac[2][0] = ( jac[0][1]*jac[1][2] - jac[0][2]*jac[1][1] )/ DetJac;
	  jac[2][1] = ( jac[0][2]*jac[1][0] - jac[0][0]*jac[1][2] )/ DetJac;
	  jac[2][2] = ( jac[0][0]*jac[1][1] - jac[0][1]*jac[1][0] )/ DetJac; 
	}
	break;
      case 2:
	{

	  double d3 = jac[0][0]*jac[1][1] - jac[0][1] * jac[1][0];  
	  double d2 = jac[0][2]*jac[1][0] - jac[0][0] * jac[1][2];  
	  double d1 = jac[0][1]*jac[1][2] - jac[0][2] * jac[1][1];  
	
	  DetJac = sqrt ( d1*d1 + d2*d2 + d3*d3 );
	
	  jac[2][0] = d1/DetJac; 
	  jac[2][1] = d2/DetJac; 
	  jac[2][2] = d3/DetJac;
	}
	break;
      case 3:
	{
	  DetJac = jac[0][0]*jac[1][1]*jac[2][2] + jac[0][2] *jac[1][0]*jac[2][1] +
	    jac[0][1]*jac[1][2]*jac[2][0] - jac[0][2] *jac[1][1]*jac[2][0] -
	    jac[0][0]*jac[1][2]*jac[2][1] - jac[0][1] *jac[1][0]*jac[2][2];
	}
	break;
      }

    Invjac(0,0) = (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1]) / DetJac;
    Invjac(1,0) = -(jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0]) / DetJac;
    Invjac(2,0) = (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]) / DetJac;
  
    Invjac(0,1) = -(jac[0][1] * jac[2][2] - jac[0][2] * jac[2][1]) / DetJac;
    Invjac(1,1) = (jac[0][0] * jac[2][2] - jac[0][2] * jac[2][0]) / DetJac;
    Invjac(2,1) = -(jac[0][0] * jac[2][1] - jac[0][1] * jac[2][0]) / DetJac;
  
    Invjac(0,2) = (jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1]) / DetJac;
    Invjac(1,2) = -(jac[0][0] * jac[1][2] - jac[0][2] * jac[1][0]) / DetJac;
    Invjac(2,2) = (jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]) / DetJac;

    return DetJac;
  }

// P. Hu added for Level Set Integration
  /*
  double Mapping::jacInverse_way2(double u, double v, double w, mTensor2 &Invjac) const
  {
          
    double jac[3][3];
    int i;
    double DetJac;
                
    deval (u,v,w,jac[0][0],jac[0][1],jac[0][2],
           jac[1][0],jac[1][1],jac[1][2],
           jac[2][0],jac[2][1],jac[2][2]);
               
    switch(ent->getElementDim())
      {
        // may take into account non x-y elements and curved elements, curved lines,...
      case 1:
        {
            
          DetJac = sqrt ( jac[0][0]*jac[0][0] + jac[0][1]*jac[0][1]
                          + jac[0][2]*jac[0][2] );
            
          // vrai is true in french
          int vrai=0;
          for(int i=0;i<3;i++){
            if(jac[0][i]==0){   // if this component of the normal vector is zero
              vrai =1;
              for(int j=0;j<3;j++){
                if (j==i)
                  jac[1][j] = 1; // then the component of the second normal must be one,
                else
                  jac[1][j] = 0; // and the other components of the second normal are zero.
              }
              continue;
            }
          }
          // surface equation with normal vector n : n1*X + n2*Y + n3*Z = 0
          if(!vrai){ // looking for the second normal vector in the plane z = 0
            double temp = sqrt ( jac[0][0]*jac[0][0] + jac[0][1]*jac[0][1] );
            jac[1][0] = -jac[0][1]/temp;
            jac[1][1] =  jac[0][0]/temp;
            jac[1][2] =  0;
          }   
          // The third normal vector
          jac[2][0] = ( jac[0][1]*jac[1][2] - jac[0][2]*jac[1][1] )/ DetJac;
          jac[2][1] = ( jac[0][2]*jac[1][0] - jac[0][0]*jac[1][2] )/ DetJac;
          jac[2][2] = ( jac[0][0]*jac[1][1] - jac[0][1]*jac[1][0] )/ DetJac;
        }
        break;
      case 2: 
        { 
    
          double d3 = jac[0][0]*jac[1][1] - jac[0][1] * jac[1][0];
          double d2 = jac[0][2]*jac[1][0] - jac[0][0] * jac[1][2];
          double d1 = jac[0][1]*jac[1][2] - jac[0][2] * jac[1][1];
    
          DetJac = sqrt ( d1*d1 + d2*d2 + d3*d3 );
          
          jac[2][0] = d1/DetJac;
          jac[2][1] = d2/DetJac;
          jac[2][2] = d3/DetJac;
        }
        break;
      case 3:
        {
          DetJac = jac[0][0]*jac[1][1]*jac[2][2] + jac[0][2] *jac[1][0]*jac[2][1] +
            jac[0][1]*jac[1][2]*jac[2][0] - jac[0][2] *jac[1][1]*jac[2][0] -
            jac[0][0]*jac[1][2]*jac[2][1] - jac[0][1] *jac[1][0]*jac[2][2];
        }
        break;
      }
          
    Invjac(0,0) = (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1]) / DetJac;
    Invjac(1,0) = -(jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0]) / DetJac;
    Invjac(2,0) = (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]) / DetJac;
              
    Invjac(0,1) = -(jac[0][1] * jac[2][2] - jac[0][2] * jac[2][1]) / DetJac;
    Invjac(1,1) = (jac[0][0] * jac[2][2] - jac[0][2] * jac[2][0]) / DetJac;
    Invjac(2,1) = -(jac[0][0] * jac[2][1] - jac[0][1] * jac[2][0]) / DetJac;
                
    Invjac(0,2) = (jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1]) / DetJac;
    Invjac(1,2) = -(jac[0][0] * jac[1][2] - jac[0][2] * jac[1][0]) / DetJac;
    Invjac(2,2) = (jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]) / DetJac;
             
    return DetJac;
  }
  */

  bool Mapping::inReferenceElement(double u, double v, double w) const
  {
    const double eps = 1.e-6;

    switch(entityType)
      {
      case TRI :
	if(u < -eps || v < -eps || 1.-u-v < -eps)return false;
	break;
      case QUAD :
	if(u < (-1.-eps) || u > (1+eps) || v < (-1.0-eps) || v > (1.0+eps))return false;
	break;
      case EDGE :
	if(u < -1.-eps || u > 1+eps)return false;
	break;
      case TET :
	if(u < -eps || v < -eps || w < -eps || 1.-u-v-w < -eps)return false;
	break;
      case HEX :
	if(u < (-1.-eps) || u > (1+eps) || v < (-1.-eps) || v > (1+eps) || w < (-1.-eps) || w > (1+eps))return false;
	break;
      case PRISM :
	if(u < -eps || v < -eps || 1.-u-v < -eps || w < (-1.-eps) || w > (1+eps))return false;
	break;
      default: break; // need to throw an exception here
      }
    return true;
  }

  bool Mapping :: interiorCheck (pEntity theMeshEntity, const mPoint &p , 
				 double &u, double &v, double &w) const
  {
    mPoint pMin,pMax;
    double eps = 1.e-5;
  
    boundingBox(pMin,pMax);

    if(p(0) > pMax(0)+eps ||
       p(1) > pMax(1)+eps ||
       p(2) > pMax(2)+eps ||
       p(0) < pMin(0)-eps ||
       p(1) < pMin(1)-eps ||
       p(2) < pMin(2)-eps)return false;
    if(!invert(p(0),p(1),p(2),u,v,w)) return false;
    if(!inReferenceElement(u,v,w))return false;
    return true;
  }


  void Mapping::COG(double &u, double &v, double &w) const
  {
    switch(entityType)
      {
      case PRISM :
      case TRI :
	u = v = 0.3333333333333;
	w = 0;
	break;
      case HEX :
      case QUAD :
      case PYRAMID :
	u = v = w = 0.;
	break;		
      case EDGE :
	u = v = w = 0.0;
	break;
      case TET :
	u = v = w = 0.25;
	break;
      default: break; // need to throw an exception here
      }
  }

  int Mapping :: simplexation (  int &dim, int &nbPnt, double *u, double *v, double *w, int combn[][4] )const
  {

    const double pts [15][3] = { 
      { -1,-1,-1 },{ 1,-1,-1 },{ 1,1,-1 },{-1,1,-1}, // vertices
      { -1,-1,1 },{ 1,-1,1 },{ 1,1,1 },{-1,1,1} ,
      {0,0,0}, // center
      {-1,0,0},{1,0,0}, //faces centers
      {0,-1,0},{0,1,0},
      {0,0,-1},{0,0,1}
    };
    int tets [24][4] = { {8,4,11,5},{8,5,11,1},{8,1,11,0},{8,0,11,4}, // face y = -1
			 {8,7,12,6},{8,6,12,2},{8,2,12,3},{8,3,12,7}, // face y = 1
			 {8,5,10,6},{8,6,10,2},{8,2,10,1},{8,1,10,5},     // face x = -1
			 {8,4,9,7},{8,7,9,3},{8,3,9,0},{8,0,9,4},     // face x =  1
			 {8,0,13,1},{8,1,13,2},{8,2,13,3},{8,3,13,0}, // face z = -1
			 {8,4,14,5},{8,5,14,6},{8,6,14,7},{8,7,14,4} }; // face z =  1

    switch(entityType)
      {
      case TET :
	u[0] = 0.0; v[0] = 0.0; w[0] = 0.0;
	u[1] = 1.0; v[1] = 0.0; w[1] = 0.0;
	u[2] = 0.0; v[2] = 1.0; w[2] = 0.0;
	u[3] = 0.0; v[3] = 0.0; w[3] = 1.0;
	dim = 4;
	combn[0][0] = 0;combn[0][1] = 1;combn[0][2] = 2;combn[0][3] = 3;
	nbPnt = 4;
	return 1;	
	break;
      case HEX :
        {
	for (int i=0;i<15;i++)
	  {
	    u[i] = pts[i][0];
	    v[i] = pts[i][1];
	    w[i] = pts[i][2];
	  }
	}
	nbPnt = 15;
	{
	for (int i=0;i<24;i++)
	  {
	    combn[i][0] = tets[i][0];
	    combn[i][1] = tets[i][1];
	    combn[i][2] = tets[i][2];
	    combn[i][3] = tets[i][3];
	  }
	}
	dim = 4;
	return 24;
	break;
      case TRI :
	u[0] = 0.0; v[0] = 0.0; w[0] = 0.0;
	u[1] = 1.0; v[1] = 0.0; w[1] = 0.0;
	u[2] = 0.0; v[2] = 1.0; w[2] = 0.0;
	nbPnt = 3;
	combn[0][0] = 0;combn[0][1] = 1;combn[0][2] = 2;
	dim = 3;
	return 1;
	break;
      case QUAD :
	u[0] = 0.0;  v[0] = 0.0;  w[0] = 0.0;
	u[1] =-1.0;  v[1] = 1.0;  w[1] = 0.0;
	u[2] = 1.0;  v[2] = 1.0;  w[2] = 0.0;
	u[3] = 1.0;  v[3] = -1.0; w[3] = 0.0;
	u[4] =-1.0;  v[4] = -1.0; w[3] = 0.0;
	nbPnt = 5;
	combn[0][0] = 0;combn[0][1] = 1;combn[0][2] = 2;
	combn[1][0] = 0;combn[1][1] = 2;combn[1][2] = 3;
	combn[2][0] = 0;combn[2][1] = 3;combn[2][2] = 4;
	combn[3][0] = 0;combn[3][1] = 4;combn[3][2] = 1;
	dim = 3;
	return 4;
	break;
      default:
	throw 1;
      }
  }

  void Mapping::normalVector(double u , double v , double w, mVector &n) const
  {
    double jac[3][3];
    deval (u,v,w,
	   jac[0][0],jac[0][1],jac[0][2],
	   jac[1][0],jac[1][1],jac[1][2],
	   jac[2][0],jac[2][1],jac[2][2]);

    switch(entityDim)
      {
      case 0:
	throw 1;
      case 1:
	{
	  mVector t (jac[0][0],jac[0][1],jac[0][2]);
	  mVector z (0,0,1);
	  n = t % z;
	  n.norm();
	}
	break;
      case 2:
	{
	  mVector t1 (jac[0][0],jac[0][1],jac[0][2]);
	  mVector t2 (jac[1][0],jac[1][1],jac[1][2]);
	  n = t1 % t2;
	  n.norm();
	}
	break;
      default:
	throw 1;
      }	
  }
