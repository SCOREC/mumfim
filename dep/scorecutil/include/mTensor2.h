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
// mTensor2.h: interface for the mTensor2 class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _SCOREC_UTIL_MTensor2_H_
#define _SCOREC_UTIL_MTensor2_H_

#include <stdio.h>
#include <cmath>
#include <algorithm>
#include "mVector.h"
#include "scorec_function_objects.h"

  /**
     Simple class for a 2nd order tensor in 3d which
     has a 3x3 matrix representation.
  */

  class mTensor2 {
  public:
    static const int NSD = 3;
    typedef double value_type;
    
    explicit mTensor2(const mVector &c1, 
		      const mVector &c2, 
		      const mVector &c3);
    explicit mTensor2(double init = 0.0);
    virtual ~mTensor2();


// These operators should be free functions!
    inline const mVector operator*(const mVector &other) const {
      mVector m(0,0,0);
      for (int i=0; i < NSD; ++i) {
	m[i] += pos[i] * other;
      }
      return m;
    }
    inline const mTensor2 operator*(const mTensor2 &other) const
    {
      mTensor2 m(0);
      for(int i=0;i<NSD;i++) {
	for(int j=0;j<NSD;j++) {
	  for(int l=0;l<NSD;l++) {
	    m(i,j) += (*this)(i,l) * other(l,j);
	  }
	}
      }
      return m;
    }
    inline mTensor2& operator*=(const double &scalar)
    {
      for(int i=0;i<NSD;i++) {
	for(int j=0;j<NSD;j++) {
	  (*this)(i,j) *= scalar;
	}
      }
      return *this;
    }

    inline mTensor2& operator/=(const double &scalar)
    {
      for(int i=0;i<NSD;i++)
	{
	  for(int j=0;j<NSD;j++)
	    {
	      (*this)(i,j) /= scalar;
	    }
	}
      return *this;
    }

    inline mTensor2& operator += (const mTensor2 &other)
    {
      for(int i=0;i<NSD;i++)
	{
	  for(int j=0;j<NSD;j++)
	    {
	      (*this)(i,j) += other.pos[i][j];
	    }
	}
      return *this;
    }

    inline mTensor2 &operator -= (const mTensor2 &other)
    {
      for(int i=0;i<NSD;i++) {
	for(int j=0;j<NSD;j++) {
	  (*this)(i,j) -= other.pos[i][j];
	}
      }
      return *this;
    }

    inline const mTensor2 operator*(const double &scalar) const
    {
      mTensor2 other;
      for(int i=0;i<NSD;i++)
	{
	  for(int j=0;j<NSD;j++)
	    {
	      other.pos[i][j] = scalar * (*this)(i,j);
	    }
	}
      return other;
    }
    double& operator()( int, int );
    double  operator()( int, int ) const;

    const mVector& operator[]( int i ) const;       
    mVector& operator[]( int i );                  
    const mVector col( int j ) const;                   

    const mTensor2 operator+(const mTensor2 &other) const;
    /// solve the system
    const mVector operator/(const mVector & other) const;
    const mTensor2 invert() const;
    /// determinant
    long eigen   (mVector e[NSD],double v[NSD]) const;
    long eigen2d (mVector e[NSD],double v[NSD]) const;
    inline void symmetrize();
    const mTensor2 intersect2d (mTensor2 &t) const;
    inline void set_all(const double &scalar);

    void transpose();
    const mTensor2 operator!() const; 
  private:
    mVector pos[NSD];
  };

  inline const mTensor2 transpose(const mTensor2& t) { return !t; }
  double det(const mTensor2& t);
  double trace(const mTensor2& t);
  double trace2(const mTensor2& t);


  inline void mTensor2::symmetrize ()
  {
    pos[0][1] = pos[1][0] = 0.5 * (pos[0][1] + pos[1][0]);
    pos[0][2] = pos[2][0] = 0.5 * (pos[0][2] + pos[2][0]);
    pos[2][1] = pos[1][2] = 0.5 * (pos[2][1] + pos[1][2]);
  }

  inline mTensor2::mTensor2(const mVector &c0, const mVector &c1, const mVector &c2)
  {
    pos[0] = c0; pos[1] = c1; pos[2] = c2;
  }

  inline double& mTensor2::operator()(int i,int j) {
    return (*this)[i][j];
  }

  inline double mTensor2::operator()(int i,int j) const
  {
    return (*this)[i][j];
  }

  inline const mVector& mTensor2::operator[]( int i ) const
  {
  //    return the ith row
                return pos[i];
  }

  inline mVector& mTensor2::operator[]( int i )
  {
  //    assign the ith row
                return pos[i];
  }

  inline const mVector mTensor2::col( int j ) const
  {
          return mVector(pos[0][j], pos[1][j], pos[2][j]);
  }

  inline const mTensor2 mTensor2::operator+(const mTensor2 &other) const
  {
    mTensor2 m(0);
    for (int i = 0; i < NSD; ++i)
      m.pos[i] += pos[i] + other[i];
    
    return m;
  }
  inline bool operator==(const mTensor2& a, const mTensor2& b) {
    for(int i = 0; i < mTensor2::NSD; ++i) {
      for(int j = 0; j < mTensor2::NSD; ++j) {
	if(a(i,j) != b(i,j)) return false;
      }
    }
    return true;
  }

  long FindCubicRoots(const double coeff[4], double x[3]);
  long NullSpace(const double *a, double *result, double eps, long n);


  inline const mVector operator*(const mVector &other,const mTensor2 &t)
  {
    mVector m(0,0,0);
    for(int i=0;i < mVector::NSD;i++) {
      for(int j=0;j < mVector::NSD;j++) {
	m[i] += t(i,j) * other[j];
      }
    }
    return m;
  }


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

  inline mTensor2::mTensor2(double init)
  {
    for(int i=0;i<NSD;i++) {
      for(int j=0;j<NSD;j++) {
	operator()(i,j) = init;
      }
    }
  }
  
  inline mTensor2::~mTensor2()
  {
  }
  
  inline double trace(const mTensor2& t)
  {
    return t[0][0] + t[1][1] + t[2][2];
  }

  inline double trace2(const mTensor2& t)
  {
    double a00 = (t[0][0] * t[0][0] + 
		  t[1][0] * t[0][1] + 
		  t[2][0] * t[0][2]); 
    double a11 = (t[1][0] * t[0][1] + 
		  t[1][1] * t[1][1] + 
		  t[1][2] * t[2][1]); 
    double a22 = (t[2][0] * t[0][2] + 
		  t[2][1] * t[1][2] + 
		  t[2][2] * t[2][2]);

    return a00 + a11 + a22;
  }

  inline double det(const mTensor2& t) {
    return (t[0][0] * (t[1][1] * t[2][2] - t[1][2] * t[2][1]) -
	    t[0][1] * (t[1][0] * t[2][2] - t[1][2] * t[2][0]) +
	    t[0][2] * (t[1][0] * t[2][1] - t[1][1] * t[2][0]) );
  }
  
  inline const mVector mTensor2::operator/(const mVector & b) const
  {
    mVector res;
    
    double detm = det(*this);
    
    if (detm == 0.0) {
      throw; // Singular!
    }
    
    double ud = 1. / (detm);
    
    res[0] = b[0] * (pos[1][1] * pos[2][2] - pos[1][2] * pos[2][1]) -
      pos[0][1] * (b[1] * pos[2][2] - pos[1][2] * b[2]) +
      pos[0][2] * (b[1] * pos[2][1] - pos[1][1] * b[2]);
    
    res[1] = pos[0][0] * (b[1] * pos[2][2] - pos[1][2] * b[2]) -
      b[0] * (pos[1][0] * pos[2][2] - pos[1][2] * pos[2][0]) +
      pos[0][2] * (pos[1][0] * b[2] - b[1] * pos[2][0]);
    
    res[2] = pos[0][0] * (pos[1][1] * b[2] - b[1] * pos[2][1]) -
      pos[0][1] * (pos[1][0] * b[2] - b[1] * pos[2][0]) +
      b[0] * (pos[1][0] * pos[2][1] - pos[1][1] * pos[2][0]);
    
    for (int i = 0; i < NSD; i++)
      res[i] *= ud;
    
    return res;
  }

  inline const mTensor2 mTensor2::invert() const
  {
    mTensor2 inv;
    
    double detm = det(*this);
    
    if (detm == 0.0) {
      throw; // Singular!
    }
    
    double ud = 1. / (detm);
    inv.pos[0][0] = ud * (pos[1][1] * pos[2][2] - pos[1][2] * pos[2][1]);
    inv.pos[0][1] = -ud * (pos[1][0] * pos[2][2] - pos[1][2] * pos[2][0]);
    inv.pos[0][2] = ud * (pos[1][0] * pos[2][1] - pos[1][1] * pos[2][0]);
    inv.pos[1][0] = -ud * (pos[0][1] * pos[2][2] - pos[0][2] * pos[2][1]);
    inv.pos[1][1] = ud * (pos[0][0] * pos[2][2] - pos[0][2] * pos[2][0]);
    inv.pos[1][2] = -ud * (pos[0][0] * pos[2][1] - pos[0][1] * pos[2][0]);
    inv.pos[2][0] = ud * (pos[0][1] * pos[1][2] - pos[0][2] * pos[1][1]);
    inv.pos[2][1] = -ud * (pos[0][0] * pos[1][2] - pos[0][2] * pos[1][0]);
    inv.pos[2][2] = ud * (pos[0][0] * pos[1][1] - pos[0][1] * pos[1][0]);
    return inv;
  }
  
  inline long FindCubicRoots(const double coeff[4], double x[3]) {
    using std::abs;
    double a1 = coeff[2] / coeff[3];
    double a2 = coeff[1] / coeff[3];
    double a3 = coeff[0] / coeff[3];
    
    double Q = (a1 * a1 - 3 * a2) / 9.;
    double R = (2. * a1 * a1 * a1 - 9. * a1 * a2 + 27. * a3) / 54.;
    double Qcubed = Q * Q * Q;
    double d = Qcubed - R * R;

    //    printf ("d = %22.15e Q = %12.5E R = %12.5E Qcubed %12.5E\n",d,Q,R,Qcubed);

    /// three roots, 2 equal 
    if (Qcubed == 0.0 || abs(Qcubed - R * R) < 	1.e-8 * (abs(Qcubed) + abs(R * R)) ) {
	double theta;
	if (Qcubed <= 0.0)theta = acos(1.0);
	else if (R / sqrt(Qcubed) > 1.0)theta = acos(1.0); 
	else if (R / sqrt(Qcubed) < -1.0)theta = acos(-1.0); 
	else theta = acos(R / sqrt(Qcubed));
	double sqrtQ = sqrt(Q);
	//	printf("sqrtQ = %12.5E teta=%12.5E a1=%12.5E\n",sqrt(Q),theta,a1);
	x[0] = -2 * sqrtQ * cos( theta           / 3) - a1 / 3;
	x[1] = -2 * sqrtQ * cos((theta + 2 * M_PI) / 3) - a1 / 3;
	x[2] = -2 * sqrtQ * cos((theta + 4 * M_PI) / 3) - a1 / 3;
      return 3;
      }

    /* Three real roots */
    if (d >= 0.0) {
      double theta = acos(R / sqrt(Qcubed));
      double sqrtQ = sqrt(Q);
      x[0] = -2 * sqrtQ * cos( theta           / 3) - a1 / 3;
      x[1] = -2 * sqrtQ * cos((theta + 2 * M_PI) / 3) - a1 / 3;
      x[2] = -2 * sqrtQ * cos((theta + 4 * M_PI) / 3) - a1 / 3;
      return 3;
    }
    
    /* One real root */
    else {
      printf("IMPOSSIBLE !!!\n");
      // If it's impossible, why not throw an exception? --BFD 10/18/2006

      double e = pow(sqrt(-d) + abs(R), 1. / 3.);
      if (R > 0)
	e = -e;
      x[0] = (e + Q / e) - a1 / 3.;
      return 1;
    }
  }

  
  inline long NullSpace(const double * const a, double * const result, const double eps, const long n) {
    const array_ref_2d<double> R(result, n);
    using std::abs;
    const int MAXN = 32;
    int r[MAXN], c[MAXN];
    register long i, j, k;
    int jj, kk, t;
    double max, temp;
    int ec;
    
    for (i = 0; i < n; i++)
      r[i] = c[i] = -1;			/* Reset row and column pivot indices */
    
    // copy the input matrix if not in place
    if (result != a) 
      for (i = 0; i < n*n; i++)  
	result[i] = a[i];
    // rest of algorithm is in place wrt result[]
    
    for (i = 0; i < n; i++) {
      /* Find the biggest element in the remaining submatrix
       * for the next full pivot.
       */
      max = 0.0;
      for (k = 0; k < n; k++) {
	if (r[k] < 0) {
	  for (j = 0; j < n; j++) {
	    if ((c[j] < 0) && ((temp = abs(R(k, j))) > max)) {
	      kk = k;
	      jj = j;
	      max = temp;
	    }
	  }
	}
      }
      if (max < eps)
	break;		/* Consider this and all subsequent pivots to be zero */
      
      c[jj] = kk;					/* The row */
      r[kk] = jj;					/*	      and column of the next pivot */
      
      temp = 1.0 / R(kk, jj);
      R(kk, jj) = 1.0;
      for (j = 0; j < n; j++)		/* Should this be for j != jj ? */
	R(kk, j) *= temp;		/* Row equilibration */
      
      for (k = 0; k < n; k++) {	/* Row elimination */
	if (k == kk)
	  continue;			/* Don't do a thing to the pivot row */
	temp = R(k, jj);
	R(k, jj) = 0.0;
	for (j = 0; j < n; j++) {
	  R(k, j) -= temp * R(kk, j);	/* Subtract row kk from row k */
	  if (abs(R(k, j)) < eps)
	    R(k, j) = 0.0;	/* Flush to zero if too small */
	}
      }
    }
    
    /* Sort into a truncated triangular matrix */
    for (j = 0; j < n; j++) {		/* For all columns... */
      while ((c[j] >= 0) && (j != c[j])) {
	for (k = 0; k < n; k++) {
	  if (r[k] < 0) {
	    /* Aha! a null column vector */
	    temp = R(k, j);	/* Get it on top */
	    R(k, j) = R(k, c[j]);
	    R(k, c[j]) = temp;
	  }
	}
	t = c[j];				/* Twiddle until pivots are on the diagonal */
	c[j] = c[t];
	c[t] = t;
      }
    }
    
    /* Copy the null space vectors into the top of the A matrix */
    ec = 0;
    for (k = 0; k < n; k++) {
      if (r[k] < 0) {
	R(k, k) = 1.0;			/* Set the pivot equal to 1 */
	if (ec != k) {
	  for (j = 0; j < n; j++) {
	    R(ec, j) = R(k, j);
	  }
	}
	ec++;
      }
    }
    /* The first  ec  rows of the matrix  a  are the vectors which are
     * orthogonal to the columns of the matrix  a.
     */
    return ec;
  }    

  inline long mTensor2::eigen(mVector e[NSD], double v[NSD]) const
  {            
    /// characteristic polynomial of T : find v root of
    /// v^3 - I1 v^2 + I2 T + I3 = 0
    /// I1 : first invariant , trace(T)
    /// I2 : second invariant , 1/2 (I1^2 -trace(T^2))
    /// I3 : third invariant , det T
    double I[4];
    I[3] = 1.0;
    I[2] = -trace(*this);
    I[1] = 0.5 * (I[2]*I[2] - trace2(*this));
    I[0] = -det(*this);

    //    printf (" %lf x^3 +  %lf x^2 + %lf x + %lf = 0\n",
    //    	  I[3],I[2],I[1],I[0]);

    const long nbEigen = FindCubicRoots(I,v);

    std::sort(v, v+NSD, greater_abs());
    
    //    printf ("nbEigen = %d %12.5E %12.5E %12.5E\n",nbEigen,v[0],v[1],v[2]);

    double result[12];
    int nb_vec=0;
    while (true) {
	double a[9] = {pos[0][0]-v[nb_vec],pos[0][1],pos[0][2],
		       pos[1][0],pos[1][1]-v[nb_vec],pos[1][2],
		       pos[2][0],pos[2][1],pos[2][2]-v[nb_vec]};
	
	double eps = 1.e-3;
	int nb = 0;
	while (1)
	  {
	    nb = NullSpace(a,result,eps,NSD);
	    if (nb != 0) break;
	    eps *= 2.0;
	  }
	int kk=0;
	for (int i=nb_vec;i<nb+nb_vec;i++)
	  {
	    e[i] = mVector (result[0+kk*NSD],result[1+kk*NSD],result[2+kk*NSD]);
	    e[i].norm();
	    kk++;
	    if (i == 2) return nbEigen;
	  }
	nb_vec += nb;
	if (nb_vec == NSD) return nbEigen;
	if (nb > NSD)throw;
      }
    //    printf (" %lf x^3 +  %lf x^2 + %lf x + %22.15E = 0\n",
    //	  I[3],I[2],I[1],I[0]);
    throw;
  }


  inline long mTensor2::eigen2d(mVector e[NSD], double v[NSD]) const
  {            
    using std::abs;
    double a = 1.0;
    double b = - pos[0][0] - pos[1][1];
    double c = (pos[0][0] * pos[1][1] - pos[1][0] * pos[0][1]);

    e[2] = mVector(0,0,1);

    double delta = b*b - 4 * a * c;

    if (delta < 0) return 0;
    

    v[0] = (-b+sqrt(delta))/(2.*a);
    v[1] = (-b-sqrt(delta))/(2.*a);
    v[2] = 1.0;

    long nbEigen = 2;

    if (abs(v[1]) > abs(v[0]))
      {
	double temp = v[0];
	v[0] = v[1];
	v[1] = temp;
      }
    
    //    printf ("nbEigen = %d %12.5E %12.5E %12.5E\n",nbEigen,v[0],v[1],v[2]);

    double result[4];
    int nb_vec=0;
    while(1)
      {
	double a[4] = {pos[0][0]-v[nb_vec],pos[0][1],
		       pos[1][0],pos[1][1]-v[nb_vec]};	
	double eps = 1.e-8;
	int nb = 0;
	while (1)
	  {
	    nb = NullSpace (a,result,eps,2);
	    if (nb != 0)break;
	    eps *= 2.0;
	    //	    printf ("esp = %12.5E\n",eps);	
	  }
	int kk=0;
	for (int i=nb_vec;i<nb+nb_vec;i++)
	  {
	    e[i] = mVector (result[0+kk*2],result[1+kk*2],0.0);
	    e[i].norm();
	    kk++;
	  }
	nb_vec += nb;
	if (nb_vec == 2)return nbEigen;
	if (nb > 2)throw;
      }
    throw;
  }

/*
  mTensor2 mTensor2::transpose() const
  {
    mTensor2 tt;
    for (int i=0;i<NSD;i++)
      for (int j=0;j<NSD;j++)
	{
	  tt(i,j) = pos[j][i];
	}
    return tt;
  }
*/

  inline void mTensor2::transpose() {
    for (int i=0;i<NSD;i++)
      for (int j=i+1;j<NSD;j++) {
                double temp = pos[j][i];
                pos[j][i] = pos[i][j];
                pos[i][j] = temp;
          }
  }

  inline const mTensor2 mTensor2::operator!() const
  {
    mTensor2 B(*this);
    B.transpose();
    return B;
  }

  inline const mTensor2 mTensor2::intersect2d(mTensor2 &t) const
  {
    using std::abs;
    mVector e1[NSD],e2[NSD];
    double v1[NSD],v2[NSD];

    eigen2d(e1,v1);
    t.eigen2d(e2,v2);
        
    if (abs(v1[0]) > abs(v2[0]))
      {
	mTensor2 TR (e1[0],e1[1],e1[2]);
	mTensor2 TRT = TR;
	TRT.transpose();
	mTensor2 D (0.0);
	D(0,0) = v1[0];
	D(2,2) = 1.0;
	double x = (t * e1[1]) * e1[1];
	if(abs(x) > abs(v1[1]))
	  D(1,1) = x;
	else
	  D(1,1) = v1[1];	
	return (TR*D)*TRT;
      }
    else
      {
	mTensor2 TR (e2[0],e2[1],e2[2]);
	mTensor2 TRT = TR;
	TRT.transpose();
	mTensor2 D (0.0);
	D(0,0) = v2[0];
	D(2,2) = 1.0;
	double x = ((*this) * e2[1]) * e2[1];
	if(abs(x) > abs(v2[1]))
	  D(1,1) = x;
	else
	  D(1,1) = v2[1];	
	return (TR*D)*TRT;
      }
    
  }
  inline void mTensor2::set_all(const double &scalar)
    {
      for(int i=0;i<NSD;i++) 
	pos[i].set_all(scalar);
    }
#endif
