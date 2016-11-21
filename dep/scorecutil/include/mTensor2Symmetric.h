/* 
   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of SCOREC written and maintained by the 
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
// mTensor2Symmetric.h: interface for the mTensor2Symmetric class.
//
//////////////////////////////////////////////////////////////////////
#ifndef _MTensor2Simmetric_H_
#define _MTensor2Simmetric_H_

#include "mVector.h"
#include <algorithm>

  class mTensor4;
  /**
     Simple class for a 2nd order tensor in 3d which
     has a 3x3 matrix representation.
  */

  class mTensor2Symmetric
  {
    mVector pos[3]; // Just use the upper triangle of these
  public:
    mTensor2Symmetric(const mVector &c1, const mVector &c2, const mVector &c3);
    mTensor2Symmetric(double init = 0.0);
    const mTensor2Symmetric operator*(const mTensor4 &other) const;
    virtual ~mTensor2Symmetric();

    inline const mVector operator*(const mVector &other) const {
      mVector m(0,0,0);
      for (int i=0; i < other.size(); ++i) {
	m[i] += pos[i] * other;
      }
      return m;
    }
    inline const mTensor2Symmetric operator*(const mTensor2Symmetric &other) const
    {
      mTensor2Symmetric m(0);
      for(int i=0;i<3;i++) {
	for(int j=0;j<3;j++) {
	  for(int l=0;l<3;l++) {
	    m(i,j) += (*this)(i,l) * other(l,j);
	  }
	}
      }
      return m;
    }
    inline mTensor2Symmetric &operator*=(const double &scalar)
    {
      for(int i=0;i<3;i++) {
	for(int j=0;j<3;j++) {
	  (*this)(i,j) *= scalar;
	}
      }
      return *this;
    }

    inline mTensor2Symmetric &operator/=(const double& scalar) {
      for(int i=0;i<3;i++) {
	for(int j=0;j<3;j++) {
	  (*this)(i,j) /= scalar;
	}
      }
      return *this;
    }

    /*    inline mTensor2Symmetric operator = (const double &scalar) const
    {
      mTensor2Symmetric other(scalar);
      return other;
      }*/

    inline mTensor2Symmetric &operator+=(const mTensor2Symmetric &other)
    {
      for(int i=0;i<3;i++) {
	for(int j=0;j<3;j++) {
	  (*this)(i,j) += other.pos[i][j];
	}
      }
      return *this;
    }

    inline mTensor2Symmetric& operator-=(const mTensor2Symmetric &other)
    {
      for(int i=0;i<3;i++) {
	  for(int j=0;j<3;j++) {
	    (*this)(i,j) -= other.pos[i][j];
	  }
      }
      return *this;
    }

    inline const mTensor2Symmetric operator*(const double &scalar) const
    {
      mTensor2Symmetric other;
      for(int i=0;i<3;i++) {
	for(int j=0;j<3;j++) {
	  other.pos[i][j] = scalar * (*this)(i,j);
	    }
	}
      return other;
    }
    double & operator() ( int, int );
    double operator() ( int, int ) const;

    const mVector operator[]( int i ) const;       
    const mVector col( int j ) const;                   

    const mTensor2Symmetric operator+(const mTensor2Symmetric &other) const;
    /// solve the system
    const mVector operator/(const mVector & other) const;
    const mTensor2Symmetric invert() const;
    /// determinant
    long eigen   (mVector e[3],double v[3]) const;
    long eigen2d (mVector e[3],double v[3]) const;
    inline void symmetrize ();
    mTensor2Symmetric intersect2d (mTensor2Symmetric &t) const;
    // mTensor2Symmetric transpose() const;
    inline void set_all(const double & scalar);

    void transpose();
    mTensor2Symmetric operator!() const; 
  };


  double det(const mTensor2Symmetric& t);
  double trace(const mTensor2Symmetric& t);
  double trace2(const mTensor2Symmetric& t);


  inline void mTensor2Symmetric::symmetrize()
  {
    pos[0][1] = pos[1][0] = 0.5 * (pos[0][1] + pos[1][0]);
    pos[0][2] = pos[2][0] = 0.5 * (pos[0][2] + pos[2][0]);
    pos[2][1] = pos[1][2] = 0.5 * (pos[2][1] + pos[1][2]);
  }

  inline void mTensor2Symmetric::set_all(const double & scalar) {
    for(int i=0;i<NSD;i++)
      pos[i].set_all(scalar);
  }

  inline mTensor2Symmetric::mTensor2Symmetric(const mVector &c0, 
					      const mVector &c1, 
					      const mVector &c2)
  {
    pos[0] = c0; pos[1] = c1; pos[2] = c2;
  }

  inline double &mTensor2Symmetric::operator () (int i,int j)
  {
    return pos[std::min(i,j)][std::max(i,j)];
  }

  inline double mTensor2Symmetric::operator () (int i,int j) const
  {
    return pos[std::min(i,j)][std::max(i,j)];
  }

  inline const mVector mTensor2Symmetric::operator[]( int i ) const
  {
    //    return the ith row
    mVector theRow = pos[i];
    for(int j = 0; j < i; ++j) {
      theRow[j] = (*this)(i,j);
    }
    return theRow;
  }

  inline const mVector mTensor2Symmetric::col( int j ) const
  {
    return mVector(pos[0][j], pos[1][j], pos[2][j]);
  }

  inline const mTensor2Symmetric mTensor2Symmetric::operator+(const mTensor2Symmetric &other) const
  {
    mTensor2Symmetric m(0);
    for (int i = 0; i < 3; ++i) {
      m.pos[i] += pos[i] + other[i];
    }
    return m;
  }    
  inline const mVector operator*(const mVector &other,const mTensor2Symmetric &t)
  {
    mVector m(0,0,0);
    for(int i=0;i<3;i++) {
      for(int j=0;j<3;j++) {
	m[i] += t(i,j) * other[j];
      }
    }
    return m;
  }

#endif // _MTensor2Simmetric_H_
