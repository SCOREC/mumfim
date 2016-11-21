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
// mPoint.h: interface for the mPoint class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _SCOREC_UTIL_MPOINT_H_
#define _SCOREC_UTIL_MPOINT_H_
#include <algorithm>

#include "mVector.h"

/**
   Simple class for 3D points
 */

class mPoint {
  /// Position
  double pos[3];
 public:
  typedef double value_type;
  static const std::size_t NSD = 3;
  typedef value_type* iterator;
  typedef const value_type* const_iterator;
  friend iterator       begin(mPoint& p)       { return p.pos; }
  friend iterator       end(mPoint& p)         { return p.pos + NSD; }
  friend const_iterator begin(const mPoint& p) { return p.pos; }
  friend const_iterator end(const mPoint& p)   { return p.pos + NSD; }

  inline mPoint() { std::fill(pos, pos + NSD, 0.0); }
  inline explicit mPoint(double x, double y =0.0, double z = 0.0) {
    pos[0] = x;
    pos[1] = y;
    pos[2] = z;
  }
  inline explicit mPoint(const double x[NSD]) {
    std::copy(x, x+NSD, pos);
  }
  mPoint(const mPoint& p) {
    std::copy(p.pos, p.pos+NSD, pos);
  }
  mPoint& operator=(const mPoint& p) {
    std::copy(p.pos, p.pos+NSD, pos);
    return *this;
  }
//  inline mPoint& operator+=(const mPoint &);
//  inline const mPoint  operator+(const mPoint &) const;
  inline const mVector operator-(const mPoint &) const;
//  inline const mPoint  operator*(double) const;
//  inline mPoint& operator*=(double);
  inline ~mPoint() {}
  inline double & operator()(int i)       { return (*this)[i]; }
  inline double   operator()(int i) const { return (*this)[i]; }
  inline double & operator[](int);
  inline double   operator[](int) const;

  inline bool operator<(const mPoint &other) const;
  /* Lexicographic sorting i.e. sort by x then y then z. This cannot
     be used for geometrical search but only for putting points in an
     associative container.
  */
  inline bool lexicographicLessThan (const mPoint &other, double EPS) const;

  friend inline double*       data(mPoint& p)       { return p.pos; }
  friend inline const double* data(const mPoint& p) { return p.pos; }
};

//inline const mPoint operator/(const mPoint& p, double x) {
//  return mPoint(p[0]/x, p[1]/x, p[2]/x);
//}

inline double & mPoint::operator[](int i) 
{
#ifdef _DEBUG_
	if(i>=NSD)throw new mException (__LINE__,__FILE__);
#endif
	return pos[i];
}

inline double mPoint::operator[](int i) const
{
#ifdef _DEBUG_
	if(i>=NSD)throw new mException (__LINE__,__FILE__);
#endif
	return pos[i];
}


namespace vector_operations_for_point {
// We should consider removing scalar operator* and operator/ because
// a point in affine space doesn't have these operations.  We could
// add mix(mPoint, mPoint, double) which would be an affine
// combination of two points. --BFD 1/5/2007
inline 
mPoint& operator+=(mPoint& a, mPoint const& b) {
  for (size_t i = 0; i != mPoint::NSD; ++i)
    a[i] += b[i];
  return a;
}

inline 
mPoint& operator-=(mPoint& a, mPoint const& b) {
  for (size_t i = 0; i != mPoint::NSD; ++i)
    a[i] -= b[i];
  return a;
}

inline 
mPoint const operator+(mPoint const& a, mPoint const& b) {
  mPoint tmp(a);
  tmp += b;
  return tmp;
}

inline 
mPoint& operator*=(mPoint& p, double x) {
  for (size_t i = 0; i != mPoint::NSD; ++i)
    p[i] *= x;
  return p;
}

inline 
mPoint& operator/=(mPoint& p, double x) {
  for (size_t i = 0; i != mPoint::NSD; ++i)
    p[i] /= x;
  return p;
}

inline 
mPoint const operator*(mPoint const& p, double x) {
  mPoint tmp(p);
  tmp *= x;
  return tmp;
}

inline 
mPoint const operator*(double x, mPoint const& p) { return p*x; }

inline 
mPoint const operator/(mPoint const& p, double x) {
  mPoint tmp(p);
  tmp /= x;
  return tmp;
}

}

// For now, use these operations anyway.
#ifndef POINTS_ARE_NOT_VECTORS
using namespace vector_operations_for_point;
#endif // POINTS_ARE_NOT_VECTORS


// It makes more sense to add a vector to a point.
inline 
const mPoint operator+(const mPoint& p, const mVector& v) {
  return mPoint(p[0]+v[0], p[1]+v[1], p[2]+v[2]);
}

// We can add in either order.
inline 
const mPoint operator+(const mVector& v, const mPoint& p) {
  return p + v;
}

// It makes sense to take the difference of a vector to a point.
inline 
const mPoint operator-(const mPoint& p, const mVector& v) {
  return mPoint(p[0]-v[0], p[1]-v[1], p[2]-v[2]);
}

// We can add in either order.
inline 
const mPoint operator-(const mVector& v, const mPoint& p) {
  return mPoint(v[0]-p[0], v[1]-p[1], v[2]-p[2]);
}

inline 
const mVector mPoint::operator-(const mPoint &other) const 
{
  return mVector(pos[0]-other.pos[0], 
		 pos[1]-other.pos[1], 
		 pos[2]-other.pos[2]);
}


inline 
bool operator==(const mPoint& a, const mPoint& b) {
  for (std::size_t i = 0; i != mPoint::NSD; ++i) {
    if (a[i] != b[i]) return false;
  }
  return true;
}


inline bool mPoint::lexicographicLessThan(const mPoint &other, double EPS) const
{
  if(pos[0]<other(0)-EPS) return 1;
  if(pos[0]>other(0)+EPS) return 0;
  if(pos[1]<other(1)-EPS) return 1;
  if(pos[1]>other(1)+EPS) return 0;
  if(pos[2]<other(2)-EPS) return 1;
  if(pos[2]>other(2)+EPS) return 0;
  return 0;
}

inline bool mPoint::operator<(const mPoint &other) const
{
  return lexicographicLessThan(other,1.e-6);
}


inline
const mPoint::value_type abs(const mPoint& p) {
  return sqrt(p[0] * p[0] + 
	      p[1] * p[1] + 
	      p[2] * p[2]);
}


inline std::ostream& operator<<(std::ostream& os, const mPoint& p) {
  os << p[0] << "\t" << p[1] << "\t" << p[2];
  return os;
}

#endif


