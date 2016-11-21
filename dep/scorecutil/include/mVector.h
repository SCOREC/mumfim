/* 
   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of SCOREC_Util written and maintained by the 
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
// mVector.h: interface for the mVector class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _SCOREC_UTIL_MVECTOR_H_
#define _SCOREC_UTIL_MVECTOR_H_

#include <cmath>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <numeric>

  class mVector  
  {
  public:
    static const int NSD = 3;
    typedef double            value_type;
    typedef value_type*       iterator;
    typedef const value_type* const_iterator;
    
    static const mVector zeroVector() { return mVector(0.0, 0.0, 0.0); }

    explicit mVector()
      { std::fill(begin(), end(), 0.0); }
    explicit mVector(double x)
      // Fill all with x
      { std::fill(begin(), end(), x); }

    explicit mVector(double x, double y, double z)
      { pos[0]=x; pos[1]=y; pos[2]=z; }
    // We can make one from anything that has indexing up to NSD.
    template <typename ThreeVectorType>
    inline explicit mVector(const ThreeVectorType& X) {
      for (int i = 0; i < NSD; ++i)
	pos[i] = X[i];
    }
    inline ~mVector() {};

    // For compatability with things that think in terms of arrays:
    inline friend value_type*       data(mVector& v)       { return v.pos; }
    inline friend const value_type* data(const mVector& v) { return v.pos; }
    
    // For STL-like behavior.
    int            size() const  { return NSD; }
    iterator       begin()       { return pos; }
    iterator       end()         { return pos+NSD; }
    const_iterator begin() const { return pos; }
    const_iterator end()  const  { return pos+NSD; }

    // Indexing
    inline double  operator() (int i) const;
    inline double& operator() (int i);
    inline double& operator[] (int i);
    inline double  operator[] (int i) const;

    inline void set_all(const double & scalar);

    // Depricated:
    inline mVector& norm();
    inline double normValue();
    inline double mag() const;              // added by T. Bui 3/03     // Depricated
  private:
    double pos[NSD];
  };


  // Inner product
  inline double operator*(const mVector& a, const mVector &b)
  {
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
  }


  // This lets us avoid all sorts of nonsense with putting off 
  // taking the square root when doing abs(). This way you don't 
  // have to worry about performance nearly as much.
  class delayed_sqrt {
  public:
    explicit delayed_sqrt(double X) : x(X) {
      assert(x >= -0.0);
    }
    operator double() const { return sqrt(x); }
    friend inline bool operator<(const delayed_sqrt& s, const delayed_sqrt& t) {
      return s.x < t.x;
    }
    friend inline bool operator<(const delayed_sqrt& s, double X) {
      return s.x < X*X;
    }
    friend inline bool operator<(double X, const delayed_sqrt& s) {
      return X*X < s.x;
    }
    friend inline bool operator>(const delayed_sqrt& s, const delayed_sqrt& t) {
      return s.x > t.x;
    }
    friend inline bool operator>(const delayed_sqrt& s, double X) {
      return s.x > X*X;
    }
    friend inline bool operator>(double X, const delayed_sqrt& s) {
      return X*X > s.x;
    }
  private:
    double x;
  };

  inline delayed_sqrt abs(const mVector& v) {
    return delayed_sqrt(v*v);
  }

  // Scalar multiply
  inline const mVector operator*(const mVector& v, const double a) {
    //BOOST_STATIC_ASSERT(mVector::NSD == 3);
    return mVector(v[0] * a,
		   v[1] * a,
		   v[2] * a);
  }
  // Scalar left multiply.
  inline const mVector operator*(const double a, const mVector& v) { return v*a; }
  
  inline void mVector::set_all(const double & scalar) {
    for(int i=0;i<NSD;i++)
      pos[i] = scalar;
  }

  inline double mVector::mag() const { return abs(*this); }

  //jfadd
  inline double  mVector::normValue()
   {
     const double n = abs(*this);
     if(n==0.0) return n;
     pos[0]/=n;pos[1]/=n;pos[2]/=n;
     return n;
  }

  inline mVector&  mVector::norm()
   {
    double n = abs(*this);
    if(n==0.0)return *this;
    pos[0]/=n;pos[1]/=n;pos[2]/=n;
    return *this;
  }

   inline double & mVector::operator()(int i)
  {
#ifdef _DEBUG_
    if(i>=NSD)throw new mException (__LINE__,__FILE__,"wrong index");
#endif
    return pos[i];
  }

  inline double mVector::operator()(int i) const
  {
    return pos[i];
  }

  inline double& mVector::operator[](int i)
  {
#ifdef _DEBUG_
    if(i>=NSD)throw new mException (__LINE__,__FILE__,"wrong index");
#endif
    return pos[i];	    
  }

  inline double mVector::operator [] (int i) const
  {
#ifdef _DEBUG_
    if(i>=NSD)throw new mException (__LINE__,__FILE__,"wrong index");
#endif
    return pos[i];	    
  }

  inline const mVector operator-(const mVector& v) 
  {
    return mVector(-v[0], -v[1], -v[2]);
  }

  //jfadd const added
  inline const mVector operator-(const mVector& a, const mVector &b)
  {
    return mVector(a[0]-b[0],a[1]-b[1],a[2]-b[2]);
  }

  // cross product
  inline const mVector cross(const mVector& a, const mVector &b)
  {
    return mVector(a[1]*b[2]-a[2]*b[1],
		   a[2]*b[0]-a[0]*b[2],
		   a[0]*b[1]-a[1]*b[0]);
  }
  // Depricated:
  inline const mVector operator%(const mVector& a, const mVector &b) { return cross(a,b); }

  inline const mVector operator+(const mVector& a, const mVector &b)
  {
    return mVector(a[0]+b[0],a[1]+b[1],a[2]+b[2]);
  }


  inline const mVector operator/(const mVector& v, const double &other)
  {
    return v * (1.0/other);
  }

  inline mVector& operator+=(mVector& v, const mVector &other)
  {
    v = v + other;
    return v;
  }

  inline mVector& operator/=(mVector& v, const double &other)
  {
    v = v/other;
    return v;
  }

  inline mVector& operator-=(mVector& v, const mVector &other)
  {
    for(int i = 0; i < mVector::NSD; ++i)
      v[i] -= other[i];
    return v;
  }

  inline mVector& operator*=(mVector& v, const double &other)
  {
    for(int i = 0; i < mVector::NSD; ++i)
      v[i] *= other;
    return v;
  }


  inline bool operator==(const mVector& a, const mVector& b) {
    for(int i = 0; i < mVector::NSD; ++i)
      if(a[i] != b[i])
	return false;
    return true;
  }

  inline const mVector normalized(const mVector& v) {
    return v/abs(v);
  }
  inline double angleRad(const mVector& v1, const mVector &v2)
  {
    // have to copy not to modify original vectors
    mVector x(normalized(v1));
    mVector y(normalized(v2));
    double cosA = x * y;
    mVector cross = x % y;
    double sinA = cross.mag();
    return atan2(sinA,cosA);
  }

  // This isn't really necessary here.
  inline double angleDeg(const mVector& v1, const mVector &v2)
  {
    return angleRad(v1, v2) * (180.0/M_PI);
  }

  inline void swap(mVector& a, mVector& b) {
    std::swap_ranges(a.begin(), a.end(),
		     b.begin());
  }



  inline std::ostream& operator<<(std::ostream& os, const mVector& v) {
    os << v[0] << "\t" << v[1] << "\t" << v[2];
    return os;
  }

#endif 

