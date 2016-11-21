#ifndef VOLUME_REGIONS_H
#define VOLUME_REGIONS_H

#include <numeric>
#include <cmath>
#include "mPoint.h"
#include "mVector.h"
#include "mTensor2.h"


  {
namespace Util {

class EmptyRegion {};
class Tetrahedron;
class Bounding_Box;
class Sphere;
class Cylinder;

inline const mPoint  center(const Bounding_Box& b);
inline const mPoint  center(const Sphere& s);
inline const mPoint  center(const Tetrahedron& tet);
template <typename Pos> inline
const mPoint& center(const Pos& p) { return position(p);}

inline const double volume(const Bounding_Box& b);
inline const double volume(const Sphere& s);
inline const double volume(const Tetrahedron& tet);
inline const double volume(const EmptyRegion&) { return 0.0; }
inline const double volume(const mPoint&) { return 0.0; }

inline const Bounding_Box  bounding_box(const Tetrahedron& tet);
inline const Bounding_Box& bounding_box(const Bounding_Box& b);
inline const Bounding_Box  bounding_box(const Sphere& s);
inline const Bounding_Box  bounding_box(const mPoint& p);


inline const Sphere& bounding_sphere(const Sphere& s) { return s; }
inline const Sphere  bounding_sphere(const Bounding_Box& b);
inline const Sphere  bounding_sphere(const mPoint& p);




template <typename T> 
inline
bool is_empty(const T& region) {
  return volume(region) == 0.0;
}

template <typename PointIterator>
const Bounding_Box bounding_box_range(PointIterator beg,
				      PointIterator end);




inline
const mPoint componentwise_min(const mPoint& a, 
			       const mPoint& b) {
  mPoint result;
  using std::min;
  for (std::size_t dim = 0; dim != mPoint::NSD; ++dim) { 
    result[dim] = min(a[dim], b[dim]);
  }
  return result;
}

inline
const mPoint componentwise_max(const mPoint& a, 
			       const mPoint& b) {
  mPoint result;
  using std::max;
  for (std::size_t dim = 0; dim != mPoint::NSD; ++dim) { 
    result[dim] = max(a[dim], b[dim]);
  }
  return result;
}



template <typename Region1,
          typename Region2>
class Region_Union {
public:
  Region_Union(const Region1& r1,
	       const Region2& r2)
    : first(r1), second(r2) {}
  Region1 first;
  Region2 second;
};


template <typename R1,
          typename R2>
const Region_Union<R1,R2>
Union(const R1& r1, const R2& r2) { 
  return Region_Union<R1,R2>(r1, r2);
}

template <typename Region1,
	  typename Region2,
	  typename PositionType>
inline
bool region_contains(const Region_Union<Region1, Region2>& r,
		     const PositionType& p) {
  return (region_contains(r.first,  p) || 
	  region_contains(r.second, p));
}
  


template <typename Region1,
          typename Region2>
class Region_Intersection {
public:
  Region_Intersection(const Region1& r1,
		      const Region2& r2)
    : first(r1), second(r2) {}
  Region1 first;
  Region2 second;
};


template <typename R1,
          typename R2>
const Region_Intersection<R1,R2>
Intersect(const R1& r1, const R2& r2) { 
  return Region_Intersection<R1,R2>(r1, r2);
}

template <typename Region1,
	  typename Region2,
	  typename PositionType>
inline
bool region_contains(const Region_Intersection<Region1, Region2>& r,
		     const PositionType& p) {
  return (region_contains(r.first,  p) xor 
	  region_contains(r.second, p));
}
  

  
template <typename R1,
	  typename R2>
inline
const Bounding_Box
bounding_box(const Region_Union<R1,R2>& reg) {
  // The bounding box of the union of two regions
  // is the componentwise low corner of both
  // and high corner of both.
  return Bounding_Box(componentwise_min(bounding_box(reg.first).low(),
					bounding_box(reg.second).low()),
		      componentwise_max(bounding_box(reg.first).low(),
					bounding_box(reg.second).low()));
}
  




namespace unchecked {
inline
const mPoint position(const double v[mPoint::NSD]);

inline
mPoint& position(double v[mPoint::NSD]);
}




struct Cylinder {
  mPoint center;
  mVector axis_half_length;
  double radius;
  Cylinder(const mPoint& c,
	   const mVector& axis_hl,
	   double r) 
    : center(c), 
      axis_half_length(axis_hl),
      radius(r) {}
};


struct Tetrahedron {
  static const std::size_t SIZE = 4;
  std::size_t size() const { return SIZE; }
  mPoint corners[SIZE];
  mPoint& operator[](std::size_t i) {
    assert(i < SIZE);
    return corners[i];
  }
  const mPoint& operator[](std::size_t i) const { 
    assert(i < SIZE);
    return corners[i];
  }
};


struct Bounding_Box {
  typedef mPoint position_type;
  position_type minmax[2];
  Bounding_Box() { minmax[0] = minmax[1] = mPoint(); }
  template <typename Pos> 
  Bounding_Box(const Pos& p) { minmax[0] = minmax[1] = position(p); }
  Bounding_Box(const position_type& low,
	       const position_type& high) {
    minmax[0] = low;
    minmax[1] = high;
    for (size_t i = 0; i != position_type::NSD; ++i) {
      assert(low[i] <= high[i]);
    }
  }
  position_type& operator[](std::size_t i) {
    assert(i == 0 || i == 1);
    return minmax[i];
  }
  const position_type& operator[](std::size_t i) const {
    assert(i == 0 || i == 1);
    return minmax[i];
  }
  position_type&       low()        { return (*this)[0]; }
  position_type&       high()       { return (*this)[1]; }
  const position_type& low()  const { return (*this)[0]; }
  const position_type& high() const { return (*this)[1]; }
};


// These are to make it easy to make open and closed regions.
inline
double bigger_by_epsilon(const double x) {
  using std::abs;
  double guess = abs(x);
  while (x + guess / 2.0 > x)
    guess /= 2.0;
  return x + guess;
}

inline
double smaller_by_epsilon(const double x) {
  using std::abs;
  double guess = abs(x);
  while (x - guess * 0.5 < x)
    guess *= 0.5;
  return x - guess;
}


struct Bounding_Box {
  typedef mPoint position_type;
  position_type minmax[2];
  Bounding_Box() { minmax[0] = minmax[1] = mPoint(); }
  template <typename Pos> 
  Bounding_Box(const Pos& p) { minmax[0] = minmax[1] = position(p); }
  Bounding_Box(const position_type& low,
	       const position_type& high) {
    minmax[0] = low;
    minmax[1] = high;
    for (size_t i = 0; i != position_type::NSD; ++i) {
      assert(low[i] <= high[i]);
    }
  }
  position_type& operator[](std::size_t i) {
    assert(i == 0 || i == 1);
    return minmax[i];
  }
  const position_type& operator[](std::size_t i) const {
    assert(i == 0 || i == 1);
    return minmax[i];
  }
  position_type&       low()        { return (*this)[0]; }
  position_type&       high()       { return (*this)[1]; }
  const position_type& low()  const { return (*this)[0]; }
  const position_type& high() const { return (*this)[1]; }
};

inline
bool operator==(const Bounding_Box& a,
		const Bounding_Box& b) {
  return a.low() == b.low() && a.high() == b.high();
}


template <typename Pos>
inline
bool region_contains(const Bounding_Box& bb,
		     const Pos& pos) {
  for (size_t i = 0;
       i != mPoint::NSD;
       ++i) {
    if (position(pos)[i] <= bb.low()[i] ||
	position(pos)[i] >= bb.high()[i])
      return false;
  }
  return true;
}

const double signed_volume(const Tetrahedron& tet);


template <typename PointInSpace>
bool region_contains(const Tetrahedron& tet, 
		     const PointInSpace& pos) {
  struct Handedness {
    double operator()(const Tetrahedron& tet) const {
      if (cross(tet[1] - tet[0], tet[2] - tet[1]) * (tet[3] - tet[1]) > 0.0)
	return 1.0;
      return -1.0;	
    }
  } handedness;

  using std::size_t;
  const double tolerance = 0.0;
  // We are looking for particles within the (tetrahedral) region
  // defined by vertcoords
  const size_t NSD = mPoint::NSD;

  Tetrahedron xyztet;
  xyztet[0] = tet[0]; 
  xyztet[1] = tet[1]; 
  xyztet[2] = tet[2]; 
  xyztet[3] = position(pos);
  const double vol1 = signed_volume(xyztet) * handedness(tet);
  if (vol1 < tolerance) return false;
  
  xyztet[0] = tet[1];
  xyztet[1] = tet[3];
  xyztet[2] = tet[2];
  const double vol2 = signed_volume(xyztet) * handedness(tet);
  if (vol2 < tolerance) return false;
  
  xyztet[0] = tet[0]; 
  xyztet[1] = tet[2]; 
  xyztet[2] = tet[3]; 
  const double vol3 = signed_volume(xyztet) * handedness(tet);
  if (vol3 < tolerance) return false;
  
  xyztet[0] = tet[0]; 
  xyztet[1] = tet[3]; 
  xyztet[2] = tet[1]; 
  const double vol4 = signed_volume(xyztet) * handedness(tet);
  if (vol4 < tolerance) return false;
 
  return true;
}


inline
mPoint& position(mPoint& v) { return v; }

inline
const mPoint& position(const mPoint& v) { return v; }

// This is in its own namespace so 
// the user doesn't accidentally use it on just any array of doubles.
namespace unchecked {

inline
const mPoint position(const double v[mPoint::NSD]) { return mPoint(v); }

// We are about to do a cast, so we want to be very sure 
// that mPoint really is just NSD doubles in a row.
BOOST_STATIC_ASSERT((sizeof(mPoint) == mPoint::NSD * sizeof(double)));

inline
mPoint& position(double v[mPoint::NSD]) { 
  return *reinterpret_cast<mPoint*>(v); 
}

}
  

inline
const Tetrahedron 
make_tetrahedron(const mPoint& a,
		 const mPoint& b,
		 const mPoint& c,
		 const mPoint& d) {
  Tetrahedron result;
  result[0] = a;
  result[1] = b;
  result[2] = c;
  result[3] = d;
  return result;
}


template <typename Region>
class is_point_in {
public:
  typedef bool& result_type;

  is_point_in(Region& r) : reg_(r) {}

  template <typename PointInSpace>
  bool operator()(const PointInSpace& p) const {
    return region_contains(reg_, position(p));
  }
  const Region& region() const { return reg_; }
private:
  Region& reg_;
};

template <typename Region> 
is_point_in<const Region>
point_in(const Region& r) {
  return is_point_in<const Region>(r);
}


// We want to be able to ask for the bounding box of a 
// query region so a spacial data structure can 
// do a quick query.
template <typename Region>
inline
const Bounding_Box loose_bounding(const is_point_in<Region>& r) {
  return bounding_box(r.region());
}

// In general, we can use the bounding box of a thing.
template <typename T>
inline
const Bounding_Box loose_bounding_box(const T& x) {
  return bounding_box(x);
}




struct Sphere {
  mPoint center;
  double radius;
  Sphere(const mPoint& c, double rad = 0.0)
    : center(c), radius(rad) {}
};


template <typename Pos>
inline bool region_contains(const Sphere& s, const Pos& pos) {
  return (abs(position(pos) - s.center)
	  < 
	  s.radius);
}





inline const double signed_volume(const Tetrahedron& tet) {
  // Volume of arbitrary tet:
  // (1/6)·det(a-b, b-c, c-d)
  Util::mTensor2 edges(tet[1] - tet[0],
		       tet[2] - tet[0],
		       tet[3] - tet[0]);
  return det(edges)/6.0;
}

inline const double volume(const Tetrahedron& tet) {
  using std::abs;
  return abs(signed_volume(tet));
}

inline 
const mPoint
center(const Tetrahedron& tet) {
  // It's just the average position.
  return (std::accumulate(tet.corners,
			  tet.corners + tet.size(),
			  mPoint()) 
	  /
	  static_cast<double>(tet.size()));
};

inline 
const Bounding_Box
bounding_box(const Tetrahedron& tet) {
  return bounding_box_range(tet.corners, 
			    tet.corners + tet.size());
};

inline const mPoint 
center(const Bounding_Box& b) {
  return (b.minmax[0] + b.minmax[1]) / 2.0;
}

inline const double 
volume(const Bounding_Box& b) {
  const mVector diagonal = b.minmax[1] - b.minmax[0];
  
  // Return the product of the extent in each direction.
  const mVector::value_type vol =
    std::accumulate(diagonal.begin(), diagonal.end(), 
		    1.0, std::multiplies<mPoint::value_type>());
  return std::max(vol, 0.0); // Negative volume is the empty box.
}

inline 
const Bounding_Box&
bounding_box(const Bounding_Box& b) {
  return b;
}

inline const mPoint 
center(const Sphere& s) { 
  assert(s.radius >= 0.0);
  return s.center; 
}

inline const double 
volume(const Sphere& s) {
  assert(s.radius >= 0.0);
  return (4.0/3.0) * M_PI * std::pow(s.radius, 3);
}

inline 
const Bounding_Box
bounding_box(const Sphere& s) {
  return Bounding_Box(s.center - mVector(s.radius, s.radius, s.radius),
		      s.center + mVector(s.radius, s.radius, s.radius));
}




template <typename PointIterator>
const Bounding_Box bounding_box_range(PointIterator beg,
				      PointIterator end) {
  // Base case: nothing.
  if (beg == end) return Bounding_Box();

  // For starters, include the first point.
  Bounding_Box result(position(*beg), position(*beg));

  ++beg; // Continue from the second one...
  for (; beg != end; ++beg) {
    result.low()  = componentwise_min(result.low(),  position(*beg));
    result.high() = componentwise_max(result.high(), position(*beg));
  }
  return result;
}

inline
bool is_empty(const Region_Intersection<Sphere,Sphere>& r) {
  return (r.first.radius + r.second.radius
	  < 
	  abs(r.first.center - r.second.center));
}


inline
const Sphere bounding_sphere(const Bounding_Box& b) {
  return Sphere(center(b), abs(b.high() - b.low()));
}

inline const Bounding_Box  bounding_box(const mPoint& p) {
  return Bounding_Box(p,p); 
}



}} // End namespaces.



#endif // VOLUME_REGIONS_H
