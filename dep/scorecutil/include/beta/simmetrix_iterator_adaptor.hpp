#ifndef SIMMETRIX_ITERATOR_ADAPTOR_HPP
#define SIMMETRIX_ITERATOR_ADAPTOR_HPP

// Simmetrix provides many "iterators" that are not 
// C++-style iterators. This file makes it easy to work
// with these iterators in an STL-like way.

// The Simmetrix iterators have an interface like this:
// LatticeModel_EXPORT LMIter LatticeModel_Iter(pLatticeModel latModel);
// LatticeModel_EXPORT void LMIter_delete(LMIter iter);
// LatticeModel_EXPORT pLattice LMIter_next(LMIter iter);
//
// Usage:
// pLatticeModel m;
// for (iterator<LMIter> i = begin(m);
//      i != end(m);
//      ++i) {
//   // Now *i as a pLattice.
//   ...
// }
//


#include <algorithm>
#include <iterator>
#include <boost/array.hpp>
#include "SimLattice.h"
#include "SimModel.h"
#include "nested_iterator.h"

  {
namespace Util {
namespace Simmetrix {

template <typename T>
class iterator;

template <typename T>
struct simmetrix_iterator_traits;

  /*
template <>
struct simmetrix_iterator_traits<LMIter> {
  typedef LMIter const    base_iterator_type;
  typedef pLattice        value_type; // The type we iterate over
  typedef pLattice&       reference;         
  typedef pLattice        pointer;
  static pointer next(base_iterator_type i) {
    return LMIter_next(i);
  }
  static void destroy(base_iterator_type i) {
    return LMIter_delete(i);
  }
  static base_iterator_type clone(base_iterator_type i) {
    return LMIter_clone(i);
  }
  static reference deref(pointer p) {
    
  };
};

template <>
struct simmetrix_iterator_traits<LIter> {
  typedef LIter base_iterator_type;
  typedef LIter value_type;
  typedef LIter reference;
  typedef LIter pointer;
  static pointer next(base_iterator_type i) {
    return (LIter_next(i) ? i : NULL);
  }
  static void destroy(base_iterator_type i) {
    return LIter_delete(i);
  }
  static base_iterator_type clone(base_iterator_type i) {
    return LIter_clone(i);
  }
  static 
};

  */
  /**** These won't work without a clone operation.
template <>
struct simmetrix_iterator_traits<VIter> {
  typedef VIter    base_iterator_type;
  typedef pVertex    value_type;
  static value_type next(base_iterator_type i) {
    return VIter_next(i);
  }
  static void destroy(base_iterator_type i) {
    return VIter_delete(i);
  }
};

template <>
struct simmetrix_iterator_traits<EIter> {
  typedef EIter    base_iterator_type;
  typedef pEdge    value_type;
  static value_type next(base_iterator_type i) {
    return EIter_next(i);
  }
  static void destroy(base_iterator_type i) {
    return EIter_delete(i);
  }
};

template <>
struct simmetrix_iterator_traits<FIter> {
  typedef FIter    base_iterator_type;
  typedef pFace    value_type;
  static value_type next(base_iterator_type i) {
    return FIter_next(i);
  }
  static void destroy(base_iterator_type i) {
    return FIter_delete(i);
  }
};

template <>
struct simmetrix_iterator_traits<RIter> {
  typedef RIter    base_iterator_type;
  typedef pRegion  value_type;
  static value_type next(base_iterator_type i) {
    return RIter_next(i);
  }
  static void destroy(base_iterator_type i) {
    return RIter_delete(i);
  }
  }; */




template <>
class iterator<LMIter> {
  typedef LMIter T; // For shorthand.
public:
  typedef pLattice value_type;
  typedef pLattice reference; // It's already a pointer, we don't need to pass it by reference.
  typedef value_type* pointer;
  typedef std::size_t difference_type;
  typedef std::forward_iterator_tag iterator_category;


  explicit iterator(T inIter = NULL) 
    : iter_(inIter), curr_(NULL) {
    if (iter_) {
      curr_ = LMIter_next(iter_);
    }
  }

  iterator(const iterator& other)
  : iter_(LMIter_clone(other.iter_)), 
    curr_(other.curr_)
  {
  }
  iterator& operator=(const iterator& other) {
    iterator tmp(other); // Copy tmp.
    swap(*this, tmp);
    // What we owned is now automatically deleted.
  }
  ~iterator() {
    if (iter_)
     LMIter_delete(iter_);
  }
  reference operator*() const {
    return curr_;
  }
  iterator& operator++() {
    assert(curr_);
    assert(iter_);
    curr_ = LMIter_next(iter_);
    return *this;
  }
  friend 
  bool operator==(const iterator& a,
		  const iterator& b) {
    return a.curr_ == b.curr_;
  }
  friend
  void swap(iterator& a, iterator& b) {
    std::swap(a.iter_, b.iter_);
    std::swap(a.curr_, b.curr_);
  }
private:
  LMIter iter_;
  pLattice curr_;
};

template <>
class iterator<LIter> {
  typedef LIter T; // For shorthand.
public:
  typedef LIter value_type;
  typedef LIter reference; // It's already a pointer, we don't need to pass it by reference.
  typedef LIter* pointer;
  typedef std::size_t difference_type;
  typedef std::forward_iterator_tag iterator_category;


 // Takes ownership of an iterator.
  explicit iterator(LIter inIter = NULL) 
    : iter_(inIter)
  {
    if (iter_) {
      if (!LIter_next(iter_)) {
	// It an empty range.
	off_end_ = true;
      }
    }
  }

  // Copy constructor
  iterator(const iterator& other)
  : iter_(LIter_clone(other.iter_)), 
    off_end_(other.off_end_)
  {
  }
  iterator& operator=(const iterator& other) {
    iterator tmp(other); // Copy tmp.
    swap(*this, tmp);
    // What we owned is now automatically deleted.
  }
  ~iterator() {
    if (iter_)
      LIter_delete(iter_);
  }
  reference operator*() const {
    assert(!off_end_);
    return iter_;
  }
  iterator& operator++() {
    assert(iter_);
    assert(!off_end_);
    off_end_ = LIter_next(iter_);
    return *this;
  }
  friend 
  bool operator==(const iterator& a,
		  const iterator& b) {
    // If either is off the end, they are equal iff both are off the
    // end.
    if (a.off_end_ || b.off_end_)
      return a.off_end_ == b.off_end_;
    return (LIter_LatticeIndex1(a.iter_) == LIter_LatticeIndex1(b.iter_) &&
	    LIter_LatticeIndex2(a.iter_) == LIter_LatticeIndex2(b.iter_) &&
	    LIter_LatticeIndex3(a.iter_) == LIter_LatticeIndex3(b.iter_) &&
	    LIter_LatticeIndex4(a.iter_) == LIter_LatticeIndex4(b.iter_));
  }
  friend
  void swap(iterator& a, iterator& b) {
    std::swap(a.iter_, b.iter_);
    std::swap(a.off_end_, b.off_end_);
  }
private:
  LIter iter_;
  bool off_end_;
};




template <typename T, typename U>
bool operator!=(const T& a,
		const U& b) {
  return !(a == b);
}


typedef iterator<LMIter> lattice_model_iterator; // Iterate over crystals.
typedef iterator<LIter>  lattice_iterator;   // Iterate over lattice sites.

lattice_model_iterator
begin(pLatticeModel x) {
  return lattice_model_iterator(LatticeModel_Iter(x));
}

lattice_model_iterator
end(pLatticeModel x) {
  return lattice_model_iterator(NULL);
}

iterator<LIter>  //lattice_iterator
begin(pLattice l, const double min[3], const double max[3]) {
  // Simmetrix tools aren't const correct.
  double theMin[3];
  double theMax[3];
  std::copy(min, min+3, theMin);
  std::copy(max, max+3, theMax);
  return lattice_iterator(Lattice_Iter(l, theMin, theMax));
}

lattice_iterator
end(pLattice l, const double min[3], const double max[3]) {
  return lattice_iterator(NULL);
}





// This traits class tells us what atoms we are looking for.

struct lattice_model_nested_traits {
  typedef boost::array<double, 3> position_type;
  lattice_model_nested_traits() {
    // By default, iterate over everything.
    std::fill(min_pos.begin(), min_pos.end(),
	      -std::numeric_limits<double>::infinity());
    std::fill(max_pos.begin(), max_pos.end(),
	      std::numeric_limits<double>::infinity());
  }
  lattice_model_nested_traits(const position_type& inMin,
			      const position_type& inMax)
    : min_pos(inMin),
      max_pos(inMax) {}
  position_type min_pos;
  position_type max_pos;
  typedef pLattice outer_reference_type;
  typedef lattice_iterator  iterator;
  iterator begin(outer_reference_type c) const {
    using Simmetrix::begin;
    return begin(c, min_pos.data(), max_pos.data());
  }
  iterator end  (outer_reference_type c) const { 
    using Simmetrix::end;
    return end(c, min_pos.data(), max_pos.data());
  }
};

typedef nested_iterator<lattice_model_iterator, 
			lattice_model_nested_traits> all_lattice_iterator;

  /*
struct all_lattice_iterator {
  typedef LatticeIter*const value_type;
  typedef value_type const& reference;
  all_lattice_iterator(lattice_model_iterator beg,
		       lattice_model_iterator end) 
    : oi_(beg), oend_(end),
      li_(NULL), lend_(NULL) {
    if (oi_ != oend_) {
      // Set the inside range.
      double lo[] = { 0, 0, 0};
      double hi[] = { 3, 3, 3};
      li_   = Simmetrix::begin(*oi_, lo, hi);
      lend_ = Simmetrix::end(*oi_, lo, hi);
      // Make sure we find a non-empty inside range.
      while (li_ == lend_ &&
	     oi_ != oend_) {
	++oi_; // Increment outter.
	if (oi_ != oend_) {
	  li_   = Simmetrix::begin(*oi_, lo, hi);
	  lend_ = Simmetrix::end(*oi_, lo, hi);
	}
      }      
    }
  }

  reference operator*() const { return *li_; }
  void operator++() {
      double lo[] = { 0, 0, 0};
      double hi[] = { 3, 3, 3};
    if (oi_ != oend_) {
      ++li_; // Increment inner
      // If we ran out of inners, look for the next one.
      while (li_ == lend_ &&
	     oi_ != oend_) { // 
	++oi_; // Increment outter.
	if (oi_ != oend_) {
	  li_   = begin(*oi_, lo, hi);
	  lend_ = end(*oi_, lo, hi);
	}
      }      
    }
  }
private:
  lattice_model_iterator oi_;
  lattice_model_iterator oend_;
  lattice_iterator li_;
  lattice_iterator lend_;
};

// Iterate over all latice points in a lattice model.
// Use nested_iterator to do this.
// Nested_iterator needs a range getter.
struct lattice_range_getter {
  typedef lattice_model_iterator::reference outer_reference_type;
  typedef lattice_model_iterator::value_type outer_value_type;
  typedef lattice_iterator     iterator;
  static iterator begin(pLattice c) { 
    double low[3] = { 0.0, 0.0, 0.0 };
    double high[3] = { 5.0, 5.0, 5.0};
    return ::SCOREC::Util::Simmetrix::begin(c, low, high); 
  }
  static iterator end  (pLattice c) { 
    double low[3] = { 0.0, 0.0, 0.0 };
    double high[3] = { 5.0, 5.0, 5.0};
    return ::SCOREC::Util::Simmetrix::end(c, low, high);  
  }
};
  */


  /*
iterator<VIter>
begin_vertices(pMesh m) {
  return iterator<VIter>(M_vertexIter(m));
}

iterator<VIter>
end_vertices(pMesh m) {
  return iterator<VIter>(NULL);
}

iterator<EIter>
begin_edges(pMesh m) {
  return iterator<EIter>(M_edgeIter(m));
}

iterator<EIter>
end_edges(pMesh m) {
  return iterator<EIter>(NULL);
}

iterator<FIter>
begin_faces(pMesh m) {
  return iterator<FIter>(M_faceIter(m));
}

iterator<FIter>
end_faces(pMesh m) {
  return iterator<FIter>(NULL);
}

iterator<RIter>
begin_regions(pMesh m) {
  return iterator<RIter>(M_regionIter(m));
}

iterator<RIter>
end_regions(pMesh m) {
  return iterator<RIter>(NULL);
  }*/


}}} // End namespaces.

#endif // SIMMETRIX_ITERATOR_ADAPTOR_HPP
