/*
 * VolumeBuckets is a simple class to store objects in 3D with O(1) access time (assuming a limit on density). It is provided with a bounding box and breaks that rectangular region up into equal-sized buckets. When an object is inserted, we find the bucket it  belongs in and add it to that bucket.
*
* The user can search for contained objects by either requesting a std::vector of pointers to all objects within a radius of a point. This allows us to do quick volume queries for neighbors of a point.
*/


#ifndef VOLUME_BUCKETS_H
#define VOLUME_BUCKETS_H

#include <vector>
#include <algorithm>
#include <iterator>
#include <boost/multi_array.hpp>
#include "mVector.h"

template <typename IterOutter, 
	  typename IterInner = typename std::iterator_traits<IterOutter>::value_type::iterator>
class NestedIterator {
  typedef IterOutter outter_iterator_type;
  typedef typename std::iterator_traits<outter_iterator_type>::value_type outter_value_type;
  typedef IterInner inner_iterator_type;
  typedef typename std::iterator_traits<inner_iterator_type>::value_type inner_value_type;
public:
  typedef std::bidirectional_iterator_tag iterator_category;
  typedef inner_value_type value_type;
  typedef value_type* pointer;
  typedef void difference_type;
  typedef value_type& reference;

  NestedIterator(outter_iterator_type b,
		 outter_iterator_type e) : bi(b), bend(e) {
    if (bi != bend) {
      vbegin = b->begin();
      vend   = b->end();
      vi = vbegin;
      }
  }
  reference operator*() const { return *vi; }
  
  friend
  inline bool operator!=(const NestedIterator& a, const NestedIterator& b) {
    return !(a == b);
  }
  friend
  inline bool operator==(const NestedIterator& a, const NestedIterator& b) {
    if (a.bi != b.bi) return false;
    // They are on the same bucket.
    if (a.bi == a.bend) { // It's at the end.
      assert(a.bend == b.bend); // We assume they have the same end.
      return true;
    }
    // Same bucket, not at the overall end.
    assert(a.vbegin == b.vbegin);
    assert(a.vend == b.vend);
    return a.vi == b.vi;
  }
  NestedIterator& operator++() {
    ++vi;
    if (vi == vend) {
      ++bi;
      if (bi != bend) {
	vbegin = bi->begin();
	vend= bi->end();
	vi = vbegin;
      }
    }
    return *this;
  }
  const NestedIterator operator++(int) {
    NestedIterator tmp = this;
    ++(*this);
    return tmp;
  }
  
  
private:
  outter_iterator_type bi;
  outter_iterator_type bend;
  inner_iterator_type  vi;
  inner_iterator_type  vbegin;
  inner_iterator_type  vend;
};




// This class stores objects in 3D with constant-time access by position.
// It stores objects by value.
// PosType foo must be addressible by double foo[i] for i = 1, 2, 3 to find its location in 3-space.
// If you want to use this class to store something else, you will have to wrap it in another class.

template <typename PosType>
class VolumeBuckets {
 public:
  static const int NSD = 3;
  typedef PosType     value_type;
  typedef value_type& reference;
  typedef value_type* pointer;
  typedef const value_type& const_reference;
  typedef std::vector<value_type> bucket_type;

  template <typename ThreeVectorType>
  VolumeBuckets(ThreeVectorType inMin, ThreeVectorType inMax) :
    buckets_(boost::extents[side_][side_][side_]) {
    for(int dim = 0; dim < NSD; ++dim) {
      min_extent[dim] = inMin[dim];
      max_extent[dim] = inMax[dim];
    }
  }



  inline void insert(const PosType& p) {
    get_bucket(p).push_back(p);
  }
  // A version for use with std::insert_iterator:
  typedef int iterator;
  inline iterator insert(iterator unused, const PosType& p) {
    get_bucket(p).push_back(p);
    return unused;
  }
  
  inline pointer find(const PosType& p) {
    bucket_type& b = get_bucket(p);
    typename bucket_type::iterator i = std::find(b.begin(), b.end(), p);
    if (i == b.end())
      return NULL;
    return &*i;
  }

  inline void move(const PosType& from, const PosType& to) {
    const pointer found = find(from);
    assert(found);
    erase(found); // This could be optimized for the case of small motion.
    insert(to);
  }

  typename bucket_type::iterator
  erase(pointer p) {
    bucket_type& b = get_bucket(*p);
    assert(!b.empty());
    const size_t p_offset = p - &b[0];
    assert(p >= &b[0]); // p must be in the array...
    assert(p_offset < b.size()); // Assert p is there.
    const typename bucket_type::iterator pi = b.begin() + p_offset;
    return b.erase(pi);
  }

//  void remove(const PosType& p) {
//    get_bucket(p).push_back(p);
//  }

  template <typename ThreeVectorType>
  bucket_type& get_bucket(ThreeVectorType pos) {
    int theIndex[NSD] = {-1, -1, -1};
    get_index_for_point(pos, theIndex);
    return buckets_[theIndex[0]][theIndex[1]][theIndex[2]];
  }
  
  template <typename ThreeVectorType>
  std::vector<bucket_type*> 
  get_all_buckets_with_points_within_r_of_pos(double r, ThreeVectorType pos);

  template <typename ThreeVectorType>
  std::vector<value_type*> 
  get_all_objects_with_points_within_r_of_pos(double r, ThreeVectorType pos);

  template <typename ThreeVectorType>
  void get_index_for_point(ThreeVectorType pos, int* index) const {
    for(int dim = 0; dim < NSD; ++dim) {
      assert(pos[dim] >= min_extent[dim]);
      assert(pos[dim] <= max_extent[dim]);

      // The index in the dim direction should range from 0 to side-1...
      // First normalize to between 0.0 and 1.0 ...
      const double theNormalizedPosition = 
	(pos[dim] - min_extent[dim]) / (max_extent[dim] - min_extent[dim]);
      // Then scale by the number of buckets in that direction. 
      // (Also scale by 0.999 to keep the range in [0, side-1] rather than having 1.0 go to "side".)
      index[dim] = 
	static_cast<int>(theNormalizedPosition * 
			 0.999 * 
			 double(buckets_.shape()[dim]));
      
      assert(index[dim] >= 0);
      assert(index[dim] < buckets_.shape()[dim]);
    }
  }    

  typedef boost::multi_array<bucket_type, NSD> buckets_type;
  //  typedef NestedIterator<typename buckets_type::iterator, 
  //			 typename bucket_type::iterator> iterator;
  /*  iterator begin() {
    return iterator(buckets_.begin(), buckets_.end());
  }
  iterator end() {
    return iterator(buckets_.end(), buckets_.end());
    }*/
private:
  double min_extent[NSD];
  double max_extent[NSD];
  static const int side_ = 50; // Later this will want to depend on the extent and on the cutoff radius.
  buckets_type buckets_;
};



template<typename T>
template <typename ThreeVectorType>
std::vector<typename VolumeBuckets<T>::bucket_type*> 
VolumeBuckets<T>::
get_all_buckets_with_points_within_r_of_pos(double r, ThreeVectorType pos) {
  // Figure out the indicies for the buckets in the rectangular solid...
  int theMinIndex[NSD];
  int theMaxIndex[NSD];
  double theMinPos[NSD];
  for(int dim = 0; dim < NSD; ++dim) {
    theMinPos[dim] = std::max(pos[dim] - r, min_extent[dim]);
  }
  get_index_for_point(theMinPos, theMinIndex);

  double theMaxPos[NSD];
  for(int dim = 0; dim < NSD; ++dim) {
    theMaxPos[dim] = std::min(pos[dim] + r, max_extent[dim]);
  }
  get_index_for_point(theMaxPos, theMaxIndex);

  // So now theMinIndex and theMaxIndex define a bounding box of indicies for the buckets we want.
  // (Actually, there may be some buckets outside the sphere of r at the corners, but for now, who cares?)
  std::vector<bucket_type*> theResult;
  for(int zi = theMinIndex[2]; zi <= theMaxIndex[2]; ++zi) {
    for(int yi = theMinIndex[1]; yi <= theMaxIndex[1]; ++yi) {
      for(int xi = theMinIndex[0]; xi <= theMaxIndex[0]; ++xi) {
	theResult.push_back(&buckets_[xi][yi][zi]);
      }
    }
  }
  return theResult;
}

template<typename T>
template <typename ThreeVectorType>
std::vector<typename VolumeBuckets<T>::value_type*> 
VolumeBuckets<T>::
get_all_objects_with_points_within_r_of_pos(double r, ThreeVectorType pos) {
  const double r_squared = r*r;
  const std::vector<bucket_type*> theBuckets = 
    get_all_buckets_with_points_within_r_of_pos(r, pos);
  std::vector<value_type*> theResult;
  for (typename std::vector<bucket_type*>::const_iterator i = theBuckets.begin();
       i != theBuckets.end();
       ++i) {
    for (typename bucket_type::iterator j = (*i)->begin();
	 j != (*i)->end();
	 ++j) {
      using SCOREC::Util::mVector;
      const mVector theDistance(mVector(pos) - mVector(*j));
      if (theDistance*theDistance < r_squared)
	theResult.push_back(&*j); // Add a pointer to that object.
    }
  }
  return theResult;
}


// DistanceOrdering is a function object.
// It isn't used by VolumeBuckets but may be helpful to users of it.
// It is constructed with a point and its operator() tells which 
// argument is closer to the point.
//
// Usage:
//   Given a container of mVectors, calling:
//     std::sort(points.begin(), points.end(),
//               DistanceOrdering(x));
//   will sort points by their distance from point x.
struct DistanceOrdering {
  const mVector x;
  DistanceOrdering(const mVector& inX) : x(inX) {}
  inline bool operator()(const mVector& a, const mVector& b) {
    return abs(a-x) < abs(b-x);
  }
};



  /**** Usage *****

  VolumeBuckets theBuckets;
  for all particles p {
  theBuckets.add(p);
  }
  ...
  // to make neighbor lists...
for all particles p {
  vector<vector<particle>*> nearby_buckets = get_all_buckets_with_points_within_r_of_pos(neighbor_radius, p.pos());
  for(int i = 0; i < nearby_buckets.size(); ++i) {
    for(j = 0; j < nearby_buckets[i]->size(); ++j) {
      // The same old distance check between (*nearby_buckets[i])[j] and p...
    }
  }
}

*******/
#endif // VOLUME_BUCKETS_H
