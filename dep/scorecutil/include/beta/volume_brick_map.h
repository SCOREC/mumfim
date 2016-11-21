#ifndef VOLUME_BRICK_MAP_H
#define VOLUME_BRICK_MAP_H

#include <map>
#include <vector>
#include <boost/array.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "beta/nested_iterator.h"
#include "beta/volume_regions.h"

// This class is to store points in space where
// the density, rho, is is known a-priori and and where the density of a region will be either almost zero or almost rho.
// It uses a std::map to search for a coordinate and a std::vector to store them.

  {
namespace Util {

// A function object to make it easy
// to iterate over std::map keys.
template <typename T>
struct select1st {
  typedef typename boost::add_reference<T>::type argument_type;
  typedef typename boost::remove_reference<T>::type::first_type result_type;

  inline
  result_type operator()(argument_type x) const { return x.first; }

};

template <typename I, typename O, typename Op>
void copy_if(I beg, I end, O out, Op op) {
  for (; beg != end; ++beg) {
    if (op(*beg)) {
      *out++ = *beg;
    }
  }
}


template <typename Op1, typename Op2>
struct unary_compose {
  Op1 o1_;
  Op2 o2_;
  typedef typename Op2::argument_type argument_type;
  typedef typename Op1::result_type result_type;
  unary_compose(Op1 one, Op2 two) 
    : o1_(one), o2_(two) {}
  result_type operator()(argument_type arg) {
    return o1_(o2_(arg));
  }
};

template <typename Op1, typename Op2>
unary_compose<Op1, Op2> compose1(Op1 one, Op2 two) {
  return unary_compose<Op1, Op2>(one, two);
}


// This way you can iterate over keys of maps, for example like so:
// typedef pair_first_iterator<std::map<int, double>::iterator >::type map_key_iter;
// copy(map_key_iter(theMap.begin()), map_key_iter(theMap.end()),
//      some_target_iterator);

template <typename Iter>
struct pair_first_iterator {
  typedef 
  boost::transform_iterator<select1st<typename boost::iterator_reference<Iter>::type>,
			    Iter> 
  type;
};



template <typename PointKey,
          typename Value>
class volume_brick_map {
public:
  static const unsigned int  NSD = 3;
  typedef PointKey                    key_type;
  typedef std::pair<PointKey, Value>  value_type;
private:
  typedef std::pair<PointKey, Value>  internal_value_type;
  typedef std::vector<internal_value_type>   internal_type;
  typedef boost::array<int, NSD>    outer_key_type;
  typedef std::map<outer_key_type, internal_type> map_type;
  typedef std::vector<internal_type*> found_brick_vector;
  typedef typename pair_first_iterator<typename internal_type::iterator>::type key_iterator;
  map_type map_;
public:
  typedef nested_iterator<typename container_iterator<map_type>::type,
			  second_range<typename map_type::value_type,
				       typename internal_type::iterator> > iterator;

  explicit volume_brick_map(const key_type& inBrickSize) : brick_size(inBrickSize) {}
  iterator begin() { return iterator(map_.begin(), map_.end()); }
  iterator end()   { return iterator(map_.end(),   map_.end()); }

  std::size_t num_bricks() const { return map_.size(); }

  void insert(const value_type& v) {
    map_[get_map_key(v.first)].push_back(v);
  }

  // Find the point.
  iterator find(const key_type& k) {
    typename map_type::iterator brick = map_.find(get_map_key(k));
    if (brick != map_.end()) {
      key_iterator beg = key_iterator(brick->second.begin()); 
      key_iterator end = key_iterator(brick->second.end()); 

      return iterator(map_.begin(), map_.end(),
		      std::find(beg,end, k).base());
    }
    return (*this).end();
  }
 // Copy exactly the points in the region.
 template <typename T, 
           typename I> // Models output iterator.
 void copy_region(const T& region, I out) {
   found_brick_vector const bricks(get_at_least_extant_bricks_in_region(region));
   for (typename found_brick_vector::const_iterator brick = bricks.begin();
        brick != bricks.end(); ++brick) {
// Copy if the key is in the region.
//     copy_if((*brick)->begin(), (*brick)->end(), out,
//	     compose1(point_in(region), 
     //                    select1st<internal_value_type>()));
     for (typename internal_type::iterator i = (*brick)->begin();
	 i != (*brick)->end(); ++i) {
       if (region_contains(region, i->first))
	 *out++ = *i;
     }
   }
 }

private:
  template <typename T>
  found_brick_vector const
  get_at_least_extant_bricks_in_region(const T& region) {
    Bounding_Box bounds = bounding_box(region);
    outer_key_type low  = get_map_key(bounds.low());
    outer_key_type high = get_map_key(bounds.high());
    found_brick_vector result;
    outer_key_type brick_key;
    BOOST_STATIC_ASSERT(NSD == 3); // Only 3-D for now.
    for (brick_key[2] = low[2]; brick_key[2] <= high[2]; ++brick_key[2]) {
      for (brick_key[1] = low[1]; brick_key[1] <= high[1]; ++brick_key[1]) {
        for (brick_key[0] = low[0]; brick_key[0] <= high[0]; ++brick_key[0]) {
          typename map_type::iterator const it = map_.find(brick_key);
          if (it != map_.end())
            result.push_back(&it->second);
        }
      }
    }
    return result;
  }
  // get_map_key(k) gets which brick k is in.
  const outer_key_type get_map_key(const key_type& k) const {
    outer_key_type result;
    for (unsigned int dim = 0; dim != NSD; ++dim)
      result[dim] = static_cast<outer_key_type::value_type>(k[dim] / brick_size[dim]);
    return result;
  }

  const key_type brick_size;
};


}} // End namespaces SCOREC::Util.

#endif // VOLUME_BRICK_MAP_H
