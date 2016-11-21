#ifndef NESTED_ITERATOR_H
#define NESTED_ITERATOR_H

// This iterator allows you to iterate over a nested container.
// For example:
// vector<vector<int> > vvi(3, vector<int>(2, 42)); // Make a 3x2 vector of ints, all "42".
// nested_iterator<vector<vector<int> > > vbeg(vvi.begin(), vvi.end());
// nested_iterator<vector<vector<int> > > vend(vvi.begin(), vvi.end());
// std::copy(vbeg, vend, ostream_iterator<int>(cout, " ")); // Prints 42 six times.
#include <map>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/type_traits/remove_reference.hpp>

  {
namespace Util {

// We use a traits class here because some containers's value_type
// isn't const when it seems like it should.
// e.g., std::set's iterator's value_type is const. 
/*template <typename T>
struct nested_iterator_traits {
  typedef T container_type;
  typedef typename container_type::iterator::value_type
                                        inner_container_type;
  typedef typename container_type::iterator         outer_iterator;
  typedef typename inner_container_type::iterator   inner_iterator;
  typedef typename inner_container_type::iterator::value_type value_type;
};

template <typename T, typename C, typename A>
struct nested_iterator_traits<std::set<T,C,A> > {
  typedef std::set<T,C,A> container_type;
  typedef typename container_type::iterator         outer_iterator;
  typedef typename outer_iterator::value_type
                                             inner_container_type;
  typedef typename inner_container_type::const_iterator   inner_iterator;
  typedef const typename inner_iterator::value_type value_type;
  };*/


// iterator_qualified_value<I>::type has a const if the iterator is const.
template <typename Iter>
struct iterator_qualified_value {
  typedef typename ::boost::iterator_reference<Iter>::type reference;
  typedef typename ::boost::remove_reference<reference>::type type;
};


template <typename C>
struct container_iterator {
  typedef typename C::iterator type;
};

template <typename C>
struct container_iterator<const C> {
  typedef typename C::const_iterator type;
};


template <typename C>
struct container_iterator<C&> {
  typedef typename container_iterator<typename boost::remove_reference<C>::type>::type type;
};


template <typename OuterIter, typename InnerIter>
struct nested_iterator_traits {
  typedef typename boost::iterator_reference<OuterIter>::type outer_reference_type;
  typedef InnerIter     iterator;
  iterator begin(outer_reference_type c) const { return c.begin(); }
  iterator end  (outer_reference_type c) const { return c.end();   }
};

  /*template <typename K, typename T, typename C, typename A>
struct nested_iterator_traits<typename std::template map<K,T,C,A>::iterator> {
  typedef typename std::template map<K,T,C,A>::iterator OuterIter;
  typedef typename std::iterator_traits<OuterIter>::value_type outer_value_type;
  typedef typename outer_value_type::const_iterator inner_iterator;
  static inner_iterator begin(outer_value_type& c) { return c.second.begin(); }
  static inner_iterator end  (outer_value_type& c) { return c.second.end();   }
};
  */
  /*template <typename T, typename C, typename A>
struct nested_iterator_traits<typename std::template set<T,C,A>::iterator> {
  typedef typename std::template set<T,C,A>::iterator OuterIter;
  typedef const typename std::iterator_traits<OuterIter>::value_type outer_value_type;
  typedef typename outer_value_type::const_iterator inner_iterator;
  static inner_iterator begin(outer_value_type& c) { return c.begin(); }
  static inner_iterator end  (outer_value_type& c) { return c.end();   }
  };*/


template<typename T, typename Iter>
struct second_range {
  typedef Iter iterator;
  Iter begin(T& c) const { return c.second.begin(); }
  Iter end  (T& c) const { return c.second.end();   }
};



template <typename OuterIter, 
  typename InnerRangeGetter = 
  nested_iterator_traits<OuterIter, 
    typename container_iterator<typename iterator_qualified_value<OuterIter>::type>::type>
         >
class nested_iterator {
public:
  typedef typename InnerRangeGetter::iterator inner_iterator;
  typedef OuterIter outer_iterator;
  typedef typename boost::iterator_value<inner_iterator>::type           value_type;
  typedef typename boost::iterator_reference<inner_iterator>::type         reference;
  typedef value_type*           pointer;
  typedef typename boost::iterator_category<inner_iterator>::type         iterator_category;
  typedef std::size_t         difference_type;

  nested_iterator(outer_iterator oib, 
		  outer_iterator oie,
		  InnerRangeGetter getter = InnerRangeGetter())
    : oi_(oib), oiend_(oie), getter_(getter) {
    if (oi_ != oiend_) { // Only initialize the inner if we can.
      ii_    = getter_.begin(*oi_);
      iiend_ = getter_.end(*oi_);
      // Skip over empty inners.
      while (ii_ == iiend_) {
	// Increment the outer.
	if (oi_ != oiend_) { // Off the end of the inner but not the outer... 
	  ++oi_; // Next outer...
	  if (oi_ != oiend_) {
	    // If we aren't at the end, reset the inner...
	    ii_    = getter_.begin(*oi_);
	    iiend_ = getter_.end(*oi_);
	  }
	} else {
	  // Off the end of the outer too.
	  assert(oi_ == oiend_);
	  return; 
	}
      }
    }
  }
  // Initialize with a particular inner iterator.
  nested_iterator(outer_iterator oib, outer_iterator oie,
		  inner_iterator it)
    : oi_(oib), oiend_(oie) {
    assert(oi_ != oiend_);
    ii_    = it;
    iiend_ = getter_.end(*oi_);
  }
  nested_iterator& operator++() {
    ++ii_; // Increment inner.
    // Skip empty inners...
    while (ii_ == iiend_ && 
	   oi_ != oiend_) { // Off the end of the inner but not the outer... 
      ++oi_; // Next outer...
      if (oi_ != oiend_) {
	// If we aren't at the end, reset the inner...
        ii_    = getter_.begin(*oi_);
	iiend_ = getter_.end(*oi_);
      }
    }
    return *this;
  }
  reference operator*() const { return *ii_; }
  pointer operator->() const { return &*ii_; }

  friend
  bool operator==(const nested_iterator& a,
                  const nested_iterator& b) {
    if (a.ii_ == b.ii_) return true; // Same place
    if (a.oi_ == a.oiend_ &&
        b.oi_ == b.oiend_) return true; // It's at the end.
    return false; // Different.
  }
  friend
  bool operator!=(const nested_iterator& a,
                  const nested_iterator& b) {
    return !(a == b);
  }
 
  // ...
private:
  outer_iterator oi_, oiend_;
  inner_iterator ii_, iiend_;
  InnerRangeGetter getter_;
};


}} // End namespaces SCOREC::Util.

#endif // NESTED_ITERATOR_H
