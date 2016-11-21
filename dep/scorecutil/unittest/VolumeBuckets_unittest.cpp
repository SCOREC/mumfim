#include <limits>
#include <vector>
#include <iterator>
#include <iostream>
#include <boost/type_traits.hpp>
#include "VolumeBuckets.h"

using std::numeric_limits;
using std::vector;
using std::ostream_iterator;
using std::cout;
using std::endl;
using boost::is_same;

template <typename T1, typename T2, typename T3>
inline void assert_close(const T1& a, const T2& b, const T3 & e) {
  assert(a-b < e || b-a < e);
}

int main() {
  vector<int> ints(3);
  ints[0] = 0;   ints[1] = 1;   ints[2] = 2;
  vector<vector<int> > vii(4, ints);
  NestedIterator<vector<vector<int> >::iterator> iterbeg(vii.begin(), vii.end());
  NestedIterator<vector<vector<int> >::iterator> iterend(vii.end(),   vii.end());
  BOOST_STATIC_ASSERT(( is_same<int, vector<int>::value_type>::value ));
  BOOST_STATIC_ASSERT(( is_same<int, vector<vector<int> >::value_type::value_type>::value ));
  BOOST_STATIC_ASSERT(( is_same<int, NestedIterator<vector<vector<int> >::iterator>::value_type>::value ));


  std::copy(iterbeg, iterend, ostream_iterator<int>(std::cout, " "));
  VolumeBuckets<mVector> theBuckets(mVector(0,0,0), mVector(2,2,2));


  ////////////////////
  mVector v(1,1,1);
  theBuckets.insert(v);
  vector<mVector*> thePoints = 
    theBuckets.get_all_objects_with_points_within_r_of_pos(1.0, mVector(0.8, 1.2, 0.9));
  assert(thePoints.size() == 1);
  theBuckets.erase(thePoints[0]);
  thePoints = 
    theBuckets.get_all_objects_with_points_within_r_of_pos(1.0, mVector(0.8, 1.2, 0.9));
  assert(thePoints.size() == 0);
}
