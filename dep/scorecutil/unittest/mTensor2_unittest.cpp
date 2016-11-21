#include <limits>
//#include <boost/test/test_tools.hpp>
//#include <boost/test/floating_point_comparison.hpp>
#include "mTensor2.h"

using std::numeric_limits;

template <typename T1, typename T2, typename T3>
inline void assert_close(const T1& a, const T2& b, const T3 & e) {
  assert(a-b < e || b-a < e);
}

int main() {
  //  mVector v(1.0, 0.0, 0.0);
  mTensor2 t2_1(mVector(1,0,0),
		mVector(0,1,0),
		mVector(0,0,1));
  assert(t2_1 == transpose(t2_1));
  assert(trace(t2_1) == 3.0);
  assert(det(t2_1) == 1.0);

  
}
