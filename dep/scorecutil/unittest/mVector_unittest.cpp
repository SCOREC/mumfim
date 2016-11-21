#include <limits>
//#include <boost/test/test_tools.hpp>
//#include <boost/test/floating_point_comparison.hpp>
#include "mVector.h"

using std::numeric_limits;
using SCOREC::Util::mVector;

template <typename T1, typename T2, typename T3>
inline void assert_close(const T1& a, const T2& b, const T3 & e) {
  assert(a-b < e || b-a < e);
}

int main() {
  mVector e1(1.0, 0.0, 0.0);
  mVector e2(0.0, 1.0, 0.0);
  mVector e3(0.0, 0.0, 1.0);
  const double epsilon = numeric_limits<mVector::value_type>::epsilon();

  // Test orthogonal vectors.
  assert_close(e1*e2, 0.0, epsilon);
  assert_close(e1*e3, 0.0, epsilon);
  assert_close(e2*e3, 0.0, epsilon);

  assert_close(abs(e1), 1.0, epsilon);
  assert_close(abs(e2), 1.0, epsilon);
  assert_close(abs(e3), 1.0, epsilon);

  assert(data(e1) == data(e1));
  assert(data(e1) != data(e2));
  
  const mVector foo = 1.3 * e1 + 324.3 * e2 -12.2 * e3;
  const mVector bar(23, 43, 63);
  assert( abs(foo + bar) <= abs(foo) + abs(bar) ); // Triangle inequality.
  
  assert(abs(foo) > 1.0);

  
  double array[] = { 2.0, 4.1, 52.3 };
  mVector baz(array);
  mVector* vecArray = reinterpret_cast<mVector*>(array);
  assert(*vecArray == baz);
}
