
#include <boost/test/minimal.hpp>
#include <vector>
#include "mPoint.h"
#include "affine_space.h"



int test_main(int argc, char** argv) {
  using std::vector;
  using std::cout;
  using std::endl;
  const double epsilon = 1.0e-15;

  BOOST_CHECK(affine_combine(1.0, 2.0, 1.0) == 1.0);
  BOOST_CHECK(affine_combine(1.0, 2.0, 0.0) == 2.0);
  BOOST_CHECK(affine_combine(1.0, 2.0, 0.5) == 1.5);

  vector<mPoint> points(3, mPoint(1.0, 1.0, 1.0));
  BOOST_CHECK(abs(mean(points.begin(), points.end()) - points[0]) < epsilon);
  points[1] = mPoint(-2.0, -2.0, -2.0);
  BOOST_CHECK(abs(mean(points.begin(), points.end()) - mPoint(0.0, 0.0, 0.0)) < epsilon);

  vector<double> weights;
  weights.push_back(1.0);
  weights.push_back(2.0);
  weights.push_back(3.0);
  
  range_reader<vector<double>::iterator> range_reader_tester(weights.begin(), weights.end());
  BOOST_CHECK(range_reader_tester(0) == 1.0);
  BOOST_CHECK(range_reader_tester() == 2.0);
  BOOST_CHECK(range_reader_tester(0) == 3.0);
  
  const double weighted_average = weighted_mean(weights.begin(), weights.end(),
						weights.begin(), weights.end());
  cout << "Weighted average: " << weighted_average << endl;
  BOOST_CHECK((7.0/3.0) + 1.0e-10 > weighted_average);
  BOOST_CHECK((7.0/3.0) - 1.0e-10 < weighted_average);

  return 0;
}
