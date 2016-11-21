#include <boost/test/minimal.hpp>
#include <vector>
#include <set>
#include "beta/volume_brick_map.h"
#include "mPoint.h"

using std::vector;
using std::set;
using std::cout;



int test_main(int, char*[]) {
  typedef volume_brick_map<mPoint,int> pointmap;
  pointmap points(mPoint(1.0, 1.0, 1.0));

  points.insert(std::make_pair(mPoint(0,0,0), 42));
  points.insert(std::make_pair(mPoint(0.5,0.3,0.3), 123));
  points.insert(std::make_pair(mPoint(3,3,4), 123));
  points.insert(std::make_pair(mPoint(3,3,3), 123));
  points.insert(std::make_pair(mPoint(-3,3,3), -123));

  size_t count = 0;
  for (pointmap::iterator i = points.begin();
       i != points.end(); ++i) {
    ++count;
    cout << "Number " << count << " is at " << i->first << ": " <<  i->second << "\n";
  }
  BOOST_CHECK(5 == count);

  cout << "Num bricks: " << points.num_bricks() << "\n";

  pointmap::iterator found = points.find(mPoint(3,3,3));
  BOOST_CHECK(found != points.end());
  BOOST_CHECK(found->second == 123);

  vector<pointmap::value_type> match;
  points.copy_region(SCOREC::Util::Sphere(mPoint(0,0,0), 1.0), back_inserter(match));
  cout << "Num matches: " << match.size() << "\n";
  
  return 0;
}

