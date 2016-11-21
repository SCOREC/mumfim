#include <boost/test/minimal.hpp>
#include <boost/type_traits.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <vector>
#include <set>
#include "beta/nested_iterator.h"

using std::vector;
using std::set;


struct Tester {
  int foo() const { return 1234; }
};

int test_main(int, char*[]) {
  vector<vector<int> > vvi(3, vector<int>(2, 42)); // Make a 3x2 vector of ints, all "42".
  nested_iterator<vector<vector<int> >::iterator> vbeg(vvi.begin(), vvi.end());
  nested_iterator<vector<vector<int> >::iterator> vend(vvi.end(), vvi.end());
  size_t count = 0;
  for (; vbeg != vend; ++vbeg) {
    std::cout << *vbeg << "\n";
    BOOST_CHECK((*vbeg == 42));
    ++count;
  }
  BOOST_CHECK(6 == count);

  using namespace boost;
  BOOST_STATIC_ASSERT((is_same<iterator_qualified_value<vector<double>::const_iterator>::type,
		       const double>::value));
  BOOST_STATIC_ASSERT((is_same<container_iterator<const vector<double> >::type,
		       vector<double>::const_iterator>::value));
  set<vector<double> > setvd;
  BOOST_STATIC_ASSERT((is_same<std::iterator_traits<set<vector<double> >::iterator>::value_type,
		       vector<double> >::value));

  setvd.insert(vector<double>(0, -1.0)); // An empty vector.
  setvd.insert(vector<double>(4, 1.0));  // Four 1.0s
  setvd.insert(vector<double>(0, -1.0)); // Another empty one.
  setvd.insert(vector<double>(3,5.0));   // Three 5.0s.

  nested_iterator<set<vector<double> >::iterator> 
    setvdbeg(setvd.begin(), setvd.end());
  nested_iterator<set<vector<double> >::iterator> 
    setvdend(setvd.end(),   setvd.end());

  count = 0;

  for (; setvdbeg != setvdend; ++setvdbeg) {
    if (count < 4) {
      BOOST_CHECK(*setvdbeg == 1.0);
    } if (count >= 4) {
      BOOST_CHECK(*setvdbeg == 5.0);
    }
    BOOST_CHECK(*setvdbeg == 1.0 || *setvdbeg == 5.0); // It's one of these.
    ++count;
    BOOST_CHECK(count <= 7);
  }
  BOOST_CHECK(count == 7);


  // Try case of zero-sized outer.
  setvd.clear();
  setvdbeg = nested_iterator<set<vector<double> >::iterator>(setvd.begin(), setvd.end());
  setvdend = nested_iterator<set<vector<double> >::iterator>(setvd.end(),   setvd.end());
  count = 0;
  for (; setvdbeg != setvdend; ++setvdbeg) {}
  BOOST_CHECK(count == 0);


  vector<vector<Tester> > test(2, vector<Tester>(2));
  nested_iterator<vector<vector<Tester> >::iterator> tbegin(test.begin(), test.end());
  nested_iterator<vector<vector<Tester> >::iterator> tend(test.end(), test.end());
  for (; tbegin != tend;
       ++tbegin) {
    BOOST_CHECK(1234 == tbegin->foo());
  }

  // Try case of doubly-nested.

  return 0;
}

