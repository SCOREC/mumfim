#include "TestSupport.h"
namespace test
{
  void compare_dynamic_matrices(const apf::DynamicMatrix & a,
                                const apf::DynamicMatrix & b)
  {
    REQUIRE(a.getRows() == b.getRows());
    REQUIRE(a.getColumns() == b.getColumns());
    for (size_t i = 0; i < a.getRows(); ++i)
    {
      for (size_t j = 0; j < a.getColumns(); ++j)
      {
        REQUIRE(a(i, j) == Approx(b(i, j)));
      }
    }
  }
  void compare_dynamic_vectors(const apf::DynamicVector & a,
                               const apf::DynamicVector & b)
  {
    REQUIRE(a.size() == b.size());
    for(size_t i = 0; i<a.size(); ++i) {
      REQUIRE(a(i) == Approx(b(i)));
    }
  }
}  // namespace test
