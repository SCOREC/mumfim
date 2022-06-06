#ifndef MUMFIM_TEST_MACROSCALE_TESTSUPPORT_H
#define MUMFIM_TEST_MACROSCALE_TESTSUPPORT_H
#include <apfDynamicMatrix.h>
#include <catch2/catch.hpp>
namespace test
{
  void compare_dynamic_matrices(const apf::DynamicMatrix & a,
                                       const apf::DynamicMatrix & b);
  void compare_dynamic_vectors(const apf::DynamicVector & a,
                                const apf::DynamicVector & b);
}
#endif  // MUMFIM_TEST_MACROSCALE_TESTSUPPORT_H
