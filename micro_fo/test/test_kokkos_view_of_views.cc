#include <vector>
#include "bioViewOfView.h"
#include "catch2/catch.hpp"
#include <Kokkos_Core.hpp>

KOKKOS_INLINE_FUNCTION
double getHostValue(int idx) { return 10.0; }
KOKKOS_INLINE_FUNCTION
double getDeviceValue(int idx) { return idx % 2 ? 20.5 : getHostValue(idx); }
TEST_CASE("view of views", "[kokkos]")
{
  // 1. verify that the host data is th same as the input data
  //    this should also check to make sure that the sizes are correct
  // v.syncToDevice()
  // 2. verify that data round trips properly
  //    a) set every other value on device
  //    b) send modified data back to host and check it
  // v.syncToHost()
  int numOuterRows = 10;
  std::vector<std::vector<double>> data;
  data.reserve(numOuterRows);
  for (int i = 0; i < numOuterRows; ++i)
  {
    int numInnerRows = i*2;
    std::vector<double> innerData;
    innerData.reserve(numInnerRows);
    for (int j = 0; j < numInnerRows; ++j)
    {
      innerData.push_back(getHostValue(j));
    }
    data.push_back(innerData);
  }
  bio::ViewOfView<double> vofv(data);
  auto hh_view = vofv.h_view();
  auto dd_view = vofv.d_view();
  REQUIRE(hh_view.extent(0) == data.size());
  REQUIRE(dd_view.extent(0) == data.size());
  for (size_t i = 0; i < hh_view.extent(0); ++i)
  {
    // check that the inner extents are the correct sizes
    REQUIRE(hh_view(i).extent(0) == data[i].size());
    for (size_t j = 0; j < hh_view(i).extent(0); ++j)
    {
      // check that the data is correct on the host side
      REQUIRE(hh_view(i)(j) == Approx(getHostValue(j)));
    }
  }
  Kokkos::parallel_for(hh_view.extent(0), KOKKOS_LAMBDA(int i) {
    for (size_t j = 0; j < dd_view(i).extent(0); ++j)
    {
      // we include this code here becuase this will verify that
      // the data actually round trips to the device and back to the host
      if(j%2)
        dd_view(i)(j) = getDeviceValue(j);
    }
  });
  Kokkos::fence();
  vofv.syncToHost();

  for (size_t i = 0; i < hh_view.extent(0); ++i)
  {
    for (size_t j = 0; j < hh_view(i).extent(0); ++j)
    {
      REQUIRE(hh_view(i)(j) == Approx(getDeviceValue(j)));
    }
  }
}

