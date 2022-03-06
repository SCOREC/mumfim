#include <array>
#include "catch2/catch.hpp"
#include <Kokkos_Core.hpp>
#include <mumfim/microscale/batched_analysis/BatchedReorderMesh.h>
using Ordinal=int;
using Scalar=double;
using ExeSpace=Kokkos::DefaultExecutionSpace;

TEST_CASE("reorder mesh", "[mesh]")
{
  Kokkos::View<Ordinal*, ExeSpace> fixed_verts_d("fixed_verts", 3);
  auto fixed_verts_h = Kokkos::create_mirror_view(fixed_verts_d);
  std::array<Ordinal,3> fixed_verts {3,7,0};
  int nverts = 8;
  for(int i=0; i<3; ++i)
  {
    fixed_verts_h(i) = fixed_verts[i];
  }
  Kokkos::deep_copy(fixed_verts_d, fixed_verts_h);
  mumfim::ReorderMesh<Scalar, Ordinal, ExeSpace::array_layout, ExeSpace> reorder(
      nverts, fixed_verts_d);
  reorder.createPermutationArray();
  SECTION("Check connectivity")
  {
    Kokkos::View<Ordinal * [2], ExeSpace::array_layout, ExeSpace>
        connectivity_d("connectivity", 10);
    auto connectivity_h = Kokkos::create_mirror_view(connectivity_d);
    std::array<Ordinal, 10> initial_connectivity{1, 7, 5, 3, 3, 1, 2, 4, 6, 0};
    std::array<Ordinal, 10> final_connectivity{0, 6, 3, 5, 5, 0, 1, 2, 4, 7};
    for (int i = 0; i < 10; ++i)
    {
      auto col = i % 2;
      auto row = i / 2;
      connectivity_h(row, col) = initial_connectivity[i];
    }
    Kokkos::deep_copy(connectivity_d, connectivity_h);
    reorder.applyPermutationToConnectivity(connectivity_d);
    Kokkos::deep_copy(connectivity_h, connectivity_d);
    for (std::size_t i = 0; i < connectivity_h.extent(0); ++i)
    {
      auto col = i % 2;
      auto row = i / 2;
      REQUIRE(connectivity_h(row, col) == final_connectivity[i]);
    }
    Kokkos::View<Scalar * [3], ExeSpace::array_layout, ExeSpace> coordinates_d(
        "coordinates", nverts);
    std::array<Scalar, 24> final_coordinates = {3,  4,  5,  6,  7,  8,  12, 13,
                                                14, 15, 16, 17, 18, 19, 20, 9,
                                                10, 11, 21, 22, 23, 0,  1,  2};
    auto coordinates_h = Kokkos::create_mirror_view(coordinates_d);
    for (std::size_t i = 0; i < 24; ++i)
    {
      auto row = i / 3;
      auto col = i % 3;
      coordinates_h(row, col) = i;
    }
    Kokkos::deep_copy(coordinates_d, coordinates_h);
    reorder.applyPermutationToCoordinates(coordinates_d);
    Kokkos::deep_copy(coordinates_h, coordinates_d);
    for (std::size_t i=0; i<coordinates_h.extent(0); ++i)
    {
      auto row = i / 3;
      auto col = i % 3;
      REQUIRE(final_coordinates[i] == Approx(coordinates_h(row, col)));
    }
  }
}
