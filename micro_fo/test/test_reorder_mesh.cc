#include <array>
#include "catch2/catch.hpp"
#include <Kokkos_Core.hpp>
#include <bioBatchedReorderMesh.h>
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
  bio::ReorderMesh<Scalar,Ordinal,ExeSpace> reorder(nverts, fixed_verts_d);
  reorder.createPermutationArray();
  SECTION("Check connectivity")
  {
    Kokkos::View<Ordinal*, ExeSpace> connectivity_d("connectivity", 11);
    auto connectivity_h = Kokkos::create_mirror_view(connectivity_d);
    std::array<Ordinal,11> initial_connectivity {1,7,5,3,3,1,2,4,6,0,1};
    std::array<Ordinal,11> final_connectivity {0,6,3,5,5,0,1,2,4,7,0};
    for(int i=0; i<11; ++i)
    {
      connectivity_h(i) = initial_connectivity[i];
    }
    Kokkos::deep_copy(connectivity_d, connectivity_h);
    reorder.applyPermutationToConnectivity(connectivity_d);
    Kokkos::deep_copy(connectivity_h, connectivity_d);
    for(std::size_t i=0; i<connectivity_h.extent(0);++i)
    {
      REQUIRE(connectivity_h(i) == final_connectivity[i]);
    }
  }
  SECTION("Check Coordinates") {
    Kokkos::View<Scalar*, ExeSpace> coordinates_d("coordinates", nverts*3);
    std::array<Scalar,24> final_coordinates = {3,4,5,6,7,8,12,13,14,15,16,17,18,19,20,9,10,11,21,22,23,0,1,2};
    auto coordinates_h = Kokkos::create_mirror_view(coordinates_d);
    for(std::size_t i=0; i<coordinates_h.extent(0); ++i)
    {
      coordinates_h(i) = i;
    }
    Kokkos::deep_copy(coordinates_d, coordinates_h);
    reorder.applyPermutationToCoordinates(coordinates_d);
    Kokkos::deep_copy(coordinates_h, coordinates_d);
    for (std::size_t i=0; i<coordinates_h.extent(0); ++i)
    {
      REQUIRE(final_coordinates[i] == Approx(coordinates_h(i)));
    }
  }
}
