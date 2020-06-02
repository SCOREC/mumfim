#ifndef BIO_BATCHED_REORDER_MESH_H__
#define BIO_BATCHED_REORDER_MESH_H__
#include <Kokkos_Core.hpp>
namespace bio
{
  /**
   * Reorder the mesh for analysis. This places all of the fixed vertices at the
   * end of the vertex array. This has the effect of allowing us to loop over only free
   * and only fixed portions of the dofs. This function assumes that the first dof
   * is zero, and that there are no gaps in the dofs
   */
  template <typename Scalar, typename Ordinal, typename ExeSpace>
  struct ReorderMesh
  {
    class TagInitPermutation{};
    class TagScanPermutation{};
    class TagPermuteConnectivity{};
    class TagPermuteCoordinates{};
    using RWOV = Kokkos::View<Ordinal *, ExeSpace>;
    using RWSV = Kokkos::View<Scalar *, ExeSpace>;
    Ordinal nvert_;
    RWOV fixed_vert_;
    RWOV permutation_array_;
    RWOV connectivity_;
    RWSV coordinates_;
    RWSV tmp_;
    ReorderMesh(Ordinal nvert, RWOV fixed_vert) : nvert_(nvert), fixed_vert_(fixed_vert), permutation_array_(RWOV("permutation array", nvert_))
    {
    }
    KOKKOS_INLINE_FUNCTION
    void operator()(TagInitPermutation, const int i) const
    {
        permutation_array_(fixed_vert_(i)) = i+1;
    }
    KOKKOS_INLINE_FUNCTION
    void operator()(TagScanPermutation, const int i, Scalar& update, const bool final) const
    {
      Ordinal updateVal = permutation_array_(i);
      bool isFixedVert = updateVal > 0;
      update+=isFixedVert?1:0;
      if(final) 
      {
        Ordinal nfixed = fixed_vert_.extent(0);
        permutation_array_(i) = isFixedVert?(nvert_-nfixed)+(updateVal-1):
          i-update;
      }
    }
    KOKKOS_INLINE_FUNCTION
    void operator()(TagPermuteConnectivity, const int i) const
    {
      connectivity_(i) = permutation_array_(connectivity_(i));
    }
    KOKKOS_INLINE_FUNCTION
    void operator()(TagPermuteCoordinates, const int i) const
    {
      coordinates_(3*permutation_array_(i/3)+i%3) = tmp_(i);
    }
    // create the map from the old to the new verts
    void createPermutationArray() const
    {
      Kokkos::parallel_for(Kokkos::RangePolicy<TagInitPermutation,ExeSpace>(0,fixed_vert_.extent(0)), *this);
      Kokkos::parallel_scan(Kokkos::RangePolicy<TagScanPermutation,ExeSpace>(0,permutation_array_.extent(0)), *this);

    }
    void applyPermutationToConnectivity(RWOV connectivity)
    {
      connectivity_ = connectivity;
      Kokkos::parallel_for(Kokkos::RangePolicy<TagPermuteConnectivity,ExeSpace>(0,connectivity_.extent(0)),*this);
    }
    void applyPermutationToCoordinates(RWSV coordinates) 
    {
      coordinates_ = coordinates;
      tmp_ = RWSV("coordinates", coordinates_.extent(0));
      Kokkos::deep_copy(tmp_, coordinates_);
      Kokkos::parallel_for(Kokkos::RangePolicy<TagPermuteCoordinates,ExeSpace>(0,coordinates_.extent(0)), *this);
      // unallocate the temporary view
      tmp_ = RWSV();
    }
    //void applyPermuationToDof();
  };
}
#endif
