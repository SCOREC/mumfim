#ifndef MUMFIM_BATCHED_REORDER_MESH_H
#define MUMFIM_BATCHED_REORDER_MESH_H
#include <Kokkos_Core.hpp>
namespace mumfim
{
  /**
   * Reorder the mesh for analysis. This places all of the fixed vertices at the
   * end of the vertex array. This has the effect of allowing us to loop over
   * only free and only fixed portions of the dofs. This function assumes that
   * the first dof is zero, and that there are no gaps in the dofs
   */
  template <typename Scalar, typename Ordinal, typename Layout, typename ExeSpace>
  struct ReorderMesh
  {
    struct TagInitPermutation
    {
    };
    struct TagScanPermutation
    {
    };
    struct TagPermuteConnectivity
    {
    };
    struct TagPermuteCoordinates
    {
    };
    struct TagPermuteNodal
    {
    };
    struct TagPermuteFixedVert
    {
    };
    using OrdinalVertViewType = Kokkos::View<Ordinal *, Layout, ExeSpace>;
    using ScalarVertViewType = Kokkos::View<Scalar *, Layout, ExeSpace>;
    using ScalarDofViewType = Kokkos::View<Scalar *[3], Layout, ExeSpace>;
    using ConnectivityType = Kokkos::View<Ordinal*[2], Layout, ExeSpace>;

    Ordinal nvert_;
    OrdinalVertViewType fixed_vert_;
    OrdinalVertViewType permutation_array_;
    ConnectivityType connectivity_;
    ScalarDofViewType coordinates_;
    ScalarVertViewType nodal_data_;
    ScalarVertViewType tmp_nodal_;
    ScalarDofViewType tmp_dof_;

    ReorderMesh(Ordinal nvert, OrdinalVertViewType fixed_vert)
        : nvert_(nvert)
        , fixed_vert_(fixed_vert)
        , permutation_array_(OrdinalVertViewType("permutation array", nvert_))
    {
    }
    KOKKOS_INLINE_FUNCTION
    void operator()(TagInitPermutation, const int i) const
    {
      permutation_array_(fixed_vert_(i)) = i + 1;
    }
    KOKKOS_INLINE_FUNCTION
    void operator()(TagScanPermutation,
                    const int i,
                    Scalar & update,
                    const bool final) const
    {
      Ordinal updateVal = permutation_array_(i);
      bool isFixedVert = updateVal > 0;
      update += isFixedVert ? 1 : 0;
      if (final)
      {
        Ordinal nfixed = fixed_vert_.extent(0);
        permutation_array_(i) =
            isFixedVert ? (nvert_ - nfixed) + (updateVal - 1) : i - update;
      }
    }
    KOKKOS_INLINE_FUNCTION
    void operator()(TagPermuteConnectivity, const int i) const
    {
      for(int j=0; j<2; ++j)
      {
        connectivity_(i,j) = permutation_array_(connectivity_(i,j));
      }
    }
    KOKKOS_INLINE_FUNCTION
    void operator()(TagPermuteCoordinates, const int i) const
    {
      //coordinates_(3 * permutation_array_(i / 3) + i % 3) = tmp_dof_(i);
      for(int j=0; j<3; ++j)
      {
        coordinates_(permutation_array_(i),j) = tmp_dof_(i,j);
      }
    }
    KOKKOS_INLINE_FUNCTION
    void operator()(TagPermuteNodal, const int i) const
    {
      nodal_data_(permutation_array_(i)) = tmp_nodal_(i);
    }
    KOKKOS_INLINE_FUNCTION
    void operator()(TagPermuteFixedVert, const int i) const
    {
      fixed_vert_(i) = permutation_array_(fixed_vert_(i));
    }
    // create the map from the old to the new verts
    void createPermutationArray() const
    {
      Kokkos::parallel_for(Kokkos::RangePolicy<TagInitPermutation, ExeSpace>(
                               0, fixed_vert_.extent(0)),
                           *this);
      Kokkos::parallel_scan(Kokkos::RangePolicy<TagScanPermutation, ExeSpace>(
                                0, permutation_array_.extent(0)),
                            *this);
    }
    void applyPermutationToConnectivity(ConnectivityType connectivity)
    {
      connectivity_ = connectivity;
      Kokkos::parallel_for(
          Kokkos::RangePolicy<TagPermuteConnectivity, ExeSpace>(
              0, connectivity_.extent(0)),
          *this);
      connectivity_ = ConnectivityType();
    }
    void applyPermutationToCoordinates(ScalarDofViewType coordinates)
    {
      coordinates_ = coordinates;
      tmp_dof_ = ScalarDofViewType("coordinates", coordinates_.extent(0));
      Kokkos::deep_copy(tmp_dof_, coordinates_);
      Kokkos::parallel_for(Kokkos::RangePolicy<TagPermuteCoordinates, ExeSpace>(
                               0, coordinates_.extent(0)),
                           *this);
      // unallocate the temporary view
      tmp_dof_ = ScalarDofViewType();
      coordinates_ = ScalarDofViewType();
    }
    void applyPermutationToNodalValue(ScalarVertViewType nodal)
    {
      nodal_data_ = nodal;
      tmp_nodal_ = ScalarVertViewType("nodal data", nodal_data_.extent(0));
      Kokkos::deep_copy(tmp_nodal_, nodal_data_);
      Kokkos::parallel_for(Kokkos::RangePolicy<TagPermuteNodal, ExeSpace>(
                               0, nodal_data_.extent(0)),
                           *this);
      // unallocate the temporary view
      tmp_nodal_ = ScalarVertViewType();
      nodal_data_ = ScalarVertViewType();
    }
    void applyPermutationToFixedVert()
    {
      Kokkos::parallel_for(Kokkos::RangePolicy<TagPermuteFixedVert, ExeSpace>(
                               0, fixed_vert_.extent(0)),
                           *this);
    }
    ScalarVertViewType getPermutationArray() { return permutation_array_; }
  };
}  // namespace mumfim
#endif
