#ifndef BIO_PACKED_DATA_H__
#define BIO_PACKED_DATA_H__
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_Pair.hpp>
#include <type_traits>
#include <vector>
template <typename DataType, typename OrdinalType, typename ExeSpace>
class PackedData
{
  public:
  using DataViewType = Kokkos::DualView<DataType *, ExeSpace>;
  using IndexViewType = Kokkos::DualView<OrdinalType *, ExeSpace>;
  using CIndexViewType = Kokkos::DualView<const OrdinalType *, ExeSpace>;
  using DeviceViewType = typename DataViewType::t_dev;
  // these types use lowerspace names for consistency with Kokkos
  using host_mirror_space = typename DataViewType::host_mirror_space;
  using memory_space = typename DataViewType::memory_space;
  // Functions below here should be private, but cannot have a private
  // host/device function
  template <typename Device>
  IndexViewType initializeRowIndex(IndexViewType row_sizes) const
  {
    row_sizes.template sync<Device>();
    IndexViewType row_index("RowIndex", row_sizes.extent(0) + 1);
    auto rows_d = row_sizes.template view<Device>();
    auto row_index_d = row_index.template view<Device>();
    Kokkos::parallel_scan(
        Kokkos::RangePolicy<typename Device::execution_space>(0, num_rows_),
        KOKKOS_LAMBDA(const OrdinalType i, OrdinalType & update,
                      const bool final) {
          update += rows_d(i);
          if (final)
          {
            row_index_d(i + 1) = update;
          }
        });
    // sync the row indices, so they are the same on the host and the device
    using device_memory_space = typename Device::memory_space;
    using other_memory_space = typename std::conditional<
        std::is_same<memory_space, device_memory_space>::value,
        host_mirror_space, memory_space>::type;
    row_index.template modify<device_memory_space>();
    row_index.template sync<other_memory_space>();
    return row_index;
  }
  // Public Functions
  PackedData(IndexViewType row_sizes) : num_rows_(row_sizes.extent(0))
  {
    row_index_ = initializeRowIndex<ExeSpace>(row_sizes);
    //row_index_ = initializeRowIndex<Kokkos::Serial>(row_sizes);
    auto row_index_h = row_index_.h_view;
    auto num_values =
        row_index_h.extent(0) > 0 ? row_index_h(row_index_h.extent(0) - 1) : 0;
    data_ = DataViewType("Data", num_values);
  }
  PackedData() : num_rows_(0) {}
  // deep_copy()
  //{
  //  // since row index is constant after initialization multiple packed
  //  // scalar data can share this without issue.
  //  //row_index_ = other.row_index_;
  //  // for multidimensional data we could use mirror_view with exe space here?
  //  //data_ = DataViewType("Data", other.data_.extent(0));
  //  //Kokkos::deep_copy(data_, other.data_);
  //}
  template <typename Device>
  KOKKOS_INLINE_FUNCTION Kokkos::View<DataType *, Device> getRow(
      OrdinalType idx) const
  {
    // only use the device version of the row index if we are using this function
    // on a device
#ifdef __CUDA_ARCH__
    auto row_index_d = row_index_.template view<Device>();
#else
    auto row_index_d = row_index_.template view<host_mirror_space>();
#endif
    // it is the responsibilty of the user to sync and mark modified the
    // device data using the sync and modify functions. Therefore, we do
    // not make any assumptions about if they modify the data here
    return Kokkos::subview(data_.template view<Device>(), Kokkos::make_pair(row_index_d(idx), row_index_d(idx + 1)));
    //return Kokkos::View<DataType *, Device>(data_.template view<Device>(), Kokkos::make_pair(row_index_d(idx), row_index_d(idx + 1)));
  }
  template <typename Device>
  void modify()
  {
    data_.template modify<Device>();
  }
  template <typename Device>
  void sync()
  {
    data_.template sync<Device>();
  }
  template <typename Device>
  auto getAllRows() -> Kokkos::View<typename DataViewType::traits::data_type,
                                    typename DataViewType::traits::array_layout,
                                    Device>
  {
    return data_.template view<Device>();
  }
  KOKKOS_INLINE_FUNCTION
  OrdinalType getNumRows() const { return num_rows_; }
  private:
  DataViewType data_;
  CIndexViewType row_index_;
  OrdinalType num_rows_;
};
#endif
