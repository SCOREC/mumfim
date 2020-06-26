#ifndef BIO_PACKED_DATA_H__
#define BIO_PACKED_DATA_H__
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_Pair.hpp>
#include <type_traits>

// get a subview with all of the data kept the same except for the first index which gets a Kokkos::pair to
// define its bounds
template<typename VT, typename PairType,
  typename std::enable_if<VT::rank==1,int>::type = 0
  >
KOKKOS_FORCEINLINE_FUNCTION
auto getSubview(VT view, PairType pair) -> 
  Kokkos::View<typename VT::data_type, typename VT::array_layout, typename VT::device_type>
{
  return Kokkos::View<typename VT::data_type, typename VT::array_layout, typename VT::device_type>(view, pair);
}

template<typename VT, typename PairType,
  typename std::enable_if<VT::rank==2,int>::type = 0
  >
KOKKOS_FORCEINLINE_FUNCTION
auto getSubview(VT view, PairType pair) -> 
  Kokkos::View<typename VT::data_type, typename VT::array_layout, typename VT::device_type>
{
  return Kokkos::View<typename VT::data_type, typename VT::array_layout, typename VT::device_type>(view, pair, Kokkos::ALL());
}
template<typename VT, typename PairType,
  typename std::enable_if<VT::rank==3,int>::type = 0
  >
KOKKOS_FORCEINLINE_FUNCTION
auto getSubview(VT view, PairType pair) -> 
  Kokkos::View<typename VT::data_type, typename VT::array_layout, typename VT::device_type>
{
  return Kokkos::View<typename VT::data_type, typename VT::array_layout, typename VT::device_type>(view, pair, Kokkos::ALL(),Kokkos::ALL());
}
template<typename VT, typename PairType,
  typename std::enable_if<VT::rank==4,int>::type = 0
  >
KOKKOS_FORCEINLINE_FUNCTION
auto getSubview(VT view, PairType pair) -> 
  Kokkos::View<typename VT::data_type, typename VT::array_layout, typename VT::device_type>
{
  return Kokkos::View<typename VT::data_type, typename VT::array_layout, typename VT::device_type>(view, pair, Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
}
template<typename VT, typename PairType,
  typename std::enable_if<VT::rank==5,int>::type = 0
  >
KOKKOS_FORCEINLINE_FUNCTION
auto getSubview(VT view, PairType pair) -> 
  Kokkos::View<typename VT::data_type, typename VT::array_layout, typename VT::device_type>
{
  return Kokkos::View<typename VT::data_type, typename VT::array_layout, typename VT::device_type>(view, pair, Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(), Kokkos::ALL());
}
template<typename VT, typename PairType,
  typename std::enable_if<VT::rank==6,int>::type = 0
  >
KOKKOS_FORCEINLINE_FUNCTION
auto getSubview(VT view, PairType pair) -> 
  Kokkos::View<typename VT::data_type, typename VT::array_layout, typename VT::device_type>
{
  return Kokkos::View<typename VT::data_type, typename VT::array_layout, typename VT::device_type>(view, pair, Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
}
template<typename VT, typename PairType,
  typename std::enable_if<VT::rank==7,int>::type = 0
  >
KOKKOS_FORCEINLINE_FUNCTION
auto getSubview(VT view, PairType pair) -> 
  Kokkos::View<typename VT::data_type, typename VT::array_layout, typename VT::device_type>
{
  return Kokkos::View<typename VT::data_type, typename VT::array_layout, typename VT::device_type>(view, pair, Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
}




template <typename DataType, typename Arg1Type=void, typename Arg2Type=void, typename Arg3Type = void>
class PackedData : public  Kokkos::ViewTraits<DataType, Arg1Type, Arg2Type, Arg3Type>
{
  public:

  using traits = Kokkos::ViewTraits<DataType,Arg1Type, Arg2Type, Arg3Type>;
  using OrdinalType = typename traits::size_type;
  using DataViewType = Kokkos::DualView<DataType, Arg1Type, Arg2Type, Arg3Type>;
  using IndexViewType = Kokkos::DualView<OrdinalType *, Arg1Type, Arg2Type, Arg3Type>;
  using CIndexViewType = Kokkos::DualView<const OrdinalType* , Arg1Type, Arg2Type, Arg3Type>;
  using DeviceViewType = typename DataViewType::t_dev;
  using host_mirror_space = typename DataViewType::host_mirror_space;
  using memory_space = typename DataViewType::memory_space;
  using ExeSpace = typename traits::execution_space;

  static_assert(traits::rank_dynamic >= 1, "The PackedData type expects a dynamic rank for the numer of rows");
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
  PackedData() = default;
  template <typename Device>
  KOKKOS_INLINE_FUNCTION
  auto getRow(
      OrdinalType idx) const -> Kokkos::View<DataType,
                                    typename traits::array_layout,
                                    Device>
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
    return getSubview(data_.template view<Device>(),Kokkos::make_pair(row_index_d(idx), row_index_d(idx + 1)));
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
  auto getAllRows() const -> Kokkos::View<DataType,
                                    typename traits::array_layout,
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
