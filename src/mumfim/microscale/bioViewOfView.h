#ifndef BIO_VIEW_OF_VIEW_H
#define BIO_VIEW_OF_VIEW_H
#include <Kokkos_Core.hpp>
#include <vector>
#include <cassert>
#include "bioMicroTypeDefinitions.h"
namespace bio
{
  template <typename DataType,
            typename LocalOrdinal=bio::LocalOrdinal,
            typename ExeSpace = Kokkos::DefaultExecutionSpace>
  class ViewOfView
  {
    public:
    using d_inner_type = Kokkos::View<DataType *, ExeSpace>;
    using h_inner_type = typename d_inner_type::HostMirror;
    using dd_type = Kokkos::View<d_inner_type *, ExeSpace>;
    using dh_type = typename dd_type::HostMirror;
    using hh_type = typename Kokkos::View<h_inner_type *, ExeSpace>::HostMirror;

    ViewOfView(std::vector<std::vector<DataType>> data)
    {
      dd_view = dd_type("dd", data.size());
      hh_view = hh_type("hh", data.size());
      dh_view = create_mirror_view(dd_view);
      assert(hh_view.extent(0) == data.size());
      assert(dd_view.extent(0) == data.size());
      assert(dh_view.extent(0) == data.size());
      for (size_t i = 0; i < hh_view.extent(0); ++i)
      {
        dh_view(i) = d_inner_type("inner", data[i].size());
        hh_view(i) = create_mirror_view(dh_view(i));
        for (size_t j = 0; j < data[i].size(); ++j)
        {
          hh_view(i)(j) = data[i][j];
        }
      }
      deep_copy(dd_view, dh_view);
      syncToDevice();
    }
    // the default constructor creates a view
    // with zero size for the device and the host
    // this should help the user understand that
    // there is nothing there if they try to take
    // the extent
    ViewOfView()
    {
      dd_view = dd_type("dd", 0);
      hh_view = hh_type("hh", 0);
      syncToDevice();
    }
    void syncToDevice()
    {
      for (size_t i = 0; i < dh_view.extent(0); ++i)
      {
        deep_copy(dh_view(i), hh_view(i));
      }
    }
    void syncToHost()
    {
      for (size_t i = 0; i < dh_view.extent(0); ++i)
      {
        deep_copy(hh_view(i), dh_view(i));
      }
    }
    hh_type & h_view() { return hh_view; }
    dd_type & d_view() { return dd_view; }
    private:
    // warning these cannot directly be used in an internal lambda
    // because the lambda will capture the this pointer. If operations
    // are added that modify these views using a lambda, first assign
    // to a new variable and use that in the lambda
    hh_type hh_view;
    dh_type dh_view;
    dd_type dd_view;
  };
}  // namespace bio
#endif
