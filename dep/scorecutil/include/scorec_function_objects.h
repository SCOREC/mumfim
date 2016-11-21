#ifndef SCOREC_FUNCTION_OBJECTS_HPP
#define SCOREC_FUNCTION_OBJECTS_HPP

struct greater_abs {
  template <typename T>
  bool operator()(const T& a, const T& b) const {
    using std::abs;
    return abs(a) > abs(b);
  }
};


// This class makes it easy
// to index a flat buffer as though
// it is two-dimensional.
// Memory is in row-major order.
// NOTE: Indexing is by (row, column).
template <typename T>
class array_ref_2d {
private:
  T* const data_;
  const size_t row_length_;
public:
  array_ref_2d(T* data, size_t row_length)  
    : data_(data), row_length_(row_length) {}
  T& operator()(size_t row, size_t col) const {
    assert(row < row_length_);
    return data_[row_length_ * row + col];
  }
};

#endif
