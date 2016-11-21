#ifndef SCOREC_UTIL_AFFINE_SPACE_H
#define SCOREC_UTIL_AFFINE_SPACE_H

// Ben FrantzDale
// 1/5/2007
// These functions and classes are for use with affine spaces
// that is, vector spaces without a zero element. This includes
// points in space in the real world.
// Legal operations on such points are subtraction (which produces
// a difference which is a vector) and adding a difference vector.
// Difference vectors can be scalar multiplied.
// Points in affine space can combined with an affine combination
// which is a weighted sum where the weights sum to 1.

#include <numeric>
#include <iterator>

template <typename AffinePoint,
	  typename Scalar>
inline
AffinePoint const affine_combine(AffinePoint const& a, 
				 AffinePoint const& b,
				 Scalar const& aWeight) {
  return b + (a - b) * aWeight;
}


// This class makes it easy to average affine points.
// It uses Knuth's online mean algorithm. See
// http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
template <typename AffinePoint>
class AffineMean {
public:
  AffineMean() : n_(0), mean_() {}
  AffineMean(const AffineMean& m) : n_(m.n_), mean_(m.mean_) {}
  AffineMean& operator+=(AffinePoint const& x) {
    n_ += 1.0;
    // We want mean = mean + (x - mean) / n
    // i.e., mean = mean + (x - mean) * (1/n)
    // i.e., mean = affine_combine(x, mean, 1/n)
    // So...
    mean_ = affine_combine(x, mean_, 1.0/n_);
  }
  AffineMean const operator+(AffinePoint const& x) const {
    AffineMean result = *this;
    result += x;
    return result;
  }
  AffinePoint const& result() const {
    assert(n_ > 0.0);
    return mean_;
  }
  std::size_t count() const {
    return static_cast<std::size_t>(n_);
  }
private:
  double n_;
  AffinePoint mean_;
};

// TODO:
// WeightedMean takes a weighted average of affine
// points. (i.e., a convex combination) To do this it takes a function object which supplies the
// weight for each point, then divides by the sum of those weights so
// that the sum of weights is 1.0.
//
// WM_n+1 = WM_n + (w_n+1 * x_n+1 - w_n+1 * WM_n) / (sum_{i=1}^{n+1} w_i)
template <typename AffinePoint,
	  typename WeightFn>
class WeightedMean {
public:
  WeightedMean(WeightFn const& f = WeightFn())
    : weight_sum_(0.0),
      weighted_mean_(),
      weight_fn_(f)
  {
  }
  WeightedMean(WeightedMean const& other) 
    : weight_sum_(other.weight_sum_),
      weighted_mean_(other.weighted_mean_),
      weight_fn_(other.weight_fn_)
  {
  }
  WeightedMean& operator+=(AffinePoint const& x) {
    double const weight = weight_fn_(x);
    assert(weight >= 0.0);
    weight_sum_ += weight;
    weighted_mean_ = affine_combine(x, weighted_mean_, weight/weight_sum_);
  }
  WeightedMean const operator+(AffinePoint const& x) const {
    WeightedMean result = *this;
    result += x;
    return result;
  }
  AffinePoint const& result() const {
    assert(weight_sum_ > 0.0); 
    return weighted_mean_;
  }
private:
  double weight_sum_;
  AffinePoint weighted_mean_;
  WeightFn weight_fn_;
};


template <typename AffinePointIterator>
inline
typename std::iterator_traits<AffinePointIterator>::value_type
mean(AffinePointIterator beg,
     AffinePointIterator end) {
  typedef
    typename std::iterator_traits<AffinePointIterator>::value_type
    AffinePoint;
  return std::accumulate(beg, end, AffineMean<AffinePoint>()).result();
}



// A function object to read off the range of weights.
template <typename ValueIter>
struct range_reader {
  ValueIter i;
  ValueIter end;
  typedef typename std::iterator_traits<ValueIter>::value_type result_type;
  range_reader(ValueIter b, ValueIter e) 
    : i(b), end(e) {}
  result_type operator()() {
    assert(i != end);
    double const result = *i;
    ++i;
    return result;
  }
  // Allow use with other types.
  template <typename T>
  result_type operator()(const T&) {
    return operator()();
  }
};



template <typename AffinePointIterator,
	  typename WeightIterator>
inline
typename std::iterator_traits<AffinePointIterator>::value_type
weighted_mean(AffinePointIterator beg,
	      AffinePointIterator end,
	      WeightIterator beg_weight,
	      WeightIterator end_weight) {
  typedef
    typename std::iterator_traits<AffinePointIterator>::value_type
    AffinePoint;
  range_reader<WeightIterator> weight_fn(beg_weight, end_weight);
  WeightedMean<AffinePoint, range_reader<WeightIterator> > mean_computer(weight_fn);
  return std::accumulate(beg, end, mean_computer).result();
}
#endif // SCOREC_UTIL_AFFINE_SPACE_H
