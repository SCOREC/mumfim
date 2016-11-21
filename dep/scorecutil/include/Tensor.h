#ifndef TENSOR_H
#define TENSOR_H

#include <cstddef>
#include <iostream>

typedef double Scalar;

template <class SpaceType>
class Index
{
  public:
    enum { limit = SpaceType::dimensions };
    Index():value(0) {}
    const size_t& Value() const { return value; }
    bool operator++()
    {
      ++value;
      if (value == limit)
      {
        value = 0;
        return false;
      }
      return true;
    }
  private:
    size_t value;
};

template <size_t Dimensions>
class R
{
  public:
    enum { dimensions = Dimensions };
    typedef Index<R<Dimensions> > index_type;
};

template <class DefinedType>
class IndexTypeOf
{
  public:
    typedef typename DefinedType::index_type Result;
};

template <>
class IndexTypeOf<Scalar>
{
  public:
    typedef Scalar Result;
};

template <class Expression>
class TensorExpression
{
};

template <class TermType>
class Indexed : public TensorExpression<Indexed<TermType> >
{
  public:
    typedef Indexed<TermType> self_type;
    typedef typename TermType::tensor_type sub_type;
    typedef typename sub_type::index_type sub_index_type;
    typedef typename sub_type::element_type tensor_type;
    typedef typename IndexTypeOf<tensor_type>::Result index_type;
    typedef Indexed<self_type> indexed_type;
    Indexed(TermType& _term, sub_index_type& _i):term(_term),i(_i) {}
    indexed_type operator[](index_type& idx)
    {
      return indexed_type(*this,idx);
    }
    tensor_type& Eval()
    {
      return term.Eval()[i.Value()];
    }
    const tensor_type& Eval() const
    {
      return term.Eval()[i.Value()];
    }
    self_type& operator=(const self_type& other)
    {
      do {
        Eval() = other.Eval();
      } while (Increment());
    }
    template <class O>
    self_type& operator=(const TensorExpression<O>& _other)
    {
      const O& other = static_cast<const O&>(_other);
      do {
        Eval() = other.Eval();
      } while (Increment());
    }
    bool Increment()
    {
      if (++i)
      {
        return true;
      }
      return term.Increment();
    }
  private:
    TermType& term;
    sub_index_type& i;
};

template <class TermType>
class ConstIndexed : public TensorExpression<ConstIndexed<TermType> >
{
  public:
    typedef ConstIndexed<TermType> self_type;
    typedef typename TermType::tensor_type sub_type;
    typedef typename sub_type::index_type sub_index_type;
    typedef typename sub_type::element_type tensor_type;
    typedef typename IndexTypeOf<tensor_type>::Result index_type;
    typedef ConstIndexed<self_type> indexed_type;
    ConstIndexed(const TermType& _term, const sub_index_type& _i):term(_term),i(_i) {}
    indexed_type operator[](index_type& idx)
    {
      return indexed_type(*this,idx);
    }
    const tensor_type& Eval()
    {
      return term.Eval()[i.Value()];
    }
  private:
    const TermType& term;
    const sub_index_type& i;
};

template <class SpaceType, class ElementType>
class Tensor : public TensorExpression<Tensor<SpaceType,ElementType> >
{
  public:
    typedef Tensor<SpaceType,ElementType> self_type;
    typedef self_type tensor_type;
    enum { size = SpaceType::dimensions };
    typedef ElementType element_type;
    typedef element_type data_type[size];
    typedef Index<SpaceType> index_type;
    typedef Indexed<self_type> indexed_type;
    typedef ConstIndexed<self_type> const_indexed_type;
    Tensor() {}
    element_type& operator[](size_t i) { return _data[i]; }
    const element_type& operator[](size_t i) const { return _data[i]; }
    indexed_type operator[](index_type& idx)
    {
      return indexed_type(*this,idx);
    }
    const_indexed_type operator[](const index_type& idx) const
    {
      return const_indexed_type(*this,idx);
    }
  private:
    data_type _data;
};

template <class SpaceType, class ElementType>
class Indexed<Tensor<SpaceType,ElementType> > : public TensorExpression<Indexed<Tensor<SpaceType,ElementType> > >
{
  public:
    typedef Tensor<SpaceType,ElementType> TermType;
    typedef Indexed<TermType> self_type;
    typedef typename TermType::tensor_type sub_type;
    typedef typename sub_type::index_type sub_index_type;
    typedef typename sub_type::element_type tensor_type;
    typedef typename IndexTypeOf<tensor_type>::Result index_type;
    typedef Indexed<self_type> indexed_type;
    Indexed(TermType& _term, sub_index_type& _i):term(_term),i(_i) {}
    indexed_type operator[](index_type& idx)
    {
      return indexed_type(*this,idx);
    }
    tensor_type& Eval()
    {
      return term[i.Value()];
    }
    const tensor_type& Eval() const
    {
      return term[i.Value()];
    }
    self_type& operator=(const self_type& other)
    {
      do {
        Eval() = other.Eval();
      } while (Increment());
    }
    template <class O>
    self_type& operator=(const TensorExpression<O>& _other)
    {
      const O& other = static_cast<const O&>(_other);
      do {
        Eval() = other.Eval();
      } while (Increment());
    }
    bool Increment()
    {
      return ++i;
    }
  private:
    TermType& term;
    sub_index_type& i;
};

template <class SpaceType, class ElementType>
class ConstIndexed<Tensor<SpaceType,ElementType> > : TensorExpression<ConstIndexed<Tensor<SpaceType,ElementType> > >
{
  public:
    typedef Tensor<SpaceType,ElementType> TermType;
    typedef ConstIndexed<TermType> self_type;
    typedef typename TermType::tensor_type sub_type;
    typedef typename sub_type::index_type sub_index_type;
    typedef typename sub_type::element_type tensor_type;
    typedef typename IndexTypeOf<tensor_type>::Result index_type;
    typedef ConstIndexed<self_type> indexed_type;
    ConstIndexed(const TermType& _term, const sub_index_type& _i):term(_term),i(_i) {}
    indexed_type operator[](index_type& idx)
    {
      return indexed_type(*this,idx);
    }
    const tensor_type& Eval()
    {
      return term[i.Value()];
    }
  private:
    const TermType& term;
    const sub_index_type& i;
};

template <class LHS, class RHS>
class Product : public TensorExpression<Product<LHS,RHS> >
{
  public:
    typedef const LHS& lhs_t;
    typedef const RHS& rhs_t;
    Product(lhs_t _lhs, rhs_t _rhs):lhs(_lhs),rhs(_rhs) {}
    Scalar Eval() const { return (lhs.Eval() * rhs.Eval()); }
  private:
    lhs_t lhs;
    rhs_t rhs;
};

template <class LHS, class RHS>
Product<LHS,RHS> operator*(const TensorExpression<LHS>& lhs, const TensorExpression<RHS>& rhs)
{
  return Product<LHS,RHS>(static_cast<const LHS&>(lhs), static_cast<const RHS&>(rhs));
}

template <class RHS>
class Product<Scalar,RHS> : public TensorExpression<Product<Scalar,RHS> >
{
  public:
    typedef const Scalar& lhs_t;
    typedef const RHS& rhs_t;
    Product(lhs_t _lhs, rhs_t _rhs):lhs(_lhs),rhs(_rhs) {}
    Scalar Eval() const { return (lhs * rhs.Eval()); }
  private:
    lhs_t lhs;
    rhs_t rhs;
};

template <class RHS>
Product<Scalar,RHS> operator*(const Scalar& lhs, const TensorExpression<RHS>& rhs)
{
  return Product<Scalar,RHS>(lhs, static_cast<const RHS&>(rhs));
}

template <class IndexType, class TermType>
class Reduction
{
  public:
    typedef const TermType& term_t;
    typedef IndexType& index_t;
    Reduction(term_t _term, index_t _i):term(_term),i(_i) {}
    Scalar Eval()
    {
      Scalar result(0);
      do {
        result += term.Eval();
      } while(++i);
      return result;
    }
  private:
    term_t term;
    index_t i;
};

template <class IndexType, class TermType>
Reduction<IndexType,TermType> Sum(IndexType& index, const TermType& term)
{
  return Reduction<IndexType,TermType>(term,index);
}

#endif
