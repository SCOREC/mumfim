#ifndef TENSORIO_H
#define TENSORIO_H

#include "Tensor.h"
#include <iostream>

template <class SpaceType, class ElementType>
std::ostream& operator<<(std::ostream& stream, const Tensor<SpaceType,ElementType>& tensor)
{
  typedef Tensor<SpaceType,ElementType> tensor_t;
  stream << '{';
  stream << tensor[0];
  for (size_t i=1; i < tensor_t::size; ++i)
  {
    stream << ',';
    stream << tensor[i];
  }
  stream << '}';
}

#endif
