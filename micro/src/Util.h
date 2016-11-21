#ifndef AMSI_UTIL_H_
#define AMSI_UTIL_H_

#include <limits>

template <class Container>
void erase_at(Container & c, size_t pos)
{
  typename Container::iterator it = c.begin();
  std::advance(it,pos);
  c.erase(it);
}

template <typename D, template<typename T, typename All = std::allocator<T> > class Container>
  void erase_ptr_at(Container<D*> & c, size_t pos)
{
  typename Container<D*>::iterator it = c.begin();
  std::advance(it,pos);
  D* d = *it;
  c.erase(it);
  delete d;
  d = NULL;
}

template <typename T>
void clean_delete_ptr(T*& data)
{
  delete data;
  data = NULL;
}

template <typename T>
void clean_delete_arr_ptr(T*& data)
{
  delete [] data;
  data = NULL;
}

inline bool close(double x, double y)
{
  return fabs(x - y) < std::numeric_limits<double>::epsilon();
}

template <typename Container>
inline void zero(Container & c)
{
  std::fill(c.begin(),c.end(),0);
}

#endif
