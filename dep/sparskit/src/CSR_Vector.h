#ifndef _CSR_VECTOR_H
#define _CSR_VECTOR_H

#include <iostream>
#include <ostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <math.h>
#include <functional>

class CSR_Vector 
{
private : 
  std::vector<double> Array;
public :
  
  inline CSR_Vector( void ) {};
  inline CSR_Vector( const int& n )
  {
    Array.insert(Array.begin(), n, 0.0); 
  }
  inline CSR_Vector(const CSR_Vector& rhs)
    : Array (rhs.Array)
  {
  }
  inline CSR_Vector(double *vecarray_beg, double *vecarray_end)
  {
    Array.insert(Array.begin(), vecarray_beg, vecarray_end);
  }
  inline ~CSR_Vector()
  {}
  inline CSR_Vector & operator=(const CSR_Vector & rhs)
  {
    if (this != &rhs)
      Array = rhs.Array;
    return (*this);
  }  

  inline void     SetSize(const int& n)
  {
    Array.clear();
    Array.insert(Array.begin(), n, 0.0); 
  }

  inline void AddVal( const int& i, const double& val)
  {
    // FORTRAN NUMBERING ?
    Array[i-1] += val;
  }

  inline void  Set(double *vecarray_beg, double *vecarray_end)
  {
    Array.clear();
    Array.insert(Array.begin(), vecarray_beg, vecarray_end); 
  }

  inline void Append(double * vecarray_beg, double * vecarray_end)
  {
    Array.insert(Array.end(), vecarray_beg, vecarray_end);
  }
  
  inline const double*  GetArray() const
  {
    return &(*Array.begin());
  }

  inline double*  GetArray()
  {
    return &(*Array.begin());
  }
  
  inline void ZeroArray()
  {
    std::fill(Array.begin(), Array.end(), 0.0);
  }

  inline double Norm() const
  {
    return sqrt(std::inner_product(Array.begin(), Array.end(), Array.begin(), 0.0));
  }
  
  inline double Dot(const CSR_Vector& vec) const
  {
    return std::inner_product(Array.begin(), Array.end(), vec.Array.begin(), 0.0);
  }
  
  inline void PrintArray( std::ostream  & out) const
  {
    out << "###########################\n";
    out << "Vector of dimension " << Array.size() << " is:\n";
    for (unsigned int i = 0; i < Array.size() ; ++i) out << Array[i] << std::endl;
    out << "\n";
    out << "###########################\n";
  }
  
  inline double& operator[](unsigned int i)
  {
    return Array[i];
  }

  inline const  double& operator[](unsigned int i) const 
  {
    return Array[i];
  }

  inline void MultipliedBy(const double& coeff)
  {
    std::transform(Array.begin(), Array.end(),  Array.begin(), 
	      std::bind2nd(std::multiplies<double>(), coeff));
  }

  inline int size() const {return Array.size() ;}
};

#endif






