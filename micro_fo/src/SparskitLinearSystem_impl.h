#include "SparskitLinearSystem.h"

namespace bio
{
  template <typename T1, typename T2>
    void setVecValues(skVec * vec, T1 & fe, T2 & rws, int nmrws, bool add = false)
  {
    if(add)
      for(int ii = 0; ii < nmrws; ii++)
	(*vec)[rws[ii]] += fe[ii];
    else
      for(int ii = 0; ii < nmrws; ii++)
	(*vec)[rws[ii]] = fe[ii];
  }

  template <typename T1, typename T2>
    void setMatValues(skMat * mat, T1 & ke, T2 & rws, int nmrws, T2 & cls, int nmcls, bool add = false)
  {
    if(add)
      for(int ii = 0; ii < nmrws; ii++)
	for(int jj = 0; jj < nmcls; jj++)
	  (*mat)(rws[ii],cls[ii]) += ke(ii,jj);
    else
      for(int ii = 0; ii < nmrws; ii++)
	for(int jj = 0; jj < nmcls; jj++)
	  (*mat)(rws[ii],cls[ii]) = ke(ii,jj);
  }
}
