/* 
   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of Trellis written and maintained by the 
   Scientific Computation Research Center (SCOREC) at Rensselaer Polytechnic
   Intitute, Troy, NY, USA.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the Rensselaer SCOREC Public License.

   This program is distributed in the hope that it will be useful, 
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.
   
   You should have received a copy of the Rensselaer SCOREC Public License
   along with this program; if not, write to Rensselaer Polytechnic Institure,
   110 8th Street, SCOREC, Troy, NY  12180, USA
*/
#ifndef _OctreeCreate2_h
#define _Octree_Create2_h

#include "Octree.h"

template <class Iterator>
Octree* OctreeCreate2 (Iterator itbeg,
		       Iterator itend,
		       int   number)
{
  double min[3] = {1.e24,1.e24,1.e24};
  double max[3] = {-1.e24,-1.e24,-1.e24};
  double locmin[3],locmax[3];
  double size  [3];
  int i;
  Iterator it = itbeg;
  while (it != itend)
    {
      mEntityBB (*it, locmin,locmax);
      for (i=0;i<3;i++)max[i] = (max[i]>locmax[i])?max[i]:locmax[i];
      for (i=0;i<3;i++)min[i] = (min[i]<locmin[i])?min[i]:locmin[i];
      ++it;
    }
  size[0] = max[0] - min[0];
  size[1] = max[1] - min[1];
  size[2] = max[2] - min[2];
  Octree* o = Octree_Create( number , min, size, mEntityBB, 
			     mEntityCentroid, mEntityInEle); 
  it = itbeg;
  while (it != itend)
    {
      Octree_Insert (*it, o);
      ++it;
    }

  Octree_Arrange(o);

  return o;
}
#endif
