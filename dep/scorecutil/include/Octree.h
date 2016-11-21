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
#ifndef _OCTREEH_
#define _OCTREEH_

#include <list>
#include "o_internals.h"


Octree* Octree_Create(int, /* max. number of elements allowed in an octant */
		      double *origin,   /* smallest x,y, z of model's bounding box */ 
		      double *size,    /* size in x, y, z of model bounding box */
		      BBFunction BB,
		      CentroidFunction Centroid,
		      InEleFunction InEle);

void  Octree_Delete(Octree *);  
void  Octree_Insert(void *, Octree *);
void  Octree_Arrange(Octree *);
void * Octree_Search(double *, Octree *);
void  Octree_SearchAll(double *, Octree *, std::list<void*>*);

extern void mEntityBB(void *a , double* outMin, double* outMax);
extern int mEntityInEle(void *a , const double* outX);
extern void mEntityCentroid( void *a , double* outX);


#endif
