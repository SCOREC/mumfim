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
#include <stdlib.h>
#include <stdio.h>
#include <list>
#include "o_internals.h"
#include "Octree.h"
#include "mPoint.h"
#include "LagrangeMapping.h"

using std::list;

void free_buckets(octantBucket *);

Octree* Octree_Create (int maxElements, double origin[3], double size[3],   
                       BBFunction BB,
                       CentroidFunction Centroid,
                       InEleFunction InEle)
{
  Octree *myOctree = new Octree;
  initializeOctantBuckets (origin, size, maxElements,
                           &(myOctree->root), &(myOctree->info));  	
  myOctree->function_BB = BB;
  myOctree->function_centroid = Centroid;
  myOctree->function_inElement = InEle;
  return myOctree;
}


void Octree_Delete(Octree *myOctree)
{
  delete myOctree->info;
  free_buckets(myOctree->root);
  delete myOctree->root;
  delete myOctree;
  
  return;
}

void Octree_Insert(void * element, Octree *myOctree)
{
  double minBB[3], maxBB[3], centroid[3];
  octantBucket *bucket;
  (*(myOctree->function_BB))(element, minBB, maxBB);
  (*(myOctree->function_centroid))(element, centroid);

// for test 
/* printf("The centroid of the insert element is %f, %f, %f\n",centroid[0], centroid[1], centroid[2]); */
//end of test

  bucket = findElementBucket(myOctree->root, centroid);
  addElement2Bucket (bucket, element, minBB, maxBB,
                            centroid,myOctree->info);	
  // printf("Ready to return from Octree_Insert\n");
  return;
}

void Octree_Arrange(Octree *myOctree)
{  
  std::list<void *>::iterator iter;
  double minPt[3], maxPt[3];
  for(iter = myOctree->info->listAllElements.begin(); iter!= 
      myOctree->info->listAllElements.end(); iter++) {
    (*(myOctree->function_BB))(*iter, minPt, maxPt);
    insertOneBB(*iter, minPt, maxPt, myOctree->root);
  }
  myOctree->info->listAllElements.clear();
  return;
}   




void* Octree_Search(double *pt, Octree *myOctree)
{
  return searchElement(myOctree->root, pt, myOctree->info, 
		       myOctree->function_BB, myOctree->function_inElement);
}


void free_buckets(octantBucket * bucket)
{
  int i, numBuck = 8;
  ELink ptr1, ptr2;

  if(bucket->next == NULL) {
    ptr1 = bucket->lhead;
    while(ptr1 != NULL) {
      ptr2 = ptr1;
      ptr1 = ptr1->next;
      delete ptr2;
    }
    bucket->listBB.clear(); 
    return;
  }

  for(i = numBuck-1; i >= 0; i--) 
    free_buckets((bucket->next)+i);	    
  delete []bucket->next;
  return;
}


void  Octree_SearchAll(double * pt, Octree * myOctree, list<void *> * output)
{
  searchAllElements(myOctree->root, pt, myOctree->info, myOctree->function_BB,
                    myOctree->function_inElement, output);	
}

void mEntityBB(void *a , double *min, double *max)
{
  pEntity e = static_cast<pEntity>(a);

  LagrangeMapping map(e);
  mPoint p1,p2;
  map.boundingBox(p1,p2);
  for(int i=0;i<3;i++) {
    min[i] = p1(i);
    max[i] = p2(i);
  }
}

int mEntityInEle(void *a , const double *x)
{
  pEntity e = static_cast<pEntity>(a);
  LagrangeMapping map(e);
  double u,v,w;
  return map.interiorCheck(e,mPoint(x[0],x[1],x[2]),u,v,w);
}


void mEntityCentroid(void *a , double *x)
{
  pEntity e = (pEntity)a;
  LagrangeMapping map(e);
  double u,v,w;
  map.COG(u,v,w);
  map.eval(u,v,w,x[0],x[1],x[2]);
}

