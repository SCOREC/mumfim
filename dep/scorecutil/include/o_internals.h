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
#ifndef OCTCLASSIFY
#define OCTCLASSIFY
#include <list>

/* file of function prototypes and macro constants */

typedef void (*BBFunction)(void *, double*, double*);
typedef int (*InEleFunction)(void *, const double *); 
typedef void (*CentroidFunction)(void *, double *);

/* structure for list of elements in an octant */
struct Elem {
  void* region;  /* the pointer to a mesh Db region */
  double centroid[3]; /* centroid of element bounding box inside of the octant */
  double minPt[3]; /* corner of element bounding box nearest the origin */
  double maxPt[3]; /* corner of elem bound box furthest from the origin */ 
  Elem* next; /* link to next item in list, NULL if end */
};

typedef Elem *ELink;

  /* stucture for octant buckets */
  struct octantBucket {
    double minPt[3];   /*  the point with the smallest coordinates */
    double maxPt[3];   /*  the point with the biggest coordinates */
    int numElements; /* number of elements contained by bucket */
    int precision;   /* the level of precision of the bucket */
    ELink lhead; /* list of elements in bucket, if NULL -> no elements */
    std::list<void*> listBB; /* list of elements in bucket by Bounding Box */ 	
    octantBucket* next; /* link to ragged digit extensions to bucket array */
    octantBucket* parent; /* link to the parent bucket */
  };

  /* octantBucket *buckets=NULL; */

  /* structure for global information and requirment */
  struct globalInfo {
    int numBuckets; /* number of octant buckets in initial grid array */
    int maxElements; /* max. number of elements allowed in an octant */
    int maxPrecision; /* current maximum octant precision for model */
    double origin[3];   /* smallest x,y, z of model's bounding box */ 
    double size[3];    /* size in x, y, z of model bounding box */
    void* ptrToPrevElement;	
    std::list<void*> listAllElements;
  };


  struct Octree {
    globalInfo *info;
    octantBucket *root;
    BBFunction function_BB;
    InEleFunction function_inElement; 
    CentroidFunction function_centroid; 
  };
  
  void refineOctants( octantBucket *buckets,
		      globalInfo *globalPara);

  int addElement2Bucket (octantBucket *bucket, void * element, 
			 double *minBB, double *maxBB,
			 double *ele_centroid, globalInfo *globalPara);
  int subdivideOctantBucket (octantBucket *bucket, globalInfo *globalPara);
  int initializeOctantBuckets (double *orig, double *size, int maxElem,
			       octantBucket **buckets, globalInfo **globalPara);
  int checkElementInBucket (octantBucket *bucket, void * element); 
  octantBucket* findElementBucket(octantBucket *buckets, const double *pt);
  void * searchElement(octantBucket *buckets, const double *pt, 
		       globalInfo *globalPara, BBFunction BBElement, 
		       InEleFunction xyzInElement);
  int xyzInElementBB(const double *xyz, void *region, BBFunction BBElement);
  void insertOneBB(void *, double *, double *, octantBucket *);
  void * searchAllElements(octantBucket *_buckets_head, 
			   const double *_pt, 
			   globalInfo *_globalPara,
			   BBFunction BBElement, 
			   InEleFunction xyzInElement, 
			   std::list<void*> *_elements);

#endif










