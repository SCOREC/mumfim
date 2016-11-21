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
#include "scorecMD_memory.h"
#include "scorecSC_mem.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

int scorec_gmem_alloc = 0;
int scorec_gmem_chunks = 0;

#define MAX_FAST_ALLOC_SIZE 20

scorecSC_Allocator * scorecSC_a[MAX_FAST_ALLOC_SIZE];

void scorecMD_memoryInit(void )
{
  int i;

  for(i = 0; i < MAX_FAST_ALLOC_SIZE; i++)
    scorecSC_a[i] = 0;
}

void scorecMD_memoryExit(void )
{
}

extern "C"
void * scorecMD_malloc(size_t size)
{
  int which;
  void *p;
  if(size%4)
    fprintf(stderr,"allocation error 2\n");

  which = size/4;
  if(which<MAX_FAST_ALLOC_SIZE){
    if(!scorecSC_a[which])
      scorecSC_a[which] = scorecSC_newAllocator(which*4,1000);
    p = scorecSC_allocate(scorecSC_a[which]);
  }  else 
    p = malloc(size);

#ifdef INSTRUMENT
  scorec_gmem_alloc += size;
  scorec_gmem_chunks += 1;
#endif

  if(!p && size)  /* see AIX note */
    std::cerr << "scorecMD_malloc - malloc failed - out of memory\n";
  return p;
}

extern "C"
void scorecMD_free(size_t size,void *p)
{
  int which;
  which = size/4;
  if(which < MAX_FAST_ALLOC_SIZE)
    scorecSC_deallocate(scorecSC_a[which],p);
  else
    free(p);
}

extern "C"
void * scorecMD_calloc(size_t num, size_t size)
{
  int which,i;
  void *p;

  size = size*num;
  if(size%4)
    fprintf(stderr,"allocation error 2\n");

  which = size/4;
  if(which<MAX_FAST_ALLOC_SIZE){
    if(!scorecSC_a[which])
      scorecSC_a[which] = scorecSC_newAllocator(which*4,1000);
    p = scorecSC_allocate(scorecSC_a[which]);
    } else
      p = malloc(size);

  for(i=0; i < which; i++)
    ((void**)p)[i] = 0;

#ifdef INSTRUMENT
  scorec_gmem_alloc += size;
  scorec_gmem_chunks += 1;
#endif

  if(!p && num && size)  /* see AIX note */
    std::cerr << "scorecMD_calloc - calloc failed - out of memory\n";
  return p;
}

/* AIX note */
/* malloc and calloc on AIX return null when a size of zero is */
/* passed, this is non-standard behavior, so the checks above */
/* have been added to deal with this */
