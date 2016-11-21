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
#include "scorecSC_mem.h"
#include <stdlib.h>
#include <stdio.h>

/* note: allocSize must be a multiple of sizeof(void*)
   */

void *scorecSC_allocNewBlock(int allocSize, int numAlloc)
     /* allocate the new block of memory and set up pointers for each
	allocatable piece in the new block to point to the next 
	available block */
{
  int i,allocInc;
  void **block, **bp;
  allocInc = allocSize/sizeof(void*);
#ifdef DEBUG
  allocInc++; /* to put in guard block */
  allocSize+=sizeof(void*);
#endif
  block = (void**)malloc(allocSize*numAlloc);
  if(block){
    bp = block;
    for(i=0; i < numAlloc-1; i++){
      *(bp) = bp+allocInc; /* store pointer to next block */
#ifdef DEBUG
      *(int*)(bp+allocInc-1) = allocSize-sizeof(void*);
#endif
      bp += allocInc;
    }
    *(bp) = 0; /* last one*/
#ifdef DEBUG
    *(int*)(bp+allocInc-1) = allocSize-sizeof(void*);
#endif

  }
  return block;
}

scorecSC_Allocator * scorecSC_newAllocator(int allocSize, int numInBlock)
{
  scorecSC_Allocator *block;
  block = (scorecSC_Allocator*)malloc(sizeof(scorecSC_Allocator));
  block->d_allocSize = allocSize;
  block->d_numInBlock = numInBlock;
  block->d_firstBlock = (scorecSC_MemBlockList*)malloc(sizeof(scorecSC_MemBlockList));
  block->d_firstBlock->d_nextBlock = 0;
  block->d_firstBlock->d_block = scorecSC_allocNewBlock(allocSize,numInBlock);
  block->d_nextFree = block->d_firstBlock->d_block;

  block->d_lastBlock = block->d_firstBlock;
  block->d_numAlloc = numInBlock;
  block->d_numFree = numInBlock;
  return block;
}

void scorecSC_deleteAllocator(scorecSC_Allocator *b)
{
  scorecSC_MemBlockList *bl, *bln;
  bl = b->d_firstBlock;
  while(bl){
    free(bl->d_block);
    bln = (scorecSC_MemBlockList*)bl->d_nextBlock;
    free(bl);
    bl = bln;
  }
}

void * scorecSC_allocate(scorecSC_Allocator *b)
{
  void *ret;
  scorecSC_MemBlockList *bl;
  ret = b->d_nextFree;
  if(!ret){
    bl = (scorecSC_MemBlockList*)malloc(sizeof(scorecSC_MemBlockList));
    bl->d_block = scorecSC_allocNewBlock(b->d_allocSize,b->d_numInBlock);
    bl->d_nextBlock = 0;
    b->d_lastBlock->d_nextBlock = bl;
    b->d_lastBlock = (scorecSC_MemBlockList*)b->d_lastBlock->d_nextBlock;
    b->d_nextFree = bl->d_block;
    ret = b->d_nextFree;
    b->d_numAlloc += b->d_numInBlock;
    b->d_numFree += b->d_numInBlock;
  } 
  b->d_nextFree = *((void**)(b->d_nextFree));
  b->d_numFree--;
  if(b->d_numFree < 0)
    fprintf(stderr, "allocator error 1\n");
#ifdef DEBUG
  if(*( ((int*)ret)+b->d_allocSize/sizeof(int)) != b->d_allocSize)
    fprintf(stderr,"allocate error 3\n");
#endif
  return ret;
}

void scorecSC_deallocate(scorecSC_Allocator *b, void *mem)
{
#ifdef DEBUG
  if(*( ((int*)mem)+b->d_allocSize/sizeof(int)) != b->d_allocSize)
    fprintf(stderr,"allocate error 4\n");
#endif
  *((void**)(mem)) = (void*)b->d_nextFree;
  b->d_nextFree = mem;
  b->d_numFree++;
}

void scorecSC_allocatorState(scorecSC_Allocator *b)
{
  fprintf(stderr,"allocSize = %d\n",b->d_allocSize);
  fprintf(stderr,"%d allocated\n",b->d_numAlloc);
  fprintf(stderr,"%d free\n",b->d_numFree);
}
