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
#ifndef H_scorecMD_Memory
#define H_scorecMD_Memory

#include <stdlib.h>
void scorecMD_memoryInit(void );
void scorecMD_memoryExit(void );

#ifdef __cplusplus
extern "C" {
#endif

void * scorecMD_malloc(size_t size);
void * scorecMD_calloc(size_t num, size_t size);
void scorecMD_free(size_t size,void *p);

#ifdef __cplusplus
}
#endif


#endif
