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
#include <stdio.h>
#include <iostream>
#include "Mapping.h"
#include "BezierMappingBuilder.h"
#include "BezierMapping.h"
#ifndef SIM
#include "mEntity.h"
#include "mVertex.h"
#include "mAOMD.h"
#else
#include "MSMapping.h"
#endif

   /*
     change a little bit here, return 0 - as straight-sided
     ortherwise return the order of geometric approximation
   */
#ifndef SIM
  int isEdgeAttachPoint(AOMD::mEntity *me)
  {
    AOMD::mAttachablePointVector *ap;    
    ap = (AOMD::mAttachablePointVector*)me->getData(AOMD::AOMD_Util::Instance()->getEmod());
    if(!ap)
      return 0;    
    int size = (ap->v).size();   
    //printf("The size is: %d\n",size);
    return size;
  }
  
  int isTriFaceAttachPoint(AOMD::mEntity *me)
  { 
    int size = me->size(1);
    int order;
    for(int i=0; i < size; i++)
      {
	AOMD::mEntity *e = me->get(1, i);
	order = isEdgeAttachPoint(e);
	if(order!=0){
	  // printf("The Loop face order is: %d", order);
	  return order;
	}
      }    
    return 0;
  }
  
  int isTetRegionAttachPoint(AOMD::mEntity *me)
  { 
    int size = me->size(1);
    for(int i=0; i < size; i++)
      {
	AOMD::mEntity *e = me->get(1, i);
	int order = isEdgeAttachPoint(e);
	if(order!=0)
	  return order;
      }    
    return 0;
  }
  
  int isUseBezier(AOMD::mEntity *me)
  {
    
    switch(me->getType())
      {
      case AOMD::mEntity::VERTEX : return 0;
      case AOMD::mEntity::EDGE   : return isEdgeAttachPoint(me); break;
      case AOMD::mEntity::TRI    : return isTriFaceAttachPoint(me); break;
      case AOMD::mEntity::TET    : return isTetRegionAttachPoint(me); break;
      }    
    return 0;
  }
#endif
  Mapping *BezierMappingBuilder::BuildMapping(pEntity e, int SystemOfCoordinates) 
  {    
#ifndef SIM
    int gorder = isUseBezier(e);
    if(gorder == 0)
      {
	//	MappingBuilder lmb;
	//	return lmb.BuildMapping(e,SystemOfCoordinates);
      }
    //printf("The Bezier order is: %d\n", gorder);
    switch(e->getType())
      {
      case AOMD::mEntity::EDGE :
	return new BezierMapping(e, gorder+1);
	break;
      case AOMD::mEntity::TRI :
	return new BezierMapping(e, gorder+1);
	break;
      case AOMD::mEntity::TET :
	return new BezierMapping(e, gorder+1);
	break;
      default: break; 
	throw 1;
      }
#else
    printf("BezierMappingBuilder currently not working with MeshSim.\n");
    throw 1;
#endif
    return NULL;
  }
  std::string BezierMappingBuilder::getMappingBuilderName() const {
    return MType();
  }






