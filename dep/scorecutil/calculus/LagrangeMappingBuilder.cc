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
#include "LagrangeMappingBuilder.h"
#include "LagrangeMapping.h"
#ifdef SIM
#include "MSMapping.h"
#endif

Mapping* LagrangeMappingBuilder::BuildMapping(pEntity e, int SystemOfCoordinates) 
{
  if(SystemOfCoordinates != 0)
    return new CylindricalCoordinatesLagrangeMapping(e);
  
  assert(e);
  switch(M_GetElementType(e)) {
  case QUAD:
  case HEX:
    break;
  case EDGE:
  case TRI:
  case TET:
    return new ConstantLagrangeMapping(e);
  default: 
    break; 
  }
  return new LagrangeMapping(e);
}

Mapping* LagrangeMappingBuilder::BuildMappingByKnots(pEntity e, int SystemOfCoordinates, std::vector<mPoint> *knots)
{
  if(SystemOfCoordinates != 0)
    return new CylindricalCoordinatesLagrangeMapping(e);

  assert(e);
  //switch(M_GetElementType(e)) {
  //case QUAD:
  //case HEX:
    //break;
  //case EDGE:
  //case TRI:
  //case TET:
    //return new ConstantLagrangeMapping(e);
  //default:
   // break;
  //}
  return new LagrangeMapping(e, knots);
}


  
std::string LagrangeMappingBuilder::getMappingBuilderName() const {
  return MType();
}






