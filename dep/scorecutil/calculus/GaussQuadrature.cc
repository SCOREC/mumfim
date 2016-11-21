/******************************************************************************

  (c) 2004-2013 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE-SCOREC file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

*******************************************************************************/

#include "GaussQuadrature.h"
#include "IntPt.h"
#include "GaussLegendreSimplex.h"

int numIntegrationPoints(int elementType, int order)
{
  switch(elementType)
    {
    case VERTEX : return 1;
    case EDGE   : return order/2+1;
    case TRI    : return (order<1)?1:getNGQTPts(order);
    case TET    : return (order<1)?1:getNGQTetPts(order);
    case QUAD   : return getNGQQPts(order);
    case HEX    : return getNGQHPts(order);
    case PRISM  : return getNGQPrismPts(order);
    case PYRAMID : return getNGQPyramidPts(order);
    default: throw 1;
    }
}

