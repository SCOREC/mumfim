/* 
   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of MCTK written and maintained by the 
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
#ifndef _GAUSS_INTEGRATOR_H_
#define _GAUSS_INTEGRATOR_H_


#include "Integrator.h"
#include "Mapping.h"
#include <vector>

template <class integrable>
double* GaussIntegration(SCOREC::Util::Mapping& mapping, const integrable& intgb )
{
  int sizeResult = intgb.size();
  int order = intgb.order(&mapping);
  mVector pos;
  double weight;

  
  double * result = NULL;
  try {
    result = new double[sizeResult];
    std::vector<double> tempResult(sizeResult);
    
    SCOREC::Util::GaussIntegrator myIntegrator(mapping.getEntity());
    int nbPt = myIntegrator.nbIntegrationPoints(order);
    std::fill(result, result + sizeResult, 0.0);
    
    for (size_t iPt = 0; iPt != nbPt; ++iPt) {
      myIntegrator.iPoint(iPt, order, pos[0], pos[1], pos[2], weight);
      intgb.eval(pos, &mapping, &(tempResult[0]));
      double detJac = mapping.detJac(pos[0],pos[1],pos[2]);
      
      for (size_t i = 0 ; i != sizeResult; ++i) 
	result[i] += tempResult[i] * weight * detJac;
    }
  } catch(...) {
    delete [] result;
    throw;
  }
  return result;
}
#endif




