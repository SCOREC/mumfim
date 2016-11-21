#include "NonLinFibMtx.h"
#include "NonFiniteElement.h"

// femanalysis
#include <FEA.h>
#include <ConvenienceFunctions.h>

// APF
#include <apf.h>
#include <apfSIM.h>

// Simmetrix
#include <MeshSim.h>

#include <memory.h>
#include <vector>
#include <list>
#include <map>
#include <iostream>

namespace Biotissue {

// obtain the element stress at gauss point
  void NonLinFibMtx::ElementStress(apf::MeshEntity * me, 
				   apf::Matrix3x3 & pk2)
{
  // might just want to pass as apf element or entity...
  apf::MeshElement * melm = apf::createMeshElement(apf_mesh,me);
  apf::Element * e = apf::createElement(apf_primary_field,me);
  //int num_nodes = apf::countNodes(e);

  apf::NewArray<apf::Vector3> disps;
  apf::getVectorNodes(e,disps);
  
  int num_gauss_pts = apf::countIntPoints(melm,integration_order);
  for(int gauss_index = 0; gauss_index < num_gauss_pts; gauss_index++) 
  {
    apf::Vector3 gauss_pt;
    apf::getIntPoint(melm,integration_order,gauss_index,gauss_pt);
  
    apf::Matrix3x3 J;
    apf::getJacobian(melm,gauss_pt,J);
    double det_jac = apf::getJacobianDeterminant(J,apf::getDimension(melm));
    
    apf::Matrix3x3 def_grad;
    double def_grad_det;
    DeformationGradient(e,gauss_pt,def_grad,def_grad_det);

    apf::Matrix3x3 C;
    RightCauchy(def_grad,C);
    PK2Stress(C,def_grad_det,poisson_ratio,shear_modulus,pk2);
  }
}
  
}
