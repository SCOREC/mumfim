#ifndef BIO_TRUSS_INTEGRATOR_H_
#define BIO_TRUSS_INTEGRATOR_H_

#include "apfUtil.h"
#include "ElementalSystem.h"

#include <apf.h>
#include <apfDynamicMatrix.h>
#include <apfField.h>

#include <cassert>

namespace bio
{
  /**
   * An integrator class to produce elemental jacobian matrices
   *  and force vectors for trusses using the provided FiberReaction.
   */
  class TrussIntegrator : public apf::Integrator
  {
  protected:
    apf::Mesh * msh;
    apf::MeshEntity * ment;
    apf::MeshElement * crnt_elmnt;
    apf::Field * u;
    
    FiberReaction * fbr_rctn;
    
    double lngth;
    double lngth_o;

    int nedofs;
    int dim;
    ElementalSystem es;
  public:
    TrussIntegrator(int o,
		    apf::Field * f,
		    FiberReaction * r)
      : apf::Integrator(o)
      , msh(NULL)
      , crnt_elmnt(NULL)
      , u(f)
      , fbr_rctn(r)
      , lngth()
      , lngth_o()
      , nedofs()
      , dim()
      , es()
    {
      assert(fbr_rctn);
      msh = apf::getMesh(f);
      dim = msh->getDimension();
      nedofs = dim*2;
    }
    // consider making an element since that makes getting the node values a bit easier?
    void inElement(apf::MeshElement * me)
    {
      ment = apf::getMeshEntity(me);
      lngth_o = apf::measure(me);
      lngth = calcDeformedLength(u,ment);
      es.resize(nedofs);
      es.zero();
    }
    void atPoint(const apf::Vector3 & p, double w, double dV)
    {
      apf::MeshEntity * vs[2];
      msh->getDownward(ment,0,&vs[0]);
      apf::Vector3 crds[2];
      getCoords(msh,&vs[0],&crds[0],2);
      apf::Vector3 frcs = (crds[1] - crds[0]) / lngth;
      
      double f = fbr_rctn->F(lngth,lngth_o);
      double fl = f / lngth;
      double dfdl = fbr_rctn->dFdl(lngth,lngth_o);
      double dfdl_fl = dfdl - fl;

      double frc = 0.0;
      for(int ii = 0; ii < dim; ii++)
      {
	frc = frcs[ii] * f;
	es.addToVector(ii      ,-frc);
	es.addToVector(2*dim+ii, frc);
      }

      apf::Matrix3x3 rctn = apf::tensorProduct(frcs,frcs*dfdl_fl) + eye()*fl;
      double op = -1.0;
      for(int ii = 0; ii < 2; ii++)
      {
	op *= -1.0;
	for(int jj = 0; jj < 2; jj++)
	{
	  op *= -1.0;
	  for(int kk = 0; kk < dim; kk++)
	    for(int ll = 0; ll < dim; ll++)
	      es.addToMatrix(ii*dim + kk, jj*dim + ll, rctn[ii][jj] * op);
	}
      }
    }
    
    const ElementalSystem * getElementalSystem() {return &es;}
  };
}

#endif
