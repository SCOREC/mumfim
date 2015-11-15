#ifndef BIO_TRUSS_INTEGRATOR_H_
#define BIO_TRUSS_INTEGRATOR_H_

#include "apfUtil.h"

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
    apf::Field * lngths;
    apf::Mesh * msh;
    apf::MeshEntity * ment;
    apf::MeshElement * crnt_elmnt;

    FiberReaction * fbr_rctn;
    
    double lngth;
    double lngth_o;

    int nedofs;
    int dim;
    apf::DynamicMatrix ke;
    apf::DynamicMatrix fe;
  public:
    TrussIntegrator(int o,
		    apf::Field * lngth_field,
		    FiberReaction * r)
      : apf::Integrator(o)
      , lngths(lngth_field)
      , msh(NULL)
      , crnt_elmnt(NULL)
      , fbr_rctn(r)
      , lngth()
      , lngth_o()
      , nedofs()
      , dim()
      , ke()
    {
      assert(lngths);
      assert(fbr_rctn);
      msh = apf::getMesh(lngths);
      dim = msh->getDimension();
      nedofs = dim*2;
    }
    void inElement(apf::MeshElement * me)
    {
      ment = apf::getMeshEntity(me);
      lngth_o = apf::getScalar(lngths,ment,0);
      lngth = apf::measure(me);

      ke.setSize(nedofs,nedofs);
      ke.zero();

      fe.setSize(nedofs);
      fe.zero();
    }
    void atPoint(const apf::Vector3 & p, double w, double dV)
    {
      apf::Vector3 frcs;
      calcEdgeVertDiffs(msh,ment,frcs);
      frcs = frcs / lngth;
      
      double f = fbr_rctn->F(lngth,lngth_o);
      double fl = f / lngth;
      double dfdl = fbr_rctn->dFdl(lngth,lngth_o);
      double dfdl_fl = dfdl - fl;

      for(int ii = 0; ii < dim; ii++)
      {
	fe(ii        ) -= frcs(ii) * f;
	fe(2*dim + ii) += frcs(ii) * f;
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
	      ke(ii*dim + kk, jj*dim + ll) += rctn[ii][jj] * op;
	}
      }
    }

    ElementalSystem * getElementalSystem()
    {
      return new FacadeElementalSystem(ke,fe);
    }
  };
}

#endif
