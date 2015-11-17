#ifndef BIO_TRUSS_INTEGRATOR_H_
#define BIO_TRUSS_INTEGRATOR_H_
#include "apfUtil.h"
#include "FiberNetwork.h"
#include "SparskitLinearSystem.h"
#include <apf.h>
#include <apfDynamicMatrix.h>
#include <apfField.h>
#include <cassert>
namespace bio
{
  class ElementalSystem;
  /**
   * An integrator class to produce elemental jacobian matrices
   *  and force vectors for trusses using the provided FiberReaction.
   *  The elemental systems are assembled into the global system formed
   *  by the skMat k and skVec f using the supplied numbering, which
   *  should already be numbered.
   * @todo Bill : refactor and simplify, it feels like there is too much
   *              happening in this class
   */
  class TrussIntegrator : public apf::Integrator
  {
  protected:
    apf::Mesh * msh;
    apf::MeshEntity * ment;
    apf::MeshElement * crnt_elmnt;
    apf::Field * u;
    apf::Numbering * num;
    FiberReaction * fbr_rctn;
    double lngth;
    double lngth_o;
    int nedofs;
    int dim;
    apf::NewArray<int> dofs;
    ElementalSystem es;
    skMat * k;
    skVec * f;
  public:
    TrussIntegrator(int o,
		    apf::Numbering * nm,
		    FiberReaction * r,
		    skMat * K,
		    skVec * F)
      : apf::Integrator(o)
      , msh(NULL)
      , crnt_elmnt(NULL)
      , u(NULL)
      , num(nm)
      , fbr_rctn(r)
      , lngth()
      , lngth_o()
      , nedofs()
      , dim()
      , es()
      , k(K)
      , f(F)
    {
      assert(fbr_rctn);
      assert(num);
      assert(k);
      assert(f);
      u = apf::getField(num);
      msh = apf::getMesh(u);
      dim = msh->getDimension();
    }
    // consider making an element since that makes getting the node values a bit easier?
    void inElement(apf::MeshElement * me)
    {
      ment = apf::getMeshEntity(me);
      lngth_o = apf::measure(me);
      lngth = calcDeformedLength(u,ment);
      nedofs = apf::getElementNumbers(num,ment,dofs);
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
    void outElement()
    {
      assembleElementalSystem(k,f,&es,dofs);
    }
    const ElementalSystem * getElementalSystem() {return &es;}
  };
}
#endif
