#ifndef BIO_DMDNI_H_
#define BIO_DMDNI_H_

#include <apf.h>
#include <apfDynamicMatrix.h>
#include <apfField.h>

namespace bio
{
  /**
   * An integrator which gives both the measure and the differential measure of
   *  an Element.
   * The measure is simply
   * \f$ M = \int_{\Omega^e} \mathbf{d}x \mathbf{d}y \mathbf{d}z
   *       = \int_{\Omega^e} \mathbf{det}(J) \mathbf{d}\xi \mathbf{d}\eta \mathbf{d}\gamma
   *       = \sum_{i=1}^{n_{int}} w \mathbf{det}(J) \f$.
   * The differential measure is given by
   * \f$ \frac{\mathbf{d}M}{\mathbf{d}N_i}  = \int_{\Omega^e} \frac{\mathbf{det}(J)}{\mathbf{d}N_i}
   *                                          \mathbf{d}\xi \mathbf{d}\eta \mathbf{d}\gamma \f$
   */
  class dMdNi : public apf::Integrator
  {
  protected:
    double m;
    apf::DynamicMatrix dM_dNi;

    apf::MeshElement * ce;
    int nends;
    int dim;
    apf::Element * e;
    apf::Field * f;
  public:
    dMdNi(apf::Field * fi, int order) 
      : apf::Integrator(order)
      , m()
      , dM_dNi()
      , ce(NULL)
      , nends(0)
      , dim(0)
      , e(NULL)
      , f(fi)
    { }

    /**
     * Set element-dependent variables appropriately.
     */
    void inElement(apf::MeshElement * me)
    {
      m = 0.0;
      
      dim = apf::getDimension(me);
      ce = me;
      e = apf::createElement(f,me);
      nends = apf::countNodes(e);

      dM_dNi.setSize(dim,nends);
    }

    /**
     * Calculate the local contribution to the measure and the differential
     *  measure for the current element.
     */
    void atPoint(apf::Vector3 const& p, double w, double dV)
    {
      apf::Matrix3x3 J;

      apf::NewArray<apf::Vector3> dNi_dxi;
      apf::getShapeGrads(e,p,dNi_dxi);

      apf::DynamicMatrix ddetJ_dNi(dim,nends);
      ddetJ_dNi.zero();
      
      for(int n = 0; n < nends; n++)
      {
	for(int ii = 0; ii < dim; ii++) // x y z
	{
	  apf::getJacobian(ce,p,J);
	  for(int jj = 0; jj < dim; jj++) // xi eta gamma
	    J[ii][jj] = dNi_dxi[n][jj];
	  ddetJ_dNi(n,ii) = apf::getJacobianDeterminant(J,dim);
	}
      }

      // numerical integration
      ddetJ_dNi *= w;
      dM_dNi += ddetJ_dNi;
      m += w * dV;
    }

    double getMeasure() { return m; }
    const apf::DynamicMatrix & getdMdNi() { return dM_dNi; }
  };
}

#endif
