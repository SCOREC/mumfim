#ifndef BIO_VOLUME_TERMS_H_
#define BIO_VOLUME_TERMS_H_
#include <apf.h>
namespace bio
{
  // need vol, det J and d_detJ / dx_rve to convert to macroscale
  class CalcdV_dx_rve : public apf::Integrator
  {
  private:
    int dim;
    int nends;
    double V;
    apf::Mesh * msh;
    apf::DynamicVector dV_dx_rve;
    apf::Field * u;
    apf::MeshEntity * me;
    apf::MeshElement * mlm;
    apf::Element * e;
    apf::FieldShape * s;
    apf::EntityShape * es;
  public:
    CalcdV_dx_rve(int o, apf::Field * du)
      : Integrator(o)
      , dim(-1)
      , nends(0)
      , V(0.0)
      , msh(apf::getMesh(du))
      , dV_dx_rve()
      , u(du)
      , me(NULL)
      , mlm(NULL)
      , e(NULL)
      , s(apf::getShape(du))
      , es(NULL)
    {}
    virtual void inElement(apf::MeshElement * m)
    {
      me = apf::getMeshEntity(m);
      es = s->getEntityShape(msh->getType(me));
      mlm = m;
      dim = apf::getDimension(mlm);
      e = apf::createElement(u,mlm);
      nends = apf::countNodes(e);
      dV_dx_rve.resize(nends*dim);
      dV_dx_rve.zero();
    }
    virtual void atPoint(apf::Vector3 const& p, double w, double dV)
    {
      apf::NewArray<apf::Vector3> dx_dxi;
      // need local grads, not global
      es->getLocalGradients(msh,me,p,dx_dxi);
      V += w * dV;
      double ddett[24] = {0};
      for(int ii = 0; ii < nends; ++ii)
      {
        for(int dd = 0; dd < dim; ++dd)
        {
          apf::Matrix3x3 J;
          apf::getJacobian(mlm,p,J);
          for(int d2 = 0; d2 < dim; ++d2)
            J[d2][dd] = dx_dxi[ii][d2];
          ddett[ii * dim + dd] = apf::getDeterminant(J);
        }
      }
      for(int ii = 0; ii < 24; ++ii)
        dV_dx_rve[ii] += w * ddett[ii];
    }
    virtual void outElement()
    {
      apf::destroyElement(e);
    }
    void getdVdxrve(apf::DynamicVector & dVdxrve)
    {
      dVdxrve = dV_dx_rve;
    }
  };

}
#endif
