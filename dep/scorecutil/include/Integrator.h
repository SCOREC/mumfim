#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include "MeshSim.h"
#include "Mapping.h"

class Integrator
{  
  public :
    virtual int nbIntegrationPoints(int order) const = 0;
    virtual void iPoint(int i, int order , double &u, double &v, double &w, double &weight) const = 0;  
};

class GaussIntegrator: public Integrator
{
  pEntity ent;
  TRELLIS_MTYPE entityType;

 public:
  GaussIntegrator(pEntity);
  GaussIntegrator(Mapping*);
  virtual int nbIntegrationPoints(int order) const;
  virtual void iPoint(int i,int order,double &u, double &v, double &w, double &weight) const;  
};

#endif


