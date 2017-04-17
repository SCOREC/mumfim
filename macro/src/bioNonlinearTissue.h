#ifndef BIO_NONLINEAR_TISSUE_H_
#define BIO_NONLINEAR_TISSUE_H_
#include "bioLinearTissue.h"
#include "MicroFOMultiscaleTypes.h"
#include "RepresentVolElem.h"       // should be able to take this out... (needed for RVE_Info struct)
#include "bioVolumeConstraint.h"
#include "bioStiffnessVariation.h"
#include <apfFEA.h>
#include <Solvers.h>
#include <NonLinElasticity.h>
#include <amsiMultiscale.h>
#include <amsiAnalysis.h>
#include <amsiUtil.h>
#include <apfsimWrapper.h>
#include <vector>
#include <list>
#include <string>
#include <iostream>
namespace bio
{
  class NonlinearTissue : public amsi::apfSimFEA
  {
  protected:
    std::map<pGEntity,amsi::ElementalSystem*> constitutives;
    std::vector<StiffnessVariation*> stf_vrtn_cnst;
    std::vector<VolumeConstraint*> vol_cnst;
    apf::Field * delta_u;
    apf::Field * strs;
    apf::Field * rcvrd_strs;
    apf::Field * strn;
    apf::Field * previous_rve;
    apf::Field * stf_vrtn;
    apf::Field * axl_yngs_mod;
    double dv_prev;
    int load_step;
    int iteration;
  public:
    NonlinearTissue(pGModel imdl,
                 pParMesh imsh,
                 pACase pd,
                 MPI_Comm cm = AMSI_COMM_SCALE);
    virtual ~NonlinearTissue();
    void computeInitGuess(amsi::LAS * las);
    virtual void ApplyBC_Dirichlet();
    void getLoadOn(pGEntity ent, double * frc);
    void step();
    void iter();
    virtual void Assemble(amsi::LAS * las);
    virtual void UpdateDOFs(const double * );
    void UpdateLambda();
    void UpdateBeta(double);
    void setBeta(double);
    void setLambda(double);
    void computeDispL2Norm(double &);
    void recoverSecondaryVariables(int);
    void storeStress(apf::MeshElement * me, double * stress);
    void storeStrain(apf::MeshElement * me, double * strain);
    int getIteration(){return iteration;}
    apf::Numbering * getNumbering() { return apf_primary_numbering; }
    apf::Field * getdUField() { return delta_u; }
    apf::Field * getUField() {return apf_primary_field;}
    //void logCnstrntParams(int ldstp, int iteration, int rnk);
  };
}
#endif
