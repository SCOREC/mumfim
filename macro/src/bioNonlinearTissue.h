#ifndef BIO_NONLINEAR_TISSUE_H_
#define BIO_NONLINEAR_TISSUE_H_
#include "bioLinearTissue.h"
#include "bioStiffnessVariation.h"
#include "bioVolumeConstraint.h"
// micro_fo
#include <bioMicroFOMultiscale.h>
#include <bioRVE2.h>
// amsi
#include <apfFEA.h>
#include <Solvers.h>
#include <NonLinElasticity.h>
#include <Solvers.h>
#include <amsiAnalysis.h>
#include <amsiMultiscale.h>
#include <amsiUtil.h>
#include <apfFEA.h>
#include <apfsimWrapper.h>
#include <iostream>
#include <list>
#include <string>
#include <vector>
namespace bio
{
  class CurrentCoordFunc;
  class NonlinearTissue : public amsi::apfSimFEA
  {
    protected:
    std::map<pGEntity, amsi::ElementalSystem*> constitutives;
    std::vector<StiffnessVariation*> stf_vrtn_cnst;
    std::vector<VolumeConstraint*> vol_cnst;
    apf::Field * delta_u;
    apf::Field * current_coords;  // coordinates in current config
    apf::Field * strs;
    apf::Field * rcvrd_strs;
    apf::Field * strn;
    apf::Field * dfm_grd;
    apf::Field * previous_rve;
    apf::Field * stf_vrtn;
    apf::Field * axl_yngs_mod;
    double dv_prev;
    int load_step;
    int iteration;
    public:
    NonlinearTissue(pGModel imdl, pParMesh imsh, pACase pd,
                    MPI_Comm cm = AMSI_COMM_SCALE);
    virtual ~NonlinearTissue();
    void computeInitGuess(amsi::LAS* las);
    virtual void ApplyBC_Dirichlet();
    void getLoadOn(pGEntity ent, double* frc);
    void step();
    void iter();
    virtual void Assemble(amsi::LAS* las);
    virtual void UpdateDOFs(const double*);
    void UpdateLambda();
    void UpdateBeta(double);
    void setBeta(double);
    void setLambda(double);
    void computeDispL2Norm(double&);
    void recoverSecondaryVariables(int);
    void storeStress(apf::MeshElement* me, double* stress);
    void storeStress(apf::MeshElement* me, apf::Matrix3x3 eps);
    void storeStrain(apf::MeshElement* me, double* strain);
    void storeStrain(apf::MeshElement* me, apf::Matrix3x3 eps);
    int getIteration() { return iteration; }
    apf::Numbering* getNumbering() { return apf_primary_numbering; }
    apf::Field* getdUField() { return delta_u; }
    apf::Field* getUField() { return apf_primary_field; }
    apf::Mesh* getMesh() { return apf_mesh; }
    CurrentCoordFunc* currentCoordFunc;
    // void logCnstrntParams(int ldstp, int iteration, int rnk);
    // function to get mapping by summing reference coordinate with
    // displacements
  };
  class CurrentCoordFunc : public apf::Function {
  private:
    apf::Field * Xf;
    apf::Field * Uf;
  public:
    CurrentCoordFunc(apf::Field * xf, apf::Field * uf) : Xf(xf), Uf(uf) {}
    void eval(apf::MeshEntity * e, double * result)
    {
      // make sure that we are only evaluating on vertices
      assert(apf::getMesh(Xf)->getType(e) == apf::Mesh::VERTEX);
      apf::Vector3 X, U;
      // get the displacements
      apf::getVector(Uf, e, 0, U);
      // get the reference coordinates
      apf::getVector(Xf, e, 0, X);
      apf::Vector3 * r = (apf::Vector3*)result;
      // set the current coordinates to be the reference plus the displacements
      *r = X + U;
    }
  };
}
#endif
