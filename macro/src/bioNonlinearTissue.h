#ifndef BIO_NONLINEAR_TISSUE_H_
#define BIO_NONLINEAR_TISSUE_H_
#include "bioLinearTissue.h"
#include "bioStiffnessVariation.h"
#include "bioVolumeConstraint.h"
// micro_fo
#include <bioMultiscaleMicroFOParams.h>
#include <bioRVE.h>
// amsi
#include <apfFEA.h>
#include <Solvers.h>
#include <NonLinElasticity.h>
#include <Solvers.h>
#include <amsiAnalysis.h>
#include <amsiMultiscale.h>
#include <amsiUtil.h>
#include <apfFEA.h>
#include <apfFunctions.h>
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
    amsi::XpYFunc * xpyfnc;
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
    NonlinearTissue(pGModel imdl, pParMesh imsh, pACase pd, pACase ss,
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
    virtual void recoverSecondaryVariables(int);
    void storeStress(apf::MeshElement* me, double* stress);
    void storeStress(apf::MeshElement* me, apf::Matrix3x3 eps);
    void storeStrain(apf::MeshElement* me, double* strain);
    void storeStrain(apf::MeshElement* me, apf::Matrix3x3 eps);
    int getIteration() { return iteration; }
    apf::Numbering* getNumbering() { return apf_primary_numbering; }
    apf::Field* getdUField() { return delta_u; }
    apf::Field* getUField() { return apf_primary_field; }
    apf::Mesh* getMesh() { return apf_mesh; }
    // void logCnstrntParams(int ldstp, int iteration, int rnk);
    // function to get mapping by summing reference coordinate with
    // displacements
  };
}
#endif
