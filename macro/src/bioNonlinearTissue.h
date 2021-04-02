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
#include <Solvers.h>
#include <amsiAnalysis.h>
#include <amsiMultiscale.h>
#include <amsiUtil.h>
#include <apfFEA.h>
#include <apfFunctions.h>
#include <iostream>
#include <list>
#include <string>
#include <vector>
#include <model_traits/CategoryNode.h>
namespace bio
{
  class CurrentCoordFunc;
  class NonlinearTissue : public amsi::apfFEA
  {
    protected:
    amsi::XpYFunc * xpyfnc;
    std::map<apf::ModelEntity*, std::unique_ptr<amsi::ElementalSystem>> constitutives;
    std::vector<std::unique_ptr<StiffnessVariation>> stf_vrtn_cnst;
    std::vector<std::unique_ptr<VolumeConstraint>> vol_cnst;
    apf::Field * delta_u;
    apf::Field * current_coords;  // coordinates in current config
    apf::Field * strs;
    apf::Field * strn;
    apf::Field * dfm_grd;
    apf::Field * previous_rve;
    apf::Field * stf_vrtn;
    apf::Field * axl_yngs_mod;
    apf::Field * prev_coords;
    amsi::XpYFunc * prv_crd_fnc;
    double dv_prev;
    int load_step;
    int iteration;
    public:
    NonlinearTissue(apf::Mesh* mesh, const mt::CategoryNode& analysis_case,
                    MPI_Comm cm = AMSI_COMM_SCALE);
    virtual ~NonlinearTissue();
    void computeInitGuess(amsi::LAS* las);
    void getLoadOn(apf::ModelEntity* ent, double* frc);
    void step();
    void iter();
    virtual void Assemble(amsi::LAS* las) override;
    virtual void UpdateDOFs(const double*) override;
    virtual void recoverSecondaryVariables(int);
    virtual void preRun() {};
    void storeStress(apf::MeshElement* me, double* stress);
    void storeStress(apf::MeshElement* me, apf::Matrix3x3 eps);
    void storeStrain(apf::MeshElement* me, double* strain);
    void storeStrain(apf::MeshElement* me, apf::Matrix3x3 eps);
    int getIteration() { return iteration; }
    apf::Numbering* getNumbering() { return apf_primary_numbering; }
    //apf::Field* getdUField() { return delta_u; }
    apf::Field* getUField() { return apf_primary_field; }
    apf::Mesh* getMesh() { return apf_mesh; }
    // void logCnstrntParams(int ldstp, int iteration, int rnk);
    // function to get mapping by summing reference coordinate with
    // displacements
  };
}
#endif
