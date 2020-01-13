#include <amsiNeumannIntegrators.h>
#include <amsiLinearElasticConstitutive.h>
#include <apfFunctions.h>
#include <sim.h>
#include "bioLinearTissue.h"
namespace bio
{
  LinearTissue::LinearTissue(pGModel imdl, pParMesh imsh, pACase ipd, pACase iss, MPI_Comm cm)
    : FEA(cm)
    , constitutives()
    , amsi::apfSimFEA(imdl,imsh,ipd,iss,cm)
  {
    apf_primary_field = apf::createLagrangeField(apf_mesh,"linear_displacement",apf::VECTOR,1);
    apf_primary_numbering = apf::createNumbering(apf_primary_field);
    GRIter ri = GM_regionIter(imdl);
    pGEntity rgn = NULL;
    while((rgn = (pGEntity)GRIter_next(ri)))
    {
      pAttribute mm = GEN_attrib(rgn,"material model");
      pAttribute cm = Attribute_childByType(mm,"continuum model");
      // should check to make sure the continuum model is iso lin elastic
      pAttributeTensor0 yngs = (pAttributeTensor0)Attribute_childByType(cm,"youngs modulus");
      pAttributeTensor0 psn = (pAttributeTensor0)Attribute_childByType(cm,"poisson ratio");
      double E = AttributeTensor0_value(yngs);
      double v = AttributeTensor0_value(psn);
      constitutives[rgn] = new amsi::LinearElasticIntegrator(apf_primary_field,1,E,v);
    }
    GRIter_delete(ri);
    int dir_tps[] = {amsi::FieldUnit::displacement};
    amsi::buildSimBCQueries(ipd,amsi::BCType::dirichlet,&dir_tps[0],(&dir_tps[0])+1,std::back_inserter(dir_bcs));
    int neu_tps[] = {amsi::NeuBCType::traction,amsi::NeuBCType::pressure};
    amsi::buildSimBCQueries(ipd,amsi::BCType::neumann,&neu_tps[0],(&neu_tps[0])+2,std::back_inserter(neu_bcs));
    /*
      int fei_tps[] = {amsi::ISO_LIN_ELASTIC};
      amsi::buildSimFiniteElementIntegrators(pd,&feis[0],(&feis[0])+1),std::back_inserter(feis));
    */
  }
  LinearTissue::~LinearTissue()
  {
    // clean up the elemental systems stored in constitutives
    for(auto it=constitutives.begin(); it!=constitutives.end(); ++it)
    {
      delete it->second;
    }
  }
  void LinearTissue::UpdateDOFs(const double * sl)
  {
    amsi::WriteOp wrop;
    amsi::FreeApplyOp frop(apf_primary_numbering,&wrop);
    amsi::ApplyVector(apf_primary_numbering,apf_primary_field,sl,first_local_dof,&frop).run();
    apf::synchronize(apf_primary_field);
  }
  void LinearTissue::Assemble(amsi::LAS* las)
  {
    ApplyBC_Neumann(las);
    apf::MeshEntity* me = NULL;
    auto it = apf_mesh->begin(analysis_dim);
    // FIXME shouldn't we skip non-owned elements similar to
    // apf::integrator::process(apf::Mesh*)?
    while ((me = apf_mesh->iterate(it))) {
      amsi::ElementalSystem* constitutive =
          constitutives[R_whatIn((pRegion)me)];
      apf::MeshElement* mlmt = apf::createMeshElement(apf_mesh, me);
      apf::Element * elm = apf::createElement(constitutive->getField(), mlmt);
      constitutive->process(mlmt);
      apf::NewArray<apf::Vector3> dofs;
      apf::getVectorNodes(elm, dofs);
      apf::NewArray<int> ids;
      apf::getElementNumbers(apf_primary_numbering, me, ids);
      AssembleDOFs(las, constitutive->numElementalDOFs(), &ids[0], &dofs[0],
                   &constitutive->getKe()(0, 0), &constitutive->getfe()(0),
                   constitutive->includesBodyForces());
      apf::destroyElement(elm);
      apf::destroyMeshElement(mlmt);
    }
    apf_mesh->end(it);
  }
  
}
