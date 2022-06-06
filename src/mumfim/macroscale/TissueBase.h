#ifndef MUMFIM_MACROSCALE_TISSUEBASE_H
#define MUMFIM_MACROSCALE_TISSUEBASE_H
#include <apfFEA.h>
#include <memory>
#include <unordered_map>
namespace mumfim
{
  class TissueBase : public amsi::apfFEA
  {
    protected:
    TissueBase(apf::Mesh * mesh,
               const mt::CategoryNode & analysis_case,
               std::vector<amsi::DirichletBCEntry> dbc,
               std::vector<amsi::NeumannBCEntry> nbc,
               const std::string & analysis_name = "",
               MPI_Comm cm = AMSI_COMM_SCALE)
        : amsi::apfFEA(mesh,
                       analysis_case,
                       std::move(dbc),
                       std::move(nbc),
                       analysis_name,
                       cm)
    {
    }
    TissueBase(apf::Mesh * mesh,
               const amsi::ModelDefinition & problem_definition,
               const amsi::ModelDefinition & solution_strategy,
               const amsi::ModelDefinition & output,
               std::vector<amsi::DirichletBCEntry> dbc,
               std::vector<amsi::NeumannBCEntry> nbc,
               const std::string & analysis_name = "",
               MPI_Comm cm = AMSI_COMM_SCALE)
        : amsi::apfFEA(mesh,
                       problem_definition,
                       solution_strategy,
                       output,
                       dbc,
                       nbc,
                       analysis_name,
                       cm)
    {
    }
    template <typename T>
    void AssembleIntegratorIntoLAS(amsi::LAS * las,
                                   T integrator,
                                   apf::Field * coordinates = nullptr)
    {
      static_assert(std::is_invocable_r_v<amsi::ElementalSystem *, T,
                                          apf::MeshEntity *, int>);
      if (coordinates == nullptr)
      {
        coordinates = apf_mesh->getCoordinateField();
      }
      apf::MeshIterator * it = apf_mesh->begin(analysis_dim);
      apf::MeshEntity * me = nullptr;
      while ((me = apf_mesh->iterate(it)))
      {
        if (!apf_mesh->isOwned(me))
        {
          continue;
        }
        apf::MeshElement * mlm = apf::createMeshElement(coordinates, me);
        auto * sys = std::invoke(integrator, me, 0);
        apf::Element * elm = apf::createElement(sys->getField(), mlm);
        sys->process(mlm);
        apf::NewArray<apf::Vector3> dofs;
        apf::getVectorNodes(elm, dofs);
        apf::NewArray<int> ids;
        apf::getElementNumbers(apf_primary_numbering, me, ids);
        AssembleDOFs(las, sys->numElementalDOFs(), &ids[0], &dofs[0],
                                &sys->getKe()(0, 0), &sys->getfe()(0),
                                sys->includesBodyForces());
        apf::destroyElement(elm);
        apf::destroyMeshElement(mlm);
      }
      apf_mesh->end(it);
    }
  };
}  // namespace mumfim
#endif  // MUMFIM_SRC_MUMFIM_MACROSCALE_MULTISCALETISSUE_CC_TISSUEBASE_H
