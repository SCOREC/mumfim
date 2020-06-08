#ifndef BIO_BATCHED_MESH_H__
#define BIO_BATCHED_MESH_H__
#include <apf.h>           // for extractCoordinateArray
#include <apfConvert.h>    // for extractCoordinateArray
#include <apfMesh2.h>      // for extractCoordinateArray
#include <apfNumbering.h>  // for getNumbering (getFixedDOF)
#include <Kokkos_Core.hpp>
#include <iostream>
#include "bioMassIntegrator.h"
#include "bioMicroFOParams.h"
#include "bioRVE.h"
#include "bioUtility.h"
namespace bio
{
  // functions that need to interact with the apf mesh
  template <typename Scalar, typename Ordinal, typename HostSpace>
  struct BatchedApfMeshFunctions
  {
    static_assert(
        Kokkos::SpaceAccessibility<HostSpace, Kokkos::HostSpace>::accessible,
        "mesh operations must have access to the host memory space");
    using RWOV = Kokkos::View<Ordinal *, HostSpace>;
    using RWSV = Kokkos::View<Scalar *, HostSpace>;
    // using DeformationGradientView = Kokkos::View<Scalar[3][3], HostSpace>;
    static void getConnectivity(apf::Mesh * m,
                                RWOV connectivity_view,
                                int cellDim)
    {
      bool first = true;
      auto nelem = m->count(cellDim);
      apf::GlobalNumbering * global =
          apf::makeGlobal(apf::numberOwnedNodes(m, "apf_destruct"));
      apf::synchronize(global);
      apf::MeshIterator * it = m->begin(cellDim);
      apf::MeshEntity * e;
      int i = 0;
      while ((e = m->iterate(it)))
      {
        auto etype = m->getType(e);
        apf::Downward verts;
        int nverts = m->getDownward(e, 0, verts);
        if (first)
        {
          first = false;
          if (connectivity_view.extent(0) != nelem * nverts)
          {
            std::cerr << "The connectivity_view has "
                      << connectivity_view.extent(0)
                      << "elements, but it should have " << nelem * nverts
                      << "elements." << std::endl;
            std::exit(EXIT_FAILURE);
          }
        }
        for (int j = 0; j < nverts; ++j)
        {
          connectivity_view(i++) =
              apf::getNumber(global, apf::Node(verts[j], 0));
        }
      }
      m->end(it);
      apf::destroyGlobalNumbering(global);
    }
    template <int NumValues = 3>
    static void getFieldValues(apf::Mesh * m,
                               apf::Field * field,
                               RWSV field_view)
    {
      apf::MeshIterator * it = m->begin(0);
      std::array<double, NumValues> values;
      int i = 0;
      while (apf::MeshEntity * v = m->iterate(it))
      {
        if (m->isOwned(v))
        {
          apf::getComponents(field, v, 0, values.data());
          for (int v = 0; v < NumValues; ++v)
          {
            field_view(i * NumValues + v) = values[v];
          }
          i++;
        }
      }
      m->end(it);
    }
    // TODO this can be made more efficient by having the integrator load
    // direcetly into the Kokkos::View, but we will hold off on this
    // because ultimately we need device side integration routines
    static void getNodalMass(apf::Mesh * mesh,
                             double fiber_density,
                             double fiber_area,
                             RWSV mass_matrix_view)
    {
      auto nodal_mass =
          apf::createLagrangeField(mesh, "nodalMass_tmp", apf::SCALAR, 1);
      apf::zeroField(nodal_mass);
      // build the mass integrator and load the nodal mass into a temporary
      // field
      MassIntegrator massInt(nodal_mass, fiber_density, fiber_area, 3,
                             MassLumpType::RowSum);
      massInt.process(mesh, 1);
      // copy the field values to the kokkos view
      getFieldValues<1>(mesh, nodal_mass, mass_matrix_view);
      apf::destroyField(nodal_mass);
    }
    static void getFixedDof(apf::Numbering * numbering,
                            RWOV fixed_dof_view,
                            std::vector<apf::MeshEntity *> bnd_nds)
    {
      for (std::size_t i = 0; i < bnd_nds.size(); ++i)
      {
        apf::MeshEntity * nd = bnd_nds[i];
        fixed_dof_view(i * 3) = apf::getNumber(numbering, nd, 0, 0);
        fixed_dof_view(i * 3 + 1) = apf::getNumber(numbering, nd, 0, 1);
        fixed_dof_view(i * 3 + 2) = apf::getNumber(numbering, nd, 0, 2);
      }
    }
    // FIXME this is just a hack for now...assuming  naive ordering
    static void getFixedVert(apf::Numbering * numbering,
                            RWOV fixed_vert_view,
                            std::vector<apf::MeshEntity *> bnd_nds)
    {
      for (std::size_t i = 0; i < bnd_nds.size(); ++i)
      {
        apf::MeshEntity * nd = bnd_nds[i];
        fixed_vert_view(i ) = apf::getNumber(numbering, nd, 0, 0)/3;
      }
    }
    // runtime needs
    template <typename Layout, typename Device>
    static void updateRVECoords(RVE & rve,
                                Kokkos::View<Scalar[3][3], Layout, Device>
                                    incremental_deformation_gradient)
    {
      apf::Matrix3x3 F;
      for (int ei = 0; ei < 3; ++ei)
        for (int ej = 0; ej < 3; ++ej)
          F[ei][ej] = incremental_deformation_gradient(ei, ej);
      ApplyIncrementalDeformationGradient(F, rve.getMesh(), rve.getdUField(),
                                          rve.getUField());
    }
  };
}  // namespace bio
#endif
