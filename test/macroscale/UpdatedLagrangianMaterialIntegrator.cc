#include <apfBox.h>
#include <mumfim/macroscale/NeoHookeanIntegrator.h>
#include <mumfim/macroscale/UpdatedLagrangianMaterialIntegrator.h>
#include <catch2/catch.hpp>
#include "TestSupport.h"
using test::compare_dynamic_matrices;
using test::compare_dynamic_vectors;

TEST_CASE("Compare UL Integrators")
{
  auto * box_mesh = apf::makeMdsBox(10, 10, 10, 3, 3, 3, true);
  auto * displacement =
      apf::createLagrangeField(box_mesh, "displacement", apf::VECTOR, 1);
  apf::zeroField(displacement);
  auto * xpyfnc =
      new amsi::XpYFunc(box_mesh->getCoordinateField(), displacement);
  auto * current_coords =
      apf::createUserField(box_mesh, "current_coordinates", apf::VECTOR,
                           apf::getLagrange(1), xpyfnc);
  auto * strs = apf::createIPField(box_mesh, "stress", apf::MATRIX, 1);
  apf::zeroField(strs);
  auto * strn = apf::createIPField(box_mesh, "strain", apf::MATRIX, 1);
  apf::zeroField(strn);
  auto * dfm_grd = apf::createIPField(box_mesh, "F", apf::MATRIX, 1);
  apf::zeroField(dfm_grd);
  // set the mesh to have some random displacements
  auto * it = box_mesh->begin(0);
  while (auto * entity = box_mesh->iterate(it))
  {
    apf::Vector3 coords;
    apf::getVector(box_mesh->getCoordinateField(), entity, 0, coords);
    apf::Vector3 disp{sin(coords[0]), sin(coords[1]), sin(coords[2])};
    apf::setVector(displacement, entity, 0, disp);
  }
  box_mesh->end(it);
  constexpr double youngs_modulus = 100;
  constexpr double poissons_ratio = 0.33;
  // create Neohookean integrator
  mumfim::NeoHookeanIntegrator neohookean(displacement, dfm_grd, current_coords,
                                          strs, strn, youngs_modulus, poissons_ratio,
                                          1);
  // UpdatedLagrangean with NeohookeanMaterial
  mumfim::UpdatedLagrangianMaterialIntegrator updated_lagrangian(
      [=](apf::Matrix3x3 F, apf::MeshEntity *, int)
      {
        auto shear_modulus = youngs_modulus / (2.0 * (1.0 + poissons_ratio));
        return mumfim::NeohookeanMaterial(shear_modulus, poissons_ratio, F); },
      strn, strs, displacement, dfm_grd, 1);
  // process mesh
  it = box_mesh->begin(3);
  while(apf::MeshEntity * ent = box_mesh->iterate(it)) {
    if (!box_mesh->isOwned(ent))
    {
      continue;
    }
    apf::MeshElement * ccmlm = apf::createMeshElement(current_coords, ent);
    apf::MeshElement * mlm = apf::createMeshElement(box_mesh->getCoordinateField(), ent);
    neohookean.process(mlm);
    updated_lagrangian.process(ccmlm);
    compare_dynamic_matrices(neohookean.getKe(),updated_lagrangian.getKe());
    compare_dynamic_vectors(neohookean.getfe(), updated_lagrangian.getfe());
    break;
  }
  box_mesh->end(it);


  delete xpyfnc;
}