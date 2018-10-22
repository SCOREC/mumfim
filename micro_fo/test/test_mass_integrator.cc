#include <amsiAnalysis.h>
#include <apf.h>
#include <apfMesh2.h>
#include <las.h>
#include <lasCSRCore.h>
#include <lasConfig.h>
#include <cassert>
#include "bioFiberNetworkIO.h"
#include "bioMassIntegrator.h"
int main(int argc, char * argv[])
{
  amsi::initAnalysis(argc, argv, MPI_COMM_WORLD);
  // create empty mesh
  // insert edge
  std::stringstream ss;
  ss << "3 2 0\n"
     << "0 0 0\n"
     << "1 1 1\n"
     << "1 2.5 0\n"
     << "0 1\n"
     << "1 2\n";
  apf::Mesh2 * mesh = bio::loadFromStream(ss);
  // apf::Mesh2 * mesh = bio::loadFromFile(
  //    "/fasttmp/mersoj/develop/biotissue_problems/fiber_networks/"
  //    "density_100_networks/del_rho100_new_1.txt");
  std::cout << "Mesh has " << mesh->count(1) << " edges" << std::endl;
  // add lagrange shape functions for "primary" field
  apf::Field * disp_fld =
      apf::createLagrangeField(mesh, "disp", apf::VECTOR, 1);
  // add step shape functions for density
  // apf::Field * dens_fld = apf::createStepField(mesh, "density", apf::SCALAR);
  apf::Field * dens_fld =
      apf::createField(mesh, "density", apf::SCALAR, apf::getConstant(1));
  // set dummy values for "primary" field
  apf::MeshEntity * me = NULL;
  for (int dim = 1; dim >= 0; --dim)
  {
    apf::MeshIterator * itr = mesh->begin(dim);
    while ((me = mesh->iterate(itr)))
    {
      // if(mesh->isOwned(me)) {
      if (true)
      {
        if (apf::getShape(disp_fld)->hasNodesIn(dim))
        {
          int nnds = apf::getShape(disp_fld)->countNodesOn(mesh->getType(me));
          for (int nd = 0; nd < nnds; ++nd)
          {
            apf::setVector(disp_fld, me, nd, apf::Vector3(0, 0, 0));
          }
        }
        // set density values
        if (apf::getShape(dens_fld)->hasNodesIn(dim))
        {
          int nnds = apf::getShape(dens_fld)->countNodesOn(mesh->getType(me));
          for (int nd = 0; nd < nnds; ++nd)
          {
            apf::setScalar(dens_fld, me, nd, 0.5);
          }
        }
      }
    }
    mesh->end(itr);
  }
  apf::Numbering * num = apf::createNumbering(disp_fld);
  int ndofs = apf::NaiveOrder(num);
  las::Sparsity * sprs = las::createCSR(num, ndofs);
  las::Mat * mass = las::createCSRMatrix(sprs);
  apf::Integrator * massIntegrator =
      new bio::MassIntegrator(num, dens_fld, disp_fld, mass, 2, bio::MassLumpType::None);
  massIntegrator->process(mesh, 1);
  // process(massIntegrator, mesh, 1);
  las::printSparskitMat(std::cout, mass, las::PrintType::full);
  // run integrate on mesh
  // construct matrix with "correct" values
  // compare matrices
  /*
  bool close = las::SparskitMatClose(baselineMat, bioIntegratorMat);
  assert(close);
  */
  delete massIntegrator;
  las::LasCreateMat * mb = las::getMatBuilder<las::sparskit>(0);
  mb->destroy(mass);
  delete mb;
  mb = NULL;
  las::destroySparsity<las::CSR *>(sprs);
  apf::destroyNumbering(num);
  amsi::freeAnalysis();
  return 0;
}
