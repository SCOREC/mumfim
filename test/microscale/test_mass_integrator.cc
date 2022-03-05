#include <amsiAnalysis.h>
#include <apf.h>
#include <apfMesh2.h>
#include <las.h>
#include <lasCSRCore.h>
#include <lasConfig.h>
#include <cassert>
#include "bioFiberNetworkIO.h"
#include "bioMassIntegrator.h"
#include "bioFiberReactions.h"
#include <vector>
int main(int argc, char * argv[])
{
  amsi::initAnalysis(argc, argv, MPI_COMM_WORLD);
  std::vector<bio::FiberReaction *> rctns;
  // create empty mesh
  // insert edge
  std::stringstream fibers;
  std::stringstream params;
  fibers << "3 2 0\n"
     << "0 0 0\n"
     << "1 1 1\n"
     << "1 2.5 0\n"
     << "0 1\n"
     << "1 2\n";
  params <<"1 2\n"<< "0 0.5 10000 1.0\n" << "0\n" << "0\n";
  apf::Mesh2 * mesh = bio::loadFromStream(fibers);
  bio::loadParamsFromStream(mesh, params, std::back_inserter(rctns));
  std::cout << "Mesh has " << mesh->count(1) << " edges" << std::endl;
  // add lagrange shape functions for "primary" field
  apf::Field * disp_fld =
      apf::createLagrangeField(mesh, "disp", apf::VECTOR, 1);
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
      }
    }
    mesh->end(itr);
  }
  apf::Numbering * num = apf::createNumbering(disp_fld);
  int ndofs = apf::NaiveOrder(num);
  las::Sparsity * sprs = las::createCSR(num, ndofs);
  las::Mat * mass = las::createCSRMatrix(sprs);
  apf::Integrator * massIntegrator =
      new bio::MassIntegrator(num, disp_fld, mass, &(rctns[0]), 2, bio::MassLumpType::RowSum);
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
  // FIXME this test currently is set to fail because we don't actually confirm the value
  // is correct in this code...I hand chacked the values for now.
  return 1;
}
