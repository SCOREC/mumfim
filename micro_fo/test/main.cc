#include <amsiAnalysis.h>
#include <amsiDetectOscillation.h>
#include <apf.h>
#include <lasCSRCore.h>
#include <mpi.h>
#include <iostream>
#include "bioFiberNetworkIO.h"
#include "bioFiberRVEAnalysis.h"
#include "bioMicroFOParams.h"
#include "bioMultiscaleRVEAnalysis.h"
#include "bioVerbosity.h"
int main(int argc, char * argv[])
{
  amsi::initAnalysis(argc, argv, MPI_COMM_WORLD);
  std::vector<bio::MicroCase> cases;
  bio::loadMicroFOFromYamlFile(
      "/fasttmp/mersoj/develop/biotissue/micro_fo/test/fiber_only.yaml", cases);
  bio::printMicroFOCase(cases[0]);
  std::string file_name = cases[0].pd.meshFile;
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::stringstream prm_name_ss;
  prm_name_ss << file_name << ".params";
  bio::FiberNetworkReactions rctns;
  apf::Mesh2 * fn_msh = bio::loadFromFile(file_name);
  bio::loadParamsFromFile(fn_msh, prm_name_ss.str(),
                          std::back_inserter(rctns.rctns));
  apf::Field * u = apf::createLagrangeField(fn_msh, "u", apf::VECTOR, 1);
  apf::Numbering * n = apf::createNumbering(u);
  int ndofs = apf::NaiveOrder(n);
  // do we need to zero field? if this assert fails we need to zero the field.
  std::cout << "Problem has " << ndofs << " degrees of freedom" << std::endl;
  assert(ndofs > 0);
  las::Sparsity * sprs = las::createCSR(n, ndofs);
  // clean up the un-needed field and numbering
  apf::destroyField(u);
  apf::destroyNumbering(n);
  las::SparskitBuffers * bfrs = new las::SparskitBuffers(ndofs);
  bio::FiberNetwork * fn = new bio::FiberNetwork(fn_msh);
  fn->setFiberReactions(rctns.rctns);
  bio::LinearStructs * vecs = bio::createLinearStructs(ndofs, sprs, bfrs);
  bio::FiberRVEAnalysis an(fn, vecs, cases[0].ss);
  assert(an.multi == NULL);
  bool result = an.run(cases[0].pd.deformationGradient);
  if (!result)
  {
    std::cerr << "The microscale analysis failed to converge" << std::endl;
  }
  std::stringstream sout;
  sout << "rnk_" << rank << "_fn_" << an.getFn()->getRVEType();
  apf::writeVtkFiles(sout.str().c_str(), an.getFn()->getNetworkMesh(), 1);
  // las::destroySparsity<CSR *>(sprs);
  amsi::freeAnalysis();
  return 0;
}
