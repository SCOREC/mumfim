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
  bio::micro_fo_solver slvr;
  slvr.data[bio::MICRO_SOLVER_EPS] = cases[0].ss.cnvgTolerance;
  slvr.data[bio::PREV_ITER_FACTOR] = cases[0].ss.oscPrms.prevNormFactor;
  bio::micro_fo_int_solver slvr_int;
  slvr_int.data[bio::MAX_MICRO_CUT_ATTEMPT] =
      cases[0].ss.oscPrms.maxMicroCutAttempts;
  slvr_int.data[bio::MICRO_ATTEMPT_CUT_FACTOR] =
      cases[0].ss.oscPrms.microAttemptCutFactor;
  slvr_int.data[bio::MAX_MICRO_ITERS] = cases[0].ss.oscPrms.maxIterations;
  slvr_int.data[bio::DETECT_OSCILLATION_TYPE] =
      static_cast<int>(cases[0].ss.oscPrms.oscType);
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
  bio::FiberRVEAnalysis an(fn, vecs, slvr, slvr_int);
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
