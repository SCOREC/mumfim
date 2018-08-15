#include <amsiAnalysis.h>
#include <amsiDetectOscillation.h>
#include <apf.h>
#include <lasCSRCore.h>
#include <mpi.h>
#include <iostream>
#include "bioFiberNetworkIO.h"
#include "bioFiberRVEAnalysis.h"
#include "bioMultiscaleRVEAnalysis.h"
#include "bioVerbosity.h"
int main(int argc, char * argv[])
{
  amsi::initAnalysis(argc, argv, MPI_COMM_WORLD);
  std::string file_name =
      "/fasttmp/mersoj/develop/biotissue_problems/fiber_networks/"
      "density_100_networks/del_rho100_new_10.txt";
  double F[9] = {1.5, 0, 0, 0, 1 / (sqrt(1.5)), 0, 0, 0, 1 / (sqrt(1.5))};
  bio::micro_fo_solver slvr;
  slvr.data[bio::MICRO_SOLVER_EPS] = 1E-6;
  slvr.data[bio::PREV_ITER_FACTOR] = 1.0;
  bio::micro_fo_int_solver slvr_int;
  slvr_int.data[bio::MAX_MICRO_CUT_ATTEMPT] = 1;
  slvr_int.data[bio::MICRO_ATTEMPT_CUT_FACTOR] = 2.0;
  slvr_int.data[bio::MAX_MICRO_ITERS] = 20;
  slvr_int.data[bio::DETECT_OSCILLATION_TYPE] =
      static_cast<int>(amsi::DetectOscillationType::IterationOnly);
  bio::micro_fo_data deformation_gradient;
  for (int i = 0; i < 9; ++i)
  {
    deformation_gradient.data[i] = F[i];
  }
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::stringstream prm_name_ss;
  prm_name_ss << file_name << ".params";
  bio::FiberNetworkReactions rctns;
  apf::Mesh2 * fn_msh = bio::loadFromFile(file_name);
  bio::loadParamsFromFile(
      fn_msh, prm_name_ss.str(), std::back_inserter(rctns.rctns));
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
  bool result = an.run(deformation_gradient);
  if (!result)
  {
    std::cerr << "The microscale analysis failed to converge" << std::endl;
  }
  std::stringstream sout;
  sout << "rnk_" << rank << "_fn_" << an.getFn()->getRVEType();
  apf::writeVtkFiles(sout.str().c_str(), an.getFn()->getNetworkMesh(), 1);
  // las::destroySparsity<CSR *>(sprs);
  amsi::freeAnalysis();
  return 1;
}
