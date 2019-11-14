#include <amsiAnalysis.h>
#include <amsiDetectOscillation.h>
#include <apf.h>
#include <lasCSRCore.h>
#include <lasCorePETSc.h>
#include <mpi.h>
#include <iostream>
#include "bioFiberNetworkIO.h"
#include "bioFiberRVEAnalysis.h"
#include "bioMicroFOParams.h"
#include "bioMultiscaleRVEAnalysis.h"
#include "bioVerbosity.h"
#include "bioMicroFOConfig.h"
#ifdef ENABLE_KOKKOS
#include <Kokkos_Core.hpp>
#endif
void stressToMat(double * stress_arr, apf::Matrix3x3 & stress)
{
  stress[0][0] = stress_arr[0];
  stress[1][1] = stress_arr[1];
  stress[2][2] = stress_arr[2];
  stress[1][2] = stress_arr[3];
  stress[2][1] = stress_arr[3];
  stress[0][2] = stress_arr[4];
  stress[2][0] = stress_arr[4];
  stress[0][1] = stress_arr[5];
  stress[1][0] = stress_arr[5];
}
int main(int argc, char * argv[])
{
  amsi::initAnalysis(argc, argv, MPI_COMM_WORLD);
#ifdef MICRO_USING_PETSC
  las::initPETScLAS(&argc, &argv, MPI_COMM_WORLD);
#endif
#ifdef ENABLE_KOKKOS
  Kokkos::initialize(argc, argv);
  {
#endif
  if(argc != 2)
  {
    std::cerr<<"Usage: "<<argv[0]<<" job.yaml"<<std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  std::vector<bio::MicroCase> cases;
  bio::loadMicroFOFromYamlFile(argv[1], cases);
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
  las::Sparsity * massSparsity = NULL;
#if defined MICRO_USING_SPARSKIT
  las::Sparsity * sprs = las::createCSR(n, ndofs);
  massSparsity = las::createIdentityCSR(ndofs);
#elif defined MICRO_USING_PETSC
  las::Sparsity * sprs = las::createPetscSparsity(n, ndofs, MPI_COMM_SELF);
#endif
  // clean up the un-needed field and numbering
  apf::destroyField(u);
  apf::destroyNumbering(n);
#if defined MICRO_USING_SPARSKIT
  las::SparskitBuffers * bfrs = new las::SparskitBuffers(ndofs);
#elif defined MICRO_USING_PETSC
  void * bfrs = NULL;
#endif
  bio::FiberNetwork * fn = new bio::FiberNetwork(fn_msh);
  fn->setFiberReactions(rctns.rctns);
  bio::LinearStructs<las::MICRO_BACKEND> * vecs =
      bio::createLinearStructs(ndofs, cases[0].ss->slvrTolerance, sprs, bfrs, massSparsity);
  bio::FiberRVEAnalysis * an = NULL;
  if(cases[0].ss->slvrType == bio::SolverType::Implicit) {
    an = bio::createFiberRVEAnalysis(
        fn, vecs, *cases[0].ss, bio::FiberRVEAnalysisType::StaticImplicit);
  }
  else if (cases[0].ss->slvrType == bio::SolverType::Explicit) {
    an = bio::createFiberRVEAnalysis(
        fn, vecs, *cases[0].ss, bio::FiberRVEAnalysisType::Explicit);
  }
  an->computeStiffnessMatrix();
  double stress[6];
  apf::Matrix3x3 strss;
#ifdef ENABLE_KOKKOS
  Kokkos::Timer timer;
  double time = timer.seconds();
  timer.reset();
#endif
  bool result = an->run(cases[0].pd.deformationGradient, stress);
  std::cout<<"Stress from run"<<std::endl;
  stressToMat(stress, strss);
  std::cout<<strss<<std::endl;
#ifdef ENABLE_KOKKOS
  double time2 = timer.seconds();
  std::cout<<"Took: "<<time2 << " seconds."<<std::endl;
#endif
  if (!result)
  {
    std::cerr << "The microscale analysis failed to converge" << std::endl;
  }
  std::stringstream sout;
  sout << "rnk_" << rank << "_fn_" << an->getFn()->getRVEType();
  apf::writeVtkFiles(sout.str().c_str(), an->getFn()->getNetworkMesh(), 1);
  las::destroySparsity<las::MICRO_BACKEND>(sprs);
#ifdef ENABLE_KOKKOS
  }
  Kokkos::finalize();
#endif
#ifdef MICRO_USING_PETSC
  las::finalizePETScLAS();
#endif
  amsi::freeAnalysis();
  return 0;
}
