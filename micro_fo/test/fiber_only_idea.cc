#include <amsiAnalysis.h>
#include <amsiDetectOscillation.h>
#include <apf.h>
#include <lasCSRCore.h>
#include <lasCorePETSc.h>
#include <mpi.h>
#include <iostream>
#include "bioFiberNetworkIO.h"
#include "bioFiberRVEAnalysis.h"
#include "bioMicroFOConfig.h"
#include "bioMicroFOParams.h"
#include "bioMultiscaleRVEAnalysis.h"
#include "bioVerbosity.h"
// we are requireing kokkos as a dependency
#include <Kokkos_Core.hpp>
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
  // ideally this replaces amsiInitAnalysis
  amsi::ScopeGuard(argc, argv, MPI_COMM_WORLD);
  // ideally this replaces las initialization including initPETScLAS
  las::ScopeGuard(&argc, &argv, MPI_COMM_WORLD);
  // we will require kokkos as a dependency
  Kokkos::ScopeGaurd(argc, argv);
  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " job.yaml" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  // load the problem settings
  std::vector<bio::MicroCase> cases;
  bio::loadMicroFOFromYamlFile(argv[1], cases);
  bio::printMicroFOCase(cases[0]);
  std::string file_name = cases[0].pd.meshFile;
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::stringstream prm_name_ss;
  prm_name_ss << file_name << ".params";
  // load the mesh
  apf::Mesh2 * fn_msh = bio::loadFromFile(file_name);
  // load the fiber network reactions
  bio::FiberNetworkReactions rctns;
  bio::loadParamsFromFile(fn_msh, prm_name_ss.str(),
                          std::back_inserter(rctns.rctns));
  bio::FiberNetwork * fn = new bio::FiberNetwork(fn_msh);
  fn->setFiberReactions(rctns.rctns);
  std::cout << "Problem has " << ndofs << " degrees of freedom" << std::endl;
  an = bio::createFiberRVEAnalysis(fn, vecs, *cases[0].ss)
           an->computeStiffnessMatrix();
  double stress[6];
  double C[36];
  apf::Matrix3x3 strss;
  Kokkos::Timer timer;
  double time = timer.seconds();
  timer.reset();
  bool result = an->run(cases[0].pd.deformationGradient, stress);
  an->computeMaterialStiffness(C);
  std::cout << "Stress from run" << std::endl;
  stressToMat(stress, strss);
  std::cout << strss << std::endl;
  for (int i = 0; i < 6; ++i)
  {
    for (int j = 0; j < 6; ++j)
    {
      std::cout << C[i * 6 + j] << " ";
    }
    std::cout << "\n";
  }
  std::cout << std::endl;
  double time2 = timer.seconds();
  std::cout << "Took: " << time2 << " seconds." << std::endl;
  if (!result)
  {
    std::cerr << "The microscale analysis failed to converge" << std::endl;
  }

  // this replaces us directly writing any vtk files and gives us a chance
  // to write all data of interest to disk which will make our lives easier
  // w.r.t. checkpointing
  // an->writeToStream(stream);
  an->checkpoint();
  return 0;
}
