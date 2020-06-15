#include <amsiAnalysis.h>
#include <amsiDetectOscillation.h>
#include <apf.h>
#include <mpi.h>
#include <iostream>
#include "bioFiberNetworkIO.h"
#include "bioFiberRVEAnalysis.h"
#include "bioMicroFOParams.h"
#include "bioMultiscaleRVEAnalysis.h"
#include "bioVerbosity.h"
#include <bioFiberRVEAnalysisStaticImplicit.h>
int main(int argc, char * argv[])
{
  amsi::initAnalysis(argc, argv, MPI_COMM_WORLD);
  std::vector<bio::MicroCase> cases;
  bio::loadMicroFOFromYamlFile(
      "./test_global_stiffness_data/global_stiffness.yaml", cases);
  las::LasCreateMat* mb = las::getMatBuilder<las::sparskit>(0);
  for (std::size_t i = 0; i < cases.size(); ++i)
  //for (std::size_t i = 0; i < 1; ++i)
  {
    bio::printMicroFOCase(cases[i]);
    std::string file_name = cases[i].pd.meshFile;
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    bio::FiberNetworkLibrary network_library;
    network_library.load(file_name,file_name+".params",0,0);
    auto fiber_network = network_library.getUniqueCopy(0, 0);
    auto solution_strategy = std::unique_ptr<bio::MicroSolutionStrategy>{new bio::MicroSolutionStrategy};
    // set the solution strategy to give me an implicit run so that I can get
    // the stiffness matrix. All the other parameters don't matter...
    auto osc_prms = bio::DetectOscillationParams();
    osc_prms.maxIterations = 10;
    osc_prms.maxMicroCutAttempts = 2;
    osc_prms.microAttemptCutFactor  = 2;
    osc_prms.oscType = amsi::DetectOscillationType::IterationOnly;
    osc_prms.prevNormFactor = 0.2;

    solution_strategy->slvrType = bio::SolverType::Implicit;
    solution_strategy->slvrTolerance = 1E-6;
    solution_strategy->cnvgTolerance = 1E-6;
    solution_strategy->oscPrms = osc_prms;

    auto an = bio::createFiberRVEAnalysis(std::move(fiber_network), std::move(solution_strategy));
    //auto an = bio::createF
    dynamic_cast<bio::FiberRVEAnalysisSImplicit*>(an.get())->computeStiffnessMatrix();
    if(an == nullptr)
    {
      std::cerr<<"Something went wrong, the analysis doesn't exist!"<<std::endl;
      std::abort();
    }
    if( an->getK() == nullptr)
    {
      std::cerr<<"Something went wrong! The stiffness matrix doesn't exist!"<<std::endl;
      std::abort();
    }


    std::ofstream out("GlobalKMatrix.mtx");
    if (!out.is_open())
    {
      std::cerr << "Could not open file for writing!" << std::endl;
      std::abort();
    }
    // we write the matrix to output to apply output filtering
    std::cout << "Writing matrix" << std::endl;
    las::printSparskitMat(out, an->getK(), las::PrintType::mmarket, true);
    out.close();

    std::cout<<"Reading biotissue matrix"<<std::endl;
    std::ifstream in("GlobalKMatrix.mtx");
    if (!in.is_open())
    {
      std::cerr << "Could not open file for reading!" << std::endl;
      std::abort();
    }
    las::Mat * readMat = las::readSparskitMat(in, las::PrintType::mmarket);
    // multiply matrix by -1 since the sign convention for the test data is opposite
    auto smm = las::getScalarMatMult<las::MICRO_BACKEND>();
    smm->exec(-1, readMat, &readMat);
    in.close();
    std::size_t found = file_name.find_last_of("/");
    std::size_t found_dot = file_name.find_last_of(".");
    std::stringstream ss;
    ss << "./test_global_stiffness_data/" << file_name.substr(found+1, found_dot-found-1) << "_globalK.mtx";
    std::cout << "Reading Python Matrix: "<<ss.str() << std::endl;
    std::ifstream in2(ss.str());
    if (!in2.is_open())
    {
      std::cerr << "Could not open "<< ss.str() << " for reading!" << std::endl;
      std::abort();
    }
    las::Mat * readMat2 = las::readSparskitMat(in2, las::PrintType::mmarket);
    in2.close();
    std::cout << "Comparing matrix" << std::endl;
    assert(readMat);
    bool close = las::sparskitMatClose(readMat, readMat2, 1E-10, 1E-15);
    std::string comp =
         close ? "True" : "False";
    std::cout << "Matrix was close " << comp << std::endl;
    assert(close);
    mb->destroy(readMat);
    mb->destroy(readMat2);
  }
  delete mb;
  mb = NULL;
  amsi::freeAnalysis();
  return 0;
}
