#include <amsiAnalysis.h>
#include <amsiMultiscale.h>
#include <amsiUtil.h>
#include <apfMDS.h>
#include <mumfim/microscale/MultiscaleRVEAnalysis.h>
#include <getopt.h>
#include <gmi_mesh.h>
#include <las.h>
#include <model_traits/ModelTraitsIO.h>
#include <Kokkos_Core.hpp>
#include <iostream>
#include "mumfim/macroscale/MultiscaleTissue.h"
#include "mumfim/macroscale/TissueMultiscaleAnalysis.h"
#include "amsi.h"
#include "mumfim/exceptions.h"
#if not defined(__APPLE__)
#include <cfenv>
#endif
bool file_exists(const std::string & name)
{
  std::ifstream f(name);
  return f.good();
}
void display_help_string()
{
  std::cout
      << "Usage: multiscale [OPTIONS]\n"
      << "  [-h, --help]                              display this help text\n"
      << "  [-g, --model model_file]                  the model file (.smd)\n"
      << "  [-m, --mesh mesh_file]                    the mesh file (.smb)\n"
      << "  [-c, --case string]                       a string specifying the "
         "analysis case to run\n"
      << "  [-b, --traits]                         model traits filename\n";
}
std::string model_filename("");
std::string mesh_filename("");
std::string analysis_case("");
std::string model_traits_filename("");
bool parse_options(int & argc, char **& argv)
{
  bool result = true;
  bool quit_loop = false;
  opterr = 0;
  while (!quit_loop)
  {
    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"model", required_argument, 0, 'g'},
        {"mesh", required_argument, 0, 'm'},
        {"traits", required_argument, 0, 'b'},
        {"case", required_argument, 0, 'c'}};
    int option_index = 0;
    int option =
        getopt_long(argc, argv, "hl:m:g:b:c:", long_options, &option_index);
    switch (option)
    {
      case 'h':
        display_help_string();
        result = false;
        break;
      case 'g':
        model_filename = optarg;
        break;
      case 'm':
        mesh_filename = optarg;
        break;
      case 'b':
        model_traits_filename = optarg;
        break;
      case 'c':
        analysis_case = optarg;
        break;
      case -1:
        // end of options
        quit_loop = true;
        break;
      default:
        break;
    }
  }
  optind = 0;  // reset global variable used by getopt_long
  opterr = 1;
  return result;
}
int run_micro_fo(int & argc,
                 char **& argv,
                 MPI_Comm comm,
                 amsi::Multiscale & multiscale)
{
  auto analysis_options = amsi::readAmsiOptions(model_traits_filename);
  analysis_options.analysis->use_petsc = false;
  amsi::Analysis analysis(analysis_options.analysis.value(), argc, argv, comm,
                          multiscale.getMPI());
  Kokkos::ScopeGuard kokkos(argc, argv);
  las::initLAS(&argc, &argv, comm);
  int rnk = -1;
  MPI_Comm_rank(comm, &rnk);
  srand(8675309 + rnk);
  mumfim::MultiscaleRVEAnalysis rves(multiscale);
  rves.init();
  rves.run();
  return 0;
}
int run_macro(int & argc,
              char **& argv,
              MPI_Comm cm,
              amsi::Multiscale & multiscale)
{
  auto analysis_options = amsi::readAmsiOptions(model_traits_filename);
  amsi::Analysis amsi_analysis(analysis_options.analysis.value(), argc, argv,
                               cm, multiscale.getMPI());
  int rnk = -1;
  MPI_Comm_rank(cm, &rnk);
  int result = 0;
  std::cerr << "Creating distribution\n";
  amsi::createDataDistribution(multiscale.getScaleManager()->getLocalTask(),
                               "micro_fo_data");
  std::cerr << "loading mesh\n";
  gmi_register_mesh();
  apf::Mesh * mesh =
      apf::loadMdsMesh(model_filename.c_str(), mesh_filename.c_str());
  std::cerr << "loading model traits\n";
  if (!file_exists(model_traits_filename))
  {
    std::cerr << "model traits file: " << model_traits_filename
              << " doesn't exist.\n";
    MPI_Abort(AMSI_COMM_WORLD, 1);
  }
  auto model_traits = mt::ReadFromFile<mt::YAML>(model_traits_filename);
  const auto * case_traits = model_traits->FindCase(analysis_case);
  if (case_traits == nullptr)
  {
    std::cerr << "\"" << analysis_case << "\" is not a valid case name.\n";
    MPI_Abort(AMSI_COMM_WORLD, 1);
  }
  mumfim::MultiscaleTissueAnalysis an(
      mesh, std::make_unique<mt::CategoryNode>(*case_traits), cm, amsi_analysis,
      multiscale);
  an.init();
  an.run();
#ifdef LOGRUN
  amsi::Log macro_stress = amsi::activateLog("macro_stresses");
  int rank = -1;
  MPI_Comm_rank(AMSI_COMM_SCALE, &rank);
  std::stringstream fname;
  fname << amsi::fs->getResultsDir() << "/macro_stress." << rank << ".log";
  {
    std::fstream strss_fs(fname.str().c_str(), std::ios::out | std::ios::app);
    amsi::flushToStream(macro_stress, strss_fs);
  }
#endif
  return result;
}
int main(int argc, char ** argv)
{
  int result = 0;
#if not defined(__APPLE__)
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif
  if (parse_options(argc, argv))
  {
    amsi::MPI mpi{argc, argv};
    auto amsi_options = amsi::readAmsiOptions(model_traits_filename);
    if (!amsi_options.multiscale)
    {
      throw mumfim::mumfim_error{"Multiscale analysis needs to be defined"};
    }
    // define the relations since they don't change and are required for
    // analysis
    amsi_options.multiscale->relations = {{"macro", "micro_fo"},
                                          {"micro_fo", "macro"}};
    const auto & scales = amsi_options.multiscale->scales;
    if (scales.size() != 2 ||
        (scales[0].name != "macro" && scales[1].name != "macro") ||
        (scales[0].name != "micro_fo" && scales[1].name != "micro_fo"))
    {
      throw mumfim::mumfim_error{
          "processors macro and micro_fo scales must be assigned in amsi "
          "input"};
    }
    if (scales[0].nprocs + scales[1].nprocs != mpi.getWorld().size())
    {
      throw mumfim::mumfim_error{
          "number of processors in macro and micro_fo scales don't match the "
          "number of processors in MPI world communicator"};
    }
    amsi::Multiscale multiscale{amsi_options.multiscale.value(), mpi};
    MPI_Barrier(MPI_COMM_WORLD);
#ifdef LOGRUN
    amsi::Log execution_time = amsi::activateLog("execution_time");
#endif
    auto * control = multiscale.getControlService();
    control->setScaleMain("macro", &run_macro);
    control->setScaleMain("micro_fo", &run_micro_fo);
    result = control->Execute(argc, argv, multiscale);
#ifdef LOGRUN
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // write the execution time log to an output file
    amsi::log(execution_time)
        << amsi::getElapsedTime(execution_time) << std::endl;
    std::stringstream filename;
    filename << amsi::fs->getResultsDir() << "execution_time." << rank
             << ".log";
    {
      std::fstream tm_fs(filename.str().c_str(), std::ios::out | std::ios::app);
      amsi::flushToStream(execution_time, tm_fs);
    }
    amsi::deleteLog(execution_time);
#endif
  }
  else
    result++;
  return result;
}
