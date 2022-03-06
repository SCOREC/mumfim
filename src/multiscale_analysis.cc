#include <amsiAnalysis.h>
#include <amsiMultiscale.h>
#include <amsiUtil.h>
#include <apfMDS.h>
#include <mumfim/microscale/bioMultiscaleRVEAnalysis.h>
#include <getopt.h>
#include <gmi_mesh.h>
#include <las.h>
#include <model_traits/ModelTraitsIO.h>
#include <Kokkos_Core.hpp>
#include <iostream>
#include "mumfim/macroscale/bioMultiscaleTissue.h"
#include "mumfim/macroscale/bioTissueMultiscaleAnalysis.h"
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
      << "  [-m, --mesh mesh_file]                    the mesh file (.sms)\n"
      << "  [-c, --case string]                       a string specifying the "
         "analysis case to run\n"
      << "  [-b, --balancing]                         model traits filename\n";
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
        {"balancing", required_argument, 0, 'b'},
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
int run_micro_fo(int & argc, char **& argv, MPI_Comm comm)
{
  Kokkos::ScopeGuard kokkos(argc, argv);
  las::initLAS(&argc, &argv, comm);
  int rnk = -1;
  MPI_Comm_rank(comm, &rnk);
  srand(8675309 + rnk);
  mumfim::MultiscaleRVEAnalysis rves;
  rves.init();
  rves.run();
  return 0;
}
int run_micro_fm(int &, char **&, MPI_Comm) { return 0; }
int run_macro(int & argc, char **& argv, MPI_Comm cm)
{
  amsi::initAnalysis(argc, argv, cm);
  int rnk = -1;
  MPI_Comm_rank(cm, &rnk);
  int result = 0;
  amsi::createDataDistribution(amsi::getLocal(), "micro_fo_data");
  gmi_register_mesh();
  apf::Mesh * mesh =
      apf::loadMdsMesh(model_filename.c_str(), mesh_filename.c_str());
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
      mesh, std::make_unique<mt::CategoryNode>(*case_traits), cm);
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
  amsi::freeAnalysis();
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
    amsi::initMultiscale(argc, argv, MPI_COMM_WORLD);
#ifdef LOGRUN
    amsi::Log execution_time = amsi::activateLog("execution_time");
#endif
    amsi::ControlService * control = amsi::ControlService::Instance();
    control->setScaleMain("macro", &run_macro);
    control->setScaleMain("micro_fo", &run_micro_fo);
    result = control->Execute(argc, argv);
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
    amsi::freeMultiscale();
  }
  else
    result++;
  return result;
}
