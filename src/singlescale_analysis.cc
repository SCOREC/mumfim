#include <amsiAnalysis.h>
#include <amsiUtil.h>
#include <apfMDS.h>
#include <fenv.h>
#include <getopt.h>
#include <gmi_mesh.h>
#include <lionPrint.h>
#include <model_traits/ModelTraitsIO.h>
#include <mpi.h>
#include <iostream>
#include <string>
#include "mumfim/macroscale/AnalysisIO.h"
#include "mumfim/macroscale/SinglescaleTissueAnalysis.h"
#include "amsiAnalysis.h"
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
      << "Usage: singlescale [OPTIONS]\n"
      << "  [-h, --help]                              display this help text\n"
      << "  [-g, --model model_file]                  the model file (.smd)\n"
      << "  [-c, --case string]                       a string specifying the "
         "analysis case to run\n"
      << "  [-m, --mesh mesh_file]                    the mesh file (.smb)\n"
      << "  [-b, --traits]                         model traits filename\n"
      << "  [-a, --amsi]                           amsi options filename\n";
}
std::string model_filename("");
std::string mesh_filename("");
std::string analysis_case("");
std::string model_traits_filename;
std::string amsi_options_filename;
// std::string vol_log("volume");
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
        {"case", required_argument, 0, 'c'},
        {"amsi", required_argument, 0, 'a'},
    };
    int option_index = 0;
    int option =
        getopt_long(argc, argv, "hl:m:g:b:c:a:", long_options, &option_index);
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
      case 'a':
        amsi_options_filename = optarg;
        break;
      case -1:
        // end of options
        quit_loop = true;
        break;
      default:
        break;
    }
  }
  optind = 0;  // reset global variable used by getopt_long in case of later use
               // in the program
  opterr = 1;
  return result;
}
int main(int argc, char ** argv)
{
  lion_set_verbosity(1);
  int result = 0;
#if not defined(__APPLE__)
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif
  if (parse_options(argc, argv))
  {
    amsi::MPI mpi(argc, argv, MPI_COMM_WORLD);
    auto amsi_options = amsi::readAmsiOptions(amsi_options_filename);
    amsi::Analysis amsi_analysis(amsi_options.analysis.value(), argc, argv,
                                 MPI_COMM_WORLD, mpi);
    int rnk = -1;
    MPI_Comm_rank(AMSI_COMM_WORLD, &rnk);
    if (rnk > 0) amsi::suppressOutput(std::cout);
    int sz = 0;
    MPI_Comm_size(AMSI_COMM_WORLD, &sz);
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
    mumfim::SinglescaleTissueAnalysis an(mesh,
                              std::make_unique<mt::CategoryNode>(*case_traits),
                              AMSI_COMM_WORLD, amsi_analysis);
    an.run();
  }
  return result;
}
