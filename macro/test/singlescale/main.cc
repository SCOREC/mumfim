#include <amsiAnalysis.h>
#include <amsiUtil.h>
#include <apfMDS.h>
#include <fenv.h>
#include <getopt.h>
#include <gmi_mesh.h>
#include <model_traits/ModelTraitsIO.h>
#include <mpi.h>
#include <iostream>
#include <string>
#include "bioAnalysisIO.h"
#include "bioTissueAnalysis.h"
#ifdef HAVE_SIMMETRIX
#include <SimUtil.h>
#include <gmi_sim.h>
#endif
void display_help_string()
{
  std::cout
      << "Usage: singlescale [OPTIONS]\n"
      << "  [-h, --help]                              display this help text\n"
      << "  [-g, --model model_file]                  the model file (.dmg)\n"
      << "  [-c, --case string]                       a string specifying the "
         "analysis case to run\n"
      << "  [-m, --mesh mesh_file]                    the mesh file (.smb)\n"
      << "  [-b, --balancing]                         model traits filename\n";
}
std::string model_filename("");
std::string mesh_filename("");
std::string analysis_case("");
std::string model_traits_filename;
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
    };
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
  optind = 0;  // reset global variable used by getopt_long in case of later use
               // in the program
  opterr = 1;
  return result;
}
int main(int argc, char ** argv)
{
  int result = 0;
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
  if (parse_options(argc, argv))
  {
    amsi::initAnalysis(argc, argv, MPI_COMM_WORLD);
    int rnk = -1;
    MPI_Comm_rank(AMSI_COMM_WORLD, &rnk);
    if (rnk > 0) amsi::suppressOutput(std::cout);
    int sz = 0;
    MPI_Comm_size(AMSI_COMM_WORLD, &sz);
#ifdef HAVE_SIMMETRIX
    Sim_readLicenseFile(0);
    gmi_sim_start();
    gmi_register_sim();
#endif
    gmi_register_mesh();
    apf::Mesh * mesh =
        apf::loadMdsMesh(model_filename.c_str(), mesh_filename.c_str());
    auto model_traits = mt::ReadFromFile<mt::YAML>(model_traits_filename);
    const auto * case_traits = model_traits->FindCase(analysis_case);
    if (case_traits == nullptr)
    {
      std::cerr << "\"" << analysis_case << "\" is not a valid case name.\n";
      MPI_Abort(AMSI_COMM_WORLD, 1);
    }
    bio::TissueAnalysis an(mesh,
                           std::make_unique<mt::CategoryNode>(*case_traits),
                           AMSI_COMM_WORLD);
    an.init();
    an.run();
#ifdef HAVE_SIMMETRIX
    gmi_sim_stop();
    Sim_unregisterAllKeys();
#endif
    amsi::freeAnalysis();
  }
  return result;
}
